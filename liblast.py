#!/usr/bin/python
# -*- coding: utf-8 -*-

'''Library for BLASTer script.'''

import os
import pickle
from Bio import Entrez, SeqIO
from Bio.Blast import NCBIXML

class Sequence(object):
    '''Represents a known sequence.

    filepath: path to fasta file.
    ref: sequence reference number.
    description: long description.
    sequence: the sequence itself.
    gene_id: NCBI gene identification.
    gene_name: NCBI gene name.
    blast_output: path to .xml file.
    loci: list of loci with recurrency?.
    '''
    Entrez.email = 'organelas@gmail.com'
    cache_folder = 'seqs'
    #TODO add __str__ to classes.

    def __init__(self, filepath=None, ref=None):
        self.loci = []
        self.gene_name = ''
        if filepath:
            self.filepath = filepath
            self.parse_fasta(filepath)
            self.get_gene_id()
        else:
            if not os.path.isdir(self.cache_folder):
                os.mkdir(self.cache_folder)
            cache_path = os.path.join(self.cache_folder, ref)
            try:
                cache_file = open(cache_path, 'rb')
                cached = pickle.load(cache_file)
                self.load_cache(cached)
                cache_file.close()
            except:
                self.fetch(ref)
                cache_file = open(cache_path, 'wb')
                pickle.dump(self, cache_file)
                cache_file.close()

    def load_cache(self, cached):
        '''Load attributes saved in cache files.'''
        self.description = cached.description
        self.ref = cached.ref
        self.sequence = cached.sequence
        self.gene_name = cached.gene_name
        self.gene_id = cached.gene_id

    def parse_fasta(self, filepath):
        '''Parse data from FASTA file.'''
        record = SeqIO.read(filepath, 'fasta')
        # Define attributes.
        self.description = record.description.split('|')[-1]
        self.ref = record.id.split('|')[3]
        self.sequence = str(record.seq)
        self.gene_name = filepath.split('/')[-1][:-3]

    def get_gene_id(self):
        '''Get gene id from NCBI.'''
        handle = Entrez.esearch(db='gene', term=self.ref)
        record = Entrez.read(handle)
        self.gene_id = record['IdList'][0]

    def fetch(self, ref):
        '''Fetch data from NCBI.'''
        handle = Entrez.efetch(db='protein', id=ref, rettype='gb')
        record = SeqIO.read(handle, 'genbank')
        # Define attributes.
        self.description = record.description
        self.ref = record.id
        self.sequence = str(record.seq)
        self.get_gene_id()

class Locus(object):
    '''Represents a new sequence.

    candidate: instance of Candidate.
    description: long description.
    sequence: the sequence itself.
    score: from Candidate BLAST.
    evalue: from Candidate BLAST.
    reverse_blast_output: reverse BLAST to human protein db.
    '''
    reverses = []

    def __init__(self, candidate, filepath, description, score, evalue, sequence):
        self.candidate = candidate
        self.filepath = filepath
        self.description = description
        self.score = int(score)
        self.evalue = evalue
        self.sequence = sequence
        self.equivalent = False

        self.write_fasta()

    def write_fasta(self):
        '''Write FASTA file to filepath.'''
        locus_file = open(self.filepath, 'w')
        locus_file.write('>%s\n%s' % (self.description, self.sequence))
        locus_file.close()

    def parse_blast(self):
        '''Parse reverse BLAST.'''
        blast_records = NCBIXML.parse(open(self.reverse_blast_output))
        for blast_record in blast_records:
            E_VALUE_THRESH = 0.001
            for alignment in blast_record.alignments:
                for hsp in alignment.hsps:
                    if hsp.expect < E_VALUE_THRESH:
                        ref = alignment.title.split('|')[1]
                        self.reverses.append(Sequence(ref=ref))

    def process(self):
        '''Define if locus is equivalent to the gene.

        True: if it matches the same gene_id and ref from the candidate gene.
        '''
        #XXX Will this rule out high scoring sequences?
        for sequence in self.reverses:
            if sequence.gene_id == self.candidate.gene_id:
                if sequence.ref == self.candidate.ref:
                    self.equivalent = True
                    break
