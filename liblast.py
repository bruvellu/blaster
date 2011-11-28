#!/usr/bin/python
# -*- coding: utf-8 -*-

'''Library for BLASTer script.'''

import os
import pickle
import re
from Bio import Entrez, SeqIO
from Bio.Blast import NCBIXML
from Bio.Blast.Applications import NcbiblastnCommandline, NcbiblastpCommandline, NcbiblastxCommandline, NcbitblastnCommandline, NcbitblastxCommandline


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
    EVALUE_THRESH = 0.001

    # Cache folder.
    cache_folder = 'seqs'
    # Check if cache folder exists.
    if not os.path.isdir(cache_folder):
        os.mkdir(cache_folder)

    #TODO add __str__ to classes.

    def __init__(self, filepath=None, ref=None, database=None):
        #TODO Make limit a variable variable.
        self.limit = 1
        self.loci = {}
        self.gene_name = ''
        self.database = database
        if filepath:
            self.filepath = filepath
            self.parse_fasta(self.filepath)
            self.set_gene_id()
            self.set_gene_name()
        else:
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
        #XXX Better parse ref!
        self.ref = record.id.split('|')[3]
        self.sequence = str(record.seq)
        self.gene_name = filepath.split('/')[-1][:-3]

    def set_gene_id(self):
        '''Get and set gene id from NCBI.'''
        #import pdb; pdb.set_trace()
        #FIXME Handle error when not returning id.
        handle = Entrez.esearch(db='gene', term=self.ref)
        record = Entrez.read(handle)
        if record['IdList']:
            self.gene_id = record['IdList'][0]
        else:
            self.gene_id = None

    def set_gene_name(self):
        '''Get and set gene name from filepath.'''
        filepath = os.path.basename(self.filepath)
        self.gene_name = os.path.splitext(filepath)[0]

    def fetch(self, ref):
        '''Fetch data from NCBI.'''
        handle = Entrez.efetch(db='protein', id=ref, rettype='gb')
        record = SeqIO.read(handle, 'genbank')
        #TODO maybe just save file locally...
        handle = Entrez.efetch(db='protein', id=ref, rettype='gb')
        handle_string = handle.read()
        # Define attributes.
        #FIXME See if this is ok!
        self.description = record.description
        self.ref = record.id
        self.sequence = str(record.seq)
        self.set_gene_id()
        # Find gene name by regular expression.
        try:
            pattern = re.compile(r'gene="(.+)"')
            reobj = pattern.search(handle_string)
            self.gene_name = reobj.group(1)
        except:
            self.gene_name = 'Not found'
        # To extract gene description, name and other info from xml parsing.
        #>>> record['Bioseq-set_seq-set'][0]['Seq-entry_set']['Bioseq-set']['Bioseq-set_seq-set'][0]['Seq-entry_seq']['Bioseq']['Bioseq_annot'][0]['Seq-annot_data']['Seq-annot_data_ftable'][1]['Seq-feat_data']['SeqFeatData']['SeqFeatData_gene']['Gene-ref']
        #>>> {u'Gene-ref_desc': 'DEAD (Asp-Glu-Ala-Asp) box polypeptide 4', u'Gene-ref_syn': ['VASA'], u'Gene-ref_locus': 'DDX4'}

    def parse_blast(self):
        '''Parse BLAST output file.'''
        blast_records = NCBIXML.parse(open(self.blast_output))
        n = 0
        # Iterate over BLAST results from database.
        for blast_record in blast_records:
            for alignment in blast_record.alignments:
                for hsp in alignment.hsps:
                    if n < self.limit:
                        if hsp.expect < self.EVALUE_THRESH:

                            print 'Instantiating >> %s' % alignment.title
                            locus_id = alignment.title.split()[0]

                            # Save locus to dictionary if key does not exists.
                            #NOTE This means that only the first hsp found in the locus is added,
                            # while the other ones are skipped
                            if not locus_id in self.loci.keys():
                                self.loci[locus_id] = {
                                        'score': hsp.score,
                                        'evalue': hsp.expect
                                        }

                            n += 1


class Locus(object):
    '''Represents a new sequence.

    candidate: instance of Candidate.
    description: long description.
    sequence: the sequence itself.
    score: from Candidate BLAST.
    evalue: from Candidate BLAST.
    reverse_blast_output: reverse BLAST to human protein db.
    '''
    EVALUE_THRESH = 0.001

    # Reverse BLAST folder.
    reverse_folder = 'reverse'

    def __init__(self, id, candidate, database):
        # Linkback to candidates.
        self.candidates = {}

        self.reciprocals = []

        self.update_candidates(candidate)

        self.reciprocal = False
        self.rank = 0
        self.id = id
        self.database = database

        # Create filename before the rest.
        self.set_filename()

        self.reverse_results_folder = os.path.join(self.reverse_folder, 'results')

        # Check if reverse folder exists.
        if not os.path.isdir(self.reverse_folder):
            os.mkdir(self.reverse_folder)
        # Check if reverse results folder exists.
        if not os.path.isdir(self.reverse_results_folder):
            os.mkdir(self.reverse_results_folder)

        self.set_filepath()
        self.set_blast_output()
        self.set_sequence()

        self.write_fasta()

    def update_candidates(self, candidate):
        '''Update list of linked candidate genes.'''
        if not candidate.gene_name in self.candidates.keys():
            self.candidates[candidate.gene_name] = candidate

    def set_sequence(self):
        '''Get and set sequence from database.'''
        parsed = SeqIO.parse(self.database, 'fasta')
        for locus in parsed:
            if locus.id == self.id:
                self.sequence = str(locus.seq)
        if not self.sequence:
            print 'SEQUENCE NOT SET! Check the IDs.'

    def set_blast_output(self):
        '''Build filepath for reverse blast output.'''
        blast_output = os.path.join(self.reverse_results_folder, '%s.xml' % self.filename)
        self.reverse_blast_output = blast_output

    def set_filepath(self):
        '''Build filepath from folder and filename.'''
        locus_filepath = os.path.join(self.reverse_folder, '%s.fa' % self.filename)
        self.filepath = locus_filepath

    def set_filename(self):
        '''Extract filename from ID.'''
        filename = '_'.join(self.id.split('_')[:2])
        self.filename = filename

    def write_fasta(self):
        '''Write FASTA file to filepath.'''
        #TODO Insert line breaks for sequences.
        locus_file = open(self.filepath, 'w')
        locus_file.write('>%s\n%s' % (self.id, self.sequence))
        locus_file.close()

    def parse_blast(self):
        '''Parse reciprocal BLAST.'''
        print 'Parsing reciprocal BLAST: %s' % self.filename
        blast_records = NCBIXML.parse(open(self.reverse_blast_output))
        for blast_record in blast_records:
            for alignment in blast_record.alignments:
                for hsp in alignment.hsps:
                    if hsp.expect < self.EVALUE_THRESH:
                        #XXX Better parse ref!
                        ref = alignment.title.split('|')[3]
                        #print 'Creating object: %s' % ref
                        sequence = Sequence(ref=ref)
                        sequence.evalue = hsp.expect
                        sequence.score = hsp.score
                        self.reciprocals.append(sequence)
        # Process results.
        #self.process()

    def process(self):
        '''Define if locus is reciprocal to the gene.

        True: if it matches the same gene_id and ref from the candidate gene.
        '''
        #TODO Handle when gene_id is None?
        #TODO import pdb; pdb.set_trace()
        #TODO print to csv!
        rank = 0
        print '\nrank\tgene\tref\t\thit'
        for sequence in self.reciprocals:
            if sequence.gene_id == current_gene_id:
                #DOSTUFF
                continue
            else:
                break
            current_gene_id = sequence.gene_id

            rank += 1
            #if sequence.gene_id == self.candidate.gene_id:
            #    if sequence.ref == self.candidate.ref:
            #        print '%d\t%s\t%s\t%s' % (rank, sequence.gene_id, sequence.ref, '+')
            #        print '\nFound reciprocal! Gene: %s, Ref: %s\n' % (sequence.gene_id, sequence.ref)
            #        self.reciprocal = True
            #        self.rank = rank
            #        break
            print '%d\t%s\t%s\t%s' % (rank, sequence.gene_id, sequence.ref, '-')


def blast(blast_type, arguments):
    '''Execute BLAST command.'''
    # Instantiate the BLAST command.
    if blast_type == 'blastn':
        cline = NcbiblastnCommandline(**arguments)
    elif blast_type == 'blastp':
        cline = NcbiblastpCommandline(**arguments)
    elif blast_type == 'blastx':
        cline = NcbiblastxCommandline(**arguments)
    elif blast_type == 'tblastn':
        cline = NcbitblastnCommandline(**arguments)
    elif blast_type == 'tblastx':
        cline = NcbitblastnCommandline(**arguments)

    # Execute BLAST.
    stdout, stderr = cline()
    print '%s BLASTed!' % arguments['query']


