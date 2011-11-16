#!/usr/bin/python
# -*- coding: utf-8 -*-

'''BLASTing batch script.

Algorhythm:
    - Take several genes of interest.
    - BLAST them locally agains a transcriptome assembly.
    - Filter and process results to extract relevant sequences.
    - BLAST resulting sequences against GenBank.
    - Retrieve results for further analysis.

Notes:
    - Read files with sequences from a folder (FASTA or GenBank format). Could also be a list of IDs or
      an pickled object, etc...

Dependencies:
    - Python 2.7.1
    - Biopython 1.56
    - BLAST+ 2.2.25

Configuration (Ubuntu):
    - Regular packages for Python and Biopython.
    - Add BLAST+ commands to path putting this code to .bashrc:

        PATH=$PATH:/home/nelas/Downloads/ncbi-blast-2.2.25+/bin
        export PATH
'''

import os

from Bio import SeqIO
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio.Blast.Applications import NcbiblastnCommandline

# Folder with candidate-genes.
CANDIDATE_FOLDER = 'candidates'

# Folder with candidate-genes BLASTs.
CANDIDATE_BLASTS = 'local_blasts'

# Folder with candidate-genes BLASTs.
GENBANK_BLASTS = 'genbank_blasts'

# Database name.
DATABASE = 'bugula.fasta'
# The source file from Bugula neritina transcriptome was downloaded from here:
# http://www.ncbi.nlm.nih.gov/sra/SRX015704

# Note: the original FASTA file needs to be converted to a BLASTable database.
# This can be done via command line and preferably before the script is run:
# makeblastdb -in bugula.fasta -parse_seqids -dbtype nucl

# Deprecated way of doing it, just for reference:

# Format a FASTA file to a readable database.
# formatdb -i bugula.fasta -p F -o
# BLAST against nucleotide with XML output:
# blast2 -p blastn -d bugula.fasta -o blastbugula.xml -i candidates/BnFoxA.gb -m 7

# Build list with path to candidate-genes.
candidates = os.listdir(CANDIDATE_FOLDER)
# Only read GenBank files.
candidates = [file for file in candidates if file.endswith('.gb')]

print '\n%d genes to be BLASTed against %s database!' % (len(candidates),
        DATABASE)

for gene_file in candidates:
    # Define filepath.
    gene_filepath = os.path.join(CANDIDATE_FOLDER, gene_file)

    # Create record object.
    record = SeqIO.read(gene_filepath, 'genbank')

    # Create temporary FASTA file, if not available.
    try:
        fasta_file = open(gene_filepath.split('.')[0] + '.fasta', 'r')
    except IOError:
        fasta_file = open(gene_filepath.split('.')[0] + '.fasta', 'w')
        fasta_file.write(record.format('fasta'))
        fasta_file.close()
        fasta_file = open(gene_filepath.split('.')[0] + '.fasta', 'r')

    # Instantiate the BLASTn command.
    cline = NcbiblastnCommandline(
            query=fasta_file.name, db=DATABASE,
            out=os.path.join(
                CANDIDATE_BLASTS, gene_file.split('.')[0]),
                #CANDIDATE_BLASTS, gene_file.split('.')[0] + '.xml'),
            #outfmt=5 # Export in XML
            )

    # Execute BLAST.
    stdout, stderr = cline()

    # Close file.
    fasta_file.close()

    print '\t' + gene_file

#TODO Process local BLASTs

#TODO BLAST results against GenBank

#TODO Process remote BLASTs

print 'Done, bye!\n'
