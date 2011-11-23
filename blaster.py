#!/usr/bin/python
# -*- coding: utf-8 -*-

'''BLASTing batch script.

Algorithm:
    - Take several genes of interest.
    - BLAST them locally against a transcriptome assembly.
    - Filter and process results to extract relevant sequences.

Notes:
    - Read files with sequences from a folder (FASTA format).

Dependencies:
    - Python 2.7.1
    - Biopython 1.56
    - BLAST+ 2.2.25

Configuration (Linux):
    - Regular packages for Python and Biopython.
    - Add BLAST+ commands to path putting this code to .bashrc:

        PATH=$PATH:/home/nelas/Downloads/ncbi-blast-2.2.25+/bin
        export PATH

'''

import getopt
import logging
import os
import subprocess
import sys

from Bio import SeqIO
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio.Blast.Applications import NcbiblastnCommandline, NcbiblastpCommandline, NcbiblastxCommandline, NcbitblastnCommandline, NcbitblastxCommandline
from liblast import Sequence, Locus

def prepare(candidates, candidates_folder):
    '''Check candidate gene files (convert to FASTA, if needed).'''
    for gene in candidates:
        gene_filepath = os.path.join(candidates_folder, gene)

        if gene.endswith('.gb'):
            # Create record object.
            record = SeqIO.read(gene_filepath, 'genbank')
            # Create FASTA file.
            try:
                fasta_file = open(gene_filepath.split('.')[0] + '.fasta', 'r')
            except IOError:
                fasta_file = open(gene_filepath.split('.')[0] + '.fasta', 'w')
                fasta_file.write(record.format('fasta'))
                fasta_file.close()
        elif gene.endswith('.fa'):
            continue
        else:
            logger.debug('File type not supported: %s', gene)

def print_loci(genes):
    '''Print loci equivalent to candidate genes.'''
    print 'Creating final file...'
    final_file = open('results.fa', 'w')
    for name, gene in genes.iteritems():
        for seq in gene.loci:
            final_file.write('>%s\n' % seq.description)
            final_file.write('%s\n\n' % seq.sequence)
    final_file.close()

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
    logger.info('%s BLASTed!', arguments['query'])

def usage():
    '''Explanation for arguments.'''

    print
    print 'USAGE: ./blaster.py -d Membranipora.fasta -b tblastn'
    print '  -c, --candidates \n\tFolder with `.fasta` or `.gb` files of candidate genes. One gene per file.'
    print '  -d, --database \n\tLocal database with new data (eg, transcriptome).'
    print '  -r, --results \n\tFolder where BLAST results will be written.'
    print '  -b, --blast \n\tBLAST command (blastn, blastp, blastx, tblastn, tblastx).'
    print

def main(argv):
    # Folder with candidate-genes.
    candidates_folder = 'candidates'
    # Folder with candidate-genes BLASTs.
    results_folder = os.path.join(candidates_folder, 'results')

    # Reverse BLAST folder
    reverse_folder = 'reverse'
    # Folder with candidate-genes BLASTs.
    reverse_results_folder = os.path.join(reverse_folder, 'results')

    # Database name.
    database = None

    # BLAST command.
    blast_type = None

    # Available BLASTs.
    blast_types = ['blastn', 'blastp', 'blastx', 'tblastn', 'tblastx']

    # Parse arguments.
    try:
        opts, args = getopt.getopt(argv, 'hc:d:r:b:', [
            'help',
            'candidates=',
            'database=',
            'reverse=',
            'blast=',
            ])
    except getopt.GetoptError:
        usage()
        logger.critical('Argument error: "%s". Aborting...',
                ' '.join(argv))
        sys.exit(2)

    # Argument handler.
    for opt, arg in opts:
        if opt in ('-h', '--help'):
            usage()
            sys.exit()
        elif opt in ('-c', '--candidates'):
            candidates_folder = arg
        elif opt in ('-d', '--database'):
            database = arg
        elif opt in ('-b', '--blast'):
            blast_type = arg

    # Print summary of arguments.
    logger.debug('Arguments: candidates=%s, database=%s', candidates_folder, database)

    # Check if BLAST command was specified.
    if not blast_type:
        logger.critical('BLAST command was not specified (use "-b"). Aborting...')
        sys.exit(2)
    else:
        if not blast_type in blast_types:
            logger.critical('Unknown BLAST command: %s', blast_type)
            logger.info('Available BLAST commands: %s', ', '.join(blast_types))
            sys.exit(2)

    # Check if results exists.
    if not os.path.isdir(results_folder):
        os.mkdir(results_folder)
    # Check if results exists.
    if not os.path.isdir(reverse_folder):
        os.mkdir(reverse_folder)
    # Check if results exists.
    if not os.path.isdir(reverse_results_folder):
        os.mkdir(reverse_results_folder)

    # Get candidate genes.
    try:
        candidates = os.listdir(candidates_folder)
    except OSError:
        logger.critical('Folder "%s" does not exist! Aborting...', candidates_folder)
        sys.exit(2)

    # Process candidate genes.
    prepare(candidates, candidates_folder)
    # Get proper genes, now.
    candidates = os.listdir(candidates_folder)
    # Only read FASTA files.
    candidates = [file for file in candidates if file.endswith('.fa')]

    logger.info('%d genes to be BLASTed against %s database!', len(candidates),
            database)

    # Main objects to store genes.
    genes = {}

    # BLAST each gene against local database.
    for genefile in candidates:
        print genefile
        gene_name = genefile[:-3]
        gene_filepath = os.path.join(candidates_folder, genefile)
        output_filepath = os.path.join(results_folder, gene_name + '.xml')

        # Only BLAST if needed.
        try:
            blastfile = open(output_filepath)
            blastfile.close()
        except:
            arguments = {
                'query': gene_filepath,
                'db': database,
                'out': output_filepath,
                'outfmt': 5, # Export in XML
                }
            # Execute BLAST command.
            blast(blast_type, arguments)

        # Instantiate Sequence.
        candidate = Sequence(filepath=gene_filepath)
        candidate.gene_name = gene_name
        candidate.blast_output = output_filepath

        print 'Parsing xml...'
        # BLAST records.
        blast_records = NCBIXML.parse(open(candidate.blast_output))

        # Limit for the first locus.
        n = 0

        # Iterate over BLAST results from database.
        for blast_record in blast_records:
            E_VALUE_THRESH = 0.001
            for alignment in blast_record.alignments:
                    for hsp in alignment.hsps:
                        if hsp.expect < E_VALUE_THRESH:
                            locus_filename = '%s_%s_Score_%s' % (gene_name, alignment.title, int(hsp.score))
                            locus_filename = locus_filename.replace('/', '-')
                            locus_filepath = os.path.join(reverse_folder, '%s.fa' % locus_filename)

                            print 'Instantiating %s...' % locus_filename
                            # Instantiate Locus object.
                            locus = Locus(candidate, locus_filepath, locus_filename, hsp.score, hsp.expect, hsp.sbjct)

                            #locus_file = open(locus_filepath)
                            locus_output_filepath = os.path.join(reverse_results_folder, '%s' % locus_filename)
                            try:
                                blastfile = open('%s.xml' % locus_output_filepath)
                                blastfile.close()
                            except:
                                reverse_args = {
                                        'query': locus_filepath,
                                        'db': 'human_protein.fa',
                                        'out': '%s.xml' % locus_output_filepath,
                                        'outfmt': 5,
                                        }

                                blast('blastp', reverse_args)

                            locus.reverse_blast_output = '%s.xml' % locus_output_filepath

                            locus.parse_blast()
                            locus.process()
                            print
                            print candidate.ref, locus.description
                            if locus.equivalent:
                                candidate.loci.append(locus)
                                print [gene.description for gene in candidate.loci]
                            print


                            n += 1

        # Add to main dictionary.
        genes[gene_name] = candidate

    print_loci(genes)

    logger.info('Done, bye!')


if __name__ == '__main__':
    # Create logger.
    logger = logging.getLogger('blaster')
    logger.setLevel(logging.DEBUG)
    logger.propagate = False
    # Define message format.
    formatter = logging.Formatter('[%(levelname)s] %(asctime)s @ %(module)s %(funcName)s (l%(lineno)d): %(message)s')

    # Create logger console handler.
    console_handler = logging.StreamHandler()
    console_handler.setLevel(logging.DEBUG)
    # Format console output.
    console_handler.setFormatter(formatter)
    # Add console to logger.
    logger.addHandler(console_handler)

    # Create logger file handler.
    file_handler = logging.FileHandler('log_blaster.log')
    file_handler.setLevel(logging.DEBUG)
    # Format file output.
    file_handler.setFormatter(formatter)
    # Add file to logger.
    logger.addHandler(file_handler)

    # Initiates main function.
    logger.info('BLASTer is running...')
    main(sys.argv[1:])
