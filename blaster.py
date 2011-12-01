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
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from liblast import Sequence, Locus, blast

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

def usage():
    '''Explanation for arguments.'''

    print
    print 'USAGE: ./blaster.py -d Membranipora.fasta -b tblastn'
    print
    print '  -c, --candidates \n\tFolder with `.fasta` or `.gb` files of candidate genes. One gene per file.'
    print '  -d, --database \n\tLocal database with new data (eg, transcriptome).'
    print '  -b, --blast \n\tBLAST command (blastn, blastp, blastx, tblastn, tblastx).'
    print

    #TODO Ability to specify limit for candidate blast results to be analysed.
    #TODO Choose evalue threshold.

def main(argv):

    # Available BLASTs.
    blast_types = ['blastn', 'blastp', 'blastx', 'tblastn', 'tblastx']
    # Database name.
    database = None
    # BLAST command.
    blast_type = None

    # Folder with candidate-genes.
    candidates_folder = 'candidates'
    # Folder with candidate-genes BLASTs.
    candidates_results_folder = os.path.join(candidates_folder, 'results')

    # Parse arguments.
    try:
        opts, args = getopt.getopt(argv, 'hc:d:b:', [
            'help',
            'candidates=',
            'database=',
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
            database = os.path.abspath(arg)
        elif opt in ('-b', '--blast'):
            blast_type = arg

    # Print summary of arguments.
    logger.debug('Arguments: candidates=%s, database=%s, blast=%s', candidates_folder, database, blast_type)


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
    if not os.path.isdir(candidates_results_folder):
        os.mkdir(candidates_results_folder)

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

    # Print info before starting.
    logger.info('%d genes to be BLASTed against %s database!', len(candidates),
            database)


    # Main object to store genes and loci.
    genes = {}
    loci = {}

    # BLAST each gene against local database.
    for genefile in candidates:
        print
        logger.info('BLASTing %s against %s...', genefile, database)

        # Instantiate Sequence.
        gene_filepath = os.path.join(candidates_folder, genefile)
        candidate = Sequence(filepath=gene_filepath, database=database)

        # Define variables.
        output_filepath = os.path.join(candidates_results_folder, '%s.xml' % candidate.gene_name)
        candidate.blast_output = output_filepath

        # Only BLAST if needed.
        try:
            blastfile = open(candidate.blast_output)
            blastfile.close()
        except:
            arguments = {
                'query': candidate.filepath,
                'db': database,
                'out': candidate.blast_output,
                'outfmt': 5, # Export in XML
                }
            # Execute BLAST command.
            blast(blast_type, arguments)

        logger.info('Parsing XML...')
        print '\nCandidate >> Gene: %s, ID:%s, Ref: %s, Description: %s' % (
                candidate.gene_name, candidate.gene_id,
                candidate.ref, candidate.description
                )
        candidate.parse_blast()

        #DO LOCUS STUFF
        for locus_id in candidate.loci.keys():
            # If locus was already specified.
            if locus_id in loci.keys():
                locus = loci[locus_id]
                locus.update_candidates(candidate)
            else:
                # Instantiate Locus object.
                locus = Locus(locus_id, candidate, database)

                # Reciprocal BLAST querying locus sequence against human database.
                try:
                    blastfile = open(locus.reverse_blast_output)
                    blastfile.close()
                except:
                    #XXX Find a better way to overcome problem when path has a "|".
                    reverse_args = {
                            'query': locus.filepath.replace('|', '\|'),
                            'db': '/home/nelas/Biologia/Doutorado/genomic_databases/human_protein.fa',
                            'out': locus.reverse_blast_output.replace('|', '\|'),
                            'outfmt': 5,
                            }
                    if blast_type == 'tblastn':
                        blast('blastx', reverse_args)
                    elif blast_type == 'blastp':
                        blast('blastp', reverse_args)

                # Parse reverse BLAST output.
                locus.parse_blast()

                # Add locus to main list.
                loci[locus_id] = locus

        # Add to main dictionary.
        genes[candidate.gene_name] = candidate


    # Print loci equivalent to candidate genes.
    logger.info('Creating result files...')

    fasta_loci = []
    output_file = open('results.txt', 'w')

    #TODO Print date and variables used to results.txt.
    for locus_id, locus in loci.iteritems():
        fasta_loci.append(SeqRecord(Seq(locus.sequence), id=locus_id, 
            description='related to: %s' % ', '.join(locus.candidates.keys())))
        print
        print locus_id, locus.candidates.keys()
        output_file.write('%s\t\t(candidates: %s)\n\n' % (locus_id, 
            ', '.join(locus.candidates.keys())))
        try:
            first = locus.reciprocals[0]
        except:
            print '\tNo reciprocals.'
        current_gene_id = first.gene_id
        n_genes = 0
        output_file.write('\tgene\t\tid\t\taccession\t\te-value\n')
        for sequence in locus.reciprocals:
            # Limits to 5 most similar genes.
            if sequence.gene_id == current_gene_id or n_genes < 5:
                print '\t%s, %s, %s, %s' % (
                        sequence.gene_name, sequence.gene_id,
                        sequence.ref, sequence.evalue
                        )

                output_file.write('\t%s\t\t%s\t\t%s\t\t%.2e' % (
                        sequence.gene_name, sequence.gene_id,
                        sequence.ref, sequence.evalue
                        ))
                accession_numbers = [seq.ref for seq in locus.candidates.values()]
                if sequence.ref in accession_numbers:
                    output_file.write(' <<')
                output_file.write('\n')
            else:
                break
            current_gene_id = sequence.gene_id
            n_genes += 1
        output_file.write('\n\n')
    print
    output_file.close()

    wrote_fasta = SeqIO.write(fasta_loci, 'results.fa', 'fasta')

    logger.info('%d records written to results.fa', wrote_fasta)
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
