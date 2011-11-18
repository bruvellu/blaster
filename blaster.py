#!/usr/bin/python
# -*- coding: utf-8 -*-

'''BLASTing batch script.

Algorithm:
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

import getopt
import logging
import os
import subprocess
import sys

from Bio import SeqIO
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio.Blast.Applications import NcbiblastnCommandline




class Gene:
    '''Sequence instance of candidate gene.'''
    def __init__(self, path):
        self.path = path

class Database:
    '''Database instance for local BLAST.

    Input file should be a FASTA (for now, yes)?
    '''
    #TODO Check if input file is FASTA.
    #TODO Check if makeblastdb is properly installed.
    #TODO Check if conversion went ok.

    def __init__(self, path):
        #XXX Need to define path to makeblastdb for it to work...
        self.MAKEBLASTDB = '/home/nelas/Downloads/ncbi-blast-2.2.25+/bin/makeblastdb'

        # Database file path.
        self.path = path

        # If database is not ready, convert it.
        if not self.check_db():
            self.convert()

    def check_db(self):
        '''Check if database is ready or not.'''
        flagger = self.path + '.nsq'
        try:
            flagger_file = open(flagger)
            return True
        except IOError:
            logger.debug('Database %s needs to be converted.', self.path)
            return False

    def convert(self):
        '''Call makeblastdb to convert database.'''
        db_call = [self.MAKEBLASTDB, '-in', self.path, '-parse_seqids', '-dbtype', 'nucl']
        try:
            subprocess.call(db_call)
            logger.debug('Database %s was converted successfully!', 
                    self.path)
        except:
            logger.debug('We failed to convert the database %s :(', 
                    self.path)


def usage():
    '''Explanation for arguments.'''

    print
    print '  -c, --candidates \n\tFolder with `.fasta` or `.gb` files of candidate genes. One gene per file.'
    print '  -d, --database \n\tLocal database with new data (eg, transcriptome). If in FASTA format it will be converted to a database BLAST+ can understand using "makeblastdb" command, example: "makeblastdb -in bugula.fasta -parse_seqids -dbtype nucl"'
    print '  -b, --local-blasts \n\tFolder where local BLAST results will be written.'
    print '  -g, --genbank-blasts \n\tFolder where GenBank BLAST results will be written.'
    print

    #TODO Parameters for BLAST analysis.


def main(argv):
    '''
    '''
    # Default values.

    # Folder with candidate-genes.
    candidate_folder = 'candidates'

    # Folder with candidate-genes BLASTs.
    local_blasts = 'local_blasts'

    # Folder with candidate-genes GenBank BLASTs.
    genbank_blasts = 'genbank_blasts'

    # Database name.
    database = 'bugula.fasta'
    # The source file from Bugula neritina transcriptome from:
    # http://www.ncbi.nlm.nih.gov/sra/SRX015704

    # Parse arguments.
    try:
        opts, args = getopt.getopt(argv, 'hc:d:b:g:', [
            'help',
            'candidates=',
            'database=',
            'local-blasts=',
            'genbank-blasts=',
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
            candidate_folder = arg
        elif opt in ('-d', '--database'):
            database = arg
        elif opt in ('-b', '--local-blasts'):
            local_blasts = arg
        elif opt in ('-g', '--genbank-blasts'):
            genbank_blasts = arg

    # Print summary of arguments.
    logger.debug(
            'Arguments: candidates=%s, database=%s, local-blasts=%s, genbank-blasts=%s',
            candidate_folder, database, local_blasts, genbank_blasts
            )

    # Instantiate database.
    blast_db = Database(database)

    # Get candidate genes.
    try:
        candidates = os.listdir(candidate_folder)
    except OSError:
        logger.debug('Folder "%s" does not exist! Aborting...', candidate_folder)
        sys.exit(2)

    # Process candidate genes.
    for gene in candidates:
        gene_filepath = os.path.join(candidate_folder, gene)
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
        elif gene.endswith('.fasta'):
            continue
        else:
            logger.debug('File type not supported: %s', gene)

    # Get proper genes, now.
    candidates = os.listdir(candidate_folder)
    # Only read GenBank files.
    candidates = [file for file in candidates if file.endswith('.fasta')]

    # Check if local_blasts exists.
    if not os.path.isdir(local_blasts):
        os.mkdir(local_blasts)

    logger.info('%d genes to be BLASTed against %s database!', len(candidates),
            blast_db.path)

    # BLAST each gene against local database.
    for gene in candidates:
        gene_filepath = os.path.join(candidate_folder, gene)
        # Instantiate the BLASTn command.
        cline = NcbiblastnCommandline(
                query=gene_filepath, db=blast_db.path,
                out=os.path.join(local_blasts, gene.split('.')[0] + '.xml'),
                outfmt=5, # Export in XML
                #out=os.path.join(local_blasts, gene.split('.')[0]), # Export txt
                )

        # Execute BLAST.
        stdout, stderr = cline()

        logger.info('%s BLASTed!', gene)

    # Check if genbank_blasts exists.
    if not os.path.isdir(genbank_blasts):
        os.mkdir(genbank_blasts)

    # Get alignments to be now BLASTed against GenBank.
    sequences = os.listdir(local_blasts)
    sequences = [seq for seq in sequences if seq.endswith('.xml')]

    # Each gene might have a bunch of alignments so iterate.
    #TODO find a proper name... it is a blast output in xml.
    logger.info('Preparing to BLAST against GenBank...')
    for gene_name in sequences:
        parse_me = open(os.path.join(local_blasts, gene_name))
        genbank_dir = os.path.join(genbank_blasts, gene_name[:-4]) # strips extension
        # Check if each gene dir exists.
        if not os.path.isdir(genbank_dir):
            os.mkdir(genbank_dir)
        blast_records = NCBIXML.parse(parse_me)
        for blast_record in blast_records:
            E_VALUE_THRESH = 0.04
            for alignment in blast_record.alignments:
                for hsp in alignment.hsps:
                    if hsp.expect < E_VALUE_THRESH:
                        blast_file = os.path.join(genbank_dir, '%s.xml' % alignment.title)
                        try:
                            genbank_output = open(blast_file)
                            logger.debug('%s already BLASTed!', alignment.title)
                        except:
                            logger.info('Consulting NCBI for %s: %s', gene_name[:-4], hsp.sbjct)
                            handle = NCBIWWW.qblast('blastn', 'nr', hsp.sbjct)
                            handle_string = handle.read()
                            genbank_output = open(blast_file, 'w')
                            genbank_output.write(handle_string)
                            genbank_output.close()

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
