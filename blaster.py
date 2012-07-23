#!/usr/bin/python
# -*- coding: utf-8 -*-

'''BLASTing batch script.

Dependencies:
    - Python 2.7.1
    - Biopython 1.56
    - BLAST+ 2.2.25

Configuration (Linux):
    - Regular packages for Python and Biopython.
    - Add BLAST+ commands to path putting this code to .bashrc:

        PATH=$PATH:/home/nelas/Downloads/ncbi-blast-2.2.25+/bin
        export PATH

See README file for details or execute the command for usage and arguments:

    python blaster.py -h

'''

import getopt
import logging
import os
import pickle
import re
import sys

from Bio import Entrez, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Blast import NCBIWWW, NCBIXML
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
    EVALUE_THRESH = 0.001

    # Cache folder.
    cache_folder = 'seqs'
    # Check if cache folder exists.
    if not os.path.isdir(cache_folder):
        os.mkdir(cache_folder)

    #TODO add __str__ to classes.

    def __init__(self, filepath=None, ref=None, database=None):
        self.limit = 3
        self.loci = {}
        self.gene_name = ''
        self.database = database
        if filepath:
            self.filepath = filepath
            self.parse_genbank(self.filepath)
        else:
            cache_filename = ref + '.gb'
            cache_path = os.path.join(self.cache_folder, cache_filename)
            # Use efetcher function here.
            try:
                self.parse_genbank(cache_path)
            except:
                handle_string = efetcher(ref)
                cache_file = open(cache_path, 'wb')
                cache_file.write(handle_string)
                cache_file.close()
                self.parse_genbank(cache_path)

    def load_cache(self, cached):
        '''Load attributes saved in cache files.'''
        self.description = cached.description
        self.ref = cached.ref
        self.sequence = cached.sequence
        self.gene_name = cached.gene_name
        self.gene_id = cached.gene_id

    def parse_genbank(self, filepath):
        '''Parse data from GENBANK file.'''
        record = SeqIO.read(filepath, 'genbank')
        # Define class attributes.
        self.description = record.description
        self.ref = record.id
        self.sequence = str(record.seq)
        self.gene_name, self.gene_id = self.get_gene_name_id(record)
        self.organism = record.description.split('[')[1].split(']')[0]
        #self.set_reciprocal_db(self.organism)

    def parse_fasta(self, filepath):
        '''Parse data from FASTA file.'''
        record = SeqIO.read(filepath, 'fasta')
        # Define attributes.
        self.description = record.description.split('|')[-1]
        self.ref = record.id.split('|')[3]
        self.sequence = str(record.seq)
        self.gene_name = filepath.split('/')[-1][:-3]
        self.organism = record.description.split('[')[1].split(']')[0]
        #self.set_reciprocal_db(self.organism)

    def get_gene_name_id(self, record):
        '''Return the gene name and id from a SeqRecord.'''
        cds = [feature for feature in record.features if feature.type == 'CDS'][0]
        try:
            gene_name = cds.qualifiers.get('gene')[0]
            gene_db_xref = cds.qualifiers.get('db_xref')
            for xref in gene_db_xref:
                if xref.startswith('GeneID'):
                    gene_id = xref.split(':')[-1]
            return gene_name, gene_id
        except:
            # No CDS in the record.
            logger.critical('No CDS was found in %s', record.id)
            return None, None

    def set_reciprocal_db(self):
        '''Set the reciprocal database according to candidate gene organism.'''
        #FIXME Instead of doing a local blast, try to do a online one
        # restricting to the organisms refseq.
        # handle = NCBIWWW.qblast("tblastn", "nr", "NP_077726.1", entrez_query="drosophila melanogaster[Organism]")
        #XXX Keep an option to run a local reciprocal BLAST for the databases
        # available (it is much faster).

        # Locally define organism.
        organism = self.organism
        base_path = '/home/nelas/Biologia/Doutorado/databases/'
        reciprocals = {
                'Homo sapiens': 'human_protein.fa',
                'Drosophila melanogaster': 'drosophila.fa',
                'Danio rerio': 'zebrafish.fa',
                'Mus musculus': 'mouse.fa',
                'Strongylocentrotus purpuratus': 'urchin.fa',
                'Ciona intestinalis': 'ciona.fa',
                'Prostheceraeus vittatus': 'prostheceraeus.fa',
                }
        try:
            self.reciprocal_db = os.path.join(base_path, reciprocals[organism])
        except:
            #XXX Set a notice to install a local database.
            self.reciprocal_db = None

    def set_gene_id(self):
        '''Get and set gene id from NCBI.'''
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
        # Fetch entry and save handle to string.
        handle = Entrez.efetch(db='protein', id=ref, rettype='gp', retmode='txt')
        handle_string = handle.read()

        # Create temporary file and write data.
        temp_gene = '/tmp/temp_gene'
        f = open(temp_gene, 'w')
        f.write(handle_string)
        f.close()

        # Create record by reading temp file.
        record = SeqIO.read(open(temp_gene, 'rU'), 'genbank')

        # Define attributes.
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
                                        'evalue': hsp.expect,
                                        'frame': hsp.frame[1],  # hit frame
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

    def __init__(self, id, frame, candidate, database):
        # Linkback to candidates.
        self.candidates = {}

        # Store BLASTs by reciprocal databases (ie, species name).
        self.reciprocal_blasts = {}

        self.update_candidates(candidate)

        self.reciprocal = False
        self.rank = 0
        self.id = id
        if frame > 0:
            self.frame = '+%d' % frame
        else:
            self.frame = str(frame)
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
        self.set_sequence()

        self.write_fasta()

    def update_candidates(self, candidate):
        '''Update list of linked candidate genes.'''
        if not candidate.gene_name in self.candidates.keys():
            self.candidates[candidate.gene_name] = candidate

    def set_sequence(self):
        '''Get and set sequence from database.'''
        #XXX parse_seqids option is needed for parsing the sequence id.
        parsed = SeqIO.parse(self.database, 'fasta')
        for locus in parsed:
            if locus.id == self.id:
                self.sequence = str(locus.seq)
        try:
            dummy_seq = self.sequence
        except:
            print 'SEQUENCE NOT SET! Check the IDs.'
            sys.exit(2)

    def set_blast_output(self, organism):
        '''Build and return filepath for reverse blast output with organism name.'''
        if not organism in self.reciprocal_blasts.keys():
            o = organism.replace(' ', '_')
            blast_output = os.path.join(
                    self.reverse_results_folder, '%s-%s.xml' % (self.filename, o))
            self.reciprocal_blasts[organism] = {
                    'path': blast_output, 'sequences': []
                    }
        return self.reciprocal_blasts[organism]['path']

    def set_filepath(self):
        '''Build filepath from folder and filename.'''
        locus_filepath = os.path.join(self.reverse_folder, '%s.fa' % self.filename)
        self.filepath = locus_filepath

    def set_filename(self):
        '''Extract filename from ID.'''
        self.filename = self.parse_filename(self.id)

    def parse_filename(self, filename):
        '''Parse id to a valid filepath.'''
        # Remove spaces and slashes to avoid filepath issues.
        filename = filename.replace(' ', '_').replace('/', '-')
        return filename

    def write_fasta(self):
        '''Write FASTA file to filepath.'''
        wrote_fasta = SeqIO.write(SeqRecord(Seq(self.sequence), id=self.id, description=''), self.filepath, 'fasta')

    def parse_blast(self, organism):
        '''Parse reciprocal BLAST.'''
        print 'Parsing reciprocal BLAST: %s' % self.filename
        #XXX Does it parse everytime?
        reverse_blast_output = self.reciprocal_blasts[organism]['path']
        blast_records = NCBIXML.parse(open(reverse_blast_output))
        for blast_record in blast_records:
            for alignment in blast_record.alignments:
                for hsp in alignment.hsps:
                    if hsp.expect < self.EVALUE_THRESH:
                        try:
                            ref = alignment.title.split('|')[3]
                        except:
                            ref = alignment.title
                        sequence = Sequence(ref=ref)
                        sequence.evalue = hsp.expect
                        sequence.score = hsp.score
                        self.reciprocal_blasts[organism]['sequences'].append(sequence)

    def process(self):
        '''Define if locus is reciprocal to the gene.

        True: if it matches the same gene_id and ref from the candidate gene.
        '''
        #TODO Handle when gene_id is None?
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


def local_blast(blast_type, arguments):
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
        cline = NcbitblastxCommandline(**arguments)

    # Execute BLAST.
    stdout, stderr = cline()
    print '%s BLASTed!' % arguments['query']


def prepare(candidates, candidates_folder):
    '''Check candidate gene files (convert to FASTA, if needed).'''
    for gene in candidates:
        gene_filepath = os.path.join(candidates_folder, gene)
        gene_gb = gene_filepath.split('.')[0] + '.gb'
        gene_fa = gene_filepath.split('.')[0] + '.fa'

        if gene.endswith(('.gb', '.genbank')):
            continue
            # # Create record object.
            # record = SeqIO.read(gene_filepath, 'genbank')
            # # Create FASTA file.
            # try:
            #     fasta_file = open(gene_filepath.split('.')[0] + '.fa', 'r')
            # except IOError:
            #     fasta_file = open(gene_filepath.split('.')[0] + '.fa', 'w')
            #     fasta_file.write(record.format('fasta'))
            #     fasta_file.close()
        elif gene.endswith(('.fa', '.txt')):
            # Only fetch if file does not exist.
            try:
                open(gene_gb, 'r')
                pass
            except IOError:
                # 1. Parse ref
                record = SeqIO.read(gene_filepath, 'fasta')
                gene_ref = record.id.split('|')[3]

                # 2. Using ref to fetch genbank
                # Fetch entry and save handle to string.
                handle_string = efetcher(ref=gene_ref)

                # 3. Write genbank file
                if handle_string:
                    f = open(gene_gb, 'w')
                    f.write(handle_string)
                    f.close()
        else:
            logger.debug('File type not supported: %s', gene)


def efetcher(ref, db='protein', rettype='gp', retmode='txt'):
    '''Container function for Entrez.efetch.

    Tries to access Entrez, retry in case of failure.

    Return handle string ready to be written into file.
    '''
    logger.debug('Initiating Entrez via efetch...')
    logger.debug('ref=%s, db=%s, rettype=%s, retmode=%s', ref, db, rettype, retmode)

    # Trying to contact Entrez.
    try:
        handle = Entrez.efetch(db=db, id=ref, rettype=rettype, retmode=retmode)
        handle_string = handle.read()
        logger.debug('Fetched %s!', ref)
        return handle_string
    except:
        logger.debug('Entrez connection failed. Retrying...')
        handle = Entrez.efetch(db=db, id=ref, rettype=rettype, retmode=retmode)
        handle_string = handle.read()
        logger.debug('Fetched %s!', ref)
        return handle_string
    else:
        logger.critical('Entrez fetching failed for %s, giving up.', ref)
        return None


def usage():
    '''Explanation for arguments.'''

    print
    print 'USAGE: ./blaster.py -d Membranipora.fasta -b tblastn'
    print
    print '  -c, --candidates \n\tFolder with `.fa` or `.gb` files of candidate genes. One gene per file.'
    print '  -d, --database \n\tLocal database with new data (eg, transcriptome).'
    print '  -b, --blast \n\tBLAST command (blastn, blastp, blastx, tblastn, tblastx).'
    print '  -e, --email \n\tYour email is required for Entrez.'
    print
    #TODO Limit number of loci from candidate blast results.
    #TODO Specify threshold evalue.
    #TODO Limit the number of genes looked up during reciprocal blast.


def main(argv):

    # Available BLASTs.
    blast_types = ['blastn', 'blastp', 'blastx', 'tblastn', 'tblastx']
    # Database name.
    database = None
    # Default BLAST command.
    blast_type = 'tblastn'

    # Folder with candidate-genes.
    candidates_folder = 'candidates'
    # Folder with candidate-genes' BLASTs.
    candidates_results_folder = 'blasts'

    # Default email.
    Entrez.email = 'your@email.com'

    # Parse arguments.
    try:
        opts, args = getopt.getopt(argv, 'hc:d:b:e:', [
            'help',
            'candidates=',
            'database=',
            'blast=',
            'email=',
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
        elif opt in ('-e', '--email'):
            Entrez.email = arg

    # Print summary of arguments.
    logger.debug('Arguments: candidates=%s, database=%s, blast=%s, email=%s', candidates_folder, database, blast_type, Entrez.email)

    # Check if BLAST command was specified.
    if not blast_type:
        logger.critical('BLAST command was not specified (use "-b"). Aborting...')
        sys.exit(2)
    else:
        if not blast_type in blast_types:
            logger.critical('Unknown BLAST command: %s', blast_type)
            logger.info('Available BLAST commands: %s', ', '.join(blast_types))
            sys.exit(2)

    # Get candidate genes.
    #TODO Make it recursively search sub-directories.
    try:
        candidates = os.listdir(candidates_folder)
    except OSError:
        logger.critical('Folder "%s" does not exist! Aborting...', candidates_folder)
        sys.exit(2)

    # Issue error if there are no candidate genes.
    if not candidates:
        logger.critical('No candidate genes at %s folder! Aborting...',
                candidates_folder)
        sys.exit(2)

    # Check if results folder exists.
    if not os.path.isdir(candidates_results_folder):
        os.mkdir(candidates_results_folder)

    # Process candidate genes.
    prepare(candidates, candidates_folder)

    # Get proper genes, now.
    candidates = os.listdir(candidates_folder)

    # Only read GENBANK files.
    candidates = [file for file in candidates if file.endswith(('.gb', '.genbank'))]

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
                'outfmt': 5,  # Export in XML
                }
            # Execute BLAST command.
            local_blast(blast_type, arguments)

        logger.info('Parsing XML...')
        print '\nCandidate >> Gene: %s, ID:%s, Ref: %s, Description: %s' % (
                candidate.gene_name, candidate.gene_id,
                candidate.ref, candidate.description
                )
        candidate.parse_blast()

        #DO LOCUS STUFF
        for locus_id in candidate.loci.keys():
            # Instantiate Locus object.
            if locus_id in loci.keys():
                locus = loci[locus_id]
                # Update candidate list.
                locus.update_candidates(candidate)
            else:
                frame = candidate.loci[locus_id]['frame']
                locus = Locus(locus_id, frame, candidate, database)

            # Generate file for BLAST output.
            reverse_blast_output = locus.set_blast_output(candidate.organism)

            # Reciprocal BLAST querying locus sequence against human database.
            try:
                blastfile = open(reverse_blast_output)
                blastfile.close()
            except:
                # 1. Check available databases using set_reciprocal_db.
                candidate.set_reciprocal_db()

                if candidate.reciprocal_db:
                    # Sanitize filepaths for shell (ie handling "|" in locus IDs).
                    reverse_args = {
                            'query': '"%s"' % locus.filepath,
                            'db': candidate.reciprocal_db,
                            'out': '"%s"' % reverse_blast_output,
                            'outfmt': 5,
                            }
                    #XXX What to do when recirocal database is not protein?
                    if blast_type == 'tblastn':
                        local_blast('blastx', reverse_args)
                    elif blast_type == 'blastp':
                        local_blast('blastp', reverse_args)
                else:
                    # 2. If None, prepare for BLAST online based on the organism.
                    #XXX Handle blastp as above, also.
                    print 'BLASTing over NCBI... (may take a while).'
                    logger.info('BLASTing %s over NCBI (restricted to %s).',
                            locus.filepath, candidate.organism)
                    handle = NCBIWWW.qblast('blastx', 'refseq_protein', open(locus.filepath).read(), entrez_query='%s[Organism]' % candidate.organism)
                    reverse_file = open(reverse_blast_output, 'w')
                    reverse_file.write(handle.read())
                    reverse_file.close()
                    handle.close()

            # Parse reverse BLAST output.
            locus.parse_blast(candidate.organism)

            # Add locus to main list.
            loci[locus_id] = locus

        # Add to main dictionary.
        genes[candidate.gene_name] = candidate

    # Print loci equivalent to candidate genes.
    logger.info('Creating result files...')

    fasta_loci = []
    output_file = open('results.txt', 'w')

    #TODO Print date and variables used to results.txt.
    #TODO Translate sequence using the correct frame/strand and output it.
    for locus_id, locus in loci.iteritems():
        locus_description = '| frame: %s | candidates: %s' % (locus.frame, ', '.join(locus.candidates.keys()))
        fasta_loci.append(SeqRecord(Seq(locus.sequence), id=locus_id, description=locus_description))
        output_file.write('%s %s\n\n' % (locus_id, locus_description))

        output_file.write('\t\t\t\t{0:20}{1:20}{2:20}{3:1}\n'.format('gene', 'id', 'accession', 'e-value'))
        #output_file.write('\t\t\t\tgene\t\t\tid\t\taccession\t\te-value\n')

        for organism, reciprocal_blast in locus.reciprocal_blasts.iteritems():
            output_file.write('  {0:>2}\n'.format(organism))
            sequences = reciprocal_blast['sequences']

            # Instantiate first sequence to get gene id.
            try:
                first = sequences[0]
            except:
                output_file.write('\t\t\t\t{0}'.format('No reciprocal sequences found.\n'))
                continue

            current_gene_id = first.gene_id
            n_genes = 0
            for sequence in sequences:
                # Limits to 5 most similar genes.
                if sequence.gene_id == current_gene_id or n_genes < 5:
                    output_file.write('\t\t\t\t{0:20}{1:20}{2:20}{3:<.2e}'.format(sequence.gene_name, sequence.gene_id, sequence.ref, sequence.evalue))
                    accession_numbers = [seq.ref for seq in locus.candidates.values()]
                    if sequence.ref in accession_numbers:
                        output_file.write(' {0:>}'.format('<<'))
                    output_file.write('\n')
                else:
                    break
                current_gene_id = sequence.gene_id
                n_genes += 1
        output_file.write('\n\n')

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
