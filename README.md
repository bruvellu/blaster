# BLASTing batch script

This script takes a series of candidate genes, BLAST them against a local 
database (eg, transcriptome from an organism) to extract related loci, 
BLAST these sequences against the reciprocal database from candidate genes (eg, human protein), retrieve and print the results for analysis.

## Installation

1. Install dependencies.

    - Python 2.7 (not Python 3)
    - Biopython 1.56
    - BLAST+ 2.2.25

2. Copy the `blaster.py` to a new folder where it will run.

3. Make sure you know where the local and reciprocal databases are (or copy 
   them to the folder).

## Usage

1. Convert databases to be BLASTable (if needed).

        makeblastdb -in membranipora.fa -parse_seqids -dbtype nucl
        makeblastdb -in human_protein.fa -parse_seqids -dbtype prot

    NOTE: `parse_seqids` argument is required for parsing the locus id from the 
    database using the Biopython SeqIO module.

2. Copy FASTA files with 1 gene sequence per file into any folder (eg, `candidates`).

3. Run the BLASTer command (database could also be in another folder).

        python blaster.py -c folder -d membranipora.fa -b tblastn -e your@email.com

4. Results (loci) are printed to 2 files. `results.txt` shows each locus with reciprocal BLAST output and e-values. `results.fa` aggregates the loci and their sequences.

## Arguments

`-c, --candidates` -- Folder with `.fa`, `.gb`, or `.txt` files of candidate genes. 
One gene per file.

`-d, --database` -- Local database with new data (eg, transcriptome).

`-b, --blast` -- BLAST command (`blastn`, `blastp`, `blastx`, `tblastn`, `tblastx`).

`-e, --email` -- Your email is required for Entrez.

## Internals

1. Each candidate is BLASTed against the local database using the specified BLAST program (default is `tblastn`).

2. Best 3 loci hits are copied to `reverse` folder as FASTA files.

3. Each locus is BLASTed against the reciprocal database to check if it returns the same gene.

4. Sequences from the best 5 genes are printed in the main output file.

## Default values

EVALUE_THRESH = 0.001

LOCI_HITS = 3

RECIPROCAL_GENES = 5
