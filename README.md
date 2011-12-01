# BLASTing batch script

This script takes a series of candidate genes, BLAST them against a local 
database (eg, transcriptome from an organism) to extract relevant sequences, 
BLAST these sequences against the reciprocal database from candidate genes (eg, human protein), retrieve and print the results for analysis.

## Arguments

`-c, --candidates` -- Folder with `.fasta` or `.gb` files of candidate genes. 
One gene per file.

`-d, --database` -- Local database with new data (eg, transcriptome). If in 
`FASTA` format it will be converted to a database `BLAST+` can understand using
`makeblastdb` command, example: `makeblastdb -in bugula.fasta -parse_seqids -dbtype nucl`

`-b, -blast` -- Folder where local BLAST results will be written.

## Usage

1. Convert databases to be BLASTable (if needed).

    makeblastdb -in membranipora.fa -parse_seqids -dbtype nucl
    makeblastdb -in human_protein.fa -parse_seqids -dbtype prot

2. Copy FASTA files with 1 gene sequence per file into the `candidates` folder.

3. Run the BLASTer command (database could also be in another folder).

    python blaster.py -d Membranipora.fasta -b tblastn

4. Results (loci) are printed to 2 files. `results.txt` shows each locus with reciprocal BLAST output and e-values. `results.fa` aggregates the loci and their sequences.

## Internals

1. Each candidate is BLASTed against the local database using the specified BLAST program.

2. Best 3 loci hits are copied to `reverse` folder as FASTA files.

3. Each locus is BLASTed against the reciprocal database to check if it returns the same gene.

4. Sequences from the best 5 genes are printed in the main output file.

## Default values

EVALUE_THRESH = 0.001
