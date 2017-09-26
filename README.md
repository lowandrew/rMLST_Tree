# rMLST Tree

This script will extract rMLST genes from FASTA files and then concatenate them and
make a tree (or at least, it will eventually. It's a work in progress!)

### Program Requirements
- MAFFT (>=7.271) installed and present on your $PATH
- FastTree (>= 2.1.8) installed and on your $PATH
- Python 3.5 (2.7 should also work, though that isn't tested.)
- NCBI blast+ (>=2.2.31) 

### Python Package Requirements
- biopython >= 1.70

### Usage
`python rMLST_tree.py fasta_directory outfile_name rMLST_database`
- fasta_directory is a path to a directory containing the assemblies you want to make a tree out of.
- outfile_name is the name of your output files. If outfile_name is example, you'll end up with example.fasta containing
concatenated genes, example_aligned.fasta containing the aligned concatenated genes, and example.nwk containing your
tree in newick format.
- rMLST_database is the path to a directory containing all fastas of all 53 rMLST genes. If you're connected to OLC NAS,
use /mnt/nas/Adam/assemblypipeline/rMLST/2017-03-29

#### Options (More to come.)
- threads (-t): Number of threads to run analysis on. Default is number of cores on your system.

