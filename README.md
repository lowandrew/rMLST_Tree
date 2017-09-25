# rMLST Tree

This script will extract rMLST genes from FASTA files and then concatenate them and
make a tree (or at least, it will eventually. It's a work in progress!)

### Program Requirements
- Muscle (>=3.8.31) installed and present on your $PATH
- FastTree (>= 2.1.8) installed and on your $PATH
- Python 3.5 (2.7 should also work, though that isn't tested.)
- NCBI blast+ (>=2.2.31) 

### Python Package Requirements
- biopython >= 1.70

### Usage
- Program takes a folder with fasta files in it as input. Will end up writing a tree to tree.nwk (this will change in future versions.)
- Runs through new_detector.py

#### Options (More to come.)
- threads (-t): Number of threads to run analysis on. Default is number of cores on your system.

