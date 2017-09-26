from Bio.Blast.Applications import NcbiblastnCommandline
import argparse
from io import StringIO
from Bio.Blast import NCBIXML
import os
import glob
from Bio.Seq import Seq
import multiprocessing


def extract_rmlst_genes(assembly):
    print('Extracting rMLST genes from {}...'.format(assembly))
    rmlst_gene_path = '/mnt/nas/Adam/assemblypipeline/rMLST/2017-03-29/'
    total_seq = '>' + assembly + '\n'
    for i in range(1, 10):
        fasta_file = rmlst_gene_path + 'BACT00000' + str(i) + '.tfa'
        if os.path.exists(fasta_file):
            seq = get_top_blast_hit(assembly, fasta_file)
            if seq != 'NA':
                total_seq += seq
    for i in range(10, 66):
        fasta_file = rmlst_gene_path + 'BACT0000' + str(i) + '.tfa'
        if os.path.exists(fasta_file):
            seq = get_top_blast_hit(assembly, fasta_file)
            if seq != 'NA':
                total_seq += seq
    return total_seq


def get_top_blast_hit(query, database):
    db_files = ['.nhr', '.nin', '.nsq']
    db_present = True
    for db_file in db_files:
        if not os.path.exists(database + db_file):
            db_present = False
    if not db_present:
        cmd = 'makeblastdb -dbtype nucl -in ' + database
        os.system(cmd)
    blastn = NcbiblastnCommandline(db=database, outfmt=5, query=query)
    stdout, stderr = blastn()
    count = 0
    j = 0
    sequence = 'NA'
    for record in NCBIXML.parse(StringIO(stdout)):
        if count > 0:
            break
        for alignment in record.alignments:
            for hsp in alignment.hsps:
                if j > 0:
                    break
                if hsp.align_length - hsp.sbjct_start < 1:
                    sequence = str(Seq(hsp.query).reverse_complement())
                else:
                    sequence = str(hsp.query)
                # print(sequence)
                j += 1
                count += 1

    return sequence

if __name__ == '__main__':
    cpu_count = multiprocessing.cpu_count()
    parser = argparse.ArgumentParser()
    parser.add_argument('fasta_directory', help='Folder containing fasta files you want to create an rMLST tree with.')
    parser.add_argument('-t', '--threads', type=int, default=cpu_count, help='Number of threads to run analysis on.'
                                                                             'Default is number of cores on machine.')
    arguments = parser.parse_args()
    fasta_files = glob.glob(arguments.fasta_directory + '/*.f*a')
    pool = multiprocessing.Pool(processes=12)
    sequences = pool.map(extract_rmlst_genes, fasta_files)
    pool.close()
    pool.join()
    f = open('sequences.txt', 'w')
    for seq in sequences:
        f.write(seq + '\n')
    f.close()
    print('Aligning sequences...')
    cmd = 'muscle -in sequences.txt -out sequences_aligned.txt'
    os.system(cmd)
    print('Making tree...')
    cmd = 'fasttree -nt -gtr sequences_aligned.txt > tree.nwk'
    os.system(cmd)
    print('DONE! :D')
