from Bio.Blast.Applications import NcbiblastnCommandline
import subprocess
import argparse
from io import StringIO
from Bio.Blast import NCBIXML
import os
import glob
from Bio.Seq import Seq
import multiprocessing


class RMLSTGeneTree(object):
    def __init__(self, arguments):
        self.rmlst_gene_path = arguments.database
        self.fasta_directory = arguments.fasta_directory
        self.threads = arguments.threads
        self.outfile = arguments.outfile
        fasta_files = glob.glob(self.fasta_directory + '/*.f*a')
        self.check_for_databases()
        pool = multiprocessing.Pool(processes=self.threads)
        sequences = pool.map(RMLSTGeneTree.extract_rmlst_genes, fasta_files)
        pool.close()
        pool.join()
        f = open(self.outfile + '.fasta', 'w')
        for seq in sequences:
            f.write(seq + '\n')
        f.close()
        # Align sequences with MAFFT, and create tree with FASTTREE.
        print('Aligning sequences...')
        with open(arguments.outfile + '.log', 'w') as logfile:
            cmd = 'mafft --auto {} > {}'.format(self.outfile + '.fasta', self.outfile + '_aligned.fasta')
            subprocess.call(cmd, shell=True, stderr=logfile)
            print('Making tree...')
            cmd = 'fasttree -nt -gtr {} > {}'.format(self.outfile + '_aligned.fasta', self.outfile + '.nwk')
            subprocess.call(cmd, shell=True, stderr=logfile)
            print('DONE! Tree can be found in {} :D'.format(arguments.outfile + '.nwk'))

    def check_for_databases(self):
        """
        Checks that the BLAST databases that we want to use have been created.
        If they haven't, they'll get created.
        """
        db_files = ['.nhr', '.nin', '.nsq']
        for i in range(1, 10):
            database = self.rmlst_gene_path + 'BACT00000' + str(i) + '.tfa'
            db_present = True
            for db_file in db_files:
                if os.path.exists(database):
                    if not os.path.exists(database + db_file):
                        db_present = False
                if not db_present:
                    cmd = 'makeblastdb -dbtype nucl -in ' + database
                    print('Creating DB for ' + database)
                    subprocess.call(cmd, shell=True)
        for i in range(10, 66):
            database = self.rmlst_gene_path + 'BACT0000' + str(i) + '.tfa'
            db_present = True
            for db_file in db_files:
                if os.path.exists(db_file):
                    if not os.path.exists(database + db_file):
                        db_present = False
                if not db_present:
                    cmd = 'makeblastdb -dbtype nucl -in ' + database
                    print('Creating DB for ' + database)
                    subprocess.call(cmd, shell=True)

    @staticmethod
    def extract_rmlst_genes(assembly):
        """
        Extracts the top hit for each rMLST gene from an assembly.
        :param assembly: Path for fasta file you want rMLST genes from. Must be uncompressed.
        :return: total_seq: String containing concatenation of all rMLST genes for the assembly.
        """
        print('Extracting rMLST genes from {}...'.format(assembly))
        rmlst_gene_path = '/mnt/nas/Adam/assemblypipeline/rMLST/2017-03-29/'
        total_seq = '>' + assembly + '\n'
        for i in range(1, 10):
            fasta_file = rmlst_gene_path + 'BACT00000' + str(i) + '.tfa'
            if os.path.exists(fasta_file):
                seq = RMLSTGeneTree.get_top_blast_hit(assembly, fasta_file)
                if seq != 'NA':
                    total_seq += seq
        for i in range(10, 66):
            fasta_file = rmlst_gene_path + 'BACT0000' + str(i) + '.tfa'
            if os.path.exists(fasta_file):
                seq = RMLSTGeneTree.get_top_blast_hit(assembly, fasta_file)
                if seq != 'NA':
                    total_seq += seq
        return total_seq

    @staticmethod
    def get_top_blast_hit(query, database):
        """
        :param query: Query sequence. Expects an assembly in FASTA format.
        :param database: Database, in fasta format.
        :return: sequence: Top hit from query to the database, reverse complemented if necessary.
        """
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
    parser.add_argument('outfile', help='Base name for your output files.')
    parser.add_argument('database', help='Path to rMLST database.')
    parser.add_argument('-t', '--threads', type=int, default=cpu_count, help='Number of threads to run analysis on.'
                                                                             'Default is number of cores on machine.')
    arguments = parser.parse_args()
    # Figure out what fasta files to do analysis on.
    RMLSTGeneTree(arguments)

