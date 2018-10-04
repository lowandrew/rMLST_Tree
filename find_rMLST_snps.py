#!/usr/bin/env python

import os
import argparse
from io import StringIO
from Bio.Seq import Seq
from Bio import pairwise2
from Bio.Blast import NCBIXML
from Bio.Blast.Applications import NcbiblastnCommandline


def extract_rmlst_genes(assembly_one, assembly_two, rmlst_folder, genes_to_ignore):
    """
    Extracts the top hit for each rMLST gene from an assembly.
    """
    total_num_snps = 0
    for i in range(1, 10):
        fasta_file = os.path.join(rmlst_folder, 'BACT00000' + str(i) + '.tfa')
        if os.path.exists(fasta_file) and os.path.split(fasta_file)[1] not in genes_to_ignore:
            seq_one = get_top_blast_hit(assembly_one, fasta_file)
            seq_two = get_top_blast_hit(assembly_two, fasta_file)
            num_identities = pairwise2.align.globalxx(seq_one, seq_two, score_only=True)
            sequence_length = min(len(seq_one), len(seq_two))
            num_snps = sequence_length - num_identities
            total_num_snps += num_snps
    for i in range(10, 66):
        fasta_file = os.path.join(rmlst_folder, 'BACT0000' + str(i) + '.tfa')
        if os.path.exists(fasta_file) and os.path.split(fasta_file)[1] not in genes_to_ignore:
            seq_one = get_top_blast_hit(assembly_one, fasta_file)
            seq_two = get_top_blast_hit(assembly_two, fasta_file)
            num_identities = pairwise2.align.globalxx(seq_one, seq_two, score_only=True)
            sequence_length = min(len(seq_one), len(seq_two))
            num_snps = sequence_length - num_identities
            total_num_snps += num_snps
    return total_num_snps


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
    parser = argparse.ArgumentParser()
    parser.add_argument('-1', '--assembly_one',
                        type=str,
                        required=True,
                        help='Full path to first fasta formatted (uncompressed) file you want to compare.')
    parser.add_argument('-2', '--assembly_two',
                        type=str,
                        required=True,
                        help='Full path to second fasta formatted (uncompressed) file you want to compare.')
    parser.add_argument('-r', '--rmlst_folder',
                        type=str,
                        required=True,
                        help='Path to rMLST gene folder - should contain BACT000001 through BACT0000065.tfa - '
                             'it\'s assumed that each .tfa file has already had a BLAST nucleotide database '
                             'created for it.')
    parser.add_argument('-i', '--ignore',
                        nargs='+',
                        default=[],
                        help='Names of any genes to ignore, since some genes appear twice in some species. Must be '
                             'the full name of the FASTA file - aka BACT000065.tfa')
    args = parser.parse_args()
    num_snps = extract_rmlst_genes(assembly_one=args.assembly_one,
                                   assembly_two=args.assembly_two,
                                   rmlst_folder=args.rmlst_folder,
                                   genes_to_ignore=args.ignore)
    print('{},{},{}'.format(os.path.split(args.assembly_one)[1], os.path.split(args.assembly_two)[1], int(num_snps)))
