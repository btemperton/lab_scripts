from Bio import SeqIO
import pandas as pd
import os
import argparse


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--clstr_file', type=argparse.FileType('r'), required=True)
    parser.add_argument('--fna_file', type=argparse.FileType('r'), required=True)
    parser.add_argument('--dir_out', required=True)
    parser.add_argument('--pc_out', required=True)
    args = parser.parse_args()
    seqs = SeqIO.index(args.fna_file, 'fasta')

    if not os.path.exists(args.dir_out):
        os.makedirs(args.dir_out)

    pc_counter = 1

    membership = []

    for line in args.clstr_file.readlines():
        pc_id = 'PC_%06d' % pc_counter
        ids = line.strip().split()
        to_keep = [seqs[x] for x in ids]
        membership += [(pc_id, x) for x in ids]
        pc_counter += 1
        if len(to_keep) > 1:
            SeqIO.write(to_keep, '%s/%s.fa' % (args.dir_out, pc_id), 'fasta')

    df = pd.DataFrame(membership, columns=['pc_id', 'gene_id'])
    df.to_csv(args.pc_out, index=False)


if __name__ == '__main__':
    main()
