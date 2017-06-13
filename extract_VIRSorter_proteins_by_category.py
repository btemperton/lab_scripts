import argparse
from Bio import SeqIO
import pandas as pd
import logging

logger = logging.getLogger('extract VIRSorter')
logger.setLevel(logging.INFO)


def main():
    logging.basicConfig()
    parser = argparse.ArgumentParser()
    parser.add_argument('--VIRSorter_prots_file', dest='virsorter_prots', type=argparse.FileType('r'), required=True)
    parser.add_argument('--VIRSorter_global_phage_signal_file', dest='virsorter_signal', metavar='STRING', required=True)
    parser.add_argument('--faa_output', dest='faa_output', type=argparse.FileType('w'), required=True)
    parser.add_argument('--protein_csv', dest='protein_csv', type=argparse.FileType('w'), required=True)
    parser.add_argument('--sample_name', dest='sample_name', metavar='STRING', default='Virsorter')
    parser.add_argument('--categories', dest='categories', metavar='STRING', default='1,2')
    args = parser.parse_args()

    virsorter_output_csv = pd.read_csv(args.virsorter_signal, comment='#', header=None, index_col=False,
                                       names = 'contig_id,nb_genes_contig,fragment,nb_genes,category,nb.hallmark,enrichment,non-caudovirales_enrichment,pfam_depletion,unchar_enrichment,strand_switch,short_genes'.split(','))

    df = virsorter_output_csv.loc[virsorter_output_csv['category'].isin([int(x) for x in args.categories.split(',')])]

    contigs = df.contig_id.values

    output = []

    args.protein_csv.write('node,contig,keywords\n')

    for record in SeqIO.parse(args.virsorter_prots, 'fasta'):
        contig_id = record.id.split('-gene')[0]
        if contig_id in contigs:
            output.append(record)
            args.protein_csv.write('%s,%s,%s\n' % (record.id, contig_id, args.sample_name))

    SeqIO.write(output, args.faa_output, 'fasta')

if __name__ == '__main__':
    main()
