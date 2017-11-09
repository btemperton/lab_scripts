import pandas as pd
from Bio import SeqIO
import multiprocessing as mp
import logging
import tqdm
import argparse


parser = argparse.ArgumentParser()
parser.add_argument('--scaffolds', metavar='STRING', required=True)
parser.add_argument('--virfinder', metavar='STRING', required=True)
parser.add_argument('--nucmer_csv', metavar='STRING', required=True)
parser.add_argument('--clstr', metavar='STRING', required=True)
parser.add_argument('--short_cov', metavar='STRING', required=True)
parser.add_argument('--hybrid_cov', metavar='STRING', required=True)
parser.add_argument('--output', metavar='STRING', required=True)
args = parser.parse_args()

logging.basicConfig(level=logging.DEBUG,
                    format='%(asctime)s %(name)s.%(funcName)s +%(lineno)s: %(levelname)-8s [%(process)d] %(message)s',
                    )

with open(args.scaffolds, 'r') as handle:
    seq_df = pd.DataFrame([(record.id, len(record.seq), record.id[0]) for record in SeqIO.parse(handle, 'fasta')],
                          columns=['contig_id', 'contig_len', 'contig_type'])


virfinder_df = pd.read_csv(args.virfinder,
                           usecols=[0,4], names=['contig_id', 'virfinder_qvalue'], header=0)
virfinder_df['contig_id'] = virfinder_df['contig_id'].str.split().str.get(0)


merged_df = pd.merge(seq_df, virfinder_df, on='contig_id', how='left')


nucmer_coords_df = pd.read_csv(args.nucmer_csv)


def parse_otu(lines):

    representative = lines[0][1]
    cluster_name = lines[0][0].replace('>', '')
    is_rep = 'Y'
    rep_sim = 'NA'

    records = []

    records.append((representative,
                           cluster_name,
                           is_rep,
                           rep_sim, 'NA', 'NA'))
    for line in lines[1:]:
        contig_name = line[0]
        contig_sim = float(line[1])

        slice = nucmer_coords_df[(nucmer_coords_df['ref_name'] == representative) &
                               (nucmer_coords_df['query_name'] == contig_name) & (nucmer_coords_df['ref_len_aln'] >=1000)]

        ref_start = '|'.join(list(slice['ref_start']))
        ref_finish = '|'.join(list(slice['ref_finish']))
        records.append((contig_name,
                               cluster_name,
                               'N',
                               contig_sim, ref_start, ref_finish))
    #logging.info('Completed processing of %s' % cluster_name)
    return records


clusters = []


with open(args.clstr, 'r') as handle:
    lines = []
    cluster_count = 0
    for line in handle.readlines():
        if line.startswith('>') and len(lines) > 0:
            clusters.append(lines)
            lines = []
            lines.append(line.strip().split())
        else:
            lines.append(line.strip().split())
    clusters.append(lines)


with mp.Pool(4) as p:
    total_records = list(tqdm.tqdm(p.imap(parse_otu, clusters), total=len(clusters)))

contig_records = [record for sublist in total_records for record in sublist]

cluster_df = pd.DataFrame(contig_records, columns=['contig_id', 'cluster_name', 'is_cluster_representative', 'similarity_to_cluster_representative', 'ref_start', 'ref_finish'])

merged_df = pd.merge(merged_df, cluster_df, on='contig_id', how='left')

##coverage

short_coverage_df = pd.read_csv(args.short_cov,
                                sep='\t', header=0, usecols=[0,2,3,4],
                                names=['contig_id', 'GC', 'mean_coverage', 'percent_covered'])
pat = r'(.*)'
repl = lambda m: "S_%s" % m.group(1)
short_coverage_df['contig_id'] = short_coverage_df['contig_id'].str.replace(pat, repl)


hybrid_coverage_df = pd.read_csv(args.hybrid_cov,
                                sep='\t', header=0, usecols=[0,2,3,4],
                                names=['contig_id', 'GC', 'mean_coverage', 'percent_covered'])
pat = r'(.*)'
repl = lambda m: "H_%s" % m.group(1)
hybrid_coverage_df['contig_id'] = hybrid_coverage_df['contig_id'].str.replace(pat, repl)

coverage_df = pd.concat([short_coverage_df, hybrid_coverage_df], ignore_index=True)

merged_df = pd.merge(merged_df, coverage_df, on='contig_id', how='left')

print(merged_df.head())
merged_df.sort_values(['cluster_name', 'is_cluster_representative'], ascending=[True, False], inplace=True)

merged_df.to_csv(args.output, index=False)



