import pandas as pd
import numpy as np
from Bio import SeqIO

with open('/Users/bt273/Downloads/short.hybrid/2k5/short.hybrid.scaffolds.2k5.fa', 'r') as handle:
    seq_df = pd.DataFrame([(record.id, len(record.seq), record.id[0]) for record in SeqIO.parse(handle, 'fasta')],
                          columns=['contig_id', 'contig_len', 'contig_type'])


virfinder_df = pd.read_csv('/Users/bt273/Downloads/short.hybrid/2k5/virfinder.combined.csv',
                           usecols=[0,4], names=['contig_id', 'virfinder_qvalue'], header=0)
virfinder_df['contig_id'] = virfinder_df['contig_id'].str.split().str.get(0)


merged_df = pd.merge(seq_df, virfinder_df, on='contig_id', how='left')


nucmer_coords_df = pd.read_csv('/Users/bt273/Downloads/short.hybrid/2k5/clstr.c80.i95/nucmer.csv')

contig_records = []


def parse_otu(lines):
    representative = lines[0][1]
    cluster_name = lines[0][0].replace('>', '')
    is_rep = 'Y'
    rep_sim = 'NA'

    contig_records.append((representative,
                           cluster_name,
                           is_rep,
                           rep_sim, 'NA', 'NA'))
    for line in lines[1:]:
        contig_name = line[0]
        contig_sim = float(line[1])

        slice = nucmer_coords_df[(nucmer_coords_df['ref_name'] == representative) &
                               (nucmer_coords_df['query_name'] == contig_name)]
        max_aln_row = slice.ix[slice['ref_len_aln'].idxmax()]
        ref_start, ref_finish = max_aln_row[['ref_start', 'ref_finish']]
        contig_records.append((contig_name,
                               cluster_name,
                               'N',
                               contig_sim, ref_start, ref_finish))





with open('/Users/bt273/Downloads/short.hybrid/2k5/clstr.c80.i95/reads_95-80.clstr', 'r') as handle:
    lines = []
    cluster_count = 0
    for line in handle.readlines():
        if line.startswith('>') and len(lines) > 0:
            parse_otu(lines)
            if cluster_count % 10 == 0:
                print('completed %i clusters' % cluster_count)
            cluster_count += 1
            lines = []
            lines.append(line.strip().split())
        else:
            lines.append(line.strip().split())
    parse_otu(lines)

cluster_df = pd.DataFrame(contig_records, columns=['contig_id', 'cluster_name', 'is_cluster_representative', 'similarity_to_cluster_representative', 'ref_start', 'ref_finish'])

merged_df = pd.merge(merged_df, cluster_df, on='contig_id', how='left')

##coverage

short_coverage_df = pd.read_csv('/Users/bt273/Downloads/short.hybrid/2k5/short.coverage.tbl',
                                sep='\t', header=0, usecols=[0,2,3,4],
                                names=['contig_id', 'GC', 'mean_coverage', 'percent_covered'])
pat = r'(.*)'
repl = lambda m: "S_%s" % m.group(1)
short_coverage_df['contig_id'] = short_coverage_df['contig_id'].str.replace(pat, repl)


hybrid_coverage_df = pd.read_csv('/Users/bt273/Downloads/short.hybrid/2k5/hybrid.coverage.tbl',
                                sep='\t', header=0, usecols=[0,2,3,4],
                                names=['contig_id', 'GC', 'mean_coverage', 'percent_covered'])
pat = r'(.*)'
repl = lambda m: "H_%s" % m.group(1)
hybrid_coverage_df['contig_id'] = hybrid_coverage_df['contig_id'].str.replace(pat, repl)

coverage_df = pd.concat([short_coverage_df, hybrid_coverage_df], ignore_index=True)

merged_df = pd.merge(merged_df, coverage_df, on='contig_id', how='left')

print(merged_df.head())
merged_df.sort_values(['cluster_name', 'is_cluster_representative'], ascending=[True, False], inplace=True)

merged_df.to_csv('/Users/bt273/Downloads/short.hybrid/2k5/short.hybrid.2k5.analysis.csv', index=False)



