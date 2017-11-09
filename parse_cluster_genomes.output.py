import pandas as pd
import numpy as np
in_file = '/Users/bt273/Downloads/reads_95-80.clstr'

contig_records = []
cluster_records = []

def parse_otu(lines):
    representative = lines[0][1]
    rep_len = int(lines[0][2])
    cluster_name = lines[0][0].replace('>', '')
    rep_type = representative[0]
    is_rep = 'Y'
    rep_sim = 'N/A'
    cluster_size = len(lines)

    cluster_type = 'unclassified'

    if cluster_size == 1:
        cluster_type = 'singleton - %s' % rep_type

    contig_records.append((cluster_name,
                    cluster_size,
                    representative,
                    rep_len,
                    rep_type,
                    is_rep,
                    rep_sim))

    for line in lines[1:]:
        contig_name = line[0]
        contig_sim = float(line[1])
        contig_len = int(line[2])
        contig_type = contig_name[0]

        if contig_type == 'S' and rep_type == 'H':
            cluster_type = 'hybrid extension'
        elif contig_type == 'H' and rep_type == 'S':
            cluster_type = 'short longer than hybrid'

        contig_records.append((cluster_name,
                               cluster_size,
                               contig_name,
                               contig_len,
                               contig_type,
                               'N',
                               contig_sim))

    #we only want to classify stuff as identical if the types are different
    if cluster_size ==2 and rep_type != lines[1][0][0]:
        if float(lines[1][2])/rep_len > .99:
            cluster_type = '99  % identical'
        if rep_len == int(lines[1][2]):
            cluster_type = 'identical'


    cluster_records.append((cluster_name, cluster_size, cluster_type))




with open(in_file, 'r') as handle:
    lines = []
    for line in handle.readlines():
        if line.startswith('>') and len(lines) > 0:
            parse_otu(lines)
            lines = []
            lines.append(line.strip().split())
        else:
            lines.append(line.strip().split())


cluster_df = pd.DataFrame(cluster_records, columns = ['cluster_name', 'cluster_size', 'cluster_type'])
cluster_df.to_csv('/Users/bt273/Downloads/clusters.csv', index=False)

contig_df = pd.DataFrame(contig_records, columns = ['cluster_name', 'cluster_size', 'contig_name', 'contig_len', 'contig_type', 'is_rep', 'contig_similarity'])
contig_df.to_csv('/Users/bt273/Downloads/contigs.csv', index=False)

print('A total of %i contigs were organised into %i clusters, which broke down as follows:' % (contig_df.shape[0], cluster_df.shape[0]))
print(cluster_df.groupby('cluster_type').size().reset_index())

print('**** Cluster 24 is weird & needs to be investigated further')
