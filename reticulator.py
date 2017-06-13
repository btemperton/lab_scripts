import argparse
import numpy as np
import pandas as pd
import scipy.stats as stats
import re
import os
import itertools
import logging
import subprocess
from multiprocessing import Pool  # Process pool
from multiprocessing import sharedctypes

logger = logging.getLogger('The Reticulator')
logger.setLevel(logging.INFO)


def parse_mcl_dump_file(dump_file):
	"""
	Reads in a dump file from MCL and converts it into a dictionary of
	protein clusters
	:param dump_file: The file to be read in
	:return: a panda dataframe
	"""
	logger.info('Parsing the MCL dump file.....')
	rows = []
	counter = 1
	id_rgx = re.compile('ref\|(.*)\|')
	for line in dump_file.readlines():
		bits = line.strip().split()
		for i in bits:
			m = id_rgx.search(i)  ## This removes any irritating ref| bits
			if m:
				i = m.group(1)
			rows.append(('PC_%06d' % counter, i))
		counter += 1
	df = pd.DataFrame.from_records(rows, columns=['PCid', 'node'])

	df.to_csv('%s/PC.clusters.txt' % OUTPUT, sep='\t', index=False)
	return df


def load_gene_map(map_file):
	"""
	Reads in a file mapping genes to contigs
	:param map_file: The mapping file with gene ids and contigs comma-separated
	:return: a pandas dataframe
	"""
	logger.info('Loading gene map file.....')
	return pd.read_csv(map_file)


def add_pc_labels(pc_dict, gm_df):
	"""
	Merge the PC dictionary with the gene map data frame
	:param pc_dict: The dictionary of PCs created with inverted mcl parse
	:param gm_df: The dataframe gene map
	:return: a pandas dataframe
	"""
	logger.info('Appending Protein Cluster labels.....')
	df = pd.merge(gm_df, pc_dict, how='left', on='node')
	df.to_csv('%s/PC.dictionary.txt' % OUTPUT, sep='\t', index=False)
	return df


def create_composition_mtx(pc_df, presence_absence=False):
	"""
	Create a composition matrix where M_ij if contig i encodes at least 1 PC j
	:param pc_df: a pandas dataframe of contigs, nodes and PC.ids
	:param presence_absence: If True, rescales values >1 to 1
	:return: an i x j matrix (rows are contigs, columns are PCs
	"""
	logger.info('Creating composition matrix.....')
	v = pc_df.groupby(['contig', 'PCid']).size().reset_index()
	v.columns = ['contig', 'PCid', 'count']
	rtnValue = v.pivot(index='contig', columns='PCid', values='count').fillna(0)
	if presence_absence:
		rtnValue[rtnValue > 1] = 1

	rtnValue.to_csv('%s/composition.mtx.txt' % OUTPUT, sep='\t')

	return rtnValue


def calculate_shared_content(args):
	"""
	Calculates the shared content based on the n x n blocks of BLOCK_SIZE
	:param args: the block windows
	:return: Nothing - it updates the SHARED_ARRAY block
	"""
	window_x, window_y = args
	tmp = np.ctypeslib.as_array(SHARED_ARRAY)
	for idx_x in range(window_x, window_x + BLOCK_SIZE):
		for idx_y in range(window_y, window_y + BLOCK_SIZE):

			if idx_x == idx_y:
				continue
			if idx_x >= len(PC_DF.index) or idx_y >= len(PC_DF.index):
				continue
			try:
				shared = PC_DF.ix[idx_x] + PC_DF.ix[idx_y]
				tmp[idx_x, idx_y] = shared.value_counts()[2.0]
			except KeyError:
				pass
			except IndexError:
				print "idx_x: %i, idx_y: %i" % (idx_x, idx_y)


def calculate_shared_matrix():
	"""
	Calculates the number of shared PCs between two contigs
	:return: A DataFrame containing the shared content between each contig
	"""
	window_idxs = [(i, j) for i, j in itertools.product(range(0, len(PC_DF.index), BLOCK_SIZE),
														range(0, len(PC_DF.index), BLOCK_SIZE))]

	p = Pool(processes=THREADS)
	p.map(calculate_shared_content, window_idxs)

	arr = np.ctypeslib.as_array(SHARED_ARRAY)
	rtnValue = pd.DataFrame(arr, index=PC_DF.index, columns=PC_DF.index)
	rtnValue.to_csv('%s/shared.mtx.txt' % OUTPUT, sep='\t')
	return rtnValue


def calculate_hypergeometric_survival(pc_df, shared_mtx):
	"""
	Creates a symmetrical matrix containing the hypergeometric survival function (1-cdf)
	for each contig pair
	:param pc_df: The dataframe of PCs - to work out the number of PCs per contig
	:param shared_mtx: The matrix with the number of shared PCs between contigs
	:return: a symmetrical matrix of probabilities
	"""
	logger.info('Calculating hypergeometric survival.....')
	pc_counts = pc_df.groupby(['contig']).PCid.nunique().reset_index()
	pc_counts.columns = ['contig', 'count']
	pc_counts.index = pc_counts['contig']
	contigs = list(shared_mtx.index)

	rtnValue = pd.DataFrame(0, index=contigs, columns=contigs)

	for i in itertools.combinations(contigs, 2):
		number_common_pcs = shared_mtx.loc[i[0], i[1]]
		a_pc_count = pc_counts.at[i[0], 'count']
		b_pc_count = pc_counts.at[i[1], 'count']
		a, b = sorted((a_pc_count, b_pc_count))
		total_pcs = pc_df.PCid.nunique()
		T = 0.5 * total_pcs * (total_pcs - 1)
		logT = np.log10(T)

		# The -1 is needed here to calculate the inverse cumlutative
		# density function
		pval = stats.hypergeom.sf(number_common_pcs - 1, total_pcs, a, b)
		sig = min(300, np.nan_to_num(-np.log10(pval) - logT))

		# If sig > 1, it is considered a match
		if sig <= 1:
			sig = 0
		rtnValue.at[i[0], i[1]] = sig

	# this matrix needs to be symmetrical for accurate calculations of weighting so
	rtnValue = rtnValue + rtnValue.T

	return rtnValue


def cluster_viruses_by_mcl(cluster_file, inflation=2.0):
	"""
	Uses the hypergeometric survival matrix to cluster the viruses using MCL
	:param cluster_file: The hypergeometric survival matrix file in long format
	:param inflation: A parameter to control inflation for MCL.
	:return: the name of the output file
	"""
	logger.info('Clustering viruses with MCL at inflation %.1f......' % inflation)
	cmd = "mcxload -abc %s --stream-mirror -o %s/viral.clusters.mci -write-tab %s/viral.clusters.tab" % (
		cluster_file, OUTPUT, OUTPUT)
	logger.info(cmd)
	fail = subprocess.check_output(cmd.split())

	if not fail:
		cmd = "mcl %s/viral.clusters.mci -I %.1f -use-tab %s/viral.clusters.tab -o %s/viral.clusters.mcl.I%i" % (
			OUTPUT, inflation, OUTPUT, OUTPUT, inflation * 10)
		logger.info(cmd)
		subprocess.check_output(cmd.split())

	return "%s/viral.clusters.mcl.I%i" % (OUTPUT, inflation * 10)


def calculate_membership_values(classification_df, hypergeometric_survival_df):
	logger.info('Calculating group membership........')
	membership_df = classification_df.pivot(index='node', columns='VCid', values='conservative_membership').fillna(0)
	membership_df[membership_df > 0] = 0

	connectivity_df = classification_df.pivot(index='node', columns='VCid', values='conservative_membership').fillna(0)
	connectivity_df[membership_df > 0] = 0

	edge_sums = hypergeometric_survival_df.sum(axis=1, numeric_only=True).astype(float)
	clusters = classification_df['node'].groupby(classification_df['VCid'])

	# for each cluster, work out the sum of edges of that cluster to each contig
	# then work out the proportion of those compared to all a contigs edges

	# if a cluster only has a single member, then the membership of that cluster must be 1
	for group_name, group in clusters:
		if len(group.values) == 1:
			membership_df.loc[group.values[0], group_name] = 1

	for row_id, data in hypergeometric_survival_df.iterrows():
		for group_name, group in clusters:

			cluster_sum = np.sum(data[group.values])
			cluster_mean = np.mean(data[group.values])
			logger.debug("summing edges between %s and %s" % (row_id, group.values))
			logger.debug("which would be %s,%s" % (row_id, str(data[group.values].values)))
			logger.debug("The sum of which is: %i " % cluster_sum)
			if edge_sums[row_id] > 0:
				membership_df.loc[row_id, group_name] = cluster_sum / edge_sums[row_id]
			connectivity_df.loc[row_id, group_name] = cluster_mean

	membership_df.fillna(0).to_csv('%s/VC.membership.txt' % OUTPUT, sep='\t')
	connectivity_df.fillna(0).to_csv('%s/VC.connectivity.txt' % OUTPUT, sep='\t')
	return membership_df, connectivity_df


def output_classification(mcl_cluster_file, cluster_name_prefix):
	logging.info('Outputting conservative membership.....')

	rows = []
	counter = 1
	with open(mcl_cluster_file, 'r') as handle:
		for line in handle.readlines():
			bits = line.strip().split()
			for i in bits:
				rows.append(("%s_%05d" % (cluster_name_prefix, counter), i, 1))
			counter += 1
	df = pd.DataFrame.from_records(rows, columns=['VCid', 'node', 'conservative_membership'])
	df.to_csv('%s/VC.clusters.txt' % OUTPUT, sep='\t', index=False)
	return df


logging.basicConfig()
parser = argparse.ArgumentParser()
parser.add_argument('--mcl_dump', dest='mcl_dump', type=argparse.FileType('r'), required=True)
parser.add_argument('--gene_map', dest='gene_map', metavar='STRING', required=True)
parser.add_argument('--output', dest='output', metavar='STRING', required=True)
parser.add_argument('--cluster_prefix', dest='cluster_prefix', metavar='STRING', default='VC')
parser.add_argument('--threads', dest='threads', type=int, default=1)
args = parser.parse_args()

if not os.path.exists(args.output):
	os.makedirs(args.output)

OUTPUT = args.output
THREADS = args.threads
BLOCK_SIZE=4

s = parse_mcl_dump_file(args.mcl_dump)
t = load_gene_map(args.gene_map)
u = add_pc_labels(s, t)
PC_DF = create_composition_mtx(u, presence_absence=True)
result = np.ctypeslib.as_ctypes(np.zeros((PC_DF.shape[0], PC_DF.shape[0])))
SHARED_ARRAY = sharedctypes.RawArray(result._type_, result)
w = calculate_shared_matrix()
hypergeometric_df = calculate_hypergeometric_survival(u, w)
hypergeometric_df.to_csv('%s/hypergeometric.survival.txt' % OUTPUT, sep='\t')

df = hypergeometric_df.where(np.triu(np.ones(hypergeometric_df.shape)).astype(np.bool))
df = df.stack().reset_index()
df.columns = ['Row', 'Column', 'Value']
df.to_csv('%s/hypergeometric.survival.long.txt' % OUTPUT, sep='\t', index=False, header=False)
viral_cluster_file = cluster_viruses_by_mcl('%s/hypergeometric.survival.long.txt' % OUTPUT)
classification_df = output_classification(viral_cluster_file, args.cluster_prefix)
calculate_membership_values(classification_df, hypergeometric_df)