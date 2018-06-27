import re
import pandas as pd
from itertools import islice
import argparse

parser = argparse.ArgumentParser(description='Process file')
parser.add_argument('--coverage', type=float)
parser.add_argument('--seed', type=int)
parser.add_argument('--in_file', type=str)
parser.add_argument('--out_file', type=str)
args = parser.parse_args()

coverage = args.coverage
seed = args.seed


in_file = args.in_file

contig_rgx = re.compile('(\S+):(\d+)-(\d+) log:')
tot_reads_rgx = re.compile('Total Reads: (\d+), Coverage: (\d+), minDepth: (\d+)')
confirmed_rgx = re.compile('Confirmed (\d+) of (\d+) bases')
corrected_rgx = re.compile('Corrected (\d+) snps; (\d+) ambiguous bases; corrected (\d+) small insertions totaling (\d+) bases, (\d+) small deletions totaling (\d+) bases')
bases_rgx = re.compile('In total, there were (\d+) and (\d+) bases in forward and reverse libraries, respectively')
chunk = []
rows = []
tot_bases = 0

def parse_chunk(chunk):
	if len(chunk)> 0:
		contig_name = chunk[0]
		start_locus = chunk[1]
		end_locus = chunk[2]
		tot_reads = chunk[4]
		confirmed = chunk[5]
		corrected = chunk[6]
		
		m = tot_reads_rgx.search(chunk[4])
		n = confirmed_rgx.search(chunk[5])
		o = corrected_rgx.search(chunk[6])
		rows.append([contig_name, tot_bases, seed, start_locus, end_locus, 
			m.group(1), m.group(2), m.group(3),
		 	n.group(1), n.group(2), 
		 	o.group(1), o.group(2), o.group(3), o.group(4), o.group(5), o.group(6)])
		
	

with open(in_file, 'r') as handle:
	skip = True
	for line in handle.readlines():
		bases = bases_rgx.search(line.strip())
		if bases:
			tot_bases = (int(bases.group(1)) + int(bases.group(2)))/1e9

		if not skip:
			chunk.append(line.strip())
		m = contig_rgx.search(line.strip())
		if m:
			parse_chunk(chunk)
			contig_name = m.group(1)
			start_locus = int(m.group(2))
			end_locus = int(m.group(3))
			chunk = [contig_name, start_locus, end_locus]
			skip = False
			
		if line.startswith('#'):
			skip = True
		
			
parse_chunk(chunk)
df = pd.DataFrame(rows,  columns=['contig', 'tot_gigabases', 'seed', 'start', 'end', 'total_reads', 'coverage', 'minDepth', 'confirmed_bases', 'tot_bases', 'snps', 'ambig_bases', 'fixed_insertions',
	'fixed_insertions_bases', 'fixed_deletions', 'fixed_deletions_bases'])
df.to_csv(args.out_file, index=False)
		
