import argparse
import sys
import numpy as np
import pandas as pd
import scipy.stats as stats
import re
from Bio import SeqIO

def main():
	parser = argparse.ArgumentParser()
	parser.add_argument('--faa_file', dest='faa_file', type=argparse.FileType('r'), required=True)
	parser.add_argument('--type', dest='type', metavar='STRING', required=True)
	parser.add_argument('--outfile', dest='outfile', type=argparse.FileType('w'), required=True)
	args = parser.parse_args()
	
	args.outfile.write('node,contig\n')
	
	desc_rgx = re.compile('.*\[(.*)\]')
        id_rgx = re.compile('ref\|(.*)\|')
        for record in SeqIO.parse(args.faa_file, 'fasta'):
                m = desc_rgx.search(record.description)
                contig = 'ERROR'
                if m:
                        contig = m.group(1)
                
                idstr = record.id
                m = id_rgx.search(record.id)
                if m:
                        idstr = m.group(1)           
                                                 
                args.outfile.write('%s,%s\n' % (idstr, contig))
        args.outfile.close()

if __name__ == '__main__':
    	main()