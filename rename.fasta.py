import argparse
import sys
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import pandas as pd

def main():
	parser = argparse.ArgumentParser()
	
	parser.add_argument('--infile', nargs='?', type=argparse.FileType('r'), default=sys.stdin)
	parser.add_argument('--outfile', nargs='?', type=argparse.FileType('w'), default=sys.stdout)
	parser.add_argument('--mapping_file', type=argparse.FileType('w'))
	parser.add_argument('--prefix', type=str, required=True)
	parser.add_argument('--no_desc', dest='no_desc', default=True, action='store_true')
	parser.set_defaults(no_desc=False)
	args = parser.parse_args()
		
	out_seqs = []
	records = []
	counter = 1
	for record in SeqIO.parse(args.infile, 'fasta'):
		new_id = "%s_%05d" % (args.prefix, counter)
		if (args.mapping_file):
			records.append((record.id, new_id))

		record.id = new_id
		if (args.no_desc):
			record.description = ''
		out_seqs.append(record)
		counter +=1
		if (len(out_seqs)>1000):
			SeqIO.write(out_seqs, args.outfile, 'fasta')
			out_seqs = []
	
	SeqIO.write(out_seqs, args.outfile, 'fasta')

	if(args.mapping_file):
		df = pd.DataFrame(records, columns=['old_id', 'new_id'])
		df.to_csv(args.mapping_file, index=False)
	
if __name__ == '__main__':
    main()	
