import argparse
import sys
from Bio import SeqIO
import pandas as pd
import string
import random


def id_generator(size=6):
	suffix = ''.join(random.SystemRandom().choice(string.ascii_uppercase +
												  string.ascii_lowercase +
												  string.digits) for _ in range(size - 1))
	return random.SystemRandom().choice(string.ascii_uppercase + string.ascii_lowercase) + suffix


def main():
	parser = argparse.ArgumentParser()

	parser.add_argument('--infile', nargs='?', type=argparse.FileType('r'), default=sys.stdin)
	parser.add_argument('--outfile', nargs='?', type=argparse.FileType('w'), default=sys.stdout)
	parser.add_argument('--mapping_file', type=argparse.FileType('w'))
	parser.add_argument('--prefix', type=str)
	parser.add_argument('--random', default=False, action='store_true')
	parser.add_argument('--sample', type=str, default='NA')
	parser.add_argument('--no_desc', dest='no_desc', default=True, action='store_true')
	parser.set_defaults(no_desc=False)
	args = parser.parse_args()

	out_seqs = []
	records = []
	counter = 1
	for record in SeqIO.parse(args.infile, 'fasta'):
		if args.random:
			new_id = id_generator(8)
		else:
			new_id = "%s_%06d" % (args.prefix, counter)
		if args.mapping_file:
			records.append((new_id, record.id, args.sample))

		record.id = new_id
		if args.no_desc:
			record.description = ''
		out_seqs.append(record)
		counter += 1

	SeqIO.write(out_seqs, args.outfile, 'fasta')

	if args.mapping_file:
		df = pd.DataFrame(records, columns=['old_id', 'new_id', 'sample'])
		df.to_csv(args.mapping_file, index=False, sep='\t', header=False)


if __name__ == '__main__':
	main()
