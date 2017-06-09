import argparse
import sys
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


def main():
	parser = argparse.ArgumentParser()
	parser.add_argument('--infile', nargs='?', type=argparse.FileType('r'), default=sys.stdin)
	parser.add_argument('--outfile', nargs='?', type=argparse.FileType('w'), default=sys.stdout)
	parser.add_argument('--window', type=int)
	parser.add_argument('--step', type=int)
	args = parser.parse_args()


	chunks = createFragments(args.infile, args.window, args.step)
	SeqIO.write(chunks, args.outfile, 'fasta')		

def createFragments(in_reads, window_size=20, step_size=1):
    chunks = []
    for record in SeqIO.parse(in_reads, 'fasta'):
    	try:
		x = list(slidingWindow(record, window_size, step_size))
		for i in x:
			chunks.append(i)
	except Exception:
		continue
    return chunks
    
def slidingWindow(record, winSize, step=1):
    """Returns a generator that will iterate through
    the defined chunks of input sequence.  Input sequence
    must be iterable."""

    sequence = record.seq
    id = record.id
    # Verify the inputs
    try:
        it = iter(sequence)
    except TypeError:
        raise Exception("**ERROR** sequence must be iterable.")
    if not ((type(winSize) == type(0)) and (type(step) == type(0))):
        raise Exception("**ERROR** type(winSize) and type(step) must be int.")
    if step > winSize:
        raise Exception("**ERROR** step must not be larger than winSize.")
    if winSize > len(sequence):
        raise Exception("**ERROR** winSize must not be larger than sequence length.")

    # Pre-compute number of chunks to emit
    numOfChunks = ((len(sequence) - winSize) / step) + 1

    # Do the work
    for i in range(0, numOfChunks * step, step):
        fragment = sequence[i:i + winSize]
        frag_id = '%s__%i-%i' % (id, i, i + winSize)
        frag_record = SeqRecord(fragment, id=frag_id, name=record.name, description='')
        yield frag_record
        
if __name__ == '__main__':
    main()
