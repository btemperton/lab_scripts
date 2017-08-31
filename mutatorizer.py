import argparse
import logging
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
import random as rnd
import sys

logger = logging.getLogger('Mutatorizer')
logger.setLevel(logging.INFO)

bases = set(list('ATCG'))

def main():
    logging.basicConfig()
    parser = argparse.ArgumentParser()
    parser.add_argument('--fasta', dest='fasta_file', metavar='STRING', required=True, type=str)
    parser.add_argument('--num_fragments', dest='num_fragments', metavar='int', required=True, type=int)
    parser.add_argument('--mean_frag_size', dest='frag_size_mu', metavar='int', required=True, type=int)
    parser.add_argument('--frag_size_std', dest='frag_size_sigma', metavar='int', required=True, type=int)
    parser.add_argument('--mean_mutation_rate', dest='mutation_rate_mu', metavar='float', required=True, type=float)
    parser.add_argument('--mutation_rate_std', dest='mutation_rate_sigma', metavar='float', required=True, type=float)
    parser.add_argument('--output', dest='output_file', metavar='string', required=True, type=str)
    args = parser.parse_args()

    outhandle = open(args.output_file, 'w')

    generated_seqs = []
    for record in SeqIO.parse(args.fasta_file, 'fasta'):
        base_id = record.id
        base_seq = str(record.seq)
        while len(generated_seqs) < args.num_fragments:
            try:
                mutation_rate = rnd.gauss(args.mutation_rate_mu, args.mutation_rate_sigma)
                subsequence = subselect_sequence(base_seq, args.frag_size_mu, args.frag_size_sigma)
                mutated_subsequence = mutate(subsequence, mutation_rate)
                new_id = '%s__mut_%.2f__len_%i' % (base_id, mutation_rate, len(subsequence))
                generated_seqs.append(SeqRecord(Seq(mutated_subsequence, IUPAC.IUPACAmbiguousDNA()),
                                            id=new_id, name=new_id, description=''))
            except Exception as e:
                print(e)

    SeqIO.write(generated_seqs, outhandle, 'fasta')

    outhandle.close()

def subselect_sequence(sequence, mean_size, sigma):
    length_of_sub = int(rnd.gauss(mean_size, sigma))
    start = rnd.randint(0, len(sequence) -length_of_sub)
    selection = list(sequence)[start:start+length_of_sub]
    return ''.join(selection)





def mutate(sequence, mutation_rate):
    list_seq = list(sequence)
    for i,b in enumerate(list_seq):
        val = rnd.uniform(0,1)
        if val <= mutation_rate:
            #convert to a set
            set_b = set(list(b))
            #get the options of mutation
            mutation_options = list((bases - set_b))
            list_seq[i] = rnd.choice(mutation_options)
    return ''.join(list_seq)



if __name__ == '__main__':
    main()
