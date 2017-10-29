from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
import pandas as pd
import argparse
import os

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--gbk', type=argparse.FileType('r'), required=True, nargs='+')
    parser.add_argument('--fna_out')
    parser.add_argument('--genomes_out')
    parser.add_argument('--faa_out')
    parser.add_argument('--tbl_out', required=True)
    parser.add_argument('--sample_id')
    args = parser.parse_args()
    dna_seqs = []
    genome_seqs = []
    amino_seqs = []
    contigs = []
    feature_ids = []
    keywords = []
    for f in args.gbk:
        for record in SeqIO.parse(f, 'genbank'):
            record.description = ''
            genome_seqs.append(record)
            for feature in record.features:
                if feature.type == 'CDS':
                    dna_record, aa_record, feature_id, keyword = extract_gene(record, feature, args.sample_id)
                    dna_seqs.append(dna_record)
                    amino_seqs.append(aa_record)
                    contig_id = record.id
                    if args.sample_id:
                        contig_id = '%s__%s' % (args.sample_id, contig_id)
                    contigs.append(contig_id)
                    feature_ids.append(feature_id)
                    keywords.append(keyword)

    if args.fna_out:
        SeqIO.write(dna_seqs, args.fna_out, 'fasta')
    if args.faa_out:
        SeqIO.write(amino_seqs, args.faa_out, 'fasta')
    if args.genomes_out:
        SeqIO.write(genome_seqs, args.genomes_out, 'fasta')

    df = pd.DataFrame({'contig': contigs, 'id': feature_ids, 'keywords': keywords})
    df = df[['id', 'contig', 'keywords']]
    df.to_csv(args.tbl_out, index=False)


def extract_gene(record, feature, sample_id):

    start = feature.location.start
    end = feature.location.end
    strand = feature.location.strand
    dna_seq = record.seq[start:end]

    feature_id = feature.qualifiers['locus_tag'][0]
    feature_id = feature_id.replace(' ', '')

    if sample_id:
        feature_id = '%s__%s' % (sample_id, feature_id)

    if strand == -1:
        dna_seq = dna_seq.reverse_complement()

    dna_record = SeqRecord(dna_seq, id=feature_id, description='')

    assert len(feature.qualifiers['translation']) == 1
    aa_seq = feature.qualifiers['translation'][0]

    keyword = feature.qualifiers['product'][0]
    if len(keyword) ==0:
        keyword = 'None'

    aa_record = SeqRecord(Seq(aa_seq, IUPAC.protein), id=feature_id, description='')
    return dna_record, aa_record, feature_id, keyword

if __name__ == '__main__':
    main()