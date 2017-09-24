import nanomath
import pysam
import pandas as pd
import re
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC

names = {'Lambda_NEB':('lambda', 0),
         'NC_021790.1':('phi_18:1', 39189),
         'NC_021796.1':('phi_38:1',72534),
         'KC821629.1':('phi_38:2',54012),
         'KF302034.1':('HM1',129439),
         'KF302037.1':('HP1',45035),
         'KF302036.1':('HS2', 38208)}


seq_dict = {'NC_021790.1':[],
         'NC_021796.1':[],
         'KC821629.1':[],
         'KF302034.1':[],
         'KF302037.1':[],
         'KF302036.1':[]}

pID_cutoff = 80

def parseMD(MDlist):
    return sum([len(item) for item in re.split('[0-9^]', MDlist)])


def parseCIGAR(cigartuples):
    return sum([item[1] for item in cigartuples if item[0] == 1])


def parseContig(bam, contig_name):
    query_names = []
    reference_names = []
    reference_lengths = []
    lengths = []
    alignedLengths = []
    quals = []
    alignedQuals = []
    mapQ = []
    pIDs = []

    for read in samfile.fetch(reference=contig_name):
        if not read.is_secondary:
            query_names.append(read.query_name)
            reference_names.append(names[read.reference_name][0])
            reference_lengths.append(names[read.reference_name][1])
            quals.append(nanomath.aveQual(read.query_qualities))
            alignedQuals.append(nanomath.aveQual(read.query_alignment_qualities))
            lengths.append(read.query_length)
            alignedLengths.append(read.query_alignment_length)
            mapQ.append(read.mapping_quality)
            try:
                pID = (1 - read.get_tag("NM") / read.query_alignment_length) * 100
            except KeyError:
                pID = (1 -
                            (parseMD(read.get_tag("MD")) + parseCIGAR(read.cigartuples))
                            / read.query_alignment_length) * 100
            pIDs.append(pID)

            if(pID >= pID_cutoff):
                new_read = SeqRecord(Seq(read.query_sequence, IUPAC.ambiguous_dna), id=read.query_name, name='')
                seq_dict[contig_name].append(new_read)


    datadf = pd.DataFrame(data={
        'query_name' : query_names,
        'reference_name' : reference_names,
        'reference_length': reference_lengths,
        'query_length' : lengths,
        'aligned_length' : alignedLengths,
        'average_read_quality' : quals,
        'average_alignment_quality' : alignedQuals,
        'mapping_quality' : mapQ,
        'percent_identity' : pIDs })

    datadf = datadf[['query_name',
                    'reference_name',
                    'reference_length',
                    'query_length',
                    'aligned_length',
                    'average_read_quality',
                    'average_alignment_quality',
                    'mapping_quality',
                    'percent_identity']]

    return (datadf)

bam='/Users/bt273/Downloads/carrier.seq.mockB.nanofilt.l1k.h50.t50.q6.align.genomes.bam'
samfile = pysam.AlignmentFile(bam, "rb")

df = pd.DataFrame()

ids = [x for x in samfile.references if x != 'Lambda_NEB']
for i in ids:
    print('Sorting out %s' % i)
    df = df.append(parseContig(samfile, i))

df.to_csv('/Users/bt273/Downloads/carrier.seq.mockB.nanofilt.l1k.h50.t50.q6.align.genomes.tbl.gz',
          index=False, compression='gzip')

for k,v in seq_dict.items():
    outhandle = open('/Users/bt273/Downloads/%s.reads.80pc.fa' % k, 'w')
    SeqIO.write(v, outhandle, 'fasta')
    outhandle.close()

