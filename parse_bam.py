import nanomath
import pysam
import pandas as pd
import re

names = {'Lambda_NEB':'lambda',
         'NC_021790.1':'phi_18:1',
         'NC_021796.1':'phi_38:1',
         'KC821629.1':'phi_38:2',
         'KF302034.1':'HM1',
         'KF302037.1':'HP1',
         'KF302036.1':'HS2'}

def parseMD(MDlist):
    return sum([len(item) for item in re.split('[0-9^]', MDlist)])


def parseCIGAR(cigartuples):
    return sum([item[1] for item in cigartuples if item[0] == 1])


def parseContig(bam, contig_name):
    query_names = []
    reference_names = []
    lengths = []
    alignedLengths = []
    quals = []
    alignedQuals = []
    mapQ = []
    pID = []

    for read in samfile.fetch(reference=contig_name):
        if not read.is_secondary:
            query_names.append(read.query_name)
            reference_names.append(names[read.reference_name])
            quals.append(nanomath.aveQual(read.query_qualities))
            alignedQuals.append(nanomath.aveQual(read.query_alignment_qualities))
            lengths.append(read.query_length)
            alignedLengths.append(read.query_alignment_length)
            mapQ.append(read.mapping_quality)
            try:
                pID.append((1 - read.get_tag("NM") / read.query_alignment_length) * 100)
            except KeyError:
                pID.append((1 -
                            (parseMD(read.get_tag("MD")) + parseCIGAR(read.cigartuples))
                            / read.query_alignment_length) * 100)


    datadf = pd.DataFrame(data={
        'query_name' : query_names,
        'reference_name' : reference_names,
        'query_length' : lengths,
        'aligned_length' : alignedLengths,
        'average_read_quality' : quals,
        'average_alignment_quality' : alignedQuals,
        'mapping_quality' : mapQ,
        'percent_identity' : pID })

    datadf = datadf[['query_name',
                    'reference_name',
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

