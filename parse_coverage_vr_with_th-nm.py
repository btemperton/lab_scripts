import os
import glob
import os
import argparse
import numpy as np
import pandas
import pysam


def get_seq_len_from_bam(samfile):
    temp = []
    for i, dic in enumerate(samfile.header['SQ']):
        try:
            # if args.v:
            # print "index: {0} key : {1} value: {2}".format(i,dic['SN'],dic['LN'])
            temp.append({'Seq': dic['SN'], 'Length': dic['LN']})
        except Exception:
            print("i: {0} d: {1}".format(i, dic))
            raise
    return temp


def coverage_vectors(contigs_size):
    coverage = {}
    for i in contigs_size:
        temp = {}
        temp["positions"] = np.zeros(i["Length"])
        temp["nb_reads"] = 0
        temp["nb_bp"] = 0
        coverage[i["Seq"]] = temp
    return coverage


def parse(samfile, coverage, th):
    for l in samfile.fetch():
        try:
            length = l.reference_length
            # print "length: {0}".format(length)
            if length > 0:
                nm_index = l.get_tag("NM")
                # print "NM index: {0}".format(nm_index)
                mismatch = int(nm_index)
                # m=re.match("(\d+)$",nm_index)
                # if m:
                # print "\t{0} has the number {1}".format(nm_index,m.group(1))
                # mismatch+=int(m.group(1))
                # else:
                # print "!!! I can't deal with this {0}".format(nm_index)
                pcent_id = float(length - mismatch) / float(length) * 100
                # print "pcent_id: {0}".format(pcent_id)
                # print "final summary: tag was {0}, which is convered into {1} mismatches for {2} length, so {3} % identity".format(l.get_tag("NM"),mismatch,length,pcent_id)
                if pcent_id >= th:
                    coverage[samfile.getrname(l.tid)]["nb_reads"] += 1
                    if args.pv == "0.7.5":
                        # pysam 0.7.5 (for BamM compatibility on OSC)
                        begin = l.pos + 1
                        end = l.pos + l.alen
                        coverage[samfile.getrname(l.tid)]["positions"][begin:end] = 1
                        coverage[samfile.getrname(l.tid)]["nb_bp"] += l.alen
                    else:
                        # pysam 0.8.4
                        begin = l.reference_start + 1
                        end = l.reference_start + l.reference_length
                        # print "begin {0} to end {1} since the length is {2}".format(begin, end, l.reference_length);
                        coverage[samfile.getrname(l.reference_id)]["positions"][begin:end] = 1
                        coverage[samfile.getrname(l.tid)]["nb_bp"] += l.reference_length
        except Exception:
            print("line: {0}".format(l))
            raise
    return coverage


def process(inputdir, code, outputfile, pv, th):
    list_files = glob.glob(inputdir + '/' + code + '_*_sorted.bam')
    print("Listing all the relevant files: {0} ").format(list_files)
    if pv == "0.7.5":
        samfile = pysam.Samfile(list_files[0], "rb")  # for older pysam version on OSC
    else:
        samfile = pysam.AlignmentFile(list_files[0], "rb")
    print("Reading the length of sequences from bam file header")
    contigs_size = get_seq_len_from_bam(samfile)
    # for index in contigs_size:
    # print "{0} - {1}".format(index["Seq"],index["Length"])
    print("Getting the coverage tables ready")
    coverage = coverage_vectors(contigs_size)
    print("Getting coverage values with th {0}".format(th))
    for bamfile in list_files:
        print("Reading {0}".format(bamfile))
        if pv == "0.7.5":
            samfile = pysam.Samfile(bamfile, "rb")  # for older pysam version on OSC
        else:
            samfile = pysam.AlignmentFile(bamfile, "rb")
        coverage = parse(samfile, coverage, th)
    print("Calculating sequence coverage ratios and writing the output")
    coverage_prop = {}
    for contig, vector in coverage.items():
        temp = {}
        for i in contigs_size:
            if contig == i["Seq"]:
                temp["length"] = i["Length"]
        temp["length_covered"] = np.sum(vector["positions"]) / float(len(vector["positions"])) * 100
        temp["nb_reads"] = vector["nb_reads"]
        temp["nb_bp"] = vector["nb_bp"]
        if vector["nb_reads"] > 1:
            coverage_prop[contig] = temp
    print("Now writing the output")
    output = pandas.DataFrame(coverage_prop).transpose()
    if (len(coverage_prop) > 0):
        output = output.sort(['nb_bp', 'length_covered'], ascending=[0, 0])
    # Only sort if there is something, otherwise no nb_bp so not happy
    # output.columns=["Length coverage (%)","Nb bp","Nb reads"]
    output.to_csv(outputfile)
    samfile.close()


def file_is_empty(path):
    return os.stat(path).st_size == 0


def main():
    parser = argparse.ArgumentParser(
        description='Parse a directory of bam files all mapped to the same reference database.')
    parser.add_argument('--bam_dir', '-b', dest='bam_dir', required=True,
                        help='directory of input bam files (sorted and indexed)')
    parser.add_argument('--code', '-c', dest='code', required=True,
                        help='code of the sample')
    parser.add_argument('--out', '-o', dest='outdir', required=True,
                        help='output directory')
    parser.add_argument('--th', '-t', dest='th', type=int, choices=range(0, 101), required=False, default=0,
                        help='minimum coverage of reference sequence (0-100), default no threshold', metavar="[0-100]")
    parser.add_argument("--verbose", '-v', dest='v', help="verbose mode",
                        action="store_true")
    parser.add_argument("--pysam_version", '-pv', dest='pv',
                        help="used to specify the version of pysam available (specific variable names for version 0.7.5)",
                        default="0.8.4")
    global args
    args = parser.parse_args()
    # print "toto"
    outfi = args.outdir + "/" + args.code + "_coverage.csv";
    if os.path.isfile(outfi):
        print("{0} is already there, we don't recompute".format(outfi))
    else:
        process(args.bam_dir, args.code, outfi, args.pv, args.th)


if __name__ == "__main__":
    output = main()
