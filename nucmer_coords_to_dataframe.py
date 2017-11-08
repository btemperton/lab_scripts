import argparse
import logging
import pandas as pd


def main():
    logging.basicConfig()
    parser = argparse.ArgumentParser()
    parser.add_argument('--coords', metavar='STRING', required=True)
    parser.add_argument('--output', metavar='STRING', required=True)
    args = parser.parse_args()

    records = []
    with open(args.coords, 'r') as handle:
        for _ in range(5):
            next(handle)
        for line in handle.readlines():
            bits = line.strip().split('|')
            ref_start, ref_finish = bits[0].split()
            query_start, query_finish = bits[1].split()
            ref_len_aln, query_len_aln = bits[2].split()
            identity = bits[3]
            ref_name, query_name = bits[4].split()
            record = [x.strip() for x in [ref_name, query_name, ref_start, ref_finish, ref_len_aln, identity, query_start, query_finish]]
            if record[0] != record[1]:
                records.append(record)

    df  = pd.DataFrame(records, columns=['ref_name', 'query_name', 'ref_start', 'ref_finish', 'ref_len_aln', 'identity', 'query_start', 'query_finish'])
    df.to_csv(args.output, index=False)

if __name__ == '__main__':
    main()