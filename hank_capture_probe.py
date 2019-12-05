# -*- coding: utf-8 -*-
'''
#-------------------------------------------------
#    File Name:      hank_capture_probe
#    Description:
#    Author:         Dengwei Zhang
#    Date:           2019/11/6
#    E-mail:         scdzzdw@163.com
#-------------------------------------------------
'''

'''
The input fasta file must be the following format:
name    proble_sequence start_position_in_read end_postion
'''

import re
import argparse
import gzip


def fq(file):
    if re.search('.gz$', file):
        fastq = gzip.open(file, 'rb')
    else:
        fastq = open(file, 'r')
    with fastq as f:
        while True:
            l1 = f.readline()
            if not l1:
                break
            l2 = f.readline()
            l3 = f.readline()
            l4 = f.readline()
            yield [l1, l2, l3, l4]

def extract(fasta, fastq, report):
    probe = dict()
    fastq_content = fq(fastq)
    with open(fasta, 'r') as fa:
        for line in fa:
            items = line.rstrip('\n').split()
            name = items[0]
            seq = items[1]
            start = items[2]
            end = items[3]
            probe[name] = [seq, start, end]

    with open(report, 'w') as fo:
        for i in iter(fastq_content):
            read_name = i[0].decode().rstrip('\n')
            seq = i[1].decode().rstrip('\n')
            for i, (j,k,l) in probe.items():
                pattern = re.compile(j)
                search = pattern.search(seq)
                if search and search.start() >= int(k) and search.end() <= int(l):
                    print(read_name, '\t', seq, '\t', name, file=fo)
                else:
                    pass


def main():
    parse = argparse.ArgumentParser(description="probe caputer for RenHan")
    parse.add_argument("--fasta", "-a", help="input sample file with FASTA format", required=True)
    parse.add_argument("--fastq", "-q", help="input fastq file", required=True)
    parse.add_argument("--output", "-o", help="output file", required=True)

    args = parse.parse_args()

    extract(args.fasta, args.fastq, args.output)


if __name__ == "__main__":
    main()

