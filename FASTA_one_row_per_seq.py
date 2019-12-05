# -*- coding: utf-8 -*-

'''
-------------------------------------------------
    File Name:      FASTA_one_seq_per_row
    Description:    This script is designed to transfer these FASTA files with multiple rows for one sequence into the
                    format that one sequence only contains one row.
    Author:         Dengwei Zhang
    Date:           2019/11/4
    E-mail:         scdzzdw@163.com
-------------------------------------------------
'''

import re
import argparse

def read_file(file_in, file_out):
    if re.search('.gz$', file_in):
        fastq = gzip.open(file_in, 'rb')
        with fastq as fin, open(file_out, 'w') as fo:
            first_line = fin.readline().decode()
            print(first_line, file=fo)
            for line in fin:
                if re.search("^>", line.decode()):
                    content = '\n' + line.decode().rstrip('\n')
                    print(content, file=fo)
                else:
                    content = line.decode().rstrip('\n')
                    print(content, file=fo)
    else:
        fastq = open(file_in, 'r')
        with fastq as fin, open(file_out, 'w') as fo:
            first_line = fin.readline()
            print(first_line, file=fo)
            for line in fin:
                if re.search("^>", line):
                    content = '\n' + line.rstrip('\n')
                    print(content, file=fo)
                else:
                    content = line.rstrip('\n')
                    print(content, file=fo)

def main():
    parse = argparse.ArgumentParser(description="tranfer FASTA to one-seq-per-row")
    parse.add_argument("--file_1", "-i", help="input FASTA filename", required=True)
    parse.add_argument("--file_2", "-o", help="output FASTA filename", required=True)

    args = parse.parse_args()

    read_file(args.file_1, args.file_2)


if __name__ == "__main__":
    main()


