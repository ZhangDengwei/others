# -*- coding: utf-8 -*-
'''
#-------------------------------------------------
#    File Name:      convert_into_methylKit
#    Description:    the coverage file much omit strand information
#    Author:         Dengwei Zhang
#    Date:           2019/11/6
#    E-mail:         scdzzdw@163.com
#-------------------------------------------------
'''

import argparse

def convert(file_in, file_out):
    with open(file_in, 'r') as fi, open(file_out, 'w') as fo:
        print('chrBase\tchr\tbase\tstrand\tcoverage\tfreqC\tfreqT', file=fo)
        for line in fi:
            items = line.rstrip('\n').split()
            chrBase = items[0] + '.' + items[1]
            chr = items[0]
            base = items[1]
            strand = 'R'
            coverage = str(int(items[4]) + int(items[5]))
            freqC = str(round((int(items[4]) / int(coverage))*100, 2))
            freqT = str(round((int(items[5]) / int(coverage))*100, 2))
            print(chrBase, '\t', chr, '\t', base, '\t', strand, '\t', coverage, '\t', freqC, '\t', freqT, file=fo)

def main():
    parse = argparse.ArgumentParser(description="convert cov file from Bismark into the file can be loaded into methylKit")
    parse.add_argument("--file_1", "-i", help="input filename", required=True)
    parse.add_argument("--file_2", "-o", help="output filename", required=True)

    args = parse.parse_args()

    convert(args.file_1, args.file_2)


if __name__ == "__main__":
    main()

