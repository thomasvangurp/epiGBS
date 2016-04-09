#!/usr/bin/env python
from itertools import izip
import argparse

def parse_args():
    "Pass command line arguments"
    parser = argparse.ArgumentParser(description='Combine watson and crick consensus and recreate original seq')
    parser.add_argument('-w','--watson',
                        help='watson consensus')
    parser.add_argument('-c','--crick',
                    help='crick consensus')
    parser.add_argument('-o','--consensus',
                        help='consensus out')
    parser.add_argument('-b','--bam',
                    help='bamfile in')
    args = parser.parse_args()
    return args


def make_ref2(args):
    """"Make a concensus from watson and crick"""
    watson = open(args.watson,'r')
    crick = open(args.crick,'r')
    out = open(args.consensus, 'w')
    count = 0
    try:
        while True:
            w = [watson.readline()[1:-1],watson.readline()[:-1].upper()]
            c = [crick.readline()[1:-1],crick.readline()[:-1].upper()]
            if w[0] == '':
                break
            while True:
                if w[0] == c[0]:
                    break
                else:
                    if int(w[0]) > int(c[0]):
                        c = [crick.readline()[1:-1],crick.readline()[:-1]]
                    else:
                        w = [watson.readline()[1:-1],watson.readline()[:-1]]
                if w == []:
                    break
            out_seq = ''
            count +=1
            if not count%1000:
                print count
            for wb, cb in izip(w[1], c[1]):
                if wb == cb:
                    out_seq += wb
                elif wb == 'T' and cb == 'C':
                    out_seq += 'C'
                elif wb  == 'G' and cb == 'A':
                    out_seq += 'G'
                elif wb  == 'N':
                    out_seq += cb
                elif cb == 'N':
                    out_seq += wb
                else:
                    out_seq += 'N'
            record_out = '>%s\n%s\n'%(w[0],out_seq)
            out.write(record_out)
    except ValueError:
        return None


def main():
    "main function loop"
    args = parse_args()
    make_ref2(args)

if __name__ == '__main__':
    main()
