#!/usr/bin/env python
from Bio import SeqIO
from itertools import izip
import pysam
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




def watson_crick_count(args, region):
    """Gets total count of watson and crick"""
    bam = pysam.Samfile(args.bam, 'rb')
    reads = bam.fetch(region)
    count_dict = {}
    for read in reads:
        count = int(read.qname.rstrip(';').split('=')[-1])
        type = read.tags[-1][-1]
        try:
            count_dict[type]+=count
        except KeyError:
             count_dict[type]=count
    str_out = ";wats_cov:%s;crick_cov:%s"%(count_dict['watson'], count_dict['crick'])
    return str_out

def make_ref2(args):
    """"Make a concensus from watson and crick"""
    watson = SeqIO.parse(open(args.watson,'r'),'fasta')
    crick = SeqIO.parse(open(args.crick,'r'),'fasta')
    out = open(args.consensus, 'w')
    try:
        while True:
            try:
                w = watson.next()
                c = crick.next()
            except StopIteration:
                break
            while True:
                if w.name == c.name:
                    break
                else:
                    try:
                        if int(w.name) > int(c.name):
                            c = crick.next()
                        else:
                            w = watson.next()
                    except StopIteration:
                        return 0
            out_seq = ''
            for wb, cb in izip(w.seq, c.seq):
                if wb.lower() == cb.lower():
                    out_seq += wb
                elif wb in 'tT' and cb in 'cC':
                    out_seq += 'C'
                elif wb in 'gG' and cb in 'aA':
                    out_seq += 'G'
                elif wb in 'nN':
                    out_seq += cb
                elif cb in 'nN':
                    out_seq += wb
                else:
                    print wb, cb
            str_out = watson_crick_count(args, w.name)
            record_out = '>%s\n%s\n'%(w.name+str_out, out_seq.upper())
            if out_seq:
                out.write(record_out)
    except ValueError:
        return None


def main():
    "main function loop"
    args = parse_args()
    make_ref2(args)

if __name__ == '__main__':
    main()
