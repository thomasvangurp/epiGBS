#!/usr/bin/env python
"""pypy only merge watson and crick calls to custom format"""
import argparse
import gzip
import os

def parse_args():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(description='Process input files')
    parser.add_argument('-w', '--watson', type=str, default=None,
                        help='watson (top strand) .vcf file input.')
    parser.add_argument('-c', '--crick', type=str, default=None,
                        help='crick (bottom strand) .vcf file input.')
    parser.add_argument('-o', '--output', type=str, default=None,
                        help='output custom tabular format')
    args = parser.parse_args()
    return args

def merge_line(watson_line,crick_line):
    """merge watson and crick output"""
    if 'DP=0;' in watson_line[7] and 'DP=0;' in crick_line[7]:
        return None
    out_line = watson_line[:2]
    out_line += [watson_line[3]]
    out_line += [watson_line[4].replace('<*>','').rstrip(',')]
    out_line += [crick_line[4].replace('<*>','').rstrip(',')]
    if out_line[-2:] == ['',''] and out_line[2] in 'AT':
        return None
    out_line += ['']*len(watson_line[9:])
    AD_index = watson_line[8].split(':').index('AD')
    watson_nt_index = out_line[3].split(',')
    crick_nt_index = out_line[4].split(',')
    for nt in 'ACGT':
        nt_pos_watson = None
        nt_pos_crick = None
        if nt == out_line[2]:
            nt_pos_watson = 0
            nt_pos_crick = 0
        if nt in watson_nt_index:
            nt_pos_watson = watson_nt_index.index(nt) + 1
        if nt in crick_nt_index:
            nt_pos_crick = crick_nt_index.index(nt) + 1
        for index,(w_value,c_value) in enumerate(zip(watson_line[9:],crick_line[9:])):
            try:
                watson_obs = w_value.split(':')[AD_index].split(',')[nt_pos_watson]
            except TypeError:
                watson_obs = 0
            try:
                crick_obs = c_value.split(':')[AD_index].split(',')[nt_pos_crick]
            except TypeError:
                crick_obs = 0
            if nt == 'T':
                out_line[index+5] += '%s,%s' % (watson_obs, crick_obs)
            else:
                out_line[index + 5] += '%s,%s:' % (watson_obs, crick_obs)
    return '\t'.join(out_line) + '\n'

def merge(args):
    """"merge watson and crick calls"""
    watson_handle = os.popen("pigz -cd %s" % args.watson)
    crick_handle = os.popen("pigz -cd %s" % args.crick)
    read_watson = True
    read_crick = True
    #TODO: define output header, should include sample names
    output = open(args.output, 'w')
    count = 0
    chrom_pos = []
    while True:
        if read_watson:
            watson_line = watson_handle.readline()
        if read_crick:
            crick_line = crick_handle.readline()
        if watson_line.startswith('#CHROM'):
            watson_header = watson_line[:-1].split('\t')
            output.write(watson_line)
            read_watson = False
        if crick_line.startswith('#CHROM'):
            crick_header = crick_line[:-1].split('\t')
            read_crick = False
        if read_watson == False and read_crick == False:
            break
    read_watson = True
    read_crick = True
    while True:
        if read_watson:
            while True:
                watson_line = watson_handle.readline()[:-1].split('\t')
                if watson_line[0] not in chrom_pos:
                    chrom_pos.append(watson_line[0])
                if 'INDEL' not in watson_line:
                    break
        if read_crick:
            while True:
                crick_line = crick_handle.readline()[:-1].split('\t')
                if crick_line[0] not in chrom_pos:
                    chrom_pos.append(crick_line[0])
                if 'INDEL' not in crick_line:
                    break
        if len(watson_line) < 2 or len(crick_line) < 2:
            break
        elif watson_line[0] == crick_line[0]:
            if watson_line[1] == crick_line[1]:
                count += 1
                if not count % 1000000:
                    print('processed %s lines' % count)
                output_line = merge_line(watson_line, crick_line)
                if output_line:
                    output.write(output_line)
                #start reading both lines again.
                read_watson = True
                read_crick = True
                continue
            elif int(watson_line[1]) > int(crick_line[1]):
                read_watson = False
                read_crick = True
            else:
                read_crick = False
                read_watson = True
        elif chrom_pos.index(watson_line[0]) > chrom_pos.index(crick_line[0]):
            read_watson = False
            read_crick = True
        else:
            read_crick = False
            read_watson = True
    output.close()
    os.popen('pigz -f %s' % args.output)
def main():
    """Main function loop"""
    args = parse_args()
    merge(args)
    
    
    
if __name__ == '__main__':
    main()