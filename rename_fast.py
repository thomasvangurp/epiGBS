#!/usr/bin/python2.6
import sys
import re
import os
from optparse import OptionParser
"""Illumina sequence identifiers
Sequences from the Illumina software use a systematic identifier:
@HWUSI-EAS100R:6:73:941:1973#0/1
HWUSI-EAS100R	the unique instrument name
6	flowcell lane
73	tile number within the flowcell lane
941	'x'-coordinate of the cluster within the tile
1973	'y'-coordinate of the cluster within the tile
#0	index number for a multiplexed sample (0 for no indexing)
/1	the member of a pair, /1 or /2 (paired-end or mate-pair reads only)
Versions of the Illumina pipeline since 1.4 appear to use #NNNNNN instead of #0 for the multiplex ID, where NNNNNN is the sequence of the multiplex tag.

In This script names can be either a simple count or the 'numbers' referring to the flowcell etc.
"""

def main():
    parser = OptionParser()
    parser.add_option("-i", "--input",  metavar = "input",  action = "store",
                      type="string",  dest = "input", 
                      help = "input fastq default is STDIN")
    parser.add_option("-o", "--output",  metavar = "output",  action = "store",
                  type="string",  dest = "output", 
                  help = "Specify output file, default is STDOUT")
    parser.add_option("-g",  "--genotype",  action = "store", 
                      type="string",  default = "" ,  dest = "genotype", 
                      help = "Optional genotype code preceding read names GTYPE_Read1")
    parser.add_option("-n",  "--number",  metavar="number",  action="store_true", 
                      default="0",   dest="number", help="Rename sequence identifiers to numbers")
    parser.add_option("-k",  "--keep",  metavar="number",  action="store_true"
                          ,  dest="id", 
                          help="Rename sequence identifiers keeping Illumina flowcell,tile and coordinates")                      
    return parser

def parse_fastq(line, index, end):
    if not index:
        while line:
            try:
                if line.startswith('@'):
                    index += 1
                    file_out.write('@%s%s%s\n'%(genotype, index, end))
                    file_out.write(file_in.next())
                    line = file_in.next()
                else:
                    pass
                    #raise ValueError
                if line.startswith('+'):
                    file_out.write('+%s%s%s\n'%(genotype, index, end))
                    file_out.write(file_in.next())
                    line = file_in.next()
                else:
                    pass
                    #raise ValueError
            except StopIteration:
                break
    else:
        while line:
            try:
                if line.startswith('@'):
                    index = ':'.join(line.split('#')[0].split(':')[1:-1])
                    file_out.write('@%s%s%s\n'%(genotype, index, end))
                    file_out.write(file_in.next())
                    line = file_in.next()
                else:
                    pass
                    #raise ValueError
                if line.startswith('+'):
                    file_out.write('+%s%s%s\n'%(genotype, index, end))
                    file_out.write(file_in.next())
                    line = file_in.next()
                else:
                    pass
                    #raise ValueError
            except StopIteration:
                break
def parse_fasta(line, index):
    while line:
        try:
            if line.startswith('>'):
                index += 1
                file_out.write('>%s%s\n'%(genotype, index))
                file_out.write(file_in.next())
                line = file_in.next()
            else:
                file_out.write(line)
                line = file_in.next()
                pass
                #raise ValueError
        except StopIteration:
            break

if __name__ == "__main__":
    parser = main()
    opts, args = parser.parse_args()
    if opts.input:
        file_in = open(opts.input, 'r')
    else:
        file_in = sys.stdin
    if opts.output:
        file_out = open(opts.output, 'w')
    else:
        file_out = sys.stdout
    if opts.genotype == 'chr':
        genotype = "%s"%(opts.genotype)
    elif opts.genotype:
        genotype = "%s_"%(opts.genotype)
    else:
        genotype  = ''
    if int(opts.number):
        index = 0
    elif opts.id:
        index = 1
    try:
        line = file_in.next()
        if line.startswith('>'):
            index = 0 # For fasta no Illumina names are available
            parse_fasta(line, index)
        elif line.startswith('@') and line.endswith('/1\n') or line.endswith('/2\n'):
            end = line[-3:-1]
            parse_fastq(line, index, end)
        elif line.startswith('@') :
            end = ''
            parse_fastq(line, index, end)
        else:
            print "File not recognized, is it a value FASTA/Q file"
            
    except StopIteration:
        print "Please provide a non-empty file"
  

