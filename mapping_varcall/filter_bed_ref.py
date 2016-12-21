__author__ = 'thomasvangurp'
"""filter bed file by pct coverage"""
import sys
import argparse
from Bio import SeqIO

def parse_args():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(description='Process input files')
    parser.add_argument('-r', '--reference', type=str, help='reference genome input.')
    parser.add_argument('-b', '--bed', help='bed file input')
    parser.add_argument('-o','--bedout', help='bed file output')
    parser.add_argument('--refout', help='subsampled reference')
    parser.add_argument('--context',default="CG,CHG",help='subsampled reference')
    parser.add_argument('-p','--percent',default="0.8",help='percent of samples called')
    parser.add_argument('-c','--mincover',default="3",help='minimum coverage for valid call')
    args = parser.parse_args()
    return args

def check_context(cluster,position,ref,context):
    """Check if context matches specified position"""
    nt = ref[str(int(cluster))].upper()[int(position)-1]
    assert nt in 'NACGT'
    if nt == 'C':
        try:
            up_2 = ref[str(int(cluster))].upper()[int(position):int(position)+2]
            if up_2[0] == 'G':
                actual_context = 'CG'
            elif up_2[1] == 'G':
                actual_context = 'CHG'
            else:
                actual_context = 'CHH'
        except IndexError:
            # TODO: replace by expected nucleotides from enz recognition site
            return 0
    elif nt == 'G':
        try:
            down_2 = ref[str(int(cluster))].upper()[int(position)-3:int(position)-1]
            if down_2[1] == 'C':
                actual_context = 'CG'
            elif down_2[0] == 'C':
                actual_context = 'CHG'
            else:
                actual_context = 'CHH'
        except IndexError:
            return 0
    else:
        print 'Not a C nucleotide!'
        return 1
    if actual_context != context:
        # print "Context mismatch!"
        return 1
    else:
        return 0

def filter_bed(args, ref_dict):
    """Filter bed file on percentage of sites being called and context"""
    bed_in_handle = open(args.bed,'r')
    bed_out_handle = open(args.bedout,'w')
    header = bed_in_handle.readline()
    bed_out_handle.write(header)
    header = header.split('\t')
    min_call = round(len(header[4:])/2 * float(args.percent))
    min_depth = int(args.mincover)
    chrom_out = set()
    for line in bed_in_handle:
        split_line = line.split('\t')
        context = split_line[2]
        cluster, position = split_line[:2]
        cluster = float(cluster)
        position = float(position)
        if check_context(cluster, position, ref_dict, split_line[2]) != 0:
            # print "Skipping cluster %s for position %s" % (cluster, position)
            continue
        if context not in args.context:
            continue
        # samples_called = int(split_line[3])
        # if samples_called >= min_call:
        #     continue
        called = 0
        for i in range(4,len(header),2):
            try:
                count = split_line[i+1].rstrip('\n')
                total_calls = int(count)
                if total_calls >= min_depth:
                    called += 1
            except ValueError:
                pass
            except IndexError:
                called = 0
                break
        if called >= min_call:
            chrom_out.update(['chr%s'%split_line[0]])
            if not line.startswith('chr'):
                line = 'chr' + line
            bed_out_handle.write(line)
    return chrom_out


def clean_fasta(fasta_input,fasta_output,seqs):
    """clean fasta file to get only reference sequences which are called"""
    fasta_input = SeqIO.parse(open(fasta_input,'r'),'fasta')
    out_handle = open(fasta_output,'w')
    for seq in fasta_input:
        name = 'chr%s'%seq.name
        if name in seqs:
            out = '>%s\n%s\n'%(name,str(seq.seq).upper())
            out_handle.write(out)
    out_handle.flush()
    out_handle.close()



def main():
    """main function loop"""
    args = parse_args()
    ref_dict = SeqIO.to_dict(SeqIO.parse(open(args.reference), 'fasta'))
    chrom_out = filter_bed(args,ref_dict)
    clean_fasta(args.reference,args.refout,chrom_out)


    return 0

if __name__ == '__main__':
    main()