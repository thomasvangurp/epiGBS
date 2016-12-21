#!/usr/bin/env python
"""Make sequence context specific bigbed per sample"""
import argparse
import os
from nested_dict import nested_dict
from Bio import SeqIO

def parse_args():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(description='Concatenate reference for display in jbrowse')
    parser.add_argument('-f', '--fasta', help='concatenated reference sequence')
    parser.add_argument('-i', '--input', help='methylation.bed input')
    parser.add_argument('-b', '--bed', help='bed file with contig to concatenated mapping')
    parser.add_argument('-o', '--outputdir', help='output directory for bigwig files')
    parser.add_argument('-s', '--samples', help='make group statistics given sample distribution')
    args = parser.parse_args()
    if not os.path.exists(args.outputdir):
        os.mkdir(args.outputdir)
    return args

def get_mapping(args):
    """parse bed file for mapping contigs to concatenated ref"""
    mapping_dict = {}
    with open(args.bed,'r') as handle:
        for line in handle:
            target,target_start,target_end,name = line.rstrip('\n').split('\t')
            mapping_dict[name] = (target,target_start)
    return mapping_dict

def get_groups(args):
    """get groups defined in sample file"""
    group_dict = nested_dict()
    with open(args.samples) as handle:
        header = handle.readline().rstrip('\n').split(',')
        for line in handle:
            split_line = line.rstrip('\n').split(',')
            sample = split_line[0]
            for name,item in zip(header[2:],split_line[2:]):
                try:
                    group_dict[name][item].append(sample)
                except AttributeError:
                    group_dict[name][item] = [sample]
    return group_dict

def check_context(cluster,position,ref,context):
    """Check if context matches specified position"""
    position = int(position)
    nt = ref[str(cluster)][int(position)]
    if nt.upper() == 'C':
        try:
            up_2 = ref[cluster][int(position)+1:int(position)+3]
            if up_2[0] == 'G':
                actual_context = 'CG'
            elif up_2[1] == 'G':
                actual_context = 'CHG'
            else:
                actual_context = 'CHH'
        except IndexError:
            # TODO: replace by expected nucleotides from enz recognition site
            return 0
    elif nt.upper() == 'G':
        try:
            down_2 = ref[cluster][int(position)-2:int(position)]
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

def make_bed_graph(mapping_dict, groups, args):
    """get methylation ratio of invididuals and groups"""
    count = 0
    with open(args.input) as input_handle:
        header = input_handle.readline().rstrip('\n').split('\t')
        #make file handles
        file_handles = {}
        ref = SeqIO.to_dict(SeqIO.parse(args.fasta,'fasta'))
        for context in ['CG','CHG','CHH']:
            for i in range(4, len(header), 2):
                sample = header[i][:-11]
                handle = open(os.path.join(args.outputdir,'%s.%s.bedGraph'%(sample,context)),'w')
                # handle.write('track type=bedGraph\n')
                file_handles['%s_%s' % (sample,context)] = handle
            for category, sub_group in groups.items():
                for value, list in sub_group.items():
                    group_name = '%s_%s' % (category, value)
                    handle = open(os.path.join(args.outputdir, '%s.%s.bedGraph' % (group_name,context)),'w')
                    # handle.write('track type=bedGraph\n')
                    file_handles['%s_%s' % (group_name,context)] = handle
        for line in input_handle:
            count += 1
            if not count%1000000:
                print '%s lines processes' % count
            split_line = line.rstrip('\n').split('\t')
            pos_in_contig = split_line[1]
            context = split_line[2]
            concat_contig, start_contig = mapping_dict[split_line[0].replace('chr', '')]
            final_position = str(int(start_contig) + int(pos_in_contig) - 1)
            nt = str(ref[concat_contig].seq)[int(final_position)]
            if check_context(concat_contig, final_position, ref, context) != 0:
                # print "Skipping cluster %s for position %s" % (cluster, position)
                continue
            meth_dict = {}
            for i in range(4,len(header),2):
                sample = header[i][:-11]
                try:
                    meth_ratio = float(split_line[i]) / float(split_line[i+1])
                except ValueError:
                    meth_ratio = None
                meth_dict[sample] = meth_ratio
            for category,sub_group in groups.items():
                for value,list in sub_group.items():
                    try:
                        meth_ratio = sum([meth_dict[e] for e in list if meth_dict[e]])/\
                                 float(len([e for e in list if meth_dict[e]]))
                    except ZeroDivisionError:
                        meth_ratio = None
                    meth_dict['%s_%s' % (category, value)] = meth_ratio
            # concat_contig,position = mapping_dict[split_line[0].replace('chr','')]
            # position = str(int(position) + int(split_line[1]) - 1)
            context = split_line[2]
            for key,value in meth_dict.items():
                if not value:
                    continue
                if nt == 'G':
                    value *= -1
                handle = file_handles['%s_%s' % (key,context)]
                out = [concat_contig,final_position,str(int(final_position)+1),'%.4f'%value]
                handle.write('\t'.join(out) + '\n')
        for name,handle in file_handles.items():
            handle.close
            print '%s closed' % handle

def bed_graph_to_bigwig(args):
    """convert bed graphs to bigwig"""
    #make chrom.sizes!
    chrom_sizes = os.path.join(args.outputdir,'chrom.sizes')
    out_handle = open(chrom_sizes,'w')
    with open(args.fasta) as handle:
        while True:
            seq_name = handle.readline().rstrip('\n')[1:]
            if not seq_name:
                break
            seq_len = len(handle.readline().rstrip('\n'))
            out_handle.write('%s\t%s\n' % (seq_name,seq_len))
        out_handle.close()
    for file in os.listdir(args.outputdir):
        if not file.endswith('bedGraph'):
            continue
        file_in = os.path.join(args.outputdir, file)
        file_sorted = os.path.join(args.outputdir, file+'.sorted')
        name = '.'.join(file.split('.')[:-2])
        context = file.split('.')[-2]
        file_out = os.path.join(args.outputdir, '%s.bw.%s' % (name, context.lower()))
        os.system('sort -k1,1 -k2,2n %s > %s' % (file_in, file_sorted))
        os.system('bedGraphToBigWig %s %s %s' % (file_sorted, chrom_sizes, file_out)) #-itemsPerSlot=1
        os.system('rm %s %s' % (file_in, file_sorted))
    return 0

def main():
    """main function"""
    args = parse_args()
    mapping_dict = get_mapping(args)
    groups = get_groups(args)
    make_bed_graph(mapping_dict, groups, args)
    bed_graph_to_bigwig(args)
    return 0


if __name__ == '__main__':
    return_code = main()