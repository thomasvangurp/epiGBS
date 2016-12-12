#!/usr/bin/env python
"""concatenate reference sequence"""
import argparse
import os
import pysam
from Bio import SeqIO
from nested_dict import nested_dict

def parse_args():
    """Pass command line arguments"""

    parser = argparse.ArgumentParser(description='Concatenate reference for display in jbrowse')
    parser.add_argument('-o', '--output_dir', help='output directory')
    parser.add_argument('-i', '--input_dir', help='input directory with bam files')
    parser.add_argument('-g', '--genes', help='xml for genes with hits')
    parser.add_argument('-r', '--ref', help='reference sequence in')
    args = parser.parse_args()
    if 'input_dir' in args:
        if os.path.exists(os.path.join(args.input_dir, 'watson.dedup.bam')):
            args.watson = os.path.join(args.input_dir, 'watson.dedup.bam')
        if os.path.exists(os.path.join(args.input_dir, 'crick.dedup.bam')):
            args.crick = os.path.join(args.input_dir, 'crick.dedup.bam')
    if not os.path.exists(args.output_dir):
        os.mkdir(args.output_dir)
    return args

def new_ref(coverage, args):
    """Get coverage of reference sequence"""
    #load current ref
    currenct_seq = SeqIO.to_dict(SeqIO.parse(args.ref,'fasta'))
    #get genes
    gene_list  = []
    ref_mapping = {'gene':{},'non-gene':{}}
    gene_ref = ''
    non_gene_ref = ''
    count = 0
    with open(args.genes) as handle:
        while True:
            # if count > 1000:
            #     break
            line = handle.readline()
            if not line:
                break
            if 'query-def' in line:
                count += 1
                query = line[line.index('>')+1:line.rindex('<')]
                gene_list.append(query)
    for key,contigs in sorted(coverage.items())[::-1]:
        for contig in contigs:
            contig_seq = str(currenct_seq[contig].seq)
            if contig in gene_list:
                ref_mapping[contig] = ('gene', len(gene_ref))
                gene_ref += contig_seq + 'N'*10
            else:
                ref_mapping[contig] = ('non_gene', len(non_gene_ref))
                non_gene_ref += contig_seq + 'N' * 10
    with open(os.path.join(args.output_dir,'ref.fa'),'w') as refout:
        if gene_ref != '':
            refout.write('>%s\n%s\n' % ('gene-related', gene_ref))
        refout.write('>%s\n%s' % ('other', non_gene_ref))
    #store length of sequences of ref for genes and non-genes
    ref_mapping['gene']['length'] = len(gene_ref)
    ref_mapping['non-gene']['length'] = len(non_gene_ref)
    return ref_mapping

def parse_bam(mapping_dict, args):
    """parse bam file to get mapping per contig and individual"""
    for bam in [args.watson, args.crick]:
        print 'start processing %s' % bam
        file_handle = pysam.AlignmentFile(bam)
        count = 0
        for read in file_handle:
            count += 1
            if not count%1000000:
                print '%s reads processed'%count
            #     break
            if read.is_qcfail or read.is_duplicate or read.is_supplementary or read.is_unmapped or \
                    (read.is_paired and not read.is_proper_pair):
                continue
            read_dict = dict(read.tags)
            try:
                mapping_dict[read.reference_name][read_dict['ST']][read_dict['RG']] += 1
            except TypeError:
                mapping_dict[read.reference_name][read_dict['ST']][read_dict['RG']] = 1
        print 'finished processing %s' % bam
    for contig in mapping_dict.keys():
        if 'Watson' not in mapping_dict[contig]:
            continue
        for sample in mapping_dict[contig]['Watson']:
            watson_coverage = mapping_dict[contig]['Watson'][sample]
            if sample in mapping_dict[contig]['Crick']:
                crick_coverage = mapping_dict[contig]['Crick'][sample]
            if min(watson_coverage,crick_coverage) >= 5:
                try:
                    mapping_dict['coverage'][contig]['coverage'] += watson_coverage + crick_coverage
                except TypeError:
                    mapping_dict['coverage'][contig]['coverage'] = watson_coverage + crick_coverage
                try:
                    mapping_dict['coverage'][contig]['samples'] += 1
                except TypeError:
                    mapping_dict['coverage'][contig]['samples'] = 1
                try:
                    mapping_dict['coverage'][contig][sample] += watson_coverage + crick_coverage
                except TypeError:
                    mapping_dict['coverage'][contig][sample] = watson_coverage + crick_coverage
    return mapping_dict




def coverage_ref(args):
    """Get coverage of reference sequence"""
    #parse watson
    mapping_dict = nested_dict()
    mapping_dict = parse_bam(mapping_dict, args)
    coverage_dict = {}
    for contig in mapping_dict['coverage'].keys():
        samples_covered = mapping_dict['coverage'][contig]['samples']
        try:
            coverage_dict[samples_covered].append(contig)
        except KeyError:
            coverage_dict[samples_covered] = [contig]
    return coverage_dict

def rewrite_bam(ref_mapping,args):
    """split bam file using pysam"""
    bam_watson_handle = pysam.AlignmentFile(args.watson)
    bam_crick_handle = pysam.AlignmentFile(args.crick)
    out_handles = {}
    #create new template
    header = bam_watson_handle.header
    header['SQ']=[{'LN':ref_mapping['gene']['length'],'SN':'gene'},
                  {'LN':ref_mapping['non-gene']['length'],'SN':'other'}]
    contig_index = {'gene':0,'non_gene':1}
    for item in bam_watson_handle.header['RG']:
        watson_path = os.path.join(args.output_dir, '%s.watson.bam' % item['SM'])
        crick_path = os.path.join(args.output_dir, '%s.crick.bam' % item['SM'])
        watson_handle = pysam.AlignmentFile(watson_path, "wb", header=header)
        crick_handle = pysam.AlignmentFile(crick_path, "wb", header=header)
        out_handles[item['SM']] = {'watson': watson_handle, 'crick': crick_handle}
    i = 0
    for read in bam_watson_handle:
        i += 1
        if not i % 100000:
            print 'processed %s reads' % i
        sample = '_'.join(dict(read.tags)['RG'].split('_')[2:])
        handle = out_handles[sample]['watson']
        #change read parameters depending on contig
        contig_name, contig_pos = ref_mapping[read.reference_name]
        read.reference_id = contig_index[contig_name]
        read.pos += contig_pos
        if read.is_paired:
            if read.is_proper_pair:
                read.pnext += contig_pos
        handle.write(read)
    for subdict in out_handles.values():
        subdict['watson'].close()
    i = 0
    for read in bam_crick_handle:
        i += 1
        if not i % 100000:
            print 'processed %s reads' % i
        sample = '_'.join(dict(read.tags)['RG'].split('_')[2:])
        handle = out_handles[sample]['crick']
        # change read parameters depending on contig
        contig_name, contig_pos = ref_mapping[read.reference_name]
        read.reference_id = contig_index[contig_name]
        read.pos += contig_pos
        if read.is_paired:
            if read.is_proper_pair:
                read.pnext += contig_pos
        handle.write(read)
    for subdict in out_handles.values():
        subdict['crick'].close()
    for item in bam_watson_handle.header['RG']:
        watson_path = os.path.join(args.output_dir, '%s.watson.bam' % item['SM'])
        crick_path = os.path.join(args.output_dir, '%s.crick.bam' % item['SM'])
        pysam.index(watson_path)
        pysam.index(crick_path)


def main():
    """main function
    :rtype: int
    """
    args = parse_args()
    coverage = coverage_ref(args)
    ref_mapping = new_ref(coverage, args)
    rewrite_bam(ref_mapping, args)
    return 0


if __name__ == '__main__':
    return_code = main()