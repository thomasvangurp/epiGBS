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
        os.mkdir(args.output_dir + '/bam/')
    return args

def new_ref(coverage, args):
    """Get coverage of reference sequence"""
    #load current ref
    currenct_seq = SeqIO.to_dict(SeqIO.parse(args.ref,'fasta'))
    #get genes
    gene_list  = []
    ref_mapping = {'gene':{},'non_gene':{}}
    gene_ref = ''
    non_gene_ref = ''
    count = 0
    with open(args.genes,'r') as handle:
        while True:
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
                ref_mapping[contig] = ('gene', len(gene_ref), len(contig_seq))
                gene_ref += contig_seq + 'N'*10
            else:
                ref_mapping[contig] = ('non_gene', len(non_gene_ref), len(contig_seq))
                non_gene_ref += contig_seq + 'N' * 10
    with open(os.path.join(args.output_dir,'ref.fa'),'w') as refout:
        if gene_ref != '':
            refout.write('>%s\n%s\n' % ('gene', gene_ref))
        refout.write('>%s\n%s' % ('non_gene', non_gene_ref))
    #store length of sequences of ref for genes and non-genes
    ref_mapping['gene']['length'] = len(gene_ref)
    ref_mapping['non_gene']['length'] = len(non_gene_ref)
    return ref_mapping

def make_bed(ref_mapping,args):
    """"make bed for positions of contigs vs concatenated ref for jbrowse"""
    bed_out = os.path.join(args.output_dir,'contigs_vs_concatenated.bed')
    with open(bed_out,'w') as bed_out_handle:
        keys = sorted([int(k) for k in ref_mapping.keys() if k.isdigit()])
        for key in keys:
            contig = str(key)
            concat_contig,concat_ref_pos,contig_len = ref_mapping[contig]
            output = [concat_contig,str(concat_ref_pos),str(concat_ref_pos + contig_len),contig]
            bed_out_handle.write('\t'.join(output)+'\n')
    return 0


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
            mapping_dict['coverage'][contig]['samples'] = 0
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
            else:
                if 'samples' in mapping_dict['coverage'][contig]:
                    continue
                else:
                    mapping_dict['coverage'][contig]['samples'] = 0
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
                  {'LN':ref_mapping['non_gene']['length'],'SN':'non_gene'}]
    contig_index = {'gene':0,'non_gene':1}
    for item in bam_watson_handle.header['RG']:
        watson_path = os.path.join(args.output_dir, '%s.watson.tmp' % item['SM'])
        crick_path = os.path.join(args.output_dir, '%s.crick.tmp' % item['SM'])
        watson_handle = pysam.AlignmentFile(watson_path, "wb", header=header)
        crick_handle = pysam.AlignmentFile(crick_path, "wb", header=header)
        out_handles[item['SM']] = {'watson': watson_handle, 'crick': crick_handle}
    i = 0
    print 'start splitting Watson reads'
    for read in bam_watson_handle:
        i += 1
        if not i % 1000000:
            print 'processed %s reads' % i
        sample = '_'.join(dict(read.tags)['RG'].split('_')[2:])
        handle = out_handles[sample]['watson']
        #change read parameters depending on contig
        try:
            contig_name, contig_pos, contig_len = ref_mapping[read.reference_name]
        except KeyError:
            print '%s not found, continue nevertheless' % read.reference_name
            continue
        read.rname = contig_index[contig_name]
        read.mrnm = contig_index[contig_name]
        read.pos += contig_pos
        if read.is_paired:
            if read.is_proper_pair:
                read.pnext += contig_pos
        handle.write(read)
    for subdict in out_handles.values():
        subdict['watson'].close()
    i = 0
    print 'start splitting Crick reads'
    for read in bam_crick_handle:
        i += 1
        if not i % 1000000:
            print 'processed %s reads' % i
        sample = '_'.join(dict(read.tags)['RG'].split('_')[2:])
        handle = out_handles[sample]['crick']
        # change read parameters depending on contig
        try:
            contig_name, contig_pos, contig_len = ref_mapping[read.reference_name]
        except KeyError:
            print '%s not found, continue nevertheless' % read.reference_name
            continue
        read.rname = contig_index[contig_name]
        read.mrnm = contig_index[contig_name]
        read.pos += contig_pos
        if read.is_paired:
            if read.is_proper_pair:
                read.pnext += contig_pos
        handle.write(read)
    for subdict in out_handles.values():
        subdict['crick'].close()
    if not os.path.exists(os.path.join(args.output_dir,'bam')):
        os.mkdir(os.path.join(args.output_dir,'bam'))
    for item in bam_watson_handle.header['RG']:
        watson_tmp = os.path.join(args.output_dir, '%s.watson.tmp' % item['SM'])
        watson_tmp2 = os.path.join(args.output_dir, '%s.watson.tmp2' % item['SM'])
        watson_path = os.path.join(args.output_dir,'bam', '%s.watson.bam' % item['SM'])
        crick_tmp = os.path.join(args.output_dir, '%s.crick.tmp' % item['SM'])
        crick_tmp2 = os.path.join(args.output_dir, '%s.crick.tmp2' % item['SM'])
        crick_path = os.path.join(args.output_dir,'bam', '%s.crick.bam' % item['SM'])
        pysam.sort(watson_tmp,'-o', watson_tmp2)
        pysam.sort(crick_tmp,'-o', crick_tmp2)
        os.system('samtools calmd -b  %s %s > %s 2>/dev/null'
                  % (watson_tmp2, os.path.join(args.output_dir,'ref.fa'), watson_path))
        os.system('samtools calmd -b  %s %s > %s 2>/dev/null'
                  % (crick_tmp2, os.path.join(args.output_dir,'ref.fa'), crick_path))
        pysam.index(watson_path)
        pysam.index(crick_path)
    os.system('rm %s/*.tmp*'%(args.output_dir))


def main():
    """main function
    :rtype: int
    """
    args = parse_args()
    coverage = coverage_ref(args)
    #determine mapping and make new reference
    ref_mapping = new_ref(coverage, args)
    #make bed file with positions of contigs
    make_bed(ref_mapping,args)
    rewrite_bam(ref_mapping, args)
    return 0


if __name__ == '__main__':
    return_code = main()