#!/usr/bin/env python
# """make coverage stats for bam files in directory"""
import sys
import pysam
import os

input_dir = sys.argv[1]

watson_bam = pysam.AlignmentFile(os.path.join(input_dir,'watson.reheader.bam'),'rb')
crick_bam  = pysam.AlignmentFile(os.path.join(input_dir,'crick.reheader.bam'),'rb')

coverage_dict = {'sample_coverage':{'watson':{},'crick':{}},'contig_coverage':{}}
contig_list = [entry['SN'] for entry in watson_bam.header['SQ']]

count = 0
for read in watson_bam:
    count += 1
    if not count % 1000000:
        print "%s reads in watson parsed"%count
    if read.is_paired and read.is_read2:
        continue
    sample = dict(read.tags)['RG']
    if sample in ['watson', 'crick']:
        continue
    contig = watson_bam.get_reference_name(read.tid)
    try:
        coverage_dict['sample_coverage'][sample] += 1
    except KeyError:
        coverage_dict['sample_coverage'][sample] = 1
    try:
        coverage_dict['contig_coverage'][contig][sample]['watson'] += 1
    except KeyError:
        if contig not in coverage_dict['contig_coverage']:
            coverage_dict['contig_coverage'][contig] = {}
        if sample not in coverage_dict['contig_coverage'][contig]:
            coverage_dict['contig_coverage'][contig][sample] = {'watson':1}
for read in crick_bam:
    count += 1
    if not count % 1000000:
        print "%s reads in crick parsed"%count
    if read.is_paired and read.is_read2:
        continue
    sample = dict(read.tags)['RG']
    if sample in ['watson','crick']:
        continue
    contig = crick_bam.get_reference_name(read.tid)
    try:
        coverage_dict['sample_coverage'][sample] += 1
    except KeyError:
        coverage_dict['sample_coverage'][sample] = 1
    try:
        coverage_dict['contig_coverage'][contig][sample]['crick'] += 1
    except KeyError:
        if contig not in coverage_dict['contig_coverage']:
            coverage_dict['contig_coverage'][contig] = {}
        if sample not in coverage_dict['contig_coverage'][contig]:
            coverage_dict['contig_coverage'][contig][sample] = {'crick':1}
        else:
            coverage_dict['contig_coverage'][contig][sample]['crick'] =  1


#get list of coverage per contig
contig_cover = {}
for contig,subdict in coverage_dict['contig_coverage'].items():
    contig_cover_watson = sum([v['watson'] for v in subdict.values() if 'watson' in v])
    contig_cover_crick = sum([ v['crick'] for v in subdict.values() if 'crick' in v])
    contig_cover[contig] = min(contig_cover_watson, contig_cover_crick)

#sorted contig list by coverage
contig_list = [i[0] for i in sorted(contig_cover.items(), key=lambda x:x[1])[::-1]]


#make list of sample from high to low coverage. Make groups
sorted_samples = [i[0] for i in sorted(coverage_dict['sample_coverage'].items(),key=lambda x:x[1])[::-1] if i[0] not in ['watson','crick']]

min_contig_cover = {}


file_out = open('/tmp/out_stats.txt','w')
for contig in contig_list:
    min_contig_cover[contig] = {}
    out_line = [contig]
    sorted_coverage = [i for i in sorted(coverage_dict['contig_coverage'][contig].items(),key=lambda x:min(x[1].values()))[::-1]]
    for quartile in [0.25,0.5,0.75,1]:
        try:
            min_contig_cover[contig][quartile] = min(sorted_coverage[int(len(sorted_samples)*quartile)-1][1].values())
        except IndexError:
            min_contig_cover[contig][quartile] = min(sorted_coverage[-1][1].values())
        out_line.append('%s' % min_contig_cover[contig][quartile])
    file_out.write('\t'.join(out_line)+'\n')






