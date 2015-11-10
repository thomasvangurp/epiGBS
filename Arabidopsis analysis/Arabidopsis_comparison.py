__author__ = 'thomasvangurp'
description = """Compare CG, CHG and CHH methylation between epiGBS data and WBGS data"""
import pysam
import os
from Bio import Restriction
from Bio import SeqIO
import vcf
import subprocess
import copy

WGBS_dir = '/Users/thomasvangurp/epiGBS/becker_nature_2011/analysis/'
epiGBS_dir = '/Users/thomasvangurp/epiGBS/Baseclear/unfiltered_sequences/seqNNAtlE/Athal/ref_mapping/'
ref_genome = '/Users/thomasvangurp/epiGBS/Baseclear/unfiltered_sequences/seqNNAtlE/Athal/TAIR10_ref/Arabidopsis_thaliana.TAIR10.22.dna.genome.fa'
cov_treshold = 10

def get_sites(ref_genome):
    """get list of all PstI locations in ref genome, put in dict(set)"""
    with open(ref_genome) as handle:
        seq_handle = SeqIO.parse(handle,'fasta')
        res_sites = {}
        for seq in seq_handle:
            import re
            res_sites[seq.id] = [m.start() for m in re.finditer('CTGCAG', str(seq.seq))]
            # to_add = set()
            # for site in res_sites[seq.id]:
            #     to_add.update(range(site+1,site+7))
            # res_sites[seq.id].update(to_add)
        #now exclude all sites that show evidence of >5% methylation in Arabidopsis
        # file_handle = open('/Users/thomasvangurp/epiGBS/becker_nature_2011/analysis/new_vcf/heatmap.igv','r')
        # header = file_handle.readline()
        # header += file_handle.readline()
        # to_add = {}
        # for line in file_handle:
        #     split_line = line[:-1].split('\t')
        #     pos = int(split_line[1])
        #     contig = split_line[0]
        #     if split_line[3] == 'CHG' and pos in res_sites[contig]:
        #         try:
        #             if max([float(v) for v in [k for n,k in enumerate(split_line) if n in [4,5,6,9] and k != "."] ]) > 0.05:
        #                 try:
        #                     index = res_sites[contig].index(pos)
        #                     pos2 = res_sites[contig][index+1]
        #                     for i in range(pos,pos2):
        #                         try:
        #                             to_add[contig].append(i)
        #                         except KeyError:
        #                             to_add[contig] = [i]
        #                 except IndexError:
        #                     pass
        #                 except ValueError:
        #                     continue
        #         except ValueError:
        #             pass
        # for k,v in to_add.items():
        #     res_sites[k]+=v
        handle.seek(0)
        seq_handle = SeqIO.parse(handle,'fasta')
        for seq in seq_handle:
            to_add = set()
            for site in [m.start() for m in re.finditer('CTGCAG', str(seq.seq))]:
                to_add.update(range(site+1,site+7))
            res_sites[seq.id] = set(res_sites[seq.id])
            res_sites[seq.id].update(to_add)
    return res_sites

def get_genome_dict(ref_genome):
    """Return genome as strings"""
    seq_dict = {}
    seq_parser = SeqIO.parse(open(ref_genome,'r'),'fasta')
    for seq in seq_parser:
        seq_dict[seq.id] = str(seq.seq)
    return seq_dict

def split_sample(header,dir,type):
    """Split watson and crick bam file into sample specific reads groups"""
    for name in header[4:]:
        if name.endswith('_methylated'):
            sample_name = name.replace('_methylated','')
        else:
            continue
        for strand in ['watson','crick']:
            split_cmd = 'samtools view -h %s/%s.bam | grep "^@SQ\|^@PG\|%s"|samtools view -Shb - > /tmp/%s_%s.bam'%\
                        (dir,strand,sample_name,sample_name+'_'+strand,type)
            index_cmd = 'samtools index /tmp/%s_%s.bam'%(sample_name+'_'+strand,type)
            file_check = os.path.exists('/tmp/%s_%s.bam'%(sample_name+'_'+strand,type))
            if file_check:
                continue
            for cmd in [split_cmd,index_cmd]:
                p = subprocess.Popen(cmd,stdout=subprocess.PIPE,stderr=subprocess.PIPE,
                             shell=True,executable='/bin/bash')
                exit_code = p.wait()
    return 0


def retain_symmetric(dict_in):
    """correct count based on positions for which there is evidence on both strands"""
    meth_dict = copy.deepcopy(dict_in)
    for ind in meth_dict.keys():
        for context,values in sorted(meth_dict[ind].items()):
            if context == 'CHG':
                interval = 2
            elif context == 'CG':
                interval = 1
            else:
                continue
            for type in values.keys():
                for chr,sub_dict in values[type].items():
                    for k,v in sub_dict.items():
                        if str(int(k)+interval) not in sub_dict \
                                and str(int(k)-interval) not in sub_dict:
                            sub_dict.pop(k)
            for type in values.keys():
                if type == 'meth':
                    continue
                correct_count = 0
                contigs = values[type].keys()
                out_dict = {}
                for contig in contigs:
                    i = 0
                    values_list = sorted([int(key) for key in values[type][contig].keys()])
                    try:
                        while True:
                            if values_list[i+1] == values_list[i] + interval:# and \
                                # values[type][contig][i][2] == values[type][contig][i+1][2]:
                                ratio_plus = values[type][contig][str(values_list[i])]
                                ratio_neg  = values[type][contig][str(values_list[i+1])]
                                position = values_list[i]
                                avg_ratio = (ratio_neg + ratio_plus) / 2
                                try:
                                    out_dict[contig][position] = avg_ratio
                                except KeyError:
                                    out_dict[contig] = {position:avg_ratio}
                                i+=2
                            else:
                                i+=1
                            if i >= len(values_list)-1:
                                break
                    except IndexError:
                        pass
            meth_dict[ind][context][type] = out_dict
    return meth_dict
def count_strand(contig,position,min_ratio,seq_dict,sample,type):
    """Determine if methylation is still above treshold given strand"""
    #1. if not both strands are represented with a minimum of 25 reads discard position
    #2. Both strands need to have a minimum of 2 "methylated" reads.
    #3. The minimum methylation ratio on both strands need to be > min_ratio
    nt_count = {'c':0,'C':0,'t':0,'T':0,'g':0,'G':0,'a':0,'A':0}
    ratio_plus = 0
    ratio_min = 0
    #get nucleotide
    nt = seq_dict[contig][int(position)-1]
    if nt == 'C':
        strand = 'watson'
    else:
        strand = 'crick'
    region = '%s:%s-%s'%(contig,position,position)
    bam = '/tmp/'+sample.rstrip('_')+'_'+strand+'_%s.bam'%type
    cmd = 'samtools mpileup -r %s %s'%(region,bam)
    p = subprocess.Popen(cmd,stdout=subprocess.PIPE,stderr=subprocess.PIPE,
                         shell=True,executable='/bin/bash')
    exit_code = p.wait()
    stdout = p.stdout.read().replace('\r','\n')
    try:
        nts = stdout.split('\t')[4]
    except IndexError:
        return ''
    # stderr = p.stderr.read()
    # if stderr:
    #     print stderr
    for nt in nt_count.keys():
        nt_count[nt]+= nts.count(nt)
    if strand == 'crick':
        plus_count = nt_count['A']+nt_count['G']
        min_count = nt_count['a']+nt_count['g']
        try:
            plus_ratio = nt_count['G'] / float(plus_count)
        except ZeroDivisionError:
            plus_ratio = ''
        try:
            _min_ratio = nt_count['g'] / float(min_count)
        except ZeroDivisionError:
            _min_ratio = ''
    else:
        plus_count = nt_count['T']+nt_count['C']
        min_count = nt_count['t']+nt_count['c']
        try:
            plus_ratio = nt_count['C'] / float(plus_count)
        except ZeroDivisionError:
            plus_ratio = ''
        try:
            _min_ratio = nt_count['c'] / float(min_count)
        except ZeroDivisionError:
            _min_ratio = ''
    if min(plus_count,min_count) > cov_treshold:
        #Here we will take the average of both estimates from forward and reverse reads
        return (_min_ratio + plus_ratio) / 2.0
    else:
        return ''

#get list of all PstI locations in ref genome, put in dict(set)
def main(input_dir,lib,res_sites):
    #loop over positions from methylation.bed based on epiGBS data
    min_meth = 1
    min_ratio = 0.05
    # res_sites = get_sites(ref_genome)
    meth_dict = {}
    vcf_location = input_dir + 'snp.vcf.gz'
    vcf_input = vcf.Reader(open(vcf_location))
    file_handle = open(input_dir+'methylation.bed','r')
    header = file_handle.readline()[:-1].split('\t')
    #split samples in watson and crick
    split_sample(header,input_dir,lib)
    #get sequences to know where watson and crick are
    seq_dict = get_genome_dict(ref_genome)
    for line in file_handle.readlines():
        split_line = line[:-1].split('\t')
        contig,position,context = split_line[:3]
        if context == '.':
            continue
        if int(position) in res_sites[contig]:
            #do not look at positions in restriction sites as they cannot be measured accurately with epiGBS
            continue
        try:
            var = vcf_input.fetch(int(contig),int(position))
            #Do not look at positions overlapping with SNPs
            if var:
                continue
            # if var:
            #     if var.QUAL > 50:
            #         continue
        except ValueError:
            pass
        for ind,i in zip(header[5:],range(5,len(split_line))):
            if ind.endswith('methylated'):
                continue
            ind_name = '_'.join(ind.split('_')[:2])
            if ind_name not in meth_dict:
                meth_dict[ind_name] = {'CG':{'meth':{},'total':{}},'CHG':{'meth':{},'total':{}},'CHH':{'meth':{},'total':{}}}
            try:
                if int(split_line[i]) > cov_treshold:
                    #only consider positions with coverage above treshold
                    context = split_line[2]
                    try:
                        ratio = float(split_line[i-1]) / float(split_line[i])
                        coverage = int(split_line[i])
                    except ValueError:
                        #no methylated site or one of the values equals none
                        ratio = 0
                    if lib == '!!epiGBS':
                        #for epiGBS libraries we have stringent criteria on meeting coverage criteria
                        # on both fw and reverse strand
                        #for WGBS data these criteria don't hold.
                        #get average ratio for both strands
                        if int(split_line[i]) >= (cov_treshold*2):
                            #we need a coverage higher or equal to coverage treshold on both strands
                            ratio = count_strand(contig,position,min_ratio,seq_dict,ind_name,lib)
                        else:
                            continue
                        if ratio == '':
                            #ratio of '' is returned if min coverage for both forward and reverse reads is not achieved
                            #so we do not count positions with insufficient coverage.
                            continue
                        # elif ratio > min_ratio and type(ratio) == type(1.0):
                        try:
                            meth_dict[ind_name][context]['meth'][contig][position] = ratio
                        except KeyError:
                            if meth_dict[ind_name][context]['meth'] != {}:
                                meth_dict[ind_name][context]['meth'][contig] = {position:ratio}
                            else:
                                meth_dict[ind_name][context]['meth'] = {contig:{position:ratio}}

                        try:
                            meth_dict[ind_name][context]['total'][contig][position]  = ratio
                        except KeyError:
                            if meth_dict[ind_name][context]['total'] != {}:
                                meth_dict[ind_name][context]['total'][contig] = {position:ratio}
                            else:
                                meth_dict[ind_name][context]['total'] = {contig:{position:ratio}}
                    else:
                        try:
                            meth_dict[ind_name][context]['total'][contig][position] = (ratio,coverage,)
                        except KeyError:
                            if meth_dict[ind_name][context]['total'] != {}:
                                meth_dict[ind_name][context]['total'][contig] = {position:(ratio,coverage,)}
                            else:
                                meth_dict[ind_name][context]['total'] = {contig:{position:(ratio,coverage,)}}
                else:
                    continue
            except ValueError:
                pass
    return meth_dict

def output_comparison(wgbs_meth,epi_meth,dir):
    """Generate output csv files for comparison between CG,CHG and CHH methylation"""
    for ind in epi_meth.keys():
        for context in epi_meth[ind].keys():
            out_file = open(output_dir+'%s_%s_comparison.csv'%(ind,context),'w')
            out_file.write('%s_epiGBS\t%s_wgbs\n'%(ind,ind))
            for chrom,positions in sorted(epi_meth[ind][context]['total'].items()):
                for pos in sorted([int(p) for p in positions]):
                    if pos not in epi_meth[ind][context]['total'][chrom]:
                        pos = str(pos)
                    epiGBS_ratio,epiGBS_total = epi_meth[ind][context]['total'][chrom][pos]
                    try:
                        wgbs_name = '30-%s_total'%ind.split('_')[1][1:]
                        wgbs_ratio,wgbs_total = wgbs_meth[wgbs_name][context]['total'][chrom][pos]
                    except KeyError:
                        #wgbs_ratio is not available, do not write output
                        continue
                    out_file.write('%s\t%.6f\t%s\t%.6f\t%s\n'%('%s\t%s'%(chrom,pos),epiGBS_ratio,epiGBS_total,wgbs_ratio,wgbs_total))
            out_file.close()

def comparison_30_31(wgbs_meth,output_dir):
    """Compare line 30 with line 31 meth"""
    for ind in wgbs_meth.keys():
        if '31' in ind:
            gen0 = ind.replace('31','30')
        else:
            continue
        for context in wgbs_meth[ind].keys():
            out_file = open(output_dir+'%s_%s_%s_comparison.csv'%(gen0,ind,context),'w')
            out_file.write('%s_epiGBS\t%s_wgbs\n'%(gen0,ind))
            for chrom,positions in sorted(wgbs_meth[ind][context]['total'].items()):
                for pos in sorted([int(p) for p in positions]):
                    if pos not in wgbs_meth[ind][context]['total'][chrom]:
                        pos = str(pos)
                    ind_ratio,ind_count = wgbs_meth[ind][context]['total'][chrom][pos]
                    try:
                        ind0_ratio,ind0_count = wgbs_meth[gen0][context]['total'][chrom][pos]
                    except KeyError:
                        #wgbs_ratio is not available, do not write output
                        continue
                    out_file.write('%s\t%.3f\t%s\t%.3f\t%s\n'%('%s\t%s'%(chrom,pos),ind0_ratio,ind0_count,ind_ratio,ind_count))
            out_file.close()

if __name__ == '__main__':
    input_dir = '/Users/thomasvangurp/epiGBS/Baseclear/unfiltered_sequences/seqNNAtlE/Athal/ref_mapping/'
    input_dir_becker = '/Users/thomasvangurp/epiGBS/becker_nature_2011/analysis/new_vcf/'
    res_sites = get_sites(ref_genome)
    wgbs_meth = main(input_dir_becker,'wgbs',res_sites)
    #output non-symmetric
    output_dir  = '/Users/thomasvangurp/epiGBS/Baseclear/unfiltered_sequences/seqNNAtlE/Athal/comparison/'
    comparison_30_31(wgbs_meth,output_dir)
    epi_meth = main(input_dir,'epiGBS',res_sites)
    output_comparison(wgbs_meth,epi_meth,output_dir)
    #output symmetric
    # output_dir  = '/Users/thomasvangurp/epiGBS/Baseclear/unfiltered_sequences/seqNNAtlE/Athal/comparison/symmetric/'
    # wgbs_meth_symmetric = retain_symmetric(wgbs_meth)
    # epi_meth_symmetric = retain_symmetric(epi_meth)
    # comparison_30_31(wgbs_meth_symmetric,output_dir)
    # output_comparison(wgbs_meth_symmetric,epi_meth_symmetric,output_dir)
