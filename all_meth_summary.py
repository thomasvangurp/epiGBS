__author__ = 'thomasvangurp'
import os
import vcf
import subprocess
from Bio import SeqIO
leak	=	"/Users/thomasvangurp/epiGBS/Baseclear/unfiltered_sequences/seqNNAtlE/Leak/output_mapping/methylation.bed"
Arabidopsis	=	"/Users/thomasvangurp/epiGBS/Baseclear/unfiltered_sequences/seqNNAtlE/Athal/output_mapping/methylation.bed"
Human	=	"/Users/thomasvangurp/epiGBS/Baseclear/unfiltered_sequences/seqNNAtlE/Human/output_mapping/methylation.bed"
Fallopia	=	"/Users/thomasvangurp/epiGBS/Baseclear/unfiltered_sequences/seqNNAtlE/Fallopia/output_mapping/methylation.bed"
Carrot	=	"/Users/thomasvangurp/epiGBS/Baseclear/unfiltered_sequences/seqNNAtlE/Carrot/output_mapping/methylation.bed"
Daphnia	=	"/Users/thomasvangurp/epiGBS/Baseclear/unfiltered_sequences/seqNNAtlE/Daphnia/output_mapping/methylation.bed"
Scabiosa	=	"/Users/thomasvangurp/epiGBS/Baseclear/unfiltered_sequences/seqNNAtlE/Scabiosa/output_mapping/methylation.bed"
Mimulus	=	"/Users/thomasvangurp/epiGBS/Baseclear/unfiltered_sequences/seqNNAtlE/Mimulus/output_mapping/methylation.bed"
Dandelion	=	"/Users/thomasvangurp/epiGBS/Baseclear/unfiltered_sequences /seqNNAtlE/dandelion/output_mapping/methylation.bed"
Lambda	=	"/Users/thomasvangurp/epiGBS/Baseclear/unfiltered_sequences/seqNNAtlE/Lambda/output_mapping/methylation.bed"
species = [Human,Arabidopsis,leak,Fallopia,Carrot,Daphnia,Scabiosa,Mimulus,Dandelion,Lambda]
sp_names = ['Homo sapiens','Arabidopsis thaliana','Allium porrum','Fallopia japonica','Daucus carrota','Daphnia magna',
            'Scabiosa columbaria','Mimulus guttatus','Taraxacum officinale','Phage lambda']

meth_dict = {}
cov_treshold = 50
min_meth = 3
min_ratio = 0.05


def count_strand(contig,position,min_ratio,dir,seq_dict,sample):
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
    bam = '/tmp/'+sample.rstrip('_')+'_'+strand+'.bam'
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
    for nt in nts:
        try:
            nt_count[nt]+=1
        except KeyError:
            pass
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
    if min(plus_count,min_count) > 25:
        return min(_min_ratio,plus_ratio)
    else:
        return ''

def split_sample(header,dir):
    """Split watson and crick bam file into sample specific reads groups"""
    for name in header[4:]:
        if name.endswith('_methylated'):
            sample_name = name.replace('_methylated','')
        for strand in ['watson','crick']:
            split_cmd = 'samtools view -h %s/%s.bam | grep "^@SQ\|^@PG\|%s"|samtools view -Shb - > /tmp/%s.bam'%\
                        (dir,strand,sample_name,sample_name+'_'+strand)
            index_cmd = 'samtools index /tmp/%s.bam'%(sample_name+'_'+strand)
            file_check = os.path.exists('/tmp/%s.bam'%(sample_name+'_'+strand))
            if file_check:
                continue
            for cmd in [split_cmd,index_cmd]:
                p = subprocess.Popen(cmd,stdout=subprocess.PIPE,stderr=subprocess.PIPE,
                             shell=True,executable='/bin/bash')
                exit_code = p.wait()
    return 0


for i,file in enumerate(species):
    name = sp_names[i]
    meth_dict[name] = {'CG':{'meth':0,'total':0},'CHG':{'meth':0,'total':0},'CHH':{'meth':0,'total':0}}
    vcf_location = '/'.join(file.split('/')[:-1]+['snp.vcf.gz'])
    vcf_input = vcf.Reader(open(vcf_location))
    file_handle = open(file,'r')
    header = file_handle.readline().split('\t')
    #Get sequence dict
    split_sample(header,'/'.join(file.split('/')[:-1]))
    fasta = '/'.join(file.split('/')[:-2])+'/output_denovo/'+'consensus_cluster.renamed.fa'
    seq_dict = {}
    seq_parser = SeqIO.parse(open(fasta,'r'),'fasta')
    for seq in seq_parser:
        seq_dict[seq.id] = str(seq.seq)
    for line in file_handle:
        split_line = line[:-1].split('\t')
        contig,position = split_line[:2]
        if int(position) < 6 or int(position) > (len(seq_dict[contig])-6):
            continue
        try:
            if vcf_input.fetch(int(contig),int(position)):
                # print "SNP position found in %s at %s"%(contig,position)
                continue
        except ValueError:
            pass
        for i in range(5,100,2):
            try:
                if i < len(split_line) and int(split_line[i]) > cov_treshold:
                    context = split_line[2]
                    sample = '_'.join(header[i].split('_')[:-1])
                    if int(split_line[i-1]) >= min_meth:
                        try:
                            # ratio = float(split_line[i-1]) / float(split_line[i])
                            # if ratio > min_ratio:
                            #     #check if observation should be used.
                            ratio = count_strand(contig,position,min_ratio,'/'.join(file.split('/')[:-1]),\
                                                 seq_dict,sample)
                            if ratio > min_ratio and type(ratio) == type(1.0):
                                meth_dict[name][context]['meth'] += 1
                            elif type(ratio) == type(1.0):
                                pass
                            else:
                                continue
                            meth_dict[name][context]['total'] += 1
                        except KeyError:
                            pass
                    else:
                        try:
                            meth_dict[name][context]['total'] += 1
                        except KeyError:
                            pass
                else:
                    break
            except ValueError:
                pass
    for k,v in sorted(meth_dict[name].items()):
        print '\t'.join([str(v) for v in [name,k,v['meth'],v['total']]])