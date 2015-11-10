__author__ = 'thomasvangurp'
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pysam
import gzip
import copy

merged_file = '/Users/thomasvangurp/epiGBS/WUR/epiGBS_2015/join.assembled.trimmed.fastq'
bam_handle = pysam.AlignmentFile('/Users/thomasvangurp/epiGBS/WUR/epiGBS_2015/sorted.bam','rb')
watson = pysam.AlignmentFile('/Users/thomasvangurp/epiGBS/WUR/epiGBS_2015/watson.bam','wb',template=bam_handle)
crick = pysam.AlignmentFile('/Users/thomasvangurp/epiGBS/WUR/epiGBS_2015/crick.bam','wb',template=bam_handle)

reads = SeqIO.index(merged_file,'fastq')

def process_reads(read_list,read_index):
    """process reads and return correct bam record"""
    nm = 1000
    best_read = None
    for read in read_list:
        for tag in read.tags:
            if tag[0] == 'NM':
                if tag[1] < nm:
                    best_read = read
                    nm = tag[1]
                break
    if nm != 1000:
        if best_read.query_alignment_length/ float(best_read.query_length ) < 0.8:
            return 0
        correct_seq = str(read_index.get(best_read.qname)._seq)
        record = SeqRecord(Seq(correct_seq))
        if best_read.is_reverse:
            read_seq = str(record.reverse_complement().seq)
        else:
            read_seq = str(record.seq)
        #change sequence of bam record
        best_read.seq = read_seq
        if 'GA' in best_read.tags[-1][1]:
            crick.write(best_read)
        else:
            watson.write(best_read)
        return best_read

# out_list = []
# for read in bam_handle:
#     if read.flag > 16:
#         continue
#     name =  read.qname
#     if out_list == []:
#         out_list+=[read]
#     elif read.qname in [item.qname for item in out_list]:
#         out_list+=[read]
#     else:
#         #process reads
#         process_reads(out_list,reads)
#         out_list = [read]


def get_ratio(nts,strand,cov_treshold):
    """Calculate corrected ratio for position"""
    nt_count = {'c':0,'C':0,'t':0,'T':0,'g':0,'G':0,'a':0,'A':0}
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
        meth = nt_count['g'] + nt_count['G']
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
        meth = nt_count['c'] + nt_count['C']

    nometh = plus_count + min_count - meth
    try:
        ratio = meth/float(nometh+meth)
    except ZeroDivisionError:
        ratio = 0.0
    return meth,nometh,ratio

watson_WGBS = gzip.open('/Users/thomasvangurp/epiGBS/WUR/epiGBS_2015/29_csp6i_watson.pileup.gz')
watson_epiGBS = gzip.open('/Users/thomasvangurp/epiGBS/WUR/epiGBS_2015/watson.pileup.gz')
crick_WGBS = gzip.open('/Users/thomasvangurp/epiGBS/WUR/epiGBS_2015/29_csp6i_crick.pileup.gz')
crick_epiGBS = gzip.open('/Users/thomasvangurp/epiGBS/WUR/epiGBS_2015/crick.pileup.gz')

ref = '/Users/thomasvangurp/epiGBS/Baseclear/unfiltered_sequences/seqNNAtlE/Athal/TAIR10_ref/Arabidopsis_thaliana.TAIR10.22.dna.genome.fa'
genome_dict = SeqIO.to_dict(SeqIO.parse(open(ref),'fasta'))

def get_context(nts):
    """return context of nucleotide"""
    if len(nts) < 5:
        return None
    if nts[2] == 'C':
        if nts[3] == 'G':
            context = 'CG'
        elif nts[4] == 'G':
            context = 'CHG'
        else:
            context = 'CHH'
    else:
        if nts[1] == 'C':
            context = 'CG'
        elif nts[0] == 'C':
            context = 'CHG'
        else:
            context = 'CHH'
    return context


def get_positions(genome_dict,meth_dict,file,mincov):
    """make a methylation_dict for a file"""

    if 'watson' in file.filename:
        nt_seq = 'C'
    else:
        nt_seq = 'G'

    for i,line in enumerate(file):
        split_line = line.rstrip('\n').split('\t')
        chrom,pos,ref,cov,nts,qual = split_line
        pos = int(pos)
        cov = float(cov)
        ref = genome_dict[chrom][pos-1]
        if ref == nt_seq and cov > mincov:
            context = get_context(genome_dict[chrom][pos-3:pos+2])
            if not context:
                continue
            try:
                # plus_count = len([i for i in nts if i.upper() == i])
                # min_count = len([i for i in nts if i.lower() == i])
                # if min([plus_count,min_count]) < (mincov):
                #     continue
                # plus_ratio = nts.count(nt_seq.upper()) / float(plus_count)
                # min_ratio = nts.count(nt_seq.lower()) / float(min_count)
                # meth_ratio = (plus_ratio + min_ratio )/ 2
                meth_ratio = nts.upper().count(nt_seq) / cov
            except ZeroDivisionError:
                continue
            try:
                meth_dict[chrom][context][pos] = (meth_ratio,cov,)
            except KeyError:
                if chrom not in meth_dict:
                    meth_dict[chrom] = {context:{pos:(meth_ratio,cov,)}}
                else:
                    meth_dict[chrom][context] = {pos:(meth_ratio,cov,)}
    return meth_dict

def make_comparison(WGBS_meth,epiGBS_meth):
    """Get xy cover for WGBS vs epiGBS methylation"""
    cg_out = open('/Users/thomasvangurp/epiGBS/WUR/epiGBS_2015/cg.csv','w')
    chg_out = open('/Users/thomasvangurp/epiGBS/WUR/epiGBS_2015/chg.csv','w')
    chh_out = open('/Users/thomasvangurp/epiGBS/WUR/epiGBS_2015/chh.csv','w')
    for chrom,context_dict in sorted(WGBS_meth.items()):
        for context,pos_dict in sorted(context_dict.items()):
            for pos,value in sorted(pos_dict.items()):
                try:
                    epiGBS_value,cov =  epiGBS_meth[chrom][context][pos]
                except KeyError:
                    continue
                output = '\t'.join([str(v) for v in [chrom,pos,value[0],value[1],epiGBS_value,cov]]) + '\n'
                if context == 'CG':
                    cg_out.write(output)
                elif context == 'CHG':
                    chg_out.write(output)
                else:
                    chh_out.write(output)

def retain_symmetric(dict_in):
    """correct count based on positions for which there is evidence on both strands"""
    meth_dict = copy.deepcopy(dict_in)
    out_dict = {}
    for chrom in meth_dict.keys():
        for context,values in sorted(meth_dict[chrom].items()):
            if context == 'CHG':
                interval = 2
            elif context == 'CG':
                interval = 1
            else:
                continue
            for type in values.keys():
                for k,v in values.items():
                    if k+interval not in values \
                            and k-interval not in values:
                        values.pop(k)
            i = 0
            values_list = sorted(values.keys())
            try:
                while True:
                    if values_list[i+1] == values_list[i] + interval:# and \
                        # values[type][contig][i][2] == values[type][contig][i+1][2]:
                        ratio_plus,count_plus = values[values_list[i]]
                        ratio_neg,count_neg = values[values_list[i+1]]
                        count_sum = count_neg + count_plus
                        position = values_list[i]
                        avg_ratio = (ratio_neg + ratio_plus) / 2
                        try:
                            out_dict[chrom][context][position] = (avg_ratio,count_sum,)
                        except KeyError:
                            if chrom not in out_dict:
                                out_dict[chrom] = {}
                            out_dict[chrom][context] = {position:(avg_ratio,count_sum,)}
                        i+=2
                    else:
                        i+=1
                    if i >= len(values_list)-1:
                        break
            except IndexError:
                pass
    meth_dict = out_dict
    return meth_dict

mincov = 10
#process epiGBS pileup
epiGBS_meth = {}
epiGBS_meth = get_positions(genome_dict,epiGBS_meth,watson_epiGBS,mincov)
epiGBS_meth = get_positions(genome_dict,epiGBS_meth,crick_epiGBS,mincov)
WGBS_meth = {}
WGBS_meth = get_positions(genome_dict,WGBS_meth,watson_WGBS,mincov)
WGBS_meth = get_positions(genome_dict,WGBS_meth,crick_WGBS,mincov)
# WGBS_meth_symmetric = retain_symmetric(WGBS_meth)
# epiGBS_meth_symmetric = retain_symmetric(epiGBS_meth)
make_comparison(WGBS_meth,epiGBS_meth)



