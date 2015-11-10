__author__ = 'thomasvangurp'
import os
from Bio import SeqIO
import subprocess
from scipy import stats
from operator import itemgetter
"""Detect DMPs between all samples, calculate their abundance and distribution"""
# Choose 2 most abundant individuals for DMP detection
# How many DMP's are there per context?
#

def get_context(file_in):
    """Get genomic context from bed file"""
    output_dict = {}
    handle = open(file_in,'r')
    header = handle.readline().rstrip('\n').split('\t')
    for line in handle:
        split_line = line.rstrip('\n').split('\t')
        contig,pos,context = split_line[:3]
        if context == '.':
            continue
        try:
            output_dict[contig][pos] = context
        except KeyError:
            output_dict[contig] = {pos:context}
    return output_dict
def make_pileup(header,dir,sample_dict):
    """make pileup"""
    for name in header[4:]:
        if name.endswith('_methylated'):
            sample_name = name.replace('_methylated','')
            if sample_name not in sample_dict['watson']:
                continue
        else:
            continue
        for strand in ['watson','crick']:
            output_loc = {}
            output_loc['bamfile'] = '/tmp/%s.bam'%(sample_name+'_'+strand)
            output_loc['output'] = os.path.join(dir,'%s.pileup'%(sample_name+'_'+strand))
            output_loc['regions'] = os.path.join(dir,'%s_%s.bed'%(sample_name,strand))
            pileup_cmd = 'samtools mpileup -l %(regions)s %(bamfile)s > %(output)s'%output_loc
            #only run command if file does not yet exist
            file_check = os.path.exists(output_loc['output'])
            if file_check:
                continue
            p = subprocess.Popen(pileup_cmd,stdout=subprocess.PIPE,stderr=subprocess.PIPE,
                             shell=True,executable='/bin/bash')
            exit_code = p.wait()
            stdout = p.stdout.read()
            stderr = p.stderr.read()
            if stderr:
                print stderr,pileup_cmd
    return 0
def make_bed_files(sample_dict,out_dir):
    """generate bed file with positions to check"""
    for strand,subdict in sample_dict.items():
        for sample,contig_dict in subdict.items():
            file_name = '%s_%s.bed'%(sample,strand)
            out_handle = os.path.join(out_dir,file_name)
            with open(out_handle,'w') as file_out:
                for contig,positions in sorted(contig_dict.items()):
                    for pos in sorted(positions):
                        out_items = (contig,pos,)
                        file_out.write('%s\t%s\n'%out_items)
    return 0

def filter_header(header,sample_total):
    """"returns index of header given sample list"""
    include_samples = open('/Users/thomasvangurp/epiGBS/Baseclear/unfiltered_sequences/seqNNAtlE/selection_DMP.txt','r')
    samples = include_samples.read().split('\n')
    max_cover = {}
    for ind,coverage in sample_total.items():
        if ind in samples:
            max_cover[coverage] = ind
    samples_out = []
    for i in range(4):
        try:
            key = max(max_cover.keys())
        except ValueError:
            break
        max_sample = max_cover[key]
        samples_out.append(max_cover[key])
        max_cover.pop(key)
    return samples_out

def process_bed(input_bed,genome_dict,mincov):
    """Process input methylation.bed and generate sample specific watson and crick bed files"""
    #could be limited to samples of interest if required.
    #limited by filter
    handle = open(input_bed,'r')
    header = [i.replace('_total','') for i in handle.readline()[:-1].split('\t')]
    sample_dict = {'watson':{},'crick':{}}
    sample_total = {}
    for linenumber,line in enumerate(handle):
        if not linenumber%100000 and linenumber:
            print 'Processed %s lines'%linenumber
        split_line = line.rstrip('\n').split('\t')
        contig = split_line[0]
        pos = int(split_line[1])
        if genome_dict[contig].seq[pos-1] == 'C':
            strand = 'watson'
        elif genome_dict[contig].seq[pos-1] == 'G':
            strand = 'crick'
        else:
            print "Nucleotide incorrect, %s"%genome_dict[contig].seq[pos-1]
        #totals are separated by 2 positions.
        for i in range(5,len(split_line),2):
            total = split_line[i]
            if total == 'None':
                continue
            name = header[i]
            try:
                sample_total[name] += int(total)
            except KeyError:
                sample_total[name] = int(total)
            try:
                if int(total) > (mincov * 2):
                    try:
                        sample_dict[strand][name][contig].append(pos)
                    except KeyError:
                        if name not in sample_dict[strand]:
                            sample_dict[strand][name] = {contig:[pos]}
                        elif contig not in sample_dict[strand][name]:
                            sample_dict[strand][name][contig] = [pos]
            except ValueError:
                pass
    #filter to only include the samples that we would like to get
    samples_to_include = filter_header(header,sample_total)
    for strand,subdict in sample_dict.items():
        for sample in subdict.keys():
            if sample not in samples_to_include:
                subdict.pop(sample)
    return sample_dict

def split_sample(header,dir,type,sample_dict):
    """Split watson and crick bam file into sample specific reads groups"""
    #to limit what is processed modify header
    for name in header[4:]:
        if name.endswith('_methylated'):
            sample_name = name.replace('_methylated','')
            if sample_name not in sample_dict['watson']:
                continue
        else:
            continue
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


def filter_pileup(header,input_dir,output_dir,libtype,cov_treshold,sample_dict):
    """generate output for methylation rate according to values in pileup"""
    #output 1 is methylation.filtered.bed
    #output 2 is methylation_ratio.filtered.igv
    output_dict = {}
    out_2_handle = open(os.path.join(output_dir,'methylation_ratio.filtered.igv'),'w')
    for name in header[4:]:
        if name.endswith('_methylated'):
            sample_name = name.replace('_methylated','')
            if sample_name not in sample_dict['watson']:
                continue
        else:
            continue
        for strand in ['watson','crick']:
            handle = open(os.path.join(output_dir,'%s.pileup'%(sample_name+'_'+strand)))
            for line in handle:
                split_line = line.rstrip('\n').split('\t')
                meth,nometh,ratio = get_ratio(split_line[4],strand,cov_treshold,libtype)
                contig,pos = split_line[:2]
                if ratio != '':
                    try:
                        output_dict[contig][int(pos)][sample_name] = {'ratio':ratio,'meth':meth,'nometh':nometh}
                    except KeyError:
                        if contig not in output_dict:
                            output_dict[contig] = {int(pos):{sample_name:{'ratio':ratio,'meth':meth,'nometh':nometh}}}
                        elif int(pos) not in output_dict[contig]:
                            output_dict[contig][int(pos)] = {sample_name:{'ratio':ratio,'meth':meth,'nometh':nometh}}
    for contig,positions in sorted(output_dict.items()):
        for pos in sorted(positions):
            out_2 = [contig,str(pos)]
            for name in header[4:]:
                if name.endswith('_methylated'):
                    sample_name = name.replace('_methylated','')
                else:
                    continue
                try:
                    out_2.append('%.3f'%positions[pos][sample_name]['ratio'])
                except KeyError:
                    out_2.append('.')
            out_2_handle.write('\t'.join(out_2)+'\n')
    return output_dict


def get_ratio(nts,strand,cov_treshold,libtype):
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
    if min(plus_count,min_count) > cov_treshold:
        #Here we will take the average of both estimates from forward and reverse reads
        nometh = plus_count + min_count - meth
        #if one of the ratios below 0.05, return 0 as ratio, this to avoid dubious methylation calls.
        if min(_min_ratio,plus_ratio) < 0.05:
            return 0,nometh+meth,0.0
        else:
            return meth,nometh,(_min_ratio + plus_ratio) / 2.0
    elif libtype == 'WGBS':
        #for WGBS libraries we do not want to be as stringent as for epiGBS
        #This exercise serves to compare becker et al with epiGBS.
        nometh = plus_count + min_count - meth
        try:
            ratio = meth/float(nometh+meth)
        except ZeroDivisionError:
            ratio = 0.0
        return meth,nometh,ratio
    else:
        return None,None,''

def fisher_exact(output_dict,out_dir):
    """Do fisher exact test for pairwise DMP's between loci"""
    two_by_two = {}
    for contig,positions in sorted(output_dict.items()):
        for pos in sorted(positions.keys()):
            for i,(ind1,counts1) in enumerate(sorted(positions[pos].items())):
                for j,(ind2,counts2) in enumerate(sorted(positions[pos].items())):
                    if ind1 == ind2 or j<i:
                        continue
                    table = [[counts1['meth'],counts1['nometh']],
                             [counts2['meth'],counts2['nometh']]]
                    total_1 = float(counts1['nometh'] + counts1['meth'])
                    total_2 = float(counts2['nometh'] + counts2['meth'])
                    ratio_1 = counts1['meth'] / total_1
                    ratio_2 = counts2['meth'] / total_2

                    oddsratio, pvalue = stats.fisher_exact(table)
                    try:
                        two_by_two['%s_%s'%(ind1,ind2)][contig].append([pos,pvalue,ratio_1,ratio_2,total_1,total_2])
                    except KeyError:
                        if '%s_%s'%(ind1,ind2) not in two_by_two:
                            two_by_two['%s_%s'%(ind1,ind2)] = {}
                        two_by_two['%s_%s'%(ind1,ind2)][contig] = [[pos,pvalue,ratio_1,ratio_2,total_1,total_2]]
    #now apply bonferroni correction to table
    for comparison in two_by_two.keys():
        handle = open(os.path.join(out_dir,'%s.csv'%comparison),'w')
        output = {}
        total_length = 0
        for contig,values in two_by_two[comparison].items():
            total_length += len(values)
        min_p = 0.05 / total_length
        handle.write('%s\t%s\n'%(comparison,min_p))
        for contig,values in two_by_two[comparison].items():
            for item in values:
            #     try:
            #         output[contig].append(item)
            #     except KeyError:
            #         output[contig] = [item]
            #         # two_by_two[comparison] = output
            # if contig in output:
            #     for item in output[contig]:
                handle.write('%s\t%s\t%s'%(comparison,contig,'\t'.join([str(i) for i in item])))
                if item[1] < min_p:
                    handle.write('\t1\n')
                else:
                    handle.write('\t0\n')
    return two_by_two


def try_to_add(entry,dict,keys,position):
    """"try to add position to dict and create key if it fails"""
    try:
        dict[keys[0]]['%s_%s'%(keys[1],entry)] += [position]
    except KeyError:
        if keys[0] not in dict:
            dict[keys[0]] = {'%s_%s'%(keys[1],entry):[position]}
        else:
            dict[keys[0]]['%s_%s'%(keys[1],entry)] = [position]
    return dict


def check_symmetric_meth(items_in):
    """Check whether methylation is present or absent for both cytosines in symmetric contexts"""
    #[position,p-value,meth_ratio_0,meth_ratio_1,total_0,total_1]
    items_out = []
    min_meth = 0.05
    if len(items_in)%2:
        print "wrong length for items_in!"
    for i in range(0,len(items_in)-1,2):
        watson,crick = items_in[i:i+2]
        obs = 0
        for pos in [2,3]:
            #if methylation is present it should be symmetric and be present in both watson and crick for both strands
            if max(watson[pos],crick[pos]) > min_meth:
                if min(watson[pos],crick[pos]) > min_meth:
                    obs += 1
            #if no evidence for methylation is found, good, retain positions as well
            elif max(watson[pos],crick[pos]) < min_meth:
                if min(watson[pos],crick[pos]) < min_meth:
                    obs += 1
        #only if presence / absence for both individuals is congruent add to output
        if obs == 2:
            items_out.append(watson)
            items_out.append(crick)
    return items_out


def retain_symmetric(dict,context_dic):
    """retain only symmetric positions in case of methylation both strand should be methylated!"""
    for comp,subdict in dict.items():
        for contig,items in subdict.items():
            #reset the content of subdict[contig] so that it is empty
            subdict[contig] = []
            context_items = {}
            in_items = []
            for item in items:
                try:
                    context = context_dic[contig][str(item[0])]
                except KeyError:
                    continue
                try:
                    context_items[context].append(item)
                except KeyError:
                    context_items[context] = [item]
            #now retain only symmetric positions in comparison for CG and CHG
            for context,items in sorted(context_items.items()):
                if context == 'CHH':
                    subdict[contig] += items
                    continue
                elif 'CHG' in context:
                    offset = 2
                elif 'CG' in context:
                    offset =1
                for i,item in enumerate(items[:-1]):
                    if items[i+1][0] == item[0] + offset:
                        if item not in in_items:
                            in_items.append(item)
                        else:
                            continue
                        if items[i+1] not in in_items:
                            in_items.append(items[i+1])
                #now check if methylated positions in CG and CHG are symmetric if methylated, otherwise remove
                out_items = check_symmetric_meth(in_items)
                subdict[contig] += out_items
                in_items = []
            dict[comp][contig] = subdict[contig]
    return dict

def retain_symmetric_DMP(summary_dict):
    """only retain symmetric DMPs"""
    out_dict = {}
    for cluster_type,subdict in summary_dict.items():
        for type,items in subdict.items():
            if type.split('_')[0] == 'CG':
                offset = 1
            elif type.split('_')[0] == 'CHG':
                offset = 2
            else:
                continue
            items_out = []
            items = sorted(items, key=itemgetter(0,1,2))
            for i,item in enumerate(items[:-1]):
                #items should be from same comparison on same contig:
                if items[i+1][0] == item[0] and items[i+1][1] == item[1]:
                    if items[i+1][2] == item[2] + offset:
                        if item not in items_out:
                            items_out.append(items[i])
                            items_out.append(items[i+1])
                        else:
                            continue
            try:
                out_dict[cluster_type][type] = items_out
            except KeyError:
                if cluster_type not in out_dict:
                    out_dict[cluster_type] = {type:items_out}
    return out_dict


def get_dmp_stats(combined_entries,output_dir,context):
    """Write the combined entries and generate tables"""
    positions = {}
    #make a dictionary with dict[cluster] = {'cg_count':12,'chg_count':30,'chh_count':20,'cg_DMP':}
    dmp_dict = {}
    for comparison,subdict in combined_entries.items():
        item_len = sum([len(v) for v in subdict.values()])
        p_min = 0.05 / float(item_len)
        for contig,list in subdict.items():
            for item in list:
                position,pvalue,ratio0,ratio1,total1,total2 = item
                ratio_diff = max(ratio0,ratio1) - min(ratio0,ratio1)
                total_sum = int(total1) + int(total2)
                try:
                    ctype = context[contig][str(position)]
                except KeyError:
                    continue
                #maybe we should have uniform representation here:
                output = (comparison,contig,position,pvalue,ratio0,ratio1,total1,total2)
                dmp_dict = try_to_add('count',dmp_dict,[contig,ctype],output)
                if pvalue < p_min:# and ratio_diff > 0.7:#
                    dmp_dict = try_to_add('DMP',dmp_dict,[contig,ctype],output)
    return dmp_dict


def write_dmp_stats(dmp_dict,output_dir,input_dir):
    """Get DMP statistics per species"""
    gene_list = open(os.path.join(input_dir,'output_denovo') + '/genes.txt','r')
    gene_list = gene_list.read().split('\n')
    te_list = open(os.path.join(input_dir,'output_denovo') + '/te.txt','r')
    te_list = te_list.read().split('\n')
    summary_dict = {}
    median_dict = {}
    diff_dict = {}
    out_path = os.path.join(output_dir,'dmp.csv')
    output = open(out_path,'w')
    for contig,subdict in dmp_dict.items():
        if contig in gene_list:
            cluster_type = 'gene'
        elif contig in te_list:
            cluster_type = 'te'
        else:
            cluster_type = 'other'
        for type,items in subdict.items():
            median_add = [sum(v[-2:]) for v in items]
            diff_add = ['%.4f'%(v[5]-v[4]) for v in items]
            if items == []:
                continue
            try:
                median_dict[cluster_type][type]+= median_add
                diff_dict[cluster_type][type]+= diff_add
            except KeyError:
                if cluster_type not in median_dict:
                    median_dict[cluster_type] = {type:median_add}
                    diff_dict[cluster_type ] = {type:diff_add}
                else:
                    median_dict[cluster_type][type] = median_add
                    diff_dict[cluster_type][type] = diff_add
            #now we go for writing DMP statistics..
            try:
                summary_dict[cluster_type][type]+=items
            except KeyError:
                if cluster_type not in summary_dict:
                    summary_dict[cluster_type] = {type:items}
                else:
                    summary_dict[cluster_type][type] = items
    # for cluster_type,subdict in summary_dict.items():
    #     for key,value in sorted(subdict.items()):
    #         median = sorted(median_dict[type][key])[value/2]
    #         print type,key,value,median
    symmetric_summary_dict = retain_symmetric_DMP(summary_dict)
    for cluster_type,subdict in symmetric_summary_dict.items():
        for key,values in sorted(subdict.items()):
            if 'CHH' not in key:
                for v in values:
                    if cluster_type != 'te':
                        out_line = [cluster_type,key] + [str(i) for i in v]
                        output.write('\t'.join(out_line)+'\n')
                if values == []:
                    continue
            else:
                continue
            out_line = [cluster_type] + key.split('_') + ['%s'%(len(values))]
            print '\t'.join(out_line)

def summary_methylation(dmp_dict,output_dir,input_dir):
    """Get methylation statistics per species"""
    gene_list = open(os.path.join(input_dir,'output_denovo') + '/genes.txt','r')
    gene_list = gene_list.read().split('\n')
    te_list = open(os.path.join(input_dir,'output_denovo') + '/te.txt','r')
    te_list = te_list.read().split('\n')
    summary_dict = {}
    for contig,subdict in dmp_dict.items():
        if contig in gene_list:
            type = 'gene'
        elif contig in te_list:
            type = 'te'
        else:
            type = 'other'
        for key,value in subdict.items():
            if value == [] or 'count' not in key:
                continue
            if 'CHH' not in key:
                methylated = 0
                methylated = len([v[-4] for i,v in enumerate(value[:-1]) if v[-4]>0.05 and value[i+1][-4]>0.05 and not i%2])
                methylated += len([v[-3] for i,v in enumerate(value[:-1]) if v[-3]>0.05 and value[i+1][-3]>0.05 and not i%2])
            else:
                methylated = len([v[-3] for v in value if v[-3]>0.05])
                methylated += len([v[-4] for v in value if v[-4]>0.05])
            non_meth = len(value) - methylated
            try:
                summary_dict[type][key]['methylated'] += methylated
                summary_dict[type][key]['total'] += methylated+non_meth
            except KeyError:
                if type not in summary_dict:
                    summary_dict[type] = {key:{'methylated':methylated,'total':methylated+non_meth}}
                elif key not in summary_dict[type]:
                    summary_dict[type][key] = {'methylated':methylated,'total':methylated+non_meth}
    for type,subdict in summary_dict.items():
        for key,value in sorted(subdict.items()):
            output = '\t'.join([type,key.split('_')[0],str(value['methylated']),str(value['total'])])
            print output


def run_epiGBS(dir):
    """Run analysis for epiGBS data in directory"""
    ref = os.path.join(dir,'output_denovo','consensus_cluster.renamed.fa')
    input_bed = os.path.join(dir,'output_mapping','methylation.bed')
    input_dir = os.path.join(dir,'output_mapping')
    out_dir = os.path.join(dir,'dmp')
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)
    with open(input_bed,'r') as file_handle:
        #get global header
        header = file_handle.readline()[:-1].split('\t')
    #set minimum coverage per strand
    mincov = 10
    #Get genome dictionary to determine strand of positions
    genome_dict = SeqIO.to_dict(SeqIO.parse(open(ref),'fasta'))
    #process bed file to get per sample list of positions that are candidates
    sample_dict = process_bed(input_bed,genome_dict,mincov)
    #make sample specific bed filed
    make_bed_files(sample_dict,out_dir)
    #split bam file in sample specific bam files
    split_sample(header,input_dir,'epiGBS',sample_dict)
    #make pileup per bam file
    make_pileup(header,out_dir,sample_dict)
    # filter pileup, calculate per position methylation rate
    libtype = 'epiGBS'
    output_dict = filter_pileup(header,input_dir,out_dir,libtype,mincov,sample_dict)
    #Do Fisher exact test for differentially methylated positions between loci
    diff_meth_epiGBS = fisher_exact(output_dict,out_dir)
    # combined_entries = compare_epi_WGBS()
    context = get_context(input_bed)
    #get counts per context
    #only retain symmetric positions
    diff_meth_epiGBS = retain_symmetric(diff_meth_epiGBS,context)
    dmp_stats = get_dmp_stats(diff_meth_epiGBS,out_dir,context)
    write_dmp_stats(dmp_stats,out_dir,dir)
    summary_methylation(dmp_stats,out_dir,dir)

run_epiGBS('/Users/thomasvangurp/epiGBS/Baseclear/unfiltered_sequences/seqNNAtlE/Fallopia/')