__author__ = 'thomasvangurp'
"""Filter methylation bed file with coverage criterion on both strands"""
from Bio import SeqIO
import os
import subprocess
from scipy import stats

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


def get_dmp_stats(combined_entries,output_dir,context):
    """Write the combined entries and generate tables"""
    positions = {}
    out_path = os.path.join(output_dir,'dmp.csv')
    output = open(out_path,'w')
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
                dmp_dict = try_to_add('count',dmp_dict,[contig,ctype],(position,total_sum))
                if pvalue < p_min and ratio_diff > 0.5:#
                    dmp_dict = try_to_add('DMP',dmp_dict,[contig,ctype],(position,total_sum))
    return dmp_dict

def process_bed(input_bed,genome_dict,mincov):
    """Process input methylation.bed and generate sample specific watson and crick bed files"""
    handle = open(input_bed,'r')
    header = handle.readline()[:-1].split('\t')
    sample_dict = {'watson':{},'crick':{}}
    for line in handle:
        split_line = line.rstrip('\n').split('\t')
        for i in range(5,len(split_line),2):
            total = split_line[i]
            contig = split_line[0]
            pos = int(split_line[1])
            name = header[i].replace('_total','')
            try:
                if int(total) > (mincov * 2):
                    try:
                        if genome_dict[contig].seq[pos-1] == 'C':
                            strand = 'watson'
                        elif genome_dict[contig].seq[pos-1] == 'G':
                            strand = 'crick'
                        else:
                            print "Nucleotide incorrect, %s"%genome_dict[contig].seq[pos-1]
                        sample_dict[strand][name][contig].append(pos)
                    except KeyError:
                        if name not in sample_dict[strand]:
                            sample_dict[strand][name] = {contig:[pos]}
                        elif contig not in sample_dict[strand][name]:
                            sample_dict[strand][name][contig] = [pos]
            except ValueError:
                pass
    return sample_dict

def make_pileup(header,dir,type):
    """make pileup"""
    for name in header[4:]:
        if name.endswith('_methylated'):
            sample_name = name.replace('_methylated','')
        else:
            continue
        for strand in ['watson','crick']:
            output_loc = {}
            output_loc['bamfile'] = '/tmp/%s_%s.bam'%(sample_name+'_'+strand,type)
            output_loc['output'] = os.path.join(dir,'filtering','%s_%s.pileup'%(sample_name+'_'+strand,type))
            output_loc['regions'] = os.path.join(dir,'filtering','%s_%s.bed'%(sample_name,strand))
            pileup_cmd = 'samtools mpileup -l %(regions)s %(bamfile)s > %(output)s'%output_loc
            #only run command if file does not yet exist
            file_check = os.path.exists(output_loc['output'])
            if file_check:
                continue
            p = subprocess.Popen(pileup_cmd,stdout=subprocess.PIPE,stderr=subprocess.PIPE,
                             shell=True,executable='/bin/bash')
            exit_code = p.wait()
    return 0

def filter_pileup(header,dir,type,cov_treshold):
    """generate output for methylation rate according to values in pileup"""
    #output 1 is methylation.filtered.bed
    #output 2 is methylation_ratio.filtered.igv
    output_dict = {}
    out_2_handle = open(os.path.join(dir,'filtering','methylation_ratio.filtered.igv'),'w')
    for name in header[4:]:
        if name.endswith('_methylated'):
            sample_name = name.replace('_methylated','')
        else:
            continue
        for strand in ['watson','crick']:
            handle = open(os.path.join(dir,'filtering','%s_%s.pileup'%(sample_name+'_'+strand,type)))
            for line in handle:
                split_line = line.rstrip('\n').split('\t')
                meth,nometh,ratio = get_ratio(split_line[4],strand,cov_treshold,type)
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

def run_epiGBS():
    """Run analysis for epiGBS data"""
    ref = '/Users/thomasvangurp/epiGBS/Baseclear/unfiltered_sequences/seqNNAtlE/Athal/TAIR10_ref/Arabidopsis_thaliana.TAIR10.22.dna.genome.fa'
    input_bed = '/Users/thomasvangurp/epiGBS/Baseclear/unfiltered_sequences/seqNNAtlE/Athal/ref_mapping/methylation.bed'
    input_dir = '/Users/thomasvangurp/epiGBS/Baseclear/unfiltered_sequences/seqNNAtlE/Athal/ref_mapping/'
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
    out_dir = '/Users/thomasvangurp/epiGBS/Baseclear/unfiltered_sequences/seqNNAtlE/Athal/ref_mapping/filtering'
    make_bed_files(sample_dict,out_dir)
    #split bam file in sample specific bam files
    split_sample(header,input_dir,'epiGBS')
    #make pileup per bam file
    make_pileup(header,input_dir,'epiGBS')
    #filter pileup, calculate per position methylation rate
    output_dict = filter_pileup(header,input_dir,'epiGBS',mincov)
    #Do Fisher exact test for differentially methylated positions between loci
    diff_meth_epiGBS = fisher_exact(output_dict,out_dir)
    context = get_context(input_bed)
    dmp_stats = get_dmp_stats(diff_meth_epiGBS,out_dir,context)
    # dmp_stats = retain_symmetric(dmp_stats)


def retain_symmetric(entries):
    """retain only symmetric positions"""
    for type in entries.keys():
        if type == 'neither':
            continue
        count_dict = {}
        for i in entries[type]:
            contig,position,context = i[:3]
            try:
                count_dict[context][contig][int(position)] = i
            except KeyError:
                if context not in count_dict:
                    count_dict[context] = {contig:{int(position):i}}
                elif contig not in count_dict[context]:
                    count_dict[context][contig] = {int(position):i}
        output = []
        for (context,offset) in [('CG',1),('CHG',2)]:
            if context not in count_dict:
                continue
            for contig in count_dict[context]:
                pos_out = []
                pos_sorted = sorted(count_dict[context][contig].keys())
                for i,pos in enumerate(pos_sorted[:-1]):
                    if pos_sorted[i+1] == pos + offset:
                        if pos not in pos_out:
                            pos_out.append(str(pos))
                        if pos_sorted[i+1] not in pos_out:
                            pos_out.append(str(pos_sorted[i+1]))
                output += [count_dict[context][contig][int(i)] for i in pos_out]
        entries[type] =  output
    return entries

def run_WGBS():
    """run  analysis WGBS data"""
    ref = '/Users/thomasvangurp/epiGBS/Baseclear/unfiltered_sequences/seqNNAtlE/Athal/TAIR10_ref/Arabidopsis_thaliana.TAIR10.22.dna.genome.fa'
    input_dir = '/Users/thomasvangurp/epiGBS/becker_nature_2011/analysis/new_vcf'
    input_bed = '/Users/thomasvangurp/epiGBS/becker_nature_2011/analysis/new_vcf/methylation.bed'
    with open(input_bed,'r') as file_handle:
        #get global header
        header = file_handle.readline()[:-1].split('\t')
    #set minimum coverage per strand
    mincov = 5
    #Get genome dictionary to determine strand of positions
    genome_dict = SeqIO.to_dict(SeqIO.parse(open(ref),'fasta'))
    #process bed file to get per sample list of positions that are candidates
    sample_dict = process_bed(input_bed,genome_dict,mincov)
    #make sample specific bed filed
    out_dir = '/Users/thomasvangurp/epiGBS/becker_nature_2011/analysis/new_vcf/filtering'
    make_bed_files(sample_dict,out_dir)
    #split bam file in sample specific bam files
    split_sample(header,input_dir,'WGBS')
    #make pileup per bam file
    make_pileup(header,input_dir,'WGBS')
    #filter pileup, calculate per position methylation rate
    output_dict = filter_pileup(header,input_dir,'WGBS',mincov)
    #Do Fisher exact test for differentially methylated positions between loci
    diff_meth_WGBS = fisher_exact(output_dict,out_dir)
    context = get_context(input_bed)
    dmp_stats = get_dmp_stats(diff_meth_WGBS,out_dir,context)
    # dmp_stats = retain_symmetric(dmp_stats)

def compare_epi_WGBS():
    """Compare the results of WGBS and epiGBS methylation rates"""
    #XY scatter plot of p-values?
    #context
    #only retain positions called in both studies
    dir_WGBS = '/Users/thomasvangurp/epiGBS/becker_nature_2011/analysis/new_vcf/filtering/'
    dir_epiGBS = '/Users/thomasvangurp/epiGBS/Baseclear/unfiltered_sequences/seqNNAtlE/Athal/ref_mapping/filtering/'
    combined_entries = {}
    for content in os.walk(dir_epiGBS):
        files = content[-1]
        for file in files:
            if file.endswith('.csv'):
                with open(dir_epiGBS+file,'r') as epiGBS_handle:
                    header = epiGBS_handle.readline()
                    epiGBS_values = {}
                    for line in epiGBS_handle:
                        split_line = line.rstrip('\n').split('\t')
                        comp,contig,pos,pvalue,ratio_1,ratio_2,total_1,total_2,incl = split_line
                        try:
                            epiGBS_values[contig][pos] = (float(pvalue),ratio_1,ratio_2)
                        except KeyError:
                            if contig not in epiGBS_values:
                                epiGBS_values[contig] = {pos:(float(pvalue),ratio_1,ratio_2)}

                try:
                    file_WGBS = open(dir_WGBS+file.replace('Athal_A','30-'),'r')
                except KeyError:
                    continue
                header = file_WGBS.readline()
                WGBS_values = {}
                for line in file_WGBS:
                    split_line = line.rstrip('\n').split('\t')
                    comp,contig,pos,pvalue,ratio_1,ratio_2,total_1,total_2,incl = split_line
                    try:
                        WGBS_values[contig][pos] = (float(pvalue),ratio_1,ratio_2)
                    except KeyError:
                        if contig not in WGBS_values:
                            WGBS_values[contig] = {pos:(float(pvalue),ratio_1,ratio_2)}
                #combined entries
                for k,v in epiGBS_values.items():
                    for pos in v:
                        if pos in WGBS_values[k]:
                            try:
                                combined_entries[file][k][pos] = {'epiGBS':epiGBS_values[k][pos],'WGBS':WGBS_values[k][pos]}
                            except KeyError:
                                if file not in combined_entries:
                                    combined_entries[file] = {k:{pos:{'epiGBS':epiGBS_values[k][pos],'WGBS':WGBS_values[k][pos]}}}
                                elif k not in combined_entries[file]:
                                    combined_entries[file][k] = {pos:{'epiGBS':epiGBS_values[k][pos],'WGBS':WGBS_values[k][pos]}}
    return combined_entries

def compare_generations():
    """Compare the results of WGBS between generations"""
    #XY scatter plot of p-values?
    #context
    #only retain positions called in both studies
    dir_WGBS = '/Users/thomasvangurp/epiGBS/becker_nature_2011/analysis/new_vcf/filtering/'
    combined_entries = {}
    for file in ['30-109_30-39.csv','30-109_31-39.csv','30-109_30-49.csv','30-109_31-49.csv']:
        if file in ['30-109_30-39.csv','30-109_30-49.csv']:
            gen30_values = {}
            with open(dir_WGBS+file,'r') as gen30_handle:
                header = gen30_handle.readline()
                for line in gen30_handle:
                    split_line = line.rstrip('\n').split('\t')
                    comp,contig,pos,pvalue,ratio_1,ratio_2,total_1,total_2,incl = split_line
                    try:
                        gen30_values[contig][pos] = (float(pvalue),ratio_1,ratio_2)
                    except KeyError:
                        if contig not in gen30_values:
                            gen30_values[contig] = {pos:(float(pvalue),ratio_1,ratio_2)}
        if file in ['30-109_31-39.csv','30-109_31-49.csv']:
            gen31_values = {}
            with open(dir_WGBS+file,'r') as gen31_handle:
                header = gen31_handle.readline()
                for line in gen31_handle:
                    split_line = line.rstrip('\n').split('\t')
                    comp,contig,pos,pvalue,ratio_1,ratio_2,total_1,total_2,incl = split_line
                    try:
                        gen31_values[contig][pos] = (float(pvalue),ratio_1,ratio_2)
                    except KeyError:
                        if contig not in gen31_values:
                            gen31_values[contig] = {pos:(float(pvalue),ratio_1,ratio_2)}
            for k,v in gen31_values.items():
                for pos in v:
                    if pos in gen30_values[k]:
                        try:
                            combined_entries[file][k][pos] = {'gen30':gen30_values[k][pos],'gen31':gen31_values[k][pos]}
                        except KeyError:
                            if file not in combined_entries:
                                combined_entries[file] = {k:{pos:{'gen30':gen30_values[k][pos],'gen31':gen31_values[k][pos]}}}
                            elif k not in combined_entries[file]:
                                combined_entries[file][k] = {pos:{'gen30':gen30_values[k][pos],'gen31':gen31_values[k][pos]}}
    return combined_entries


def write_combined_gen(combined_entries,context):
    """Write the combined entries and generate tables"""
    output = open('/Users/thomasvangurp/epiGBS/becker_nature_2011/analysis/new_vcf/filtering/comparison_DMP_39_49.csv','w')
    positions = {}
    for comparison,subdict in combined_entries.items():
        item_len = sum([len(v) for v in subdict.values()])
        p_min = 0.05 / float(item_len)
        count = {'gen30':[],'gen31':[],'neither':[],'both':[]}
        for contig,list in subdict.items():
            for key,value in list.items():
                gen30dif = float(value['gen30'][1]) - float(value['gen30'][2])
                gen31dif = float(value['gen31'][1]) - float(value['gen31'][2])
                gen30p = float(value['gen30'][0])
                gen31p = float(value['gen31'][0])
                # if abs(gen30dif) < 0.7:
                #     if abs(gen31dif) < 0.7:
                #         continue
                ctype = context[contig][key]
                if gen31p < p_min:
                    if gen30p < p_min:
                        count['both'].append((contig,key,ctype,gen30dif,gen31dif,gen30p,gen31p,))
                    else:
                        count['gen31'].append((contig,key,ctype,gen30dif,gen31dif,gen30p,gen31p,))
                elif gen30p < p_min:
                    count['gen30'].append((contig,key,ctype,gen30dif,gen31dif,gen30p,gen31p,))
                else:
                    count['neither'].append((gen30dif,gen31dif,gen30p,gen31p,))
        #set function to only retain symmetric with min_count
        # count = retain_symmetric(count)
        for k,v in sorted(count.items()):
            if k!='neither':
                for ratios in v:
                    output.write('%s\t%s\t'%(comparison,k) + '\t'.join([str(I).replace('.',',') for I in ratios]) + '\n')
            print comparison,k,len(v)
            if '39' in comparison:
                title = 'line 39'
            else:
                title = 'line 49'
        make_venn((len(count['gen30']),len(count['gen31']),len(count['both'])),title,p_min)
        get_coverage_WGBS(combined_entries)
def get_context(file_in):
    """Get genomic context from bed file"""
    output_dict = {}
    handle = open(file_in,'r')
    header = handle.readline().rstrip('\n').split('\t')
    for line in handle:
        split_line = line.rstrip('\n').split('\t')
        contig,pos,context = split_line[:3]
        try:
            output_dict[contig][pos] = context
        except KeyError:
            output_dict[contig] = {pos:context}
    return output_dict



def write_combined_entries(combined_entries,context):
    """Write the combined entries and generate tables"""
    output = open('/Users/thomasvangurp/epiGBS/becker_nature_2011/analysis/new_vcf/filtering/comparison_DMP.csv','w')
    positions = {}
    for comparison,subdict in combined_entries.items():
        item_len = sum([len(v) for v in subdict.values()])
        p_min = 0.05 / float(item_len)
        count = {'WGBS':[],'epiGBS':[],'neither':[],'both':[]}
        for contig,list in subdict.items():
            for key,value in list.items():
                epidif = float(value['epiGBS'][1]) - float(value['epiGBS'][2])
                WGBSdif = float(value['WGBS'][1]) - float(value['WGBS'][2])
                ctype = context[contig][key]
                if abs(epidif) < 0.7:
                    if abs(WGBSdif) < 0.7:
                        continue
                if value['WGBS'][0] < p_min:
                    if value['epiGBS'][0] < p_min:
                        count['both'].append((contig,key,ctype,epidif,WGBSdif,))
                    else:
                        count['WGBS'].append((contig,key,ctype,epidif,WGBSdif,))
                elif value['epiGBS'][0] < p_min:
                    count['epiGBS'].append((contig,key,ctype,epidif,WGBSdif,))
                else:
                    count['neither'].append((epidif,WGBSdif,))
        #only retain symmetric counts
        # count = retain_symmetric(count)
        for k,v in sorted(count.items()):
            if k!='neither':
                for ratios in v:
                    output.write('%s\t%s\t'%(comparison,k) + '\t'.join([str(I) for I in ratios]) + '\n')
            print comparison,k,len(v)
        # get_coverage(combined_entries)
        make_venn((len(count['WGBS']),len(count['epiGBS']),len(count['both'])),comparison,p_min)

def get_coverage(combined_entries):
    """return average coveraged for positions covered"""
    cov_WGBS = []
    cov_epiGBS = []
    bed_WGBS = open('/Users/thomasvangurp/epiGBS/becker_nature_2011/analysis/new_vcf/methylation.bed','r')
    bed_epiGBS = open('/Users/thomasvangurp/epiGBS/Baseclear/unfiltered_sequences/seqNNAtlE/Athal/ref_mapping/methylation.bed','r')
    positions = {}
    for comparison,subdict in combined_entries.items():
         for contig,list in subdict.items():
            for key,value in list.items():
                try:
                    if int(key) not in positions[contig]:
                        positions[contig].update([int(key)])
                except KeyError:
                    positions[contig] = set([int(key)])
    header = bed_WGBS.readline().rstrip('\n').split('\t')
    for line in bed_WGBS:
        split_line = line.rstrip('\n').split('\t')
        values = []
        contig,pos = split_line[:2]
        if int(pos) not in positions[contig] or split_line[2]=='.':
            continue
        for value,name in zip(split_line,header):
            if name in ['30-109_total','30-119_total','30-29_total','30-89_total']:
                try:
                    values.append(int(value))
                except ValueError:
                    pass
        if values != []:
            cov_WGBS.append(sum(values)/float(len(values)))
    avg_cov_WGBS = sum(cov_WGBS)/float(len(cov_WGBS))
    header = bed_epiGBS.readline().rstrip('\n').split('\t')
    for line in bed_epiGBS:
        split_line = line.rstrip('\n').split('\t')
        values = []
        contig,pos = split_line[:2]
        if int(pos) not in positions[contig] or split_line[2]=='.':
            continue
        for value,name in zip(split_line,header):
            if 'total' in name:
                try:
                    values.append(int(value))
                except ValueError:
                    pass
        if values != []:
            cov_epiGBS.append(sum(values)/float(len(values)))
    avg_cov_epiGBS = sum(cov_epiGBS)/float(len(cov_epiGBS))
    print 'average coverage WGBS over %s\t%s\n'%(len(cov_WGBS),avg_cov_WGBS)
    print 'average coverage epiGBS over %s\t%s\n'%(len(cov_epiGBS),avg_cov_epiGBS)


def get_coverage_WGBS(combined_entries):
    """return average coveraged for positions covered"""
    cov_30 = []
    cov_31 = []
    bed_WGBS = open('/Users/thomasvangurp/epiGBS/becker_nature_2011/analysis/new_vcf/methylation.bed','r')
    positions = {}
    for comparison,subdict in combined_entries.items():
         for contig,list in subdict.items():
            for key,value in list.items():
                try:
                    if int(key) not in positions[contig]:
                        positions[contig].update([int(key)])
                except KeyError:
                    positions[contig] = set([int(key)])
    header = bed_WGBS.readline().rstrip('\n').split('\t')
    for line in bed_WGBS:
        split_line = line.rstrip('\n').split('\t')
        values_30 = []
        values_31 = []
        contig,pos = split_line[:2]
        if int(pos) not in positions[contig] or split_line[2]=='.':
            continue
        for value,name in zip(split_line,header):
            if name in ['30-39_total','30-49_total']:
                try:
                    values_30.append(int(value))
                except ValueError:
                    pass
            if name in ['31-39_total','31-49_total']:
                try:
                    values_31.append(int(value))
                except ValueError:
                    pass
        if values_30 != []:
            cov_30.append(sum(values_30)/float(len(values_30)))
        if values_31 != []:
            cov_31.append(sum(values_31)/float(len(values_31)))
    avg_cov_30 = sum(cov_30)/float(len(cov_30))
    avg_cov_31 = sum(cov_31)/float(len(cov_31))
    print 'average coverage gen 30 over %s\t%s\n'%(len(cov_30),avg_cov_30)
    print 'average coverage gen 31 over %s\t%s\n'%(len(cov_31),avg_cov_31)


def make_venn(sizes,title,minp):
    """Write venn diagram"""
    from matplotlib import pyplot as plt
    from matplotlib_venn import venn2
    try:
        elements = title.split('.')[0].split('_')
        title = 'DMPs with p < %.2E between %s and %s'%(minp,elements[1],elements[3])
        c = venn2(subsets = sizes, set_labels = ('WGBS','epiGBS',''))
    except IndexError:
        elements = title.split('.')[0].split('_')
        title = 'DMPs with p < %.2E in %s'%(minp,title)
        c = venn2(subsets = sizes, set_labels = ('WGBS gen30', 'WGBS gen31'))
    # c.set_lw(1.0)
    # c.set_ls('dotted')
    plt.title(title)
    # # plt.annotate('Unknown set', xy=v.get_label_by_id('100').get_position() - np.array([0, 0.05]), xytext=(-70,-70),
    #              ha='center', textcoords='offset points', bbox=dict(boxstyle='round,pad=0.5', fc='gray', alpha=0.1),
    #              arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=0.5',color='gray'))
    plt.show()
    # plt.savefig('/Users/thomasvangurp/Dropbox/Thomas_NIOO/epiGBS paper/Rebuttal/figures_additional_DMP/%s.pdf'%title,dpi=gcf().dpi)


def main():
    # run_epiGBS()
    # run_WGBS()
    # make_venn()
    # combined_entries = compare_epi_WGBS()
    context = get_context('/Users/thomasvangurp/epiGBS/becker_nature_2011/analysis/new_vcf/methylation.bed')
    entries_generations = compare_generations()
    # write_combined_entries(combined_entries,context)
    write_combined_gen(entries_generations,context)
if __name__ == '__main__':
    main()