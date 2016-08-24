#!/usr/bin/env pypy
import argparse
import os
from Bio import SeqIO,Restriction
import re
import gzip

def parse_args():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(description='Process input files')
    parser.add_argument('-m', '--mergedcalls', type=str, default=None,
                        help='Merged watson and crick calls')
    parser.add_argument('-s', '--SNP_input', type=str, default=None,
                        help='SNP input file, disable')
    parser.add_argument('-r', '--reference', type=str, default=None,
                        help='reference genome')
    parser.add_argument('-b', '--barcodes', type=str, default=None,
                        help='barcodes and enzymes used')
    parser.add_argument('-o', '--methylation_output', type=str, nargs='?', default=None,
                        help='methylation.bed output')
    parser.add_argument('-heat', '--heatmap_output', type=str, nargs='?', default=None,
                        help='Heatmap igv file output name')
    args = parser.parse_args()
    return args

def make_header(args,handle,split_line):
    """create header for bed and IGV file"""
    header = ['chr', 'pos', 'context', 'samples_called']
    for element in split_line[9:]:
        header.append('%s_methylated' % element)
        header.append('%s_total' % element)
    output =  '\t'.join(header) + '\n'
    return output

def check_contect_SNP(ref_base, chrom, pos, context, line, SNP_nearby):
    """return dictionary with lines changes per context"""
    context_dict = {context:line}
    if context == 'CG':
        #TODO: methylation polymorphis changed induced by SNP leading to absence of G should result in CHG or CHH!
        return context_dict
    if ref_base == 'C':
        if context == 'CHH':
            to_check = {pos + 1: 'CG', pos + 2: 'CHG'}
        elif context == 'CHG':
            to_check = {pos + 1: 'CG'}
    else:
        if context == 'CHH':
            to_check = {pos - 1: 'CG', pos - 2: 'CHG'}
        elif context == 'CHG':
            to_check = {pos - 1: 'CG'}
    for SNP in SNP_nearby:
        if SNP[0] == chrom:
            if int(SNP[1]) > pos + 2:
                break
            if int(SNP[1]) in to_check:
                new_context = to_check[int(SNP[1])]
                if (ref_base == 'C' and 'G' in SNP[4]) or (ref_base == 'G' and 'C' in SNP[4]):
                    for n,call in enumerate(SNP[9:]):
                        if call.split(':')[0] not in ['./.','0/0']:
                            if new_context == 'CHG' and 'CG' in context_dict:
                                #A C in CGG context is CG, not CHG. if CG was called previously skip this call.
                                if context_dict['CG'][n + 5] != '':
                                    continue
                            try:
                                context_dict[new_context][n + 5] = call
                                context_dict[context][n + 5] = ':'.join(['0,0']*4)
                            except KeyError:
                                context_dict[new_context] = context_dict[context][:5] + [':'.join(['0,0']*4)] * (len(line) - 5 )
                                context_dict[new_context][n + 5] = context_dict[context][n + 5]
                                context_dict[context][n + 5] = ':'.join(['0,0']*4)

        else:
            break
    return context_dict

def calc_context(split_line, genome, SNP_nearby):
    """
    If there are no SNP's neighbouring the methylation call, the context can be called
    by using the reference genome.
    """
    chrom, pos, ref_base = split_line[:3]
    pos = int(pos)
    slice_start = max(int(pos) - 3, 0)  # negative positions excluded
    slice_end = int(pos) + 2
    reference_bases = genome[chrom].seq[slice_start:slice_end]
    if ref_base == 'G':
        ref_context = reference_bases[0:2][::-1]
        if re.match('C.', str(ref_context)):
            context = 'CG'
        elif re.match('[ATG]C', str(ref_context)):
            context = 'CHG'
        else:
            context = 'CHH'
    elif ref_base == 'C':
        ref_context = reference_bases[3:5]
        if re.match('G.', str(ref_context)):
            context = 'CG'
        elif re.match('[ATC]G', str(ref_context)):
            context = 'CHG'
        else:
            context = 'CHH'
    else:
        context = '.'
    #check if context is affected by neighbouring SNPs
    context = check_contect_SNP(ref_base, chrom, pos, context, split_line, SNP_nearby)
    return context

def get_range(args):
    """Get range for which methylation polymorphisms can be called given enzyme overhang"""
    #parse barcodes for enzymes being used
    with open(args.barcodes,'r') as barcode_handle:
        header = barcode_handle.readline().rstrip('\n').split('\t')
        split_line =  barcode_handle.readline().rstrip('\n').split('\t')
        enzyme_left = split_line[header.index('ENZ_R1')]
        enzyme_right = split_line[header.index('ENZ_R2')]
        for enzyme in Restriction.AllEnzymes:
            if "%s"%(enzyme) == enzyme_left:
                left_start = len(enzyme.ovhgseq)
            elif "%s"%(enzyme) == enzyme_right:
                right_end = -1 *len(enzyme.ovhgseq)
    return left_start,right_end

def is_SNP(chrom, pos, SNP_nearby):
    """determine if variant is also in SNPs"""
    # Determine if there are any potential nearby SNPs
    for potential_SNP in SNP_nearby:
        if potential_SNP[:2] == [chrom, pos]:
            SNP = potential_SNP
            break
        else:
            SNP = None
    try:
        return SNP
    except UnboundLocalError:
        return None

def methylation_calling(split_line, context, SNP_nearby):
    """
    Main base calling algorithm. Determines methylation/SNP status for each sample having a watson and crick record.
    """
    chrom, pos, ref_base, watson_ALT, crick_ALT = split_line[:5]
    if ref_base.upper() in 'ATN':
        return None
    if len(ref_base) > 1:
        return None
    if ref_base.upper() == "G":
        out_line = [chrom, pos, context, '']
        #Note that nucleotides order of A and G are inverted (e.g. G/A) as obs is inverted
        #obs.split(',')[1] is for the Crick position
        [[out_line.append(int(obs.split(',')[0])) for nt, obs in zip('TGCA', obs[::-1].split(':'))
          if nt in 'GA'] for obs in split_line[5:]]
        for n in range(5, len(out_line), 2):
            out_line[n] += out_line[n - 1]
            if out_line[n] == 0:
                out_line[n - 1] = 'None'
                out_line[n] = 'None'
        out_line[3] = (len(out_line[4:]) - out_line.count('None')) / 2
        return out_line
    elif ref_base.upper() == "C":
        out_line = [chrom, pos, context,'']
        [[out_line.append(int(obs.split(',')[0])) for nt,obs in zip('ACGT',obs.split(':'))
          if nt in 'CT'] for obs in split_line[5:]]
        for n in range(5,len(out_line),2):
            out_line[n] += out_line[n-1]
            if out_line[n] == 0:
                out_line[n - 1] = 'None'
                out_line[n] = 'None'

        out_line[3] = (len(out_line[4:]) - out_line.count('None')) / 2
    return out_line

    return None
    # print ''
    # # If the sample is methylated, the processed_samples will be filled under the methylated key, or when it's
    # # Determine the reference base at the position of the VCF call.
    # ref_base = self.watson_record.REF
    # # loop watson and crick record for combined samples.
    # if min([self.watson_record.INFO['DP'], self.crick_record.INFO['DP']]) == 0:
    #     return None
    # # if self.watson_record.REF not in ['C', 'G']:
    # #     watson_alt = sum(self.watson_record.INFO['AD'][1:])/float(self.watson_record.INFO['DP'])
    # #     crick_alt = sum(self.crick_record.INFO['AD'][1:])/float(self.crick_record.INFO['DP'])
    # #     if max(watson_alt, crick_alt) < 0.05:
    # #         return None
    # self.call_genotypes()
    # for watson_sample, crick_sample in izip(self.watson_record, self.crick_record):
    #     # If there are is no call for both the watson and crick record sample, continue as we can not determine
    #     # whether polymorphism is a SNP/methylation polymorphism.
    #     if not watson_sample.called or not crick_sample.called:
    #         continue
    #     # Assigning the right alt base to the records.
    #     alt_watson = watson_sample.gt_bases.split('/')[1]
    #     alt_crick = crick_sample.gt_bases.split('/')[1]
    #
    #     sample_name = watson_sample.sample
    #     # TODO: move SNP calling to separate algorithm.
    #     if ref_base == 'C':
    #         if alt_crick == 'C' and alt_watson in 'CT':
    #             # Methylation in watson C/T No polymorphism in Crick: methylation
    #             self.processed_samples[sample_name]['methylated'] = watson_sample
    #             # try:
    #             #     #'AC' is not present in homozygous situations.
    #             #     if sum(crick_sample.site.INFO['AD'][1:]) > 0:
    #             #         #The Alternate alleles need to be called in at least one sample to be valid!
    #             #         # self.processed_samples[sample_name]['snp'] = crick_sample
    #             # except KeyError:
    #             #     pass
    #         elif alt_crick == 'A':
    #             # alt_crick contains another base, this SNP is valid if present in watson.
    #             if alt_watson in crick_sample.site.ALT:
    #                 # Both watson and crick contain the same alternate allele, the SNP is real?
    #                 alt_index = crick_sample.site.ALT.index(alt_watson)
    #                 try:
    #                     crick_alt_pct = crick_sample.data.AO[alt_index] / float(crick_sample.data.DP)
    #                 except TypeError:
    #                     crick_alt_pct = crick_sample.data.AO / float(crick_sample.data.DP)
    #                 alt_index = watson_sample.site.ALT.index(alt_watson)
    #                 try:
    #                     watson_alt_pct = watson_sample.data.AO[alt_index] / float(watson_sample.data.DP)
    #                 except TypeError:
    #                     watson_alt_pct = watson_sample.data.AO / float(watson_sample.data.DP)
    #                 if watson_alt_pct != 0.0 and crick_alt_pct != 0.0:
    #                     if max(crick_alt_pct, watson_alt_pct) / min(crick_alt_pct, watson_alt_pct) < 1.5:
    #                         continue
    #                         # self.processed_samples[sample_name]['snp'] = crick_sample
    #                         # TODO: merge alt counts for watson and crick here
    #                 elif crick_alt_pct == 0.0:
    #                     # REF:C watson C/T/A called as C/T crick C/A
    #                     # We can call both SNP and methylation. SNP from crick reliable
    #                     # Watson information on C/T ratio informative for methylation call
    #                     self.processed_samples[sample_name]['methylated'] = watson_sample
    #                 else:
    #                     pass
    #                     # Can this occur? TODO: check if this can be true
    #
    #         elif alt_watson == 'G' and alt_crick in 'AG':
    #             # C/G polymorphism in watson, C/G or C/A in Crick
    #             # Methylation information from Crick
    #             self.processed_samples[sample_name]['methylated'] = crick_sample
    #         elif alt_watson == 'T':
    #             if alt_crick == 'T':
    #                 # C/T variant in both watson and crick: SNP ==> only information from crick is reliable
    #                 # self.processed_samples[sample_name]['snp'] = crick_sample
    #                 if set(watson_sample.gt_bases.replace('/', '')) == set(['T']):
    #                     self.processed_samples[sample_name]['methylated'] = watson_sample
    #             if alt_crick != 'T' and alt_crick != 'C':
    #                 # Watson contains C/T methylation polymorphism and potentially other information
    #                 # Step 1: We can call Methylation polymorphism, crick is not C/T!
    #                 self.processed_samples[sample_name]['methylated'] = watson_sample
    #
    #     elif ref_base == 'G':
    #         # Watson is homozygous reference (i.e. no SNP) and crick has Methylation variation
    #         if alt_watson == 'G' and alt_crick in 'GA':
    #             self.processed_samples[sample_name]['methylated'] = crick_sample
    #         elif alt_crick == 'A' and alt_watson == 'A':
    #             # The crick allele can only be queried for methylation variation if it is fully converted
    #             # this means that no G can be in the Genotype.
    #             if set(crick_sample.gt_bases.replace('/', '')) == set(['A']):
    #                 self.processed_samples[sample_name]['methylated'] = crick_sample
    #         elif alt_watson == 'C' and alt_crick in 'CT':
    #             self.processed_samples[sample_name]['methylated'] = watson_sample
    # return 1

def get_SNP(handle, split_line, SNP_nearby):
    """return SNPs 2 basepairs up and downstream from methylation polymorphism"""
    chrom, pos = split_line[:2]
    for n,SNP in enumerate(SNP_nearby):
        chrom_SNP,pos_SNP = SNP[:2]
        if int(chrom_SNP) < int(chrom):
            SNP_nearby.pop(n)
        elif int(pos_SNP) < int(pos) -2 and chrom_SNP == chrom:
            SNP_nearby.pop(n)
        else:
            break
    while True:
        if len(SNP_nearby) < 5:
            next_SNP = handle.readline().rstrip('\n').split('\t')
            if next_SNP[0] == '':
                break
            SNP_nearby.append(next_SNP)
        else:
            break
    return SNP_nearby

def remove_SNP(split_line, SNP_nearby):
    """Eliminate Allele calls caused by a SNP from methylation call"""
    chrom, pos = split_line[:2]
    SNP = is_SNP(chrom, pos, SNP_nearby)
    if not SNP:
        return split_line
    else:
        out_line = split_line[:5]
        for SNP_call,observations in zip(SNP[9:],split_line[5:]):
            if SNP_call.split(':')[0] == '0/0':
                #TODO: revisit stringency of this rule. Should depend on which alleles are present in SNP per ind
                out_line.append(observations)
            else:
                out_line.append(':'.join(['0,0']*4))
        return out_line
        #remove SNP variation from methylation variant position.

def make_IGV_output(meth_record):
    """return IGV-style output of quantitative methylation value output"""
    out_record = [meth_record[0],str(int(meth_record[1])-1),meth_record[1], meth_record[2]]
    for i in range(4,len(meth_record),2):
        try:
            meth_float = '%.3f'% (meth_record[i] / float(meth_record[i+1]))
        except ZeroDivisionError:
            meth_float = '0.00'
        except ValueError:
            meth_float = '.'
        out_record.append(meth_float)
    return '\t'.join(out_record) + '\n'

def get_IGV_header(samples):
    """get header for IGV output file"""
    header = '#type=DNA_METHYLATION\nChromosome\tStart\tEnd\tFeature'
    for sample in samples[9:]:
        header += '\t%s' % sample
    return header + '\n'

def main():
    """main function loop"""
    args = parse_args()
    # make_header(args,handle,line)
    count = 0
    # enz_range = get_range(args)
    reference_genome = SeqIO.to_dict(SeqIO.parse(args.reference, "fasta"))
    SNP_nearby = []
    igv_handle = open(args.heatmap_output,'w')
    with open(args.methylation_output, 'w') as handle:
        with os.popen("pigz -cd %s" % args.SNP_input) as SNP_handle:
            #read all lines with comments
            while True:
                SNP = SNP_handle.readline().rstrip('\n').split('\t')
                if not SNP[0].startswith('#'):
                    SNP_nearby.append(SNP)
                    break
            # for line in os.popen("pigz -cd %s " % args.mergedcalls):
            for line in gzip.open(args.mergedcalls,'r'):
                split_line = line.rstrip('\n').split('\t')
                if not count and line.startswith('#'):
                    header = make_header(args,handle,split_line)
                    igv_header = get_IGV_header(split_line)
                    igv_handle.write(igv_header)
                    handle.write(header)
                    continue
                if split_line[2] not in ['C','G']:
                    #skip ref position with no CG
                    continue
                SNP_nearby = get_SNP(SNP_handle, split_line, SNP_nearby)
                #remove methylation variation observations distorted by SNP observations
                split_line = remove_SNP(split_line, SNP_nearby)
                context_dict = calc_context(split_line, reference_genome, SNP_nearby)
                count += 1
                if not count % 1000000:
                    print 'processed %s lines ' % count
                #if multiple contexts results, give key-value pairs for context/line
                for context, split_line in sorted(context_dict.items()):
                    #TODO: change field code for SNP induced different meth contexts
                    meth_record = methylation_calling(split_line, context, SNP_nearby)
                    if meth_record:
                        igv_record = make_IGV_output(meth_record)
                        igv_handle.write(igv_record)
                    if meth_record:
                        handle.write('\t'.join([str(s) for s in meth_record]) + '\n')
    handle.close()
    igv_handle.close()
    # os.popen('bgzip -f %s' % args.SNP_output)
    # os.popen('tabix -p vcf %s.gz' % args.SNP_output)

if __name__ == '__main__':
    main()