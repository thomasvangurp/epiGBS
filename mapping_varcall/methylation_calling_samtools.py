#!/usr/bin/env python
# __author__ = 'Bjorn Wouters'
# Date created: 16/09/2014 (europe date)
# Function: Base calling for methylated nucleotides and SNP's
#Python version: 2.7.3
#External modules: vcf, HTSeq
#Known bugs: None
#Modifications: None

import vcf
from vcf import utils
from itertools import izip
from Bio import SeqIO
import argparse
import re
import subprocess
import heapq

# argparse used for commandline interpretation.
def parse_args():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(description='Process input files')
    parser.add_argument('-r', '--reference', type=str, nargs='?', default=None,
                        help='reference genome input.')
    parser.add_argument('-w', '--watson', type=str, nargs='?', default=None,
                        help='watson (top strand) .vcf file input.')
    parser.add_argument('-c', '--crick', type=str, nargs='?', default=None,
                        help='crick (bottom strand) .vcf file input.')
    parser.add_argument('-m', '--methylation_output', type=str, nargs='?', default=None,
                        help='Methylation vcf file output name')
    parser.add_argument('-s', '--SNP_output', type=str, nargs='?', default=None,
                        help='SNP vcf file output name')
    parser.add_argument('-heat', '--heatmap_output', type=str, nargs='?', default=None,
                        help='Heatmap igv file output name')
    parser.add_argument('-methylation_called', '-methylation_called', type=str, nargs='?', default=None,
                        help='Called sample information .txt file output name')
    parser.add_argument('-snp_called', '-snp_called', type=str, nargs='?', default=None,
                        help='Called sample information .txt file output name')
    parser.add_argument('-qual', '--min_quality', type=int, nargs='?', default=0,
                            help='Minimum Freebayes call quality before processing, default: 0')
    args = parser.parse_args()
    return args


def main():
    """Main function loop"""
    args = parse_args()
    files = parse_vcf(args)
    zip_tabix(args)


def filter_records(records):
    """filter variant calls on strand placement"""
    out_records = []
    for record in records:
        if not record.is_monomorphic:
            print 'Number of reference observations on the forward strand:\t%s'%record.INFO['SRF']
            print 'Number of reference observations on the reverse strand:\t%s'%record.INFO['SRR']
            print 'Strand balance probability for the reference allele:\t%s'%record.INFO['SRP']
            print 'Number of alternate observations on the forward strand:\t%s'%record.INFO['SAF']
            print 'Number of alternate observations on the reverse strand:\t%s'%record.INFO['SAR']
            print 'Strand balance probability for the alternate allele:\t%s'%record.INFO['SAP']

            print record


def parse_vcf(args):
    """
    Iterates through the given two vcf files. For each unique variant
    the scripts determines if it's a methylation or a SNP call or both.
    """
    methyl_called_file = open(args.methylation_called, 'w')
    snp_called_file = open(args.snp_called, 'w')
    watson_file = vcf.Reader(open(args.watson, 'r'))  # Watson input .vcf file
    crick_file = vcf.Reader(open(args.crick, 'r'))  # Crick input .vcf file
    reference_genome = SeqIO.to_dict(SeqIO.parse(args.reference, "fasta"))  # reference genome .fasta file
    methylation_file = vcf.Writer(open(args.methylation_output, 'w'),  # methylation output file
                                  watson_file, )
    snp_file = vcf.Writer(open(args.SNP_output, 'w'),  # snp output file
                          watson_file, )
    igv_file = open(args.heatmap_output, 'w')  # heatmap (igv) output file

    # Creates the header of the .igv file
    igv_file.write('#type=DNA_METHYLATION\n')
    samples = '\t'.join(watson_file.samples)
    igv_file.write('Chromosome\tStart\tEnd\tFeature\t' + samples + '\n')
    #Iterates through the watson file. If the variant is methylated (True); the script will write
    #the record in the methylation file. If the variant is not methylated (False); the script will write the
    #record in the snp file
    call_base = CallBase(watson_file, crick_file,
                         reference_genome, methylation_file,
                         snp_file, igv_file, methyl_called_file, snp_called_file)
    #use vcf.utils.walk_together
    combined_records = vcf.utils.walk_together(watson_file,crick_file)
    old_chrom = None
    for records in combined_records:
        if None not in records:
            #both records need to be present and valid
            # qsum = sum(record.QUAL for record in records)
            # filter_records(records)
            call_base.watson_record,call_base.crick_record = records

            #process IGV records when done with CHROM
            if records[0].CHROM != old_chrom:
               for site in call_base.methylated_records:
                   write_igv_file(call_base, site)
               call_base.methylated_records = []

            #Call methylation / SNPs: method of callbase class
            #TODO: check quality parameters elsewhere
            if call_base.watson_record.REF in ['C', 'G']:
                return_code = call_base.methylation_calling()
                if not return_code:
                    continue
            return_code = call_base.filter_snps()

            if not return_code:
                continue

            call_base.write_records()

            # if int(call_base.watson_record.CHROM) > 1000:
            #     break

            #TODO If there are no SNP's in the cluster/chromosome, the igv file needs to be written without a sliding window.
            old_chrom = records[0].CHROM



def make_empty_sample(sample):
    """
    Returns an empty pyVCF sample for the given sample site.
    """
    return vcf.model._Call(sample.site,
                           sample.sample,
                           tuple([None]*len(sample.site.FORMAT.split(':'))))


def make_sample(sample, genotype):
    # TODO: CG context
    """
    Returns the current given sample with the corrected genotype.
    """
    return vcf.model._Call(sample.site,
                           sample.sample,
                           (genotype,sample.data.DP, sample.data.AD, sample.data.RO,
                            sample.data.AO))


def combine_record_samples(sample1, sample2):
    """
    Returns two sample records values together.
    """
    #TODO:implement merging function for SNP calling from watson/crick with methylation polymorphisms
    #TODO: Create records given that QR, QA, RO, AO are not available, only use AD and merge here. Check how AD records
    #TODO: are generated when multiple calls exist.
    if not sample1.called:
        depth = sample2.data.DP
        ref_observations = sample2.data.RO
        alt_observations = sample2.data.AO
    elif not sample2.called:
        depth = sample1.data.DP
        ref_observations = sample1.data.RO
        alt_observations = sample1.data.AO
    else:
        #determine on which sample we should base output record
        alt_observations = []
        depth = sample1.data.DP + sample2.data.DP
        ref_observations = sample1.data.RO + sample2.data.RO
        if len(sample1.site.ALT) > len(sample2.site.ALT):
            out_prim = sample1
            out_sec = sample2
        elif len(sample2.site.ALT) > len(sample1.site.ALT):
            out_prim = sample2
            out_sec = sample1
        elif sample1.site.ALT == sample2.site.ALT:
            out_prim = sample2
            out_sec = sample1
        elif sample1.site.ALT in [None,[None]]:
            out_prim = sample2
            out_sec = sample1
        elif sample2.site.ALT in [None,[None]]:
            out_prim = sample1
            out_sec = sample2
        else:
            #TODO: make SNP calls when methylation changes alt alleles
            out_prim = sample1
            out_sec = sample2
        try:
            for n,nt in enumerate(out_prim.site.ALT):
                if nt not in out_sec.site.ALT:
                    if str(nt) == 'C' and 'T' in [str(i) for i in out_sec.site.ALT]:
                        print 'SNP masked by meth polymorphism'
                    if str(nt) == 'T' and 'C' in [str(i) for i in out_sec.site.ALT]:
                        print 'SNP masked by meth polymorphism'
                    alt_observations.append(out_prim.data.AO[n])
                else:
                    sec_index = out_sec.site.ALT.index(nt)
                    sec_count = out_sec.data.AO[sec_index]
                    prim_count = out_prim.data.AO[n]
                    out_count = 0
                    for add in [prim_count,sec_count]:
                        if type(add) == type(1):
                            out_count += add
                    alt_observations.append(out_count)
            header = sample1.site.FORMAT
            header_list = header.split(':')
            call_data = vcf.model.make_calldata_tuple(header_list)
            values = [out_prim.data.GT,depth,out_prim.data.AD,ref_observations,
                      alt_observations]
            model = vcf.model._Call(out_prim.site,
                                out_prim.sample,
                                call_data(*values))
            return model
        except TypeError:
            return out_prim


def write_igv_file(call_base, methyl_record):
    """
    Object function to write the record values to a heatmap (.igv file)
    """
    total_samples = 0
    processed_samples = {'CG': dict(), 'CHG': dict(), 'CHH': dict(), '.': dict()}
    methyl_pos = methyl_record.POS
    # If there are no SNP's at the given chromosome/cluster there is no dictionary key, so there are also no snp
    # positions in the chromosome/cluster
    if call_base.snp_record_dict.has_key(methyl_record.CHROM):
        snp_pos = set(pos.POS for pos in call_base.snp_record_dict[methyl_record.CHROM])
    else:
        snp_pos = set()

    def write_to_file(dataset, context):
        """
        Writes the given dataset to an .igv file with the right context.
        Dataset is a list with ratio's.
        """
        chr = methyl_record.CHROM
        # call_base.igv_file.write(str(chr) + '\t' + str(methyl_pos - 1) + '\t' + str(methyl_pos) + '\t' + context)
        out = '%s\t%s\t%s\t%s\t'%(chr,methyl_pos-1,methyl_pos,context)
        out += '\t'.join([str(v) for v in dataset]) + '\n'
        call_base.igv_file.write(out)
        # for value in dataset:
        #     if value == '.': # If there's no call, the value is set to '.'
        #         call_base.igv_file.write('\t.')
        #     else:
        #         call_base.igv_file.write('\t'+str(value))
        # call_base.igv_file.write('\n')

    def write_sample_data(sample_data, context, total_samples):
        """
        Writes a file with for each called sample the total calls and methylated calls.
        """
        chr = methyl_record.CHROM
        call_base.samples_called.write(str(chr) + '\t' + str(methyl_pos) + '\t' + context + '\t' + str(total_samples))
        for value in sample_data:
            if not isinstance(value, basestring):
                call_base.samples_called.write('\t'+'\t'.join(map(str, value)))
            else:
                call_base.samples_called.write(value)
        call_base.samples_called.write('\n')

    def calc_context(ref,pos):
        """
        If there are no SNP's neighbouring the methylation call, the context can be called
        by using the reference genome.
        """
        slice_start = max(pos-3,0) #negative positions excluded
        slice_end = pos+2
        reference_bases = call_base.reference_genome[methyl_record.CHROM].seq[slice_start:slice_end]
        if ref == 'G':
            ref_context = reference_bases[0:2][::-1]
            if re.match('C.', str(ref_context)):
                context = 'CG'
            elif re.match('[ATG]C', str(ref_context)):
                context = 'CHG'
            else:
                context = 'CHH'
        elif ref == 'C':
            ref_context = reference_bases[3:5]
            if re.match('G.', str(ref_context)):
                context = 'CG'
            elif re.match('[ATC]G', str(ref_context)):
                context = 'CHG'
            else:
                context = 'CHH'
        else:
            context = '.'
        return context

    def calc_methylated_observations(sample):
        if methyl_record.REF in 'CG':
            return sample.data.RO
        else:
            if sample.gt_bases[-1] in 'CG':
                if not isinstance(sample.data.AO, int):
                    float_index = int(sample.gt_alleles[1])-1
                    return sample.data.AO[float_index]
                else:
                    sample.data.AO
            else:
                return 0

    def calc_alt_observations(sample):
        try:
            if sample.site.REF == 'C':
                if 'T' in sample.site.ALT:
                    if not isinstance(sample.data.AO, int):
                        float_index = sample.site.ALT.index('T')
                        return sample.data.RO+sample.data.AO[float_index]
                    else:
                        return sample.data.RO+sample.data.AO
                else:
                    return sample.data.RO
            elif sample.site.REF == 'G':
                if 'A' in sample.site.ALT:
                    if not isinstance(sample.data.AO, int):
                        float_index = sample.site.ALT.index('A')
                        return sample.data.RO+sample.data.AO[float_index]
                    else:
                        return sample.data.RO+sample.data.AO
                else:
                    return sample.data.RO
            else:
                if sample.gt_bases[-1] == 'C':
                    if not isinstance(sample.data.AO, int):
                        C_index = sample.site.ALT.index('C')
                        T_index = sample.site.ALT.index('T')
                        return sample.data.AO[C_index]+sample.data.AO[T_index]
                    else:
                        return sample.data.AO
                elif sample.gt_bases[-1] == 'G':
                    if not isinstance(sample.data.AO, int):
                        G_index = sample.site.ALT.index('G')
                        A_index = sample.site.ALT.index('A')
                        return sample.data.AO[G_index]+sample.data.AO[A_index]
                    else:
                        return sample.data.AO
                else:
                    return 0
        except TypeError:
            return sample.data.RO
        except IndexError:
            return sample.data.RO

    def calc_ratio(sample):
        """
        Calculated the methylation ratio by dividing the methylated counts by total (methylated+unmethylated) counts.
        In case of a SNP the count of the alternate allele will not be taken into account for determining the ratio!
        """
        try:
            if methyl_record.REF == 'C':
                if sample.data.RO == 0:
                    return 0
                if 'T' in methyl_record.ALT:
                    if not isinstance(sample.data.AO, int):
                        float_index = methyl_record.ALT.index('T')
                        float_number = float(sample.data.RO) / float((sample.data.RO+sample.data.AO[float_index]))
                    else:
                        float_number = float(sample.data.RO) / float((sample.data.RO+sample.data.AO))
                else:
                    float_number = 1
            elif methyl_record.REF == 'G':
                if sample.data.RO == 0:
                    return 0
                if 'A' in methyl_record.ALT:
                    if not isinstance(sample.data.AO, int):
                        float_index = methyl_record.ALT.index('A')
                        float_number = float(sample.data.RO) / float((sample.data.RO+sample.data.AO[float_index]))
                    else:
                        float_number = float(sample.data.RO) / float((sample.data.RO+sample.data.AO))
                else:
                    float_number = 1
            else:
                if sample.gt_bases[-1] == 'C':
                    if not isinstance(sample.data.AO, int):
                        C_index = methyl_record.ALT.index('C')
                        T_index = methyl_record.ALT.index('T')
                        float_number = float(sample.data.AO[C_index]) / float((sample.data.AO[C_index]+sample.data.AO[T_index]))
                    else:
                        float_number = 1
                elif sample.gt_bases[-1] == 'G':
                    if not isinstance(sample.data.AO, int):
                        G_index = methyl_record.ALT.index('G')
                        A_index = methyl_record.ALT.index('A')
                        float_number = float(sample.data.AO[G_index]) / float((sample.data.AO[G_index]+sample.data.AO[A_index]))
                    else:
                        float_number = 1
                else:
                    float_number = 0
            ratio = "%.3f" % float_number
        except IndexError:
            ratio = "0.00"
        except TypeError:
            ratio = "0.00"
        return ratio

    def get_snp_sample_genotype(current_pos, sample_name):
        """
        Returns the sample of the SNP at the given position.
        """
        record = next((record for record in call_base.snp_record_dict[methyl_record.CHROM]
                       if record.POS == current_pos), None)
        sample = next((sample for sample in record.samples if sample.sample == sample_name), None)
        if not sample.called:
            return None
        try:
            return sample.gt_bases
        except IndexError:
            #TODO: break here and see what causes this!
            return None

    for sample in methyl_record.samples:
        sample_name = sample.sample
        if sample.called:
            total_samples += 1
            if methyl_record.REF in 'CG':
                ref = sample.site.REF
            else:
                ref = sample.gt_alleles[-1]

            if ref == 'C':
                p1 = 1
                p2 = 2
                alt_chk = 'G'
                ref_chk = 'A'
            elif ref == 'G':
                p1 = -1
                p2 = -2
                alt_chk = 'C'
                ref_chk = 'T'
            else:
                continue
            context = calc_context(ref,methyl_record.POS)
            if methyl_pos+p1 not in snp_pos and methyl_pos+p2 not in snp_pos:
                #No SNPs are found for this sample in 2 downstream adjacent positions
                context = calc_context(ref,methyl_record.POS)
                #TODO: define output for raw numbers instead of ratio!
                ratio = calc_ratio(sample)
                processed_samples[context].update({sample_name: ratio})
                continue
            if methyl_pos+p1 in snp_pos:
                gt = get_snp_sample_genotype(methyl_pos+p1, sample_name)
                if gt:
                    alt = gt[-1]
                    ref = gt[0]
                    if alt == alt_chk:
                        if ref == ref_chk:
                            context = '.'
                        else:
                            #Assume that SNP is CG as this is most common. TODO: determine if valid!
                            context = 'CG'
                    ratio = calc_ratio(sample)
                    processed_samples[context].update({sample_name: ratio})
                    continue
            if methyl_pos+p2 in snp_pos:
                gt = get_snp_sample_genotype(methyl_pos+p2, sample_name)
                if gt:
                    alt = gt[-1]
                    ref = gt[0]
                    if alt == alt_chk:
                        if ref == ref_chk:
                            context = '.'
                        else:
                            context = 'CHG'
                ratio = calc_ratio(sample)
                processed_samples[context].update({sample_name: ratio})

    for context in processed_samples:
        #Context can be CG,CHG,CHH or unknown:.
        if processed_samples[context]: #false if empty
            ratio_dataset = []
            sample_data = []
            for sample in methyl_record.samples:
                sample_name = sample.sample
                if sample_name in processed_samples[context]:
                    sample_ratio = processed_samples[context][sample_name]
                    ratio_dataset.append(sample_ratio)
                    sample_data.append([calc_methylated_observations(sample), calc_alt_observations(sample)])
                else:
                    ratio_dataset.append('.')
                    sample_data.append('\tNone\tNone')
            write_to_file(ratio_dataset, context)
            write_sample_data(sample_data, context, total_samples)


def write_snp_file(call_base, snp_record):
    snp_dict = {'A': dict(), 'T': dict(), 'G': dict(), 'C': dict()}

    total_samples = 0

    for sample in snp_record.samples:
        if sample.called:
            try:
                #new field is sample.AD sample.data.AO does no longer exist
                if isinstance(sample.data.AO,list) and len(sample.data.AO) > 1:
                    if sum(sample.data.AO) > 0:
                        total_samples += 1
                    for i, alt_base in enumerate(sample.site.ALT):
                        snp_dict[str(alt_base)].update({sample.sample: sample.data.AO[i]})
                #TODO: make sure sample.site contains relevant alt alleles
                elif sample.data.AO in [None,[None]] or sample.site in [None,[None]]:
                    continue
                elif len(sample.data.AO) == 1:
                    snp_dict[str(sample.site.ALT[0])].update({sample.sample:sample.data.AO[0]})
                elif isinstance(sample.data.AO, int):
                    if sample.data.AO > 0:
                        total_samples += 1
                        alt_base = str(sample.site.ALT[0])
                        snp_dict[alt_base].update({sample.sample: sample.data.AO})
            except KeyError:
                continue


    call_base.snp_output_file.write(str(snp_record.CHROM) + '\t' + str(snp_record.POS) +
                                    '\t' + str(total_samples))

    for sample in snp_record.samples:
        sample_name = sample.sample
        for base in snp_dict.keys():
            if sample_name in snp_dict[base]:
                call_base.snp_output_file.write('\t'+str(snp_dict[base][sample_name]))
            else:
                call_base.snp_output_file.write('\t'+str(0))
    call_base.snp_output_file.write('\n')


class CallBase(object):
    def __init__(self, watson_file, crick_file, reference_genome,
                 methylation_file, snp_file, igv_file, called_file, snp_output):
        self.watson_file = watson_file
        self.crick_file = crick_file
        self.samples_called = called_file
        self.snp_output_file = snp_output
        self.samples_called.write('chr\tpos\tcontext\tsamples_called')
        self.snp_output_file.write('chr\tpos\tsamples_called')
        for sample in watson_file.samples:
            self.samples_called.write('\t'+sample+'_methylated')
            self.samples_called.write('\t'+sample+'_total')
            self.snp_output_file.write('\t'+sample+'_A')
            self.snp_output_file.write('\t'+sample+'_T')
            self.snp_output_file.write('\t'+sample+'_G')
            self.snp_output_file.write('\t'+sample+'_C')
        self.samples_called.write('\n')
        self.snp_output_file.write('\n')
        self.watson_record = None
        self.crick_record = None
        self.reference_genome = reference_genome
        self.methylation_file = methylation_file
        self.snp_file = snp_file
        self.igv_file = igv_file
        self.methylation_calls = {'C': set(['T/C', 'C/C', 'T/T']),
                                  'G': set(['A/G', 'G/G', 'A/A'])}
        self.snp_record_dict = dict()
        self.methylated_records = list()

    def set_offsets(self, qual_offset, min_alt_observations):
        """
        Sets the offset parameters of the call_base object.
        """
        self.qual_offset = qual_offset
        self.min_alt_observations = min_alt_observations
        return self.qual_offset

    def check_change_samtools_call(self, sample):
        """Check variant calling in samtools provided sample"""
        empty_sample = make_empty_sample(sample)
        GT = None
        #if genotype call is homozygous wheras the percentage of alt counts is higher than 5 %, change genotype call!
        if not sample.is_het and type(sample.data.AD) == type([]):
            if sample.data.AD[0] == sum(sample.data.AD):
                pass
            else:
                try:
                    if sample.data.AD[0] / float(sample.data.DP) > 0.05:
                        max_alt = max(sample.data.AD[1:])
                        alt_pos = sample.data.AD.index(max_alt)
                        GT = '0/%s'%alt_pos
                        values = [GT]
                except ZeroDivisionError:
                    pass
        if not GT:
            values = [sample.data.GT]
        if 'PL' in sample.data._fields:
            header = ['GT','DP','AD','RO','AO']
            call_data = vcf.model.make_calldata_tuple(header)
            values += [sample.data.DP,
                      sample.data.AD,
                      ]
            if type(sample.data.AD) == type(1):
                values += [sample.data.AD,0]
            elif sample.site.ALT[0].type == '*':
                values += [sample.data.AD[0],0]
            else:
                values += [sample.data.AD[0]]
                if len(sample.data.AD[1:]) == 1:
                    values += [sample.data.AD[1]]
                else:
                    values += [sample.data.AD[1:]]
        elif sample.data.GT == '0/0' and type(sample.data.AD) == type(1):
            header = ['GT','DP','AD','RO','AO']
            call_data = vcf.model.make_calldata_tuple(header)
            values += [sample.data.DP,
                      sample.data.AD,
                      #RO is position 0 of AD
                      sample.data.AD,
                      #AO should be 0
                      0
                      ]
        else:
            pass
        #change sample.site so that it contains all the FORMAT fields in use
        sample.site.FORMAT = ':'.join(header)
        model = vcf.model._Call(sample.site,
                                    sample.sample,
                                    call_data(*values))
        return model

    def call_samples(self, record):
        samples_out = []
        alleles_observed = []
        header = ['GT','DP','AD','RO','AO']
        call_data = vcf.model.make_calldata_tuple(header)
        for sample in record.samples:
            out_count = {}
            for pos,nt in enumerate(record.ALT):
                count = sample.data.AD[pos+1]
                if count:
                    #ADF is Allelic Depth on forward strand
                    if 'ADFR' in sample.data._fields:
                        #Get minimum count for both forward and reverse mapped reads
                        min_count = float(min([sample.data.ADF[pos + 1], sample.data.ADR[pos + 1]]))
                        #if this count is 0 or lower than 5%, do not take the allele into account.
                        min_strand_cover = min([sum(sample.data.ADF), sum(sample.data.ADR)])
                        #min cover of one of the strands is too low. take total
                        if min_strand_cover < 10:
                            #not enough data on one strand, use both
                            #TODO: check if this could be done differently!
                            if count / float(sample.data.DP) > 0.05:
                                out_count[str(nt)] = count
                            continue
                        if min_count / min_strand_cover > 0.05:
                            out_count[str(nt)] = count
                    else:
                        #method for when no forward or reverse allelic depth is known.
                        if count / float(sample.data.DP) > 0.05:
                            out_count[str(nt)] = count
            if out_count == {}:
                if sum(sample.data.AD) == 0:
                    #sample is not called as the sum of all calls is 0
                    #TODO: implement treshold here based on number of observations for valid allele call?
                    GT = './.'
                else:
                    #only reference bases were found, sample is homozygous reference
                    GT = '0/0'
            elif len(out_count) == 1:
                #only one alternate allele found
                alt_pos = [str(nt) for nt in record.ALT].index(out_count.keys()[0]) + 1
                if sample.data.AD[0] / float(sample.data.DP) > 0.05:
                    #reference allele is present more than 5%
                    if 'ADFR' not in sample.data._fields:
                        GT = '0/%s'%alt_pos
                    else:
                        #Check if reference allele is indeed observed in both forward and reverse mapping reads
                        GT = []
                        for i in [0,1]:

                            try:
                                if min([float(sample.data.ADF[i])/sum(sample.data.ADF),
                                        float(sample.data.ADR[i])/sum(sample.data.ADR)]) > 0.05:
                                    GT.append(str(i))
                            except ZeroDivisionError:
                                pass
                        if GT == []:
                            GT = './.'
                        else:
                            GT = '%s/%s'%(GT[0],GT[-1])
                else:
                    #only alternate allele is present
                    GT = '%s/%s'%(alt_pos,alt_pos)
            else:
                if sample.data.AD[0] / float(sample.data.DP) > 0.05:
                    GT = ['0']
                    i = 1
                else:
                    i = 2
                    GT = []
                #more than one alternate allele found determine top 2 alleles
                max_values = sorted(out_count.values())[:i]
                for value in max_values:
                    for nt,count in out_count.items():
                        if count == value:
                            alt_pos = [str(v) for v in record.ALT].index(str(nt)) + 1
                            GT.append(str(alt_pos))
                GT = '/'.join(GT)
            if 'ADFR' in sample.data._fields and min([sum(sample.data.ADF),sum(sample.data.ADR)]) > 10:
                AD = []
                for obs_fw,obs_rev in zip(sample.data.ADF,sample.data.ADR):
                    try:
                        fw_ratio =  obs_fw / float(sum(sample.data.ADF))
                    except ZeroDivisionError:
                        fw_ratio  = 0
                    try:
                        rev_ratio = obs_rev / float(sum(sample.data.ADR))
                    except ZeroDivisionError:
                        rev_ratio = 0
                    if min([fw_ratio,rev_ratio]) > 0.05:
                        AD.append(obs_fw + obs_rev)
                    else:
                        AD.append(0)
                AO = [count for (count,allele) in zip(AD[1:-1],sample.site.ALT) if allele in record.ALT]
            elif 'ADFR' in sample.data._fields:
                AD = sample.data.AD
                AO = AD[1:]
                AO = [count for (count,allele) in zip(AD[1:-1],sample.site.ALT) if allele in record.ALT]
            else:
                AO = [count for (count,allele) in zip(sample.data.AD[1:-1],sample.site.ALT) if allele in record.ALT]
                AD = sample.data.AD
            if AO == []:
                AO = None
            values = [GT,
                      sample.data.DP,
                      sample.data.AD,
                      AD[0],
                      AO
                      ]
            sample.site.FORMAT = ':'.join(header)
            model = vcf.model._Call(sample.site,
                                        sample.sample,
                                        call_data(*values))
            model.site = record
            samples_out.append(model)
        record.samples = samples_out
        try:
            assert len(record.samples[0].data.AO) == len(record.samples[0].site.ALT)
        except TypeError:
            assert len(record.samples[0].site.ALT) == 1


        return record

    def call_genotypes(self):
        """"call genotypes for watson and crick"""
        # determine which alt records to keep for watson
        # only alt records that are present with more than 5% in any sample are preserved
        keep_nt = list()
        for pos, nt in enumerate(self.watson_record.ALT[:-1]):
            max_pct = max([sample.data.AD[pos+1] / float(sample.data.DP) for sample in self.watson_record.samples if
                           sample.data.DP != 0])
            if max_pct > 0.05:
                keep_nt.append(nt)
        if not keep_nt:
            keep_nt = [None]
        self.watson_record.ALT = keep_nt
        self.watson_record.alleles = [self.watson_record.REF] + keep_nt

        self.watson_record = self.call_samples(self.watson_record)

        keep_nt = list()
        for pos, nt in enumerate(self.crick_record.ALT[:-1]):
            max_pct = max([sample.data.AD[pos+1] / float(sample.data.DP) for sample in self.crick_record.samples if
                           sample.data.DP != 0])
            if max_pct > 0.05:
                keep_nt.append(nt)
        if not keep_nt:
            keep_nt = [None]
        self.crick_record.ALT = keep_nt
        self.crick_record.alleles = [self.crick_record.REF] + keep_nt

        self.crick_record = self.call_samples(self.crick_record)

    def methylation_calling(self):
        """
        Main base calling algorithm. Determines methylation/SNP status for each sample having a watson and crick record.
        """

        # If the sample is methylated, the processed_samples will be filled under the methylated key, or when it's
        # a normal SNP, it will be filed under the snp key or when it's both, under both.

        self.processed_samples = {key: {'methylated': None} for key in self.watson_file.samples}
        # Determine the reference base at the position of the VCF call.
        ref_base = self.watson_record.REF
        # loop watson and crick record for combined samples.
        if min([self.watson_record.INFO['DP'], self.crick_record.INFO['DP']]) == 0:
            return None
        # if self.watson_record.REF not in ['C', 'G']:
        #     watson_alt = sum(self.watson_record.INFO['AD'][1:])/float(self.watson_record.INFO['DP'])
        #     crick_alt = sum(self.crick_record.INFO['AD'][1:])/float(self.crick_record.INFO['DP'])
        #     if max(watson_alt, crick_alt) < 0.05:
        #         return None
        self.call_genotypes()
        for watson_sample, crick_sample in izip(self.watson_record, self.crick_record):
            # If there are is no call for both the watson and crick record sample, continue as we can not determine
            # whether polymorphism is a SNP/methylation polymorphism.
            if not watson_sample.called or not crick_sample.called:
                continue
            # Assigning the right alt base to the records.
            alt_watson = watson_sample.gt_bases.split('/')[1]
            alt_crick = crick_sample.gt_bases.split('/')[1]
        
            sample_name = watson_sample.sample
            #TODO: move SNP calling to separate algorithm.
            if ref_base == 'C':
                if alt_crick == 'C' and alt_watson in 'CT':
                    #Methylation in watson C/T No polymorphism in Crick: methylation
                    self.processed_samples[sample_name]['methylated'] = watson_sample
                    # try:
                    #     #'AC' is not present in homozygous situations.
                    #     if sum(crick_sample.site.INFO['AD'][1:]) > 0:
                    #         #The Alternate alleles need to be called in at least one sample to be valid!
                    #         # self.processed_samples[sample_name]['snp'] = crick_sample
                    # except KeyError:
                    #     pass
                elif alt_crick == 'A':
                    if alt_watson == 'A':
                        pass
                        #Both watson and crick records should contain information on this alternate allele
                        #records are combined and written
                        # combined_record = combine_record_samples(watson_sample, crick_sample)
                        # self.processed_samples[sample_name]['snp'] = combined_record
                    else:
                        #alt_crick contains another base, this SNP is valid if present in watson.
                        if alt_watson in crick_sample.site.ALT:
                            #Both watson and crick contain the same alternate allele, the SNP is real?
                            alt_index = crick_sample.site.ALT.index(alt_watson)
                            try:
                                crick_alt_pct = crick_sample.data.AO[alt_index] / float(crick_sample.data.DP)
                            except TypeError:
                                crick_alt_pct = crick_sample.data.AO / float(crick_sample.data.DP)
                            alt_index = watson_sample.site.ALT.index(alt_watson)
                            try:
                                watson_alt_pct = watson_sample.data.AO[alt_index] / float(watson_sample.data.DP)
                            except TypeError:
                                watson_alt_pct = watson_sample.data.AO / float(watson_sample.data.DP)
                            if watson_alt_pct != 0.0 and crick_alt_pct != 0.0:
                                if max(crick_alt_pct,watson_alt_pct)/min(crick_alt_pct,watson_alt_pct)< 1.5:
                                    continue
                                    # self.processed_samples[sample_name]['snp'] = crick_sample
                                    #TODO: merge alt counts for watson and crick here
                            elif crick_alt_pct == 0.0:
                                #REF:C watson C/T/A called as C/T crick C/A
                                #We can call both SNP and methylation. SNP from crick reliable
                                #Watson information on C/T ratio informative for methylation call
                                # self.processed_samples[sample_name]['snp'] = crick_sample
                                self.processed_samples[sample_name]['methylated'] = watson_sample
                            else:
                                pass
                                #Can this occur? TODO: check if this can be true
                        else:
                            #alt_crick = 'A' ref == 'C' and alt_watson not seen in alt_watson.
                            #Likely a C/G SNP for which one or both alleles got converted watson C>T crick G>A
                            #TODO: make new C/G SNP record combining watson and crick variation.

                            pass
                elif alt_watson == 'G' and alt_crick in 'AG':
                    #C/G polymorphism in watson, C/G or C/A in Crick
                    #SNP information from watson, Methylation information from Crick
                     #TODO: After SNP, methylation > SNP
                    # self.processed_samples[sample_name]['snp'] = watson_sample
                    self.processed_samples[sample_name]['methylated'] = crick_sample
                # elif crick_sample.gt_bases == 'C/A' and watson_sample.gt_bases == 'G/T':
                #     #Situation: C/G SNP, both alleles unmethylated.on watson C gets converted to T on watson G to A..
                #     #TODO: create combined record for both Watson and crick replacing alt in watson A=> G
                #     # and crick C=>T
                #     pass
                elif alt_watson == 'T':
                    if alt_crick == 'T':
                    #C/T variant in both watson and crick: SNP ==> only information from crick is reliable
                        # self.processed_samples[sample_name]['snp'] = crick_sample
                        if set(watson_sample.gt_bases.replace('/','')) == set(['T']):
                            self.processed_samples[sample_name]['methylated'] = watson_sample
                    if alt_crick != 'T' and alt_crick != 'C':
                    #Watson contains C/T methylation polymorphism and potentially other information
                    #Step 1: We can call Methylation polymorphism, crick is not C/T!
                        self.processed_samples[sample_name]['methylated'] = watson_sample
                        #Step 2. Is the SNP supported in watson strand?
                        #Can we call the SNP?
                        # if alt_crick in watson_sample.site.ALT:
                        #     alt_index = watson_sample.site.ALT.index(alt_crick)
                        #     try:
                        #         alt_count = watson_sample.data.AO[alt_index]
                        #     except TypeError:
                        #         alt_count = watson_sample.data.AO
                        #     try:
                        #         t_count = watson_sample.data.AO[watson_sample.site.ALT.index('T')]
                        #     except TypeError:
                        #         t_count = watson_sample.data.AO
                        #     if alt_count > t_count:
                        #         self.processed_samples[sample_name]['snp'] = crick_sample
                        #         #TODO: merge alt counts for watson and crick here
            elif ref_base == 'G':
                #Watson is homozygous reference (i.e. no SNP) and crick has Methylation variation
                if alt_watson == 'G' and alt_crick in 'GA':
                    self.processed_samples[sample_name]['methylated'] = crick_sample
                    # try:
                    #     #If one or more sample has an alternate allele call it for this individual as well
                    #     if sum(watson_sample.site.INFO['AD'][1:]) > 0:
                    #          self.processed_samples[sample_name]['snp'] = watson_sample
                    # except KeyError:
                    #     pass
                elif alt_crick == 'A' and alt_watson == 'A':
                    #The watson sample contains information on a SNP, G/A crick not reliable for SNP
                    # self.processed_samples[sample_name]['snp'] = watson_sample
                    #The crick allele can only be queried for methylation variation if it is fully converted
                    #this means that no G can be in the Genotype.
                    if set(crick_sample.gt_bases.replace('/','')) == set(['A']):
                        self.processed_samples[sample_name]['methylated'] = crick_sample
                # elif alt_crick == 'T' and alt_watson == 'T':
                #     combined_record = combine_record_samples(watson_sample,crick_sample)
                #     self.processed_samples[sample_name]['snp'] = combined_record
                #     #TODO: After SNP, methylation > SNP
                elif alt_watson == 'C' and alt_crick in 'CT':
                    # self.processed_samples[sample_name]['snp'] = crick_sample
                    self.processed_samples[sample_name]['methylated'] = watson_sample
            # elif ref_base == 'T':
            #     if alt_watson == 'T' and alt_crick == 'T':
            #         #both samples have the reference allele
            #         #Determine the allele count for the Alternate allele for the site
            #         #If sum of alternate allele count > 0 ?? add relevant Allele to SNP output
            #         # if 'AC' in watson_sample.site.INFO and 'AC' in crick_sample.site.INFO:
            #         #     #Only Crick contains the SNP, TODO: check if valid!
            #         #     if sum(crick_sample.site.INFO['AC']) == 0 and sum(watson_sample.site.INFO['AC']) > 0:
            #         #         self.processed_samples[sample_name]['snp'] = watson_sample
            #         #     #Only Watson contains the SNP, TODO: check if valid!
            #         #     elif sum(crick_sample.site.INFO['AC']) > 0 and sum(watson_sample.site.INFO['AC']) == 0:
            #         #         self.processed_samples[sample_name]['snp'] = crick_sample
            #         #Both contain the SNP, combine records, TODO: check if combination is warranted!
            #         #Both watson and crick should contain the same type of non reference allele
            #         combined_record = combine_record_samples(watson_sample,crick_sample)
            #         self.processed_samples[sample_name]['snp'] = combined_record
            #         # elif 'AC' in watson_sample.site.INFO:
            #         #     #Only watson contains the SNP,
            #         #     #TODO: assert if valid. Only true if SNP is T/C!
            #         #     self.processed_samples[sample_name]['snp'] = watson_sample
            #         # elif 'AC' in crick_sample.site.INFO:
            #         #     #Only watson contains the SNP,
            #         #     #TODO: assert if valid. Only true if SNP is T/G!
            #         #     self.processed_samples[sample_name]['snp'] = crick_sample
            #     elif alt_watson == 'A' and alt_crick == 'A':
            #         #TODO: check what is causing the error in combined_record it's not writing now.
            #         combined_record = combine_record_samples(watson_sample,crick_sample)
            #         # self.processed_samples[sample_name]['snp'] = combined_record
            #         self.processed_samples[sample_name]['snp'] = watson_sample
            #     elif alt_watson == 'G' and alt_crick in 'GA':
            #         self.processed_samples[sample_name]['snp'] = watson_sample
            #         # self.processed_samples[sample_name]['methylated'] = crick_sample
            #     elif alt_watson == 'C' and alt_crick in 'CT':
            #         self.processed_samples[sample_name]['snp'] = crick_sample
            #     elif 'A' in alt_watson + alt_crick and ref_base in alt_watson + alt_crick:
            #         #one of the samples
            #         combined_record = combine_record_samples(watson_sample,crick_sample)
            #         self.processed_samples[sample_name]['snp'] = combined_record
            # elif ref_base == 'A':
            #     if alt_watson == 'A' and alt_crick == 'A':
            #         #both samples have the reference allele
            #         #Determine the allele count for the Alternate allele for the site
            #         #If sum of alternate allele count > 0 ?? add relevant Allele to SNP output
            #         if 'AC' in watson_sample.site.INFO and 'AC' in crick_sample.site.INFO:
            #             #Only Crick contains the SNP, TODO: check if valid!
            #             if sum(crick_sample.site.INFO['AC']) == 0 and sum(watson_sample.site.INFO['AC']) > 0:
            #                 self.processed_samples[sample_name]['snp'] = watson_sample
            #             #Only Watson contains the SNP, TODO: check if valid!
            #             elif sum(crick_sample.site.INFO['AC']) > 0 and sum(watson_sample.site.INFO['AC']) == 0:
            #                 self.processed_samples[sample_name]['snp'] = crick_sample
            #             #Both contain the SNP, combine records, TODO: check if combination is warranted!
            #             #Both watson and crick should contain the same type of non reference allele
            #             elif sum(watson_sample.site.INFO['AC']) > 0 and sum(crick_sample.site.INFO['AC']) > 0:
            #                 combined_record = combine_record_samples(watson_sample, crick_sample)
            #                 self.processed_samples[sample_name]['snp'] = combined_record
            #         elif 'AC' in watson_sample.site.INFO:
            #             #Only watson contains the SNP,
            #             #TODO: assert if valid. Only true if SNP is A/G!
            #             self.processed_samples[sample_name]['snp'] = watson_sample
            #         elif 'AC' in crick_sample.site.INFO:
            #             #Only watson contains the SNP,
            #             #TODO: assert if valid. Only true if SNP is A/C!
            #             self.processed_samples[sample_name]['snp'] = crick_sample
            #     elif alt_watson == 'T' and alt_crick == 'T':
            #         combined_record = combine_record_samples(watson_sample, crick_sample)
            #         self.processed_samples[sample_name]['snp'] = combined_record
            #     elif alt_crick == 'C' and alt_watson in 'CT':
            #         self.processed_samples[sample_name]['snp'] = crick_sample
            #         #self.processed_samples[sample_name]['methylated'] = watson_sample
            #     elif alt_watson == 'G' and alt_crick in 'GA':
            #         self.processed_samples[sample_name]['snp'] = watson_sample
            #         #self.processed_samples[sample_name]['methylated'] = crick_sample
            #     elif 'T' in alt_watson + alt_crick and ref_base in alt_watson + alt_crick:
            #         #one of the samples
            #         combined_record = combine_record_samples(watson_sample, crick_sample)
            #         self.processed_samples[sample_name]['snp'] = combined_record
        return 1

    def combine_snp_record(self,watson_record, crick_record, convert_dict):
        """Combined SNP calls into one record taking into account expected bisulfite conversions"""
        if not watson_record.called:
            return None
            #TODO: make SNP calling rules for when watson or SNP record is uncalled
            depth = crick_record.data.DP
            ref_observations = crick_record.data.RO
            alt_observations = crick_record.data.AO
        elif not crick_record.called:
            #TODO: make SNP calling rules for when watson or SNP record is uncalled
            return None
            depth = watson_record.data.DP
            ref_observations = watson_record.data.RO
            alt_observations = watson_record.data.AO
        else:
            #determine on which sample we should base output record
            depth = watson_record.data.DP + crick_record.data.DP
            ref_observations = 0
            ref_base = watson_record.site.REF
            nt_counts = {'C': 0, 'T': 0, 'G': 0, 'A': 0}
            alt_records_watson, alt_records_crick = ([], [])
            if type(watson_record.data.AO) == type([]):
                alt_records_watson += [str(r) for c, r in zip(watson_record.data.AO,watson_record.site.ALT) if c > 0]
            if type(crick_record.data.AO) == type([]):
                alt_records_crick += [str(r) for c, r in zip(crick_record.data.AO,crick_record.site.ALT) if c > 0]

            # account reference base observations for C
            if ref_base == 'C':
                nt_counts['C'] += crick_record.data.RO
                # Methylated C's are evidence of C. This will create an inbalance in the allele count as T's in watson cannot be taken into accoun
                nt_counts['C'] += watson_record.data.RO
                try:
                    crick_index_T = [str(r) for r in crick_record.site.ALT].index('T')
                except ValueError:
                    crick_index_T = None
                # We can only add all C and T watson observations if there is no evidence of a C/T SNP in Crick
                if crick_index_T == None or crick_record.data.AO[crick_index_T] / float(crick_record.data.DP) < 0.05:
                    # all T and C counts for the watson allele are stored as
                    # C observations as no evidence of a T alt allele is present on Crick
                    try:
                        alt_index = [str(r) for r in watson_record.site.ALT].index('T')
                        nt_counts['C'] += watson_record.data.AO[alt_index]
                    except ValueError:
                        pass
            if ref_base == 'G':
                #Add watson record reference observations as these are never disputed.
                nt_counts['G'] += watson_record.data.RO
                nt_counts['G'] += crick_record.data.RO
                try:
                    watson_index_A = [str(r) for r in watson_record.site.ALT].index('A')
                except ValueError:
                    watson_index_A = None
                #TODO 1/2: check if we should use implied evidence for absence of SNP to proceed with calling
                #TODO 2/2: converted reference allele. for now, leave intact. Evaluate!
                if watson_index_A == None or watson_record.data.AO[watson_index_A] / float(watson_record.data.DP) < 0.05:
                    #all A counts for the watson allele are stored as G observations
                    try:
                        alt_index = [str(r) for r in crick_record.site.ALT].index('A')
                        nt_counts['G'] += crick_record.data.AO[alt_index]
                    except ValueError:
                        pass
            if ref_base == 'A':
                #Add watson record reference observations as these are never disputed.
                nt_counts['A'] += watson_record.data.RO
                if 'G' not in alt_records_watson:
                    #no evidenve for G presence is available
                    nt_counts['A'] += crick_record.data.RO

            if ref_base == 'T':
                #Add Crick record reference observations as these are never disputed.
                nt_counts['T'] += crick_record.data.RO
                if 'C' not in alt_records_crick:
                    #if there is no C allele called on the crick allele than All T observations are legit
                    nt_counts['T'] += watson_record.data.RO

            alt_records = []
            if type(watson_record.data.AO) == type([]):
                alt_records += alt_records_watson
            if type(crick_record.data.AO) == type([]):
                alt_records += alt_records_crick

            for nt in convert_dict['watson'].keys():
                #only process records that exist in either the watson or crick alt site
                if nt not in alt_records:
                    continue
                try:
                    watson_alt_index = [str(r) for r in watson_record.site.ALT].index(nt)
                except ValueError:
                    watson_alt_index = None
                try:
                    crick_alt_index = [str(r) for r in crick_record.site.ALT].index(nt)
                except ValueError:
                    crick_alt_index = None
                watson_process = convert_dict['watson'][nt]
                crick_process = convert_dict['crick'][nt]
                if nt == 'G' and crick_process == 'NA' and watson_record.data.RO != 0:
                    crick_process = 'NU'
                #NU = non-convert and use. No conversion is needed
                #NA = Non available for combined call
                #CU = convert and use. This means that if we are looking for evidence of C in watson we might
                #actually have to be looking for C and T combined. This is of course dependent on C not being in crick.
                if watson_process == 'NA':
                    if nt in [str(r) for r in crick_record.site.ALT]:
                        nt_counts[nt] += crick_record.data.AO[crick_alt_index]
                    continue
                elif crick_process == 'NA':
                    if nt in [str(r) for r in watson_record.site.ALT]:
                        nt_counts[nt] += watson_record.data.AO[watson_alt_index]
                    continue
                elif watson_process == crick_process and watson_process == 'NU':
                    #both allelic observations can be used without a problem
                    if nt in [str(r) for r in watson_record.site.ALT]:
                        nt_counts[nt] += watson_record.data.AO[watson_alt_index]
                    if nt in [str(r) for r in crick_record.site.ALT]:
                        nt_counts[nt] += crick_record.data.AO[crick_alt_index]
                elif watson_process == 'CU' and crick_process == 'NU':
                    if watson_alt_index != None:
                        #we have an observation in watson which could be a SNP or a methylation polymorphism.
                        #check if alternate allele is in crick.
                        if nt == 'T':
                            if 'C' in crick_record.site.ALT:
                                c_index = [str(r) for r in crick_record.site.ALT].index('C')
                                c_count = crick_record.data.AO[c_index]
                                if c_count / float(crick_record.data.DP) > 0.05:
                                    #we cannot assume that the T observation in watson is not a converted C.
                                    #only call the T in crick
                                    if crick_alt_index:
                                        nt_counts[nt] += crick_record.data.AO[crick_alt_index]
                                    continue
                            #C Allele is not present in crick record, no evidence for a SNP on this position.
                            #Assume that T count is legible T, used counts from both watson and crick here
                            nt_counts[nt] += watson_record.data.AO[watson_alt_index]
                            if crick_alt_index:
                                nt_counts[nt] += crick_record.data.AO[crick_alt_index]
                        elif nt != 'T':
                            print nt
                            raise FloatingPointError("Seen wrong nucleotide when assuming T")
                    else:
                        #no observations for alternate allele in watson as crick is NU proceed normally
                        if crick_alt_index:
                            continue
                            nt_counts[nt] += crick_record.data.AO[crick_alt_index]
                elif watson_process == 'NU' and crick_process == 'CU':
                    if crick_alt_index != None:
                        #we have an observation in crick which could be a SNP or a methylation polymorphism.
                        #check if alternate allele is in watson.
                        if nt == 'A':
                            if 'G' in watson_record.site.ALT:
                                g_index = [str(r) for r in watson_record.site.ALT].index('G')
                                g_count = watson_record.data.AO[g_index]
                                if g_count / float(watson_record.data.DP) > 0.05:
                                    #we cannot assume that the T observation in crick is not a converted C.
                                    #only call the T in watson
                                    if watson_alt_index:
                                        nt_counts[nt] += watson_record.data.AO[watson_alt_index]
                                    continue
                            #C Allele is not present in watson record, no evidence for a SNP on this position.
                            #Assume that T count is legible T, used counts from both crick and watson here
                            nt_counts[nt] += crick_record.data.AO[crick_alt_index]
                            if watson_alt_index:
                                nt_counts[nt] += watson_record.data.AO[watson_alt_index]
                            continue
                        elif nt != 'A':
                            print nt
                            raise FloatingPointError("Seen wrong nucleotide when assuming T")
                    else:
                        if watson_alt_index:
                            #TODO: evaluate decision to not take into account single-strand SNP observations
                            continue
                            nt_counts[nt] += watson_record.data.AO[watson_alt_index]
            nt_out = {}
            DP = float(sum(nt_counts.values()))
            for nt,count in nt_counts.items():
                try:
                    #TODO: set to parsable parameter
                    if count/DP > 0.05 and count > 1:
                        nt_out[nt] = count
                except ZeroDivisionError:
                    continue

        return nt_out

            #TODO: make GT record

            # header = watson_record.site.FORMAT
            # header_list = header.split(':')
            # call_data = vcf.model.make_calldata_tuple(header_list)
            # values = [out_prim.data.GT,depth,out_prim.data.AD,ref_observations,
            #           alt_observations]
            # model = vcf.model._Call(out_prim.site,
            #                     out_prim.sample,
            #                     call_data(*values))
            # return model
        # except TypeError:
        #     return out_prim

    def filter_snps(self):
        # Sets all the "snp" calles to None, TODO: if filter_snps is finished, this can be removed
        for key in self.watson_file.samples:
            self.processed_samples[key]['snp'] = None
        # Sets all the "snp" calles to None

        # Reference, Watson and Crick record REF are both the same.
        ref_base = self.watson_record.REF

        # If one of the records contains a depth of 0: return None.
        if min([self.watson_record.INFO['DP'], self.crick_record.INFO['DP']]) == 0:
            return None
        #Check if genotypes are already called
        genotypes_called = None
        for sample in self.processed_samples:
            methylation_record = self.processed_samples[sample]['methylated']
            if methylation_record:
                if methylation_record.called:
                    genotypes_called = 1
                    break
        if not genotypes_called:
            self.call_genotypes()
        for watson_sample, crick_sample in izip(self.watson_record, self.crick_record):
            # If there is no call for both the watson and crick record sample, continue as we can not determine
            # whether the polymorphism is a SNP polymorphism.
            if not watson_sample.called or not crick_sample.called:
                continue
            sample_name = watson_sample.sample

            # NU = non-convert and use. No conversion is needed
            # NA = Non available for combined call
            # CU = convert and use. Conversion is always T>C for Watson and A>G for Crick
            if ref_base == "C":
                convert_dict = {'watson':{'A':'NU','T':'NA','G':'NU'},
                                'crick' :{'A':'CU','T':'NU','G':'NU'}}
            elif ref_base == "T":
                convert_dict = {'watson':{'A':'NU','C':'NA','G':'NU'},
                                'crick' :{'A':'CU','C':'NU','G':'NU'}}
            elif ref_base == "G":
                convert_dict = {'watson':{'A':'NU','C':'NU','T':'CU'},
                                'crick' :{'A':'NA','C':'NU','T':'NU'}}
            elif ref_base == "A":
                convert_dict = {'watson':{'C':'NU','T':'CU','G':'NU'},
                                'crick' :{'C':'NU','T':'NU','G':'NA'}}
            else:
                #Base is N
                #TODO: check if we cannot do some basecalling here.
                continue
            combined_count = self.combine_snp_record(watson_sample,crick_sample,convert_dict)
            self.processed_samples[sample_name]["snp"] = (combined_count,watson_sample,crick_sample)

        self.call_SNPs()
        return 1

    def call_SNPs(self):
        """Call SNPs considering the observations made for all individuals"""
        #Rules:
        #1. Only alleles that contain a SNP are called. homozygous ref observations only do not count
        #2. VCF record SNP alleles are independent from methylation VCF records in terms of ALT records.
        combined_allele_count = {}
        for sample in self.processed_samples:
            try:
                allele_count = self.processed_samples[sample]["snp"][0]
            except TypeError:
                continue
            #add
            for k, v in allele_count.items():
                try:
                    combined_allele_count[k] += v
                except KeyError:
                    combined_allele_count[k] = v
        if len(combined_allele_count.keys()) == 1 and combined_allele_count.keys()[0] == self.watson_record.REF:
            #all observations are reference allele. No SNP calling output is required here.
            #set all SNP records to None
            for sample in self.processed_samples:
                self.processed_samples[sample]['snp'] = None
        else:
            combined_count_tuple = sorted(combined_allele_count.items(), key=lambda x: x[1])[::-1]
            # get a valid SNP record object.
            site_obj = self.watson_record
            # Add alt alleles in order of appearance
            alt_alleles = []
            for allele in [i[0] for i in combined_count_tuple]:
                if allele != self.watson_record.REF:
                    alt_alleles.append(vcf.model._Substitution(allele))
            site_obj.ALT = alt_alleles
            site_obj.INFO['DP'] = sum(combined_allele_count.values())
            site_obj.INFO['AD'] = [int(i[1]) for i in combined_count_tuple]
            #TODO: check which other INFO objects could be made here
            #Start calling samples with site object
            samples_out = list()
            for sample in sorted(self.processed_samples):
                if not self.processed_samples[sample]["snp"]:
                    empty_model = vcf.model._Call(site_obj,
                    sample, tuple([None]*len(site_obj.FORMAT.split(':'))))
                    samples_out.append(empty_model)
                    continue
                #try:
                allele_count = self.processed_samples[sample]["snp"][0]
                AO = []
                try:
                    # Watson and Crick record references are the same.
                    RO = allele_count[self.watson_record.REF]
                except KeyError:
                    RO = 0
                for nt in [str(i) for i in alt_alleles]:
                    if nt in allele_count:
                        AO.append(allele_count[nt])
                    else:
                        AO.append(0)
                DP = sum(allele_count.values())
                #call GT
                if DP == 0:
                    GT = "./."
                elif RO == DP:
                    GT = "0/0"
                elif RO > max(AO):
                    if len(AO) == 1:
                        GT = "0/1"
                    else:
                        GT = "0/"+str(AO.index(max(AO)))
                elif max(AO) == DP:
                    GT = "1/1"
                else:
                    if len (AO) == 1:
                        GT = "0/1"
                    else:
                        allele_1, allele_2, = heapq.nlargest(2, AO)
                        GT = str(AO.index(allele_1)) + "/" + str(AO.index(allele_2))

                AD = [RO]
                for item in AO:
                    AD.append(item)

                header = ['GT', 'DP', 'AD', 'RO', 'AO']
                call_data = vcf.model.make_calldata_tuple(header)
                values = [GT, DP, AD, RO, AO]
                model = vcf.model._Call(site_obj,
                                    sample,
                                    call_data(*values))
                samples_out.append(model)
                #except TypeError:
                    #print "Print bad genotype"
                    #make empty record
                    #pass

            # Sets a called SNP record as the parent record for a new SNP call object.
            out_sites = dict()

            chr = site_obj.CHROM
            pos = site_obj.POS
            id = site_obj.ID
            ref = site_obj.REF
            alt = site_obj.ALT
            qual = site_obj.QUAL
            filter = site_obj.FILTER
            format = site_obj.FORMAT
            sample_indexes = site_obj._sample_indexes

            processed_record = vcf.model._Record(
                chr, pos, id, ref,
                alt, qual,
                filter, 0, format,
                sample_indexes, samples_out
            )

            if processed_record.ALT:
                self.snp_file.write_record(processed_record)

    def write_records(self):
        """
        Writes the samples that are saved in the self.processed_samples variable.
        """
        # If the methylation or the snp variable in the unique sample is None, the None will
        # be replaced with an empty sample made by the function: make_empty_sample.
        for sample in self.watson_file.samples:
            index = self.crick_file.samples.index(sample)
            if self.processed_samples[sample]['methylated'] == None:
                record_sample = self.watson_record.samples[index]
                empty_sample = make_empty_sample(record_sample)
                self.processed_samples[sample]['methylated'] = empty_sample
            #if self.processed_samples[sample]['snp'] == None:
            #    record_sample = self.watson_record.samples[index]
            #    empty_sample = make_empty_sample(record_sample)
            #    self.processed_samples[sample]['snp'] = empty_sample

        # Converts the dictionary to two separate lists with all the methylation and snp samples in it.
        methylated_samples = []
        snp_samples = []
        for sample in self.watson_file.samples:
            methylated_samples.append(self.processed_samples[sample]['methylated'])
            snp_samples.append(self.processed_samples[sample]['snp'])

        def define_record(self, samples):
            """
            Returns a pyVCF parsable record if a list of samples is given.
            """
            # Sets the record as a parent record depending on reference base.
            if self.watson_record and self.crick_record:

                if self.watson_record.REF == 'C':
                    vcf_file = self.watson_record
                elif self.watson_record.REF == 'G':
                    vcf_file = self.crick_record
                else:
                    vcf_file = self.watson_record
            else:
                vcf_file = self.watson_record

            chr = vcf_file.CHROM
            pos = vcf_file.POS
            id = vcf_file.ID
            ref = vcf_file.REF
            alt = vcf_file.ALT
            qual = vcf_file.QUAL
            filter = vcf_file.FILTER
            format = vcf_file.FORMAT
            sample_indexes = vcf_file._sample_indexes

            processed_record = vcf.model._Record(
                chr, pos, id, ref,
                alt, qual,
                filter, 0, format,
                sample_indexes, samples
            )

            return processed_record

        def define_record_snp(samples):
            """
            Returns a pyVCF parsable record if a list of samples is given.
            """
            # Sets a called SNP record as the parent record for a new SNP call object.
            out_sites = {}
            for sample in samples:
                if not sample.site.ALT[0] in [None,0]:
                    out_sites[len(sample.site.ALT)] = sample.site
            try:
                vcf_file = out_sites[max(out_sites.keys())]
            except ValueError:
                vcf_file = [sample for sample in samples if sample.called][0].site
            chr = vcf_file.CHROM
            pos = vcf_file.POS
            id = vcf_file.ID
            ref = vcf_file.REF
            alt = vcf_file.ALT
            qual = vcf_file.QUAL
            filter = vcf_file.FILTER
            format = vcf_file.FORMAT
            sample_indexes = vcf_file._sample_indexes

            processed_record = vcf.model._Record(
                chr, pos, id, ref,
                alt, qual,
                filter, 0, format,
                sample_indexes, samples
            )

            self.snp_file.write_record(processed_record)
            return processed_record
        #TODO: check SNP calling here
        #if any(sample.called for sample in snp_samples):
            #Use different method here to make sure the number of ALT alleles matches if crick
            #and watson are different. #TODO: determine if this method can not be used always.
            #snp_record = define_record_snp(snp_samples)
            #write_snp_file(self, snp_record)
            #self.snp_file.write_record(snp_record)
            #snp_record_dict is purely used to determine context for CG,CHG and CHH methylation.
            #It does not contain any useful information pertaining to SNPs
            #if self.snp_record_dict.has_key(snp_record.CHROM) == True:
            #    self.snp_record_dict[snp_record.CHROM].append(snp_record)
            #else:
                # self.snp_record_dict.clear()
            #    self.snp_record_dict[snp_record.CHROM] = [snp_record]
        if any(sample.called for sample in methylated_samples):
            methylation_record = define_record(self, methylated_samples)
            self.methylation_file.write_record(methylation_record)
            self.methylated_records.append(methylation_record)

        if self.snp_record_dict.has_key(self.watson_record.CHROM):
            if len(self.snp_record_dict[self.watson_record.CHROM]) >= 5:
                for record in self.methylated_records[:-2]:
                    write_igv_file(self, record)
                del self.methylated_records[:-2]
                del self.snp_record_dict[self.watson_record.CHROM][0]


def zip_tabix(args):
    """Zip and tabix index the output files"""
    cmds = ['bgzip -f %s'% args.methylation_output,
            'tabix -p vcf %s.gz'% args.methylation_output,
            'bgzip -f %s'% args.SNP_output,
            'tabix -p vcf %s.gz'% args.SNP_output
            ]
    for cmd in cmds:
        p = subprocess.Popen(cmd,stdout=subprocess.PIPE,stderr=subprocess.PIPE,shell=True,executable='/bin/bash')
        p.wait()

# If script is called; calls for the main
if __name__ == '__main__':
    main()