#!/usr/bin/env pypy
import argparse
import os

def parse_args():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(description='Process input files')
    parser.add_argument('-m', '--mergedcalls', type=str, default=None,
                        help='Merged watson and crick calls')
    parser.add_argument('-s', '--SNP_output', type=str, nargs='?', default=None,
                        help='SNP vcf file output name')
    args = parser.parse_args()
    return args


def make_header(args,handle,line):
    """make vcf header for SNP output"""
    #TODO: implement history samtools version / commands to generate VCF calls from watson and crick
    header = '\n'
    #TODO: define header properties
    with open(args.SNP_output,'w') as handle:
        handle.write(header)
    return 0


def combine_counts(observations, watson_ALT, crick_ALT, convert_dict, ref_base):
    """Combine SNP calls into one record taking into account expected bisulfite conversions"""
    #observations is a list of observations separated by ":"
    # Watson_A,Crick_A:Watson_C,Crick_C:Watson_T,Crick_T:Watson_G,Crick_G
    watson_count = {}
    crick_count = {}
    for combined_count,nt in zip(observations,'ACGT'):
        watson_count[nt], crick_count[nt] = [int(s) for s in combined_count.split(',')]
    #prevent spurious low counts to interfere with SNP calling algorithm, discard observations lower than 5%.
    for nt in 'ACGT':
        try:
            if watson_count[nt] / float(sum(watson_count.values())) < 0.05:
                watson_count[nt] = 0
        except ZeroDivisionError:
            pass
        try:
            if crick_count[nt] / float(sum(crick_count.values())) < 0.05:
                crick_count[nt] = 0
        except ZeroDivisionError:
            pass
    # determine on which sample we should base output record
    nt_counts = {'C': 0, 'T': 0, 'G': 0, 'A': 0}
    alt_records_watson, alt_records_crick = ([], [])
    alt_records_watson += [str(nt) for nt, count in watson_count.items() if count != 0 and nt != ref_base]
    alt_records_crick += [str(nt) for nt, count in crick_count.items() if count != 0 and nt != ref_base]
    # account reference base observations for C
    if ref_base == 'C':
        nt_counts['C'] += crick_count['C']
        if alt_records_crick == [] and alt_records_watson == ['T']:
            # fast routine for methylation polymorphisms only
            nt_counts['C'] = crick_count['C'] + sum(watson_count.values())
            return nt_counts
        # We can only add all C and T watson observations if there is no evidence of a C/T SNP in Crick
        if crick_count['T'] > 0:
            if crick_count['T'] / float(sum(crick_count.values())) < 0.05:
                # all T and C counts for the watson allele are stored as
                # C observations as no evidence of a T alt allele is present on Crick
                nt_counts['C'] += watson_count['T']
        # exclude situation in which reference observations are made in one but not the other strand
        if crick_count['C'] > 0 and watson_count['C'] == 0 and watson_count['T'] == 0:
            return {}
    if ref_base == 'G':
        # Add watson record reference observations as these are never disputed.
        nt_counts['G'] += watson_count['G']
        if alt_records_watson == [] and alt_records_crick == ['A']:
            # fast routine for methylation polymorphisms only
            nt_counts['G'] = sum(watson_count.values()) + sum(crick_count.values())
            return nt_counts
        # TODO 1/2: check if we should use implied evidence for absence of SNP to proceed with calling
        # TODO 2/2: converted reference allele. for now, leave intact. Evaluate!
        if watson_count['A'] > 0:
            if watson_count['A']  / float(sum(watson_count.values())) < 0.05:
                # all A counts for the watson allele are stored as G observations
                nt_counts['G'] += crick_count['A']
        # exclude situation in which reference observations are made in one but not the other strand
        if watson_count['G'] > 0 and crick_count['G'] == 0 and crick_count['A'] == 0:
            return {}
    if ref_base == 'A':
        # Add watson record reference observations as these are never disputed.
        nt_counts['A'] += watson_count['A']
        # Have very stringent conditions on absence for taking implied alleles into account
        if watson_count['G'] == 0 and crick_count['G'] == 0:
            # no evidenve for G presence is available
            if watson_count['A'] > 0 and crick_count['A'] > 0:
                # ref observations should be present for both watson and crick allele
                nt_counts['A'] += crick_count['A']
            elif min(crick_count['A'], watson_count['A']) == 0 and \
                            max(crick_count['A'], watson_count['A']) > 0:
                return {}
    if ref_base == 'T':
        # Add Crick record reference observations as these are never disputed.
        nt_counts['T'] += crick_count['T']
        # Have very stringent conditions on absence for taking implied alleles into account
        if 'C' not in alt_records_crick and 'C' not in watson_ALT:
            # if there is no C allele called on the crick allele than All T observations are legit
            if crick_count['T'] > 0 and watson_count['T'] > 0:
                nt_counts['T'] += watson_count['T']
            elif min(crick_count['T'], watson_count['T']) == 0 and \
                            max(crick_count['T'], watson_count['T']) > 0:
                return {}

    alt_records = []
    if 'A' in alt_records_crick:
        # if A in alt_records_crick A or G must be in alt_records_watson
        if 'G' in alt_records_watson and 'A' not in alt_records_watson:
            # DO not add A as alt observation as it is a converted G
            pass
        elif 'A' in alt_records_watson or ref_base in 'GA':
            # TODO: check exception to this rule. In case of watson G/T and Crick A/T with ref G
            # A does not need to be in alt_records.
            alt_records.append('A')
        else:
            return {}
    elif 'A' in alt_records_watson:
        # A was not seen in alt_record watson whereas it should have been. invalidate call
        return {}
    if 'T' in alt_records_watson:
        if 'C' in alt_records_crick and 'T' not in alt_records_crick:
            # DO not add T as alt observation as it can only be a converted C
            pass
        # if T in alt_records_watson C or T must be in alt_record_crick
        elif 'T' in alt_records_crick or ref_base in 'CT':
            alt_records.append('T')
        else:
            return {}
    elif 'T' in alt_records_crick and 'T' not in alt_records_watson:
        return {}
    if 'C' in alt_records_crick:
        # if C in alt_records_crick, T or C should be in REF or ALT records watson
        if 'T' in alt_records_watson or 'C' in alt_records_watson or ref_base in 'CT':
            alt_records.append('C')
        else:
            return {}
    elif 'C' in alt_records_watson:
        return {}
        # C was not seen in alt_record_crick whereas it should have been. invalidate call
    if 'G' in alt_records_watson:
        # if G in alt_records_watson, G or A should be in REF or ALT records watson
        if 'A' in alt_records_crick or 'G' in alt_records_crick or ref_base in 'GA':
            alt_records.append('G')
        else:
            return {}
    elif 'G' in alt_records_crick:
        # G was not seen in alt_record_crick whereas it should have been. invalidate call
        return {}

    # only non conflicting alleles are present in alt_records

    for nt in convert_dict['watson'].keys():
        # only process records that exist in either the watson or crick alt site
        if nt not in alt_records:
            continue
        watson_process = convert_dict['watson'][nt]
        crick_process = convert_dict['crick'][nt]
        if watson_process == 'NA' and crick_process == 'NU':
            if crick_count[nt] == 0:
                # allele is not found in crick
                continue
            nt_counts[nt] += crick_count[nt]
            continue
        elif watson_process == 'NU' and crick_process == 'NA':
            if watson_count[nt] == 0:
                # allele is not found in watson
                continue
            nt_counts[nt] += watson_count[nt]
            continue
        elif watson_process.startswith('NO'):
            no, nt_not, strand = watson_process.split('_')
            assert strand == 'crick'
            if nt_not not in alt_records_crick:
                nt_counts[nt] += watson_count[nt]
                nt_counts[nt] += crick_count[nt]
                continue
        elif crick_process.startswith('NO'):
            no, nt_not, strand = crick_process.split('_')
            assert strand == 'watson'
            if nt_not not in alt_records_watson:
                nt_counts[nt] += watson_count[nt]
                nt_counts[nt] += crick_count[nt]
                continue
        elif crick_process.startswith('ADD') and watson_process == 'NU':
            add, nt_search, no, nt_not, strand = crick_process.split('_')
            assert strand == 'watson'
            # 1. Nucleotide which is supposed to be absent from watson alt records is indeed absend
            # 2. Nucleotide which we are searching for is present in crick alt records
            nt_counts[nt] += watson_count[nt]
            if nt_not not in alt_records_watson:
                if crick_count[nt] != 0:
                    nt_counts[nt] += crick_count[nt]
                    nt_counts[nt] += crick_count[nt_search]
            continue

        elif watson_process.startswith('ADD') and crick_process == 'NU':
            add, nt_search, no, nt_not, strand = watson_process.split('_')
            assert strand == 'crick'
            # 1. Nucleotide which is supposed to be absent from crick alt records is indeed absend
            # 2. Nucleotide which we are searching for is present in watson alt records
            nt_counts[nt] += crick_count[nt]
            if nt_not not in alt_records_crick:
                nt_counts[nt] += watson_count[nt]
                nt_counts[nt] += watson_count[nt_search]
            continue
        else:
            print ''
    nt_out = {}
    DP = float(sum(nt_counts.values()))
    for nt, count in nt_counts.items():
        try:
            # TODO: set to parsable parameter
            if count / DP > 0.05 and count > 0:
                nt_out[nt] = count
            else:
                nt_out[nt] = 0
        except ZeroDivisionError:
            continue
    return nt_out


def call_SNP(line):
    """main SNP calling algorithm"""
    split_line = line.rstrip('\n').split('\t')
    chrom, pos, ref_base, watson_ALT, crick_ALT = split_line[:5]
    if len(ref_base) > 1:
        return None
    if ref_base.upper() == 'N':
        return None
    if ref_base.upper() == "C":
        if crick_ALT == '' and watson_ALT == 'T':
            # Only a methylation polymorphism
            return None
        convert_dict = {'watson': {'A': 'NU', 'T': 'NA', 'G': 'NU'},
                        'crick': {'A': 'NO_G_watson', 'T': 'NU', 'G': 'ADD_A_NO_A_watson'}}
    elif ref_base.upper() == "T":
        convert_dict = {'watson': {'A': 'NU', 'C': 'NA', 'G': 'NU'},
                        'crick': {'A': 'NO_G_watson', 'C': 'NU', 'G': 'ADD_A_NO_A_watson'}}
    elif ref_base.upper() == "G":
        if crick_ALT == 'A' and watson_ALT == '':
            #Only a methylation polymorphism
            return None
        convert_dict = {'watson': {'A': 'NU', 'C': 'ADD_T_NO_T_crick', 'T': 'NO_C_crick'},
                        'crick': {'A': 'NA', 'C': 'NU', 'T': 'NU'}}
    elif ref_base.upper() == "A":
        convert_dict = {'watson': {'C': 'ADD_T_NO_T_crick', 'T': 'NO_C_crick', 'G': 'NU'},
                        'crick': {'C': 'NU', 'T': 'NU', 'G': 'NA'}}
    if watson_ALT == '' and crick_ALT == '':
        return None
    else:
        counts = []
        DP = 0
        ALT = {}
        for observations in split_line[5:]:
            observations = observations.split(':')
            count = combine_counts(observations, watson_ALT, crick_ALT , convert_dict, ref_base)
            for nt,value in count.items():
                if value == 0 or nt == ref_base:
                    continue
                try:
                    ALT[nt] += value
                except KeyError:
                    ALT[nt] = value
            DP += sum(count.values())
            counts.append(count)
        ALT = [s[0] for s in sorted(ALT.items(), key=lambda x: x[1])[::-1]]
        if ALT != []:
            record = make_vcf_record(chrom, pos, ref_base, DP, ALT, counts, split_line[5:])
            return record
        else:
            return None

def get_GT(count,ALT,ref_base):
    """return numeric genotype code given counts of alleles"""
    try:
        nt_count = {ref_base:count[ref_base]}
    except KeyError:
        return './.'
    for nt in ALT:
        nt_count[nt] = count[nt]
    #TODO: determine what to do if allele has 3 nucleotides.
    order = sorted(nt_count.items(), key=lambda x: x[1])[::-1][:2]
    gt = []
    for (nt,count) in order:
        if count / float(sum(nt_count.values())) < 0.05:
            continue
        if nt == ref_base:
            gt.append('0')
        else:
            gt.append('%s' % (ALT.index(nt) + 1))
    if len(gt) == 1:
        gt = gt * 2
    return "/".join(gt)

def make_vcf_record(chrom, pos, ref_base, DP, ALT, counts, observations):
    """Call SNPs considering the observations made for all individuals, make vcf record without"""
    # Rules:
    # 1. Only alleles that contain a SNP are called. homozygous ref observations only do not count
    # 2. VCF record SNP alleles are independent from methylation VCF records in terms of ALT records.
    #TODO: add strand evidence specification
    #TODO: add num-het, num-hom, number called, other custom fields
    output = [chrom, pos, '.', ref_base, ','.join(ALT), '.', '.','DP=%s'%DP]
    RO = sum([c[ref_base] for c in counts if c != {}])
    output[-1] += ';RO=%s' % RO
    output[-1] += ';AO=%s' % (DP - RO)
    output.append('GT:DP:ADW:ADC:RO:AO')
    ADW = 0
    ADC = 0
    GT_LIST = {'uncalled':0,'HET':0,'HOM_REF':0,'HOM_ALT':0}
    ALT_FREQ = []
    for count,observation in zip(counts,observations):
        values = {'DP':sum(count.values())}
        values['ADW'] = ','.join([i.split(',')[0] for i in observation.split(':')])
        values['ADC'] = ','.join([i.split(',')[1] for i in observation.split(':')])
        values['GT'] = get_GT(count, ALT, ref_base)
        if values['GT'][0] == values['GT'][2]:
            if values['GT'][0] == '.':
                GT_LIST['uncalled'] += 1
            elif values['GT'][0] == '0':
                GT_LIST['HOM_REF'] += 1
            else:
                GT_LIST['HOM_ALT'] += 1
        else:
            GT_LIST['HET'] += 1
            ALT_FREQ.append((sum(count.values()) - count[ref_base]) / float(sum(count.values())))
        try:
            values['RO'] = '%s' % count[ref_base]
            for nt in ALT:
                if count[nt] != 0:
                    append_value = str(count[nt])
                else:
                    append_value = '0'
                try:
                    values['AO'].append(append_value)
                except KeyError:
                    values['AO'] = [append_value]

        except KeyError:
            values['RO'] = '.'
            values['AO'] = '.'
        values['AO'] = ','.join(values['AO'])
        output.append('%(GT)s:%(DP)s:%(ADW)s:%(ADC)s:%(RO)s:%(AO)s'%values)
    try:
        output[7] += ';AB=%.2f' % (sum(ALT_FREQ) / float(len([v for v in ALT_FREQ if v != 0.0])))
    except ZeroDivisionError:
        output[7] += ';AB=0'
    for k,v in GT_LIST.items():
        output[7] += ';%s=%s'% (k,v)
    #AC gets translated into allele count
    #AF gets translated into allele frequency
    #AN gets translated into total # of alleles.
    vcf_line = '\t'.join(output) + '\n'
    return vcf_line


def main():
    """main function loop"""
    args = parse_args()
    # make_header(args,handle,line)
    count = 0
    with open(args.SNP_output, 'w') as handle:
        if args.mergedcalls.endswith('.gz'):
            agent = 'pigz -cd '
        else:
            agent = 'cat '
        for line in os.popen(agent + " %s" % args.mergedcalls):
            if not count:
                make_header(args,handle,line)
            snp_record = call_SNP(line)
            count += 1
            if not count % 1000000:
                print 'processed %s lines ' % count
            if snp_record:
                handle.write(snp_record)
    handle.close()
    os.popen('bgzip -f %s' % args.SNP_output)
    os.popen('tabix -p vcf %s.gz' % args.SNP_output)


if __name__ == '__main__':
    main()