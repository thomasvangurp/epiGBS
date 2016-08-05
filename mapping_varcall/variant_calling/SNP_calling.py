#!/usr/bin/env pypy
import argparse
import gzip

def parse_args():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(description='Process input files')
    parser.add_argument('-r', '--reference', type=str, nargs='?', default=None,
                        help='reference genome input.')
    parser.add_argument('-m', '--mergedcalls', type=str, default=None,
                        help='Merged watson and crick calls')
    parser.add_argument('-s', '--SNP_output', type=str, nargs='?', default=None,
                        help='SNP vcf file output name')
    args = parser.parse_args()
    return args


def make_header(args):
    """make vcf header for SNP output"""
    header = '\n'
    with open(args.SNP_output) as handle:
        handle.write(header)
    return 0


def combine_counts(observations, watson_ALT, crick_ALT, convert_dict, ref_base):
    """Combine SNP calls into one record taking into account expected bisulfite conversions"""
    #observations is a list of observations separated by ":"
    # Watson_A,Crick_A:Watson_C,Crick_C:Watson_T,Crick_T:Watson_G,Crick_G
    watson_count = {}
    crick_count = {}
    for combined_count,nt in zip(observations,'ACGT'):
        watson_count[nt], crick_count[nt] = combined_count.split(',') 
    
    # determine on which sample we should base output record
    nt_counts = {'C': 0, 'T': 0, 'G': 0, 'A': 0}
    alt_records_watson, alt_records_crick = ([], [])
    alt_records_watson += [str(nt) for nt, count in watson_count.items() if count != 0]
    alt_records_crick += [str(nt) for nt, count in crick_count.items() if count != 0]
    # account reference base observations for C
    if ref_base == 'C':
        nt_counts['C'] += crick_count['C']
        if alt_records_crick == [] and alt_records_watson == ['T']:
            # fast routine for methylation polymorphisms only
            nt_counts['C'] = crick_count['C'] + sum(watson_count.values())
            return nt_counts
        # We can only add all C and T watson observations if there is no evidence of a C/T SNP in Crick
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
            nt_counts['G'] = sum(watson_count.values(),crick_count.values())
            return nt_counts
        # TODO 1/2: check if we should use implied evidence for absence of SNP to proceed with calling
        # TODO 2/2: converted reference allele. for now, leave intact. Evaluate!
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
    elif 'T' in alt_records_crick:
        # T was not seen in alt_record watson whereas it should have been. invalidate call
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
            if watson_count[nt] == 0:
                # allele is not found in crick
                continue
            nt_counts[nt] += crick_count[nt]
            continue
        elif watson_process == 'NU' and crick_process == 'NA':
            if watson_alt_index == None:
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
        except ZeroDivisionError:
            continue
    return nt_out


def call_SNP(line, ref_base):
    """main SNP calling algorithm"""
    split_line = line.rstrip('\n').split('\t')
    chrom, pos, ref_base, watson_ALT, crick_ALT = split_line[:5]
    if ref_base == "C":
        convert_dict = {'watson': {'A': 'NU', 'T': 'NA', 'G': 'NU'},
                        'crick': {'A': 'NO_G_watson', 'T': 'NU', 'G': 'ADD_A_NO_A_watson'}}
    elif ref_base == "T":
        convert_dict = {'watson': {'A': 'NU', 'C': 'NA', 'G': 'NU'},
                        'crick': {'A': 'NO_G_watson', 'C': 'NU', 'G': 'ADD_A_NO_A_watson'}}
    elif ref_base == "G":
        convert_dict = {'watson': {'A': 'NU', 'C': 'ADD_T_NO_T_crick', 'T': 'NO_C_crick'},
                        'crick': {'A': 'NA', 'C': 'NU', 'T': 'NU'}}
    elif ref_base == "A":
        convert_dict = {'watson': {'C': 'ADD_T_NO_T_crick', 'T': 'NO_C_crick', 'G': 'NU'},
                        'crick': {'C': 'NU', 'T': 'NU', 'G': 'NA'}}
    if watson_ALT == '' and crick_ALT == '':
        return None
    else:
        for observations in split_line[6:]:
            counts = combine_counts(observations, watson_ALT, crick_ALT , convert_dict, ref_base)



def main():
    """main function loop"""
    args = parse_args()
    make_header(args)
    with open(args.SNP_output) as handle:
        for line in gzip.open(args.mergedcalls):
            snp_record = call_SNP(line)
            if snp_record:
                handle.write(snp_record)


if __name__ == '__main__':
    main()