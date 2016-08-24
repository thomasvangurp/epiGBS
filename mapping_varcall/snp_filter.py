#/usr/bin/env pypy
import vcf
import argparse

def parse_args():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(description='Process input files')
    parser.add_argument('-i', '--input', type=str, default=None,
                        help='SNP input VCF file')
    parser.add_argument('-o', '--output', type=str, default=None,
                        help='filtered SNP output file')
    parser.add_argument('-t', '--template', type=str, default=None,
                        help='template to take VCF header from')
    args = parser.parse_args()
    return args



def filter(args):
    """filter vcf output records"""
    input_handle = open(args.input, 'r')
    output_handle = open(args.output, 'w')
    template = vcf.Reader(args.template, 'r')
    vcf_handle = vcf.Reader(input_handle)
    vcf_out = vcf.Writer(output_handle, template, )
    for record in vcf_handle:
        if record.num_called > (len(record.samples) * 0.6):
            # min_call = 3
            # if min([sample.data.DP for sample in record.samples if sample.called]) < min_call:
            #     continue
            if float(record.INFO['AB'][0]) < 0.3 or float(record.INFO['AB'][0]) > 0.8:
                continue
            if int(record.INFO['HET'][0]) > 0.9 * record.num_called or int(record.INFO['HET'][0]) < 0.1 * record.num_called:
                continue
            # if max(record.num_hom_alt,record.num_hom_ref) < (0.4 * len(record.samples)):
            #     continue
            ht_rates = []
            # for sample in record.samples:
            #     if sample.is_het:
            #         try:
            #             ht_rates.append(min(sum(sample.data.AO),int(sample.data.RO)))
            #         except TypeError:
            #             ht_rates.append(min(sample.data.AO, int(sample.data.RO)))
            # if ht_rates != []:
            #     if min(ht_rates) < 2:
            #         continue
            # if record.REF == 'C':
            #     if 'T' in record.alleles:
            #         continue
            # elif record.REF == 'T':
            #     if 'C' in record.alleles:
            #         continue
            # elif record.REF == 'G':
            #     if 'A' in record.alleles:
            #         continue
            # elif record.REF == 'A':
            #     if 'G' in record.alleles:
            #         continue
            vcf_out.write_record(record)


def main():
    """Main function loop"""
    args = parse_args()
    filter(args)


if __name__ == '__main__':
    main()