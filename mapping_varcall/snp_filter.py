#/usr/bin/env pypy
import vcf

input_handle = open('/Users/thomasvangurp/epiGBS/Zwitserland/pilot_60/seq5U7ms_/Gal_mol/output_mapping/pypy_samtools_full_snp.vcf','r')
output_handle = open('/Users/thomasvangurp/epiGBS/Zwitserland/pilot_60/seq5U7ms_/Gal_mol/output_mapping/snp.selection.vcf','w')
template = vcf.Reader(open('/Users/thomasvangurp/epiGBS/Zwitserland/pilot_60/seq5U7ms_/Gal_mol/output_mapping/watson.vcf.gz','r'))
vcf_handle = vcf.Reader(input_handle)
vcf_out = vcf.Writer(output_handle,template,)

for record in vcf_handle:
    if record.num_called > (len(record.samples) * 0.7):
        min_call = 3
        if min([sample.data.DP for sample in record.samples if sample.called]) < min_call:
            continue
        if record.aaf[0] < 0.1:
            continue
        if max(record.num_hom_alt,record.num_hom_ref) < 2:
            continue
        ht_rates = []
        for sample in record.samples:
            if sample.is_het:
                ht_rates.append(min(sum(sample.data.AD[1:]),int(sample.data.RO)))
        if ht_rates != []:
            if min(ht_rates) < 2:
                continue
        if sorted([str(s) for s in record.alleles]) in (['C','T'],['A','G']):
            continue
        vcf_out.write_record(record)