#/usr/bin/env pypy
import vcf

input_handle = open('/Volumes/3tb-deena/Scabiosa/output_mapping/Scabiosa.snp.vcf','r')
output_handle = open('/Volumes/3tb-deena/Scabiosa/output_mapping/Scabiosa.filtered2.snp.vcf','w')
template = vcf.Reader(open('/Volumes/3tb-deena/Scabiosa/output_mapping/watson.vcf.gz','r'))
vcf_handle = vcf.Reader(input_handle)
vcf_out = vcf.Writer(output_handle,template,)

for record in vcf_handle:
    if record.num_called > (len(record.samples) * 0.8):
        # min_call = 3
        # if min([sample.data.DP for sample in record.samples if sample.called]) < min_call:
        #     continue
        if record.aaf[0] < 0.3:
            continue
        if max(record.num_hom_alt,record.num_hom_ref) < (0.4 * len(record.samples)):
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