"""analyse coverage of vcf file"""
import vcf
from nested_dict import nested_dict

vcf_file = '/Volumes/Elements/C101HW17070511/raw_data/Scabi_RRBS/output_mapping/snp.vcf.gz'
vcf_handle = vcf.Reader(filename=vcf_file)
called = {}
contigs = {}
#TODO: only include samples with sufficient sequence coverage
for n,record in enumerate(vcf_handle):
    if not n%100000:
        # print n
        print(len(contigs))
    call_count = 0
    for sample in record.samples:
        if sample.called:
            if sample.data.DP >= 4:
                call_count += 1
    try:
        called[call_count] += 1
    except KeyError:
        called[call_count] = 1
    if sample.site.CHROM not in contigs:
        contigs[sample.site.CHROM] = call_count
    else:
        if call_count > contigs[sample.site.CHROM]:
            contigs[sample.site.CHROM] = call_count

for k,v in sorted(called.items()):
    print k,v

# for k,v in sorted(contigs.items()):
#     print k,v