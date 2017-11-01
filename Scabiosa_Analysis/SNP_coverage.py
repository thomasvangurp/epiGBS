"""analyse coverage of vcf file"""
import vcf
from nested_dict import nested_dict

vcf_file = '/Volumes/Elements/C101HW17070511/raw_data/Scabi_RRBS/output_mapping/snp.vcf.gz'
vcf_handle = vcf.Reader(filename=vcf_file)
called = {}
contigs = {}
for n,record in enumerate(vcf_handle):
    if not n%1000:
        print n
    call_count = 0
    for sample in record.samples:
        if sample.called:
            if sample.data.DP > 10:
                call_count += 1
    try:
        called[call_count] += 1
    except KeyError:
        called[call_count] = 1

for k,v in sorted(called.items()):
    print k,v