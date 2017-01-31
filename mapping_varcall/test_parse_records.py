import gzip
watson_file = gzip.open("/Volumes/3tb-deena/Daphnia/output_mapping/Daphnia.watson.vcf.gz", 'r')

count = 0
for n,line in enumerate(watson_file):
    split_line = line.split('\t')
    if not n%1000:
        print n