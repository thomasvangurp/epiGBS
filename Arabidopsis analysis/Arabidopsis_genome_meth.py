__author__ = 'thomasvangurp'
"""Get median coverage and methylations percentage per 1MB of arabidopsis genome"""


genome_dict = {}
bed_file = open('/Users/thomasvangurp/epiGBS/Baseclear/unfiltered_sequences/seqNNAtlE/Athal/ref_mapping/methylation.bed')
header = bed_file.readline().split('\t')
for line in bed_file:
    split_line = line[:-1].split('\t')
    chrom,pos = split_line[:2]
    MB = int(pos)/1000000
    for i in range(5,len(split_line),2):
        methylation = split_line[i-1]
        coverage = split_line[i]
        if coverage != 'None':
            try:
                ratio = int(methylation)/int(coverage)
            except ZeroDivisionError:
                ratio = 0
            try:
                genome_dict[chrom][MB].append((int(coverage),ratio))
            except KeyError:
                if chrom not in genome_dict:
                    genome_dict[chrom] = {MB:[(int(coverage),ratio)]}
                elif MB not in genome_dict[chrom]:
                    genome_dict[chrom][MB] =[(int(coverage),ratio)]
                    # print ''
for chrom,subdict in sorted(genome_dict.items()):
    for MB,coverage in sorted(subdict.items()):
        median_coverage = sorted([item[0] for item in coverage if item[0] > 10])[len([item[0] for item in coverage if item[0] > 10])/2]
        average_coverage = sum([item[0] for item in coverage if item[0] > 10])/len([item[0] for item in coverage if item[0] > 10])
        average_methylation = (sum([item[1] for item in coverage if item[0]>10])/float(len([item[1] for item in coverage if item[0]>10])))
        print chrom,MB,median_coverage,'%.4f'%average_methylation