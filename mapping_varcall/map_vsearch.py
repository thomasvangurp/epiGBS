"""mape files with vsearch"""
import gzip


#1. split ref in joined and merged portion
#2. map joined reads with C/T and G/A reads
#
# ref = open('/Users/thomasvangurp/epiGBS/Zwitserland/pilot_60/seq5U7ms_/Ver_cha/output_denovo/consensus_cluster.renamed.fa')
#
#
# ref_joined_watson = open('/tmp/joined_ref_watson.fa','w')
# ref_joined_crick = open('/tmp/joined_ref_crick.fa','w')
# ref_merged_watson = open('/tmp/merged_ref_watson.fa','w')
# ref_merged_crick = open('/tmp/merged_ref_crick.fa','w')
#
# seq = ''
# for line in ref:
#     if line.startswith('>'):
#         if seq != '':
#             if 'NNNNNNNN' in seq.upper():
#                 ref_joined_watson.write(header + seq.upper().replace('C','T')+'\n')
#                 ref_joined_crick.write(header + seq.upper().replace('G','A')+'\n')
#             else:
#                 ref_merged_watson.write(header + seq.upper().replace('C', 'T') + '\n')
#                 ref_merged_crick.write(header + seq.upper().replace('G', 'A') + '\n')
#         seq = ''
#         header = line
#     else:
#         seq += line.rstrip('\n')
#
# #make fasta for merged reads
# merged_reads = gzip.open('/Users/thomasvangurp/epiGBS/Zwitserland/pilot_60/seq5U7ms_/Ver_cha/output_denovo/Assembled.trimmed.fq.gz')
# merged_fa_watson = open('/tmp/merged_reads_watson.fa','w')
# merged_fa_crick = open('/tmp/merged_reads_crick.fa','w')
# while True:
#     read = []
#     for i in range(4):
#         try:
#             read.append(merged_reads.next())
#         except StopIteration:
#             break
#     if read == []:
#         break
#     if 'watson' in read[0].lower():
#         convert_seq = read[1].replace('C','T')
#         header = '>%s' % (read[0][1:-1].replace(' ', '|').replace('\t', '|'))
#         header += '|' + read[1]
#         assert header.count('\n') == 1
#         merged_fa_watson.write(header + convert_seq)
#     else:
#         convert_seq = read[1].replace('G', 'A')
#         header = '>%s' % (read[0][1:-1].replace(' ', '|').replace('\t', '|'))
#         header += '|' + read[1]
#         assert header.count('\n') == 1
#         merged_fa_crick.write(header + convert_seq)

#Vsearch cmd:
#vsearch --usearch_global /tmp/merged_reads_watson.fa --db /tmp/merged_ref_watson.fa --samout /tmp/merged_watson.sam --id 0.95


def parse_sam(in_file, out_file):
    """parse sam file and write correct output"""
    out_handle = open(out_file , 'w')
    for line in open(in_file, 'r'):
        if line.startswith('@'):
            out_handle.write(line)
            continue
        split_line = line.split('\t')
        header = split_line[0].split('|')
        out_line = [header[0]] + split_line[1:9] + [header[-1]] + ['*'] + header[1:-1] + split_line[11:]
        out_handle.write('\t'.join(out_line))


parse_sam('/tmp/STAR_index/watson_mergedAligned.out.sam', '/tmp/merged_watson.rewrite.sam')
#samtools view  -F 256 -@ 4 -Shu /tmp/merged_watson.rewrite.sam |sambamba sort -u -m 4GB -t 4 /dev/stdin -o /tmp/sorted
