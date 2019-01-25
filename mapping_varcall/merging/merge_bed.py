__author__ = 'thomasvangurp'
"""merge two bed filed"""
import sys
import gzip

f1,f2,output = sys.argv[1:]

def merge_lines(l1,l2,to_add):
    """merge 2 lines"""
    merged_line = l1 + l2[4:]
    merged = 0
    for item in to_add:
        value_to_add = merged_line[item[0] + len(l2[4:]) - merged]
        merged_line.pop(item[0] + len(l2[4:]) - merged)
        merged += 1
        if value_to_add == 'None':
            continue
        elif merged_line[item[1]] == 'None':
            merged_line[item[1]] = str(value_to_add)
        else:
            merged_line[item[1]] = str(int(merged_line[item[1]]) + int(value_to_add))
    merged_line[3] = str(int((len(merged_line[4:]) - merged_line.count('None'))/2))
    return '\t'.join(merged_line) + '\n'

output_handle = open(output,'wt')
with open(f1) as f1_handle:
    with open(f2) as f2_handle:
        while True:
            header_1 = f1_handle.readline()[:-1].split('\t')
            header_2 = f2_handle.readline()[:-1].split('\t')
            to_add = [(header_2.index(p),header_1.index(p),) for p in header_2[4:] if p in header_1]
            merged_header = header_1 + header_2[4:]
            popped = 0
            for item in to_add:
                assert merged_header[item[0] + len(header_1[4:]) - popped] == merged_header[item[1]]
                merged_header.pop(item[0] + len(header_1[4:]) - popped)
                popped += 1
            output_handle.write('\t'.join(merged_header) + '\n')
            break
        while True:
            try:
                line_1 = f1_handle.readline()[:-1].split('\t')
                line_2 = f2_handle.readline()[:-1].split('\t')
                while True:
                    chrom1, pos1 = line_1[:2]
                    chrom2, pos2 = line_2[:2]
                    if chrom1 == chrom2 and pos1 == pos2:
                        break
                    elif chrom1 == chrom2:
                        if pos1 > pos2:
                            line_2 = f2_handle.readline()[:-1].split('\t')
                        elif pos1 < pos2:
                            line_1 = f1_handle.readline()[:-1].split('\t')
                    elif chrom1 > chrom2:
                        line_2 = f2_handle.readline()[:-1].split('\t')
                    else:
                        line_1 = f1_handle.readline()[:-1].split('\t')
            except IndexError:
                break
            merged_line  = merge_lines(line_1,line_2,to_add)
            output_handle.write(merged_line)
output_handle.close()