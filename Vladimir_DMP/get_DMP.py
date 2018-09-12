#!/usr/bin/env pypy3
import argparse
import sys
from openpyxl import load_workbook
import os

def parse_args():
    """Pass command line arguments"""
    if not sys.argv[1:]:
        sys.argv.append('-h')
    parser = argparse.ArgumentParser(description='parse command line options for input')
    # input files
    parser.add_argument('-x', '--xlsx',
                        help='xlsx template for input with links to DMP-csv files')
    args = parser.parse_args()
    return args

def get_DMP_dict(wb, context):
    """get a list of DMP locations"""
    ws = wb[context]
    contents = []
    for i, row in enumerate(ws.rows):
        contents.append([])
        for j, cell in enumerate(row):
            contents[-1].append((cell.value, cell.hyperlink,))
    header = contents[0]
    output_dict = {}
    for row in contents[1:]:
        for position, item in enumerate(row):
            DMP_type, link_object = item
            genotype = header[position][0]
            try:
                file_path = os.path.join("/Users/thomasvangurp/Dropbox/", link_object.target.replace('../',''))
            except AttributeError:
                continue
            assert os.path.exists(file_path)
            try:
                output_dict[DMP_type][genotype] = file_path
            except KeyError:
                output_dict[DMP_type] = {genotype: file_path}
    return output_dict

def get_DMP(**kwargs):
    """get list of DMP from csv file given criteria
    for now, diffmeth_p_adj_fdr_max and the **type** can be set.
    When type is set to chrom the chrom is taken as present. When it is set to chrom_pos the
    position within the chrom is set"""
    output = set()
    file_handle = open(kwargs['file'],'r')
    header = file_handle.readline().rstrip('\n').split(',')
    for line in file_handle:
        split_line = line.rstrip('\n').split(',')
        if float(split_line[header.index('diffmeth.p.adj.fdr')]) < 0.05:
            if kwargs['type'] == 'chrom':
                output.add(split_line[1])
            else:
                output.add('%s_%s' % (split_line[1], split_line[2]))
    return output


def make_DMP_list(outputdir, location_dict):
    """makes list of DMPs considering a certain comparison type and a location_dict"""
    # limit_to_type will only produce lists of DMPs for comparisons that contain the keywords in this list
    limit_to_type = "treatment"
    for comparison, link_dict in location_dict.items():
        if limit_to_type in comparison:
            output_subdir = os.path.join(outputdir,comparison)
            if not os.path.exists(output_subdir):
                os.mkdir(output_subdir)
            for gt, file in link_dict.items():
                DMPs = get_DMP(file=file, diffmeth_p_adj_fdr_max = 0.05,type = "pos")
                with open(os.path.join(output_subdir, '%s.dmp.perpos.txt'% gt),'w') as output_handle:
                    output_handle.write('\n'.join(DMPs) + '\n')
    pass

def main():
    """main function loop"""
    wb = load_workbook('Phragmites_DMP.xlsx')
    context = "CG"
    output_dict = get_DMP_dict(wb, context)
    output_dir = "output"
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)
    make_DMP_list(output_dir, output_dict)
    pass


if __name__ == '__main__':
    main()
