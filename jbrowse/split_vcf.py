import vcf
import argparse
import os
"""split watson and crick bam file into sample specific bam files"""


def parse_args():
    """Pass command line arguments"""

    parser = argparse.ArgumentParser(description='use bwameth for mapping reads')
    parser.add_argument('-i', '--input_dir', help='input directory')
    parser.add_argument('-o', '--output_dir', help='output directory')
    parser.add_argument('-w', '--watson', help='watson vcf file')
    parser.add_argument('-c', '--crick', help='crick vcf file')
    args = parser.parse_args()
    if 'input_dir' in args:
        if os.path.exists(os.path.join(args.input_dir, 'watson.bam')) and args.watson is None:
            args.watson = os.path.join(args.input_dir, 'watson.bam')
        if os.path.exists(os.path.join(args.input_dir, 'crick.bam')) and args.crick is None:
            args.crick = os.path.join(args.input_dir, 'crick.bam')
    if not os.path.exists(args.output_dir):
        os.mkdir(args.output_dir)
    return args


def split_bam(args):
    """split bam file using pysam"""
    bam_watson_handle = pysam.AlignmentFile(args.watson)
    bam_crick_handle = pysam.AlignmentFile(args.crick)
    out_handles = {}
    for item in bam_watson_handle.header['RG']:
        watson_path = os.path.join(args.output_dir, '%s.watson.bam' % item['SM'])
        crick_path = os.path.join(args.output_dir, '%s.crick.bam' % item['SM'])
        watson_handle = pysam.AlignmentFile(watson_path, "wb", template=bam_watson_handle)
        crick_handle = pysam.AlignmentFile(crick_path, "wb", template=bam_crick_handle)
        out_handles[item['SM']] = {'watson': watson_handle, 'crick': crick_handle}
    i = 0
    for read in bam_watson_handle:
        i += 1
        if not i % 100000:
            print 'processed %s reads' % i
        sample = dict(read.tags)['RG'].split('_')[-1]
        handle = out_handles[sample]['watson']
        handle.write(read)
    i = 0
    for read in bam_crick_handle:
        i += 1
        if not i % 100000:
            print 'processed %s reads' % i
        sample = dict(read.tags)['RG'].split('_')[-1]
        handle = out_handles[sample]['crick']
        handle.write(read)
    for subdict in out_handles.values():
        subdict['watson'].close()
        pysam.index(subdict['watson'], "%s.bai" % subdict['watson'])
        subdict['crick'].close()
        pysam.index(subdict['crick'], "%s.bai" % subdict['crick'])

def main():
    """main function"""
    args = parse_args()
    split_bam(args)
    return 0


if __name__ == '__main__':
    main()