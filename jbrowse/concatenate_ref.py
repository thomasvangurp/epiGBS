"""concatenate reference sequence"""
import argparse
import os
import pysam


def parse_args():
    """Pass command line arguments"""

    parser = argparse.ArgumentParser(description='Concatenate reference for display in jbrowse')
    parser.add_argument('-o', '--output_dir', help='output directory')
    parser.add_argument('-i', '--input_dir', help='input directory with bam files')
    parser.add_argument('-r', '--ref', help='reference sequence in')
    args = parser.parse_args()
    if 'input_dir' in args:
        if os.path.exists(os.path.join(args.input_dir, 'watson.dedup.bam')) and args.watson is None:
            args.watson = os.path.join(args.input_dir, 'watson.dedup.bam')
        if os.path.exists(os.path.join(args.input_dir, 'crick.dedup.bam')) and args.crick is None:
            args.crick = os.path.join(args.input_dir, 'crick.dedup.bam')
    if not os.path.exists(args.output_dir):
        os.mkdir(args.output_dir)
    return args

def coverage_ref(args):
    """Get coverage of reference sequence"""

def parse_bam(bam):
    """parse bam file to get mapping per contig and individual"""
    mapping_dict = {}
    file_handle = pysam.AlignmentFile(bam)
    for read in file_handle:
        print ''




def coverage_ref(args):
    """Get coverage of reference sequence"""
    #parse watson
    mapping_watson = parse_bam(args.watson)
    #parse crick
    mapping_crick = parse_bam(args.crick)

def rewrite_bam(args):
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
        sample = '_'.join(dict(read.tags)['RG'].split('_')[2:])
        handle = out_handles[sample]['watson']
        handle.write(read)
    for subdict in out_handles.values():
        subdict['watson'].close()
    i = 0
    for read in bam_crick_handle:
        i += 1
        if not i % 100000:
            print 'processed %s reads' % i
        sample = '_'.join(dict(read.tags)['RG'].split('_')[2:])
        handle = out_handles[sample]['crick']
        handle.write(read)
    for subdict in out_handles.values():
        subdict['crick'].close()
    for item in bam_watson_handle.header['RG']:
        watson_path = os.path.join(args.output_dir, '%s.watson.bam' % item['SM'])
        crick_path = os.path.join(args.output_dir, '%s.crick.bam' % item['SM'])
        pysam.index(watson_path)
        pysam.index(crick_path)


def main():
    """main function
    :rtype: int
    """
    args = parse_args()
    rewrite_bam(args)
    return 0


if __name__ == '__main__':
    return_code = main()