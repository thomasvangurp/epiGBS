#!/usr/bin/env pypy
# -*- coding: utf-8 -*-
import argparse
import subprocess
import os
import math
import gzip
import tempfile
import urllib
import json
import ssl
import codecs
__author__ = 'thomasvangurp'
__description__ = "map reads orders of magnitudes faster using STAR"


def parse_args():
    """Pass command line arguments"""
    parser = argparse.ArgumentParser(description='use STAR for mapping reads')
    # input files
    parser.add_argument('-s', '--sequences',
                        help='number of sequences to take for testing',default=None)
    parser.add_argument('--tmpdir',
                        help='tmp directory', default="/tmp/")
    parser.add_argument('--input_dir',
                        help='optional: Choose input directory')
    parser.add_argument('--reads_R1',
                        help='Forward unmerged reads')
    parser.add_argument('--reads_R2',
                        help='Reverse unmerged reads')
    parser.add_argument('--merged',
                        help='merged watson and crick fastq')
    parser.add_argument('--reference',
                        help='reference clusters')
    parser.add_argument('--barcodes',
                        help='Barcodes used in output')
    parser.add_argument('--species',
                        help='Species: if selected only that species will be put in BAM RG header')
    parser.add_argument('--threads',
                        help='Number of threads to used where multithreading is possible')
    parser.add_argument('--output_dir',
                        help='Choose output directory')
    parser.add_argument('--extraflags',
                        help='extra flags for testing')
    args = parser.parse_args()
    if args.input_dir:
        args.reads_R1 = os.path.join(args.input_dir, 'all.R1.fq.gz')
        args.reads_R2 = os.path.join(args.input_dir, 'all.R2.fq.gz')
        args.merged = os.path.join(args.input_dir, 'all.merged.fq.gz')
        args.joined = os.path.join(args.input_dir, 'all.joined.fq.gz')
        args.reference = os.path.join(args.input_dir, 'ref.fa')
    if args.output_dir:
        if not os.path.exists(args.output_dir):
            os.mkdir(args.output_dir)
        if 'log' not in args:
            args.log = os.path.join(args.output_dir, 'mapping_variantcalling.log')
    args.tmpdir = tempfile.mkdtemp(suffix='STAR', prefix='tmp', dir=args.tmpdir)
    return args


def get_version():
    """get version of current script"""
    parent_dir = os.path.dirname(os.path.realpath(__file__))
    while True:
        if '.git' in os.listdir(parent_dir):
            break
        parent_dir = os.path.dirname(parent_dir)
    git_log = os.path.join(parent_dir, '.git', 'logs', 'HEAD')
    handle = open(git_log, 'r')
    log_lines = [l.split('\t') for l in handle.readlines()]
    # now get latest github commit
    url = 'https://api.github.com/repos/thomasvangurp/epiGBS/commits'
    context = ssl._create_unverified_context()
    result = json.load(urllib.urlopen(url, context=context))
    print('')


def remove_PCR_duplicates(args):
    """us subprocess to run external remove_PCR_duplicates routine"""
    cmd = "mark_PCR_duplicates.py -i %s" % os.path.join(args.output_dir,'out.bam')
    cmd += " --output %s" % os.path.join(args.output_dir,'out.dedup.bam')
    cmd += " -b %s" % args.barcodes
    cmd += " -r %s" % args.reference
    log = "Removal of PCR duplicates"
    run_subprocess([cmd], args, log)
    return args


def run_subprocess(cmd, args, log_message):
    "Run subprocess under standardized settings"
    # force the cmds to be a string.
    if len(cmd) != 1:
        cmd = [" ".join(cmd)]
    with open(args.log, 'a') as log:
        log.write("now starting:\t%s\n" % log_message)
        log.write('running:\t%s\n' % (' '.join(cmd)))
        log.flush()
        p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True, executable='/bin/bash')
        stdout, stderr = p.communicate()
        stdout = stdout.decode("utf-8").replace('\r', '\n')
        stderr = stderr.decode("utf-8").replace('\r', '\n')
        if stdout:
            log.write('stdout:\n%s\n' % stdout)
        if stderr:
            log.write('stderr:\n%s\n' % stderr)
        return_code = p.poll()
        if return_code:
            raise RuntimeError(stderr)
        log.write('finished:\t%s\n\n' % log_message)
    return 0


def process_reads_merged(args):
    """process reads and make them ready for mapping with STAR"""
    merged = tempfile.NamedTemporaryFile(suffix=".fastq", prefix='merged', dir=args.tmpdir, delete=False)
    args.merged_out = merged.name
    if args.merged.endswith('.gz'):
        cmd = "pigz -cd %s" % args.merged
    else:
        cmd = "cat %s" % args.merged
    if args.sequences != None:
        cmd += '|head -n %s' % (int(args.sequences) * 4)
    log = 'use mawk to convert spaces and tabs to | in merged reads'
    cmd+= """ |mawk '{if (NR%%4==1) gsub(" ","|") gsub(/\t/,"|");print}' > %s""" % (args.merged_out)
    run_subprocess([cmd], args, log)
    return args


def process_reads_joined(args):
    """process reads and make them ready for mapping with STAR"""
    joined_r1 = tempfile.NamedTemporaryFile(suffix=".fastq", prefix='joined_R1_', dir=args.tmpdir,
                                                   delete=False)
    joined_r2 = tempfile.NamedTemporaryFile(suffix=".fastq", prefix='joined_R2_', dir=args.tmpdir,
                                                   delete=False)
    #use mawk to write fastq files, replacing tabs and spaces for |
    args.joined_r1_out = joined_r1.name
    args.joined_r2_out = joined_r2.name
    if args.reads_R1.endswith('.gz'):
        cmd = "pigz -cd %s" % args.reads_R1
    else:
        cmd = "cat %s" % args.reads_R1
    if args.sequences != None:
        cmd += '|head -n %s' % (int(args.sequences) * 4)
    log = 'use mawk to convert spaces and tabs to | in forward reads'
    cmd+= """ |mawk '{if (NR%%4==1) gsub(" ","|") gsub(/\t/,"|");print}' > %s""" % (args.joined_r1_out)
    run_subprocess([cmd], args, log)

    #process reverse reads
    if args.reads_R2.endswith('.gz'):
        cmd = "pigz -cd %s" % args.reads_R2
    else:
        cmd = "cat %s" % args.reads_R2
    if args.sequences != None:
        cmd += '|head -n %s' % (int(args.sequences) * 4)
    log = 'use mawk to convert spaces and tabs to | in reverse reads'
    cmd+= """ |mawk '{if (NR%%4==1) gsub(" ","|") gsub(/\t/,"|");print}' > %s""" % (args.joined_r2_out)
    run_subprocess([cmd], args, log)
    return args


def index_STAR(args):
    """make STAR index for merged and joined reads"""

    # make STAR index folder for merged path
    merged_STAR_index = os.path.join(args.output_dir, 'STAR_merged')
    if not os.path.exists(merged_STAR_index):
        os.mkdir(merged_STAR_index)
    ref_merged = os.path.join(merged_STAR_index, '%s.merged.fa' % args.species)

    # make STAR index folder for joined path
    joined_STAR_index = os.path.join(args.output_dir, 'STAR_joined')
    if not os.path.exists(joined_STAR_index):
        os.mkdir(joined_STAR_index)
    ref_joined = os.path.join(joined_STAR_index, '%s.joined.fa' % args.species)

    # get file handle for input reference file
    try:
        file_handle = open(args.reference, 'r')
    except IOError:
        raise IOError('file %s does not exist' % args.reference)

    # iterate over input lines and write to references
    joined_len = 0
    merged_len = 0
    joined_count = 0
    merged_count = 0
    ref_merged_handle = open(ref_merged, 'w')
    ref_joined_handle = open(ref_joined, 'w')
    seq = ''
    for line in file_handle:
        if line.startswith('>'):
            if seq != '':
                if 'NNNNNNNN' in seq.upper():
                    joined_len += len(seq)
                    joined_count += 1
                    ref_joined_handle.write(header + seq.upper() + '\n')
                else:
                    merged_len += len(seq)
                    merged_count += 1
                    ref_merged_handle.write(header + seq.upper()+ '\n')
            seq = ''
            header = line
        else:
            seq += line.rstrip('\n')
    # write final sequence, this is always merged
    merged_len += len(seq)
    merged_count += 1
    ref_merged_handle.write(header + seq.upper() + '\n')
    # close file handles
    ref_merged_handle.close()
    ref_joined_handle.close()

    # MAKE LIST for indexes to be made
    index_list = [(joined_len, joined_count, joined_STAR_index, ref_joined),
                  (merged_len, merged_count, merged_STAR_index, ref_merged)]
    # calculate parameters for indexing reference for merged and joined reads.
    for (genome_len, no_clusters, genome_dir, ref) in index_list:
        index_cmd = 'STAR --runThreadN %s --runMode genomeGenerate --genomeDir %s' % (args.threads, genome_dir)
        fasta_file = [file for file in os.listdir(genome_dir) if file.endswith('.fa')][0]
        index_cmd += ' --genomeFastaFiles %s' % os.path.join(genome_dir, fasta_file)
        genomeSAindexNbases = min(14, math.log(genome_len, 2) / 2 - 1)
        index_cmd += ' --genomeSAindexNbases %i' % genomeSAindexNbases
        genomeChrBinNbits = min(18, math.log(genome_len / no_clusters, 2))
        index_cmd += ' --genomeChrBinNbits %i' % genomeChrBinNbits
        log = 'making STAR index of %s' % (ref)
        if 'Genome' not in os.listdir(genome_dir):
            run_subprocess([index_cmd], args, log)
    return args


def map_STAR(args):
    """map reads with STAR"""
    for type in ['joined', 'merged']:

        STAR_index_dir = os.path.join(args.output_dir, 'STAR_%s' % (type))
        cmd = "STAR --runThreadN %s --genomeDir %s" % (args.threads, STAR_index_dir)

        if type == 'joined':
            # cmd += " --readFilesIn %(joined)s" % vars(args)
            cmd += " --readFilesIn %(joined_r1_out)s %(joined_r2_out)s" % vars(args)
        else:
            cmd += " --readFilesIn %(merged_out)s" % vars(args)

        cmd += " --outSAMattributes NM MD AS --outSAMtype SAM"
        cmd += " --outFileNamePrefix %s" % (os.path.join(args.output_dir, '%s' % (type)))
        cmd += " --outReadsUnmapped Fastx"  # output of unmapped reads for inspection
        cmd += " --scoreGapATAC -2 --scoreGapNoncan -2"
        if args.joined_r1_out.endswith('.gz'):
            cmd += " --readFilesCommand pigz -cd"
        # outFilterScoreMinOverLread : float: sam as outFilterMatchNmin, but normalized to the read length (sum of mates lengths for paired-end reads)
        # outFilterMatchNminOverLread: float: same as outFilterScoreMin, but normalized to read length (sum of mates lengths for paired-end reads)

        # –outFilterMultimapNmax 1 int: maximum number of loci the read is allowed to map to. Alignments (all of
        # them) will be output only if the read maps to no more loci than this value.
        cmd += " --outFilterMismatchNoverLmax 0.90"
        cmd += "--outFilterMatchNminOverLread 0.9 --scoreGap -4 " \
               " --alignEndsType Extend5pOfReads12" \
               " --alignSoftClipAtReferenceEnds No" \
               " --outSAMorder PairedKeepInputOrder" \
               " --outFilterMultimapNmax 1" \
               " --scoreInsOpen -1" \
            # make sure we have a bam file sorted by name
        if args.sequences:
            cmd += ' --readMapNumber %s' % args.sequences
        if args.extraflags:
            cmd += ' %s' % args.extraflags
        log = "run STAR for %s reads" % (type)
        if not os.path.exists(os.path.join(args.output_dir, '%sAligned.out.sam' % type)):
            run_subprocess([cmd], args, log)
            log = "write final log of STAR to normal log"
            cmd = "cat %s " % os.path.join(args.output_dir, '%s' % (type) + 'Log.final.out')
            run_subprocess([cmd], args, log)
    return args


def parse_sam(in_file, out_file, read_type, strand):
    """parse sam file and write correct output"""
    out_handle = open(out_file, 'a')
    # if strand == 'watson':
    #     nt = ['C']
    # else:
    #     nt = ['G']
    count = 0
    # print 'Warning, only works for forward mapped reads'
    mismatch = 0
    clip_count_total = 0
    for line in open(in_file, 'r'):
        modulo_line_no = count % 2
        # alternates between 0 and 1
        if line.startswith('@'):
            continue
        split_line = line.rstrip('\n').split('\t')
        # skip read pairs with improper flags.
        # TODO: do this filtering in mark_PCR_duplicates or elsewhere with access to pysam.
        if split_line[1] not in ['0', '99', '147']:
            mismatch += 1
            count += 1
            # continue
        char_count = ''
        clip_count = 0
        for char in split_line[5]:
            if not char.isalpha():
                char_count += char
            elif char == 'S':
                clip_count += int(char_count)
            else:
                char_count = ''
        if clip_count > 6:
            clip_count_total += 1
            count += 1
            # continue
        header = split_line[0].split('|')
        out_line = [header[0]]
        out_line += split_line[1:]
        out_line += header[3:6]
        out_handle.write('\t'.join(out_line) + '\n')
        count += 1
    print('%s mismatches out of %s' % (mismatch, count))
    print('%s reads out of  %s soft clipped more than 5' % (clip_count_total, count))


def addRG(in_files, args):
    """make header for output bamfile and split in watson and crick"""
    # define readgroup header lines by combining the following

    """
    -
    read group
    ID*
    Unique read group identifier. The value of the ID field is used in the RG tags of alignment records.
    SM*
    Sample (use pool name where a pool is being sequenced)
    LB
    Library
    DS
    Description
    PU
    Platform unit (e.g. lane for Illumina or slide for SOLiD); should be a full, unambiguous identifier
    PI
    Predicted median insert size (maybe different from the actual median insert size)
    CN
    Name of sequencing center producing the read.
    DT
    Date the run was produced (ISO 8601 date or date/time).
    PL
    Platform/technology used to produce the read."""

    with open(args.barcodes, 'r') as barcodes:
        sam_out = open(in_files['header'], 'a')
        header = barcodes.readline().split('\t')
        for line in barcodes:
            RG = ['@RG']
            split_line = line.split('\t')
            if args.species and 'Species' in header:
                if split_line[(header.index('Species'))] != args.species:
                    continue
            fc = split_line[(header.index('Flowcell'))]
            lane = split_line[(header.index('Lane'))]
            sample = split_line[(header.index('Sample'))]
            RG.append('ID:%s_%s_%s' % (fc, lane, sample))
            RG.append('SM:%s' % (sample))
            RG.append('LB:%s_%s' % (fc, sample))
            RG.append('PL:ILLUMINA\n')
            sam_out.write('\t'.join(RG))
    sam_out.close()
    return in_files


def make_header(args):
    """Make header for watson and crick bam file"""
    header = os.path.join(args.output_dir, 'header.sam')
    args.header = header
    header_handle = open(header, 'w')
    header_handle.write('@HD\tVN:1.4\n')
    joined_sam = open(os.path.join(args.output_dir, 'joinedAligned.out.sam'))
    merged_sam = open(os.path.join(args.output_dir, 'mergedAligned.out.sam'))
    for line in joined_sam:
        if line.startswith('@'):
            if line.startswith('@SQ'):
                header_handle.write(line)
        else:
            break
    for line in merged_sam:
        if line.startswith('@'):
            if line.startswith('@SQ'):
                header_handle.write(line)
            elif not line.startswith('@HD'):
                header_handle.write(line)
        else:
            break
    header_handle.close()
    in_files = {'header': os.path.join(args.output_dir, 'header.sam')}
    addRG(in_files, args)
    return args


def bam_output(args):
    """Generate watson and crick output bam file"""

    merged_sam = os.path.join(args.output_dir, 'mergedAligned.out.sam')
    joined_sam = os.path.join(args.output_dir, 'joinedAligned.out.sam')
    out_sam = tempfile.NamedTemporaryFile(prefix='tmp', suffix='.sam', dir=args.output_dir)
    # rewrite sam file merged and joined for watson and crick
    parse_sam(merged_sam, out_sam.name, 'merged', 'tmp')
    # TODO: determine why joined reads have more soft-clips or single read matches
    parse_sam(joined_sam, out_sam.name, 'joined', 'tmp')
    # convert to sorted and indexed bam
    cmd = 'cat %s %s |samtools view -@ 4 -Shb |sambamba sort -m 4GB -t %s -o %s  /dev/stdin' % (args.header,
                                                                                                out_sam.name,
                                                                                                args.threads,
                                                                                                os.path.join(
                                                                                                    args.output_dir,
                                                                                                    'out.bam'))
    log = "make sorted bam file"
    run_subprocess([cmd], args, log)
    out_sam.close()
    return args


def clean(args):
    """delete non-used intermediate files"""
    log = 'removing tmp dir %s ' % (args.tmpdir)
    if args.tmpdir.endswith('STAR'):
        cmd = ['rm -rf %s' % (args.tmpdir)]
        run_subprocess(cmd, args, log)
    log = "remove tmp files from output dir"
    cmd = ['rm -rf %s/merged*' % args.output_dir]
    run_subprocess(cmd, args, log)
    cmd = ['rm -rf %s/joined*' % args.output_dir]
    run_subprocess(cmd, args, log)
    cmd = ['rm -rf %s/joined* header.sam' % args.output_dir]
    run_subprocess(cmd, args, log)


def main():
    """main function loop"""
    # 1 get command line arguments
    args = parse_args()
    # version = get_version()alignSoftClipAtReferenceEnds
    if __name__ == "__main__":
        log = open(args.log, 'w')
    else:
        log = open(args.log, 'a')
    log.write("started run\n")
    # 2 make reference genome fo STAR in appropriate directory
    args = index_STAR(args)
    # 4 map processed reads
    if not os.path.exists(os.path.join(args.output_dir,'header.sam')):
        args = process_reads_merged(args)
        args = process_reads_joined(args)
        args = map_STAR(args)
        args = make_header(args)
        args = bam_output(args)
    args = remove_PCR_duplicates(args)
    # clean(args)


if __name__ == '__main__':
    main()







