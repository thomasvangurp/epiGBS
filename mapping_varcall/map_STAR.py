#!/usr/bin/env pypy
import argparse
import subprocess
import os
import math
import gzip
import tempfile
from Bio import SeqIO

__author__ = 'thomasvangurp'
__description__ = "map reads orders of magnitudes faster using STAR"


def parse_args():
    """Pass command line arguments"""
    parser = argparse.ArgumentParser(description='use STAR for mapping reads')
    #input files
    parser.add_argument('-s', '--sequences',
                        help='number of sequences to take for testing')
    parser.add_argument('--tmpdir',
                        help='tmp directory',default="/tmp/")
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
    args = parser.parse_args()
    if args.input_dir:
        args.reads_R1 = os.path.join(args.input_dir,'Unassembled.R1.watson_trimmed.fq.gz')
        args.reads_R2 = os.path.join(args.input_dir,'Unassembled.R2.crick_trimmed.fq.gz')
        args.merged = os.path.join(args.input_dir,'Assembled.trimmed.fq.gz')
        args.reference = os.path.join(args.input_dir,'consensus_cluster.renamed.fa')
    if args.output_dir:
        if not os.path.exists(args.output_dir):
            os.mkdir(args.output_dir)
        if 'log' not in args:
            args.log = os.path.join(args.output_dir,'mapping_variantcalling.log')
        args.watson_vcf = os.path.join(args.output_dir,'watson.vcf')
        args.crick_vcf = os.path.join(args.output_dir,'crick.vcf')
        args.snp_vcf = os.path.join(args.output_dir,'snp.vcf')
        args.methylation_vcf = os.path.join(args.output_dir,'methylation.vcf')
        args.heatmap = os.path.join(args.output_dir,'heatmap.igv')
        #2 bed files should be made for subsequent analysis using Rnbeads or other software
        args.mastermeth = os.path.join(args.output_dir,'methylation.bed')
    args.tmpdir = tempfile.mkdtemp(suffix='STAR', prefix='tmp', dir=args.tmpdir)
    return args

def get_version():
    """get version of current script"""
    parent_dir = os.path.dirname(os.path.realpath(__file__))
    while True:
        if '.git' in os.listdir(parent_dir):
            break
        parent_dir = os.path.dirname(parent_dir)
    git_log = os.path.join(parent_dir,'.git','logs','HEAD')



def remove_PCR_duplicates(args):
    """us subprocess to run external remove_PCR_duplicates routine"""
    cmd =  "mark_PCR_duplicates.py --input_dir %s" % args.output_dir
    cmd += " -b %s" % args.barcodes
    cmd += " -r %s" % args.reference
    log = "Removal of PCR duplicates"
    run_subprocess([cmd], args, log)

def run_subprocess(cmd,args,log_message):
    "Run subprocess under standardized settings"
    #force the cmds to be a string.
    if len(cmd) != 1:
        cmd = [" ".join(cmd)]
    with open(args.log, 'a') as log:
        log.write("now starting:\t%s\n" % log_message)
        log.write('running:\t%s\n' % (' '.join(cmd)))
        log.flush()
        p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True, executable='/bin/bash')
        stdout, stderr = p.communicate()
        stdout = stdout.replace('\r', '\n')
        stderr = stderr.replace('\r', '\n')
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
    watson_merged = tempfile.NamedTemporaryFile(suffix=".fastq", prefix='watson_merged', dir=args.tmpdir, delete=False)
    crick_merged = tempfile.NamedTemporaryFile(suffix=".fastq", prefix='crick_merged', dir=args.tmpdir, delete=False)
    args.watson_merged = watson_merged.name
    args.crick_merged = crick_merged.name
    watson_convert = ['C', 'T']
    crick_convert = ['G', 'A']
    print 'Started processing merged reads'
    if args.merged.endswith('.gz'):
        file_in_handle = gzip.open(args.merged, 'rb')
    else:
        file_in_handle = open(args.merged, 'r')
    watson_out_handle = open(watson_merged.name, 'w')
    crick_out_handle = open(crick_merged.name, 'w')
    j = 0
    while True:
        read = []
        for i in range(4):
            try:
                read.append(file_in_handle.next())
            except StopIteration:
                break
        j += 1
        if not j % 1000000:
            print 'Processed %s reads' % j
        if not read:
            break
        if 'watson' in read[0].lower():
            convert_seq = read[1].upper().replace(watson_convert[0], watson_convert[1])
            c_pos = [str(n) for n,i in enumerate(read[1]) if i.upper() == 'C']
            header = '@%s' % (read[0][1:-1].replace(' ', '|').replace('\t', '|'))
            header += '|%s\n' % (','.join(c_pos))
            watson_out_handle.write(header + convert_seq + '+\n' + read[3])
        else:
            convert_seq = read[1].upper().replace(crick_convert[0], crick_convert[1])
            g_pos = [str(n) for n, i in enumerate(read[1]) if i.upper() == 'G']
            header = '@%s' % (read[0][1:-1].replace(' ', '|').replace('\t', '|'))
            header += '|%s\n' % (','.join(g_pos))
            crick_out_handle.write(header + convert_seq + '+\n' + read[3])
    watson_out_handle.close()
    crick_out_handle.close()
    return args

def process_reads_joined(args):
    """process reads and make them ready for mapping with STAR"""

    watson_joined_r1 = tempfile.NamedTemporaryFile(suffix=".fastq", prefix='watson_joined', dir=args.tmpdir,
                                                   delete=False)
    watson_joined_r2 = tempfile.NamedTemporaryFile(suffix=".fastq", prefix='watson_joined', dir=args.tmpdir,
                                                   delete=False)
    crick_joined_r1 = tempfile.NamedTemporaryFile(suffix=".fastq", prefix='crick_joined', dir=args.tmpdir,
                                                   delete=False)
    crick_joined_r2 = tempfile.NamedTemporaryFile(suffix=".fastq", prefix='crick_joined', dir=args.tmpdir,
                                                   delete=False)
    args.watson_joined_r1 = watson_joined_r1.name
    args.watson_joined_r2 = watson_joined_r2.name
    args.crick_joined_r1 = crick_joined_r1.name
    args.crick_joined_r2 = crick_joined_r2.name
    watson_convert = ['C', 'T']
    crick_convert = ['G', 'A']
    print 'Started processing joined reads'
    if args.reads_R1.endswith('.gz'):
        r1_handle = gzip.open(args.reads_R1, 'rb')
        r2_handle = gzip.open(args.reads_R2, 'rb')
    else:
        r1_handle = open(args.reads_R1, 'rb')
        r2_handle = open(args.reads_R2, 'rb')
    #make 4 file handles for forward and reverse watson and crick
    watson_r1_handle = open(args.watson_joined_r1, 'w')
    watson_r2_handle = open(args.watson_joined_r2, 'w')
    crick_r1_handle = open(args.crick_joined_r1, 'w')
    crick_r2_handle = open(args.crick_joined_r2, 'w')
    j = 0
    while True:
        read_r1 = []
        read_r2 = []
        for i in range(4):
            try:
                read_r1.append(r1_handle.next())
                read_r2.append(r2_handle.next())
            except StopIteration:
                break
        j += 1
        if not j % 1000000:
            print 'Processed %s reads' % (j)
        if not read_r1:
            break
        if 'watson' in read_r1[0].lower():
            convert_r1 = read_r1[1].upper().replace('C', 'T')
            convert_r2 = read_r2[1].upper().replace('G', 'A')
            c_pos = [str(n) for n, i in enumerate(read_r1[1]) if i.upper() == 'C']
            g_pos = [str(n) for n, i in enumerate(read_r2[1].rstrip('\n')[::-1]) if i.upper() == 'G']
            header = '@%s' % (read_r1[0][1:-1].replace(' ', '|').replace('\t', '|'))
            header += '|%s\n' % (','.join(c_pos) + '|' + ','.join(g_pos))
            watson_r1_handle.write(header + convert_r1 + '+\n' + read_r1[3])
            watson_r2_handle.write(header + convert_r2 + '+\n' + read_r2[3])
        else:
            convert_r1 = read_r1[1].upper().replace('G', 'A')
            convert_r2 = read_r2[1].upper().replace('C', 'T')
            g_pos = [str(n) for n, i in enumerate(read_r1[1]) if i.upper() == 'G']
            c_pos = [str(n) for n, i in enumerate(read_r2[1].rstrip('\n')[::-1]) if i.upper() == 'C']
            header = '@%s' % (read_r1[0][1:-1].replace(' ', '|').replace('\t', '|'))
            header += '|%s\n' % (','.join(g_pos) + '|' + ','.join(c_pos))
            crick_r1_handle.write(header + convert_r1 + '+\n' + read_r1[3])
            crick_r2_handle.write(header + convert_r2 + '+\n' + read_r2[3])
    crick_r1_handle.close()
    crick_r2_handle.close()
    watson_r1_handle.close()
    watson_r2_handle.close()
    return args


def index_STAR(args):
    """make STAR index for merged and joined reads"""

    # make STAR index folder for merged path
    merged_STAR_watson_index = os.path.join(args.output_dir,'STAR_merged_watson')
    merged_STAR_crick_index = os.path.join(args.output_dir,'STAR_merged_crick')
    if not os.path.exists(merged_STAR_watson_index):
        os.mkdir(merged_STAR_watson_index)
        os.mkdir(merged_STAR_crick_index)
    ref_merged_watson = os.path.join(merged_STAR_watson_index, '%s.merged.watson.fa' % args.species)
    ref_merged_crick = os.path.join(merged_STAR_crick_index, '%s.merged.crick.fa' % args.species)
        
    #make STAR index folder for joined path
    joined_STAR_watson_index = os.path.join(args.output_dir,'STAR_joined_watson')
    joined_STAR_crick_index = os.path.join(args.output_dir,'STAR_joined_crick')
    if not os.path.exists(joined_STAR_watson_index):
        os.mkdir(joined_STAR_watson_index)
        os.mkdir(joined_STAR_crick_index)
    ref_joined_watson = os.path.join(joined_STAR_watson_index, '%s.joined.watson.fa' % args.species)
    ref_joined_crick = os.path.join(joined_STAR_crick_index, '%s.joined.crick.fa' % args.species)

    #get file handle for input reference file
    try:
        file_handle = open(args.reference, 'r')
    except IOError:
        raise IOError('file %s does not exist' % args.reference)

    #iterate over input lines and write to references
    joined_len = 0
    merged_len = 0
    joined_count = 0
    merged_count = 0
    ref_merged_watson_handle = open(ref_merged_watson,'w')
    ref_merged_crick_handle = open(ref_merged_crick,'w')
    ref_joined_watson_handle = open(ref_joined_watson,'w')
    ref_joined_crick_handle = open(ref_joined_crick,'w')
    seq = ''
    for line in file_handle:
        if line.startswith('>'):
            if seq != '':
                if 'NNNNNNNN' in seq.upper():
                    joined_len += len(seq)
                    joined_count += 1
                    ref_joined_watson_handle.write(header + seq.upper().replace('C','T')+'\n')
                    ref_joined_crick_handle.write(header + seq.upper().replace('G','A')+'\n')
                else:
                    merged_len += len(seq)
                    merged_count += 1
                    ref_merged_watson_handle.write(header + seq.upper().replace('C', 'T') + '\n')
                    ref_merged_crick_handle.write(header + seq.upper().replace('G', 'A') + '\n')
            seq = ''
            header = line
        else:
            seq += line.rstrip('\n')
    #write final sequence, this is always merged
    merged_len += len(seq)
    merged_count += 1
    ref_merged_watson_handle.write(header + seq.upper().replace('C', 'T') + '\n')
    ref_merged_crick_handle.write(header + seq.upper().replace('G', 'A') + '\n')
    #close file handles
    ref_joined_watson_handle.close()
    ref_joined_crick_handle.close()
    ref_merged_watson_handle.close()
    ref_merged_crick_handle.close()
    #MAKE LIST for indexes to be made 
    index_list = [(joined_len, joined_count, joined_STAR_watson_index, ref_joined_watson),
                  (joined_len, joined_count, joined_STAR_crick_index, ref_joined_crick),
                  (merged_len, merged_count, merged_STAR_watson_index, ref_merged_watson),
                  (merged_len, merged_count, merged_STAR_crick_index, ref_merged_crick)]
    #calculate parameters for indexing reference for merged and joined reads.
    for (genome_len, no_clusters, genome_dir, ref) in index_list:
        index_cmd = 'STAR --runThreadN %s --runMode genomeGenerate --genomeDir %s'%(args.threads,genome_dir)
        fasta_file = [file for file in os.listdir(genome_dir) if file.endswith('.fa')][0]
        index_cmd += ' --genomeFastaFiles %s'%os.path.join(genome_dir,fasta_file)
        genomeSAindexNbases = min(14, math.log(genome_len,2)/2 - 1)
        index_cmd += ' --genomeSAindexNbases %i'%genomeSAindexNbases
        genomeChrBinNbits = min(18, math.log(genome_len/no_clusters, 2))
        index_cmd += ' --genomeChrBinNbits %i' % genomeChrBinNbits
        log = 'making STAR index of %s'%(ref)
        if 'Genome' not in os.listdir(genome_dir):
            run_subprocess([index_cmd], args, log)
    return args

def map_STAR(args):
    """map reads with STAR"""
    for type in ['joined', 'merged']:
        for strand in ['watson', 'crick']:
            if strand == 'watson':
                n = 1
            else:
                n = 3
            STAR_index_dir = os.path.join(args.output_dir,'STAR_%s_%s'%(type, strand))
            cmd = "STAR --runThreadN %s --genomeDir %s"%(args.threads, STAR_index_dir)

            if type == 'merged':
                cmd += " --readFilesIn %s" % vars(args)['%s_%s' % (strand, type)]
            else:
                cmd += " --readFilesIn %s " %   vars(args)['%s_%s_r1' % (strand, type)]
                cmd += " %s" %                  vars(args)['%s_%s_r2' % (strand, type)]

            cmd += " --outSAMattributes NM MD AS --outSAMtype SAM"
            cmd += " --outFileNamePrefix %s" % (os.path.join(args.output_dir,'%s_%s'%(type,strand)))
            #outFilterScoreMinOverLread : float: sam as outFilterMatchNmin, but normalized to the read length (sum of mates’ lengths for paired-end reads)
            #outFilterMatchNminOverLread: float: same as outFilterScoreMin, but normalized to read length (sum of mates’ lengths for paired-end reads)

            cmd += " --outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 "

            # –outFilterMultimapNmax 1 int: maximum number of loci the read is allowed to map to. Alignments (all of
            # them) will be output only if the read maps to no more loci than this value.
            cmd += " -–outFilterMultimapNmax 1 --outFilterMismatchNoverLmax 0.95"
            # TODO: implement --alignEndsType endtoend mapping after joined reads are merged
            cmd += "--outFilterMatchNminOverLread 0.9 --scoreGap -4 --seedPerReadNmax 5000" \
                   " --alignEndsType Extend5pOfRead1" \
                   " --alignSoftClipAtReferenceEnds No" \
                   " --outSAMorder PairedKeepInputOrder" \
                   " --outFilterMultimapNmax 1" \
                   " --scoreInsOpen -1" \
            #make sure we have a bam file sorted by name
            log = "run STAR for % strand on %s reads"%(strand, type)
            run_subprocess([cmd],args, log)
            log = "write final log of STAR to normal log"
            cmd = "cat %s " % os.path.join(args.output_dir, '%s_%s' % (type, strand) + 'Log.final.out')
            run_subprocess([cmd], args, log)
    return args


def parse_sam(in_file, out_file, read_type , strand):
    """parse sam file and write correct output"""
    out_handle = open(out_file , 'a')
    if strand == 'watson':
        nt = ['C']
    else:
        nt = ['G']
    count = 0
    print 'Warning, only works for forward mapped reads'
    mismatch = 0
    clip_count_total = 0
    for line in open(in_file, 'r'):
        modulo_line_no = count % 2
        #alternates between 0 and 1
        if line.startswith('@'):
            continue
        split_line = line.rstrip('\n').split('\t')
        #skip read pairs with improper flags.
        #TODO: do this filtering in mark_PCR_duplicates or elsewhere with access to pysam.
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
        meth_pos_list = header[6:]
        out_line = [header[0]]
        out_line += split_line[1:9]
        seq = list(split_line[9])
        try:
            meth_pos = [int(n) for n in meth_pos_list[-modulo_line_no].split(',')]
            for n in meth_pos:
                assert seq[n] in ['T','A']
                seq[n] = nt[-modulo_line_no]
        except ValueError:
            pass
        out_line += [''.join(seq)]
        out_line += split_line[10:]
        out_line += header[3:6]
        out_handle.write('\t'.join(out_line) + '\n')
        count += 1
    print '%s mismatches out of %s' % (mismatch, count)
    print '%s reads out of  %s soft clipped more than 5' % (clip_count_total, count)

def addRG(in_files,args):
    """make header for output bamfile and split in watson and crick"""
    #define readgroup header lines by combining the following

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

    with open(args.barcodes,'r') as barcodes:
        sam_out= open(in_files['header'],'a')
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
            RG.append('ID:%s_%s_%s'%(fc,lane,sample))
            RG.append('SM:%s'%(sample))
            RG.append('LB:%s_%s'%(fc,sample))
            RG.append('PL:ILLUMINA\n')
            sam_out.write('\t'.join(RG))
    sam_out.close()
    return in_files


def make_header(args):
    """Make header for watson and crick bam file"""
    header = os.path.join(args.output_dir,'header.sam')
    args.header = header
    header_handle = open(header,'w')
    header_handle.write('@HD\tVN:1.4\n')
    joined_sam = open(os.path.join(args.output_dir, 'joined_watsonAligned.out.sam'))
    merged_sam = open(os.path.join(args.output_dir, 'merged_watsonAligned.out.sam'))
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
    in_files = {'header':os.path.join(args.output_dir,'header.sam')}
    addRG(in_files, args)
    return args


def bam_output(args):
    """Generate watson and crick output bam file"""

    for strand in ['watson', 'crick']:
        merged_sam = os.path.join(args.output_dir, 'merged_%sAligned.out.sam' % strand)
        joined_sam = os.path.join(args.output_dir, 'joined_%sAligned.out.sam' % strand)
        out_sam = tempfile.NamedTemporaryFile(prefix=strand, suffix='.sam', dir=args.output_dir)
        #rewrite sam file merged and joined for watson and crick
        parse_sam(merged_sam, out_sam.name, 'merged', strand)
        #TODO: determine why joined reads have more soft-clips or single read matches
        parse_sam(joined_sam, out_sam.name, 'joined', strand)
        #convert to sorted and indexed bam
        cmd = 'cat %s %s |samtools view -@ 4 -Shb |sambamba sort -m 4GB -t %s -o %s  /dev/stdin'%(args.header,
                                                                            out_sam.name,args.threads,
                                                            os.path.join(args.output_dir,'%s.bam' % strand) )
        log = "make sorted bam file"
        run_subprocess([cmd], args, log)
        out_sam.close()
    return args


def clean(args):
    """delete non-used intermediate files"""
    log =  'removing tmp dir %s ' % (args.tmpdir)
    if args.tmpdir.endswith('STAR'):
        cmd = ['rm -rf %s' % (args.tmpdir)]
        run_subprocess(cmd,args,log)
    log = "remove sam file outputs from output dir"
    cmd = ['rm %s/*Aligned.out.sam' % args.output_dir]
    # run_subprocess(cmd, args, log)


def main():
    """main function loop"""
    #1 get command line arguments
    args = parse_args()
    version = get_version()
    log = open(args.log,'w')
    log.write("started run\n")
    #2 make reference genome fo STAR in appropriate directory
    args = index_STAR(args)
    #3 rewrite fastq files to contain
    args = process_reads_joined(args)
    args = process_reads_merged(args)
    #4 map processed reads
    args = map_STAR(args)
    args = make_header(args)
    args = bam_output(args)
    args = remove_PCR_duplicates(args)
    clean(args)

if __name__ == '__main__':
    main()







