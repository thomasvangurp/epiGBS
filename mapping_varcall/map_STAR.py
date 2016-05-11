#!/usr/bin/env pypy
import argparse
import subprocess
import os
import math
import gzip
import tempfile

__author__ = 'thomasvangurp'
__description__ = "map reads orders of magnitudes faster using STAR"


def parse_args():
    """Pass command line arguments"""
    parser = argparse.ArgumentParser(description='use STAR for mapping reads')
    #input files
    parser.add_argument('-s','--sequences',
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
    return args

def remove_PCR_duplicates(in_files,args):
    """Remove PCR duplicates and non-paired PE-reads per cluster"""
    #check if random tag is present in fastq file, otherwise do not perform function
    # fastq_tags = open(in_files[''])
    #TODO: implement sample specific PCR duplicate detection
    for strand,bamfile in in_files['bam_out'].items():
        clusters = SeqIO.parse(open(args.reference),'fasta')
        handle = pysam.AlignmentFile(bamfile,'rb')
        out_bam = tempfile.NamedTemporaryFile(suffix='uniq.bam',dir=args.output_dir,delete=False)
        out_handle = pysam.AlignmentFile(out_bam.name,'wb', template=handle)
        read_count = {}
        for cluster in clusters:
            enzymes = ["Csp6I","NsiI"]
            if len(cluster.seq) > 350:
                #this must be a reference genome / chromosome: look for regions with mapping reads
                regions = get_regions(cluster,enzymes)
            else:
                regions = [None]
            for region in regions:
                if region:
                    reads = handle.fetch(cluster.id,region[0],region[1])
                else:
                    reads = handle.fetch(cluster.id)
                if 'NNNNNNNN' in cluster._seq.upper() and not region:
                    cluster_is_paired = True
                elif region:
                    if region[1] - region[0] > 240:
                        cluster_is_paired = True
                    else:
                        cluster_is_paired = False
                else:
                    cluster_is_paired = False
                read_out = {}
                for read in reads:
                    tag_dict = dict(read.tags)
                    try:
                        tag = tag_dict['RN']
                        sample = tag_dict['RG']
                        AS = tag_dict['AS']
                    except KeyError:
                        break
                    if not read.is_proper_pair and cluster_is_paired:
                        continue
                    if sample not in read_out:
                        read_out[sample] = {}
                    if tag not in read_out[sample]:
                        read_out[sample][tag] = {read.qname:AS}
                    else:
                        try:
                            read_out[sample][tag][read.qname]+= AS
                        except KeyError:
                            read_out[sample][tag][read.qname] = AS
                #process read_out
                if read_out == {} and 'RN' not in tag_dict:
                    #random tag not yet implemented. return in_files and do not process further
                    return in_files
                if region:
                    reads = handle.fetch(cluster.id, region[0], region[1])
                else:
                    reads = handle.fetch(cluster.id)
                for read in reads:
                    if not read.is_proper_pair and cluster_is_paired:
                        continue
                    # if not read_count%100000:
                    #     print '%s reads processed for %s strand'%(read_count,strand)
                    tag_dict = dict(read.tags)
                    tag = tag_dict['RN']
                    sample = tag_dict['RG']
                    try:
                        read_count[sample]['count'] += 1
                    except KeyError:
                        if sample not in read_count:
                            read_count[sample] = {'count':1}
                        else:
                            read_count[sample]['count'] =  1
                    max_AS = max(read_out[sample][tag].values())
                    qname = [name for name,AS in read_out[sample][tag].items() if AS == max_AS][0]
                    if read.qname == qname:
                        out_handle.write(read)
                    else:
                        try:
                            read_count[sample]['dup_count'] += 1
                        except KeyError:
                            read_count[sample]['dup_count'] = 1
        for key , subdict in sorted(read_count.items()):
            count = subdict['count']
            if 'dup_count' in subdict:
                dup_count = subdict['dup_count']
                dup_pct = dup_count / float(count)
                print '%s has %s reads and %s duplicates. Duplicate rate: %.2f%%'%(key,count,dup_count,100*dup_pct)
            else:
                print '%s has %s reads and 0 duplicates. Duplicate rate: 0%%' % (key, count)
        # out_handle.flush()
        out_handle.close()
        old_bam = in_files['bam_out'][strand]
        log = "move old bam file %s to %s"%(old_bam,old_bam.replace('.bam','.old.bam'))
        cmd = ["mv %s %s"%(old_bam,old_bam.replace('.bam','.old.bam'))]
        run_subprocess(cmd,args,log)
        log = "move uniq bam file %s to %s"%(out_bam.name,old_bam)
        cmd = ["mv %s %s"%(out_bam.name,in_files['bam_out'][strand])]
        run_subprocess(cmd,args,log)
        log = "index bam file %s"%(old_bam)
        cmd = ["samtools index %s"%(in_files['bam_out'][strand])]
        run_subprocess(cmd,args,log)
    return in_files

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


def process_reads(args):
    """process reads and make them ready for mapping with STAR"""
    watson_merged = tempfile.NamedTemporaryFile(suffix=".fastq", prefix='watson_merged', dir=args.tmpdir, delete=False)
    crick_merged = tempfile.NamedTemporaryFile(suffix=".fastq", prefix='crick_merged', dir=args.tmpdir, delete=False)
    watson_r1 = tempfile.NamedTemporaryFile(suffix=".fastq", prefix='watson_r1', dir=args.tmpdir, delete=False)
    crick_r1 = tempfile.NamedTemporaryFile(suffix=".fastq", prefix='crick_r1', dir=args.tmpdir, delete=False)
    watson_r2 = tempfile.NamedTemporaryFile(suffix=".fastq", prefix='watson_r2', dir=args.tmpdir, delete=False)
    crick_r2 = tempfile.NamedTemporaryFile(suffix=".fastq", prefix='crick_r2', dir=args.tmpdir, delete=False)
    files = {'merged': (args.merged, watson_merged.name, ['C','T'], crick_merged.name, ['G','A']),
             'reads_forward': (args.reads_R1, watson_r1.name, ['C', 'T'], crick_r1.name, ['G', 'A']),
             'reads_reverse': (args.reads_R2, watson_r2.name, ['G', 'A'], crick_r2.name, ['C', 'T'])}
    for name, (file_in, watson_out, watson_convert, crick_out, crick_convert) in files.items():
        print 'Started processing %s' % name
        if file_in.endswith('.gz'):
            file_in_handle = gzip.open(file_in, 'rb')
        else:
            file_in_handle = open(file_in, 'r')
        watson_out_handle = open(watson_out, 'w')
        crick_out_handle = open(crick_out, 'w')
        j = 0
        while True:
            read = []
            for i in range(4):
                try:
                    read.append(file_in_handle.next())
                except StopIteration:
                    break
            j+=1
            if not j % 10000:
                print 'Processed %s reads in %s' % (j, name)
            if not read:
                break
            if 'watson' in read[0].lower():
                convert_seq = read[1].replace(watson_convert[0], watson_convert[1])
                # header = '@%s' % (read[0][1:-1].replace(' ', '|').replace('\t', '|'))
                # header += '|' + read[1]
                watson_out_handle.write(read[0] + convert_seq + '+\n' + read[3])
            else:
                convert_seq = read[1].replace(crick_convert[0], crick_convert[1])
                # header = '@%s' % (read[0][1:-1].replace(' ', '|').replace('\t', '|'))
                # header += '|' + read[1]
                crick_out_handle.write(read[0] + convert_seq + '+\n' + read[3])
        watson_out_handle.close()
        crick_out_handle.close()
    return files

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
    reads = process_reads(args)
    for type in ['joined', 'merged']:
        for strand in ['watson', 'crick']:
            if strand == 'watson':
                n = 1
            else:
                n = 3
            STAR_index_dir = os.path.join(args.output_dir,'STAR_%s_%s'%(type, strand))
            cmd = "STAR --runThreadN %s --genomeDir %s"%(args.threads, STAR_index_dir)
            if type == 'merged':
                if args.merged.endswith('.gz'):
                    cmd += " --readFilesIn %s" % reads['merged'][n]
            else:
                if args.reads_R1.endswith('.gz'):
                    cmd += " --readFilesIn %s %s "% \
                           (reads['reads_forward'][n], reads['reads_reverse'][n])
            cmd += " --outSAMattributes NH HI NM MD AS --outSAMtype BAM Unsorted"
            cmd += " --outFileNamePrefix %s" % (os.path.join(args.output_dir,'%s_%s'%(type,strand)))
            cmd += " --outStd BAM_Unsorted"
            cmd += " --outSAMreadID Number"
            #outFilterScoreMinOverLread : float: sam as outFilterMatchNmin, but normalized to the read length (sum of mates’ lengths for paired-end reads)
            #outFilterMatchNminOverLread: float: same as outFilterScoreMin, but normalized to read length (sum of mates’ lengths for paired-end reads)

            cmd += " --outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 "

            # –outFilterMultimapNmax 1 int: maximum number of loci the read is allowed to map to. Alignments (all of
            # them) will be output only if the read maps to no more loci than this value.
            cmd += " -–outFilterMultimapNmax 1 --outFilterMismatchNoverLmax 0.95"
            #make sure we have a bam file sorted by name
            cmd += "|sambamba sort -N -m 4GB -t 6 -o %s.namesorted.bam /dev/stdin" % \
                   (os.path.join(args.output_dir,'%s_%s'%(type,strand)))
            log = "run STAR for % strand on %s reads"%(strand, type)
            run_subprocess([cmd],args, log)
            log = "write final log of STAR to normal log"
            cmd = "cat %s " % os.path.join(args.output_dir, '%s_%s' % (type, strand) + 'Log.final.out')
            run_subprocess([cmd], args, log)
    return args

def make_bam(args):
    for strand in ['watson', 'crick']:
        for type in ['merged', 'joined']:
            bam_file = os.path.join(args.output_dir,'%s_%s.namesorted.bam' % (type, strand))
            if type == 'merged':
                fastq_file  = open('/tmp/watson_mergedfNqeem.fastq','r')
            cmd = "samtools view -@ 4 %s" % bam_file

            p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True,
                                 executable='/bin/bash')
            out_file = open('/tmp/test.sam','w')
            record, previous_record = None, None
            records_processed = 0
            for line in p.stdout:
                split_line = line.rstrip('\n').split('\t')
                if split_line[0] == previous_record:
                    #multimapping read is present, ignore
                    record = None
                    continue
                else:
                    if record:
                        out_file.write(record)
                while True:
                    fastq_entry = [fastq_file.readline() for i in range(0, 4)]
                    records_processed += 1
                    try:
                        assert records_processed == int(split_line[0])
                        break
                    except AssertionError:
                        pass
                read_name = fastq_entry[0].replace(' ', '\t').split('\t')[0]
                read_tags = fastq_entry[0].replace(' ', '\t').split('\t')[1:]

                line_out = [read_name] + split_line[1:9] + [fastq_entry[1].rstrip('\n')] + [split_line[10]]
                line_out += split_line[11:]
                line_out += read_tags
                previous_record = split_line[0]
                record = '\t'.join(line_out)
            p.communicate()
            out_file.close()
            break


def parse_sam(in_file, out_file):
    """parse sam file and write correct output"""
    out_handle = open(out_file , 'a')
    for line in open(in_file, 'r'):
        if line.startswith('@'):
            continue
        split_line = line.rstrip('\n').split('\t')
        header = split_line[0].split('|')
        out_line = [header[0]]
        out_line += split_line[1:9]
        out_line += [header[-1]]
        out_line += split_line[10:]
        out_line += header[1:-1]
        out_handle.write('\t'.join(out_line) + '\n')


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
        parse_sam(merged_sam, out_sam.name)
        #TODO: fix joined output
        # parse_sam(joined_sam, out_sam.name)
        #convert to sorted and indexed bam
        out_bam = tempfile.NamedTemporaryFile(prefix=strand, suffix='.sorted.bam', dir=args.output_dir,delete=False)
        cmd = 'cat %s %s |samtools view -F 256 -@ 4 -Shb |sambamba sort -m 4GB -o %s -t %s /dev/stdin'%(args.header,
                                                                            out_sam.name, out_bam.name, args.threads)
        log = "make sorted bam file"
        run_subprocess([cmd], args, log)




def main():
    """main function loop"""
    #1 get command line arguments
    args = parse_args()
    log = open(args.log,'w')
    log.write("started run\n")
    make_bam(args)
    #2 make reference genome fo STAR in appropriate directory
    # args = index_STAR(args)
    # args = map_STAR(args)
    # args = make_header(args)
    # bam_output(args)
if __name__ == '__main__':
    main()







