#!/usr/bin/env python
__author__ = 'thomasvangurp'
# Date created: 22/11/2014 (europe date)
# Function: Pipeline for mapping reads to reference
#Python version: 2.7.3
#External dependencies: samtools,pysam,methylation_calling.py
#Known bugs: None
#Modifications: None
import pysam
import argparse
import subprocess
import tempfile
import os
import shutil
import sys
from Bio import SeqIO
from Bio import Restriction



def getScriptPath():
    return os.path.dirname(__file__)

def parse_args():
    "Pass command line arguments"
    if not sys.argv[1:]:
        sys.argv.append('-h')
    parser = argparse.ArgumentParser(description='use bwameth for mapping reads')
    #input files
    parser.add_argument('-s','--sequences',
                        help='number of sequences to take for testing')
    parser.add_argument('--subsample_treshold',
                        help='Subsample treshold',default='100000')
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
    parser.add_argument('--refgenome',
                        help='reference genome instead of clusters')
    parser.add_argument('-b','--barcodes',
                    help='Barcodes used in output')
    parser.add_argument('--species',
                        help='Species: if selected only that species will be putin BAM RG header')
    parser.add_argument('--bamout',
                        help='output for bam file with RGs')
    parser.add_argument('--threads',
                        help='Number of threads to used where multithreading is possible')
    parser.add_argument('--log',
                        help='log of output operation')
    parser.add_argument('--output_dir',
                        help='optional: Choose output directory')
    parser.add_argument('--watson_vcf',
                        help='watson vcf output')
    parser.add_argument('--crick_vcf',
                        help='crick vcf output')
    parser.add_argument('--snp_vcf',
                        help='vcf output snp')
    parser.add_argument('--methylation_vcf',
                        help='Methylation vcf output')
    parser.add_argument('--heatmap',
                        help='heatmap output methylation')
    args = parser.parse_args()
    if args.input_dir:
        if not args.reads_R1:
            args.reads_R1 = os.path.join(args.input_dir,'Unassembled.R1.watson.fq.gz')
        if not args.reads_R2:
            args.reads_R2 = os.path.join(args.input_dir,'Unassembled.R2.crick.fq.gz')
        if not args.merged:
            args.merged = os.path.join(args.input_dir,'Assembled.fq.gz')
        if args.reference == None and args.refgenome == None:
            args.reference = os.path.join(args.input_dir,'consensus_cluster.renamed.fa')
        if args.barcodes == None:
            args.barcodes = os.path.join(args.input_dir,'barcodes.csv')
    if args.output_dir:
        if not os.path.exists(args.output_dir):
            os.mkdir(args.output_dir)
        if not args.log:
            args.log = os.path.join(args.output_dir,'mapping_variantcalling.log')
        args.watson_vcf = os.path.join(args.output_dir,'watson.vcf.gz')
        args.crick_vcf = os.path.join(args.output_dir,'crick.vcf.gz')
        args.snp_vcf = os.path.join(args.output_dir,'snp.vcf.gz')
        args.methylation_vcf = os.path.join(args.output_dir,'methylation.vcf.gz')
        args.heatmap = os.path.join(args.output_dir,'heatmap.igv')
        args.mastermeth = os.path.join(args.output_dir,'methylation.bed')
    return args

def run_subprocess(cmd,args,log_message):
    "Run subprocess under standardized settings"
    #force the cmds to be a string.
    if len(cmd) != 1:
        cmd = [" ".join(cmd)]
    with open(args.log,'a') as log:
        log.write("now starting:\t%s\n"%log_message)
        log.write('running:\t%s\n'%(' '.join(cmd)))
        log.flush()
        p = subprocess.Popen(cmd,stdout=subprocess.PIPE,stderr=subprocess.PIPE,shell=True,executable='/bin/bash')
        stdout, stderr = p.communicate()
        stdout = stdout.decode().replace('\r','\n')
        stderr = stderr.decode().replace('\r','\n')
        if stdout:
            log.write('stdout:\n%s\n'%stdout)
        if stderr:
            log.write('stderr:\n%s\n'%stderr)
        return_code = p.poll()
        if return_code:
            raise RuntimeError(stderr)
        log.write('finished:\t%s\n\n'%log_message)
    return 0

def make_header(in_files,args):
    """Make sam header given input file species and"""
    #parse input barcode file and make list of individuals

    in_files['header'] = 'location of header'
    return in_files

def run_bwameth(in_files,args):
    "run bwa_meth for mapping"

    in_files['bam_out'] = {}
    in_files['bam_out']['watson'] = os.path.join(args.output_dir,'watson.bam')
    in_files['bam_out']['crick'] = os.path.join(args.output_dir,'crick.bam')
    in_files['header'] = os.path.join(args.output_dir,'header.sam')

    log = "index renamed reference using bwameth"
    ref = args.reference
    if not os.path.exists('%s.bwameth.c2t'%ref):
        cmd = ['bwameth.py index %s'%ref]
        run_subprocess(cmd,args,log)

    log = "run bwameth for merged reads"
    if args.sequences:
        add = '|head -n %s'%(4*int(args.sequences))
    else:
        add = ''
    if args.merged:
        cmd = ['bwameth.py -t %s -p %s --reference %s <(pigz -cd %s %s) NA'%
               (args.threads,
                os.path.join(args.output_dir,'merged'),
                ref,
                args.merged,add
                )]
        run_subprocess(cmd,args,log)

    log = "run bwameth for non-merged reads"
    cmd = ['bwameth.py -t %s -p %s --reference %s <(pigz -cd %s %s) <(pigz -cd %s %s)'%
           (args.threads,
            os.path.join(args.output_dir,'pe'),
            ref,
            args.reads_R1,add,
            args.reads_R2,add
            )]
    run_subprocess(cmd,args,log)

    log = "get start of header from pe.bam, add RG information using addRG function"
    cmd = ["samtools view -H %s > %s"%
           ((os.path.join(args.output_dir,'pe.bam')),
            (os.path.join(args.output_dir,'header.sam')))]
    run_subprocess(cmd,args,log)


    log = "Append RG to header"
    in_files = addRG(in_files,args)

    log = "merge bam files"
    cmd = ["samtools merge -h %s -fcp -@ %s %s <(samtools reheader %s %s) <(samtools reheader %s %s)"%
           (in_files['header'], args.threads, os.path.join(args.output_dir,'combined.bam'), in_files['header'],
           os.path.join(args.output_dir, 'pe.bam'), in_files['header'],os.path.join(args.output_dir, 'merged.bam'))]

    run_subprocess(cmd,args,log)
    
    log = "index combined bam file"
    cmd = ["samtools index %s"%(os.path.join(args.output_dir,'combined.bam'))]
    run_subprocess(cmd, args, log)



    log = "split in watson and crick bam file"
    bam_input = pysam.AlignmentFile(os.path.join(args.output_dir,'combined.bam'),'rb')
    watson_output = pysam.AlignmentFile(os.path.join(args.output_dir,'watson.bam'),'wb', template=bam_input)
    crick_output = pysam.AlignmentFile(os.path.join(args.output_dir,'crick.bam'),'wb', template=bam_input)
    for record in bam_input:
        tag_dict = dict(record.tags)
        #remove reads with alternate mapping positions
        # if 'XA' in tag_dict:
        #     continue
        try:
            if (record.is_reverse and record.is_paired == False) or \
                (record.is_paired and record.is_read1 and record.is_reverse == True) or \
                (record.is_paired and record.is_read2 and record.is_reverse == False):
                    if tag_dict['ST'].lower() == 'crick':
                        watson_output.write(record)
                    elif tag_dict['ST'].lower() == 'watson':
                        crick_output.write(record)
            else:
                if tag_dict['ST'].lower() == 'watson':
                    watson_output.write(record)
                elif tag_dict['ST'].lower() == 'crick':
                    crick_output.write(record)
        except KeyError:
            continue
    watson_output.close()
    crick_output.close()

    in_files['bam_out'] = {}
    in_files['bam_out']['watson'] = os.path.join(args.output_dir,'watson.bam')
    in_files['bam_out']['crick'] = os.path.join(args.output_dir,'crick.bam')


    log = "index watson bam file"
    cmd = ["samtools index %s"%in_files['bam_out']['watson']]
    run_subprocess(cmd,args,log)
    log = "index crick bam file"
    cmd = ["samtools index %s"%in_files['bam_out']['crick']]
    run_subprocess(cmd,args,log)
    return in_files

def run_STAR(in_files, args):
    "run bwa_meth for mapping"

    in_files['bam_out'] = {}
    in_files['bam_out']['watson'] = os.path.join(args.output_dir, 'watson.dedup.bam')
    in_files['bam_out']['crick'] = os.path.join(args.output_dir, 'crick.dedup.bam')
    in_files['header'] = os.path.join(args.output_dir, 'header.sam')
    cmd = ["map_STAR.py",
           '--reads_R1 %s' % args.reads_R1,
           '--reads_R2 %s' % args.reads_R2,
           '--merged %s' % args.merged,
           "--barcodes %s" % args.barcodes,
           "--tmpdir %s" % args.tmpdir,
           "--threads %s" % args.threads,
           "--output_dir %s" % args.output_dir]
    if not args.reference:
        cmd += ['--refgenome %s' % args.refgenome]
    else:
        cmd += ['--reference %s' % args.reference]
    if args.sequences != None:
        cmd += ['--sequences %s' % args.sequences]
    log = "Map reads using STAR"
    run_subprocess(cmd, args, log)

    return in_files

def addRG(in_files,args):
    "make header for output bamfile and split in watson and crick"
    #define readgroup header lines by combining the following
    species_name = ''
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
            if args.species:
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

def get_enz(enz):
    """Get enzyme from biopython restriction library"""
    for enzyme in Restriction.AllEnzymes:
        if "%s"%(enzyme) == enz:
            return enzyme


def get_regions(contig,enzymes):
    """return loci with start and end locations"""
    out_sites = []
    enz_1 = get_enz(enzymes[0])
    enz_2 = get_enz(enzymes[1])
    enz_1_sites = enz_1.search(contig.seq)
    enz_2_sites = enz_2.search(contig.seq)
    combined_sites = sorted(enz_1_sites + enz_2_sites)
    for i in range(len(combined_sites)):
        site_A = combined_sites[i]
        try:
            site_B = combined_sites[i+1]
        except IndexError:
            break
        if site_B - site_A < 30:
            continue
        if site_A in enz_1_sites and site_B in enz_2_sites:
            out_sites.append((site_A + 1, site_B - len(enz_2.site)))
        elif site_A in enz_2_sites and site_B in enz_1_sites:
            out_sites.append((site_A + 1, site_B - len(enz_1.site)))
    return out_sites


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
                print('%s has %s reads and %s duplicates. Duplicate rate: %.2f%%'%(key,count,dup_count,100*dup_pct))
            else:
                print('%s has %s reads and 0 duplicates. Duplicate rate: 0%%' % (key, count))
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


def run_Freebayes(in_files,args):
    "run freebayes on watson and crick bam file with threadpool"
    in_files['variants'] = {}
    log = open(args.log,'a')
    for strand in ['watson','crick']:
        processes = set()
        max_processes = int(args.threads)
        outdir = tempfile.mkdtemp(prefix='vcf', dir=args.tmpdir)
        outlist = []
        in_files['variants'][strand] = outdir
        n = 0
        skipped = 0
        with open(in_files['header']) as header:
            for line in header:
                if line.startswith('@SQ'):
                    contig = line.split('\t')[1].split(':')[1]
                    length = line[:-1].split('\t')[2].split(':')[1]
                else:
                    continue
                #determine coverage on contig in bam file
                #set depth at 100.000.000
                n+=1
                if not n%1000:
                    print('Done processing %s contigs on %s,skipped %s'%(n,strand,skipped))
                bamfile = in_files['bam_out'][strand]
                cmd = ['samtools mpileup -d 10000000 %s -r %s:10-10'%(bamfile,contig)]
                # if int(contig) > 1000:
                #     break
                if int(length) <500:
                    p = subprocess.Popen(cmd,stdout=subprocess.PIPE,
                                               stderr=subprocess.PIPE,shell=True,executable='/bin/bash')
                    stdout, stderr = p.communicate()
                    return_code = p.poll()
                    stderr = stderr.replace('\r','\n')
                    if return_code:
                        raise RuntimeError(stderr)
                    try:
                        out = stdout.split('\t')
                        depth = int(out[3])
                        n_samples = int(stderr.split(' ')[1])
                    except IndexError:
                        #no reads for this contig, skip
                        continue
                else:
                    # it does not make sense to calculate depth here!
                    depth = 0
                if depth < n_samples * 10:
                    skipped += 1
                    continue
                outlist.append('%s.vcf'%contig)
                #freebayes --bam  <(samtools view -hs 0.89552238806 /tmp/watson.bam 1|samtools view -Shb -)
                # --fasta-reference /Volumes/data/epiGBS/Baseclear/unfiltered_sequences/seqNNAtlE/Carrot/consensus.clustered.renamed.fa
                # -F 0 -E 1 -C 0 -G 0 --haplotype-length 1 -k -K -X -u -i -q 21 -w -a
                # --report-all-haplotype-alleles --report-monomorphic --report-genotype-likelihood-max
                cmd = """freebayes -f %s -F 0 -E 1 \
                -C 0 -G 0 --haplotype-length 1 \
                --report-all-haplotype-alleles --report-monomorphic\
                 --report-genotype-likelihood-max \
                --haplotype-length 1 -KkXuiwaq 21 """%\
                      (args.reference)
                if depth > int(args.subsample_treshold):
                    factor = int(args.subsample_treshold) / float(depth)
                    log = "Subsample bam file with high coverage"
                    bamfile = "--bam <(samtools view -hs %s %s %s|samtools view -Shb -)"%\
                              (factor,bamfile,contig)
                    cmd += " %s > %s"%\
                    (bamfile, os.path.join(outdir,'%s.vcf'%contig))
                # elif depth == 0:
                #     cmd += " --bam %s > %s"%\
                #     (bamfile, os.path.join(outdir,'%s.vcf'%contig))
                else:
                    cmd += " -r %s:0-%s --bam %s > %s"%\
                    (contig, int(length)-1,bamfile, os.path.join(outdir,'%s.vcf'%contig))
                processes.add(subprocess.Popen(cmd,stdout=subprocess.PIPE,
                                               stderr=subprocess.PIPE,shell=True,executable='/bin/bash'))
                while len(processes) >= max_processes:
                    os.wait()
                    processes.difference_update([
                        p for p in processes if p.poll() is not None])
        #Make sure that all Freebayes processes are done before continuing to next st   ep.
        while len(processes):
            os.wait()
            processes.difference_update([
                        p for p in processes if p.poll() is not None])
        if strand == 'watson':
            target = args.watson_vcf
        else:
            target = args.crick_vcf
        print(outdir,outlist[0],target)
        shutil.move(os.path.join(outdir,outlist[0]),target)
        for vcf_file in outlist[1:]:
            file_in = os.path.join(outdir,vcf_file)
            cmd = ['grep -v "^#" %s >> %s'%(file_in,target)]
            p = subprocess.Popen(cmd,stdout=subprocess.PIPE,
                                       stderr=subprocess.PIPE,shell=True,executable='/bin/bash')
            stdout, stderr = p.communicate()
    return in_files

def variant_calling_samtools(in_files,args):
    """Do variant calling with freebayes"""
    #run mpileup on watson bam file
    in_files['vcf_out'] = {}
    in_files['vcf_out']['watson'] = os.path.join(args.output_dir,'watson.vcf.gz')
    in_files['vcf_out']['crick'] = os.path.join(args.output_dir,'crick.vcf.gz')

    cmd = ["samtools mpileup --reference %s -gt DP,AD,INFO/AD" % (args.reference) +
           " --max-depth  999999999 " +  # call at max-depth of 10.000.000
           "-q 0 " +  # Do not skip alignments with low mapQ
           "-Q 15 " +  # Skip bases with baseQ/BAQ smaller than 15
           "--skip-indels " +  # skip indels
           "-vu %s" % (
           in_files['bam_out']['watson']) +  # v = generate genotype likelihoods in VCF format u = uncompressed
           "|grep -v '^##contig='|bgzip -c > %s" % (in_files['vcf_out']['watson'])]

    log = "use samtools mpileup to get variant observation counts for watson"
    run_subprocess(cmd, args, log)


    cmd = ["samtools mpileup --reference %s -gt DP,AD,INFO/AD" % (args.reference) +
           " --max-depth  999999999 " + #call at max-depth of 10.000.000
           "-q 0 " + #Do not skip alignments with low mapQ #TODO: investigate option
           "-Q 15 " + #Skip bases with baseQ/BAQ smaller than 15
           "--skip-indels " + #skip indels
           "-vu %s" % (in_files['bam_out']['crick']) + #v = generate genotype likelihoods in VCF format u = uncompressed
           "|grep -v '^##contig='|bgzip -c > %s" % (in_files['vcf_out']['crick'])]

    log = "use samtools mpileup to get variant observation counts for crick"
    run_subprocess(cmd, args, log)
    return in_files

def merge_watson_crick(in_files, args):
    """create merged.tsv.gz with watson and crick calls merged"""
    if 'vcf_out' not in in_files:
        in_files['vcf_out'] = {}
        in_files['vcf_out']['watson'] = os.path.join(args.output_dir, 'watson.vcf.gz')
        in_files['vcf_out']['crick'] = os.path.join(args.output_dir, 'crick.vcf.gz')
    in_files['vcf_out']['merged'] = os.path.join(args.output_dir,'merged.tsv')
    cmd = ["merge_watson_crick.py",
           "-w %s" % in_files['vcf_out']['watson'],
           "-c %s" % in_files['vcf_out']['crick'],
           "-o %s" % in_files['vcf_out']['merged']]

    log = "Create custom tsv file for combining watson and crick observation counts per individual"
    run_subprocess(cmd, args, log)
    in_files['vcf_out']['merged'] = os.path.join(args.output_dir, 'merged.tsv.gz')
    return in_files

def SNP_calling(in_files, args):
    """run SNP calling"""
    if 'vcf_out' not in in_files:
        in_files['vcf_out'] = {}
    in_files['vcf_out']['SNP'] = os.path.join(args.output_dir, 'snp.vcf')
    in_files['vcf_out']['merged'] = os.path.join(args.output_dir, 'merged.tsv.gz')
    cmd = ["SNP_calling.py",
           "-m %s" % in_files['vcf_out']['merged'],
           "-s %s" % in_files['vcf_out']['SNP'],
           "-w %s" % os.path.join(args.output_dir, 'watson.vcf.gz')]
    log = "perform SNP calling"
    run_subprocess(cmd, args, log)

    return in_files


def methylation_calling(in_files,args):
    "run methylation calling script."
    log = ["Run methylation calling script"]
    in_files['vcf_out']['SNP'] = os.path.join(args.output_dir, 'snp.vcf.gz')
    in_files['vcf_out']['merged'] = os.path.join(args.output_dir, 'merged.tsv.gz')
    cmd = ["methylation_calling.py",
           " -r %s"%(args.reference),
           " -m %s"%(in_files['vcf_out']['merged']),
           " -s %s"%(in_files['vcf_out']['SNP']),
           " -o %s"%(os.path.join(args.output_dir,'methylation.bed')),
           " -heat %s"%(os.path.join(args.output_dir,'heatmap.igv'))
           ]



    run_subprocess(cmd,args,log)
    return in_files

def main():
    "Main function loop"
    args = parse_args()
    #Make sure log is empty at start
    if os.path.isfile(args.log):
        os.remove(args.log)
    #Step 1: discover files in input #todo
    files = {}
    #Step 2: map reads using STAR
    #TODO: replace for running map_STAR
    files = run_STAR(files,args)
    if args.refgenome:
        args.reference = args.refgenome
    files = variant_calling_samtools(files, args)
    files = merge_watson_crick(files,args)
    files = SNP_calling(files, args)
    files = methylation_calling(files,args)
    print('done')
if __name__ == '__main__':
    main()
