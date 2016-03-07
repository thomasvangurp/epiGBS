#!/usr/bin/env python
__author__ = 'thomasvangurp'
# Date created: 20/11/2014 (europe date)
# Function: Pipeline for creation of reference
#Python version: 2.7.3
#External dependencies: vsearch,samtools,vcfutils.pl,rename_fast.py
#Known bugs: None
#Modifications: None
import argparse
import subprocess
import tempfile
import os
import shutil
from distutils.spawn import find_executable
from Bio import SeqIO
import gzip
# usearch = "usearch"
# seqtk = "seqtk"
# pear = "/mnt/data/tools/pear_merge/default/pear"
# mergeBSv3 = "~/epiGBS/mergeBSv3.py"
# create_consensus = "~/epiGBS/create_consensus.py"
# vcfutils = "~/epiGBS/vcfutils.pl"
origWD = os.getcwd()
os.chdir(origWD)


dependencies = {}
dependencies['samtools'] = '<0.1.18'
dependencies['vcfutils.pl'] = ''
dependencies['usearch'] = 'usearch_8.0.1409'
dependencies['seqtk'] = '1.0-r31'
dependencies['pear'] = 'v0.9.7'
dependencies['pigz'] = ''

usearch = "usearch_8.0.1409_i86osx32"
vsearch = "vsearch"
seqtk = "/opt/bin/seqtk"
pear = "pear"
mergeBSv3 = "mergeBSv3.py"
create_consensus = "create_consensus.py"
vcfutils = "vcfutils.pl"

def parse_args():
    "Pass command line arguments"
    parser = argparse.ArgumentParser(description='Process input files')
    #input files
    parser.add_argument('-s','--sequences',
                        help='number of sequences to take for testing, useful for debugging')
    parser.add_argument('--forward',
                        help='forward reads fastq')
    parser.add_argument('--reverse',
                    help='reverse reads fastq')
    parser.add_argument('--min_unique_size',default="2",
                    help='Minimum unique cluster size')
    parser.add_argument('--clustering_treshold',default="0.95",
                    help='Clustering treshold for final clustering step')
    parser.add_argument('-t','--tmpdir',
                        help='tmp directory',default='/tmp')
    parser.add_argument('--threads',
                        help='Number of threads to used where multithreading is possible')
    parser.add_argument('--outputdir',
                        help='Optional: output directory')
    parser.add_argument('--log',
                        help='log of output operation')
    parser.add_argument('--samout',#TODO: describe what is the purpose of this file
                        help='sam output of clustering process')
    parser.add_argument('--consensus',#TODO: describe what is the purpose of this file
                        help='consensus output')
    parser.add_argument('--consensus_cluster',#TODO: describe what is the purpose of this file
                    help='consensus clustering output')
    args = parser.parse_args()
    if args.outputdir:
        if not os.path.exists(args.outputdir):
            try:
                os.mkdir(args.outputdir)
            except OSError:
                raise
        args.log = os.path.join(args.outputdir,'make_reference.log')
        args.samout = os.path.join(args.outputdir,'clustering.bam')
        args.consensus = os.path.join(args.outputdir,'consensus.fa')
        args.consensus_cluster = os.path.join(args.outputdir,'consensus_cluster.fa')
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
        while True:
            nextline = p.stdout.readline().replace('\r','\n')
            if nextline == '' and p.poll() != None:
                break
            log.write(nextline)
            log.flush()
        return_code = p.poll()
        stderr = p.stderr.read().replace('\r','\n')
        if return_code:
            raise RuntimeError(stderr)
        if stderr:
            log.write('stderr:\n%s\n'%stderr)
        log.write('finished:\t%s\n\n'%log_message)
    return 0

def merge_reads(args):
    "Unzip / Merge Watson and crick reads using pear"
    out_files = {}
    for strand in ['watson','crick']:
        fwd_out = tempfile.NamedTemporaryFile(suffix=".fastq.gz",prefix=strand,dir=args.tmpdir)
        rev_out = tempfile.NamedTemporaryFile(suffix=".fastq.gz",prefix=strand,dir=args.tmpdir)
        join_out = tempfile.NamedTemporaryFile(prefix="join_%s"%strand,dir=args.tmpdir)
        if args.sequences:
            head = '|head -n %s'%(int(args.sequences)*4)
        else:
            head = ''
        if strand == 'watson':
            grep_watson = "|grep Watson -A 3 |sed '/^--$/d'"
            cmd1 = ['gzcat '+args.forward + head + grep_watson + '|pigz -c >'+fwd_out.name]
            cmd2 = ['gzcat '+args.reverse + head + grep_watson + '|pigz -c >'+rev_out.name]
        else:
            grep_crick = "|grep Crick -A 3 |sed '/^--$/d'"
            cmd1 = ['gzcat '+args.forward + head + grep_crick + '|pigz -c  >'+fwd_out.name]
            cmd2 = ['gzcat '+args.reverse + head + grep_crick + '|pigz -c  >'+rev_out.name]
        log = "Write input files to tmpdir using gzcat"
        run_subprocess(cmd1,args,log)
        run_subprocess(cmd2,args,log)
        #todo: check if pear is on path
        pear = find_executable('pear')
        if not pear:
            raise EnvironmentError("Pear was not found, please make sure it exists on your path or install it from"+
                        " \nhttp://www.exelixis-lab.org/web/software/pear")
        cmd = [pear]
        #set reads
        cmd+=['-f',fwd_out.name]
        cmd+=['-r',rev_out.name]
        #set number of threads for pear
        cmd+=['-j','%s'%args.threads]
        cmd+=['-p','0.001']
        #set minimum overlap
        cmd+=['-v','10']
        #set minimum assembly length to 0
        cmd+=['-n','0']
        #specify output
        cmd+=['-o',join_out.name]
        log = "run pear for merging reads"
        run_subprocess(cmd,args,log)
        #Delete input files and output file name that are no longer needed??
        fwd_out.close()
        rev_out.close()
        fwd_out.close()
        #append output files as dictionary
        out_files[strand] = {'merged':join_out.name+".assembled.fastq",
                             'single_R1':join_out.name+".unassembled.forward.fastq",
                             'single_R2':join_out.name+".unassembled.reverse.fastq"}
    return out_files

def remove_methylation(in_files,args):
    """Remove methylation in watson and crick using sed"""
    for strand in ['watson','crick']:
        for key,value in in_files[strand].items():
            name_out = key + "_demethylated"
            file_out = value.split('.')[0]+ "_demethylated." + '.'.join(value.split('.')[1:])
            file_out += '.gz'
            if strand == 'watson':
                #sed only replaces values in lines that only contain valid nucleotides
                cmd = ["sed '/^[A,C,G,T,N]*$/s/C/T/g' "+value + "| pigz -c >"+file_out]
            else:
                #sed only replaces values in lines that only contain valid nucleotides
                cmd = ["sed '/^[A,C,G,T,N]*$/s/G/A/g' "+value + "| pigz -c >"+file_out]
            log = "use sed to remove methylation variation in %s_%s"%(strand,key)
            run_subprocess(cmd,args,log)
            in_files[strand][name_out] = file_out
    return in_files

def reverse_complement(read):
    """fast reverse complement"""
    nts = {'A':'T','C':'G','G':'C','T':'A','N':'N'}
    output = [nts[nt] for nt in read[::-1]]
    return ''.join(output)

def join_fastq(r1,r2,outfile):
    """join fastq files with 'NNNN' between forward and reverse complemented reverse read"""
    if r1.endswith('gz'):
        f1_handle = gzip.open(r1)
        f2_handle = gzip.open(r2)
    else:
        f1_handle = open(r1)
        f2_handle = open(r2)
    if outfile.endswith('gz'):
        out_handle = gzip.open(outfile,'w')
    else:
        out_handle = open(outfile,'w')
    while True:
        read1 = []
        read2 = []
        try:
            for i in range(4):
                read1.append(f1_handle.next())
                read2.append(f2_handle.next())
        except StopIteration:
            break
        out = '>%s\n%s'%(read1[0][1:-1],read1[1][:-1])
        #add N nucleotides to read s1
        out += 'N'*8
        #add reverse complement of read 2 to read 1
        #TODO: check if output file is correct in terms of strandedness
        out += read2[1][:-1]
        out += '\n'
        out_handle.write(out)
    out_handle.close()
    return True

def join_non_overlapping(in_files,args):
    """join non overlapping PE reads"""
    for strand in ['watson','crick']:
        #TODO: Trim end of R1 read length for demultiplexing? Different bc lengths
        out_file = in_files[strand]['single_R1_demethylated'].replace('unassembled.forward.fastq','joined.fa')
        join_fastq(in_files[strand]['single_R1_demethylated'],in_files[strand]['single_R2_demethylated'],out_file)
        #store output files in dictionary
        out_name = 'demethylated_joined'
        in_files[strand][out_name] = out_file
    return in_files



def trim_and_zip(in_files,args):
    """Trim fastq files and return using pigz"""
    in_files['trimmed'] = {}
    #Process single files
    #TODO: determine number of nucleotide trimmed from reads dependent on restriction enzyme being used.
    log = 'Process single watson reads: Trim first %s bases of R1'%'3'
    file_in = in_files['watson']['single_R1']
    if args.outputdir:
        file_out = os.path.join(args.outputdir, 'Unassembled.R1.watson_trimmed.fq.gz')
    else:
        file_out = '_'.join(args.watson_forward.split('_')[:-1])+'.unassembled.watson.R1.trimmed.fq.gz'
    in_files['trimmed']['watson_R1'] = file_out
    cmd = [seqtk + ' trimfq -b 4 %s |pigz -c > %s'%(file_in,file_out)]
    run_subprocess(cmd,args,log)

    log = 'Process single watson reads: reverse complement but no trimming required for R2 '
    file_in = in_files['watson']['single_R2']
    if args.outputdir:
        file_out = os.path.join(args.outputdir, 'Unassembled.R2.crick_trimmed.fq.gz')
    else:
        file_out = '_'.join(args.watson_reverse.split('_')[:-1])+'.unassembled.watson.R2.fq.gz'
    in_files['trimmed']['watson_R2'] = file_out
    #Take reverse complement as pear outputs R2 in reverse complement
    cmd = [seqtk + ' seq %s |%s seq -r - |pigz -c > %s'%(file_in,seqtk,file_out)]
    run_subprocess(cmd,args,log)

    log = 'Process single crick reads: no trimming required for R1'
    file_in = in_files['crick']['single_R1']
    if args.outputdir:
        file_out = os.path.join(args.outputdir, 'Unassembled.R1.watson_trimmed.fq.gz')
    else:
        file_out = '_'.join(args.crick_forward.split('_')[:-1])+'.unassembled.crick.R1.fq.gz'
    in_files['trimmed']['crick_R1'] = file_out
    cmd = [seqtk + ' seq %s |pigz -c >> %s'%(file_in,file_out)]
    run_subprocess(cmd,args,log)

    log = 'Process single crick reads: reverse complement an trim first 4 of R2'
    file_in = in_files['crick']['single_R2']
    if args.outputdir:
        file_out = os.path.join(args.outputdir, 'Unassembled.R2.crick_trimmed.fq.gz')
    else:
        file_out = '_'.join(args.crick_reverse.split('_')[:-1])+'.unassembled.crick.R2.trimmed.fq.gz'
    in_files['trimmed']['crick_R2'] = file_out
    #Take reverse complement as pear outputs R2 in reverse complement
    cmd = [seqtk + ' trimfq -e 4 %s |%s seq -r - |pigz -c >> %s'%(file_in,seqtk,file_out)]
    run_subprocess(cmd,args,log)

    #Process merged files
    log = 'Process merged watson reads: Trim first 4 bases of R1'
    file_in = in_files['watson']['merged']
    if args.outputdir:
        file_out = os.path.join(args.outputdir, 'Assembled.trimmed.fq.gz')
    else:
        file_out = '_'.join(args.watson_forward.split('_')[:-1])+'.assembled.watson.trimmedR1.fq.gz'
    in_files['trimmed']['watson_merged'] = file_out
    cmd = [seqtk + ' trimfq -b 4 %s |pigz -c > %s'%(file_in,file_out)]
    run_subprocess(cmd,args,log)

    log = 'Process merged crick reads: Trim first 4 bases of R2'
    file_in = in_files['crick']['merged']
    if args.outputdir:
        file_out = os.path.join(args.outputdir, 'Assembled.trimmed.fq.gz')
    else:
        file_out = '_'.join(args.crick_forward.split('_')[:-1])+'.assembled.crick.trimmedR2.fq.gz'
    in_files['trimmed']['crick_merged'] = file_out
    cmd = [seqtk + ' trimfq -e 4 %s |pigz -c >> %s'%(file_in,file_out)]
    run_subprocess(cmd,args,log)

    return in_files



def dereplicate_reads(in_files,args):
    """dereplicate reads using vsearch"""
    for strand in ['watson','crick']:
        for name in ["demethylated_joined","merged_demethylated"]:
            file_in = in_files[strand][name]
            file_out = '.'.join(file_in.split('.')[:-2]) + '.derep.' + file_in.split('.')[-2]
            in_files[strand][name] = [file_in]
            cmd = [vsearch+' -derep_fulllength %s -sizeout -minuniquesize 2 -output %s'%(file_in,file_out)]
            log = "Dereplicate full_length of %s using vsearch"%(strand)
            run_subprocess(cmd,args,log)
            name_out = name + "_derep"
            in_files[strand][name_out] = file_out
    return in_files

def combine_and_convert(in_files,args):
    """Combine and convert watson and crick files using sed"""
    in_files['combined'] = {}
    #First process merged files
    #"cat <(sed 's/>/>c/' atha_crick_derep.fa ) <(sed 's/>/>w/' atha_watson_derep.fa ) | tee >(sed 's/C/T/g;s/G/A/g' > atha_combined.CTGA.fa) > atha_combined.fa"
    #'/tmp/join_crickeBvj46_demethylated.assembled.derep.fastq'
    out_normal = in_files['crick']['merged_demethylated_derep'].replace('crick','combined')
    out_CTGA = out_normal.replace('derep','derep.CTGA')
    cmd = ["cat <(sed 's/>/>c/' %s ) <(sed 's/>/>w/' %s ) |"%(in_files['crick']['merged_demethylated_derep'],
            in_files['watson']['merged_demethylated_derep'])+
           "tee >(sed '/^[A,C,G,T,N]*$/s/C/T/g;/^[A,C,G,T,N]*$/s/G/A/g'"+
           "> %s) > %s"%(out_CTGA,out_normal)
            ]
    log = "Combine and convert merged watson and crick files using sed"
    run_subprocess(cmd,args,log)
    in_files['combined']["merged_combined"] = out_normal
    in_files['combined']["merged_combined_CTGA"] = out_CTGA

    #Now Process the joined files
    out_normal = in_files['crick']['demethylated_joined_derep'].replace('crick','combined')
    out_CTGA = out_normal.replace('derep','derep.CTGA')
    cmd = ["cat <(sed 's/>/>c/' %s ) <(sed 's/>/>w/' %s ) |"%(in_files['crick']['demethylated_joined_derep'],
            in_files['watson']['demethylated_joined_derep'])+
           "tee >(sed '/^[A,C,G,T,N]*$/s/C/T/g;/^[A,C,G,T,N]*$/s/G/A/g'"+
           "> %s) > %s"%(out_CTGA,out_normal)
            ]
    log = "Combine and convert joined watson and crick files using sed"
    run_subprocess(cmd,args,log)
    in_files['combined']["joined_combined"] = out_normal
    in_files['combined']["joined_combined_CTGA"] = out_CTGA
    return in_files

def make_uc(in_files,args):
    """run vsearch to create uc file"""
    #run vsearch to create uc file
    out_CTGA = in_files['combined']["merged_combined_CTGA"]
    merged_uc = tempfile.NamedTemporaryFile(suffix=".uc",prefix='derep_merged',dir=args.tmpdir,delete=False)
    joined_uc = tempfile.NamedTemporaryFile(suffix=".uc",prefix='derep_joined',dir=args.tmpdir,delete=False)
    cmd = [vsearch + ' -derep_fulllength %s -strand both -uc %s'%
    (out_CTGA,merged_uc.name)]
    log = "Make uc output for merged reads"
    run_subprocess(cmd,args,log)
    in_files['combined']["derep_merged"] = merged_uc.name
    #TODO: do the same processing on the joined files here.
    out_CTGA = in_files['combined']["joined_combined_CTGA"]
    cmd = [vsearch + ' -derep_fulllength %s -strand both -uc %s'%
    (out_CTGA,joined_uc.name)]
    log = "Make uc output for joined reads"
    run_subprocess(cmd,args,log)
    in_files['combined']["derep_joined"] = joined_uc.name
    return in_files

def make_sam(in_files,args):
    """Run mergeBSv3.py to make a sam file of clusters"""
    merged_sam = tempfile.NamedTemporaryFile(suffix=".sam",prefix='merged',dir=args.tmpdir)
    joined_sam = tempfile.NamedTemporaryFile(suffix=".sam",prefix='joined',dir=args.tmpdir)
    cmd1 = ['%s -s %s -c %s --clusters %s --samout %s'%
           (mergeBSv3,
            in_files['combined']['merged_combined'],
           in_files['crick']['merged_demethylated_derep'],
           in_files['combined']['derep_merged'],
           merged_sam.name)]
    cmd2 = ['%s -s %s -c %s --clusters %s --samout %s'%
           (mergeBSv3,
            in_files['combined']['joined_combined'],
           in_files['crick']['demethylated_joined_derep'],
           in_files['combined']['derep_joined'],
           joined_sam.name)]
    log = "Make sam output for merged reads"
    run_subprocess(cmd1,args,log)
    log = "Make uc output for joined reads"
    run_subprocess(cmd2,args,log)
    in_files['sam_out'] = {}
    in_files['sam_out']['merged'] = merged_sam.name
    in_files['sam_out']['joined'] = joined_sam.name
    log = "create sorted and indexed bam for next step"
    cmd1 = ["samtools view -Shb %s|samtools sort - %s;samtools index %s"%
           (merged_sam.name,merged_sam.name.replace('.sam',''),merged_sam.name.replace('.sam','.bam'))]
    cmd2 = ["samtools view -Shb %s|samtools sort - %s;samtools index %s"%
           (joined_sam.name,joined_sam.name.replace('.sam',''),joined_sam.name.replace('.sam','.bam'))]
    run_subprocess(cmd1,args,log)
    run_subprocess(cmd2,args,log)
    in_files['sam_out']['merged'] = merged_sam.name.replace('.sam','.bam')
    in_files['sam_out']['joined'] = joined_sam.name.replace('.sam','.bam')
    return in_files

def split_output(in_files,args):
    """Split output in indexed Watson and crick  bam files"""
    #Create bam file for watson and crick sam for both sam files
    watson_merged = tempfile.NamedTemporaryFile(suffix=".bam",prefix='watson_merged',dir=args.tmpdir,delete=False)
    crick_merged = tempfile.NamedTemporaryFile(suffix=".bam",prefix='crick_merged',dir=args.tmpdir,delete=False)
    watson_joined = tempfile.NamedTemporaryFile(suffix=".bam",prefix='watson_joined',dir=args.tmpdir,delete=False)
    crick_joined = tempfile.NamedTemporaryFile(suffix=".bam",prefix='crick_joined',dir=args.tmpdir,delete=False)
    cmd1 = ["samtools view -h %s | tee >(grep '^@\|watson' |"%(in_files['sam_out']['merged'])+
            "samtools view -Shb - |samtools sort -o - %s > %s )|"%
            (watson_merged.name.replace('.bam',''),watson_merged.name)+
            "grep '^@\|crick' |samtools view -Shb - |samtools sort -o - %s > %s"%
            (crick_merged.name.replace('.bam',''),crick_merged.name)]
    cmd2 = ["samtools view -h %s | tee >(grep '^@\|watson' |"%(in_files['sam_out']['joined'])+
            "samtools view -Shb - |samtools sort -o - %s > %s )|"%
            (watson_joined.name.replace('.bam',''),watson_joined.name)+
            "grep '^@\|crick' |samtools view -Shb - |samtools sort -o - %s > %s"%
            (crick_joined.name.replace('.bam',''),crick_joined.name)]
    log = "make watson sorted and indexed bam"
    run_subprocess(cmd1,args,log)
    log = "make crick sorted and indexed bam"
    run_subprocess(cmd2,args,log)
    in_files['sam_out']['watson_joined'] = watson_joined.name
    in_files['sam_out']['crick_joined'] = crick_joined.name
    in_files['sam_out']['watson_merged'] = watson_merged.name
    in_files['sam_out']['crick_merged'] = crick_merged.name
    return in_files

def get_consensus(in_files,args):
    """ 1. Make consensus using vcfutils.pl for both watson and crick
        2. Combine consensus to recreate original ref
        3. sort by length for subsequent processing"""
    in_files['consensus'] = {}
    cons_out = tempfile.NamedTemporaryFile(suffix=".fq",prefix='cons_out',dir=args.tmpdir,delete=False)
    for type in ['merged','joined']:
        crick_out = tempfile.NamedTemporaryFile(suffix=".fq",prefix='cns_crick',dir=args.tmpdir)
        watson_out = tempfile.NamedTemporaryFile(suffix=".fq",prefix='cns_watson',dir=args.tmpdir)
        consensus = tempfile.NamedTemporaryFile(suffix=".fq",prefix='consensus',dir=args.tmpdir)
        cons_sort = tempfile.NamedTemporaryFile(suffix=".fq",prefix='cons_sort',dir=args.tmpdir)
        #TODO: find way to make sure new samtools versions play well with this command trick..
        cmd1 = ['samtools_old mpileup -u %s | bcftools view -cg - | %s vcf2fq |%s seq -A - > %s'%
                (in_files['sam_out']['crick_%s'%type],vcfutils,seqtk,crick_out.name)]
        cmd2 = ['samtools_old mpileup -u %s | bcftools view -cg - | %s vcf2fq |%s seq -A -  > %s'%
               (in_files['sam_out']['watson_%s'%type],vcfutils,seqtk,watson_out.name)]
        log = "get consensus crick %s"%type
        run_subprocess(cmd1,args,log)
        log = "get consensus watson %s"%type
        run_subprocess(cmd2,args,log)
        in_files['consensus']['watson_%s'%(type)]=watson_out.name,
        in_files['consensus']['crick_%s'%(type)]=crick_out.name,
        in_files['consensus']['consensus_%s'%(type)]=cons_sort.name

        #combine watson and crick
        log = "combining watson and crick fastq for %s reads"%type
        cmd = ['%s -w %s -c %s -b %s -o %s'%
               (create_consensus,
                watson_out.name,
                crick_out.name,
                in_files['sam_out'][type],
                consensus.name
               )]
        run_subprocess(cmd,args,log)
        log = "Append consensus %s out to output"%type
        cmd = ["cat %s |tee >(wc -l) >> %s"%(consensus.name,cons_out.name)]
        run_subprocess(cmd,args,log)

    log = "sort combined consensus %s on length for vsearch"%type
    cmd = [vsearch + " -sortbylength %s --output %s"%(cons_out.name,args.consensus)]
    try:
        result = run_subprocess(cmd,args,log)
    except Exception:
        print Exception
        pass

    # log = "Write combined output to final destination"
    # cmd = ["cp %s %s"%(cons_out.name,args.consensus)]
    # run_subprocess(cmd,args,log)
    in_files['consensus']['consensus'] = args.consensus
    return in_files

def cluster_consensus(in_files,args):
    "Cluster concensus with preset id"
    cluster_cons = args.consensus_cluster
    cmd = [vsearch+" -cluster_smallmem %s -id 0.95 -centroids %s -sizeout -strand both"%
           (in_files['consensus']['consensus'],
           cluster_cons)]
    log = "Clustering consensus with 95% identity"
    run_subprocess(cmd,args,log)
    in_files['consensus']['consensus_clustered'] = cluster_cons
    log = "rename cluster_cons for bwa_meth compatibility"
    cluster_renamed = cluster_cons.replace('.fa','.renamed.fa')
    cmd = ['cat %s | rename_fast.py -n > %s'%(cluster_cons, cluster_renamed)]
    run_subprocess(cmd,args,log)
    return in_files

def check_dependencies():
    """check for presence of dependencies and if not present say where they can be installed"""

    return 0


def main():
    "Main function loop"
    #Check if os is windows, if so, throw error and tell the user that the software cannot run on windows.
    if os.name == 'nt':
        raise OSError("This pipeline relies on unix/linux with a bash shell.")
    check_dependencies()
    args = parse_args()
    #Make sure log is empty at start
    if os.path.isfile(args.log):
        os.remove(args.log)
    #Step 1: Merge Watson and crick reads returns dictionary
    files = merge_reads(args)
    #Step 2: Remove methylation variation from both watson and crick using sed
    files = remove_methylation(files,args)
    #Step 3: join the non overlapping PE reads from watson and crick using vsearch
    files = join_non_overlapping(files,args)
    #Step 3a use seqtk to trim merged and joined reads from enzyme recognition site
    files = trim_and_zip(files,args)
    #Step 4: Dereplicate all watson and crick reads
    files = dereplicate_reads(files,args)
    #Step 5 Combine crick and watson and create a AT only file
    files = combine_and_convert(files,args)
    #Step 6: create uc file
    files = make_uc(files,args)
    #Step 7: run mergeBSv3.py
    files = make_sam(files,args)
    #Step 8: Split output in indexed Watson and crick  bam files
    files = split_output(files,args)
    #Step 9:  Get consensus sequence from both bam files using mpileup,bcfview and vcfutils.pl
    #Step 9b: Combine sequences to get consensus
    files = get_consensus(files,args)
    #Step 10: cluster consensus with 95% identity and rename for bwameth compatibility.
    files = cluster_consensus(files,args)
    print 'boe'
if __name__ == '__main__':
    main()