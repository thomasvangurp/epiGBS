#!/usr/bin/env pypy
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
import operator
import unittest
import gzip

#unittest for robust modules, for in installation mode help to create robust function
#integration testing
origWD = os.getcwd()
os.chdir(origWD)


dependencies = {}
dependencies['samtools'] = '<0.1.18'
dependencies['vcfutils.pl'] = ''
dependencies['usearch'] = 'usearch_8.0.1409'
dependencies['seqtk'] = '1.0-r31'
dependencies['pear'] = 'v0.9.7'
dependencies['pigz'] = ''

usearch = "usearch"#_8.0.1409_i86osx32"
vsearch = "vsearch"
seqtk = "seqtk"
pear = "pear"
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
    parser.add_argument('--barcodes',
                        help='max barcode length used to trim joined reads')
    parser.add_argument('--cycles',default='126',
                        help='Number of sequencing cycles / read length')
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
    assert os.path.exists(args.barcodes)
    try:
        assert os.path.exists(args.forward)
    except AssertionError:
        raise AssertionError("%s does not exist" % args.forward)
    try:
        assert os.path.exists(args.reverse)
    except AssertionError:
        raise AssertionError("%s does not exist" % args.reverse)
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
        p = subprocess.Popen(cmd,stdout=subprocess.PIPE,stderr=subprocess.PIPE,shell=True,executable='bash')
        stdout, stderr = p.communicate()
        stdout = stdout.replace('\r','\n')
        stderr = stderr.replace('\r','\n')
        if stdout:
            log.write('stdout:\n%s\n'%stdout)
        if stderr:
            log.write('stderr:\n%s\n'%stderr)
        log.write('finished:\t%s\n\n'%log_message)
    return 0

def merge_reads(args):
    "Unzip / Merge Watson and crick reads using pear"
    #TODO: run only once for both watson and crick at same time
    out_files = {}
    for strand in ['watson','crick']:
        fwd_out = tempfile.NamedTemporaryFile(suffix=".fastq.gz",prefix=strand,dir=args.tmpdir)
        rev_out = tempfile.NamedTemporaryFile(suffix=".fastq.gz",prefix=strand,dir=args.tmpdir)
        join_out = tempfile.NamedTemporaryFile(prefix="join_%s"%strand,dir=args.tmpdir)
        if args.sequences:
            head = '|head -n %s'%(int(args.sequences)*4)
        else:
            head = ''
        if args.forward.endswith('.gz'):
            cat = 'pigz -p %s -cd '%args.threads
        else:
            cat = 'cat '
        if strand == 'watson':
            grep_watson = "|grep 'Watson\|watson' -A 3 |sed '/^--$/d'"
            cmd1 = [ cat + args.forward + head + grep_watson + '|pigz -p %s -c >'%(args.threads)+fwd_out.name]
            cmd2 = [ cat  + args.reverse + head + grep_watson + '|pigz -p %s -c >'%(args.threads)+rev_out.name]
        else:
            grep_crick = "|grep 'Crick\|crick' -A 3 |sed '/^--$/d'"
            cmd1 = [cat + args.forward + head + grep_crick + '|pigz -p %s -c  >'%(args.threads)+fwd_out.name]
            cmd2 = [cat + args.reverse + head + grep_crick + '|pigz -p %s -c  >'%(args.threads)+rev_out.name]
        log = "Write input files to tmpdir using gzcat"
        run_subprocess(cmd1,args,log)
        run_subprocess(cmd2,args,log)
        #todo: check if pear is on path
        cmd = ['pear']
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
                cmd = ["sed '/^[A,C,G,T,N]*$/s/C/T/g' "+value + "| pigz -p %s -c >"%(args.threads)+file_out]
            else:
                #sed only replaces values in lines that only contain valid nucleotides
                cmd = ["sed '/^[A,C,G,T,N]*$/s/G/A/g' "+value + "| pigz -p %s -c >"%(args.threads)+file_out]
            log = "use sed to remove methylation variation in %s_%s"%(strand,key)
            run_subprocess(cmd,args,log)
            in_files[strand][name_out] = file_out
    return in_files

def reverse_complement(read):
    """fast reverse complement"""
    nts = {'A':'T','C':'G','G':'C','T':'A','N':'N'}
    output = [nts[nt] for nt in read[::-1]]
    return ''.join(output)

def join_fastq(r1,r2,outfile,args):
    """join fastq files with 'NNNN' between forward and reverse complemented reverse read"""
    #get max length of forward and reverse barcodes
    if args.barcodes:
        with open(args.barcodes) as bc_handle:
            header = bc_handle.readline()[:-1].split('\t')
            barcode_1_index = header.index('Barcode_R1')
            barcode_2_index = header.index('Barcode_R2')
            try:
                wobble_R1_index = header.index('Wobble_R1')
                wobble_R2_index = header.index('Wobble_R2')
            except ValueError:
                wobble_R1_index = None
                wobble_R2_index = None
            barcode_1_max_len = 0
            barcode_2_max_len = 0
            for line in bc_handle:
                split_line = line.rstrip('\n').split('\t')
                try:
                    #TODO: make control nucleotide explicit option in barcode file, now harccoded!
                    wobble_R1_len = int(split_line[wobble_R1_index]) + 1
                    wobble_R2_len = int(split_line[wobble_R2_index]) + 1
                except TypeError:
                    wobble_R1_len = 0
                    wobble_R2_len = 0
                if len(split_line[barcode_1_index]) > barcode_1_max_len:
                    barcode_1_max_len = len(split_line[barcode_1_index])
                if len(split_line[barcode_2_index]) > barcode_2_max_len:
                    barcode_2_max_len = len(split_line[barcode_2_index])
        max_len_R1 = int(args.cycles) - barcode_1_max_len - wobble_R1_len
        max_len_R2 = int(args.cycles) - barcode_2_max_len - wobble_R2_len
    else:
        #no trimming required
        max_len_R1 = 200
        max_len_R2 = 200
    #Trim the reads up to the min expected length to improve de novo reference creation for joined reads
    cmd = ["paste <(pigz -p %s -cd %s |seqtk seq -A - | cut -c1-%s) " % (args.threads, r1, max_len_R1) +
           "<(pigz -p %s -cd %s |seqtk seq -A -|cut -c1-%s)|cut -f1-5" % (args.threads, r2, max_len_R2)+
           "|sed '/^>/!s/\t/NNNNNNNN/g' |pigz -p %s -c > %s" % (args.threads, outfile)]
    log = "Combine joined fastq file"
    run_subprocess(cmd,args,log)
    return True

def join_non_overlapping(in_files,args):
    """join non overlapping PE reads"""
    for strand in ['watson','crick']:
        #TODO: Trim end of R1 read length for demultiplexing? Different bc lengths
        out_file = in_files[strand]['single_R1_demethylated'].replace('unassembled.forward.fastq','joined.fa')
        join_fastq(in_files[strand]['single_R1_demethylated'],in_files[strand]['single_R2_demethylated'],out_file,args)
        #store output files in dictionary
        out_name = 'demethylated_joined'
        in_files[strand][out_name] = out_file
    return in_files



def trim_and_zip(in_files,args):
    """Trim fastq files and return using pigz"""
    in_files['trimmed'] = {}

    log = 'Zip single watson reads: '
    file_in = in_files['watson']['single_R1']
    if args.outputdir:
        file_out = os.path.join(args.outputdir, 'Unassembled.R1.watson.fq.gz')
    else:
        file_out = '_'.join(args.watson_forward.split('_')[:-1])+'.Unassembled.watson.R1.fq.gz'
    in_files['trimmed']['watson_R1'] = file_out
    cmd = [seqtk + ' seq %s |pigz -p %s -c > %s'%(file_in, args.threads, file_out)]
    run_subprocess(cmd,args,log)

    log = 'Process single watson reads: reverse complement required for R2 '
    file_in = in_files['watson']['single_R2']
    if args.outputdir:
        file_out = os.path.join(args.outputdir, 'Unassembled.R2.crick.fq.gz')
    else:
        file_out = '_'.join(args.watson_reverse.split('_')[:-1])+'.Unassembled.watson.R2.fq.gz'
    in_files['trimmed']['watson_R2'] = file_out
    #Take reverse complement as pear outputs R2 in reverse complement
    cmd = [seqtk + ' seq -r %s |pigz -p %s -c > %s'%(file_in, args.threads, file_out)]
    run_subprocess(cmd,args,log)

    log = 'Process single crick reads: no trimming required for R1'
    file_in = in_files['crick']['single_R1']
    if args.outputdir:
        file_out = os.path.join(args.outputdir, 'Unassembled.R1.watson.fq.gz')
    else:
        file_out = '_'.join(args.crick_forward.split('_')[:-1])+'.Unassembled.crick.R1.fq.gz'
    in_files['trimmed']['crick_R1'] = file_out
    cmd = [seqtk + ' seq %s |pigz -p %s -c >> %s'%(file_in, args.threads, file_out)]
    run_subprocess(cmd,args,log)

    log = 'Process single crick reads: reverse complement for R2'
    file_in = in_files['crick']['single_R2']
    if args.outputdir:
        file_out = os.path.join(args.outputdir, 'Unassembled.R2.crick.fq.gz')
    else:
        file_out = '_'.join(args.crick_reverse.split('_')[:-1])+'.Unassembled.crick.R2.fq.gz'
    in_files['trimmed']['crick_R2'] = file_out
    #Take reverse complement as pear outputs R2 in reverse complement
    cmd = [seqtk + ' seq -r %s |pigz -p %s -c >> %s'%(file_in, args.threads, file_out)]
    run_subprocess(cmd,args,log)

    #Process merged files
    log = 'Process merged watson reads:'
    file_in = in_files['watson']['merged']
    if args.outputdir:
        file_out = os.path.join(args.outputdir, 'Assembled.fq.gz')
    else:
        file_out = '_'.join(args.watson_forward.split('_')[:-1])+'.Assembled.R1.fq.gz'
    in_files['trimmed']['watson_merged'] = file_out
    cmd = [seqtk + ' seq %s |pigz -p %s -c > %s'%(file_in, args.threads, file_out)]
    run_subprocess(cmd,args,log)

    log = 'Process merged crick reads:'
    file_in = in_files['crick']['merged']
    if args.outputdir:
        file_out = os.path.join(args.outputdir, 'Assembled.fq.gz')
    else:
        file_out = '_'.join(args.crick_forward.split('_')[:-1])+'.Assembled.fq.gz'
    in_files['trimmed']['crick_merged'] = file_out
    cmd = [seqtk + ' seq %s |pigz -p %s -c >> %s'%(file_in, args.threads, file_out)]
    run_subprocess(cmd,args,log)

    return in_files



def dereplicate_reads(in_files,args):
    """dereplicate reads using vsearch"""
    for strand in ['watson','crick']:
        for name in ["demethylated_joined","merged_demethylated"]:
            file_in = in_files[strand][name]
            file_out = '.'.join(file_in.split('.')[:-2]) + '.derep.' + file_in.split('.')[-2]
            in_files[strand][name] = [file_in]
            cmd = [vsearch +' -derep_fulllength %s -sizeout -minuniquesize 2 -output %s'%(file_in,file_out)]
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

def make_binary_output(in_files,args):
    """make binary output sequence for uc generation"""
    #TODO: implement this routine instead of combine_and_convert
    watson_handle = open(in_files['watson']['merged_demethylated_derep'],'r')
    crick_handle = open(in_files['crick']['merged_demethylated_derep'],'r')
    #joined entries
    watson_handle_join = open(in_files['watson']['demethylated_joined_derep'],'r')
    crick_handle_join = open(in_files['crick']['demethylated_joined_derep'],'r')
    uc_input_fa = tempfile.NamedTemporaryFile(suffix=".fa", prefix='uc_input', dir=args.tmpdir, delete=False)
    in_files['combined'] = {'uc_input_fa':uc_input_fa.name}
    output = open(in_files['combined']['uc_input_fa'],'w')
    for i,handle in enumerate([watson_handle,crick_handle,watson_handle_join,crick_handle_join]):
        if i == 0 or i == 2:
            name_start = '>w_'
        else:
            name_start = '>c_'
        for line in handle:
            if line.startswith('>'):
                try:
                    seq = seq.replace('C','T').replace('G','A')
                    out = name + '\n' + seq + '\n'
                    output.write(out)
                except NameError:
                    pass
                name = '%s'%name_start
                seq = ''
            else:
                name += line.rstrip('\n')
                seq += line.rstrip('\n')
    return in_files

def make_uc(in_files,args):
    """run usearch to create uc file"""
    #run usearch to create uc file
    fasta_input = in_files['combined']['uc_input_fa']
    uc_out = tempfile.NamedTemporaryFile(suffix=".uc",prefix='combined_out',dir=args.tmpdir,delete=False)
    cmd = [vsearch + ' -derep_fulllength %s -strand both -uc %s'%(fasta_input,uc_out.name)]
    log = "Make uc output for reads with original sequences in read header"
    run_subprocess(cmd,args,log)
    in_files['combined']["uc_out"] = uc_out.name
    return in_files


def get_ref(clusters):
    """Generate reference from sequences that cluster together"""
    output_count = {}
    id = int(clusters[0][1]) + 1
    transform = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    for cluster in clusters:
        direction = cluster[4]
        seq = cluster[8][2:]
        strand = cluster[8][0]
        if direction != '-':
            pass
        else:
            #get reverse complement
            seq = [transform[nt] for nt in seq][::-1]
            if strand == 'w':
                strand = 'c'
            else:
                strand = 'w'
        for i,nt in enumerate(seq):
            try:
                output_count[i][strand][nt] += 1
            except KeyError:
                if i not in output_count:
                    output_count[i] = {strand:{nt:1}}
                if strand not in output_count[i]:
                    output_count[i][strand] = {nt:1}
                if nt not in output_count[i][strand]:
                    output_count[i][strand][nt] = 1
    output_fasta = '>%s\n'%id
    for i in sorted(output_count.keys()):
        try:
            watson_nt = max(output_count[i]['w'].iteritems(), key=operator.itemgetter(1))[0]
            crick_nt = max(output_count[i]['c'].iteritems(), key=operator.itemgetter(1))[0]
        except KeyError:
            return None
        if watson_nt == crick_nt:
            output_fasta += watson_nt
        elif watson_nt == 'G' and crick_nt == 'A':
            output_fasta += watson_nt
        elif crick_nt == 'C' and watson_nt == 'T':
            output_fasta += crick_nt
        elif watson_nt == 'N' and crick_nt != 'N':
            output_fasta += crick_nt
        elif watson_nt != 'N' and crick_nt == 'N':
            output_fasta += watson_nt
        else:
            output_fasta += 'N'
    output_fasta += '\n'
    return output_fasta
    return output_fasta



def make_ref_from_uc(in_files,args):
    """make reference directly from uc output"""
    #ref_handle contains non-clustered output sequences
    uc_handle = open(in_files["combined"]["uc_out"],'r')
    ref_output = tempfile.NamedTemporaryFile(suffix=".fa", prefix='tmp_ref', dir=args.tmpdir, delete=True)
    ref_handle = open(ref_output.name,'w')
    clusters = []
    for line in uc_handle:
        if line.startswith('H') or line.startswith('S'):
            split_line = line.split('\t')
            cluster_id = int(split_line[1])
            try:
                if cluster_id != int(clusters[-1][1]):
                    ref = get_ref(clusters)
                    if ref:
                        ref_handle.write(ref)
                    clusters = []
            except NameError:
                pass
            except IndexError:
                pass
            clusters.append(split_line)
    ref = get_ref(clusters)
    if ref:
        ref_handle.write(ref)
    ref_handle.close()

    #sort output by length
    cmd = ['vsearch -sortbylength %s --output %s'%(ref_output.name,args.consensus)]
    log = "sort sequences by length"
    run_subprocess(cmd, args, log)
    return in_files



def cluster_consensus(in_files,args):
    "Cluster concensus with preset id"
    #TODO: vsearch ignores joined sequences, temporarily revert back to vsearch to prevent this from happening.
    cmd = [vsearch + "-qmask none -cluster_smallmem %s -id 0.95 -centroids %s -sizeout -strand both"%
           (args.consensus,
            args.consensus_cluster)]
    log = "Clustering consensus with 95% identity"
    run_subprocess(cmd,args,log)
    # in_files['consensus']['consensus_clustered'] = args.consensus_cluster
    log = "rename cluster_cons for bwa_meth compatibility"
    cluster_renamed = args.consensus_cluster.replace('.fa','.renamed.fa')
    cmd = ['cat %s | rename_fast.py -n > %s'%(args.consensus_cluster, cluster_renamed)]
    run_subprocess(cmd,args,log)
    log = "faidx index %s" % cluster_renamed
    cmd = ["samtools faidx %s" % cluster_renamed]
    run_subprocess(cmd, args, log)
    return in_files

def check_dependencies():
    """check for presence of dependencies and if not present say where they can be installed"""

    return 0

def add_numbers(a,b):
    """"return sum of a and b"""



def clear_tmp(file_dict):
    """clear tmp files"""
    purge_list = []
    for v in file_dict.keys():
        for key,value in file_dict[v].items():
            try:
                if value.startswith('/tmp'):
                    purge_list.append(value)
            except AttributeError:
                if type(value) == type([]):
                    purge_list.append(value[0])
    for item in purge_list:
        print "removing %s" % item
        os.remove(item)
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
    #Step 5 create binary A/T only fasta output, with headers containing original sequence
    files = make_binary_output(files,args)
    #Step 6: create uc file
    files = make_uc(files,args)
    #Step 7: make reference directly from single uc file
    files = make_ref_from_uc(files,args)
    #step 8: Cluster consensus
    files = cluster_consensus(files,args)
    # step 8: Clean tmp files
    files = clear_tmp(files)

if __name__ == '__main__':
    main()