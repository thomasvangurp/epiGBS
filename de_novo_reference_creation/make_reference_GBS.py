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
import gzip

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
    "Merged reads using PEAR"
    out_files = {}
    if args.sequences:
        head = '|head -n %s'%(int(args.sequences)*4)
    else:
        head = ''
    if args.forward.endswith('.gz'):
        cat = 'pigz -p %s -cd '%args.threads
    else:
        cat = 'cat '
    #todo: check if pear is on path
    cmd = ['pear']
    #set reads
    cmd+=['-f',args.forward]
    cmd+=['-r',args.reverse]
    #set number of threads for pear
    cmd+=['-j','%s'%args.threads]
    cmd+=['-p','0.001']
    #set minimum overlap
    cmd+=['-v','10']
    #set minimum assembly length to 0
    cmd+=['-n','0']
    #specify output
    cmd+=['-o','%s/merged'%args.outputdir]
    log = "run pear for merging reads"

    run_subprocess(cmd,args,log)
    #Delete input files and output file name that are no longer needed??
    #append output files as dictionary
    out_files = {'merged':'%smerged'%args.outputdir + ".assembled.fastq",
                         'single_R1':'%smerged'%args.outputdir + ".unassembled.forward.fastq",
                         'single_R2':'%smerged'%args.outputdir + ".unassembled.reverse.fastq"}
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
    cmd = ["paste <(seqtk seq -A %s | cut -c1-%s) " % (r1, max_len_R1) +
           "<(seqtk seq  -A %s |cut -c1-%s)|cut -f1-5" % (r2, max_len_R2)+
           "|sed '/^>/!s/\t/NNNNNNNN/g' |pigz -p %s -c > %s" % (args.threads, outfile)]
    log = "Combine joined fastq file into single fasta file"
    if not os.path.exists(outfile):
        run_subprocess(cmd,args,log)
    return True

def join_non_overlapping(in_files,mapping_dict, args):
    """join non overlapping PE reads"""
    for key in mapping_dict:
        key = key.replace(' ','_')
        if 'mono' not in key and key != 'all':
            continue
        left = os.path.join(args.outputdir, '%s.R1.fq.gz' % key)
        right = os.path.join(args.outputdir, '%s.R2.fq.gz' % key)
        joined = os.path.join(args.outputdir, '%s.joined.fq.gz' % key)
        if not os.path.exists(joined):
            join_fastq(left, right, joined, args)
    return in_files

def get_mapping_dict(args):
    """get mapping dict for mono's"""
    mapping_dict = {'all':None}
    with open(args.barcodes) as barcode_handle:
        header = barcode_handle.readline()[:-1].split('\t')
        for line in barcode_handle:
            split_line = line[:-1].split('\t')
            try:
                id = split_line[header.index('Sample')]
            except KeyError:
                raise KeyError('barcode file does not contain Sample, please revise your file %s' % args.barcodes)
            mapping_dict[id] = {}
            for k, v in zip(header, split_line):
                if k != 'Sample':
                    mapping_dict[id][k] = v
    return mapping_dict

def trim_split_and_zip(in_files, mapping_dict, args):
    """Trim , split and zip fastq files for mono species"""
    in_files['trimmed'] = {}

    # log = 'Zip and split forward reads'
    # file_in = in_files['single_R1']
    # cmd = seqtk + ' seq %s |tee >(pigz -c > %s)' % (file_in, os.path.join(args.outputdir,'all.R1.fq.gz'))
    # mono_list = [k for k in mapping_dict.keys() if 'mono' in k.lower()]
    # for k in mono_list[:-1]:
    #     cmd += "|tee >(grep '%s' -A 3|sed '/^--$/d' |pigz -c > %s)" % (k, os.path.join(args.outputdir,'%s.R1.fq.gz'%(k.replace(' ','_'))))
    # #add the latest key to mono 1.
    # k = mono_list[-1]
    # cmd += "|grep '%s' -A 3|sed '/^--$/d' |pigz -c > %s" % (k, os.path.join(args.outputdir,'%s.R1.fq.gz'%k.replace(' ','_')))
    # if not os.path.exists(os.path.join(args.outputdir,'all.R1.fq.gz')):
    #     run_subprocess([cmd],args,log)
    #
    # log = 'Zip and split reverse reads'
    # file_in = in_files['single_R2']
    # cmd = seqtk + ' seq -r %s |tee >(pigz -c > %s)' % (file_in, os.path.join(args.outputdir, 'all.R2.fq.gz'))
    # mono_list = [k for k in mapping_dict.keys() if 'mono' in k.lower()]
    # for k in mono_list[:-1]:
    #     cmd += "|tee >(grep '%s' -A 3|sed '/^--$/d' |pigz -c > %s)" % (
    #     k, os.path.join(args.outputdir, '%s.R2.fq.gz' % (k.replace(' ', '_'))))
    # # add the latest key to mono 1.
    # k = mono_list[-1]
    # cmd += "|grep '%s' -A 3|sed '/^--$/d' |pigz -c > %s" % (
    # k, os.path.join(args.outputdir, '%s.R2.fq.gz' % k.replace(' ', '_')))
    # if not os.path.exists(os.path.join(args.outputdir, 'all.R2.fq.gz')):
    #     run_subprocess([cmd], args, log)

    #Process merged files
    log = 'Zip and split merged reads'
    file_in = in_files['merged']
    cmd = seqtk + ' seq %s |pigz -c > %s' % (file_in, os.path.join(args.outputdir, 'all.merged.fq.gz'))
    # mono_list = [k for k in mapping_dict.keys() if 'mono' in k.lower()]
    # for k in mono_list[:-1]:
    #     cmd += "|tee >(grep '%s' -A 3|sed '/^--$/d' |pigz -c > %s)" % (
    #     k, os.path.join(args.outputdir, '%s.merged.fq.gz' % (k.replace(' ', '_'))))
    # add the latest key to mono 1.
    # k = mono_list[-1]
    # cmd += "|grep '%s' -A 3|sed '/^--$/d' |pigz -c > %s" % (
    # k, os.path.join(args.outputdir, '%s.merged.fq.gz' % k.replace(' ', '_')))
    if not os.path.exists(os.path.join(args.outputdir, 'all.merged.fq.gz')):
        run_subprocess([cmd], args, log)

    return in_files



def dereplicate_reads(in_files, mapping_dict, args):
    """dereplicate reads using vsearch"""
    in_files['derep'] = {}
    for type in ['joined','merged']:
        for name in [k for k in mapping_dict.keys() if 'mono' in k]:
            name = name.replace(' ','_')
            file_in = os.path.join(args.outputdir, '%s.%s.fq.gz' % (name, type))
            file_out = '.'.join(file_in.split('.')[:-2]) + '.derep.fa'
            cmd = [vsearch +' -derep_fulllength %s -sizeout -minuniquesize 2 -output %s'%(file_in,file_out)]
            log = "Dereplicate full_length of %s using vsearch"%(name)
            if not os.path.exists(file_out):
                run_subprocess(cmd, args, log)
            in_files['derep']['%s_%s_derep' % (name, type)] = file_out
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



def cluster_consensus(in_files, mapping_dict, args):
    "Cluster concensus with preset id"
    #concatenate and order derep output
    for name in [k for k in mapping_dict.keys() if 'mono' in k]:
        #concatenate and order
        name = name.replace(' ','_')
        joined = os.path.join(args.outputdir,'%s.joined.derep.fa' % name)
        merged = os.path.join(args.outputdir,'%s.merged.derep.fa' % name)
        combined = os.path.join(args.outputdir,'%s.combined.derep.fa' % name)
        combined_ordered = os.path.join(args.outputdir,'%s.combined.ordered.derep.fa' % name)
        output_derep = os.path.join(args.outputdir,'%s.clustered.fa' % name)
        output_derep_renamed = os.path.join(args.outputdir,'%s.clustered.renamed.fa' % name)
        cmd = ['cat %s %s > %s' % (joined, merged, combined)]
        log = "combine joined and merged dereplicates reads of %s for sorting" % name
        if not os.path.exists(combined):
            run_subprocess(cmd, args, log)


        cmd = ['vsearch -sortbylength %s --output %s' % (combined, combined_ordered)]
        log = 'Sort derep sequences by length.'
        if not os.path.exists(combined_ordered):
            run_subprocess(cmd, args, log)

        cmd = [vsearch+" -cluster_smallmem %s -id 0.95 -centroids %s -sizeout -strand both"%
               (combined_ordered,output_derep
                )]
        log = "Clustering consensus with 95% identity"
        if not os.path.exists(output_derep):
            run_subprocess(cmd,args,log)
        # in_files['consensus']['consensus_clustered'] = args.consensus_cluster

        cmd = ['cat %s |rename_fast.py -g %s -n > %s'%(output_derep, name, output_derep_renamed)]
        log = "rename resulting clusters with number and origin name"
        if not os.path.exists(output_derep_renamed):
            run_subprocess(cmd,args,log)
            cmd = ["samtools faidx %s"%output_derep_renamed]
            log = "faidx index %s" % output_derep_renamed
            run_subprocess(cmd, args, log)
    #concatenate all renamed clusters
    ref = os.path.join(args.outputdir, 'ref.fa')
    cmd = ['cat %s/*.renamed.fa > %s' % (args.outputdir, ref)]
    log = 'concatenate all renamed clusters'
    run_subprocess(cmd, args, log)
    #index this new reference sequence
    cmd = ["samtools faidx %s" % ref]
    log = "faidx index %s" % ref
    run_subprocess(cmd, args, log)
    return in_files

def check_dependencies():
    """check for presence of dependencies and if not present say where they can be installed"""

    return 0

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
    check_dependencies()
    args = parse_args()
    #Make sure log is empty at start
    if os.path.isfile(args.log):
        os.remove(args.log)
    #Step 1: get mapping dict for mono's
    mapping_dict = get_mapping_dict(args)
    #Step 1: Merge Watson and crick reads returns dictionary
    #Step 3a use seqtk to trim merged and joined reads from enzyme recognition site
    files = merge_reads(args)
    #Step 4: Dereplicate all watson and crick reads
    files = trim_split_and_zip(files, mapping_dict, args)
    #Step 3: join the non overlapping PE reads from watson and crick using vsearch
    files = join_non_overlapping(files, mapping_dict, args)
    files = dereplicate_reads(files, mapping_dict, args)
    #Step 5 create binary A/T only fasta output, with headers containing original sequence
    #step 8: Cluster consensus
    files = cluster_consensus(files, mapping_dict, args)
    cmd = 'rm %s/mono*' % args.outputdir
    log = 'removing all intermeditate files that are no longer needed'
    run_subprocess(cmd, args, log)

if __name__ == '__main__':
    main()