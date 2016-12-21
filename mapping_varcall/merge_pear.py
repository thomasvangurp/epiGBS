#/usr/bin/env pypy
import subprocess
import tempfile
import argparse
import os
seqtk = "seqtk"

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
    parser.add_argument('-t','--tmpdir',
                        help='tmp directory',default='/tmp')
    parser.add_argument('--threads',
                        help='Number of threads to used where multithreading is possible')
    parser.add_argument('--outputdir',
                        help='Optional: output directory')
    parser.add_argument('--log',
                        help='log of output operation')
    args = parser.parse_args()
    if args.outputdir:
        if not os.path.exists(args.outputdir):
            try:
                os.mkdir(args.outputdir)
            except OSError:
                raise
        args.log = os.path.join(args.outputdir,'merge_reads.log')
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
    args = parse_args()
    #Make sure log is empty at start
    if not args.log:
        return 0
    if os.path.isfile(args.log):
        os.remove(args.log)
    #Step 1: Merge Watson and crick reads returns dictionary
    files = merge_reads(args)
    files = trim_and_zip(files, args)
    clear_tmp(files)
if __name__ == '__main__':
    main()