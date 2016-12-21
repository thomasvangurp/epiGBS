#!/usr/bin/env pypy
import argparse
import subprocess
import tempfile
import os
import gzip

def parse_args():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(description='Process input files')
    parser.add_argument('-s', '--sequences', type=str, default=None,
                        help='Input sequences to annotate')
    parser.add_argument('-d', '--db', type=str, default='database',
                        help='watson vcf for header')
    parser.add_argument('-x', '--xmloutput', type=str, default=None,
                        help='xml output with diamond results')
    parser.add_argument('-t', '--threads', type=str, default=None,
                        help='number of threads to use simultaneously')
    parser.add_argument('-s', '--sensitive', action='store_true',
                        help='number of threads to use simultaneously')
    parser.add_argument('-m', '--maxtargetseqs', type=str, default='20',
                        help='The maximum number of target sequences per query to keep alignments for')
    parser.add_argument('-e', '--evalue', type=str, default='0.0001',
                        help='Maximum expected value to keep an alignment.')
    parser.add_argument('-l', '--log', type=str, default=None,
                        help='log file')
    parser.add_argument('--tmpdir', type=str, default=tempfile.gettempdir(),
                        help='tmp dir, defaults to system tmp dir')
    args = parser.parse_args()
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

def diamond_command(args):
    """build and run diamond command"""
    #diamond blastx -k 20  --salltitles --sensitive -f 6 -p 8
    # --db /Volumes/BACKUP_THOM/sequenceserver/blast_database/plant.refseq.03.11.2016.prot.fa.diamond.db.dmnd
    # -t /tmp
    # -q /Volumes/BACKUP_THOM/Zwitserland/Ono_vic/output_denovo/consensus_cluster.renamed.fa
    # -o /Volumes/BACKUP_THOM/Zwitserland/Ono_vic/output_denovo/consensus_cluster.renamed.diamond.xml
    cmd = ['diamond','blastx','-k %s' % args.maxtargetseqs]
    cmd.append('--salltitles') #add full description to blast results, for pretty blast2go output
    if args.sensitive:
        cmd.append('--sensitive')
    cmd.append('-j %s' % args.threads)
    cmd.append('-f 5') #xml output
    cmd.append('-db %s' % args.db) #xml output
    cmd.append('-t %s' % args.tmpdir) #temp directory
    cmd.append('-q %s' % args.sequences) #query, input fasta
    cmd.append('-o %s' % args.xmloutput) #xml output
    if args.xmloutput.endswith('.gz'):
        cmd.append('-compress 1') #enable gzip compression
    log = 'run diamond'


def main():

    return 0


if __name__ == '__main__':
    main()
