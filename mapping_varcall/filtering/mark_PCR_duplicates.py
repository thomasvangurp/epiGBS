#!/usr/bin/env python
import pysam
import os
import sys
import tempfile
import subprocess
import argparse
from Bio import SeqIO
from Bio import Restriction

def parse_args():
    """Pass command line arguments"""
    parser = argparse.ArgumentParser(description='Remove PCR duplicates from bam file')
    #input files
    parser.add_argument('-i', '--input',
                        help='bam input file')
    parser.add_argument('-o', '--output',
                        help='output bam file')
    parser.add_argument('-b', '--barcodes',
                        help='barcodes to detect enzyme combination used')
    parser.add_argument('--input_dir',
                        help='input directory, all bam file will be filtered for PCR duplicates')
    parser.add_argument('-l', '--log',
                        help='log, default is stdout',default=sys.stdout)
    parser.add_argument('--stat',
                        help='statistics on PCR duplicate removal')
    parser.add_argument('-r','--reference',
                    help='reference clusters')
    args = parser.parse_args()
    if args.input_dir:
        args.watson = os.path.join(args.input_dir,'watson.bam')
        if not os.path.exists(args.watson):
            raise OSError("%s does not exist" % args.watson)
        args.crick = os.path.join(args.input_dir,'crick.bam')
        if not os.path.exists(args.crick):
            raise OSError("%s does not exist" % args.crick)
    return args

#TODO: make statistics output per contig. (Do lowly covered contigs get less PCR duplicate rates?

def get_enz(enz):
    """Get enzyme from biopython restriction library"""
    for enzyme in Restriction.AllEnzymes:
        if "%s"%(enzyme) == enz:
            return enzyme

def get_regions(contig,enzymes):
    """return loci with start and end locations"""
    out_sites = []
    enz_1 = [enz for enz in Restriction.AllEnzymes if "%s"%enz == enzymes[0]][0]
    enz_2 = [enz for enz in Restriction.AllEnzymes if "%s"%enz == enzymes[1]][0]
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



def remove_PCR_duplicates(bam_in, bam_out, ref):
    """Remove PCR duplicates and non-paired PE-reads per cluster"""
    #TODO: implement as module to be ran concurrently with appropriate statistics.
    clusters = SeqIO.parse(open(ref),'fasta')
    handle = pysam.AlignmentFile(bam_in,'rb')
    out_handle = pysam.AlignmentFile(bam_out,'wb', template=handle)
    read_count = {}
    for cluster in clusters:
        #TODO: detect enzyme from barcode file to enable other combinations.
        enzymes = ["Csp6I","NsiI"]
        cluster_seq_len = len(str(cluster.seq).upper().replace('N',''))
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
            # if 'NNNNNNNN' in cluster._seq.upper() and not region:
            #     cluster_is_paired = True
            # elif region:
            #     if region[1] - region[0] > 240:
            #         cluster_is_paired = True
            #     else:
            #         cluster_is_paired = False
            # else:
            #     cluster_is_paired = False
            read_out = {}
            read = None
            qc_fail = set()
            for read in reads:
                tag_dict = dict(read.tags)
                try:
                    tag = tag_dict['RN']
                    sample = tag_dict['RG']
                    AS = tag_dict['AS']
                except KeyError:
                    break
                if sample not in read_out:
                    read_out[sample] = {}
                if tag not in read_out[sample]:
                    read_out[sample][tag] = {read.qname:AS}
                else:
                    try:
                        read_out[sample][tag][read.qname]+= AS
                    except KeyError:
                        read_out[sample][tag][read.qname] = AS
                if tag_dict['AS'] < 0.9 * cluster_seq_len:
                    #alignment score lower than read length, subset read
                    read.is_qcfail = True
                    #update the 'failed' dictionary to exclude this read
                    qc_fail.update([read.qname])
                    try:
                        read_count[sample]['qc_fail'] += 1
                    except KeyError:
                        if sample not in read_count:
                            read_count[sample] = {}
                        read_count[sample]['qc_fail'] = 1
                if read.is_paired and not read.is_proper_pair:
                    read.is_qcfail = True
                    qc_fail.update([read.qname])
                    try:
                        read_count[sample]['qc_fail'] += 1
                    except KeyError:
                        read_count[sample]['qc_fail'] = 1
            if read == None:
                #no reads were found for this contig, continue
                continue
            #process read_out
            if read_out == {} and 'RN' not in tag_dict:
                #random tag not implemented for this library. return in_files and do not process further
                return 0
            if region:
                reads = handle.fetch(cluster.id, region[0], region[1])
            else:
                reads = handle.fetch(cluster.id)
            for read in reads:
                tag_dict = dict(read.tags)
                tag = tag_dict['RN']
                sample = tag_dict['RG']
                max_AS = max(read_out[sample][tag].values())
                qname = [name for name, AS in read_out[sample][tag].items() if AS == max_AS][0]
                try:
                    read_count[sample]['count'] += 1
                except KeyError:
                    if sample not in read_count:
                        read_count[sample] = {'count':1}
                    else:
                        read_count[sample]['count'] =  1
                if read.qname in qc_fail:
                    read.is_qcfail = True
                if read.qname != qname:
                    #Read with duplicate tag with a lower quality score is marked as a PCR duplicate.
                    read.is_duplicate = True
                    try:
                        read_count[sample]['dup_count'] += 1
                    except KeyError:
                        read_count[sample]['dup_count'] = 1
                out_handle.write(read)
    print 'Sample:\tReads:\tDuplicates:\tDuplicate rate:'
    for key , subdict in sorted(read_count.items()):
        count = subdict['count']
        if 'dup_count' in subdict:
            dup_count = subdict['dup_count']
            dup_pct = dup_count / float(count)
            print '%s\t%s\t%s\t%.2f%%'%(key,count,dup_count,100*dup_pct)
        else:
            print '%s \t%s \t0\t0%%' % (key, count)
    print 'Sample:\tReads:\tQC-fail:\tQC-fail rate:'
    for key , subdict in sorted(read_count.items()):
        count = subdict['count']
        if 'qc_fail' in subdict:
            qc_fail = subdict['qc_fail']
            qc_fail_pct = qc_fail / float(count)
            print '%s\t%s\t%s\t%.2f%%'%(key,count,qc_fail,100*qc_fail_pct)
        else:
            print '%s \t%s \t0\t0%%' % (key, count)
    out_handle.close()
    cmd = ['samtools sort %s > %s' % (bam_out, bam_out.replace('.bam','.sorted.bam') )]
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True, executable='/bin/bash')
    stdout, stderr = p.communicate()
    stdout = stdout.replace('\r', '\n')
    stderr = stderr.replace('\r', '\n')
    print stdout
    print stderr
    os.rename(bam_out.replace('.bam','.sorted.bam'), bam_out)
    cmd = ['samtools index %s' % (bam_out)]
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True, executable='/bin/bash')
    stdout, stderr = p.communicate()
    stdout = stdout.replace('\r', '\n')
    stderr = stderr.replace('\r', '\n')
    print stdout
    print stderr
    return 0

def main():
    """main function loop"""
    #1 get command line arguments
    args = parse_args()
    if args.log == sys.stdout:
        log = sys.stdout
    else:
        log = open(args.log,'w')
    # 2 Remove PCR duplicates.
    if not args.input and args.input_dir:
        for strand in ['watson', 'crick']:
            print "started PCR duplicate removal for %s strand" % (strand)
            in_bam = os.path.join(args.input_dir,'%s.bam' % strand)
            out_bam = os.path.join(args.input_dir,'%s.dedup.bam' % strand)
            remove_PCR_duplicates(in_bam, out_bam, args.reference)
    else:
        remove_PCR_duplicates(args.input, args.output, args.reference)

if __name__ == '__main__':
    main()

