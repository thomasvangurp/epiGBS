__author__ = 'thomasvangurp'
"""filter bed file by pct coverage"""
import argparse
from Bio import SeqIO
import numpy

def parse_args():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(description='Process input files')
    parser.add_argument('-r', '--reference', type=str, help='reference genome input.')
    parser.add_argument('-i', '--info', type=str, help='Sample information files')
    parser.add_argument('-b', '--bed', help='bed file input')
    parser.add_argument('-o','--bedout', help='bed file output')
    parser.add_argument('-p','--percent',default="0.8",help='percent of samples called')
    parser.add_argument('-c','--mincover',default="3",help='minimum coverage for valid call')
    args = parser.parse_args()
    return args

def parse_sample_inf(args):
    """Parse sample informations and return as dictionary"""
    sample_groups = {}
    group_handle = open(args.info)
    header = group_handle.readline()
    for line in group_handle:
        split_line = line[:-1].split('\t')
        sample_groups[split_line[0]] = split_line[1]
    groups = sorted(set(sample_groups.values()))
    header =  4 * ['']
    for group in groups:
        header += [group, '']
    header += ['\n'] + (3 * ['']) + ['Average_methylation_%', 'SD_methylation_%'] * len(groups)
    return sample_groups, '\t'.join(header) + '\n'

def calculate_group(header, group_info, split_line):
    """Calculates average methylation and SD per group"""
    meth_ratios = {}
    for i in range(4, len(header),2):
        meth_count, total_count = split_line[i:i+2]
        try:
            meth_ratio = int(meth_count) / float(total_count)
        except ValueError:
            meth_ratio = None
        except ZeroDivisionError:
            meth_ratio = 0
        if meth_ratio != None:
            sample_name = header[i]
            group = group_info[sample_name.split('_')[0]]
            try:
                meth_ratios[group].append(meth_ratio)
            except KeyError:
                meth_ratios[group] = [meth_ratio]
    return meth_ratios

def process_chrom_stat(chrom_stat, bed_out_handle, split_line):
    """Process per contig statistics of average methylation rates"""
    if chrom_stat['CHH'] == {}:
        return 0
    for context, subdict in sorted(chrom_stat.items()):
        out_line = [split_line[0],'0' ,context,'contig_avg']
        for group, meth_values in sorted(subdict.items()):
            avg_meth = sum(meth_values) / float(len(meth_values))
            out_line[1] = str(len(meth_values))
            sd_meth = numpy.std(meth_values , ddof=1)
            out_line += ['%.4f' % avg_meth, '%.4f' % sd_meth]
        bed_out_handle.write('\t'.join(out_line) + '\n')


def summary_bed(args):
    """Filter bed file on percentage of sites being called and context"""
    bed_in_handle = open(args.bed,'r')
    bed_out_handle = open(args.bedout,'w')
    header = bed_in_handle.readline().split('\t')
    group_info, header_out = parse_sample_inf(args)
    bed_out_handle.write(header_out)
    min_call = round(len(header[4:])/2 * float(args.percent))
    min_depth = int(args.mincover)
    chrom_out = set()
    current_chrom = None
    chrom_stat = {'CG':{},'CHG':{},'CHH':{}}
    for line in bed_in_handle:
        split_line = line.split('\t')
        chrom = split_line[0]
        if current_chrom == None:
            current_chrom = split_line[0]
        elif chrom != current_chrom:
            current_chrom = split_line[0]
            process_chrom_stat(chrom_stat, bed_out_handle, split_line)
            chrom_stat = {'CG': {}, 'CHG': {}, 'CHH': {}}
        context = split_line[2]
        called = 0
        for i in range(4,len(header),2):
            try:
                count = split_line[i+1].rstrip('\n')
                total_calls = int(count)
                if total_calls >= min_depth:
                    called += 1
            except ValueError:
                pass
            except IndexError:
                called = 0
                break
        if called >= min_call:
            chrom_out.update(['chr%s'%split_line[0]])
            if not line.startswith('chr'):
                line = 'chr' + line
            summary = calculate_group(header, group_info, split_line)
            out_line = split_line[:4]
            for key , meth_values in sorted(summary.items()):
                avg_meth = sum(meth_values) / float(len(meth_values))
                try:
                    chrom_stat[context][key].append(avg_meth)
                except KeyError:
                    chrom_stat[context][key] = [avg_meth]
                sd_meth = numpy.std(meth_values , ddof=1)
                out_line += ['%.4f'%avg_meth, '%.4f'%sd_meth]
            bed_out_handle.write('\t'.join(out_line)+'\n')
    return chrom_out


def clean_fasta(fasta_input,fasta_output,seqs):
    """clean fasta file to get only reference sequences which are called"""
    fasta_input = SeqIO.parse(open(fasta_input,'r'),'fasta')
    out_handle = open(fasta_output,'w')
    for seq in fasta_input:
        name = 'chr%s'%seq.name
        if name in seqs:
            out = '>%s\n%s\n'%(name,seq.seq.tostring().upper())
            out_handle.write(out)
    out_handle.flush()
    out_handle.close()



def main():
    """main function loop"""
    args = parse_args()
    chrom_out = summary_bed(args)
    # clean_fasta(args.reference,args.refout,chrom_out)


    # return 0

if __name__ == '__main__':
    main()