#!/usr/bin/env python
__author__ = 'Bjorn Wouters'
__email__ = "bjorn-wouters@hotmail.com"

"""
Description: Preparing analysis files for the RnBeads R package.
Version: 1.0.0.
"""

import argparse
import os

def main():
    args = ParseArgs()
    input, output_dict, samples = ParseFiles(args.input, args.output)
    output = args.output
    IgvToRnBeads(input, output_dict, samples, output)


def ParseArgs():
    """
    Parses the arguments from the argparse module to the args variable.
    """
    parser = argparse.ArgumentParser(description='Process input files')
    parser.add_argument('-i', '--input', type=str, help='Input needs to be in .bed  or .igv format.')
    parser.add_argument('-o', '--output', type=str, help='Output file path')
    return parser.parse_args()


def ParseFiles(input, output):
    """
    Makes a bed file out of each sample that is in the header of the input .bed file.
    """
    input_file = open(input, 'r')
    output_file_dict = dict()
    header = input_file.next()
    samples = header.split()[4:]

    # Iterates through the sample in the header.
    for i in range(0, len(samples), 2):
        sample = samples[i].strip('_methylated')
        output_file_dict.update({sample: open(output+'/'+sample+'.bed', 'w')})

    return [input_file, output_file_dict, samples]


def IgvToRnBeads(input_file, output_dict, samples, output, given_samples=None, min_reads=5):
    """
    Fills each sample.bed file with the methylation info that is obtained from the input .bed file.
    """
    # Creates a file where every chromosome is put in if there is methylation data from with a coverage > 5.
    chr_file = open(os.path.join(output, "chromosomes.txt"), "w")
    # Valid chromosome = chromosome with atleast 1 sample that has a coverage higher than 5 on it.
    valid_chromosomes = set()
    for line in input_file:
        chr, pos, context = line.split()[:3]
        if context == 'CG':
            try:
                # Need info of next line because the context also, needs to be "CG".
                next_line = input_file.next()
            except StopIteration:
                pass
            if next_line.split()[2] != 'CG':  # Needs to be "CG" as well.
                continue
            else:
                # None means actually that there is no read data but data with 0 reads
                # but data with zero (0) won't be written to the .bed file's anyway.
                normalised_line = line.replace('None', '0')
                normalised_next_line = next_line.replace('None', '0')
                values_C = normalised_line.split()[4:]
                values_G = normalised_next_line.split()[4:]
                for i in range(0, len(samples), 2):
                    total = int(values_C[i+1])+int(values_G[i+1])
                    # Needs to be more than the offset.
                    if total <= int(min_reads):
                        continue
                    # Appends the chromosome list if a new valid chromosome is passed the total offset.
                    if chr not in valid_chromosomes:
                        valid_chromosomes.add(chr)
                        chr_file.write("".join([chr.strip("chr"), "\n"]))
                    sample = samples[i].strip('_methylated')
                    sample_file = output_dict[sample]
                    out_str = chr+"\t%s\t%s\t"%(int(pos)-1, int(pos)+1)
                    sample_file.write(out_str)
                    #Methylation in CG sites is determined for both watson and crick strand combined!
                    methylated = int(values_C[i])+int(values_G[i])
                    if total == 0:
                        ratio = 0
                    else:
                        ratio = float(methylated)/total*1000
                    strand = '+'
                    # Sample file needs to be in EPP format: see RnBeads vignette for a more detailed description.
                    sample_file.write("'"+str(methylated)+'/'+str(total)+"'"+'\t'+str("%.0f" % ratio)+'\t'+strand+'\n')

    # Closes all files in the dictionary.
    for file in output_dict.values():
        file.close()

    invalid_samples = list()
    valid = True
    # Makes a dictionary of the samples and the number of lines the files`s files. The number of line represent
    # The number of CG sites each sample has with a coverage >= 5.
    content = {sample: sum(1 for line in open(output_dict[sample].name)) for sample in output_dict.keys()}
    #mean = sum(content.values())/len(content.values())
    offset = max(content.values())*0.05
    # If given samples (list). Than only the samples given in that list will be analysed on sites differences.
    if given_samples:
        for sample in given_samples:
            if content[sample] < offset:
                valid = False
                invalid_samples.append(sample)
    else:
        for sample in content.keys():
            if content[sample] < offset:
                valid = False
                invalid_samples.append(sample)

    input_file.close()
    chr_file.close()

    if not valid:
        return invalid_samples
    else:
        return None


def chrom_sizes(fasta, tmp, ac):
    """Prepare chromsizes file for RNbeads"""
    from Bio import SeqIO
    output_file = open(os.path.join(tmp, ac+".chrom.sizes"), 'w')
    for record in SeqIO.parse(fasta, "fasta"):
        record_size = len(record.seq._data)
        output_file.write("chr"+record.name+"\t"+str(record_size)+"\n")
    output_file.close()
    return output_file


class SampleFileError(Exception):
    def __init__(self, value):
        self.value = value

    def __str__(self):
        return repr(self.value)

#Calls the main function.
if __name__ == '__main__':
    main()