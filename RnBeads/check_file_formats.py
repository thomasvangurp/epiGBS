
import argparse
from Bio import SeqIO
import re


def main():
    # Parses the given arguments with the argparse module.
    args = parse_args()
    # If the fasta file is given, check on file syntax:
    if args.fasta:
        check_fasta(args.fasta)
    if args.bed:
        check_bed(args.bed)
    if args.sample:
        check_sample_file(args.sample)


def check_fasta(fasta, valid_nucl="[ATGCN]"):
    """
    Checks the given .fasta file on specific syntax:
    if sequence name is numeric: it has to begin with "chr".
    Sequence only consist out of: A, T, G, C and N nucleotides.
    """
    for record in SeqIO.parse(fasta, "fasta"):
        seq = record.seq._data
        if not re.match(valid_nucl, seq):
            raise FormatError("Non-valid nucleotides in the given fasta file.")
        name = record.name
        if not name.startswith("chr"):
            if not name.isalpha():
                raise FormatError("Non-valid chromosome for the fasta name (only chr+[0-9] or [a-z][A-Z] possible).")
    return True


def check_bed(bed):
    bed_file = open(bed, "r")
    header = bed_file.next().split()
    header_format = ["chr", "pos", "context", "samples_called"]
    for i, item in enumerate(header):
        if i > 3:
            if not (item.endswith("_total") or item.endswith("_methylated")):
                raise FormatError("Non-valid sample header (has to end with _total or _methylated).")
        else:
            if item != header_format[i]:
                raise FormatError("Non-valid header names: format = chr\\tpos\\tcontext\\tsamples_called")
    bed_file.close()
    return True


def check_sample_file(sample_input):
    sample_file = open(sample_input, "r")
    header = sample_file.next().rstrip('\n').split(",")
    header_format = ["sampleID", "filename_bed"]
    for i, item in enumerate(header[0:2]):
        if item != header_format[i]:
            raise FormatError("Wrong header format. (Needs to be sampleID\\tfilename_bed)")
    for item in header[2:]:
        split_item = item.split("_")
        for word in split_item:
            if not word.strip().isalpha():
                raise FormatError("Wrong header name (may only consist out of [a-z][A-Z] and underscore ( _ )")
    for sample in sample_file:
        split_line = sample.rstrip('\n').split(",")
        if len(split_line) != len(header):
            raise FormatError("Length not the sample as the header, you're missing a value with a sample")
        if split_line[1] != split_line[0]+".bed":
            raise FormatError("Bed column incorrect: second column = first column name + .bed")
    return True


def parse_args():
    parser = argparse.ArgumentParser(description="Check RnBeads input file formats")
    parser.add_argument('-f', '--fasta', help='Input fasta file', default=None)
    parser.add_argument('-b', '--bed', help='Input bed file', default=None)
    parser.add_argument('-s', '--sample', help='Input sample (csv) file', default=None)
    args = parser.parse_args()
    return args


class FormatError(Exception):
    def __init__(self, value):
        self.value = value

    def __str__(self):
        return repr(self.value)


# Calls the main function if the script is executed.
if __name__ == "__main__":
    main()