#!/usr/bin/env python
"""Python script for converting PE fastq files for a epiGBS run
This script takes as input
--r1 left-hand fastq file /1
--r2 right-hand fastq file /2

For the right-hand file the reverse complement is generated.
Additionally, the correct barcode is appended at the start  to the read.
"""

import Levenshtein
import re
import sys
import os
import shutil
import subprocess
import shutil
from optparse import OptionParser
from Bio import SeqIO, Seq, Alphabet, Restriction
from itertools import product
from Bio.Seq import Seq
from Bio.Data.IUPACData import *
from StringIO import StringIO
from Bio.SeqRecord import SeqRecord
import operator
import tempfile
import gzip, bz2
import time


def parse_options():
    """Parses command line options"""
    parser = OptionParser()
    parser.add_option("--r1_in", metavar="reads1", action="store",
                      type="string", dest="reads1", help="left-hand fastq file")
    parser.add_option("--r2_in", metavar="reads2", action="store",
                      type="string", dest="reads2",
                      help="right-hand fastq file")
    parser.add_option("--mode", metavar="mode", action="store", type="string",
                      dest="mode", help="pe or single end mode", \
                      default="pe")
    parser.add_option("-b", "--barcodes", metavar="input", action="store",
                      type="string", dest="barcode", \
                      default="barcodes.csv",
                      help="input tab separated barcode file")
    parser.add_option("--output-dir", metavar="outputdir", action="store",
                      type="string", dest="outputdir", default="",
                      help="Specify output directory, only for galaxy")
    parser.add_option("-s", "--split", action="store_true",
                      default=0, dest="split",
                      help="Create multiple output files")
    parser.add_option("--addRG", action="store_true",
                      default=True, dest="addRG",
                      help="""Append append FASTA/Q comment to SAM output.
                      This option can be used to transfer read meta information
                      (e.g. barcode) to the SAM output. Note that the FASTA/Q
                      comment (the string after a space in the header line)
                      must conform the SAM spec (e.g. BC:Z:CGTAC).
                      Malformated comments lead to incorrect SAM output.""")
    parser.add_option("--match1", action="store", metavar="match1",
                      type="string", default="matching-R1",
                      dest="match1")
    parser.add_option("--match2", action="store", metavar="match2",
                      type="string", default="matching-R2",
                      dest="match2")
    parser.add_option("--stat", action="store", metavar="stat",
                      type="string", default="stats.txt",
                      dest="stat", help="statistics of read_nr per barcode")
    parser.add_option("--nomatch1", action="store", metavar="nomatch1",
                      type="string", default="/tmp/non-matching-R1.fastq",
                      dest="nomatch1", help="statistics of read_nr per barcode")
    parser.add_option("--nomatch2", action="store", metavar="nomatch2",
                      type="string", default="/tmp/non-matching-R2.fastq",
                      dest="nomatch2", help="statistics of read_nr per barcode")
    parser.add_option("-m", "--mismatches", action="store",
                      type="int", default=2, dest="mismatch",
                      help="Number of mismatches allowed")
    parser.add_option("-d", "--delete", action="store_true",
                      default=1, dest="delete",
                      help="Remove the barcode from the sequence, default is TRUE")
    parser.add_option("--control-nt", action="store_true",
                      default=1, dest="control_nucleotide",
                      help="implement barcode design with control nucleotide")
    return parser


def sort_barcodes(barcodes):
    """Sorts barcodes by length"""
    barcodes_out = []
    barcode_dic = {}
    for i in range(2):
        for barcode in barcodes:
            try:
                barcode_dic[len(barcode)] += [barcode]
            except KeyError:
                barcode_dic[len(barcode)] = [barcode]
        for key, list in sorted(barcode_dic.items()):
            for barcode in list:
                barcodes_out.append(barcode)
    return barcodes_out


def search_fast(sequence, barcodes, mismatch, position, enz_sites, max_bc_length):
    """Faster implementation of levenshtein"""
    sequence = sequence[1]
    max_total_len = max_bc_length + len(enz_sites[0])
    bc_len = min(sequence[:max_total_len].index(enz_sites[0]), sequence[:max_total_len].index(enz_sites[1]))
    if sequence[:bc_len] in barcodes:
        return sequence[:bc_len]


def get_strand(control_nt):
    """give strand given control nucleotide"""
    if control_nt == 'T':
        strand = 'Watson'
    elif control_nt == 'C':
        strand = 'Crick'
    else:
        strand = 'NA'
    return strand


def levenshtein(read, bc_set, enz_sites, mismatch, max_total_len):
    """Calculates the levenshtein distance between a sequence and a set of
    Barcodes. If the longest barcode with a perfect match
    is found this is returned and the script automatically quits"""
    sequence = read[1][:-1]
    # Process read 1 barcode
    short_sequence = read[1][
                     :max_total_len + 1]  # TODO: remove 1 and add conversion control nucleotides to barcodes.csv
    # for enz in enz_sites:
    #     if enz in short_sequence:
    #         short_sequence
    for (start, bc) in sorted(bc_set, key=lambda i: len(i[1]), reverse=True):
        wobble = short_sequence[:start]
        if short_sequence[start:start + len(bc)] == bc:
            control_nt = short_sequence[start + len(bc)]
            strand = get_strand(control_nt)
            return bc, wobble, start, strand
    matches = {}
    for barcode in sorted(bc_set, key=lambda i: len(i[1]), reverse=True):
        dist = []
        if barcode[1] in short_sequence[1:]:
            # this can happen if the wobble is shorter than it should be. minimum  wobble length > 1
            index = short_sequence.index(barcode[1])
            if short_sequence[index + len(barcode[1]):index + len(barcode[1]) + len(enz_sites[0])] not in enz_sites:
                break
            wobble = short_sequence[:index]
            control_nt = short_sequence[index + len(barcode[1])]
            strand = get_strand(control_nt)
            return barcode[1], wobble, index, strand
        for enz_site in enz_sites:
            try:
                part1 = short_sequence[barcode[0]:short_sequence.rindex(enz_site) + len(enz_site)]
            except ValueError:
                dist.append(100)
                continue
            part2 = barcode[1] + enz_site
            dist.append(Levenshtein.distance(part1, part2))
        try:
            # get the enz site with the min distance
            matches[min(dist)] += [barcode, enz_sites[dist.index(min(dist))]]
        except KeyError:
            matches[min(dist)] = [barcode, enz_sites[dist.index(min(dist))]]
    try:
        best_match, enz_site = matches[min(matches.keys())]
    except ValueError:
        # there are multiple matches do not return these conflicting values.
        return None, None, None, None
    if min(matches.keys()) <= mismatch:
        # the first nucleotide of the enz_site is the control nucleotide
        control_nt = matches[min(matches.keys())][1][0]
        strand = get_strand(control_nt)
        return best_match[1], short_sequence[:3], start, strand
        # return  best_match, enz_site
        # return  left_bc,wobble_left,left_start,control_left
    else:
        return None, None, None, None


def get_cutrem(enzyme):
    """Returns a list with the sequences of possible cut site remnants"""
    cutrem = []
    if enzyme.is_3overhang():
        site = enzyme.elucidate().split('_')[1].replace('^', '')
    if enzyme.is_5overhang():
        site = enzyme.elucidate().split('^')[1].replace('_', '')
    # compute all possible combination given one or more ambiguous bases
    pos_list = []
    #    site = site[:-1] turn on for epiGBS
    for base in site:
        pos_list += [ambiguous_dna_values[base]]
    cutrem = []
    for comb in product(*pos_list):
        cutrem.append(''.join(comb))
    return cutrem


class Barcode(object):
    """CLass to hold Barcode and enzyme informations """

    def __init__(self):
        self.Sample = None
        self.Flowcell = None
        self.Lane = None
        self.Barcode_R1 = None
        self.Barcode_R2 = None
        self.ENZ_R1 = None
        self.ENZ_R2 = None
        self.Wobble_R1 = 0
        self.Wobble_R2 = 0
        self.enz_remnant_R1 = ''
        self.enz_remnant_R2 = ''

    def get_seq(self):
        """Return sequence to search on left and right read"""
        R1_start = (self.Wobble_R1, self.Barcode_R1)
        R2_start = (self.Wobble_R2, self.Barcode_R2)
        return (R1_start, R2_start)


def get_enz_remnant(enz):
    """Get enzyme recognition site remnant sequence"""
    if enz.ovhg > 0:
        remnant = enz.site[enz.fst3:]
        return remnant
    else:
        remnant = enz.site[enz.fst5:]
        return remnant


def parse_bc(barcodes, fc, ln):
    """Parses barcode file and matches barcodes for specified flowcell and lane"""
    file_in = open(barcodes, 'r')
    bc_dict = {}
    header_index = {}
    for line in file_in.readlines():
        if line.startswith("#") or line.startswith("Flowcell"):
            for n, item in enumerate(line.rstrip('\n').split('\t')):
                header_index[n] = item
        else:
            bc_instance = Barcode()
            for n, item in enumerate(line.rstrip('\n').split('\t')):
                if header_index[n] in bc_instance.__dict__:
                    bc_instance.__setattr__(header_index[n], item)
            if bc_instance.ENZ_R1 != None:
                bc_instance.ENZ_R1 = get_enz(bc_instance.ENZ_R1)
            else:
                bc_instance.ENZ_R1 = get_enz('PstI')
            if bc_instance.ENZ_R2 != None:
                bc_instance.ENZ_R2 = get_enz(bc_instance.ENZ_R2)
            else:
                bc_instance.ENZ_R2 = get_enz('PstI')
            bc_instance.enz_remnant_R1 = get_enz_remnant(bc_instance.ENZ_R1)
            bc_instance.enz_remnant_R2 = get_enz_remnant(bc_instance.ENZ_R2)
            bc_instance.Wobble_R1 = int(bc_instance.Wobble_R1)
            bc_instance.Wobble_R2 = int(bc_instance.Wobble_R2)
            if Flowcell == fc and Lane == ln:
                bc_dict[bc_instance.get_seq()] = bc_instance
    return bc_dict


def read_type(left_read, right_read, left_enzsite, right_enzsite, left_bc, right_bc):
    """Determine if bisulfite read is watson or crick"""
    lr_enz_left = left_read[1][len(left_bc):len(left_bc) + 5]
    rr_enz_right = right_read[1][len(right_bc):len(right_bc) + 5]
    if left_enzsite == 'TACAA' and right_enzsite == 'TGCAG':
        return 'crick'
    elif right_enzsite == 'TACAA' and left_enzsite == 'TGCAG':
        return 'watson'
    elif right_enzsite == 'TGCAG' and left_enzsite == 'TGCAG':
        return 'gbs'
    else:
        # enzyme sites have not been establshed correctly, establish read
        # type based on closest matching enz site and CG count.
        watson_count = left_read[1].count('G') + right_read[1].count('C') + 0.001
        crick_count = left_read[1].count('C') + right_read[1].count('G') + 0.001
        left_distance = Levenshtein.distance(lr_enz_left, left_enzsite)
        right_distance = Levenshtein.distance(rr_enz_right, right_enzsite)
        if left_distance < right_distance:
            # left enz_site should be leading since it has fewer mismatches.
            if left_enzsite == 'TACAA' and crick_count / float(watson_count) > 2:
                return 'crick'
            else:
                return 'nodet'
        else:
            if left_enzsite == 'TGCAG' and watson_count / float(crick_count) > 2:
                return 'watson'
            else:
                return 'nodet'


def parse_seq_pe(opts, bc_dict, Flowcell, Lane):
    """Fastq/a-parser for PE-reads"""
    if opts.reads1.endswith('.gz'):
        seq1_handle = gzip.open(opts.reads1, "rb")
        seq2_handle = gzip.open(opts.reads2, "rb")
    elif opts.reads1.endswith('.bz2'):
        seq1_handle = bz2.open(opts.reads1, "rb")
        seq2_handle = bz2.open(opts.reads2, "rb")
    else:
        try:
            seq1_handle = open(opts.reads1, "r")
            seq2_handle = open(opts.reads2, "r")
        except IOError:
            seq1_handle = gzip.open(opts.reads1 + '.gz', "rb")
            seq2_handle = gzip.open(opts.reads2 + '.gz', "rb")
            opts.reads1 += '.gz'

    if not opts.split:
        seq1_name = '%(code)s_%(Flowcell)s_s_%(lane)s_fastq.txt' % \
                    ({'code': 'R1_%s' % opts.output.split('/')[-2], 'Flowcell': Flowcell, 'lane': Lane})
        seq2_name = '%(code)s_%(Flowcell)s_s_%(lane)s_fastq.txt' % \
                    ({'code': 'R2_%s' % opts.output.split('/')[-2], 'Flowcell': Flowcell, 'lane': Lane})
        if opts.reads1.endswith('.gz'):
            seq1_name += '.gz'
            seq2_name += '.gz'
            seq1_out = gzip.open(os.path.join(opts.output, seq1_name), 'a')
            seq2_out = gzip.open(os.path.join(opts.output, seq2_name), 'a')
        else:
            seq1_out = open(os.path.join(opts.output, seq1_name), 'a')
            seq2_out = open(os.path.join(opts.output, seq2_name), 'a')
    if opts.reads1.endswith('.gz'):
        nomatch1_out = gzip.open(opts.nomatch1, "w")
        nomatch2_out = gzip.open(opts.nomatch2, "w")
    else:
        nomatch1_out = open(opts.nomatch1, "w")
        nomatch2_out = open(opts.nomatch2, "w")
    seq = 0
    bc_set_left = set(k[0] for k in bc_dict.keys())
    bc_set_right = set(k[1] for k in bc_dict.keys())
    elements_1 = set(entry.enz_remnant_R1 for entry in bc_dict.values())
    elements_2 = set(entry.enz_remnant_R2 for entry in bc_dict.values())
    enz_sites_left = []
    enz_sites_right = []
    if opts.control_nucleotide:
        nt = 'C'
        for element in elements_1:
            if nt + element not in enz_sites_left:
                # implement search which includes control nucleotide
                enz_sites_left += [nt + element]
        for element in elements_2:
            if nt + element not in enz_sites_right:
                enz_sites_right += [nt + element]
    else:
        for element in elements_1[0]:
            if element[0] not in enz_sites_left:
                # implement search which includes control nucleotide
                enz_sites_left += [element]
        for element in elements_2[0]:
            if element[0] not in enz_sites_right:
                enz_sites_right += [element]
    max_bc_len_left = max(k[0][0] + len(k[0][1]) for k in bc_dict.keys()) + max(len(k) for k in enz_sites_left)
    max_bc_len_right = max(k[1][0] + len(k[1][1]) for k in bc_dict.keys()) + max(len(k) for k in enz_sites_right)
    left_read = [True]
    while left_read[0]:
        seq += 1
        left_read = []
        right_read = []
        for i in range(4):
            try:
                left_read += [seq1_handle.readline()]
                right_read += [seq2_handle.readline()]
            except StopIteration:
                break
        left_bc, wobble_left, left_start, control_left = levenshtein(left_read, bc_set_left, enz_sites_left,
                                                                     opts.mismatch, max_bc_len_left)
        right_bc, wobble_right, right_start, control_right = levenshtein(right_read, bc_set_right, enz_sites_right,
                                                                         opts.mismatch, max_bc_len_right)
        if left_bc and right_bc:
            # Put the correct sequence of the barcode
            try:
                bc_dict['%s_%s' % (left_bc, right_bc) + '_count'] += 1
            except KeyError:
                bc_dict['%s_%s' % (left_bc, right_bc) + '_count'] = 1
            if opts.addRG:
                # determine if read is watson or crick.
                try:
                    SM_id = bc_dict[((3, left_bc), (3, right_bc))].Sample
                except KeyError:
                    # This can only happen if the barcode is incorrectly read
                    try:
                        SM_id = bc_dict[((0, left_bc), (0, right_bc))].Sample
                    except KeyError:
                        continue
                # one control nucleotide should be converted the other not. If this succeeds than call read type (watson,crick)
                # based on left nucleotide. if this is
                if control_left != control_right:
                    strand = control_left
                else:
                    strand = control_left
                RG_id = '%s_%s_%s' % (Flowcell, Lane, SM_id)
                if wobble_left == '':
                    wobble_left = 'NNN'
                if wobble_right == '':
                    wobble_right = 'NNN'
                wobble = wobble_left + "_" + wobble_right
                left_read[0] = left_read[0].split(' ')[0].rstrip('\n') \
                               + '\tBC:Z:%s\tBC:Z:%s\tRG:Z:%s\tST:Z:%s\n' % (left_bc, right_bc, RG_id, strand)

                right_read[0] = right_read[0].split(' ')[0].rstrip('\n') \
                                + '\tBL:Z:%s\tBR:Z:%s\tRG:Z:%s\tST:Z:%s\n' % (left_bc, right_bc, RG_id, strand)
                if opts.control_nucleotide:
                    left_read[0] = left_read[0][:-1] + '\tRN:Z:%s\n' % wobble
                    right_read[0] = right_read[0][:-1] + '\tRN:Z:%s\n' % wobble
            else:
                id = left_read[0][:-1]
            if opts.delete:
                # +1 because of control nucleotide after barcode
                if opts.control_nucleotide:
                    control_NT = 'C'
                else:
                    control_NT = ''
                left_read[1] = left_read[1][left_start + len(left_bc + control_NT):]
                left_read[3] = left_read[3][left_start + len(left_bc + control_NT):]
                right_read[1] = right_read[1][right_start + len(right_bc + control_NT):]
                right_read[3] = right_read[3][right_start + len(right_bc + control_NT):]
            if not opts.split:
                seq1_out.write(''.join(left_read))
                seq2_out.write(''.join(right_read))
            else:
                # If splitting is activated, compression takes too long, disable!
                output_location_1 = os.path.join(opts.output,
                                                 "%s_%s_1.fastq" % (bc_dict[((3, left_bc), (3, right_bc))].Sample))
                output_location_2 = os.path.join(opts.output,
                                                 "%s_%s_2.fastq" % (bc_dict[((3, left_bc), (3, right_bc))].Sample))
                output_handle_1 = open(output_location_1, 'a')
                output_handle_2 = open(output_location_2, 'a')
                output_handle_1.write(''.join(left_read))
                output_handle_2.write(''.join(right_read))
        else:
            # Barcode sequence was not recognized
            nomatch1_out.write(''.join(left_read))
            nomatch2_out.write(''.join(right_read))
    seq1_out.close()
    seq2_out.close()
    nomatch1_out.close()
    nomatch2_out.close()
    return bc_dict


def parse_seq(opts, bc_sorted, bc_dict, Flowcell, Lane):
    """Fastq/a-parser for se-reads"""
    seq1_handle = open(opts.reads1, "rb")
    left_read = [1]
    enz_sites = get_cutrem(opts.enzyme)
    max_bc_length = len(bc_sorted[-1])
    if not opts.split:
        seq1_name = '%(code)s_%(Flowcell)s_s_%(lane)s_fastq.txt' % \
                    ({'code': 'R1samplecode123', 'Flowcell': Flowcell, 'lane': Lane})
        seq1_out = open(os.path.join(opts.output, seq1_name), 'a')
    nomatch_out = open(opts.nomatch1, "a")
    while left_read[0]:
        left_read = []
        for i in range(4):
            try:
                left_read += [seq1_handle.readline()]
            except Error:
                brake
        start_position = 0  # Position to start searching for barcode match.
        barcode = levenshtein(left_read, bc_sorted, opts.mismatch, \
                              start_position, enz_sites, max_bc_length)
        if barcode:
            # Put the correct sequence of the barcode
            try:
                bc_dict[barcode + '_count'] += 1
            except KeyError:
                bc_dict[barcode + '_count'] = 1
            if opts.rename:
                id = '@' + bc_dict[barcode] + '_%s' % bc_dict[barcode + '_count']
                left_read[0] = id + '\n'
            else:
                id = 0
            if opts.delete:
                left_read[1] = left_read[1][len(barcode):]
                left_read[3] = left_read[3][len(barcode):]
                right_read[1] = right_read[1][len(barcode):]
                right_read[3] = right_read[3][len(barcode):]
            if not opts.split:
                seq1_out.write(''.join(left_read))
            else:
                # galaxy required output: "%s_%s_%s_%s_%s" % ( 'primary', output1.id, name, 'visible', file_type,)
                output_location_1 = os.path.join(opts.output, "%s_%s.fastq" % (bc_dict[barcode], barcode))
                output_handle_1 = open(output_location_1, 'a')
                output_handle_1.write(''.join(left_read))
                output_handle_1.close()
        else:
            # Barcode sequence was not recognized
            nomatch_out.write(''.join(left_read))
    if not opts.split:
        seq1_out.close()
    nomatch_out.close()
    return bc_dict


def get_details_flow(opts):
    """Returns Flowcell and Lanes basef on fastq input"""
    if opts.reads1.endswith('.gz'):
        seq1_handle = gzip.open(opts.reads1, "rb")
    elif opts.reads1.endswith('.bz2'):
        seq1_handle = bz2.open(opts.reads1, "rb")
    else:
        try:
            seq1_handle = open(opts.reads1, "rb")
        except IOError:
            seq1_handle = gzip.open(opts.reads1 + '.gz', "rb")
    illumina_id = seq1_handle.readline().split(':')
    Flowcell, Lane = illumina_id[2:4]
    return Flowcell, Lane


def make_bc_record(record, barcode, id):
    """Returns a new SeqRecord with barcode plus sequence."""
    id += '/2'
    str_record = record
    if str_record[1][0] == "N":
        str_record[1] = opts.enzyme.ovhgseq[0] + str_record[1][1:]
        str_record[3] = "J" + str_record[3][1:]
    str_record[1] = barcode + str_record[1]
    str_record[3] = "J" * len(barcode) + str_record[3]
    if id:
        str_record[0] = id + '\n'
    return str_record


def parse_dir(opts):
    """Parse directory and return"""
    return 0


def get_enz(enz):
    """Get enzyme from biopython restriction library"""
    for enzyme in Restriction.AllEnzymes:
        if "%s" % (enzyme) == enz:
            return enzyme


def put_output(dir_in, opts, Flowcell, Lane):
    """Uses shutil to move the output into galaxy directory"""
    seq1_name = '%(code)s_%(Flowcell)s_s_%(lane)s_fastq.txt' % \
                ({'code': 'R1samplecode123', 'Flowcell': Flowcell, 'lane': Lane})
    seq2_name = '%(code)s_%(Flowcell)s_s_%(lane)s_fastq.txt' % \
                ({'code': 'R2samplecode123', 'Flowcell': Flowcell, 'lane': Lane})
    if not os.exists(os.path.join(dir_in, seq1_name)):
        seq1_name += ".gz"
        seq2_name += ".gz"
    shutil.move(os.path.join(dir_in, seq1_name), opts.match1)
    shutil.move(os.path.join(dir_in, seq2_name), opts.match2)
    return 0


def write_stats(bc_dict, opts):
    """Write stats to output file"""
    barcode_in = open(opts.barcode, 'r')
    stat_out = open(opts.stat, "w")
    # Write the first line = header to output barcode field.
    for line in barcode_in.readlines():
        fc = line.split('\t')[0]
        ln = line.split('\t')[1]
        if not (fc == Flowcell and ln == Lane):
            continue
        name = line.split('\t')[4]
        bc_left, bc_right = line.split('\t')[2:4]
        try:
            bc_count = bc_dict['%s_%s' % (bc_left, bc_right) + '_count']
            stat_out.write("%s\t" * 3 % (name, '%s_%s' % (bc_left, bc_right), bc_count) + '\n')
        except KeyError:
            stat_out.write("%s\t" * 3 % (name, '%s_%s' % (bc_left, bc_right), 0) + '\n')
    stat_out.close()


if __name__ == "__main__":
    parser = parse_options()
    opts, args = parser.parse_args()

    # Make sure we identify the  flowcell and lane records in the barcodefile that correspond to our fastq file
    Flowcell, Lane = get_details_flow(opts)
    bc_dict = parse_bc(opts.barcode, Flowcell, Lane)
    opts.output = tempfile.mkdtemp(prefix='seq', dir=opts.outputdir)
    if os.path.exists(opts.output):
        # TODO: check content to see if deletion is warranted
        pass
        # shutil.rmtree(opts.output)
        # os.mkdir(opts.output)
    else:
        os.mkdir(opts.output)
    # if opts.outputdir:
    #     try:
    #        file_out = open(opts.outputdir, 'w')
    #        file_out.write('%s'%opts.output)
    #        file_out.close()
    #     except OSError:
    #         #TODO: determine error type
    #         pass
    if opts.mode == 'pe':
        parse_seq_pe(opts, bc_dict, Flowcell, Lane)
    else:
        parse_seq(opts, bc_sorted, bc_dict, Flowcell, Lane)
    write_stats(bc_dict, opts)
    if opts.match1 != 'matching-R1':
        # match1 is the default variable name.
        put_output(opts.output, opts, Flowcell, Lane)
    print "Done."


