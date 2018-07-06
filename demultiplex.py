#!/usr/bin/env python
# use http://rst.ninjs.org for viewing these figures
# see http://docutils.sourceforge.net/docs/user/rst/quickref.html for a rst reference
"""
Python module for barcode deconvolution of paired-end fastq files using a barcode file

purpose
-------
The *purpose* of this module is to assign sequences to samples and remove barcodes.
For this, the following steps are taken:

1. Identify the barcode and enzyme recognition site in the forward and reverse read allowing for some mismatches
2. Remove the nucleotides and quality score letters of the barcode and adapter-derived oligonucleotides
3. Match the identified barcode from forward and reverse read to a sample and Flowcell using a barcodes.tsv file
4. Add a sample and read group identifier based on the sample and the flowcell number and lane

**Mandatory** parameter to be defined for this module are

Input files
===========
*gzip compressed files are allowed and even encouraged!*

--r1_in  left-hand or **forward** fastq file or **R1** reads e.g. *R1.fastq.gz*
--r2_in  right-hand or **reverse** fastq file or **R2** reads e.g. *R2.fastq.gz*

Configuration options
=====================

--addRG  Append sample and read group tags in SAM format to the read name. BWA-mem and other
aligners can include this information in the SAM/BAM output.


Output files
============
*will be compressed if the .gz extension is present*

--match1  left-hand or **forward** fastq file or **R1** output e.g. *R1.out.fastq.gz*
--match2  right-hand or **reverse** fastq file or **R2** output e.g. *R2.out.fastq.gz*
--nomatch1  non-matching left-hand or **forward** fastq file or **R1** output e.g. *R1.out.fastq.gz*
--nomatch2  non-matching right-hand or **reverse** fastq file or **R2** output e.g. *R2.out.fastq.gz*

V3.0 accepts different enzyme combinations and barcodes simultaneously

:Authors:
     Thomas P. van Gurp
:Version: 3.0 of 2017/10/2
"""

__version__ = 3.0

import Levenshtein
import os
import shutil
from optparse import OptionParser
from Bio import Restriction
from itertools import product
from Bio.Data.IUPACData import *
import tempfile
import gzip
import bz2
# from mem_top import mem_top



def parse_options():
    """Parses command line options"""
    parser = OptionParser()
    parser.add_option("--r1_in", metavar="reads1", action="store",
                      type="string", dest="reads1", help="left-hand fastq file")
    parser.add_option("--r2_in", metavar="reads2", action="store",
                      type="string", dest="reads2",
                      help="right-hand fastq file")
    parser.add_option("-b", "--barcodes", metavar="input", action="store",
                      type="string", dest="barcode", default="barcodes.tsv",
                      help="input tab separated barcode file")
    parser.add_option("--output-dir", metavar="outputdir", action="store",
                      type="string", dest="outputdir", default="",
                      help="Specify output directory, only for galaxy")
    parser.add_option("-s", "--split", action="store_true",
                      default=False, dest="split",
                      help="Create multiple output files *NOT* recommended")
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
                      type="string", default=None,
                      dest="stat", help="statistics of read_nr per barcode")
    parser.add_option("--nonconversion", action="store", metavar="nonconversion",
                      type="string", default=None,
                      dest="nonconversion", help="Non conversion statistics")
    parser.add_option("--nomatch1", action="store", metavar="nomatch1",
                      type="string", default="/tmp/non-matching-R1.fastq",
                      dest="nomatch1", help="statistics of read_nr per barcode")
    parser.add_option("--nomatch2", action="store", metavar="nomatch2",
                      type="string", default="/tmp/non-matching-R2.fastq",
                      dest="nomatch2", help="statistics of read_nr per barcode")
    parser.add_option("-m", "--mismatches", action="store",
                      type="int", default=2, dest="mismatch",
                      help="Number of mismatches allowed")
    parser.add_option("--mode", action="store",
                      default="pe", dest="mode",
                      help="Paired-end is default and only valid choice")
    parser.add_option("-d", "--delete", action="store_true",
                      default=1, dest="delete",
                      help="Remove the barcode from the sequence, default is TRUE")
    parser.add_option("--control-nt", action="store_true",
                      default=0, dest="control_nucleotide",
                      help="implement barcode design with control nucleotide")
    opts, args = parser.parse_args()
    if opts.stat == None:
        opts.stat = os.path.join(opts.outputdir, 'demultiplex_stats.tsv')
    if opts.nonconversion == None:
        opts.nonconversion = os.path.join(opts.outputdir, 'nonconversion.tsv')
    return opts, args


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
    if control_nt in 'CT':
        return control_nt
    else:
        return None


def get_variants(barcode):
    """return list of possible barcodes given ambiguous nucleotides"""
    pos_list = []
    for base in barcode:
        pos_list += [ambiguous_dna_values[base]]
    cutrem = []
    for comb in product(*pos_list):
        cutrem.append(''.join(comb))
    return cutrem


def levenshtein(read, bc_set, mismatch, max_total_len, control_IUPAC='Y'):
    """Calculates the levenshtein distance between a sequence and a set of
    Barcodes. If the longest barcode with a perfect match
    is found this is returned and the script automatically quits"""
    try:
        var = read * 2
    except:
        var = read * int(read)

    sequence = read[1][:-1]
    # Process read 1 barcode
    short_sequence = read[1][
                     :max_total_len]
    for (start, bc) in sorted(bc_set, key=lambda i: len(i[1]), reverse=True):
        wobble = short_sequence[:start]
        break
    matches = {}
    for barcode in sorted(bc_set, key=lambda i: len(i[1]), reverse=True):
        dist = []
        bc_variants = get_variants(barcode[1])

        # we need a function to make multiple barcodes with enzyme sites for ambiguous nucleotides
        for bc_variant in bc_variants:
            if bc_variant in short_sequence[min(1, start):]:
                # this can happen if the wobble is shorter than it should be. minimum  wobble length > 1
                wobble = short_sequence[:short_sequence.index(bc_variant)]
                control_nt_index = barcode[1].index(control_IUPAC)
                control_nt = short_sequence[len(wobble) + control_nt_index]
                # wobble = short_sequence[:index]
                return barcode[1], wobble, len(wobble), control_nt, 0
            try:
                part1 = short_sequence[barcode[0]:len(barcode[1]) + barcode[0]]
            except ValueError:
                dist.append(100)
                continue
            part2 = bc_variant
            dist.append(Levenshtein.distance(part1.decode(), part2.decode()))
        try:
            # get the enz site with the min distance
            matches[min(dist)] += [barcode]
        except KeyError:
            matches[min(dist)] = [barcode]
    if len(matches[min(matches.keys())]) == 1:
        start, barcode = matches[min(matches.keys())][0]
    else:
        # there are multiple matches do not return these conflicting values.
        return None, None, None, None, None
    if min(matches.keys()) <= mismatch:
        # the left_most Y gives the location of the control-nucleotide
        control_nt_index = matches[min(matches.keys())][0][1].index(control_IUPAC)
        control_nt = short_sequence[len(wobble) + control_nt_index]
        return barcode, short_sequence[:start], start, control_nt, min(matches.keys())
    else:
        return None, None, None, None, None


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
        # design of Read_1 is NNN|BARCODE|CONTROL-NT|ENZ-REMNANT
        # CONTROL-NT for R1  is either C or T, put Y as control nucleotide
        R1_start = (self.Wobble_R1, self.Barcode_R1 + 'Y' + self.enz_remnant_R1)
        # CONTROL-NT for R2  is either G or A, put R as control nucleotide
        R2_start = (self.Wobble_R2, self.Barcode_R2 + 'Y' + self.enz_remnant_R2)
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
    for line_number, line in enumerate(file_in.readlines()):
        if line_number == 0:
            for n, item in enumerate(line.rstrip('\n').split('\t')):
                header_index[n] = item
            try:
                assert 'Sample' in header_index.values()
            except AssertionError:
                raise KeyError('"Sample" not in header, please revise the header ' +
                               '(first lines) of the barcode file %s' % barcodes)
        else:
            bc_instance = Barcode()
            for n, item in enumerate(line.rstrip('\n').split('\t')):
                if header_index[n] in bc_instance.__dict__:
                    bc_instance.__setattr__(header_index[n], item)
            if bc_instance.ENZ_R1 != None:
                bc_instance.ENZ_R1 = get_enz(bc_instance.ENZ_R1)
            else:
                bc_instance.ENZ_R1 = get_enz('PstI')
            bc_instance.enz_remnant_R1 = get_enz_remnant(bc_instance.ENZ_R1)
            bc_instance.Wobble_R1 = int(bc_instance.Wobble_R1)

            if bc_instance.ENZ_R2 != None:
                bc_instance.ENZ_R2 = get_enz(bc_instance.ENZ_R2)
            else:
                bc_instance.ENZ_R2 = get_enz('PstI')
            bc_instance.enz_remnant_R2 = get_enz_remnant(bc_instance.ENZ_R2)
            bc_instance.Wobble_R2 = int(bc_instance.Wobble_R2)
            if bc_instance.Flowcell == fc and bc_instance.Lane == ln:
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

def write_non_conversion_stats(args, conversion):
    """write non conversion statistics on a per sample basis."""
    with open(args.nonconversion,'w') as stat_out:
        header = ['Sample','Watson-count','Crick-count','Non-conversion-count','non-conversion-rate','other-count']
        stat_out.write('\t'.join(header) + '\n')
        for sample, count_dict in sorted(conversion.items()):
            try:
                watson_count = count_dict['C_T']
                crick_count = count_dict['T_C']
                non_conv_count = count_dict['C_C']
                non_conv_rate = non_conv_count / float(non_conv_count +
                                                       watson_count +
                                                       crick_count)
            except KeyError:
                continue
            line_out = [sample, str(watson_count), str(crick_count), str(non_conv_count),
                        '%.5f' % non_conv_rate, str(sum(count_dict.values()) - watson_count - crick_count - non_conv_count)]
            stat_out.write('\t'.join(line_out) + '\n')
    print(stat_out.name)

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
            seq1_out = gzip.open(os.path.join(opts.output, seq1_name), 'wb')
            seq2_out = gzip.open(os.path.join(opts.output, seq2_name), 'wb')
        else:
            seq1_out = open(os.path.join(opts.output, seq1_name), 'w')
            seq2_out = open(os.path.join(opts.output, seq2_name), 'w')
    if opts.reads1.endswith('.gz'):
        nomatch1_out = gzip.open(opts.nomatch1, "wb")
        nomatch2_out = gzip.open(opts.nomatch2, "wb")
    else:
        nomatch1_out = open(opts.nomatch1, "w")
        nomatch2_out = open(opts.nomatch2, "w")
    seq = 0
    bc_set_left = set((k[0]) for k in bc_dict.keys())
    bc_set_right = set(k[1] for k in bc_dict.keys())
    # elements_1 = set(entry.enz_remnant_R1 for entry in bc_dict.values())
    # elements_2 = set(entry.enz_remnant_R2 for entry in bc_dict.values())
    # enz_sites_left = []
    # enz_sites_right = []
    # if opts.control_nucleotide:
    #     nt = 'C'
    #     for element in elements_1:
    #         if nt + element not in enz_sites_left:
    #             # implement search which includes control nucleotide
    #             enz_sites_left += [nt + element]
    #     for element in elements_2:
    #         if nt + element not in enz_sites_right:
    #             enz_sites_right += [nt + element]
    # else:
    #     for element in elements_1[0]:
    #         if element[0] not in enz_sites_left:
    #             # implement search which includes control nucleotide
    #             enz_sites_left += [element]
    #     for element in elements_2[0]:
    #         if element[0] not in enz_sites_right:
    #             enz_sites_right += [element]
    max_bc_len_left = max([k[0] + len(k[1]) for k in bc_set_left])
    max_bc_len_right = max([k[0] + len(k[1]) for k in bc_set_right])
    left_read = [True]
    conversion = {}
    while left_read[0]:
        seq += 1
        left_read = []
        right_read = []
        for i in range(4):
            try:
                left_read += [seq1_handle.readline().decode()]
                right_read += [seq2_handle.readline().decode()]
            except StopIteration:
                break
        left_bc, wobble_left, left_start, control_left, mismatch_left = levenshtein(left_read, bc_set_left,
                                                                                    opts.mismatch, max_bc_len_left)
        right_bc, wobble_right, right_start, control_right, mismatch_right = levenshtein(right_read, bc_set_right,
                                                                                         opts.mismatch,
                                                                                         max_bc_len_right)
        if left_bc and right_bc:
            # Put the correct sequence of the barcode
            try:
                SM_id = bc_dict[((3, left_bc), (3, right_bc))].Sample
            except KeyError:
                continue
            try:
                bc_dict[SM_id + '_count'] += 1
            except KeyError:
                bc_dict[SM_id + '_count'] = 1
            if opts.addRG:
                # determine if read is watson or crick.
                try:
                    SM_id = bc_dict[((3, left_bc), (3, right_bc))].Sample
                except KeyError:
                    # This can only happen if the barcode is incorrectly read
                    try:
                        SM_id = bc_dict[((0, left_bc), (3, right_bc))].Sample
                    except KeyError:
                        continue
                RG_id = '%s_%s_%s' % (Flowcell, Lane, SM_id)
                try:
                    conversion[SM_id][control_left + '_' + control_right] += 1
                except KeyError:
                    if SM_id not in conversion:
                        conversion[SM_id] = {}
                    conversion[SM_id][control_left + '_' + control_right] = 1
                if control_left != control_right:
                    #this is the default situation. A read is arbitrarily assigned to "Watson" or "Crick"
                    #depending on the conversion occuring on the forward or reverse read
                    if control_left == 'C':
                        #the control nucleotide on the forward read was not converted, we call this a Watson read
                        if control_right == 'T':
                            strand = 'Crick'
                        else:
                            nomatch1_out.write(''.join(left_read).encode('utf-8'))
                            nomatch2_out.write(''.join(right_read).encode('utf-8'))
                            continue
                    elif control_left == 'T':
                        # the control nucleotide on the forward read was converted, we call this a Crick read
                        if control_right == 'C':
                            strand = 'Watson'
                        else:
                            nomatch1_out.write(''.join(left_read).encode('utf-8'))
                            nomatch2_out.write(''.join(right_read).encode('utf-8'))
                            continue
                    else:
                        nomatch1_out.write(''.join(left_read).encode('utf-8'))
                        nomatch2_out.write(''.join(right_read).encode('utf-8'))
                        continue
                else:
                    nomatch1_out.write(''.join(left_read).encode('utf-8'))
                    nomatch2_out.write(''.join(right_read).encode('utf-8'))
                    continue
                if wobble_left == '':
                    wobble_left = 'NNN'
                if wobble_right == '':
                    wobble_right = 'NNN'
                wobble = wobble_left + "_" + wobble_right
                left_read[0] = left_read[0].split(' ')[0].rstrip('\n') \
                               + '\tBL:Z:%s\tBR:Z:%s\tRG:Z:%s\tML:i:%s\tMR:i:%s\tST:Z:%s\n' \
                                 % (left_bc, right_bc, RG_id, mismatch_left, mismatch_right, strand)

                right_read[0] = right_read[0].split(' ')[0].rstrip('\n') \
                                + '\tBL:Z:%s\tBR:Z:%s\tRG:Z:%s\tML:i:%s\tMR:i:%s\tST:Z:%s\n' \
                                  % (left_bc, right_bc, RG_id, mismatch_left, mismatch_right, strand)
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
                left_bc_only = bc_dict[((3, left_bc), (3, right_bc))].Barcode_R1
                right_bc_only = bc_dict[((3, left_bc), (3, right_bc))].Barcode_R2
                left_read[1] = left_read[1][left_start + len(left_bc_only + control_NT):]
                left_read[3] = left_read[3][left_start + len(left_bc_only + control_NT):]
                right_read[1] = right_read[1][right_start + len(right_bc_only + control_NT):]
                right_read[3] = right_read[3][right_start + len(right_bc_only + control_NT):]
            if not opts.split:
                seq1_out.write(''.join(left_read).encode('utf-8'))
                seq2_out.write(''.join(right_read).encode('utf-8'))
            else:
                # If splitting is activated, compression takes too long, disable!
                output_location_1 = os.path.join(opts.output,
                                                 "%s_%s_1.fastq" % (bc_dict[((3, left_bc), (3, right_bc))].Sample))
                output_location_2 = os.path.join(opts.output,
                                                 "%s_%s_2.fastq" % (bc_dict[((3, left_bc), (3, right_bc))].Sample))
                output_handle_1 = open(output_location_1, 'a')
                output_handle_2 = open(output_location_2, 'a')
                output_handle_1.write(''.join(left_read).encode('utf-8'))
                output_handle_2.write(''.join(right_read).encode('utf-8'))
        else:
            # Barcode sequence was not recognized
            nomatch1_out.write(''.join(left_read).encode('utf-8'))
            nomatch2_out.write(''.join(right_read).encode('utf-8'))
    seq1_out.close()
    seq2_out.close()
    nomatch1_out.close()
    nomatch2_out.close()
    write_non_conversion_stats(opts, conversion)
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
    illumina_id = seq1_handle.readline().decode().split(':')
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


# #module unittest
# def test_funtion():
#

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
    indexes = {}
    for line in barcode_in.readlines():
        fc = line.split('\t')[0]
        ln = line.split('\t')[1]
        if not (fc == Flowcell and ln == Lane):
            for k, n in enumerate(line.rstrip('\n').split('\t')):
                indexes[n] = k
            continue
        name = line.split('\t')[indexes['Sample']]
        bc_left, bc_right = line.split('\t')[2:4]
        try:
            bc_count = bc_dict[name + '_count']
            stat_out.write("%s\t" * 3 % (name, '%s_%s' % (bc_left, bc_right), str(bc_count)) + '\n')
        except KeyError:
            stat_out.write("%s\t" * 3 % (name, '%s_%s' % (bc_left, bc_right), '0') + '\n')
    stat_out.close()

opts, args = parse_options()
Flowcell, Lane = get_details_flow(opts)

bc_dict = parse_bc(opts.barcode, Flowcell, Lane)
if not os.path.exists(opts.outputdir):
    os.mkdir(opts.outputdir)
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
    parse_seq(opts, bc_dict, Flowcell, Lane)
write_stats(bc_dict, opts)
if opts.match1 != 'matching-R1':
    # match1 is the default variable name.
    put_output(opts.output, opts, Flowcell, Lane)
print("Done.")

