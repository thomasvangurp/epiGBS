#!/usr/bin/env python
"""add DMPs to jbrowse as gff3 and bigbed files"""
import argparse
import os
import gzip
import cPickle as pickle
import xml.etree.cElementTree as ET
from Bio import SeqIO
from nested_dict import nested_dict


def parse_args():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(description='Concatenate reference for display in jbrowse')
    parser.add_argument('--rnbeads_dir', help='Output directory from Rnbeads')
    parser.add_argument('--code', help='code to look for')
    parser.add_argument('-c', '--context', help='Nucleotide context for DMPs')
    parser.add_argument('-b', '--bed', help='bed file with contig to concatenated mapping')
    parser.add_argument('-r', '--ref', help='reference sequence for contigs')
    parser.add_argument('-o', '--outputdir', help='output directory for gff3 and bigbed files')
    parser.add_argument('--treshold', help='minimum p value for DMP',default='0.05')
    args = parser.parse_args()
    for context in ['cg','chg','chh']:
        path = os.path.join(args.rnbeads_dir,args.code + context,'analysis')
        if os.path.exists(path):
            os.chdir(path)
            dirs = filter(os.path.isdir, os.listdir(path))
            dirs = [os.path.join(path, d) for d in dirs]  # add path to each file
            dirs.sort(key=lambda x: os.path.getmtime(x))
            #take dir with latest creation date
            dir = dirs[-1]
            vars(args)['rnbeads_dir_%s' % context] = dir
    if not os.path.exists(args.outputdir):
        os.mkdir(args.outputdir)
    return args

def get_mapping(args):
    """parse bed file for mapping contigs to concatenated ref"""
    mapping_dict = {}
    with open(args.bed,'r') as handle:
        for line in handle:
            target,target_start,target_end,name = line.rstrip('\n').split('\t')
            mapping_dict[name] = (target,target_start,target_end)
    return mapping_dict


def get_DMP_list(dir, args):
    """parse HTML of Rnbeads to get DMP list"""
    handle = open(os.path.join(dir,'differential_methylation.html'),'r')
    DMP_list = {}
    for line in handle:
        if '<li>' in line and '.csv' in line:
            file_location = os.path.join(dir, line[line.index('"')+1:line.rindex('"')])
            comparison = line[line.index('csv')+5:]
            comparison = comparison[:comparison.index('</a>')]
            DMP_list[comparison] = file_location
    return DMP_list

def check_context(cluster,position,context):
    """Check if context matches specified position"""
    nt = cluster[position]
    assert nt in 'NACGT'
    if nt == 'C':
        try:
            up_2 = cluster[position+1:position+3]
            if up_2[0] == 'G':
                actual_context = 'CG'
            elif up_2[1] == 'G':
                actual_context = 'CHG'
            else:
                actual_context = 'CHH'
        except IndexError:
            # TODO: replace by expected nucleotides from enz recognition site
            return 1
    elif nt == 'G':
        try:
            down_2 = cluster[position-2:position]
            if down_2[1] == 'C':
                actual_context = 'CG'
            elif down_2[0] == 'C':
                actual_context = 'CHG'
            else:
                actual_context = 'CHH'
        except IndexError:
            return 1
    else:
        print 'Not a C nucleotide!'
        return 1
    if actual_context != context.upper():
        # print "Context mismatch!"
        return 1
    else:
        return 0


def context_mapping(contig, context):
    """reverse mapping of 'fake' RNbeads positions to 'real' positions"""
    mapping_dict = {}
    offset = 0
    if context.upper() == 'CG':
        for pos,nt in enumerate(str(contig.seq.upper())):
            if nt == 'C' or nt == 'G':
                return_code = check_context(str(contig.seq).upper(), pos, context)
                if return_code == 0: #means we found a hit
                    mapping_dict[pos + 1  + offset] = pos + 1
                    offset +=  1
                    if nt == 'G':
                        offset += 1
        return mapping_dict
    elif context.upper() == 'CHG':
        sequence = str(contig.seq).upper().replace('CG', 'XY')
        converted_sequence = str()
        mapping_dict = {}
        i = 0
        while True:
            if len(sequence[i:]) < 3:
                break
            if sequence[i] in 'CX' and sequence[i + 2] in 'GY':
                if sequence[i] != 'X':
                    mapping_dict[len(converted_sequence) + 1]  = i + 1
                    converted_sequence += 'CG'
                if sequence[i + 2] != 'Y':
                    mapping_dict[len(converted_sequence) + 1]  = i + 3
                    converted_sequence += 'CG'
            else:
                converted_sequence += sequence[i]
            i += 1
        return mapping_dict
    elif context.upper() == 'CHH':
        sequence = str(contig.seq).upper().replace('CG', 'XY')
        converted_sequence = str()
        mapping_dict = {}
        i = 0
        while True:
            if len(sequence[i:]) < 3:
                break
            if sequence[i:i + 2] == 'CG' or (sequence[i] in 'CX' and sequence[i + 2] in 'GY'):
                converted_sequence += sequence[i]
            elif sequence[i] == 'C' and sequence[i + 1] not in ['GY'] and sequence[i + 2] not in ['GY']:
                mapping_dict[len(converted_sequence) + 1]  = i + 1
                converted_sequence += 'CG'
            elif sequence[i] not in 'CX' and sequence[i + 1] not in ['CX'] and sequence[i + 2] == 'G':
                mapping_dict[len(converted_sequence) + 1] = i + 3
                converted_sequence += 'CG'
                i += 2
            else:
                converted_sequence += sequence[i]
            i += 1
        return mapping_dict

def pickle_gff3_entry(in_file, out_file, args, mapping_dict):
    """pickle gff3 file entry with DMPs for merging"""
    DMP_dict = nested_dict()
    with open(in_file,'r') as in_handle:
        header = in_handle.readline()[:-1].split('\t')
        context = in_file.split('.')[-2]
        for line in in_handle:
            split_line = line[:-1].split('\t')
            content = {'context':context}
            for k,v in zip(header,split_line):
                if ',' in v:
                    v = float(v.replace(',','.'))
                content[k] = v
            try:
                if content['diffmeth.p.val']<= float(args.treshold):
                    DMP_dict[content['concatenated_contig']][int(content['concat_contig_pos'])] = content
            except ValueError:
                continue
    with open(out_file,'wb') as out_handle:
        pickle.dump(DMP_dict.to_dict(),out_handle,2)


def rewrite_DMP(in_file, out_name, args, mapping_dict, context):
    currenct_seq = SeqIO.to_dict(SeqIO.parse(args.ref, 'fasta'))
    out_name = out_name.replace(' ','-') + '.%s.tsv' % context
    out_file = os.path.join(args.outputdir, out_name)
    with open(in_file,'r') as in_handle:
        with open(out_file,'w') as out_handle:
            header = in_handle.readline().rstrip('\n').split(',')
            header[1] = 'concatenated_contig'
            header[2] = 'concat_contig_pos'
            header = header[:3] + ['original_contig'] + ['org_contig_pos'] + header[3:]
            out_handle.write('\t'.join(header) + '\n')
            for line in in_handle:
                line = line.replace(',',';')
                line = line.replace('.',',')
                split_line = line.rstrip('\n').split(';')
                contig = split_line[1][3:]
                contig_mapping_dict = context_mapping(currenct_seq[contig], context)
                position = contig_mapping_dict[int(split_line[2])]
                #determine strand
                if currenct_seq[contig][position - 1] == 'G':
                    split_line[3] = '-'
                concat_contig, concat_pos, concat_end = mapping_dict[contig]
                split_line = split_line[:3] + [contig] + [str(position)] + split_line[3:]
                position += int(concat_pos)
                for n,item in enumerate(split_line):
                    if ',' in item:
                        split_line[n] = item[:6]
                split_line[1] = concat_contig
                split_line[2] = str(position)
                out_handle.write('\t'.join(split_line) + '\n')

def make_gff3(base_name, args, mapping_dict):
    """make gff3 output of blast hits for jbrowse display"""
    DMP_dict = nested_dict()
    for context in ['cg','chg','chh']:
        file = os.path.join(args.outputdir,'%s.%s.pickle' % (base_name, context))
        with open(file) as input_handle:
            dict_entry = nested_dict(pickle.load(input_handle))
            DMP_dict.update(dict_entry)
    gff3_dir = os.path.join(args.outputdir,'gff3')
    if not os.path.exists(gff3_dir):
        os.mkdir(gff3_dir)
    gff3_output = open(os.path.join(gff3_dir, base_name.replace('(','').replace(')','').replace('.','') + '.gff3'), 'w')
    gff3_output.write('##gff-version 3.2.1\n')
    for gene, subdict in DMP_dict.items():
        contig = None
        for position,subdict in sorted(subdict.items()):
            if subdict['original_contig'] != contig:
                contig = subdict['original_contig']
                concat_contig, concat_start_pos, concat_end_pos = mapping_dict[contig]
                gff3_output.write('##sequence-region %s %s %s\n' % (concat_contig, concat_start_pos, concat_end_pos))
            out_line = []
            out_line.append('%(concatenated_contig)s' % subdict)            #1 seqid
            out_line.append('RnBeads_%(context)s' % subdict)                #2 source
            out_line.append('5_methylcytosine')                             #3 type 5_methylcytosine  see http://www.sequenceontology.org/browser/current_svn/term/SO:0001918
            out_line.append(str(int('%(concat_contig_pos)s' % subdict)))  #4 start
            out_line.append(str(int('%(concat_contig_pos)s' % subdict)))  #5 end
            out_line.append('%(diffmeth.p.val)s' % subdict)                 #6 score
            out_line.append('%(Strand)s' % subdict)                         #6 strand
            out_line.append('0')                                            #6 Phase . or 0
            attributes = 'ID=%(context)s_%(combinedRank)s;' % subdict
            subdict['diffmeth.p.val'] = '%.2e' % float(subdict['diffmeth.p.val'])
            subdict['diffmeth.p.adj.fdr'] = '%.2e' % float(subdict['diffmeth.p.adj.fdr'])
            attributes += 'Name=%(mean.diff)s p-value:%(diffmeth.p.val)s str:(%(Strand)s);' % (subdict)       #7 Attributes, start with unique ID
            attributes += 'Description=FDR-adjusted p-value:%(diffmeth.p.adj.fdr)s;' % subdict       #7 Attributes, start with unique ID
            attributes += 'Alias=%(concatenated_contig)s_%(concat_contig_pos)s;' % subdict       #7 Attributes, start with unique ID
            attributes += 'Ontology_term=SO:0001918;' % subdict       #7 Attributes, start with unique ID
            for k,v in subdict.items():
                k = k[0].upper() + k[1:]
                attributes += '%s=%s;' % (k,v)
            attributes = attributes[:-1]
            out_line.append(attributes)                                            #6 Phase . or 0
            gff3_output.write('\t'.join(out_line) + '\n')
    os.system('bgzip -f %s' % (os.path.join(gff3_dir, base_name.replace('(','').replace(')','').replace('.','') + '.gff3')) )
    os.system('tabix -p gff %s.gz' % (os.path.join(gff3_dir, base_name.replace('(','').replace(')','').replace('.','') + '.gff3')) )
    # contig = tree._root._children[8]._children[0]._children[2].text
    # try:
    #     concatenated_contig, position = ref_mapping[contig]
    # except KeyError:
    #     return None
    # record['contig_start_pos'] = int(position)
    # record['contig_end_pos'] = int(position) + int(tree._root._children[6].text)
    # record['seqid'] = concatenated_contig
    # blast_hits = tree._root._children[8]._children[0]._children[4]._children
    # record['hits'] = {}
    # names = []
    # for hit in blast_hits:
    #     hit_name = (hit._children[2].text.split('|')[-1].lstrip(' '))
    #     hit_gi = hit._children[2].text.split('|')[1]
    #     ncbi_code = hit._children[2].text.split('|')[3]
    #     if hit_name not in names:
    #         names.append(hit_name)
    #     else:
    #         continue
    #     record['hits'][hit._children[2].text] = {'hsps':{},'gi':hit_gi,'ncbi_code':ncbi_code}
    #     hitdict = record['hits'][hit._children[2].text]
    #     for i,hsp in enumerate(hit._children[5]._children):
    #         record['hits'][hit._children[2].text]['hsps'][i] = {}
    #         dict = record['hits'][hit._children[2].text]['hsps'][i]
    #         for element in hsp._children:
    #             dict[element.tag] = element.text
    #             if element.tag == 'Hsp_query-from':
    #                 #calculate starting position of HSP
    #                 hsp_start_pos = record['contig_start_pos'] + int(element.text)
    #                 dict[element.tag] = hsp_start_pos
    #                 if 'record_start' not in hitdict:
    #                     hitdict['start'] = hsp_start_pos
    #                 else:
    #                     if hsp_start_pos < hitdict['start']:
    #                         hitdict['start'] = hsp_start_pos
    #             elif element.tag == 'Hsp_query-to':
    #                 hsp_end_pos = record['contig_start_pos'] + int(element.text)
    #                 dict[element.tag] = hsp_end_pos
    #                 if 'record_end' not in hitdict:
    #                     hitdict['end'] = hsp_end_pos
    #                 else:
    #                     if hsp_end_pos > hitdict['end']:
    #                         hitdict['end'] = hsp_end_pos
    # gff3 = gff3_record(record)
    # return gff3


def main():
    """main function"""
    args = parse_args()
    mapping_dict = get_mapping(args)
    for context in ['chg','chh','cg']:
        dir = vars(args)['rnbeads_dir_%s' % context]
        DMP_list = get_DMP_list(dir, args)
        for name,file in DMP_list.items():
            out_name = name.replace(' ', '-') + '.%s.tsv' % context
            if not os.path.exists(os.path.join(args.outputdir, out_name)):
                rewrite_DMP(file, name, args, mapping_dict, context)
    #pickle dictionaries for gff3 for merging them.
    for file in os.listdir(args.outputdir):
        if file.endswith('.tsv') and 'based-on-' in file:
            output = '.'.join(file.split('.')[:-1]) + '.pickle'
            if not os.path.exists(os.path.join(args.outputdir, output)):
                pickle_gff3_entry(os.path.join(args.outputdir,file), os.path.join(args.outputdir,output), args, mapping_dict)
    for file in os.listdir(args.outputdir):
        if file.endswith('.cg.pickle') and 'based-on-' in file:
            base_name = '.'.join(file.split('.')[:-2])
            make_gff3(base_name, args, mapping_dict)
    return 0


if __name__ == '__main__':
    return_code = main()