#!/usr/bin/env python
"""add DMPs to jbrowse as gff3 and bigbed files"""
import argparse
import os
import gzip
import xml.etree.cElementTree as ET
from Bio import SeqIO


def parse_args():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(description='Concatenate reference for display in jbrowse')
    parser.add_argument('-x', '--xml', help='xml input with annotation')
    parser.add_argument('--rnbeads_dir', help='Output directory from Rnbeads')
    parser.add_argument('-c', '--context', help='Nucleotide context for DMPs')
    parser.add_argument('-b', '--bed', help='bed file with contig to concatenated mapping')
    parser.add_argument('-r', '--ref', help='reference sequence for contigs')
    parser.add_argument('-o', '--outputdir', help='output directory for gff3 and bigbed files')
    parser.add_argument('--treshold', help='minimum p value for DMP',default='0.0001')
    args = parser.parse_args()
    if not os.path.exists(args.outputdir):
        os.mkdir(args.outputdir)
    return args

def get_mapping(args):
    """parse bed file for mapping contigs to concatenated ref"""
    mapping_dict = {}
    with open(args.bed,'r') as handle:
        for line in handle:
            target,target_start,target_end,name = line.rstrip('\n').split('\t')
            mapping_dict[name] = (target,target_start)
    return mapping_dict


def get_DMP_list(args):
    """parse HTML of Rnbeads to get DMP list"""
    handle = open(os.path.join(args.rnbeads_dir,'differential_methylation.html'),'r')
    DMP_list = {}
    for line in handle:
        if '<li>' in line and '.csv' in line:
            file_location = os.path.join(args.rnbeads_dir, line[line.index('"')+1:line.rindex('"')])
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
    if actual_context != context:
        # print "Context mismatch!"
        return 1
    else:
        return 0


def context_mapping(contig, args):
    """reverse mapping of 'fake' RNbeads positions to 'real' positions"""
    mapping_dict = {}
    offset = 0
    if args.context == 'CG':
        for pos,nt in enumerate(str(contig.seq.upper())):
            if nt == 'C' or nt == 'G':
                return_code = check_context(str(contig.seq).upper(),pos,args.context)
                if return_code == 0: #means we found a hit
                    mapping_dict[pos + 1  + offset] = pos + 1
                    offset +=  1
                    if nt == 'G':
                        offset += 1
        return mapping_dict


def make_gff3_entry(in_file, out_file, args):
    """make gff3 file entry for DMP"""
    DMP_dict = {}
    with open(in_file,'r') as in_handle:
        header = in_handle.readline()[:-1].split(',')
        with open(out_file,'w') as out_handle:
            for line in in_handle:
                split_line = line[:-1].split(',')
                content = {}
                for k,v in zip(header,split_line):
                    content[k] = v
                if float(content['diffmeth.p.val']) <= float(args.treshold):
                    try:
                        DMP_dict[content['Chromosome']][content['Start']]
                    except KeyError:
                        pass
    gff3_output = '##sequence-region %(seqid)s %(contig_start_pos)s %(contig_end_pos)s\n'


def rewrite_DMP(in_file,out_name, args, mapping_dict):
    currenct_seq = SeqIO.to_dict(SeqIO.parse(args.ref, 'fasta'))
    out_name = out_name.replace(' ','_') + '.tsv'
    out_file = os.path.join(args.outputdir,out_name)
    out_file.replace(' ','_')
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
                contig_mapping_dict = context_mapping(currenct_seq[contig], args)
                try:
                    position = contig_mapping_dict[int(split_line[2])]
                except KeyError:
                    continue
                #determine strand
                if currenct_seq[contig][position - 1] == 'G':
                    split_line[3] = '-'
                concat_contig, concat_pos = mapping_dict[contig]
                split_line = split_line[:3] + [contig] + [str(position)] + split_line[3:]
                position += int(concat_pos)
                for n,item in enumerate(split_line):
                    if ',' in item:
                        split_line[n] = item[:6]
                split_line[1] = concat_contig
                split_line[2] = str(position)
                out_handle.write('\t'.join(split_line) + '\n')

def gff3_record(record):
    """make gff3 record output"""
    gff3_output = '##sequence-region %(seqid)s %(contig_start_pos)s %(contig_end_pos)s\n' % record
    for name,hit in record['hits'].items():
        hit['ID'] = name.split('|')[1] + '.' + record['seqid'] + '.%s' % hit['start']
        hit['name'] = name.split('|')[-1].lstrip(' ')
        hit['e-value'] = 1
        hit['seqid'] = record['seqid']
        if hit['hsps'][0]['Hsp_query-frame'][0] == '-':
            hit['strand'] = '-'
        else:
            hit['strand'] = '+'
        #TODO:check if phase and score are correct
        hit['phase'] = '0'
        hit['score'] = '.'
        #now add attributes
        hsps = ''
        for hsp in hit['hsps'].values():
            hsp['Parent'] = hit['ID']
            hsp['seqid'] = hit['seqid']
            hsp['strand'] = hit['strand']
            hsp['Hsp_query-frame'] = hsp['Hsp_query-frame'].lstrip('-')
            #todo: check if this phase is correct
            if hit['phase'] == '0':
                hit['phase'] = hsp['Hsp_query-frame']
            hsps += '%(seqid)s\tdiamond\tmatch_part\t%(Hsp_query-from)s\t%(Hsp_query-to)s\t' % hsp
            hsps += '%(Hsp_bit-score)s\t%(strand)s\t%(Hsp_query-frame)s\t' % hsp
            attributes = 'Parent=%(Parent)s' % hsp
            for key,value in hsp.items():
                attributes += ';%s=%s' % (key.lower(),value)
                if key == 'Hsp_evalue':
                    if float(value) < float(hit['e-value']):
                        hit['e-value'] = value
            hsps += attributes + '\n'
        gff3_output += '%(seqid)s\tdiamond\tprotein_match\t%(start)s\t%(end)s\t%(score)s\t%(strand)s\t%(phase)s\t' % hit
        gff3_output += 'ID=%(ID)s;name=%(name)s;gi=%(gi)s;ncbi_code=%(ncbi_code)s' % hit
        gff3_output += ';e-value=%s\n' % hit['e-value']
        gff3_output += hsps
    return gff3_output

def make_gff3_genes(ref_mapping, args, tree):
    """make gff3 output of blast hits for jbrowse display"""
    record = {}
    contig = tree._root._children[8]._children[0]._children[2].text
    try:
        concatenated_contig, position = ref_mapping[contig]
    except KeyError:
        return None
    record['contig_start_pos'] = int(position)
    record['contig_end_pos'] = int(position) + int(tree._root._children[6].text)
    record['seqid'] = concatenated_contig
    blast_hits = tree._root._children[8]._children[0]._children[4]._children
    record['hits'] = {}
    names = []
    for hit in blast_hits:
        hit_name = (hit._children[2].text.split('|')[-1].lstrip(' '))
        hit_gi = hit._children[2].text.split('|')[1]
        ncbi_code = hit._children[2].text.split('|')[3]
        if hit_name not in names:
            names.append(hit_name)
        else:
            continue
        record['hits'][hit._children[2].text] = {'hsps':{},'gi':hit_gi,'ncbi_code':ncbi_code}
        hitdict = record['hits'][hit._children[2].text]
        for i,hsp in enumerate(hit._children[5]._children):
            record['hits'][hit._children[2].text]['hsps'][i] = {}
            dict = record['hits'][hit._children[2].text]['hsps'][i]
            for element in hsp._children:
                dict[element.tag] = element.text
                if element.tag == 'Hsp_query-from':
                    #calculate starting position of HSP
                    hsp_start_pos = record['contig_start_pos'] + int(element.text)
                    dict[element.tag] = hsp_start_pos
                    if 'record_start' not in hitdict:
                        hitdict['start'] = hsp_start_pos
                    else:
                        if hsp_start_pos < hitdict['start']:
                            hitdict['start'] = hsp_start_pos
                elif element.tag == 'Hsp_query-to':
                    hsp_end_pos = record['contig_start_pos'] + int(element.text)
                    dict[element.tag] = hsp_end_pos
                    if 'record_end' not in hitdict:
                        hitdict['end'] = hsp_end_pos
                    else:
                        if hsp_end_pos > hitdict['end']:
                            hitdict['end'] = hsp_end_pos
    gff3 = gff3_record(record)
    return gff3




def get_xml_record(ref_mapping, args, output_handle):
    """read xml file to yield single record"""
    record_out = ''
    with open(args.xml) as handle:
        while True:
            pos = handle.tell()
            line = handle.readline()
            if not line:
                break
            if pos and line == '<?xml version="1.0"?>\n' and record_out != '':
                tree = ET.ElementTree(ET.fromstring(record_out))
                output = make_gff3_genes(ref_mapping,args,tree)
                if output:
                    output_handle.write(output)
                record_out = ''
                handle.seek(pos)
            else:
                record_out += line
    output_handle.close()

def main():
    """main function"""
    args = parse_args()
    mapping_dict = get_mapping(args)
    DMP_list = get_DMP_list(args)
    for name,file in DMP_list.items():
        # make_gff3_entry(file, name +'.gff3', args)
        rewrite_DMP(file,name,args,mapping_dict)
    return 0


if __name__ == '__main__':
    return_code = main()