"""Do separate analysis per genetic cluster"""
import argparse
import os
import gzip
import xml.etree.cElementTree as ET

def parse_args():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(description='Concatenate reference for display in jbrowse')
    parser.add_argument('-s', '--samples', help='input sample file with genetic clusters')
    parser.add_argument('-b', '--bed', help='original bed file to filter again per cluster')
    parser.add_argument('-c', '--context', help='contexts for which to do analysis', default ='CG,CHG,CHH')
    parser.add_argument('-r', '--rnbeads script', help='annotation file from blast2go',
                        default='/Users/thomasvangurp/Dropbox/Thomas_NIOO/bioinformatics/epiGBS'+
                                '/github_repository/epiGBS/rnbeads/rnbeads_analysis.py')
    parser.add_argument('-o', '--outputdir', help='output directory for gff3 file')
    parser.add_argument('-s', '--species', help='list of species to run for',
                        )
    args = parser.parse_args()
    return args

def get_mapping(args):
    """parse bed file for mapping contigs to concatenated ref"""
    mapping_dict = {}
    with open(args.bed,'r') as handle:
        for line in handle:
            target,target_start,target_end,name = line.rstrip('\n').split('\t')
            mapping_dict[name] = (target,target_start)
    return mapping_dict


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
                tree = ET.ElementTree(ET.fromstring(record_out.replace('&','%26')))
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
    output_handle = open(os.path.join(args.outputdir,'out.gff3'),'w')
    output_handle.write('##gff-version 3\n')
    get_xml_record(mapping_dict, args, output_handle)
    return 0


if __name__ == '__main__':
    return_code = main()