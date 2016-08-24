"""split watson and crick bam file into sample specific bam files"""
import json
import argparse
import os
import subprocess
import json

def parse_args():
    """Pass command line arguments"""

    parser = argparse.ArgumentParser()
    parser.add_argument('-j', '--jbrowse_dir', help='jbrowse directory', default='/Applications/MAMP/htdocs/jbrowse/')
    parser.add_argument('-i', '--input_dir', help='input directory')
    parser.add_argument('-b', '--barcodes', help='Barcode directory for determining categories')
    parser.add_argument('-c', '--categories', help='comma separated list of categories to include')
    parser.add_argument('-r', '--reference', help='reference fasta file')
    parser.add_argument('-n', '--name', help='name of track in jbrowse',default='test jbrowse entry epiGBS')
    args = parser.parse_args()
    return args


def add_genome(args):
    """Add reference genome to jbrowse"""
    add_genome_script = os.path.join(args.jbrowse_dir,'bin','prepare-refseqs.pl')
    options = ' --fasta %s ' % args.reference
    options += ' --out %s ' % os.path.join(args.jbrowse_dir, 'data', args.name)
    options += ' --trackLabel %s ' % args.name
    cmd = add_genome_script + options
    os.system(cmd)


def make_categories(args):
    """add categories required for listing files in jbrowse output"""
    categories = {}
    try:
        barcode_handle = open(args.barcodes, 'r')
    except OSError:
        raise OSError('File %s does not exist' % args.barcode)
    header = barcode_handle.readline()[:-1].split('\t')
    for line in barcode_handle:
        split_line = line[:-1].split('\t')
        sample_id = split_line[header.index('Sample')]
        categories[sample_id] = {}
        for k, v in zip(header, split_line):
            if k not in args.categories.split(','):
                continue
            categories[sample_id][k] = v

    return categories


def parse_json(args):
    """add or modify bam files to json"""
    json_in = os.path.join(args.jbrowse_dir, 'data', args.name, 'trackList.json')
    json_out = json.load(open(json_in, 'r'))
    return json_out


def add_bam(args, jbrowse_json, categories):
    """add bam files to json"""
    #make directory for symlinking bam files
    bam_path = os.path.join(args.jbrowse_dir, 'data', args.name, 'bam')
    if not os.path.exists(bam_path):
        os.mkdir(bam_path)
    files = os.listdir(os.path.join(args.input_dir,'bam'))
    for file in files:
        src = os.path.join(args.input_dir,'bam', file)
        dst = os.path.join(bam_path, file)
        if os.path.exists(dst):
            os.remove(dst)
        try:
            os.symlink(src, dst)
        except OSError:
            break

    for strand in ['watson', 'crick']:
        for name,values in categories.items():
            metadata = dict(values.items())
            metadata['name'] = '.'.join([name,strand])
            #TODO: add proper categories
            metadata['location'] = os.path.join('bam', "%s.%s.bam" % (name, strand))
            template_SNP = {
                 "category" : "%(history)s/SNPCoverage" % metadata,
                 "key" : "%(name)s.SNP" % metadata,
                 "storeClass" : "JBrowse/Store/SeqFeature/BAM",
                 "urlTemplate" : "%(location)s" % metadata,
                 "metadata": metadata,
                 "label" : "%(name)s.SNP" % metadata,
                 "type" : "JBrowse/View/Track/SNPCoverage"
                }
            template_read = {
                "category": "%(history)s/Alignments2" % metadata,
                "key": "%(name)s.Reads" % metadata,
                "storeClass": "JBrowse/Store/SeqFeature/BAM",
                "urlTemplate": "%(location)s" % metadata,
                "metadata": metadata,
                "label": "%(name)s.Reads" % metadata,
                "type": "JBrowse/View/Track/Alignments2"
            }
            if template_SNP['key'] not in [v['key'] for v in jbrowse_json['tracks']]:
                jbrowse_json['tracks'].append(template_SNP)
                jbrowse_json['tracks'].append(template_read)
            else:
                index = [n for n,v in enumerate(jbrowse_json['tracks'])
                         if v['key'] == template_SNP['key']][0]
                jbrowse_json['tracks'][index] = template_SNP
                index = [n for n, v in enumerate(jbrowse_json['tracks'])
                         if v['key'] == template_read['key']][0]
                jbrowse_json['tracks'][index] = template_read
    return jbrowse_json


def write_json(jbrowse_json, args):
    """Write json to output location"""
    out_file = open(os.path.join(args.jbrowse_dir, 'data', args.name, 'trackList.json'),'w')
    json.dump(jbrowse_json, out_file, indent=4, sort_keys=True)
    #write dataset id to config
    in_config = open(os.path.join(args.jbrowse_dir, 'data', args.name, 'tracks.conf'),'r')
    config = "[general]\ndataset_id = %s\n" % args.name
    config += in_config.read()
    out_config = open(os.path.join(args.jbrowse_dir, 'data', args.name, 'tracks.conf'),'w')
    out_config.write(config)
    out_config.close()
    #write dataset id to jbrowse config
    jbrowse_config_handle = open(os.path.join(args.jbrowse_dir, 'jbrowse.conf'),'r')
    jbrowse_config = jbrowse_config_handle.read().rstrip('\n')
    if '[datasets.%s]' % args.name not in jbrowse_config:
        jbrowse_config += '\n[datasets.%(name)s]\nurl  = ?data=data/%(name)s/\nname = %(name)s\n' % vars(args)
        jbrowse_config_handle = open(os.path.join(args.jbrowse_dir, 'jbrowse.conf'), 'w')
        jbrowse_config_handle.write(jbrowse_config)
        jbrowse_config_handle.close()




def main():
    """main function"""
    args = parse_args()
    if not os.path.exists(os.path.join(args.jbrowse_dir,'data', args.name, 'seq')):
        add_genome(args)
    categories = make_categories(args)
    json = parse_json(args)
    add_bam(args, json, categories)
    write_json(json,args)
    return 0


if __name__ == '__main__':
    main()