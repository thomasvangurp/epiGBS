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
    parser.add_argument('--bedmapping', help='bed file mapping contigs to ')
    parser.add_argument('--genemapping', help='gff3 file with gene mapping')
    parser.add_argument('-n', '--name', help='name of track in jbrowse',default='test jbrowse entry epiGBS')
    parser.add_argument('--remove',action='store_true', help='remove entry from jbrowse')
    args = parser.parse_args()
    return args


def add_genome(args):
    """Add reference genome to jbrowse"""
    add_genome_script = os.path.join(args.jbrowse_dir,'bin','prepare-refseqs.pl')
    options = ' --fasta %s ' % args.reference
    options += ' --out %s ' % os.path.join(args.jbrowse_dir, 'data', args.name)
    # options += ' --compress ' #compress reference sequences for faster loading
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
    header = barcode_handle.readline()[:-1].split(',')
    for line in barcode_handle:
        split_line = line[:-1].split(',')
        sample_id = split_line[header.index('sampleID')]
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
    # jbrowse_json['tracks'] = jbrowse_json['tracks'][:1]
    for strand in ['watson', 'crick']:
        for name,values in categories.items():
            metadata = dict(values.items())
            for category, value in metadata.items():
                metadata['key'] = '.'.join([name,strand])
                metadata['label'] = '%s_%s_%s' % (category, name, strand)
                metadata['location'] = os.path.join('bam', "%s.%s.bam" % (name, strand))
                template_SNP = {
                     "category": "%s/%s/SNPCoverage" % (category, value),
                     "key" : "%(key)s.SNP" % metadata,
                     "storeClass" : "JBrowse/Store/SeqFeature/BAM",
                     "urlTemplate" : "%(location)s" % metadata,
                     "metadata": metadata,
                     "label" : "%(label)s.SNP" % metadata,
                     "type" : "JBrowse/View/Track/SNPCoverage"
                    }
                template_read = {
                    "category": "%s/%s/Alignments2" % (category, value),
                    "key": "%(key)s.Reads" % metadata,
                    "storeClass": "JBrowse/Store/SeqFeature/BAM",
                    "urlTemplate": "%(location)s" % metadata,
                    "metadata": metadata,
                    "label": "%(label)s.Reads" % metadata,
                    "type": "JBrowse/View/Track/Alignments2"
                }
                if template_SNP['label'] not in [v['label'] for v in jbrowse_json['tracks']]:
                    jbrowse_json['tracks'].append(template_SNP)
                    jbrowse_json['tracks'].append(template_read)
                else:
                    index = [n for n,v in enumerate(jbrowse_json['tracks'])
                             if v['label'] == template_SNP['label']][0]
                    jbrowse_json['tracks'][index] = template_SNP
                    index = [n for n, v in enumerate(jbrowse_json['tracks'])
                             if v['label'] == template_read['label']][0]
                    jbrowse_json['tracks'][index] = template_read
    return jbrowse_json

def add_bigwig(args, jbrowse_json, categories):
    """add bigwig files to json"""
    input_dir = os.path.join(args.input_dir,'bigwig')
    bigwig_outputdir = os.path.join(args.jbrowse_dir, 'data', args.name, 'bigwig')
    if not os.path.exists(bigwig_outputdir):
        os.mkdir(bigwig_outputdir)

    files = os.listdir(input_dir)
    for file in files:
        if '.bw.' not in file:
            continue
        src = os.path.join(args.input_dir, 'bigwig', file)
        dst = os.path.join(bigwig_outputdir, file)
        if os.path.exists(dst):
            os.remove(dst)
        try:
            os.symlink(src, dst)
        except OSError:
            break
    for name, values in categories.items():
        metadata = dict(values.items())
        for category,value in metadata.items():
            metadata['label'] = '%s_%s.bw' % (category,name)
            metadata['key'] = '%s.bw' % (name)
            metadata['location'] = os.path.join('bigwig','%s.bw'%name)
            template_bigwig = {
            "category": "%s/%s/bigwig" % (category,value),
            "key" : "%(key)s" % metadata,
            "label" : "%(label)s" % metadata,
            "style" : { "height" : 50 },
            "storeClass" : "MethylationPlugin/Store/SeqFeature/MethylBigWig",
            "urlTemplate" : "%(location)s" % metadata,
            "metadata" : metadata,
            "type" : "MethylationPlugin/View/Track/Wiggle/MethylPlot"
            }
            if template_bigwig['label'] not in [v['label'] for v in jbrowse_json['tracks']]:
                jbrowse_json['tracks'].append(template_bigwig)
            else:
                index = [n for n, v in enumerate(jbrowse_json['tracks'])
                         if v['label'] == template_bigwig['label']][0]
                jbrowse_json['tracks'][index] = template_bigwig
            #add group specific bw track if this does not yet exist
            metadata_group = {}
            metadata_group['label'] = '%s_%s.bw' % (category, value)
            metadata_group['key'] = '%s_%s.group.bw' % (category, value)
            metadata_group['description'] = 'bigwig for average methylation in group %s_%s' % (category, value)
            metadata_group['location'] = os.path.join('bigwig', '%s_%s.bw' % (category,value))
            template_bigwig = {
                "category": "%s/%s/bigwig" % (category, value),
                "key": "%(key)s" % metadata_group,
                "label": "%(label)s" % metadata_group,
                "style": {"height": 50},
                "storeClass": "MethylationPlugin/Store/SeqFeature/MethylBigWig",
                "urlTemplate": "%(location)s" % metadata_group,
                "metadata": metadata_group,
                "type": "MethylationPlugin/View/Track/Wiggle/MethylPlot"
            }
            if template_bigwig['label'] not in [v['label'] for v in jbrowse_json['tracks']]:
                jbrowse_json['tracks'].append(template_bigwig)
            else:
                index = [n for n, v in enumerate(jbrowse_json['tracks'])
                         if v['label'] == template_bigwig['label']][0]
                jbrowse_json['tracks'][index] = template_bigwig
    return jbrowse_json

def add_bed(args, jbrowse_json):
    """add bed file with mappings of contigs to concatenated refs"""
    add_bed_script = os.path.join(args.jbrowse_dir, 'bin', 'flatfile-to-json.pl')
    options = ' --bed  %s ' % args.bedmapping
    options += ' --trackLabel %s.contig_mapping' % args.name
    options += ' --category Annotation'
    options += ' --out %s ' % os.path.join(args.jbrowse_dir, 'data', args.name)
    # options += ' --compress '  # compress reference sequences for faster loading
    cmd = add_bed_script + options
    if  '%s.contig_mapping' % args.name not in [v['label'] for v in jbrowse_json['tracks']]:
        os.system(cmd)
    jbrowse_json = parse_json(args)
    bed_track_index = [n for n,track in enumerate(jbrowse_json['tracks']) \
                       if track['key'] == args.name + '.contig_mapping'][0]
    bed_track = jbrowse_json['tracks'][bed_track_index]
    #change bed track annotation record
    bed_track['category'] = 'Annotation'
    jbrowse_json['tracks'][bed_track_index] = bed_track
    return jbrowse_json

def add_gff3(args, jbrowse_json):
    """add bed file with mappings of contigs to concatenated refs"""
    add_gff3_script = os.path.join(args.jbrowse_dir, 'bin', 'flatfile-to-json.pl')
    options = ' --gff  %s ' % args.genemapping
    options += ' --trackLabel %s.gene_mapping' % args.name
    options += ' --out %s ' % os.path.join(args.jbrowse_dir, 'data', args.name)
    #options += ' --compress '  # compress reference sequences for faster loading
    cmd = add_gff3_script + options
    if  '%s.gene_mapping' % args.name not in [v['label'] for v in jbrowse_json['tracks']]:
        os.system(cmd)
    jbrowse_json = parse_json(args)
    track_index = [n for n, track in enumerate(jbrowse_json['tracks']) \
                       if track['key'] == args.name + '.gene_mapping'][0]
    track = jbrowse_json['tracks'][track_index]
    track['category'] = 'Annotation'
    track['style']['className'] = 'cds'
    track['style']['description'] = 'e-value'
    track['style']['label'] = 'name'
    track['track'] = "CDS"
    track['type'] = "JBrowse/View/Track/CanvasFeatures"
    track['trackType'] = "CanvasFeatures"
    track["feature"] = ["match_part:bare_predicted","protein_match:predicted"]
    test = {
        "maxFeatureScreenDensity": 400,
        "maxFeatureGlyphExpansion": 500,
        "maxHeight": 600,
        "style": {
            "_defaultHistScale": 4,
            "_defaultLabelScale": 30,
            "_defaultDescriptionScale": 120,
            "showLabels": 'true',
            "showTooltips": 'true',
            "className": "cds",
            "linkTemplate": "http://www.ncbi.nlm.nih.gov/gquery/?term={name}-{start}-{end}"
        },
        "displayMode": "normal",
        "events": {},
        "menuTemplate": [
            {
                "label": "View details",
                "title": "{type} {name}",
                "action": "contentDialog",
                "iconClass": "dijitIconTask"
            },
            {
                "iconClass": "dijitIconFilter"
            }
        ],
        "autocomplete": "all",
        "track": "CDS",
        "key": "CanvasFeatures - mixed",
        "feature": [
            "CDS:bare_predicted",
            "mRNA:exonerate",
            "mRNA:predicted"
        ],
        "storeClass": "JBrowse/Store/SeqFeature/NCList",
        "trackType": "CanvasFeatures",
        "urlTemplate": "tracks/CDS/{refseq}/trackData.json",
        "compress": 0,
        "category": "Genes",
        "type": "JBrowse/View/Track/CanvasFeatures",
        "label": "CDS",
        "baseUrl": "http://www.pristionchus.org/jbrowse/JBrowse-1.10.9/sample_data/json/volvox/",
        "metadata": {}
    }
    #change the json entry, add appropriate category etc.
    jbrowse_json['tracks'][track_index] = track
    return  jbrowse_json


def add_DMP(args, categories):
    """Add bgzipped gff3 files containing DMPs"""
    gff3_dir = os.path.join(args.input_dir,'gff3')
    existing_tracks = []
    #make output dir for symlinks
    output_dir = os.path.join(args.jbrowse_dir, 'data', args.name, 'gff3')
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)
    with open(os.path.join(args.jbrowse_dir,'data',args.name, 'tracks.conf')) as track_handle:
        for line in track_handle:
            if '[tracks.' in line.replace(' ',''):
                line = line.replace(' ','')
                track_name = line[8:line.rindex(']')]
                if track_name not in existing_tracks:
                    existing_tracks.append(track_name)
    with open(os.path.join(args.jbrowse_dir, 'data', args.name, 'tracks.conf'),'a') as track_handle:
        for file in os.listdir(gff3_dir):
            if file.endswith(('.gz')):
                file = file.replace('.','').replace('gff3gz','.gff3.gz')
                g1, vs, g2,based,on,category = file.split('-')
                category = category[:-len('.gff3.gz')]
                dict = {'g1':g1,'g2':g2,'category':category,'file':file}
                dict['trackname'] = 'DMP-' + file[:-len('.gff3.gz')]
                dict['key'] = 'DMP-' + file[:file.index('-based')]
                template = """\n[ tracks . %(trackname)s ]
                storeClass = JBrowse/Store/SeqFeature/GFF3Tabix
                urlTemplate = gff3/%(file)s
                tbiUrlTemplate = gff3/%(file)s.tbi
                type = CanvasFeatures
                style.color = {variation_color}
                metadata.description = DMPs between %(g1)s and %(g2)s based on %(category)s
                category = %(category)s/%(g1)s/DMP
                key = %(key)s\n""" % dict
                template += """\n[ tracks . %(trackname)s_DMP ]
                storeClass = JBrowse/Store/SeqFeature/GFF3Tabix
                urlTemplate = gff3/%(file)s
                tbiUrlTemplate = gff3/%(file)s.tbi
                type = CanvasFeatures
                style.color = {variation_color}
                metadata.description = DMPs between %(g1)s and %(g2)s based on %(category)s
                category = DMP/%(category)s
                key = %(key)s\n""" % dict
                if dict['trackname'] not in existing_tracks:
                    track_handle.write(template.replace('                ',''))
                #now make or update symlinks
                file_in = os.path.join(gff3_dir, file)
                file_out = os.path.join(output_dir,file)
                if not os.path.exists(file_out):
                    os.system('ln -s %s %s' % (file_in, file_out))
                    os.system('ln -s %s.tbi %s.tbi' % (file_in, file_out))



def index_name(args):
    """index names in jbrowse directory"""
    index_names_script = os.path.join(args.jbrowse_dir, 'bin', 'generate-names.pl')
    options = ' --out  %s ' % os.path.join(args.jbrowse_dir, 'data', args.name)
    options += ' --incremental'
    options += ' --mem 2560000000' #2.5 gig allowed
    cmd = index_names_script + options
    os.system(cmd)
    return 0

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
    json = add_bed(args, json)
    json = add_gff3(args, json)
    # index_name(args)
    add_bam(args, json, categories)
    add_bigwig(args, json, categories)
    add_DMP(args, categories)
    write_json(json,args)
    return 0


if __name__ == '__main__':
    main()