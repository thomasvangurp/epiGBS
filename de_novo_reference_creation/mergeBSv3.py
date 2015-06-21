#!/usr/bin/env python
"Merges Bisulphite reads by finding reciprocal pairs for PE watson / crick sequencing reads"
import sys
import Bio.SeqIO as SeqIO
import pysam
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from pyfaidx import Fasta
import argparse
#seq_in = SeqIO.parse(open('/Volumes/data/epiGBS/test_scabi/join_all.fq', 'r'), 'fastq')


def parse_options():
    """Parses command line options"""
    parser = argparse.ArgumentParser(description='Process input files')
    # parser.add_argument('-r', '--reference', type=str, nargs='?', default=None,
    #                 help='reference genome input.')
    parser.add_argument("-s","--seqin",
                      help = "Combined watson and crick file")
    parser.add_argument("-c","--crickin",dest = "crick",
                      help = "Crick fasta for CRICK_MAX")
    parser.add_argument("--clusters",dest = "clusters",
                      help = "uc input file to make sam output")
    parser.add_argument("--samout",dest = "samout",
                      help = "Sam output file")
    args = parser.parse_args()
    return args


def count_crickMAX(args):
    """Count the number of sequences in the Crick fasta file"""
    with open(args.crick, 'r') as crick_in:
        count = 0
        for line in crick_in:
            if line.startswith('>'):
                count +=1
    return count


if __name__ == '__main__':
    args = parse_options()
    CRICK_MAX =  count_crickMAX(args)
    print "now starting Fasta import"
    seq_in = Fasta(args.seqin)
    print "done with Fasta import"
    clusters = open(args.clusters)
    outsam = args.samout


# path = '/Volumes/data/epiGBS/Baseclear/Athal/'
# path = '/Volumes/data/epiGBS/DNAVISION/Project_DNA11032___140919_SN170_0407_AC52R6ACXX/Sample_DNA11032-001-L1/output/seqykJJfz/scabiosa/'
# path = '/tmp/'
# path = '/Volumes/data/epiGBS/FINAL/Scabiosa/BASECLEAR/'
# seq_in = Fasta(path+'Scabiosa_combined.fa')
#fasta_in = SeqIO.parse(open('/tmp/test.fa', 'r'), 'fasta')
seq_in_keymap = {}
for key in seq_in.keys():
    seq_in_keymap[key.split(';')[0]] = key
    faidx_rec = seq_in[key]

# clusters = open(path +'derep.uc', 'r')
# outsam = path+'derep_out.sam'

#clusters = open('/Volumes/data/epiGBS/test_scabi/cluster_sorted_a.uc', 'r')
#out_fa = open('/Volumes/data/epiGBS/test_scabi/output3.fa', 'w')
#outsam = '/Volumes/data/epiGBS/test_scabi/output3.sam'

#seq_in = SeqIO.parse(open('/Volumes/data/galaxy/database/files/009/dataset_9152.dat', 'r'), 'fastq')
#clusters = open('/Volumes/data/epiGBS/test_scabi/cluster_923.uc', 'r')
#cluster_records = pickle.load(open( "/tmp/save.p", "rb" ))
#
#print 'boe'


def make_bio_object(id, seq_in):
    """Recreates seq-record from pyfaidx record"""
    key = seq_in_keymap[id.split(';')[0]]
    faidx_rec = seq_in[key]
    seq = faidx_rec[0:faidx_rec.__len__()].seq
    record = SeqRecord(Seq(seq), id=key, name= faidx_rec.name, description=faidx_rec.name)
    return record
    

def create_SAM_header(clusters):
    """Create SAM header for bisulphite read mapping"""
    header = { 'HD': {'VN': '1.0'},
            'SQ': [], 'RG':[]}
    header_SQ = {}
    count = 0
    for i in ['crick', 'watson']:
        header['RG'].append({'ID':i, 'SM':i})
    for line in clusters:
        if line.startswith('S'):
            name, length = line.split('\t')[1:3]
            cluster_dict = {'LN':int(length), 'SN':name}
            header_SQ[int(cluster_dict['SN'])] = cluster_dict['LN']
    for k, v in sorted(header_SQ.items()):
        header['SQ'].append({'LN':v, 'SN':str(k)})
    clusters.seek(0)
    return header


#seq_in, clusters, out_fa, outsam = sys.argv[1:]
#clusters = open(clusters, 'r')
#seq_in = SeqIO.parse(open(seq_in, 'r'), 'fastq')
#out_fa = open(out_fa, 'w')



#cluster_records = pickle.load(open( "/tmp/save.p", "rb" ))
#name = outsam.getrname(int(name))
#print 'boe'

count = 0
cluster_records = {}


def add_to_cluster(cluster_instance, cluster_records):
    """Add a read to cluster_records"""
    cluster = cluster_instance.cluster
    seq = cluster_instance.spaced_READ()
    if cluster in cluster_records:
        for pos, nt in enumerate(seq):
            #spaces represent non-covered bases.
            if nt == ' ': 
                continue
            try:
                try:
                    cluster_records[cluster][cluster_instance.read_type][pos][nt]+=1
                except KeyError:
                    cluster_records[cluster][cluster_instance.read_type][pos][nt]=1
            except KeyError:
                try:
                    cluster_records[cluster][cluster_instance.read_type][pos]={nt:1}
                except KeyError:
                    cluster_records[cluster][cluster_instance.read_type] = {}
                    cluster_records[cluster][cluster_instance.read_type][pos]={nt:1}
    else:
        #Add a dictionary in which the nucleotide composition of the cluster is saved
        cluster_records[cluster] = {}
        for pos, nt in enumerate(seq):
            if nt != ' ':
                try:
                    cluster_records[cluster][cluster_instance.read_type][pos] = {nt:1}
                except KeyError:
                    cluster_records[cluster][cluster_instance.read_type] = {}
                    cluster_records[cluster][cluster_instance.read_type][pos] = {nt:1}
            else:
                try:
                    cluster_records[cluster][cluster_instance.read_type][pos] = {}
                except KeyError:
                    cluster_records[cluster][cluster_instance.read_type] = {}
                    cluster_records[cluster][cluster_instance.read_type][pos] = {}
    return cluster_records


def  parse_USEARCH_cigar(cigar, seq_len):
    """Parses usearch specific cigar and return appropriate
    tuple output
    USEARCH uses an alternative output format option in which I and D
    are interchanged
    
    0123456
    MDINSHP
    
    """
    DECODE = 'MDINSHP'
    _ENCODE = dict( (c,i) for (i, c) in enumerate(DECODE) )
    result = []
    n = ''
    total_M = 0 #to keep track of total matching characters to not exceed limit.
    pos = 0
    for c in cigar:
        if c.isdigit():
            n += c
        elif c in _ENCODE:
            pos += 1
            if n == '':
                n='1'
            if pos > 1:
                if total_M +int(n) > seq_len:
                    n = seq_len - total_M
                    total_M += int(n)
                else:
                    if c != 'I':
                        total_M += int(n)
            elif c in 'DM':
                if total_M +int(n) > seq_len:
                    n = seq_len - total_M
                    total_M += int(n)
                else:
                    total_M += int(n)
            result.append( (_ENCODE[c], int(n)) )
            if int(n) == seq_len and c == 'M':
                return result
            n = ''
    return result

# cigar = '91M2I63M'
# seq_len = 154
# parse_USEARCH_cigar(cigar, seq_len)

def create_SAM_line(cluster_obj):
    """Create SAM record for read mapping."""   
    read = cluster_obj.read
    cigar = cluster_obj.cigar_usearch
    record = pysam.AlignedRead()
    record.qname = read.id
    record.tid = int(cluster_obj.cluster)
    if cluster_obj.read_type == 'watson' and cluster_obj.direction == '-':
        record.seq = str(read.seq.reverse_complement())
        record.flag = long(16)
    elif cluster_obj.read_type == 'crick' and cluster_obj.direction == '+':
        record.seq = str(read.seq)
    elif cluster_obj.read_type == 'crick' and cluster_obj.direction == '-':
        record.seq = str(read.seq.reverse_complement())
        record.flag = long(16)
    elif cluster_obj.read_type == 'watson' and cluster_obj.direction == '+':
        record.seq = str(read.seq)
    elif cluster_obj.direction == '*':
        record.seq = str(read.seq)
#    record.qual = read.format('fastq')[:-1].split('\n')[-1]
    record.qual = 'H'*len(record.seq)
    #TODO: properly call mapq using new algorithm.
    record.mapq = len(cluster_obj.read.seq) 
    if cigar in '*=':
        record.cigar = ((0, len(read.seq), ), )
    else:
        cigar_tuple = parse_USEARCH_cigar(cluster_obj.cigar_usearch, len(read.seq))
        if cigar_tuple[0][0] == 2:
            record.pos = cigar_tuple[0][1]
            cigar_tuple = cigar_tuple[1:]
        if cigar_tuple[0][0] == 0:
            if int(cigar_tuple[0][1]) > len(read.seq):
                record.cigar = ((0, len(read.seq), ), )
                record.tags = ( ("NM", 1),("RG", "BS1") )
                return record
        record.cigar = cigar_tuple
#    if int(record.qname) < CRICK_MAX:
#        rg_type = "crick"
#    else:
#        rg_type = "watson"
    tags = []
    for tag in read.description.split('\t')[1:]:
        name, type, content = tag.split(':')
        tags.append(tuple((name, content)))
#    record.tags = tuple(tags)
    if cluster_obj.read_type == 'watson':
        record.tags = (('RG', 'watson', ), )
    else:
        record.tags = (('RG', 'crick', ), )
    return record
    


count = 0

class Cluster_obj:
    """Contains usearch cluster objects"""
    def __init__(self, line):
        """Initializes a class instance"""
        split_line = line.split()
        self.type = split_line[0] #type of cluster object
        cluster_tid_nr = str(outsam.gettid(split_line[1]))
        self.cluster = cluster_tid_nr #split_line[1] #target cluster
        self.direction = split_line[4]#forward or reverse match to centroid.
        if self.direction == '*':
            self.direction = '+'
        self.cigar_usearch = split_line[7] #usearch type cigar
        self.readname = split_line[8]
        self.read = None #seq read object
        #Determines if the read is watson or crick.
        self.read_type = None #watson or crick read?
    def SAM_record(self):
        """Returns a SAM record of cluster_object"""
        sam_record = create_SAM_line(self)
        return sam_record
    def spaced_READ(self):
        """"Returns a string with contig aligned read values
        the string contains spaces if there is no coverage for a certain
        contig position"""
        if self.direction == '-':
            rv_read = read.reverse_complement()
            rv_read.description = self.read.description
            rv_read.id = self.read.id
            rv_read.name = self.read.name
            self.read = rv_read
#            self.read.seq = self.read.seq.reverse_complement()
        if self.cigar_usearch in ['*', '=']:
            return str(self.read.seq)
        else:
            pos = 0
            out = ''
            cigar = self.SAM_record().cigar
            #Add the appropriate number of leading spaces.
            out += self.SAM_record().pos * ' '
            for item in cigar:
                operation, length = item
                #Match
                if operation == 0:
                    out+=str(self.read.seq)[pos:pos+length]
                    pos+=length
                #Insertion
                if operation == 2:
                    out+=' '*length
                #deletion
                if operation == 1:
                    pos += length
            return out



#parse fasta input of valid clusters with min_size = 2!
cluster_ids = set()
total_lines = 0
for cluster in clusters:
    if cluster.startswith('C'):
        split_line = cluster.split('\t')
        if int(split_line[2])>1:
            cluster_ids.add(split_line[1])
    else:
        total_lines+=1
clusters.seek(0)
print len(cluster_ids)

header = create_SAM_header(clusters)
outsam = pysam.Samfile(outsam,'wh', header=header)
n = 0
cluster_list = seq_in.keys()

for line in clusters:
    #clusters contains an ordered list of sequences which should be processed
    #Forward reads /1 or merged fastq reads always represent crick reads. 
    #These have number below CRICK_MAX
    if not n%1000 and n:
        print "processed %s out of %s lines"%(n, total_lines)
#        make_ref(cluster_records)
    if line.startswith('S') or line.startswith('H'):
        n+=1
    if line.split('\t')[1] not in cluster_ids:
        continue
    cluster_instance = Cluster_obj(line)
    #get read with name count.
    read = make_bio_object(cluster_instance.readname, seq_in)
    cluster_pos = cluster_list.index(read.name)
    cluster_instance.read = read
    #Centroids of clusters are always in direction +, so the type depends only on 
    if cluster_instance.type == 'C':
        cluster_instance.direction = '+'
    #Now we will only encounter crick reads, as they are >crick first
    if cluster_pos  > CRICK_MAX:
       if cluster_instance.direction == '-':
            cluster_instance.read_type = "crick"
       elif cluster_instance.direction == '+':
            cluster_instance.read_type = "watson"
    #Now we will only encounter Watson reads, as they are >crick first
    elif cluster_pos <= CRICK_MAX:
        if cluster_instance.direction == '-':
            cluster_instance.read_type = "watson"
        elif cluster_instance.direction == '+':
            cluster_instance.read_type = "crick"

#    if cluster_instance.type not in ['C', 'N']: #see http://www.drive5.com/usearch/manual/ucout.html
    sam_record = create_SAM_line(cluster_instance)
    outsam.write(sam_record)
outsam.close()
