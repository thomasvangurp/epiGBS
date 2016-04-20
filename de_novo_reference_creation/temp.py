"""Temporary file for testing new uc creation"""
import operator
import subprocess


# watson_input = "/Volumes/3tb-deena/tmp/join_watsonraW6ZM_demethylated.assembled.derep.fastq"
watson_input = "/Volumes/3tb-deena/tmp/join_watsonraW6ZM_demethylated.joined.derep.fa"
# crick_input = "/Volumes/3tb-deena/tmp/join_crickB52Fbr_demethylated.assembled.derep.fastq"
crick_input = "/Volumes/3tb-deena/tmp/join_crickB52Fbr_demethylated.joined.derep.fa"
# output = "/Volumes/3tb-deena/tmp/test_merged_daphnia.fa"
output = "/Volumes/3tb-deena/tmp/test_joined_daphnia.fa"
#inputs are generated after vsearch -derep_fulllength


def make_binary_output(watson_handle,crick_handle,output):
    """make binary output sequence for uc generation"""
    for i,handle in enumerate([watson_handle,crick_handle]):
        if i == 0:
            name_start = '>w_'
        else:
            name_start = '>c_'
        for line in handle:
            if line.startswith('>'):
                try:
                    seq = seq.replace('C','T').replace('G','A')
                    out = name + '\n' + seq + '\n'
                    output.write(out)
                except NameError:
                    pass
                name = '%s'%name_start
                seq = ''
            else:
                name += line.rstrip('\n')
                seq += line.rstrip('\n')


def get_ref(clusters):
    """Generate reference from sequences that cluster together"""
    output_count = {}
    id = int(clusters[0][1]) + 1
    transform = {'A':'T','C':'G','G':'C','T':'A','N':'N'}
    for cluster in clusters:
        direction = cluster[4]
        seq = cluster[8][2:]
        strand = cluster[8][0]
        if direction != '-':
            pass
        else:
            #get reverse complement
            seq = [transform[nt] for nt in seq][::-1]
            if strand == 'w':
                strand = 'c'
            else:
                strand = 'w'
        for i,nt in enumerate(seq):
            try:
                output_count[i][strand][nt] += 1
            except KeyError:
                if i not in output_count:
                    output_count[i] = {strand:{nt:1}}
                if strand not in output_count[i]:
                    output_count[i][strand] = {nt:1}
                if nt not in output_count[i][strand]:
                    output_count[i][strand][nt] = 1
    output_fasta = '>%s\n'%id
    for i in sorted(output_count.keys()):
        try:
            watson_nt = max(output_count[i]['w'].iteritems(), key=operator.itemgetter(1))[0]
            crick_nt  = max(output_count[i]['c'].iteritems(), key=operator.itemgetter(1))[0]
        except KeyError:
            return None
        if watson_nt == crick_nt:
            output_fasta += watson_nt
        elif watson_nt == 'G' and crick_nt == 'A':
            output_fasta += watson_nt
        elif crick_nt == 'C' and watson_nt == 'T':
            output_fasta += crick_nt
        elif watson_nt == 'N' and crick_nt != 'N':
            output_fasta += crick_nt
        elif watson_nt != 'N' and crick_nt == 'N':
            output_fasta += watson_nt
        else:
            output_fasta += 'N'
    output_fasta += '\n'
    return output_fasta



def make_ref_from_uc(uc,ref):
    """make reference directly from uc output"""
    uc_handle = open(uc,'r')
    ref_handle = open(ref,'w')
    clusters = []
    for line in uc_handle:
        if line.startswith('H') or line.startswith('S'):
            split_line = line.split('\t')
            cluster_id = int(split_line[1])
            try:
                if cluster_id != int(clusters[-1][1]):
                    ref = get_ref(clusters)
                    if ref:
                        ref_handle.write(ref)
                    clusters = []
            except NameError:
                pass
            except IndexError:
                pass
            clusters.append(split_line)

watson_handle = open(watson_input,'r')
crick_handle = open(crick_input,'r')
output = open(output,'w')
make_binary_output(watson_handle,crick_handle,output)

output_ref = "/Volumes/3tb-deena/tmp/ref_joined_daphnia.fa"
output = "/Volumes/3tb-deena/tmp/test_joined_daphnia.fa"
#get reference sequence directly from uc output
uc = "/Volumes/3tb-deena/tmp/test_joined_daphnia.uc"
cmd = "vsearch -derep_fulllength %s -strand both -uc %s"%(output,uc)
p = subprocess.Popen(cmd,stdout=subprocess.PIPE,stderr=subprocess.PIPE,shell=True,executable='bash')
stdout, stderr = p.communicate()
make_ref_from_uc(uc,output_ref)

# now starting:   Combine and convert merged watson and crick files using sed
# running:        cat <(sed 's/>/>c/' /Volumes/3tb-deena/tmp/join_crickb85Kba_demethylated.assembled.derep.fastq ) <(sed 's/>/>w/' /Volumes/3tb-deena/tmp/join_watsoncH1pg2_demethylated.assembled.derep.fastq ) |tee >(sed '/^[A,C,G,T,N]*$/s/C/T/g;/^[A,C,G,T,N]*$/s/G/A/g'> /Volumes/3tb-deena/tmp/join_combinedb85Kba_demethylated.assembled.derep.CTGA.fastq) > /Volumes/3tb-deena/tmp/join_combinedb85Kba_demethylated.assembled.derep.fastq
# finished:       Combine and convert merged watson and crick files using sed

output_derep_CTGA = "/Volumes/3tb-deena/tmp/join_combinedb85Kba_demethylated.assembled.derep.CTGA.fastq"