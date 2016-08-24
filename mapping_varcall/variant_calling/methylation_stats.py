"""Decribe methylation stats per species"""
import sys

# handle = open(sys.argv[1])
handle = open('/Volumes/3tb-deena/Daphnia/output_DMP/Daphnia.methylation.filtered.bed','r')

header = handle.readline()[:-1].split('\t')

individuals = ['_'.join(i.split('_')[:-1]) for i in header if i.endswith('_total')]
meth_count = {}
count = 0
for line in handle:
    count+=1
    if not count%100000:
        print count
    split_line = line[:-1].split('\t')
    for ind,i in zip(individuals,range(4,len(header),2)):
        try:
            methylated = int(split_line[i])
            total = int(split_line[i+1])
        except ValueError:
            continue
        try:
            meth_perct = methylated / float(total)
        except ZeroDivisionError:
            meth_perct = 0
        if meth_perct > 0.5:
            try:
                meth_count[ind]['meth'] += 1
                meth_count[ind]['total'] += 1
            except KeyError:
                if ind not in meth_count:
                    meth_count[ind] = {'meth':1}
                    meth_count[ind] = {'total':1}
                else:
                    meth_count[ind]['meth'] = 1
                    meth_count[ind]['total'] = 1
        else:
            try:
                meth_count[ind]['total'] += 1
            except KeyError:
                if ind not in meth_count:
                    meth_count[ind] = {'total': 1}
                else:
                    meth_count[ind]['total'] = 1
for k,v in sorted(meth_count.items()):
    print k,'\t','\t'.join(['%s' % (b) for a,b in v.items()])

