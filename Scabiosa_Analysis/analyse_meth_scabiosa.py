from nested_dict import nested_dict

groups = {}
with open('groups.tsv','r') as handle:
    header = handle.readline().rstrip('\n').split('\t')
    for line in handle:
        split_line = line.rstrip('\n').split('\t')
        groups[split_line[0]] = split_line[-1]
print ''


def wur_buxton_diff(header, split_line):
    """calculates average difference in meth between WUR and buxton"""
    meth_values = nested_dict()
    for k,v in zip(header,split_line):
        try:
            location,sample,type = k.split('_')
            if sample in groups:
                meth_values[groups[sample]][location][sample][type] = v
        except ValueError:
            pass
    for group in meth_values.keys():
        buxton_meth_values = []
        for ind, meth_dict in meth_values[group]['BUXTON'].items():
            try:
                buxton_meth_values.append(int(meth_dict['methylated'])
                                          / float(meth_dict['total']))
            except TypeError:
                buxton_meth_values.append(None)
            except ValueError:
                buxton_meth_values.append(None)
        wur_meth_values = []
        for ind, meth_dict in meth_values[group]['WUR'].items():
            try:
                wur_meth_values.append(int(meth_dict['methylated'])
                                       / float(meth_dict['total']))
            except TypeError:
                wur_meth_values.append(None)
            except ValueError:
                wur_meth_values.append(None)
        diff = [a - b for a,b in zip(buxton_meth_values, wur_meth_values)
                if a and b]
        if diff == []:
            abs_diff = None
            rel_diff = None
        else:
            abs_diff = sum([abs(v) for v in diff]) / float(len(diff))
            rel_diff = sum(diff) / float(len(diff))
        meth_values[group]['abs_diff'] = abs_diff
        meth_values[group]['rel_diff'] = rel_diff
    return meth_values

output = open('/Users/thomasvangurp/epiGBS/scabi_aug2017/buxton_wur_trans.bed','w')

with open('/Users/thomasvangurp/epiGBS/scabi_aug2017/valid_pairs.filtered.bed') as handle:
    header = handle.readline().rstrip('\n').split('\t')
    output_header = header[:3] + ['C_abs','SD_abs','SR_abs','WW_abs',
                                  'WWSD_abs','WWSR_abs']
    output_header += ['C_rel','SD_rel','SR_rel','WW_rel',
                                  'WWSD_rel','WWSR_rel']
    output.write('\t'.join(output_header) + '\n')
    for line in handle:
        split_line = line[:-1].split('\t')
        diff = wur_buxton_diff(header, split_line)
        output_line = split_line[:3]
        for k,v in sorted(diff.items()):
            try:
                output_line.append('%.3f'%v['abs_diff'])
            except TypeError:
                output_line.append('None')
        for k,v in sorted(diff.items()):
            try:
                output_line.append('%.3f'%v['rel_diff'])
            except TypeError:
                output_line.append('None')
        output.write('\t'.join(output_line) + '\n')
output.close()