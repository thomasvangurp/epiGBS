#!/usr/bin/env python3
import sys
import gzip

input_file = sys.argv[1]
input_handle = gzip.open(input_file)
output_handle = sys.stdout

count_dict = {}
while True:
    lines = []
    for i in range(4):
        lines.append(input_handle.readline().decode())
    try:
        bc = lines[0].split('\t')[3]
    except IndexError:
        break
    try:
        count_dict[bc] += 1
    except KeyError:
        count_dict[bc] = 1
    if count_dict[bc] <= 100000:
        [output_handle.write(l) for l in (lines)]