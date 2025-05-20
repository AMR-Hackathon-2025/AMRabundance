#!/usr/bin/env python

# python merge_presence_absence.py /shared/team/efrat/example_pres_abs test_combined_pres_abs.tsv

import os
import glob
import sys

if len(sys.argv) != 3:
    print(f"Usage: merge_presence_absence.py <input_folder> <output_file>")
    sys.exit(1)

input_folder = sys.argv[1]
output_file = sys.argv[2]

tsv_files = glob.glob(os.path.join(input_folder, "*.tsv"))

with open(output_file, 'w') as outfile:
    for i, fname in enumerate(tsv_files):
        with open(fname) as infile:
            for j, line in enumerate(infile):
                if i > 0 and j == 0:
                    continue
                outfile.write(line)
