#!/usr/bin/env python
# This scripts takes an input file (assumed to be an output of samtools depth command),
#  and it also reads a file with reference sequence lengths. Then we calculate the % of the reference genome 
#  covered by mapped reads. If it is above X% it considers that gene as present.

# For testing (on CLIMB):
# cd /shared/team/efrat/
# python calc_presence_absence.py example_depth_files/ERR2749674.tsv /shared/team/databases/card-data/nucleotide_fasta_protein_homolog_model_contig_lengths.tsv test_calc_presence_absence_V2.tsv

import os
import sys

# Arguments
if len(sys.argv) < 4:
    print("Usage: python calc_presence_absence.py <input_file> <lengths_file> <output_file> [threshold]")
    sys.exit(1)

depth_file = sys.argv[1]
lengths_file = sys.argv[2]
output_file = sys.argv[3]
threshold = float(sys.argv[4]) if len(sys.argv) > 4 else 10.0

# Extract run ID from file name (without .tsv suffix)
run_id = os.path.basename(depth_file).replace(".tsv", "")

# Load sequence lengths
lengths = {}
with open(lengths_file) as f:
    for line in f:
        if line.strip():
            raw_id, length = line.strip().split('\t')
            seq_id = raw_id.split()[0]  # Take only the first part before any whitespace
            lengths[seq_id] = int(length)

# Parse depth file and count covered positions
coverage_counts = {}
with open(depth_file) as f:
    for line in f:
        seq_id, pos, depth = line.strip().split('\t')
        if seq_id not in lengths:
            print(f"Warning: sequence ID '{seq_id}' not found in lengths file. Skipping.")
            sys.exit(1)
        if int(depth) > 0:
            coverage_counts.setdefault(seq_id, 0) # Default: 0
            coverage_counts[seq_id] += 1

# Calculate coverage and determine presence/absence
results = {}
for seq_id in lengths:
    covered = coverage_counts.get(seq_id, 0)
    length = lengths[seq_id]
    percent_covered = (covered / length) * 100
    results[seq_id] = 1 if percent_covered >= threshold else 0

# Write output
with open(output_file, 'w') as out:
    out.write("run_id\tcard_id\tpresence\n")
    for seq_id in sorted(results):
        out.write(f"{run_id}\t{seq_id}\t{results[seq_id]}\n")

print("Done.")
