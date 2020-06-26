#!/usr/bin/env/ python
import sys
from Bio.Seq import Seq

def get_gRNA_lib(gRNA_file):
    f = open(gRNA_file)
    gRNA_lib = {}
    i = 0
    for line in f:
        if line[0] != ">":
            gRNA_seq = line.rstrip()[0:-3]  # remove NGG PAM sequence (present in genome but not in plasmid gRNA)

            if gRNA_seq in gRNA_lib:
                raise AssertionError("gRNA library file contains duplicate entries!")
            else:
                gRNA_lib[gRNA_seq] = 0
    f.close()
    return gRNA_lib


def count_gRNA(gRNA_lib, seq_file):
    f = open(seq_file)

    gRNA_count = gRNA_lib.copy()

    for gRNA in gRNA_count:
        gRNA_revcomp = str(Seq(gRNA).reverse_complement())
        f.seek(0)  # go to the begginng of the file
        for line in f:
            # we assume that R1 and R2 reads are always overlapping, hence we only need to look at
            # the R1 file and look for both the oriented and the reverse complement.
            if gRNA in line or gRNA_revcomp in line:
                gRNA_count[gRNA] += 1
            # how could we count sequences with no gRNA?
            # how can we count the sequences from the Mismatch in the Perf files?
            # Probably reusing this script but imputing the mismatch library! 
    f.close
    return gRNA_count


# get filename from command line arguments
# could be improved with getopt
try:
    plasmid_seq_data = sys.argv[1]
    gRNA_lib_file = sys.argv[2]
    csv_out_file = sys.argv[3]
except IndexError:
    raise SystemExit(f"Usage: gRNA_dist_count.py <plasmid_seq_file> <gRNA_library_file> <output_file_name.csv>")

# Create dictionary from gRNAs
gRNA_lib = get_gRNA_lib(gRNA_lib_file)
# Count each gRNA
gRNA_counts = count_gRNA(gRNA_lib, plasmid_seq_data)

# Export to file
with open(csv_out_file, 'w') as f:
    for key in gRNA_counts.keys():
        f.write("%s,%s\n"%(key,gRNA_counts[key]))