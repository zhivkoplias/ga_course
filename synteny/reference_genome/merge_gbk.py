#!/usr/bin/env python
# Copyright 2016 Peter Cock, James Hutton Institute
# This example Biopython script is released for free re-use under
# the MIT license https://opensource.org/licenses/MIT
# For backgound see http://seqanswers.com/forums/showthread.php?t=28924
import sys
from Bio.Alphabet import generic_dna
from Bio import SeqIO

if len(sys.argv) < 2:
    print("Usage: Takes as arguments filenames for input DNA GenBank files,")
    print("outputs merged record to stdout using 50bp runs of N as linker.")
    print("")
    print("python merge_gbk.py input*.gbk > combined.gbk")

# argv[0] = python script
# argv[1] = first record
# argv[2] onwards = subsequence records
filename = sys.argv[1]
sys.stderr.write("Loading %s\n" % filename)
combined_record = SeqIO.read(filename, "genbank")
# Forcing DNA alphabet in case of improper input files
combined_record.seq.alphabet = generic_dna
for filename in sys.argv[2:]:
    sys.stderr.write("Loading %s\n" % filename)
    record = SeqIO.read(filename, "genbank")
    record.seq.alphabet = generic_dna
    combined_record = combined_record + ("N" * 50) + record
sys.stderr.write("Merging into one record\n")
count = SeqIO.write(combined_record, sys.stdout, "genbank")
assert count == 1, count
