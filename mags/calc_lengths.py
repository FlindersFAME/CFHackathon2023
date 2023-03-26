#!/usr/bin/env python3
import argparse
from argparse import RawTextHelpFormatter
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import pandas as pd
from Bio.Seq import Seq

def get_input():
    usage = 'python3 calc_lengths.py ...'
    parser = argparse.ArgumentParser(description='script to calculate lengths from multiFASTA files ', formatter_class=RawTextHelpFormatter)
    parser.add_argument('-i', '--infile', action="store", help='input FASTA file',  required=True)
    parser.add_argument('-o', '--outfile', action="store", help='outfile file in csv format',  required=True)
    args = parser.parse_args()

    return args

args = get_input()


contigs = []
lengths = []


# write the fasta

for dna_record in SeqIO.parse(args.infile, "fasta"):

    l = len(dna_record.seq)
    contigs.append(dna_record.id)
    lengths.append(l)

# contig and length df

length_df = pd.DataFrame(
{'contig': contigs,
    'length': lengths
})

# save the length_gc.tsv also
length_df.to_csv(args.outfile, sep=",", index=False)

