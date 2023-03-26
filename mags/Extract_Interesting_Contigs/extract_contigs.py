#!/usr/bin/env python3
import argparse
from argparse import RawTextHelpFormatter
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import pandas as pd
from Bio.Seq import Seq

def get_input():
    usage = 'python3 extract_contigs.py ...'
    parser = argparse.ArgumentParser(description='script to get FASTA sequences for certain contigs input in a csv ', formatter_class=RawTextHelpFormatter)
    parser.add_argument('-i', '--infile', action="store", help='input FASTA file',  required=True)
    parser.add_argument('-c', '--contig', action="store", help='input file of desired contigs from gene_data.csv in csv (with header row)',  required=True)
    parser.add_argument('-o', '--outfile', action="store", help='outfile file in fasta format',  required=True)
    args = parser.parse_args()

    return args

args = get_input()

# read the csv in 

cols = ['contigs']

try:
    contigs_df = pd.read_csv(args.contig, delimiter= ',\t', index_col=False, names = cols, skiprows=1)
except pd.errors.EmptyDataError:
    print('csv is empty')

# write the fasta
with open(args.outfile, 'w') as out_fasta:
    for dna_record in SeqIO.parse(args.infile, "fasta"):
        # if contig is in contigs save it
        # https://stackoverflow.com/questions/30944577/check-if-string-is-in-a-pandas-dataframe
        if contigs_df['contigs'].str.contains(dna_record.id ).any():
            SeqIO.write(dna_record, out_fasta, 'fasta')


