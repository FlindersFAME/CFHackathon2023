"""
Find all the reads that are in a bam file that are NOT in one or more other bam files.

We make a bam file with Pseudomonas hits, and then other bam files with 16S, 5S, 23S (well, in theory. At the moment
just 16S).
"""

import os
import sys
import argparse
import pysam

__author__ = 'Rob Edwards'

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='What unique reads are in my bam file?')
    parser.add_argument('-b', '--bam', help='bam file where all reads should be', required=True)
    parser.add_argument('-n', '--notin', help='bam file where reads should NOT be', required=True, action='append')
    parser.add_argument('-f', '--fasta', help='output reads to a fasta file')
    parser.add_argument('-r', '--reads', help='list the reads')
    grp = parser.add_mutually_exclusive_group(required=True)
    grp.add_argument('-m', '--mapped', help='Each read is mapped', action='store_true')
    grp.add_argument('-p', '--proper', help='read is mapped in a proper pair', action='store_true')
    parser.add_argument('-v', '--verbose', help='verbose output', action='store_true')
    args = parser.parse_args()

    dontwant = set()
    for bamfile in args.notin:
        if args.verbose:
            print(f"Excluding reads from {bamfile}", file=sys.stderr)
        bampysam = pysam.AlignmentFile(bamfile, "rb")
        for read in bampysam.fetch(until_eof=True):
            if args.mapped and read.is_mapped:
                dontwant.add(read.query_name)
            elif args.proper and read.is_proper_pair:
                dontwant.add(read.query_name)


    if args.verbose:
        print(f"Found {len(dontwant)} reads to exclude", file=sys.stderr)
        print(f"Checking {args.bam}", file=sys.stderr)
    rc = 0
    trc = 0
    bampysam = pysam.AlignmentFile(args.bam, "rb")
    fasta = None
    if args.fasta:
        fasta = open(args.fasta, 'w')
    for read in bampysam.fetch(until_eof=True):
        if (args.mapped and read.is_mapped) or (args.proper and read.is_proper_pair):
            trc += 1
            if read.query_name not in dontwant:
                rc += 1
                if args.reads:
                    print(read.query_name)
                if fasta:
                    print(f">{read.query_name}\n{read.query_sequence}", file=fasta)


    if fasta:
        fasta.close()

    if args.verbose:
        print(f"In {args.bam}, we kept a {rc} reads from {trc} reads in total", file=sys.stderr)
    print(f"Bamfile: {args.bam}\tReads: {rc}")
