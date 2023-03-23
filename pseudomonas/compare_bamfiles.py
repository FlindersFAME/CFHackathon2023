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
    parser.add_argument('-v', '--verbose', help='verbose output', action='store_true')
    args = parser.parse_args()

    dontwant = set()
    for bamfile in args.notin:
        if args.verbose:
            print(f"Excluding reads from {bamfile}", file=sys.stderr)
        bampysam = pysam.AlignmentFile(bamfile, "rb")
        for read in bampysam.fetch(until_eof=True):
            dontwant.add(read.query_name)

    bampysam = pysam.AlignmentFile(args.b, "rb")
    for read in bampysam.fetch(until_eof=True):
        if read.query_name not in dontwant:
            print(read.query_name)
