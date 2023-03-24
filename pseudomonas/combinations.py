"""
Print the possible combinations of files
"""

import os
import sys
import argparse
from itertools import permutations

__author__ = 'Rob Edwards'

file_endings = ['P_aeruginosa.bam', '16S.bam', 'Myco100.bam', 'Staph.bam', 'Strep.bam', 'Pseudomonas_aeruginosa.fasta']

with open("IDs.tsv", "r") as ids:
    for sample in ids:
        joined = [f"{s}.{x}" for x in file_endings]
        for i in range(len(joined)):
            p = joined.pop(0)
            t = " -n ".join(joined)
            if "16S" not in p:
                print(f"python ~/GitHubs/CFHackathon2023/pseudomonas/compare_bamfiles.py -p -b {p} -n {t} ")
            joined.append(p)
