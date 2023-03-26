"""
Parse mmseqs hits and report those that don't overlap the hit with the lowest e-value
"""

import os
import sys
import argparse
from hackathon_lib import stream_blast_results

__author__ = 'Rob Edwards'


def non_overlapping_mmseqs(m8file: str, mineval: float, evalue_check: bool = False, verbose: bool = False):
    """
    Find non-overlapping hits
    :param m8file: the mmseqs2 m8-format output file
    :param mineval: the minimum evalue
    :param evalue_check: check whether discarded hits have a better e value than stored hits
    :param verbose: more output
    """

    hits = []
    last_query =  None
    positions = {}
    for r in stream_blast_results(m8file, verbose):
        if r.evalue > mineval:
            continue
        if r.query != last_query:
            if verbose and last_query and not hits:
                print(f"There were no hits for {last_query}", file=sys.stderr)
            for h in hits:
                print(str(h))
            hits = []
            positions = {}
            last_query = r.query

        qsta = r.query_start
        qend = r.query_end
        if qsta > qend:
            (qend, qsta) = (qsta, qend)
        keep = True
        for i in range(qsta, qend + 1):
            if i in positions:
                keep = False
                if evalue_check:
                    if r.evalue < positions[i].evalue:
                        print(f"We had a hit for {positions[i].query} -> {positions[i].db} : {positions[i].evalue} and now",
                              f" {r.query} -> {r.db} : {r.evalue}", file=sys.stderr)
        if keep:
            for i in range(qsta, qend + 1):
                positions[i] = r
            hits.append(r)
    print(list(map(str, hits)))







if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Parse m8 format mmseqs results')
    parser.add_argument('-f', help='mmseqs m8 format output file', required=True)
    parser.add_argument('-e', help=f'minimum e value [Default: %(default)f]', default=1e-5, type=float)
    parser.add_argument('-c', help='Check whether better e values are discarded', action='store_true')
    parser.add_argument('-v', help='verbose output', action='store_true')
    args = parser.parse_args()

    non_overlapping_mmseqs(args.f, args.e, args.c, args.v)
