#!/usr/bin/python

import argparse
import sys
import re

# define args
parser = argparse.ArgumentParser(description='Filter bwa NM tag <= n from a bam file')
parser.add_argument('-n', '--nm', type = int, help='max edit distance allowed')

# parse args
args = parser.parse_args()
n = int(args.nm)
for line in sys.stdin:
    NMsearch = re.search('NM:i:([0-9]+)', line)
    if NMsearch:
        NM = NMsearch.group(1)
        if int(NM) <= n:
            sys.stdout.write(line)
