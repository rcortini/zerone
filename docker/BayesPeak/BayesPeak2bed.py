#!/usr/bin/env python
#-*- coding:utf-8 -*-

import sys

with open(sys.argv[1]) as fin:
    for line in fin:
        chrom, start, end = line.split()
        start, end = int(start) + 1, int(end) + 1
        print '\t'.join([chrom, str(start), str(end)])
