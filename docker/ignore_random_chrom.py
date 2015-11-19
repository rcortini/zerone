#!/usr/bin/env python
#-*- coding:utf-8 -*-

import sys

with open(sys.argv[1]) as f:
    for l in f:
        if len(l.split()[4].split(':')[0]) <= 5:
            print l.rstrip()
