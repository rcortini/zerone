#!/usr/bin/env python
#-*- coding:utf-8 -*-

import sys

for line in sys.stdin:
    if len(line.split()[4].split(':')[0]) <= 5:
        print line.rstrip()
