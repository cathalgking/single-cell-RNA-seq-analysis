#!/usr/bin/env python

import sys
for line in sys.stdin:
     seq = line.strip()
     seq = f"{seq[:10]}{seq[22:31]}{seq[44:]}"
     print(seq)
