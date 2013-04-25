#!/usr/bin/env python
# -*- coding:utf-8 -*-

import argparse
import sys
import qmpy

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('-o', nargs='*')
    parser.add_argument('-dc', nargs='*')
    parser.add_argument('-v', '--version',
                        action='version',
                        version='%(prog)s 0.0.1') # version
    args = parser.parse_args()
    
    minterm_true = map((lambda x: int(x)), args.o) if args.o else []
    minterm_dc = map((lambda x: int(x)), args.dc) if args.dc else []

    result = qmpy.QM(minterm_true, minterm_dc)

    min_num   = sys.maxint
    min_array = []
    for r in result:    
        if min_num > len(r):
            min_num = len(r)
            min_array = r

    print "Number of Equations : "+str(min_num)
    for i in min_array:
        print i,
    print ""

    qmpy.QM_validation(min_array, minterm_true, minterm_dc)

