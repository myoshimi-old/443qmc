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
    
    minterm_true = list(map((lambda x: int(x)), args.o)) if args.o else []
    minterm_dc = list(map((lambda x: int(x)), args.dc)) if args.dc else []

    result = qmpy.QM(minterm_true, minterm_dc)

    print("Number of Equations : "+str(len(result)))
    for r in result:
        print(r)
    print("")
    
    qmpy.QM_validation(result, minterm_true, minterm_dc)

