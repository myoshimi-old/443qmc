#!/usr/bin/env python
# -*- coding:utf-8 -*-

import argparse
import sys

def count_bit(x):
    i = 0
    while(x):
        x &= x-1
        i=i+1
    return i

def hamming_distance_mask(a, b, ma, mb):
    d = (a|ma|mb)^(b|ma|mb)
    dm = ma ^ mb
    return count_bit(d) + count_bit(dm)

def format_ba(a, m, d):
    x = ""
    while(d):
        if m%2:
            x = 'X'+x
        else:
            x = str(a%2)+x
        a /= 2
        m /= 2
        d -= 1
    return x

def compute_stage(in_hamming, minterm):
    out_hamming = [[] for i in xrange(len(in_hamming))]
    primes = []
    for idx in xrange(len(in_hamming)-1):
        print "idx",idx
        for s in in_hamming[idx]:
            for d in in_hamming[idx+1]:
                # hd: ハミング距離
                hd = hamming_distance_mask(minterm[s['idx'][0]],
                                           minterm[d['idx'][0]],
                                           s['mask'], d['mask'])
                if hd == 1:
                    # dc : dont care位置
                    dc = ((minterm[s['idx'][0]]|s['mask']) \
                              ^ (minterm[d['idx'][0]]|d['mask'])) \
                              | s['mask'] | d['mask']
                    # cb : dont careを除くビットの数
                    cb = count_bit(minterm[s['idx'][0]]|dc) - count_bit(dc)
                    
                    print " s",format(minterm[s['idx'][0]],'0'+str(blength)+'b'),
                    print " d",format(minterm[d['idx'][0]],'0'+str(blength)+'b'),
                    print " hamming:",hd,
                    print " dc",format(dc,'0'+str(blength)+'b'),"cb:",cb

                    nidx = s['idx']+d['idx']
                    nidx.sort()
                    ndic = {"idx": nidx, "mask":dc, "flag":0}
                    s['flag'] += 1
                    d['flag'] += 1
                    if ndic not in out_hamming[cb]:
                        print "append",ndic
                        out_hamming[cb].append(ndic)
    for idx in range(len(in_hamming)):
        for s in in_hamming[idx]:
            if s['flag'] == 0:
                print "result",
                primes.append(s)
                print s, "*"
    return (out_hamming, primes)



parser = argparse.ArgumentParser(description='')
parser.add_argument('-o', nargs='*')
parser.add_argument('-dc', nargs='*')
parser.add_argument('-v', '--version',
                    action='version',
                    version='%(prog)s 0.0.1') # version
args = parser.parse_args()

print args

minterm_true = []
minterm_dc   = []

if args.o:
    for n in args.o:
        minterm_true.append(int(n))

if args.dc:
    for n in args.dc:
        minterm_dc.append(int(n))

print "======================== Computing ========================="
minterm = minterm_true + minterm_dc
blength= len(format(max(minterm),'0b')) # ビット長の計算
print "blength:", blength

"""
i = 0
for n in minterm_true:
    print format(i, '4d')+" : "+ \
        format(n, '0'+str(blength)+'b') + " : " + str(n)
    i = i+1
for n in minterm_dc:
    print format(i,'4d')+" : "+ \
        format(n, '0'+str(blength)+'b') + " : " + str(n)
    i = i+1
"""
print "======================== Input Data ========================"

primes = []
hamming = [[] for i in range(blength+1)]

i=0
for n in minterm:
    # bc : ビットカウント
    bc = count_bit(n)
    hamming[bc].append({"idx":[i], "mask":0, "flag":0})
    i=i+1

n=0
for h in hamming:
    print "["+str(n)+"]"
    for m in h:
        print "\t m"+str(minterm[m['idx'][0]])+"\t "+\
            format(minterm[m['idx'][0]], '0'+str(blength)+'b')
    print "------------------------------------------------------------"
    n += 1


print "=============== Mid Process Computation ===================="
pp = 0
while 1:
    hamming_tuple = compute_stage(hamming, minterm)
    hamming2      = hamming_tuple[0]
    primes.extend(hamming_tuple[1])
    n = 0
    for h in hamming2:
        print h
        n += len(h)
        i=i+1
    print "remains : ", n
    hamming = hamming2
    pp += 1
    if n == 0:
        break

print "======================== Prime Expression ========================"
i = 0
for p in primes:
    print format(i,'4d'),\
        format_ba(minterm[p['idx'][0]], p['mask'], blength), \
        format(minterm[p['idx'][0]], '5d'), format(p['mask'], '08b')
    i+=1

btable = [0 for i in xrange(len(minterm_true))]
for j in range(len(minterm_true)):
    btable[j] = 0
    for p in primes:
        btable[j] *= 2
        if j in p['idx']:
            btable[j] += 1

# 被覆問題
g = [btable[0]]
r = []

for bt in btable:
    print "bt : "+format(bt,'0'+str(len(primes))+'b')
    r = []
    for ib in reversed(xrange(len(primes))):
        eidx = bt&(2**ib)
        if eidx:
            #print "  "+format(eidx,'0'+str(blength)+'b')
            for p in g:
                r.append(p|eidx)

    g = []
    for p in r:
        if p not in g:
            g.append(p)
    """
    print " r[",
    for i in r:
        print format(i,'0'+str(blength)+'b')+",",
    print "]"
    print " g[",
    for i in g:
        print format(i,'0'+str(blength)+'b')+",",
    print "]"
    """
 
print " g[",
for i in g:
    print format(i,'0'+str(len(primes))+'b')+","+str(count_bit(i))
print "]"

# 解の生成
result = []
for r in g:
    rtmp = []
    for ib in reversed(xrange(len(primes))):
        if r & 2**ib :
            pidx = len(primes)-1-ib
            midx = primes[pidx]['idx'][0]
            mask = primes[pidx]['mask']
            rtmp.append(format_ba(minterm[midx],mask, blength))
    result.append(rtmp)

min_num   = sys.maxint
min_array = []
for r in result:    
    if min_num > len(r):
        min_num = len(r)
        min_array = r

print "Number of Equations : "+str(len(r))
print r

