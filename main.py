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

def hamming_output(in_hamming, minterm, blength):
    for idx in xrange(len(in_hamming)):
        print "["+str(idx)+"]"
        for m in in_hamming[idx]:
            print "\t "+format_ba(minterm[m['idx'][0]], m['mask'], blength)
        print "------------------------------------------------------------"
        

def compute_stage(in_hamming, minterm):
    out_hamming = [[] for i in xrange(len(in_hamming))]
    primes = []
    for idx in xrange(len(in_hamming)-1):
        # print "idx",idx
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

                    nidx = s['idx']+d['idx']
                    nidx.sort()
                    ndic = {"idx": nidx, "mask":dc, "flag":0}
                    s['flag'] += 1
                    d['flag'] += 1
                    if ndic not in out_hamming[cb]:
                        out_hamming[cb].append(ndic)
    for idx in xrange(len(in_hamming)):
        for s in in_hamming[idx]:
            if s['flag'] == 0:
                primes.append(s)
    hamming_output(out_hamming, minterm, blength)
    print "Primes"
    for p in primes:
        print "\t "+format_ba(minterm[p['idx'][0]], p['mask'], blength)
    print "------------------------------------------------------------"
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
print "# of 1  :"+str(len(minterm_true))
print "# of DC :"+str(len(minterm_dc))

print "======================== Input Data ========================"

primes = []
hamming = [[] for i in range(blength+1)]

# Generating dictionary for each input
# 
for n in xrange(len(minterm)):
    hamming[count_bit(minterm[n])].append({"idx":[n], "mask":0, "flag":0})

for h in xrange(len(hamming)):
    print "["+str(h)+"]"
    for m in hamming[h]:
        print "\t m"+str(minterm[m['idx'][0]])+"\t "+\
            format(minterm[m['idx'][0]], '0'+str(blength)+'b')
    print "------------------------------------------------------------"

print "=============== Mid Process Computation ===================="
pp = 0
while 1:
    hamming_tuple = compute_stage(hamming, minterm)
    hamming2      = hamming_tuple[0]
    primes.extend(hamming_tuple[1])
    n = 0
    for h in hamming2:
        n += len(h)
    print "remains : ", n
    hamming = hamming2
    pp += 1

    if n == 0:
        break

print "======================== Prime Expression ========================"
i = 0
for p in primes:
    print format(i,'4d'),\
        format_ba(minterm[p['idx'][0]], p['mask'], blength)
    i+=1


# Generating btable;    
print "======================== Tables ========================"
btable = [0 for i in xrange(len(minterm_true))]

for j in xrange(len(minterm_true)):
    btable[j] = 0
    for p in primes:
        btable[j] *= 2
        if j in p['idx']:
            btable[j] += 1

for bt in xrange(len(btable)):
    print "p["+str(bt)+"]\t"+str(minterm_true[bt])+"\t"+\
        format(btable[bt],'0'+str(len(primes))+'b')
            
# 被覆問題
g = [btable[0]]
r = []

for bt in xrange(len(btable)):
    r = []
    for ib in reversed(xrange(len(primes))):
        eidx = btable[bt]&(2**ib)
        if eidx:
            for p in g:
                r.append(p|eidx)
    g = list(set(r))


g.sort(cmp=lambda x,y: cmp(count_bit(x), count_bit(y)))

print "========= Prime Combinations for Logical equations ========="
for i in g:
    print format(i,'0'+str(len(primes))+'b')+","+str(count_bit(i))

# 解の生成
result = []
for r in g:
    rtmp = []
    for ib in reversed(xrange(len(primes))):
        if r & (2**ib) :
            pidx = len(primes)-ib-1
            rtmp.append(format_ba(primes[pidx]['idx'][0],
                                  primes[pidx]['mask'],
                                  blength))
    result.append(rtmp)

min_num   = sys.maxint
min_array = []
for r in result:    
    if min_num > len(r):
        min_num = len(r)
        min_array = r

print "Number of Equations : "+str(min_num)
for i in r:
    print i,
print ""


#if __name__ == '__main__':



