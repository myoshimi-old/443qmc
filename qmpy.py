#!/usr/bin/env python
# -*- coding:utf-8 -*-

import argparse
import sys
from functools import reduce

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
        a = int(a / 2)
        m = int(m / 2)
        d -= 1
    return x

def hamming_output(in_hamming, minterm, blength):
    for idx in range(len(in_hamming)):
        print(("["+str(idx)+"]"))
        for m in in_hamming[idx]:
            print(("\t "+(format_ba(minterm[m['idx'][0]], m['mask'], blength))))
        print("------------------------------------------------------------")
        

def compute_stage(in_hamming, minterm, blength):
    out_hamming = [[] for i in range(len(in_hamming))]
    primes = []
    for idx in range(len(in_hamming)-1):
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
    for idx in range(len(in_hamming)):
        for s in in_hamming[idx]:
            if s['flag'] == 0:
                primes.append(s)

    hamming_output(out_hamming, minterm, blength)
    print("Primes")
    for p in primes:
        print(("\t "+format_ba(minterm[p['idx'][0]], p['mask'], blength)))
    print("------------------------------------------------------------")
    return (out_hamming, primes)

def QM(minterm_true, minterm_dc):
    print("======================== Computing =========================")
    minterm = minterm_true + minterm_dc
    blength= len(format(max(minterm),'0b')) # ビット長の計算
    print(("blength : "+str(blength)))
    print(("# of 1  : "+str(len(minterm_true))))
    print(("# of DC : "+str(len(minterm_dc))))

    print("======================== Input Data ========================")

    primes = []
    hamming = [[] for i in range(blength+1)]

    # Generating dictionary for each input
    # 
    for n in range(len(minterm)):
        hamming[count_bit(minterm[n])].append({"idx":[n], "mask":0, "flag":0})

    for h in range(len(hamming)):
        print(("["+str(h)+"]"))
        for m in hamming[h]:
            print(("\t m"+str(minterm[m['idx'][0]])+"\t "+\
                       format(minterm[m['idx'][0]], '0'+str(blength)+'b')))
        print("------------------------------------------------------------")

    print("=============== Mid Process Computation ====================")
    n = sys.maxsize
    while n != 0:
        hamming_tuple = compute_stage(hamming, minterm, blength)
        hamming2      = hamming_tuple[0]
        primes.extend(hamming_tuple[1])
        hamming = hamming2
        n = reduce(lambda a,b: a+b, list(map((lambda x: len(x)), hamming)))
        print(("remains : "+str(n)))

    print("======================== Prime Expression ========================")
    i = 0
    for p in primes:
        print(format(i,'4d')+" "+\
                  format_ba(minterm[p['idx'][0]], p['mask'], blength))
        i+=1

    # Generating btable;    
    print("======================== Tables ========================")
    btable = [0 for i in range(len(minterm_true))]

    for j in range(len(minterm_true)):
        btable[j] = 0
        for p in primes:
            btable[j] *= 2
            if j in p['idx']:
                btable[j] += 1

    for bt in range(len(btable)):
        print(("p["+str(bt)+"]\t"+str(minterm_true[bt])+"\t"+\
                   format(btable[bt],'0'+str(len(primes))+'b')))
            
    # 被覆問題
    g = [btable[0]]
    r = []

    for bt in range(len(btable)):
        r = []
        for ib in reversed(list(range(len(primes)))):
            eidx = btable[bt]&(2**ib)
            if eidx:
                for p in g:
                    r.append(p|eidx)
        g = list(set(r))

    #g.sort(cmp=lambda x,y: cmp(count_bit(x), count_bit(y)))
    g.sort(key=lambda x: count_bit(x))

    print("========= Prime Combinations for Logical equations =========")
    for i in g:
        print((format(i,'0'+str(len(primes))+'b')+","+str(count_bit(i))))

    # 解の生成
    result = []
    for r in g:
        nr = format(r, '0'+str(len(primes))+'b')
        rtmp = []
        for j in range(len(nr)):
            if nr[j] == '1':
                rtmp.append(format_ba(minterm[primes[j]['idx'][0]],
                                      primes[j]['mask'],
                                      blength))
        #print rtmp
        result.append(rtmp)
    
    print("============== Computation Finished ========================")

    return result

def QM_validation(logeq, minterm_true, minterm_dc):
    minterm = minterm_true + minterm_dc
    blength= len(format(max(minterm),'0b')) # ビット長の計算

    valeq = []
    
    for eq in logeq:
        t = int(''.join(map((lambda x: '0' if x == 'X' else x), eq)), 2)
        d = int(''.join(map((lambda x: '1' if x == 'X' else '0'), eq)), 2)
        valeq.append((t, d))

    print("=================== Validation =============================")
    valid = []
    for m in range(2**blength):
        for v in valeq:
            if (m | v[1]) ^ (v[0] | v[1]) == 0:
                valid.append(m)
    valid = sorted(list(set(valid)))

    mm = minterm_true
    mm.sort()
    s = ""
    for m in mm:
        s+=str(m)+" "
    print(s)
    
    s = ""
    for v in valid:
        s+=str(v)+" "
    print(s)
    
    f = 0
    if len(mm) == len(valid):
        for x in range(len(mm)):
            if mm[x] != valid[x]:
                f += 1
    if (f == 0) and (len(mm) == len(valid)):
        print("True")
    else:
        print("Failed")


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

    result = QM(minterm_true, minterm_dc)

    min_num   = sys.maxsize
    min_array = []
    for r in result:    
        if min_num > len(r):
            min_num = len(r)
            min_array = r

    print(("Number of Equations : "+str(min_num)))
    s = ""
    for i in min_array:
        s+=str(i)+" "
    print(s)

    QM_validation(min_array, minterm_true, minterm_dc)

