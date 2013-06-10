#!/usr/bin/env python
# -*- coding:utf-8 -*-

import sys
import qmpy
import random

# ========================================================================

num = 400
array = []
minvalue=1
maxvalue=1000

while len(array) < num:
    r = random.randint(minvalue, maxvalue)
    if not r in array:
        array.append(r)

print sorted(array)
sorted(array)

# ========================================================================

result = qmpy.QM(array, [])

print("Number of Equations : "+str(len(result)))
for r in result:
    print(r)
print("")
    
qmpy.QM_validation(result, array, [])

