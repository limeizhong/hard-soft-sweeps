import re
import numpy as np 
from collections import Counter
import sys

def h12(freq_list):
    freq = np.asarray(freq_list)/float(sum(freq_list))
    freq_square = [freq[i]**2 for i in range(len(freq))]
    H1 = sum(freq_square)
    H2 = sum(freq_square[1:])
    H3 = sum(freq_square[2:])
    H12 = sum(freq[:2])**2+H3
    return H12,H2/H1

f=open(sys.argv[1],"rt")
f_out=open(sys.argv[1]+'.h12',"wt")

all_ms=[]
temp=[]
for line in f.readlines():
    l=line.strip()
    if l.startswith('0'):
        temp.append(l)
    if l.startswith('1'):
        temp.append(l)
    if l.startswith('positions'):
        if temp!=[]:
            all_ms.append(temp)
        temp=[]
        
all_ms.append(temp)
print (len(all_ms))

for i in range(len(all_ms)):
    hap=Counter(all_ms[i])
    values =sorted([v for v in hap.values()],reverse = True)
    h_statistics=h12(values)
    print(i,h_statistics[0],h_statistics[1],file=f_out)

