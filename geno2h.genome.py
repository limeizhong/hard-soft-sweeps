import numpy as np 
from collections import Counter
import sys
import re

wg = [0,10,13,16,18,24,32,34,45,50,52,56,59,64,69,74,78,84,89,95,101,103,106,110,113,123,128,130,133,137,140,144,151,160,163,165,170,173,192,194,197,200,204,215,225,228,232,235,241,243,247,251,256,258,262,265,271,278,282,284,287,297]
lg = [1,3,4,6,7,8,9,15,20,21,22,26,27,30,31,35,39,40,41,43,47,48,49,51,53,54,58,61,65,66,67,71,72,76,77,81,82,86,87,88,90,94,96,97,98,100,108,109,111,114,115,121,124,125,126,127,129,136,139,141,142,145,149,150,152,153,154,155,156,157,159,162,164,167,168,171,172,174,175,177,178,179,182,184,186,187,190,191,193,196,202,206,207,212,213,216,217,218,219,220,222,224,226,227,230,231,233,240,242,244,245,249,250,252,254,255,257,259,260,264,275,283,286,288,290,292,293,298,300,301]
ig = [2,5,11,12,14,17,19,23,25,28,29,33,36,37,38,42,44,46,55,57,60,62,63,68,70,73,75,79,80,83,85,91,92,93,99,102,104,105,107,112,116,117,118,119,120,122,131,132,134,135,138,143,146,147,148,158,161,166,169,176,180,181,183,185,188,189,195,198,199,201,203,205,208,209,210,211,214,221,223,229,234,236,237,238,239,246,248,253,261,263,266,267,268,269,270,272,273,274,276,277,279,280,281,285,289,291,294,295,296,299]

def read_imputed_012_file(file):
    f= open(file,"rt")
    geno=[]
    for line in f.readlines():
        l=line.split()[2:]
        geno.append(l)
    geno=np.asarray(geno,dtype=str)
    return geno

def h12(freq_list):
    freq = np.asarray(freq_list)/float(sum(freq_list))
    freq_square = [freq[i]**2 for i in range(len(freq))]
    H1 = sum(freq_square)
    H2 = sum(freq_square[1:])
    H3 = sum(freq_square[2:])
    H12 = sum(freq[:2])**2+H3
    return H12,H2/H1

def hapkeys_freq(data,keys,group):
    hap=Counter([data[i] for i in group])
    values =sorted([v for v in hap.values()],reverse = True)
    h_statistics=h12(values)
    for key in keys:
        if key not in hap.keys():
            hap[key]=0
    return hap,h_statistics

def main_h12(file1,file2,win=50,step=20):
    geno=read_imputed_012_file(file1)
    f1_out=open(file2,"wt")
    n=(len(geno)-int(win))/int(step)
    if int(n)<n:
        n=n+1
    for i in range(int(n)):
        st=int(step)*i
        ed=int(win)+int(step)*i
        li = list(range(st,ed))
        current_hap=[]
        for ind in range(302):
            current_hap.append("".join(geno[li,ind]))
        all_hap = Counter(current_hap)
        vk = sorted([(v,k) for k,v in all_hap.items()],reverse = True)
        keys =[e[1] for e in vk]
        whap,wh=hapkeys_freq(current_hap,keys,wg)
        lhap,lh=hapkeys_freq(current_hap,keys,lg) 
        ihap,ih=hapkeys_freq(current_hap,keys,ig) 
        print(st,ed,wh[0],wh[1],lh[0],lh[1],ih[0],ih[1],file=f1_out)

if __name__ == '__main__': 
    main_h12(sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4])

# python geno2h.genome.py Gm01_012_geno Gm01.h12.50.stats 50 20
