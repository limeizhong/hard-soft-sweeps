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
      
def gene_index_file(file):
    f= open(file,"rt")
    gene_sub={}
    for line in f.readlines():
        l=line.split()
        gene_sub[l[0]]=[int(l[2]),int(l[3])]
    return gene_sub

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

def main_h12(genofile,geneindexfile,out,step=0):
    geno=read_imputed_012_file(genofile)
    gene_sub=gene_index_file(geneindexfile)
    f_out=open(out,"wt")
    for key in gene_sub.keys():
        gene_st=max(gene_sub[key][0]-int(step),0)
        gene_ed=min(gene_sub[key][1]+int(step),len(geno))
        if gene_ed != gene_st:
            li = list(range(gene_st,gene_ed))
            current_hap=[]
            for ind in range(302):
                current_hap.append("".join(geno[li,ind]))
            all_hap = Counter(current_hap)
            vk = sorted([(v,k) for k,v in all_hap.items()],reverse = True)
            keys =[e[1] for e in vk]
            whap,wh=hapkeys_freq(current_hap,keys,wg)
            lhap,lh=hapkeys_freq(current_hap,keys,lg) 
            ihap,ih=hapkeys_freq(current_hap,keys,ig) 
            lvk = sorted([(v,k) for k,v in lhap.items()],reverse = True)
            lk =[e[1] for e in lvk]
            k1=lk[0]
            k2=lk[1]
            sim_l12=len([i for i in range(len(k1)) if k1[i]==k2[i]])/len(k1)
            sim_wl1=[]
            sim_wl2=[]
            whap0=Counter([current_hap[i] for i in wg])
            for wk in whap0.keys():
                assert len(wk)==len(k1)
                sim=len([i for i in range(len(k1)) if wk[i]==k1[i]])/len(k1)
                sim_wl1.append(sim)
                assert len(wk)==len(k2)
                sim=len([i for i in range(len(k2)) if wk[i]==k2[i]])/len(k2)
                sim_wl2.append(sim)
            wl1=max(sim_wl1)
            wl2=max(sim_wl2)
            res=min([wl1,wl2])
            if res >= sim_l12:
                res1="wl"
            else:
                res1="ll"
            posl1=sim_wl1.index(wl1)
            posl2=sim_wl2.index(wl2)
            if posl1==posl2:
                res2="same"
            else:
                res2="different"
            print(key,1,sim_l12,wl1,wl2,res1,res2,whap[k1],lhap[k1],ihap[k1],file=f_out)
            print(key,2,sim_l12,wl1,wl2,res1,res2,whap[k2],lhap[k2],ihap[k2],file=f_out)
              
if __name__ == '__main__': 
    main_h12(sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4]) 

# python hap.similarity.py Gm01_012_geno Gm01_gene.index Gm01.h12.sim 100