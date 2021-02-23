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

def combine_x2( x):
    return float(x * (x-1) / 2)

def compute_homozygosity(haps):
    haps=np.asarray(haps)
    uniques = np.unique(haps) # get unique keys
    lengths = [len(np.where(haps == u)[0]) for u in uniques]
    Len = len(haps)
    return sum(combine_x2(l) for l in lengths) / combine_x2(Len)

def ehh_keys(haplotype):  
    ehh_key = []
    uniques = np.unique(haplotype)
    for u in uniques:
        kp = np.where(haplotype == u)[0]        
        kp_w = [i for i in kp if i in wg]
        kp_l = [i for i in kp if i in lg]
        kp_i = [i for i in kp if i in ig]
        if len(kp_w) >= 4 or len(kp_l) >= 7 or len(kp_i) >= 6:
            ehh_key.append([kp_w,kp_l,kp_i,u])       
    return ehh_key

def main_ehh(file1,file2,file3,stop=0.05,step=20):
    geno=read_imputed_012_file(file1)
    gene_sub=gene_index_file(file2)
    f_out=open(file3,"wt")
    for key in gene_sub.keys():
        gene_st=gene_sub[key][0]
        gene_ed=gene_sub[key][1]
        if gene_ed != gene_st:
            li = list(range(gene_st,gene_ed))
            current_hap=[]
            for ind in range(302):
                current_hap.append("".join(geno[li,ind]))
            current_hap=np.asarray(current_hap)
            ehh_key=ehh_keys(current_hap)

            for k in ehh_key:
                print(key,len(k[0]),len(k[1]),len(k[2]),end = ' ',file=f_out)
                total_res=[]
                for g in [k[0],k[1],k[2]]:                    
                    total_up=0
                    total_dw=0
                    hom_up=float(stop)
                    hom_dw=float(stop)
                    if len(g)>3:
                        for h in range(1,5000):
                            if hom_up >= float(stop):
                                li_up = list(range(max(gene_st-int(step)*h,0),gene_ed))
                                current=[]
                                for indv in g:
                                    current.append("".join(geno[li_up,indv]))
                                hom_up=compute_homozygosity(current)
                                total_up+=hom_up*int(step)
                            if max(gene_st-int(step)*h,0)==0:
                                break
                            if hom_up < float(stop):
                                break
                        for h in range(1,5000):
                            if hom_dw >= float(stop):
                                li_dw = list(range(gene_st,min(gene_ed+int(step)*h,len(geno))))
                                current=[]
                                for indv in g:
                                    current.append("".join(geno[li_dw,indv]))
                                hom_dw=compute_homozygosity(current)
                                total_dw+=hom_dw*int(step)
                            if min(gene_ed+int(step)*h,len(geno))==len(geno):
                                break
                            if hom_dw < float(stop):
                                break                            
                    total_res.append(total_up)
                    total_res.append(total_dw)
                for res in total_res:
                    print(res,end = ' ',file=f_out)
                print(k[3],file=f_out)  
                
if __name__ == '__main__': 
    main_ehh(sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4],sys.argv[5])

# python geno2ihh.py Gm01_012_geno Gm01_gene.index Gm01.gene.hap.ehh 0.05 20