import numpy as np
import sys

def process(pos012_file,gene_file,output_file):
    f1=open(pos012_file,'rt')
    f2=open(gene_file,'rt')
    f3=open(output_file,'wt')
    chr_pos = {}
    chrom="Gm00"
    pos=[]
    for line in f1.readlines():
        l = line.split()     
        if l[0]==chrom: 
            pos.append(int(l[1]))
        if l[0]!=chrom:
            chr_pos[chrom]=pos
            chrom=l[0]
            pos=[]
            pos.append(int(l[1]))
    chr_pos[chrom]=pos # for "Gm20"

    for line in f2.readlines()[1:]:
        l = line.split()
        gene_start=max(int(l[2])-1000,min(chr_pos[l[1]]))
        gene_end=min(int(l[3])+1000,max(chr_pos[l[1]])) 
        gene_len=gene_end-gene_start+1                  
        add1 = sorted([gene_start]+chr_pos[l[1]])
        start = add1.index(gene_start)
        add2 = sorted([gene_end]+chr_pos[l[1]])
        end = add2.index(gene_end)    
        print(l[0],l[1],start,end,file=f3)

if __name__ == '__main__': 
    process(sys.argv[1],sys.argv[2],sys.argv[3])