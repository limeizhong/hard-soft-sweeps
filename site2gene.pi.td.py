import numpy as np
import sys

def a1(n):
    a1 = 0.0
    for i in range(1,int(n)):
        a1 += 1.0/i
    return a1

def Theta(S,length,N):
    theta = S/(length* a1(N))*1000
    return theta 

def Pi(K,length):
    pi = K/length*1000
    return pi

def tajimaD(S, K, N):
    N=int(N)
    if S == 0:
        D = -10.0
    else:    	
        a2 = 0
        for i in range(1, N):
            a2 += 1.0/(i*i)
        b1 = (N+1)/(3*(N-1))
        c1 = b1-1/a1(N)	
        b2 = (2*(N**2+N+3))/(9*N*(N-1))
        c2 = b2-(N+2)/(a1(N)*N)+a2/(a1(N)**2)
        e1 = c1/a1(N)
        e2 = c2/(a1(N)**2+a2)
        D = (K-S/a1(N))/np.sqrt(e1*S+e2*S*(S-1))
    return D

def process(chr_pop,gene,pop_size):
    f1 = open(chr_pop+'.sites.pi.theta','rt')
    pos=[]
    div_array=[]
    for line in f1.readlines():
        l = line.split()
        pos.append(int(l[1]))
        div_array.append(l[2:])              
    div_array=np.asarray(div_array,dtype=float)

    f2=open(gene,'rt')
    f3=open(chr_pop+'.gene.div','wt')
    # f4=open(chr_pop+'.gene.index','wt')
    for line in f2.readlines():
        l = line.split()
        gene_start=max(int(l[2])-1000,min(pos))
        gene_end=min(int(l[3])+1000,max(pos)) 
        gene_len=gene_end-gene_start+1                  
        add1 = sorted([gene_start]+pos)
        start = add1.index(gene_start)
        add2 = sorted([gene_end]+pos)
        end = add2.index(gene_end)
        (theta,pi,td)=(0,0,0)
        if end-start==0:
            theta=0
            pi=0
            td=-10
        else:
            gene_index = list(range(start,end))
            K=sum([div_array[i][0] for i in gene_index]) 
            S=sum([div_array[i][1] for i in gene_index])
            pi=Pi(K,gene_len)
            theta =Theta(S,gene_len,pop_size)
            td=tajimaD(S,K,pop_size)
        print(l[0],l[1],K,S,pi,theta,td,file=f3)
        # print(l[0],l[1],start,end,file=f4)

if __name__ == '__main__': 
    process(sys.argv[1],sys.argv[2],sys.argv[2])
