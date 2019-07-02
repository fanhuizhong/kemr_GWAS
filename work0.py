#!/usr/bin/env python

import re,os,time
import pandas as pd
import numpy as np
from Bio import SeqIO
  
def kmer_get(aastr,k):
    rst = []
    for i in xrange(0,len(aastr)-k+1):
        getStr = str(aastr[i:i+k])
        if not re.search('Z|X|\*|U',getStr):
            rst.append(KmerOpt.aa2int(aastr[i:i+k]))
    return rst
    
    
def read_data(dir):
    m = [];maxN=1280000000
    df = pd.DataFrame(np.zeros((maxN,102)),index = xrange(maxN))
    c = 0
    for i in os.listdir(dir):
        if not re.search('.cv.txt',i):continue
        m.append(i.split('.')[0])
        a = pd.read_table(dir + i,sep=' ',skiprows=3,header=None)
        df.iloc[a.iloc[:,0].tolist(),c] = a.iloc[:,1].tolist()
        c += 1
    df.columns = m
    return df
    
def gene_kmer(df,ref_fa,k):
    gene = {}
    for i in SeqIO.parse(ref_fa,'fasta'):
        m = re.sub('-','',str(i.seq))
        if gene.has_key(i.id):
            gene [i.id] += kmer_get(m,k)
        else:
            gene [i.id] = kmer_get(m,k)

    gene_all = {}
    for i,j in gene.items():
        gene_all[i]= df.iloc[np.unique(np.array(j)),:]# bin frequency
    for i,j in gene_all.items():
        j.to_csv('/home/wuq/projects/kmerEvol/data_version3/work0_result/' + i, encoding='utf-8', index=True)
    
if __name__ == "__main__":
    
    print time.ctime()
    df = read_data('/home/wuq/projects/kmerEvol/data_version3/frequency/')# The folder contains kmer frequency file for each species, this frequency files were produced using CVtree software 
    print time.ctime()
    gene_kmer(df,'/home/wuq/projects/kmerEvol/data_version3/Homo_sapiens.GRCh38.pep.all.fa',7)
    print time.ctime()

