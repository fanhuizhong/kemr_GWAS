#!/usr/bin/env python

import re,os,time
import pandas as pd
import numpy as np
from statsmodels.sandbox.stats import multicomp
import multiprocessing as mp


    
def filter(dir,gene):
    aa1 =  pd.read_csv(dir + gene,header = None,low_memory=False)
    a = aa1
    c = []
    for j in xrange(1,len(a)):
        m = a.iloc[j,:]
        aa = len(m[m=='0.0'])
        bb = len(m[m=='-1.0'])
        if float((aa+bb)/102.) > 0.1:
            c.append(j)
    b = a.drop(c,axis=0)
    return(b)

def yourFunc(dir,geneList,qRst):
    allRst = []
    for i in geneList:
        temp = filter(dir,i)
        if len(temp) < 6 :continue
        allRst.append( [i,temp] )
        temp.to_csv('/home/wuq/projects/kmerEvol/data_version3/work1_1_filter_result/' + i,encoding='utf-8',index = False,header=False)


if __name__ == "__main__":
    
    dir = '/home/wuq/projects/kmerEvol/data_version3/work0_result/'
    
    corN = 50
    aimLst =  os.listdir(dir)
    splitPosition = [i for i in range(0,len(aimLst),len(aimLst)/corN) ]
    plist = []
    manager = mp.Manager() 
    qRst = manager.Queue(maxsize=0)
    p = mp.Pool()
    for i,ii in enumerate(splitPosition):
        if i+1 < len(splitPosition):
            inputList = aimLst[ii:splitPosition[i+1]]
        else:
            inputList = aimLst[ii:]
        print ii,len(inputList)
        pw = p.apply_async( yourFunc, args = (dir,inputList,qRst) )
    p.close()
    p.join()
