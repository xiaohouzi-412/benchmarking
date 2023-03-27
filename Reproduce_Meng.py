#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 27 17:23:15 2023

@author: xinjia
"""


import numpy as np
import netCDF4 as nc
from netCDF4 import Dataset
from numpy import NaN
import itertools
from itertools import *

 
# y_all3 is the daily temperature data during the period 1979.01.01-2017.12.31, with the shape (365,39,105), 105 is the amount of nodes

# calculate daily anomalies, for 2 periods; the first part is 1979-1983， and the second is 1984-2017
# (1) the period 1979-1983, for each year, the climatology is the same
t_climatology1=np.empty(shape=[365,105]) 
y1=y_all3[:,0:5,:]  #daily data during 1979.01.01-1983.12.31
t_climatology1=np.average(y1,axis=1)
# calculate anomalies
t_anomalies1=np.empty(shape=[365,5,105]) 
for i in range(0,5):
    t_anomalies1[:,i,:]=y1[:,i,:]-t_climatology1
    
    
# (2) the period during 1984-2017, the climatology is changing with each year, since the data is increasing by time
t_climatology2=np.empty(shape=[365,34,105]) 
for i2 in range(5,39):
    for i3 in range(0,365):
        temp3=y_all3[i3,0:i2+1,:]
        t_climatology2[i3,i2-5,:]=np.average(temp3,axis=0)
    
#calculate anomalies
y2=y_all3[:,5:39,:]  
t_anomalies2=y2-t_climatology2

t_anomalies=np.empty(shape=[365,39,105])
t_anomalies[:,0:5,:]=t_anomalies1
t_anomalies[:,5:39,:]=t_anomalies2

T_anomalies=np.empty(shape=[365*39,105])
for m1 in range(0,39):
    p=m1*365
    q=p+365
    T_anomalies[p:q,:]=t_anomalies[:,m1,:]   


ap=np.arange(0,105)  
pij1=list(itertools.combinations(ap,2))  

ds=[0,31,59,90,120,151,181,212,243,273,304,334] 
dsa=np.empty(shape=[12,39]) 
for i in range(0,12):
    dsa[i,:]=np.arange(ds[i],ds[i]+365*38+1,365)   
    
t1=dsa.reshape(12*39,order='F') 
t=t1[12:468]   

t=t.astype(np.int64) 

#calculate Cij(-tlag)
Cij=np.empty(shape=[456,201,5460])

for i in range(0,len(pij1)):
    aa=pij1[i]   
    i_p=aa[0]  
    i_p=i_p.astype('int64')
    j_p=aa[1] 
    j_p=j_p.astype('int64')
    
    #calculate <Ti(t)>
    A=np.empty(shape=[456])
    for a in range(0,456):
        s=t[a]-364
        T_temp1=T_anomalies[np.arange(s,t[a]+1),i_p]  
        A[a]=(1/365)*np.sum(T_temp1,axis=0)
    #calculate <Tj(t-tlag)>
    B=np.empty(shape=[456,201])
    B1=np.empty(shape=[456])
    for tlag1 in range(0,201):
        for a1 in range(0,456):
            t1=t-tlag1
            s1=t1[a1]-364
            T_temp2=T_anomalies[np.arange(s1,t1[a1]+1),j_p]
            B1[a1]=(1/365)*np.sum(T_temp2,axis=0)
        B[:,tlag1]=B1
        
    Y=np.empty(shape=[456,201])
    for y22 in range(0,201):
        Y[:,y22]=A*B[:,y22]
    
    
    # calculate <Ti(t)*Tj(t-tlag)>
    C11=np.empty(shape=[456])
    C22=np.empty(shape=[456,201])
    for tlag2 in range(0,201):
        for a2 in range(0,456):
            ss=t[a2]-364
            T_temp11=T_anomalies[np.arange(ss,t[a2]+1),i_p]
            t11=t-tlag2
            s11=t11[a2]-364
            T_temp22=T_anomalies[np.arange(s11,t11[a2]+1),j_p]
            C_temp=T_temp11*T_temp22
            C11[a2]=(1/365)*np.sum(C_temp)
        C22[:,tlag2]=C11
    
    D1=np.empty(shape=[456])
    for d in range(0,456):
        d1=t[d]-364
        T_temp1D=T_anomalies[np.arange(d1,t[d]+1),i_p]    
        D1[d]=np.sqrt(np.var(T_temp1D,ddof=1,axis=0))
        
    D2=np.empty(shape=[456,201])
    D22=np.empty(shape=[456])
    for tdlag1 in range(0,201):
        for ad1 in range(0,456):
            td1=t-tdlag1    
            sd1=td1[ad1]-364
            T_tempd2=T_anomalies[np.arange(sd1,td1[ad1]+1),j_p]
            D22[ad1]=np.sqrt(np.var(T_tempd2,ddof=1,axis=0))
        D2[:,tdlag1]=D22
        
    DD_temp1=np.empty(shape=[456,201])
    for dd2 in range(0,201):
        DD_temp1[:,dd2]=D1*D2[:,dd2]
   
    Cij[:,:,i]=(C22-Y)/DD_temp1  
        
    
#calculate Cij(timelag)
ap=np.arange(0,105)   
pij1=list(itertools.combinations(ap,2))  

ds=[0,31,59,90,120,151,181,212,243,273,304,334] 
dsa=np.empty(shape=[12,39]) 
for i in range(0,12):
    dsa[i,:]=np.arange(ds[i],ds[i]+365*38+1,365)   
    
t1=dsa.reshape(12*39,order='F')  
t=t1[12:468]  
t=t.astype(np.int64)     
    
#calculate Cij   
Cij1=np.empty(shape=[456,201,5460])
for i in range(0,len(pij1)):
    aa=pij1[i]   
    i_p=aa[0]  
    i_p=i_p.astype('int64')
    j_p=aa[1]  
    j_p=j_p.astype('int64')
    
    # calculate <Tj(t)>
    A=np.empty(shape=[456])
    for a in range(0,456):
        s=t[a]-364
        T_temp1=T_anomalies[np.arange(s,t[a]+1),j_p]  
        A[a]=(1/365)*np.sum(T_temp1,axis=0)
    # calculate <Ti(t-tlag)>
    B=np.empty(shape=[456,201])
    B1=np.empty(shape=[456])
    for tlag1 in range(0,201):
        for a1 in range(0,456):
            t1=t-tlag1
            s1=t1[a1]-364
            T_temp2=T_anomalies[np.arange(s1,t1[a1]+1),i_p]
            B1[a1]=(1/365)*np.sum(T_temp2,axis=0)
        B[:,tlag1]=B1
        
    Y=np.empty(shape=[456,201])
    for y22 in range(0,201):
        Y[:,y22]=A*B[:,y22]
    
    
    #calculate <Tj(t)*Ti(t-tlag)>
    C11=np.empty(shape=[456])
    C22=np.empty(shape=[456,201])
    for tlag2 in range(0,201):
        for a2 in range(0,456):
            ss=t[a2]-364
            T_temp11=T_anomalies[np.arange(ss,t[a2]+1),j_p]
            t11=t-tlag2
            s11=t11[a2]-364
            T_temp22=T_anomalies[np.arange(s11,t11[a2]+1),i_p]
            C_temp=T_temp11*T_temp22
            C11[a2]=(1/365)*np.sum(C_temp)
        C22[:,tlag2]=C11
    
    D1=np.empty(shape=[456])
    for d in range(0,456):
        d1=t[d]-364
        T_temp1D=T_anomalies[np.arange(d1,t[d]+1),j_p]      
        D1[d]=np.sqrt(np.var(T_temp1D,ddof=1,axis=0))
        
    D2=np.empty(shape=[456,201])
    D22=np.empty(shape=[456])
    for tdlag1 in range(0,201):
        for ad1 in range(0,456):
            td1=t-tdlag1    
            sd1=td1[ad1]-364
            T_tempd2=T_anomalies[np.arange(sd1,td1[ad1]+1),i_p]
            D22[ad1]=np.sqrt(np.var(T_tempd2,ddof=1,axis=0))
        D2[:,tdlag1]=D22
        
    DD_temp1=np.empty(shape=[456,201])
    for dd2 in range(0,201):
        DD_temp1[:,dd2]=D1*D2[:,dd2]
   
    Cij1[:,:,i]=(C22-Y)/DD_temp1  
        

Cijtotal=np.concatenate((Cij, Cij1), axis = 1)
Cijmax1=Cijtotal.max(axis=1)  
Ct1=np.sum(Cijmax1,axis=1)
Ct=1/5460*Ct1

#calculate FI(t)
Ct_1=Ct[12:456]    
f1=np.empty(shape=[408])
j=0
for i in range(36,444):
    c1=Ct_1[0:i]
    c11=np.log(c1)
    f1[j]=1/(i+1)*sum(c11)
    j=j+1

Ct_2=Ct[48:456]  
f2=np.log(Ct_2)
FI=np.empty(shape=[408])
FI=f1-f2














    
   
    
    
    
    
    
    
    