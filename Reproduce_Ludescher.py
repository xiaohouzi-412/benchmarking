#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 27 16:17:44 2023

@author: xinjia

Reproduction of Ludescher's work


"""

import numpy as np
import netCDF4 as nc
from netCDF4 import Dataset

#obtain the information of grid points in the outside basin
T_outside1=np.empty(shape=[365,207,63])
T_outside11=np.empty(shape=[365,207])
   

lon1=np.arange(48,117,3)  #longitude index of j
lat1=np.arange(24,51,3)  #latitude index of j
k1=0

for i in range(1949,2012):
    for j in range(1,366):
        file_path1='D:\\reproduce\\data\\tanomalies\\tanomalies_'+str(j)+'_'+str(i)+'.nc'
        b=Dataset(file_path1)
        to=np.array(b.variables['tanomalies'])
        to_j=to[lat1]   
        to_j=to_j[:,lon1]   
        to_j1=to_j.reshape((1,207))
        T_outside11[j-1,:]=to_j1    
    T_outside1[:,:,k1]=T_outside11  
    k1=k1+1
        
T_outside2=np.empty(shape=[365*63,207])
for m1 in range(0,63):
    p=m1*365
    q=p+365
    T_outside2[p:q,:]=T_outside1[:,:,m1]     


T_elbasin=T_outside2[:,[102,103,104,105,106,107,108,109,110,111,112,113,130,136]]  # grid points in elnino basin (14 in total)， during 1949.01.01-2011.12.31
T_outside=np.delete(T_outside2,[102,103,104,105,106,107,108,109,110,111,112,113,130,136],axis=1) # grid points in outside area (193 in total)

#开始计算Cij
t=np.arange(732,22995,99)  #each 100th day during 1951.01.01-2011.12.31 
# in ordor to calculate Cij(-tao), we need to first construct Ti(t)和Tj(t-tlag),Ti is the girds of T_elbasin,Tj is T_outside
Ti=T_elbasin[t,:]   
Tj=np.empty(shape=[225,193,201])
for tlag in range(0,201):
    k=t-tlag
    Tj[:,:,tlag]=T_outside[k,:]


#calculate <Ti(t)>,<Tj(t-tlag)>,<Ti*Tj(t-tlag)>, seperately
A=np.empty(shape=[225,14])   #A is <Ti(t)>
T_temp1=np.empty(shape=[365,14])   
#calculate <Ti(t)>
for a in range(0,225):
    s=t[a]-364
    T_temp1=T_elbasin[np.arange(s,t[a]+1),:]      
    A[a,:]=(1/365)*np.sum(T_temp1,axis=0)

#calculate <Tj(t-tlag)>
B=np.empty(shape=[225,193,201])
B1=np.empty(shape=[225,193])
for tlag1 in range(0,201):
    for a1 in range(0,225):
        t1=t-tlag1    
        s1=t1[a1]-364
        T_temp2=T_outside[np.arange(s1,t1[a1]+1),:]
        B1[a1,:]=(1/365)*np.sum(T_temp2,axis=0)
    B[:,:,tlag1]=B1

#Calculate <Ti(t)*Tj(t-tlag)>
C=np.empty(shape=[225,201,193,14])
C11=np.empty(shape=[225])
C22=np.empty(shape=[225,201])
C33=np.empty(shape=[225,201,193])

for i1 in range(0,14):
    for j1 in range(0,193):
        for tlag2 in range(0,201):
            for a2 in range(0,225):
                ss=t[a2]-364
                T_temp11=T_elbasin[np.arange(ss,t[a2]+1),i1]
                t11=t-tlag2
                s11=t11[a2]-364
                T_temp22=T_outside[np.arange(s11,t11[a2]+1),j1]
                C_temp=T_temp11*T_temp22
                C11[a2]=(1/365)*np.sum(C_temp)
            C22[:,tlag2]=C11
        C33[:,:,j1]=C22
    C[:,:,:,i1]=C33
    

D1=np.empty(shape=[225,14])
for d in range(0,225):
    d1=t[d]-364
    T_temp1D=T_elbasin[np.arange(d1,t[d]+1),:]      
    D1[d,:]=np.sqrt(np.var(T_temp1D,ddof=1,axis=0))
    

D2=np.empty(shape=[225,193,201])
D22=np.empty(shape=[225,193])
for tdlag1 in range(0,201):
    for ad1 in range(0,225):
        td1=t-tdlag1    
        sd1=td1[ad1]-364
        T_tempd2=T_outside[np.arange(sd1,td1[ad1]+1),:]
        D22[ad1,:]=np.sqrt(np.var(T_tempd2,ddof=1,axis=0))
    D2[:,:,tdlag1]=D22


#calculate Cij(-tlag)
Cij=np.empty(shape=[225,201,193*14])
C_new=C.reshape((225,201,193*14))

#calculate A*B
Y=np.empty(shape=[225,201,193,14])
Y_temp1=np.empty(shape=[225,201])
Y_temp2=np.empty(shape=[225,201,193])
for y1 in range(0,14):
    for y11 in range(0,193):
        for y22 in range(0,201):
            Y_temp1[:,y22]=A[:,y1]*B[:,y11,y22]
        Y_temp2[:,:,y11]=Y_temp1
    Y[:,:,:,y1]=Y_temp2
   
Y_new=Y.reshape(225,201,193*14)

#calculate D1*D2
Dd=np.empty(shape=[225,201,193,14])
DD_temp1=np.empty(shape=[225,201])
DD_temp2=np.empty(shape=[225,201,193])
for dd in range(0,14):
    for dd1 in range(0,193):
        for dd2 in range(0,201):
            DD_temp1[:,dd2]=D1[:,dd]*D2[:,dd1,dd2]
        DD_temp2[:,:,dd1]=DD_temp1
    Dd[:,:,:,dd]=DD_temp2
            
Dd_new=Dd.reshape(225,201,193*14)

Cij=(C_new-Y_new)/Dd_new   #Cij is the matrix of (225,201,193*14)
#calculate Sij
mean_Cij=np.mean(Cij,axis=1)
max_Cij=np.max(Cij,axis=1)
sd_Cij=np.std(Cij,axis=1)
Sij=(max_Cij-mean_Cij)/sd_Cij
St=np.mean(Sij,axis=1)