#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 27 17:17:00 2023

@author: xinjia
"""


import numpy as np
import xlrd
import pandas as pd


data=pd.read_excel('C:/work/move.xlsx',sheet_name='Sheet3')
F = data['st']   # time series of S(t)
O = data['threshold']  # threshold for S(t)
F1=np.empty(shape=[408,1])
O1=np.empty(shape=[408,1])

for i in range(0,408):
    F1[i]=F[i]
    O1[i]=O[i]  

w=np.arange(12,408,12)  #window length for each movement 
F_move=np.empty(shape=[408,1])
F_all=np.empty(shape=[408,33])   # move 1 year, 2 years, 3 years...each time. Put the moved result into F_all

for i in w:
    F_move[:i]=F1[-i:]
    F_move[i:]=F1[:-i]
    j=int(i/12-1)
    F_all[:,j]=F_move[:,0]   # we get all the shifted series
    


