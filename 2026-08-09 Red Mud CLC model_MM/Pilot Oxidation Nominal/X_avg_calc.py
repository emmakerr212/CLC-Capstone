# -*- coding: utf-8 -*-
"""
Created on Wed Feb 25 15:00:05 2026

@author: d4mcbrid
"""

import numpy as np
import pandas as pd
import clc_des__param_L as dp
from scipy.integrate import simpson


def _X_avg_(filename, z, r, Rp, L):
    X = pd.read_excel(filename, index_col=[0,1,2])
    Xs = np.full([len(z),len(r)],0,dtype=float)
    for i in range(len(z)):
        for j in range(len(r)):
            Xs[i,j] = X.loc[z[i]].loc[r[j]].loc[1].val
    X_r2 = Xs*(4*np.pi*r[None,:]**2) # spherical weighting  
    # integrate over r
    int_r = simpson(X_r2,r,axis=1)#np.trapezoid(X_r2,r,axis=1)
    
    z -= 0.02
    z /= 0.966667-0.02 # this converts the red mud y coordinates to be between 0 and 1 to get the average conversion
    z = z*L
    int_total = simpson(int_r,z)#np.trapezoid(int_r,z)
    # calculate volumes
    V = 4/3*np.pi*Rp**3*L
    
    X_avg = int_total/V
    return X_avg
    
if __name__ == '__main__':

# Calculating the average conversion of the OC at the end of the reduction stage
    L = 3.08#dp.L
    Rp = dp.Dp/2
    filename = 'X_base_red.xlsx'
    X  = pd.read_excel(filename, index_col=[0,1,2]) # 0: y, 1: r, 2: t
    # z = np.array([0,0.003333,0.01,0.023333,0.05,0.066667,0.1,0.116667,
    #               0.15,0.166667,0.2,0.216667,0.25,0.266667,0.3,0.316667,
    #               0.35,0.358333,0.375,0.383333,0.4,0.416667,0.45,0.466667,
    #               0.5,0.516667,0.55,0.566667,0.6,0.616667,0.65,0.666667,
    #               0.7,0.716667,0.75,0.766667,0.8,0.816667,0.85,0.866667,
    #               0.9,0.916667,0.95,0.966667,1])
    z = np.array([  0.02    , 0.03    , 0.036667,
           0.05    , 0.058333, 0.075   , 0.083333, 0.1     , 0.108333,
           0.125   , 0.133333, 0.15    , 0.166667, 0.2     , 0.216667,
           0.25    , 0.266667, 0.3     , 0.316667, 0.35    , 0.358333,
           0.375   , 0.383333, 0.4     , 0.416667, 0.45    , 0.466667,
           0.5     , 0.516667, 0.55    , 0.566667, 0.6     , 0.616667,
           0.65    , 0.666667, 0.7     , 0.716667, 0.75    , 0.766667,
           0.8     , 0.816667, 0.85    , 0.866667, 0.9     , 0.916667,
           0.95    , 0.966667])

    #0.      , 0.005   , 0.015   ,, 1. 
    r = np.array([0,0.000833,0.001667,0.0025])
    
    Xs = np.full([len(z),len(r)],0,dtype=float)
    
    for i in range(len(z)):
        for j in range(len(r)):
            Xs[i,j] = X.loc[z[i]].loc[r[j]].loc[1].val
    
    X_r2 = Xs*(4*np.pi*r[None,:]**2) # spherical weighting
    
    # integrate over r
    int_r = simpson(X_r2,r,axis=1)#np.trapezoid(X_r2,r,axis=1)
    z -=0.02
    z/=0.966667-0.02 
    print(z)
    z *= L # add dimension back to axial domain
    # integrate over z
    int_total = simpson(int_r,z)#np.trapezoid(int_r,z)
    
    # calculate volumes
    V = 4/3*np.pi*Rp**3*L
    
    X_avg = int_total/V
   # print('end of code')