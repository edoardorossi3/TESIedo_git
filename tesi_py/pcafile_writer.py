#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 18 20:29:02 2022

@author: edoardo
"""

import astropy.io.fits as fits
import numpy as np
import os
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA

work_dir=os.getenv('HOME')+'/Desktop/TESI/models/Sandage_varZ_v4.2eq_CB16MILES_ChFall/' #working directory
f_name=work_dir+'sandage_varZ_v4.2eq_spec_dcomb0p25null_001.fits'  #file from which we compute the pc
n_pc=100  #number of pc
file_pca=work_dir+'PCA_sandage_varZ_v4.2eq_spec_dcomb0p25null.fits'  #principal components file

hdul=fits.open(f_name)
N=np.size(hdul[0].data)
y=hdul[1].data
wl=hdul[0].data
sel_wl=(wl<8500) & (wl>3500)
wl=wl[sel_wl]
y=y[...,sel_wl]
N_wl=np.size(wl)
N_models=np.size(hdul[1].data, axis=0)
norm=np.zeros(N_models)
norm2=np.zeros(N_models)

d_wl=np.zeros(np.size(wl)-1)
wl_grid=np.zeros(np.size(wl)-1)
y_interp=np.zeros(np.size(wl)-1)
y_interp2=np.zeros(np.size(wl)-1)

#normalizing 
for i_wl in range (0, np.size(wl)-1):
    d_wl[i_wl]=wl[i_wl+1]-wl[i_wl]
    wl_grid[i_wl]=(wl[i_wl+1]+wl[i_wl])/2.0
        
for i_mod in range(0, N_models):
    y_interp=np.interp(wl_grid, wl, y[i_mod, ...])
    norm[i_mod]=np.dot(y_interp, d_wl)
    
for i_model in range(0, N_models):
    y[i_model, ...]=y[i_model, ...]/norm[i_model]

y=StandardScaler().fit_transform(y)  #standardization: (y-y_mean)/std

pca=PCA(n_components=n_pc)
pca.fit(y)
#y_proj=pca.fit_transform(y)
#y_new=(np.dot(y_proj, pca.components_)*np.std(hdul[1].data, axis=0)+np.mean(hdul[1].data, axis=0))
PC=pca.components_
hdu=fits.PrimaryHDU(data=PC)
hdu.writeto(file_pca, overwrite=True)