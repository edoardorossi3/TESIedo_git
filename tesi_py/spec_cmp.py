#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 22 14:57:19 2021

@author: edoardo
"""

import numpy as np
import os
import atropy.io.fits as fits
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA

#SETTINGS
i_chunk=1
work_dir=os.getenv('HOME')+'/Desktop/TESI/models/Sandage_varZ_v4.2eq_CB16MILES_ChFall/' #working directory
f_name=work_dir+'sandage_varZ_v4.2eq_spec_dcomb0p25null_{:03d}.fits' #file to be compressed
n_pc=100  #number of pc
file_pca='PCA_sandage_varZ_v4.2eq_spec_dcomb0p25null.fits'
file_cmp='sandage_varZ_v4.2eq_spec_dcomb0p25null_{:03d}_compressed.fits'

hdul=fits.open(f_name)
N=np.size(hdul[0].data)

y=hdul[1].data
StandardScaler().fit(y)
y=StandardScaler().fit_transform(y)  #standardization: (y-y_mean)/std

pca=PCA(n_components=n_pc)
y_proj=pca.fit_transform(y)
y_new=(np.dot(y_proj, pca.components_)*np.std(hdul[1].data, axis=0)+np.mean(hdul[1].data, axis=0))
PC=pca.components_