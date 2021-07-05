#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul  5 18:34:15 2021

@author: edoardo
"""
def spectre(suffix_file, age50_min, age50_max):
    import numpy as np
    import astropy.io.fits as fits
    import os
    import matplotlib.pyplot as plt
    import matplotlib.cm as cm
    from matplotlib.collections import LineCollection
    from matplotlib.colors import ListedColormap, BoundaryNorm
    ## generate color vector to be used in plot (color=c)
    
    work_dir=os.getenv('HOME')+'/Desktop/TESI/models/Sandage_v4.1_Zfix_noburst_cb16MILES_1M/'
    #name_file='sandage_varZ_v4.1_m62fix_noburst_1M_spec_dcombnull_001.fits'
    
    file_spec=work_dir+suffix_file+'.fits'
    file_par=work_dir+suffix_file+'_physpar_wagef.fits'
    hdul=fits.open(file_spec)
    hdul_par=fits.open(file_par)
    
    wl=hdul[0].data
    spec=hdul[1].data
    age10=hdul_par[1].data['age10']
    age50=hdul_par[1].data['age50']
    age90=hdul_par[1].data['age90']
    
    N_models=np.size(age10)
    idx_array=np.arange(N_models)
    
    dage_norm=np.log10(age10-age90)/age50
    
    wl_range=[3500, 8000]
    wl_norm=[5450, 5550]
    sel_wl=((wl>wl_range[0]) & (wl<wl_range[1]))
    sel_norm=((wl>wl_norm[0]) & (wl<wl_norm[1]))
        
    sel_models=((np.log10(age50)>age50_min) & (np.log10(age50)<age50_max))
    idx_sel=idx_array[sel_models]
    for i_model in range(0, np.size(idx_sel)):
        _norm=np.mean(spec[idx_sel[i_model],sel_norm])
        c=cm.coolwarm((dage_norm[idx_sel[i_model]]-np.min(dage_norm[sel_models]))/(np.max(dage_norm[sel_models])-np.min(dage_norm[sel_models])))
        plt.plot(wl[sel_wl], spec[idx_sel[i_model], sel_wl]/_norm, alpha = 0.01,linewidth = 1, color=c)
        
        
        
    
        
        
        
        