#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul  5 18:34:15 2021

@author: edoardo
"""
def spectre(name_file):
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
    
    file_spec=work_dir+name_file
    hdul=fits.open(file_spec)
    
    wl=hdul[0].data
    spec=hdul[1].data
    
    wl_range=[3500, 8000]
    sel_wl=((wl>wl_range[0]) & (wl<wl_range[1]))
    #norm tra 5450 e 5550
    
        #c = cm.coolwarm((dage_norm-np.min(dage_norm))/(np.max(pc)-np.min(pc)))
    c=cm.coolwarm((1.0-(-1.3))/(2.0-(-1.3)))
    plt.plot(wl[sel_wl], spec[0, sel_wl],alpha = 0.5,linewidth = 1, color=c)