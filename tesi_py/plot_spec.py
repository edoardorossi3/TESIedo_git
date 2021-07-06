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
    from matplotlib.colors import ListedColormap, BoundaryNorm, Normalize
    import matplotlib.colorbar as clb
    import function_plot as f_plt
    ## generate color vector to be used in plot (color=c)
    
    work_dir=os.getenv('HOME')+'/Desktop/TESI/models/Sandage_v4.1_Zfix_noburst_cb16MILES_1M/'
    #name_file='sandage_varZ_v4.1_m62fix_noburst_1M_spec_dcombnull_001.fits'
    snr='200'
    file_spec=work_dir+suffix_file+'_001.fits'
    file_par=work_dir+suffix_file+'_001_physpar_wagef.fits'
    file_time=work_dir+'Time_resol_Zfix1M_SNR'+snr+'_tot.fits'
    file_pert=work_dir+suffix_file+'_perterr_SNR'+snr+'_001.fits'
    
    hdul=fits.open(file_spec)
    hdul_par=fits.open(file_par)
    hdul_time=fits.open(file_time)
    hdul_pert=fits.open(file_pert)
    
    wl=hdul[0].data
    spec=hdul[1].data
    age10=hdul_par[1].data['age10']
    age50=hdul_par[1].data['age50']
    age90=hdul_par[1].data['age90']
    time_res_vector=hdul_time[1].data['Log_d1090n50_min_z62']
    d4000n=hdul_pert[1].data['D4000N_SIGMA200_PERT']
    hdhg=hdul_pert[1].data['HDHG_SIGMA200_PERT']
    hb=hdul_pert[1].data['LICK_HB_SIGMA200_PERT']
    mg2fe=hdul_pert[1].data['MG2FE_SIGMA200_PERT']
    mgfep=hdul_pert[1].data['MGFE_PRIME_SIGMA200_PERT']
    u=hdul_pert[1].data['ABMAG_U_PERT']
    g=hdul_pert[1].data['ABMAG_G_PERT']
    r=hdul_pert[1].data['ABMAG_R_PERT']
    i=hdul_pert[1].data['ABMAG_I_PERT']
    z=hdul_pert[1].data['ABMAG_Z_PERT']
    
    n_bin=int((age50_min-6.0)/0.05)
    N_models=np.size(age10)
    idx_array=np.arange(N_models)
    dage_norm=np.log10((age10-age90)/age50)
    dage_min=-1.3
    dage_max=0.3
    time_res=time_res_vector[n_bin]
    
    wl_range=[3500, 8000]
    wl_norm=[5450, 5550]
    sel_wl=((wl>wl_range[0]) & (wl<wl_range[1]))
    sel_norm=((wl>wl_norm[0]) & (wl<wl_norm[1]))
        
    sel_models=((np.log10(age50)>age50_min) & (np.log10(age50)<age50_max))
    print(np.size(sel_models))
    print(np.sum(sel_models))
    sel_ref=((np.log10(age50)>age50_min) & (np.log10(age50)<age50_max)&(dage_norm<-1.0))
    idx_sel=idx_array[sel_models]
    widths=[5,0.05]
    heights=[1,1]
    
    gs=dict(width_ratios=widths, height_ratios=heights)
    fig, axs=plt.subplots(2,2,figsize=(40,10), gridspec_kw=gs)

    delta=f_plt.chi_q(dage_norm,d4000n, hdhg, hb, mg2fe, mgfep, u-r,g-r, r-i, r-z, sel_models, sel_ref,mkplot=False, sigma_obs=False)
    
    spec_sel=spec[sel_models, ...]
    sel_delta_min=(delta==np.min(delta))
    norm=np.mean(spec_sel[sel_delta_min,sel_norm])
    #axs[0,0].plot(wl[sel_wl], spec_sel[sel_delta_min, sel_wl]/norm, color='black')
    
    for i_model in range(0, np.size(idx_sel)):
        _norm=np.mean(spec[idx_sel[i_model],sel_norm])
        c=cm.rainbow_r((dage_norm[idx_sel[i_model]]-dage_min)/(dage_max-dage_min))
        axs[0,0].plot(wl[sel_wl], spec[idx_sel[i_model], sel_wl]/_norm, alpha = 0.1,linewidth = 0.1, color=c)
        axs[1,0].plot(wl[sel_wl], spec[idx_sel[i_model], sel_wl]/_norm-spec_sel[sel_delta_min, sel_wl]/norm, alpha = 0.1,linewidth = 0.1, color=c)

    
    NORM=Normalize(vmin=dage_min, vmax=dage_max)
    clb.ColorbarBase(axs[0,1], cmap=cm.rainbow_r, norm=NORM)
    axs[0,1].plot([-1,1], [time_res, time_res], linewidth=5, color='black') 
    clb.ColorbarBase(axs[1,1], cmap=cm.rainbow_r, norm=NORM)
    axs[1,1].plot([-1,1], [time_res, time_res], linewidth=5, color='black') 
    
    axs[0,0].set_facecolor('#d8dcd6')
    axs[1,0].set_facecolor('#d8dcd6')
    axs[0,0].tick_params(labelsize=30)
    axs[1,0].tick_params(labelsize=30)
    axs[0,1].tick_params(labelsize=30)
    axs[1,1].tick_params(labelsize=30)

    
    
    
        