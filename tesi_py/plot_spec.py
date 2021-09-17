#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul  5 18:34:15 2021

@author: edoardo
"""
def spectre(suffix_file, age50_min, age50_max, dage_n_max=None):
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
    
    #work_dir=os.getenv('HOME')+'/Desktop/TESI/models/Sandage_v4.1_Zfix_noburst_cb16MILES_1M/'
    work_dir='/export/home/extragal/zibetti/no_ownCloud/SteMaGE/data/SEDlibraries/Sandage_v4.1_Zfix_noburst_cb16MILES_1M/'
    #name_file='sandage_varZ_v4.1_m62fix_noburst_1M_spec_dcombnull_001.fits'
    snr='500'
    SNR=200.0
    file_spec=work_dir+suffix_file+'_001.fits'
    file_par=work_dir+suffix_file+'_001_physpar_wagef.fits'
    file_pert=work_dir+suffix_file+'_perterr_SNR'+snr+'_001.fits'
    
    file_time5=work_dir+'Time_resol_Zfix1M_SNR5_tot.fits'
    file_time10=work_dir+'Time_resol_Zfix1M_SNR10_tot.fits'
    file_time20=work_dir+'Time_resol_Zfix1M_SNR20_tot.fits'
    file_time50=work_dir+'Time_resol_Zfix1M_SNR50_tot.fits'
    file_time100=work_dir+'Time_resol_Zfix1M_SNR100_tot.fits'
    file_time200=work_dir+'Time_resol_Zfix1M_SNR200_tot.fits'
    file_time500=work_dir+'Time_resol_Zfix1M_SNR500_tot.fits'

    
    hdul=fits.open(file_spec)
    hdul_par=fits.open(file_par)
    hdul_pert=fits.open(file_pert)
    
    hdul_time5=fits.open(file_time5)
    hdul_time10=fits.open(file_time10)
    hdul_time20=fits.open(file_time20)
    hdul_time50=fits.open(file_time50)
    hdul_time100=fits.open(file_time100)
    hdul_time200=fits.open(file_time200)
    hdul_time500=fits.open(file_time500)

    
    wl=hdul[0].data
    age10=hdul_par[1].data['age10']
    age50=hdul_par[1].data['age50']
    age90=hdul_par[1].data['age90']
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
    
    time_res_vector5=hdul_time5[1].data['Log_d1090n50_min_z62']
    time_res_vector10=hdul_time10[1].data['Log_d1090n50_min_z62']
    time_res_vector20=hdul_time20[1].data['Log_d1090n50_min_z62']
    time_res_vector50=hdul_time50[1].data['Log_d1090n50_min_z62']
    time_res_vector100=hdul_time100[1].data['Log_d1090n50_min_z62']
    time_res_vector200=hdul_time200[1].data['Log_d1090n50_min_z62']
    time_res_vector500=hdul_time500[1].data['Log_d1090n50_min_z62']

    
    #starting spectra
    print('age50:', age50[0])  
    #spec_sel=np.arange(np.size(hdul[1].data[0,...]))
    #for j in range(0, np.size(age50)):
    #    if ((np.log10(age50[j])<age50_max) & (np.log10(age50[j])>age50_min)):
    #        spec_sel=np.array([spec_sel, hdul[1].data[j,...]])
    
    print('checkpoint 1')
    #print('shape spec:', np.shape(spec_sel))
    
    for i_chunk in range (2,51):
        _file_spec=work_dir+suffix_file+'_{:03d}.fits'
        _file_par=work_dir+suffix_file+'_{:03d}_physpar_wagef.fits'
        _file_pert=work_dir+suffix_file+'_perterr_SNR'+snr+'_{:03d}.fits'

        _hdul_par=fits.open(_file_par.format(i_chunk))
        _hdul_pert=fits.open(_file_pert.format(i_chunk))
        
        
        _age10=_hdul_par[1].data['age10']
        _age50=_hdul_par[1].data['age50']
        _age90=_hdul_par[1].data['age90']
        _d4000n=_hdul_pert[1].data['D4000N_SIGMA200_PERT']
        _hdhg=_hdul_pert[1].data['HDHG_SIGMA200_PERT']
        _hb=_hdul_pert[1].data['LICK_HB_SIGMA200_PERT']
        _mg2fe=_hdul_pert[1].data['MG2FE_SIGMA200_PERT']
        _mgfep=_hdul_pert[1].data['MGFE_PRIME_SIGMA200_PERT']
        _u=_hdul_pert[1].data['ABMAG_U_PERT']
        _g=_hdul_pert[1].data['ABMAG_G_PERT']
        _r=_hdul_pert[1].data['ABMAG_R_PERT']
        _i=_hdul_pert[1].data['ABMAG_I_PERT']
        _z=_hdul_pert[1].data['ABMAG_Z_PERT']
        
        #for k in range(0, np.size(_age50)):
        #    if ((np.log10(_age50[k])<age50_max) & (np.log10(_age50[k])>age50_min)):
        #        spec_sel=np.array([spec_sel, _hdul_spec[1].data[k,...]])
        
        age50=np.append(age50, _age50)
        age10=np.append(age10, _age10)
        age90=np.append(age90, _age90)
        d4000n=np.append(d4000n, _d4000n)
        hdhg=np.append(hdhg, _hdhg)
        hb=np.append(hb, _hb)
        mg2fe=np.append(mg2fe, _mg2fe)
        mgfep=np.append(mgfep, _mgfep)
        u=np.append(u, _u)
        g=np.append(g, _g)
        r=np.append(r, _r)
        i=np.append(i, _i)
        z=np.append(z, _z)
        
    dage_norm=np.log10((age10-age90)/age50)

    if (dage_n_max!=None):
        sel_models=((np.log10(age50)<age50_max) & (np.log10(age50)>age50_min) &(dage_norm<dage_n_max))
    else:
        sel_models=((np.log10(age50)<age50_max) & (np.log10(age50)>age50_min))


    
    N_row=np.sum(sel_models)
    N_col=np.size(hdul[1].data[0,...])
    matrix_spec=np.zeros(shape=(N_row, N_col))
    #spec_sel=np.delete(spec_sel, 0, 0)  #delete the first row    
    
    n_bin=int((age50_min-6.0)/0.05)
    N_models=np.size(age10)
    idx_array=np.arange(N_models)
    dage_min=-1.3
    dage_max=0.3
    time_res5=time_res_vector5[n_bin]
    time_res10=time_res_vector10[n_bin]
    time_res20=time_res_vector20[n_bin]
    time_res50=time_res_vector50[n_bin]
    time_res100=time_res_vector100[n_bin]
    time_res200=time_res_vector200[n_bin]
    time_res500=time_res_vector500[n_bin]
    
    wl_range=[3500, 8000]
    #wl_range=[4827,4892]
    wl_norm=[5450, 5550]
    #wl_norm=[4827, 4848, 4876, 4892]
    sel_wl=((wl>wl_range[0]) & (wl<wl_range[1]))
    sel_norm=((wl>wl_norm[0]) & (wl<wl_norm[1]))
    #sel_norm=(((wl>wl_norm[0]) & (wl<wl_norm[1])) | ((wl>wl_norm[2]) & (wl<wl_norm[3])))    
    
    #sel_models=((np.log10(age50)>age50_min) & (np.log10(age50)<age50_max))
    #print(np.size(sel_models))
    #print(np.sum(sel_models))
    sel_ref=((np.log10(age50)>age50_min) & (np.log10(age50)<age50_max)&(dage_norm<-1.0))
    idx_sel=idx_array[sel_models]
    widths=[5,0.15]
    #widths=[3,0.15,3]
    heights=[1,1]
    
    #gs=dict(width_ratios=widths, height_ratios=heights)
    cm2inch = 1/2.54 
    #fig, axs=plt.subplots(2,2,figsize=(30*cm2inch,20*cm2inch), gridspec_kw=gs)
    fig = plt.figure(figsize=(30*cm2inch,20*cm2inch))
    #gs = fig.add_gridspec(2, 2, height_ratios=heights, width_ratios=widths, hspace=0.2, wspace=0.3)
    gs = fig.add_gridspec(2, 2, height_ratios=heights, width_ratios=widths, hspace=0.2, wspace=0.4)
    
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[1, 0])
    ax3 = fig.add_subplot(gs[:, 1])
    #ax4 = fig.add_subplot(gs[0,2])
    #ax5 = fig.add_subplot(gs[1,2])

    #clb.ColorbarBase(ax3, cmap=cm.rainbow_r)
    
    print('size r:', np.size(r))
    print('size i:', np.size(i))

    delta=f_plt.chi_q(dage_norm,d4000n, hdhg, hb, mg2fe, mgfep, u-r,g-r, r-i, r-z, sel_models, sel_ref,mkplot=False, sigma_obs=False)
    first_row=0
    for i_chunk in range(0,50):
        _file_spec=work_dir+suffix_file+'_{:03d}.fits'
        _hdul_spec=fits.open(_file_spec.format(i_chunk+1))
        _sel_mod=sel_models[i_chunk*20000:(i_chunk+1)*20000]
        _spec=_hdul_spec[1].data
        last_row=np.sum(_sel_mod)+first_row
        matrix_spec[first_row:last_row,...]=_spec[_sel_mod,...]
        
        print('first_row:', first_row,'last_row:', last_row)
        first_row=last_row
    
    print('n_row:', N_row)
    #return
    #spec_sel=spec[sel_models, ...]
    #sel_delta_min=(delta==np.min(delta))
    arg_delta_min=np.argmin(delta)
    #print('first spec:', np.shape(hdul[1].data[0,...]))
    #norm=np.mean(matrix_spec[sel_delta_min,sel_norm])
    
    norm=np.mean(matrix_spec[arg_delta_min,sel_norm]) #fluxes
    
    for i_model in range(0, N_row):
        #_norm=np.mean(spec[idx_sel[i_model],sel_norm])
        
        _norm=np.mean(matrix_spec[i_model, sel_norm])
        
        c=cm.rainbow_r((dage_norm[idx_sel[i_model]]-dage_min)/(dage_max-dage_min))
        #axs[0,0].plot(wl[sel_wl], spec[idx_sel[i_model], sel_wl]/_norm, alpha = 0.1,linewidth = 0.1, color=c)
        #axs[1,0].plot(wl[sel_wl], spec[idx_sel[i_model], sel_wl]/_norm-spec_sel[sel_delta_min, sel_wl]/norm, alpha = 0.1,linewidth = 0.1, color=c)
        ax1.set_ylim([0.0,1.7])
        ax1.plot(wl[sel_wl], matrix_spec[i_model, sel_wl]/_norm, alpha = 0.1,linewidth = 0.05, color=c)
        ax2.plot(wl[sel_wl], (matrix_spec[i_model, sel_wl]/_norm)/(matrix_spec[arg_delta_min, sel_wl]/norm), alpha = 0.1,linewidth = 0.05, color=c)
        
    
    NORM=Normalize(vmin=dage_min, vmax=dage_max)
    clb.ColorbarBase(ax3, cmap=cm.rainbow_r, norm=NORM)
    ax3.plot([-1,1], [time_res5, time_res5], linewidth=1, color='black') 
    ax3.plot([-1,1], [time_res10, time_res10], linewidth=1, color='black') 
    ax3.plot([-1,1], [time_res20, time_res20], linewidth=1, color='black') 
    ax3.plot([-1,1], [time_res50, time_res50], linewidth=1, color='black') 
    ax3.plot([-1,1], [time_res100, time_res100], linewidth=1, color='black') 
    ax3.plot([-1,1], [time_res200, time_res200], linewidth=1, color='black') 
    ax3.plot([-1,1], [time_res500, time_res500], linewidth=1, color='black') 
    
    ax1.set_facecolor('#d8dcd6')
    ax2.set_facecolor('#d8dcd6')
    ax1.tick_params(labelsize=15)
    ax2.tick_params(labelsize=15)
    ax3.tick_params(labelsize=10)
    ax2.set_ylim([0.95,1.05])
    #axs[1,0].plot([np.min(wl[sel_wl]), np.max(wl[sel_wl])], np.array([1.0,1.0])/SNR)
    #axs[1,0].plot([np.min(wl[sel_wl]), np.max(wl[sel_wl])], -np.array([1.0,1.0])/SNR)
    ax1.minorticks_on()
    ax2.minorticks_on()
    
    ax2.set_xlabel(r'$\lambda [\AA]$', size=20)
    ax2.set_ylabel(r'$L_{\lambda}/L_{\lambda,ref}$', size=20)
    ax1.set_ylabel(r'$L_{\lambda}/L_{\lambda,cont}$', size=20)
    ax3.set_ylabel(r'$log_{10}(\Delta age_{n})$', size=20)
    
    #arrow
    ax3.annotate('S/N=5',
    xy=(-1, time_res5),
    xycoords='data',
    xytext=(-3, time_res5),
    arrowprops=
    dict(facecolor='darkred', shrink=0.05),
    horizontalalignment='right',
    verticalalignment='center', size=15)
    
    ax3.annotate('',
    xy=(-1, time_res10),
    xycoords='data',
    xytext=(-3, time_res10),
    arrowprops=
    dict(facecolor='red', shrink=0.05),
    horizontalalignment='right',
    verticalalignment='center')
    
    ax3.annotate('',
    xy=(-1, time_res20),
    xycoords='data',
    xytext=(-3, time_res20),
    arrowprops=
    dict(facecolor='orange', shrink=0.05),
    horizontalalignment='right',
    verticalalignment='center')
    
    ax3.annotate('S/N=50',
    xy=(-1, time_res50),
    xycoords='data',
    xytext=(-3, time_res50),
    arrowprops=
    dict(facecolor='yellow', shrink=0.05),
    horizontalalignment='right',
    verticalalignment='center', size=15)
    
    ax3.annotate('',
    xy=(-1, time_res100),
    xycoords='data',
    xytext=(-3, time_res100),
    arrowprops=
    dict(facecolor='green', shrink=0.05),
    horizontalalignment='right',
    verticalalignment='center')
    
    ax3.annotate('',
    xy=(-1, time_res200),
    xycoords='data',
    xytext=(-3, time_res200),
    arrowprops=
    dict(facecolor='blue', shrink=0.05),
    horizontalalignment='right',
    verticalalignment='center')
    
    ax3.annotate('S/N=500',
    xy=(-1, time_res500),
    xycoords='data',
    xytext=(-3, time_res500),
    arrowprops=
    dict(facecolor='darkblue', shrink=0.05),
    horizontalalignment='right',
    verticalalignment='center', size=15)
     
    return fig


def spec_indices(suffix_file, age50_min, age50_max, wl_range1, wl_range2, wl_norm1, wl_norm2, title1='', title2=''):
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
    
    #work_dir=os.getenv('HOME')+'/Desktop/TESI/models/Sandage_v4.1_Zfix_noburst_cb16MILES_1M/'
    work_dir='/export/home/extragal/zibetti/no_ownCloud/SteMaGE/data/SEDlibraries/Sandage_v4.1_Zfix_noburst_cb16MILES_1M/'
    #name_file='sandage_varZ_v4.1_m62fix_noburst_1M_spec_dcombnull_001.fits'
    snr='500'
    SNR=200.0
    file_spec=work_dir+suffix_file+'_001.fits'
    file_par=work_dir+suffix_file+'_001_physpar_wagef.fits'
    file_pert=work_dir+suffix_file+'_perterr_SNR'+snr+'_001.fits'
    
    file_time5=work_dir+'Time_resol_Zfix1M_SNR5_tot.fits'
    file_time10=work_dir+'Time_resol_Zfix1M_SNR10_tot.fits'
    file_time20=work_dir+'Time_resol_Zfix1M_SNR20_tot.fits'
    file_time50=work_dir+'Time_resol_Zfix1M_SNR50_tot.fits'
    file_time100=work_dir+'Time_resol_Zfix1M_SNR100_tot.fits'
    file_time200=work_dir+'Time_resol_Zfix1M_SNR200_tot.fits'
    file_time500=work_dir+'Time_resol_Zfix1M_SNR500_tot.fits'

    
    hdul=fits.open(file_spec)
    hdul_par=fits.open(file_par)
    hdul_pert=fits.open(file_pert)
    
    hdul_time5=fits.open(file_time5)
    hdul_time10=fits.open(file_time10)
    hdul_time20=fits.open(file_time20)
    hdul_time50=fits.open(file_time50)
    hdul_time100=fits.open(file_time100)
    hdul_time200=fits.open(file_time200)
    hdul_time500=fits.open(file_time500)

    
    wl=hdul[0].data
    age10=hdul_par[1].data['age10']
    age50=hdul_par[1].data['age50']
    age90=hdul_par[1].data['age90']
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
    
    time_res_vector5=hdul_time5[1].data['Log_d1090n50_min_z62']
    time_res_vector10=hdul_time10[1].data['Log_d1090n50_min_z62']
    time_res_vector20=hdul_time20[1].data['Log_d1090n50_min_z62']
    time_res_vector50=hdul_time50[1].data['Log_d1090n50_min_z62']
    time_res_vector100=hdul_time100[1].data['Log_d1090n50_min_z62']
    time_res_vector200=hdul_time200[1].data['Log_d1090n50_min_z62']
    time_res_vector500=hdul_time500[1].data['Log_d1090n50_min_z62']

    
    #starting spectra
    print('age50:', age50[0])  
    #spec_sel=np.arange(np.size(hdul[1].data[0,...]))
    #for j in range(0, np.size(age50)):
    #    if ((np.log10(age50[j])<age50_max) & (np.log10(age50[j])>age50_min)):
    #        spec_sel=np.array([spec_sel, hdul[1].data[j,...]])
    
    print('checkpoint 1')
    #print('shape spec:', np.shape(spec_sel))
    
    for i_chunk in range (2,51):
        _file_spec=work_dir+suffix_file+'_{:03d}.fits'
        _file_par=work_dir+suffix_file+'_{:03d}_physpar_wagef.fits'
        _file_pert=work_dir+suffix_file+'_perterr_SNR'+snr+'_{:03d}.fits'

        _hdul_par=fits.open(_file_par.format(i_chunk))
        _hdul_pert=fits.open(_file_pert.format(i_chunk))
        
        
        _age10=_hdul_par[1].data['age10']
        _age50=_hdul_par[1].data['age50']
        _age90=_hdul_par[1].data['age90']
        _d4000n=_hdul_pert[1].data['D4000N_SIGMA200_PERT']
        _hdhg=_hdul_pert[1].data['HDHG_SIGMA200_PERT']
        _hb=_hdul_pert[1].data['LICK_HB_SIGMA200_PERT']
        _mg2fe=_hdul_pert[1].data['MG2FE_SIGMA200_PERT']
        _mgfep=_hdul_pert[1].data['MGFE_PRIME_SIGMA200_PERT']
        _u=_hdul_pert[1].data['ABMAG_U_PERT']
        _g=_hdul_pert[1].data['ABMAG_G_PERT']
        _r=_hdul_pert[1].data['ABMAG_R_PERT']
        _i=_hdul_pert[1].data['ABMAG_I_PERT']
        _z=_hdul_pert[1].data['ABMAG_Z_PERT']
        
        #for k in range(0, np.size(_age50)):
        #    if ((np.log10(_age50[k])<age50_max) & (np.log10(_age50[k])>age50_min)):
        #        spec_sel=np.array([spec_sel, _hdul_spec[1].data[k,...]])
        
        age50=np.append(age50, _age50)
        age10=np.append(age10, _age10)
        age90=np.append(age90, _age90)
        d4000n=np.append(d4000n, _d4000n)
        hdhg=np.append(hdhg, _hdhg)
        hb=np.append(hb, _hb)
        mg2fe=np.append(mg2fe, _mg2fe)
        mgfep=np.append(mgfep, _mgfep)
        u=np.append(u, _u)
        g=np.append(g, _g)
        r=np.append(r, _r)
        i=np.append(i, _i)
        z=np.append(z, _z)
        
    
    sel_models=((np.log10(age50)<age50_max) & (np.log10(age50)>age50_min))
    N_row=np.sum(sel_models)
    N_col=np.size(hdul[1].data[0,...])
    matrix_spec=np.zeros(shape=(N_row, N_col))
    #spec_sel=np.delete(spec_sel, 0, 0)  #delete the first row    
    
    n_bin=int((age50_min-6.0)/0.05)
    N_models=np.size(age10)
    idx_array=np.arange(N_models)
    dage_norm=np.log10((age10-age90)/age50)
    dage_min=-1.3
    dage_max=0.3
    time_res5=time_res_vector5[n_bin]
    time_res10=time_res_vector10[n_bin]
    time_res20=time_res_vector20[n_bin]
    time_res50=time_res_vector50[n_bin]
    time_res100=time_res_vector100[n_bin]
    time_res200=time_res_vector200[n_bin]
    time_res500=time_res_vector500[n_bin]
    
    
    
    sel_wl1=((wl>wl_range1[0]) & (wl<wl_range1[1]))
    sel_norm1=(((wl>wl_norm1[0]) & (wl<wl_norm1[1])) | ((wl>wl_norm1[2]) & (wl<wl_norm1[3])))    
    
    sel_wl2=((wl>wl_range2[0]) & (wl<wl_range2[1]))
    sel_norm2=(((wl>wl_norm2[0]) & (wl<wl_norm2[1])) | ((wl>wl_norm2[2]) & (wl<wl_norm2[3])))    
    
    #sel_models=((np.log10(age50)>age50_min) & (np.log10(age50)<age50_max))
    #print(np.size(sel_models))
    #print(np.sum(sel_models))
    sel_ref=((np.log10(age50)>age50_min) & (np.log10(age50)<age50_max)&(dage_norm<-1.0))
    idx_sel=idx_array[sel_models]
    widths=[3,0.15,3]
    heights=[1,1]
    
    #gs=dict(width_ratios=widths, height_ratios=heights)
    cm2inch = 1/2.54 
    #fig, axs=plt.subplots(2,2,figsize=(30*cm2inch,20*cm2inch), gridspec_kw=gs)
    fig = plt.figure(figsize=(30*cm2inch,15*cm2inch))
    #gs = fig.add_gridspec(2, 2, height_ratios=heights, width_ratios=widths, hspace=0.2, wspace=0.3)
    gs = fig.add_gridspec(2, 3, height_ratios=heights, width_ratios=widths, hspace=0.2, wspace=0.8)
    
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[1, 0])
    ax3 = fig.add_subplot(gs[:, 1])
    ax4 = fig.add_subplot(gs[0,2])
    ax5 = fig.add_subplot(gs[1,2])

    #clb.ColorbarBase(ax3, cmap=cm.rainbow_r)
    
    print('size r:', np.size(r))
    print('size i:', np.size(i))

    delta=f_plt.chi_q(dage_norm,d4000n, hdhg, hb, mg2fe, mgfep, u-r,g-r, r-i, r-z, sel_models, sel_ref,mkplot=False, sigma_obs=False)
    first_row=0
    for i_chunk in range(0,50):
        _file_spec=work_dir+suffix_file+'_{:03d}.fits'
        _hdul_spec=fits.open(_file_spec.format(i_chunk+1))
        _sel_mod=sel_models[i_chunk*20000:(i_chunk+1)*20000]
        _spec=_hdul_spec[1].data
        last_row=np.sum(_sel_mod)+first_row
        matrix_spec[first_row:last_row,...]=_spec[_sel_mod,...]
        
        print('first_row:', first_row,'last_row:', last_row)
        first_row=last_row
    
    print('n_row:', N_row)
    #return
    #spec_sel=spec[sel_models, ...]
    #sel_delta_min=(delta==np.min(delta))
    arg_delta_min=np.argmin(delta)
    #print('first spec:', np.shape(hdul[1].data[0,...]))
    #norm=np.mean(matrix_spec[sel_delta_min,sel_norm])
    
    norm1=np.mean(matrix_spec[arg_delta_min,sel_norm1]) #fluxes
    norm2=np.mean(matrix_spec[arg_delta_min,sel_norm2])
    for i_model in range(0, N_row):
        #_norm=np.mean(spec[idx_sel[i_model],sel_norm])
        
        _norm1=np.mean(matrix_spec[i_model, sel_norm1])
        _norm2=np.mean(matrix_spec[i_model, sel_norm2])
        
        c=cm.rainbow_r((dage_norm[idx_sel[i_model]]-dage_min)/(dage_max-dage_min))
        #axs[0,0].plot(wl[sel_wl], spec[idx_sel[i_model], sel_wl]/_norm, alpha = 0.1,linewidth = 0.1, color=c)
        #axs[1,0].plot(wl[sel_wl], spec[idx_sel[i_model], sel_wl]/_norm-spec_sel[sel_delta_min, sel_wl]/norm, alpha = 0.1,linewidth = 0.1, color=c)
        ax1.plot(wl[sel_wl1], matrix_spec[i_model, sel_wl1]/_norm1, alpha = 0.1,linewidth = 0.05, color=c)
        ax2.plot(wl[sel_wl1], (matrix_spec[i_model, sel_wl1]/_norm1)/(matrix_spec[arg_delta_min, sel_wl1]/norm1), alpha = 0.1,linewidth = 0.05, color=c)
        
        ax4.plot(wl[sel_wl2], matrix_spec[i_model, sel_wl2]/_norm2, alpha = 0.1,linewidth = 0.05, color=c)
        ax5.plot(wl[sel_wl2], (matrix_spec[i_model, sel_wl2]/_norm2)/(matrix_spec[arg_delta_min, sel_wl2]/norm2), alpha = 0.1,linewidth = 0.05, color=c)
        
    
    NORM=Normalize(vmin=dage_min, vmax=dage_max)
    clb.ColorbarBase(ax3, cmap=cm.rainbow_r, norm=NORM)
    ax3.plot([-1,1], [time_res5, time_res5], linewidth=1, color='black') 
    ax3.plot([-1,1], [time_res10, time_res10], linewidth=1, color='black') 
    ax3.plot([-1,1], [time_res20, time_res20], linewidth=1, color='black') 
    ax3.plot([-1,1], [time_res50, time_res50], linewidth=1, color='black') 
    ax3.plot([-1,1], [time_res100, time_res100], linewidth=1, color='black') 
    ax3.plot([-1,1], [time_res200, time_res200], linewidth=1, color='black') 
    ax3.plot([-1,1], [time_res500, time_res500], linewidth=1, color='black') 
    
    ax1.set_facecolor('#d8dcd6')
    ax2.set_facecolor('#d8dcd6')
    ax4.set_facecolor('#d8dcd6')
    ax5.set_facecolor('#d8dcd6')
    ax1.tick_params(labelsize=15)
    ax2.tick_params(labelsize=15)
    ax4.tick_params(labelsize=15)
    ax5.tick_params(labelsize=15)
    ax3.tick_params(labelsize=10)
    ax2.set_ylim([0.95,1.05])
    ax5.set_ylim([0.95,1.05])

    #axs[1,0].plot([np.min(wl[sel_wl]), np.max(wl[sel_wl])], np.array([1.0,1.0])/SNR)
    #axs[1,0].plot([np.min(wl[sel_wl]), np.max(wl[sel_wl])], -np.array([1.0,1.0])/SNR)
    ax1.minorticks_on()
    ax2.minorticks_on()
    ax4.minorticks_on()
    ax5.minorticks_on()
    
    ax2.set_xlabel(r'$\lambda [\AA]$', size=20)
    ax5.set_xlabel(r'$\lambda [\AA]$', size=20)

    ax2.set_ylabel(r'$L_{\lambda}/L_{\lambda,ref}$', size=20)
    ax1.set_ylabel(r'$L_{\lambda}/L_{\lambda,cont}$', size=20)
    ax3.set_ylabel(r'$log_{10}(\Delta age_{n})$', size=20)
    
    ax1.set_title(title1)
    ax4.set_title(title2)
    
    #arrow
    ax3.annotate('S/N=5',
    xy=(-1, time_res5),
    xycoords='data',
    xytext=(-3, time_res5),
    arrowprops=
    dict(facecolor='darkred', shrink=0.05),
    horizontalalignment='right',
    verticalalignment='center', size=15)
    
    ax3.annotate('',
    xy=(-1, time_res10),
    xycoords='data',
    xytext=(-3, time_res10),
    arrowprops=
    dict(facecolor='red', shrink=0.05),
    horizontalalignment='right',
    verticalalignment='center')
    
    ax3.annotate('',
    xy=(-1, time_res20),
    xycoords='data',
    xytext=(-3, time_res20),
    arrowprops=
    dict(facecolor='orange', shrink=0.05),
    horizontalalignment='right',
    verticalalignment='center')
    
    ax3.annotate('S/N=50',
    xy=(-1, time_res50),
    xycoords='data',
    xytext=(-3, time_res50),
    arrowprops=
    dict(facecolor='yellow', shrink=0.05),
    horizontalalignment='right',
    verticalalignment='center', size=15)
    
    ax3.annotate('',
    xy=(-1, time_res100),
    xycoords='data',
    xytext=(-3, time_res100),
    arrowprops=
    dict(facecolor='green', shrink=0.05),
    horizontalalignment='right',
    verticalalignment='center')
    
    ax3.annotate('',
    xy=(-1, time_res200),
    xycoords='data',
    xytext=(-3, time_res200),
    arrowprops=
    dict(facecolor='blue', shrink=0.05),
    horizontalalignment='right',
    verticalalignment='center')
    
    ax3.annotate('S/N=500',
    xy=(-1, time_res500),
    xycoords='data',
    xytext=(-3, time_res500),
    arrowprops=
    dict(facecolor='darkblue', shrink=0.05),
    horizontalalignment='right',
    verticalalignment='center', size=15)
     
    return fig

    





