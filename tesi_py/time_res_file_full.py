#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 20 17:06:54 2022

@author: edoardo.rossi
"""
#import libraries
import numpy as np
import astropy.io.fits as fits
import scipy.ndimage as scind
import function_plot as f_plt
from tqdm import tqdm
import astropy.io.fits as fits
from astropy.table import Table
from astropy.table import Column



#files and directory
work_dir='/export/home/extragal/zibetti/no_ownCloud/SteMaGE/data/SEDlibraries/Sandage_v4.1_Zfix_noburst_cb16MILES_1M/'
f_name=work_dir+'sandage_varZ_v4.1_m62fix_noburst_1M_spec_dcombnull_001.fits'
f_par=work_dir+'sandage_varZ_v4.1_m62fix_noburst_1M_spec_dcombnull_001_physpar_wagef.fits'


t=Table()

#age 
prefix_file='sandage_varZ_v4.1_m62fix_noburst_1M_spec_dcombnull'
_file_par=work_dir+prefix_file+'_{:03d}_physpar_wagef.fits'
_hdul_par=fits.open(_file_par.format(1))
hdul=fits.open(f_name)

np.random.seed(0)

age10=_hdul_par[1].data['age10']
age50=_hdul_par[1].data['age50']
age90=_hdul_par[1].data['age90']
for i_chunk in range (2,51):

    _hdul_par=fits.open(_file_par.format(i_chunk))

    _age10=_hdul_par[1].data['age10']
    _age50=_hdul_par[1].data['age50']
    _age90=_hdul_par[1].data['age90']
    age50=np.append(age50, _age50)
    age10=np.append(age10, _age10)
    age90=np.append(age90, _age90)
    
dage_n=np.log10((age10-age90)/age50)
age10=np.log10(age10)
age90=np.log10(age90)
age50=np.log10(age50)

#bin age50
age50_bin=np.histogram(age50, bins=83, range=(6.0, 10.15))[1]
Nbin_age50=np.size(age50_bin)
dage_n_lim=np.zeros(Nbin_age50)
#analysis for each bin
wl=hdul[0].data
sel_wl=(wl<5600) & (wl>3800)
wl_sel=wl[sel_wl]
N_col=np.size(wl_sel)

snr=20
wl_0=5500
wl_1=5550
sel_wl_sub=((wl_sel<5550) & (wl_sel>5500))
wl_sub=wl_sel[sel_wl_sub]
wl_mid=np.zeros(np.size(wl_sub)-1)
dwl=np.zeros(np.size(wl_sub)-1)
for i in range(0, np.size(wl_sub)-1):
    wl_mid[i]=(wl_sub[i]+wl_sub[i+1])/2
    dwl=wl_sub[i+1]-wl_sub[i]


for i_bin in tqdm(range(1,Nbin_age50)):
    
    sel_models=(((age50)<age50_bin[i_bin]) & ((age50)>age50_bin[i_bin-1]))
    N_row=np.sum(sel_models)
    matrix_spec=np.zeros(shape=(N_row, N_col))
    first_row=0
    delta_spec=np.zeros(N_row)

    
    for i_chunk in range(0,50):
        _file_spec=work_dir+prefix_file+'_{:03d}.fits'
        _hdul_spec=fits.open(_file_spec.format(i_chunk+1))
        _sel_mod=sel_models[i_chunk*20000:(i_chunk+1)*20000]
        _y_data=_hdul_spec[1].data
        _y_data=_y_data[..., sel_wl]
        last_row=np.sum(_sel_mod)+first_row
        
        matrix_spec[first_row:last_row,...]=_y_data[_sel_mod,...]
        first_row=last_row
    
    dage_n_sub=dage_n[sel_models]
    idx_ref=np.argmin(dage_n_sub)

    y_data_interp=np.interp(wl_mid, wl_sub, matrix_spec[idx_ref, sel_wl_sub])
    err_spec=np.sum((y_data_interp*dwl)/((wl_1-wl_0)*snr))
    err_spec_norm=err_spec/scind.median_filter(matrix_spec[idx_ref, ...], size=500)
    _pert=np.random.normal(scale=err_spec_norm, size=np.size((matrix_spec[idx_ref, ...]/scind.median_filter(matrix_spec[idx_ref, ...], size=500))))
    
    spec_ref_pert=(matrix_spec[idx_ref, ...]/scind.median_filter(matrix_spec[idx_ref, ...], size=500))+_pert
    sel_ref=(dage_n_sub<-1.0)
    for i_mod in range(0, N_row):
        y_data_interp=np.interp(wl_mid, wl_sub, matrix_spec[i_mod, sel_wl_sub])
        err_spec=np.sum((y_data_interp*dwl)/((wl_1-wl_0)*snr))
        err_spec_norm=err_spec/scind.median_filter(matrix_spec[i_mod, ...], size=500)
        _pert=np.random.normal(scale=err_spec_norm, size=np.size((matrix_spec[i_mod, ...]/scind.median_filter(matrix_spec[idx_ref, ...], size=500))))
        spec_mod_pert=(matrix_spec[i_mod, ...]/scind.median_filter(matrix_spec[i_mod, ...], size=500))+_pert
        for i_wl in range(0, N_col):
            #_std_ref=np.std((matrix_spec[sel_ref, i_wl]/scind.median_filter(matrix_spec[sel_ref, i_wl], size=500))+_pert[sel_ref])
            delta_spec[i_mod]=delta_spec[i_mod]+((spec_mod_pert[i_wl]-spec_ref_pert[i_wl])/err_spec_norm[i_wl])**2
    
    
    sel_ref=(dage_n_sub<-1.0)
    std_ref=np.std(delta_spec[sel_ref])
    mean_ref=np.mean(delta_spec[sel_ref])
    perc_84_ref=np.percentile(delta_spec[sel_ref], 84)
    
    bins=50
    bin_edges=np.histogram(dage_n_sub, bins=bins)[1]
    bin_centre=np.zeros(bins)
    perc_delta_84=np.zeros(bins)
    std=np.zeros(bins)
    perc_delta_16=np.zeros(bins)
    
    for i in range(0,bins):
        bin_centre[i]=(bin_edges[i]+bin_edges[i+1])/2
        _sel=((dage_n_sub<bin_edges[i+1]) & (dage_n_sub>bin_edges[i]))
        perc_delta_84[i]=np.percentile(delta_spec[_sel], 84)
        perc_delta_16[i]=np.percentile(delta_spec[_sel], 16)
        
    inter_delta=lambda x: np.interp(x, bin_centre, perc_delta_16-(mean_ref+std_ref))
    dage_n_lim[i_bin]=f_plt.bisection(inter_delta, -1.0, 1.2, 0.01)

name_col="dage_n_lim_m62"
new_col=Column(dage_n_lim, name=name_col)
t.add_column(new_col)
t.write(work_dir+'Time_resol_full_Zfix1M_SNR'+str(snr)+'.fits',format='fits', overwrite=True)
