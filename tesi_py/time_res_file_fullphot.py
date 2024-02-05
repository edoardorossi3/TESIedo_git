#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 20 17:06:54 2022

@author: stefano.zibetti from edoardo.rossi
"""
#%% import libraries

import numpy as np
import astropy.io.fits as fits
import scipy.ndimage as scind
import function_plot as f_plt
from tqdm import tqdm
import astropy.io.fits as fits
from astropy.table import Table
from astropy.table import Column



# working parameters
work_dir='/export/home/extragal/zibetti/ownCloud_Arcetri/SteMaGE/data/SEDlibraries/Sandage_v4.1_Zfix_noburst_cb16MILES_1M/'

## spectral window
wl_min=3800. 
wl_max=5600.
# wl_min=3500.
# wl_max=7000.

prefix_file='sandage_varZ_v4.1_m62fix_noburst_1M_spec_dcombnull'
name_col="dage_n_lim_m62"

snr=100
out_fname='Time_resol_fullphot_Z62fix1M_SNR'+str(snr)+'.fits'

## normalization window
wl_0=5500
wl_1=5550

t=Table()
#%%

# read in wl array
hdul=fits.open(work_dir+prefix_file+'_001.fits')
wl=hdul[0].data
sel_wl=(wl<wl_max) & (wl>wl_min) 
wl_sel=wl[sel_wl]

# read in and define ages
age10=np.zeros(0)
age50=np.zeros(0)
age90=np.zeros(0)
pert_ur=np.zeros(0)
pert_gr=np.zeros(0)
pert_ri=np.zeros(0)
pert_rz=np.zeros(0)

_file_par=work_dir+prefix_file+'_{:03d}_physpar_wagef.fits'
_file_pert=work_dir+prefix_file+'_perterr_SNR'+str(snr)+'_{:03d}.fits'

for i_chunk in range (1,51):
    # fill in agef arrays
    _hdul_par=fits.open(_file_par.format(i_chunk))
    # _age10=_hdul_par[1].data['age10']
    # _age50=_hdul_par[1].data['age50']
    # _age90=_hdul_par[1].data['age90']
    age50=np.append(age50, _hdul_par[1].data['age50'])
    age10=np.append(age10, _hdul_par[1].data['age10'])
    age90=np.append(age90, _hdul_par[1].data['age90'])
    _hdul_par.close()
    
    # fill in color arrays
    _hdul_pert=fits.open(_file_pert.format(i_chunk))
    _mag_u=_hdul_pert[1].data['ABMAG_U_PERT']
    _mag_g=_hdul_pert[1].data['ABMAG_G_PERT']
    _mag_r=_hdul_pert[1].data['ABMAG_R_PERT']
    _mag_i=_hdul_pert[1].data['ABMAG_I_PERT']
    _mag_z=_hdul_pert[1].data['ABMAG_Z_PERT']
    pert_ur=np.append(pert_ur, _mag_u-_mag_r)
    pert_gr=np.append(pert_gr, _mag_g-_mag_r)
    pert_ri=np.append(pert_ri, _mag_r-_mag_i)
    pert_rz=np.append(pert_rz, _mag_r-_mag_z)
    _hdul_pert.close()

dage_n=np.log10((age10-age90)/age50)
lgage50=np.log10(age50)

#bin age50
lgage50_bin=np.histogram(lgage50, bins=83, range=(6.0, 10.15))[1]
Nbin_lgage50=np.size(lgage50_bin)
dage_n_lim=np.zeros(Nbin_lgage50-1)
#analysis for each bin
N_wl_sel=np.size(wl_sel)

sel_wl_sub=((wl_sel<wl_1) & (wl_sel>wl_0))
wl_sub=wl_sel[sel_wl_sub]
wl_mid=np.zeros(np.size(wl_sub)-1)
dwl=np.zeros(np.size(wl_sub)-1)
for i in range(0, np.size(wl_sub)-1):
    wl_mid[i]=(wl_sub[i]+wl_sub[i+1])/2
    dwl=wl_sub[i+1]-wl_sub[i]

#%%
np.random.seed(0)

for i_bin in tqdm(range(0,Nbin_lgage50-1)): #tqdm(range(0,3)): #
    
    sel_models=(((lgage50)<lgage50_bin[i_bin+1]) & ((lgage50)>lgage50_bin[i_bin]))
    dage_n_sub=dage_n[sel_models]
    idx_ref=np.argmin(dage_n_sub)
    sel_ref=(dage_n_sub<-1.0) #reference range
    
    # spectra
    N_row=np.sum(sel_models)
    matrix_spec=np.zeros(shape=(N_row, N_wl_sel))
    delta_spec=np.zeros(N_row)

    
    first_row=0
    for i_chunk in range(0,50):
        _file_spec=work_dir+prefix_file+'_{:03d}.fits'
        _hdul_spec=fits.open(_file_spec.format(i_chunk+1))
        _sel_mod=sel_models[i_chunk*20000:(i_chunk+1)*20000]
        _y_data=_hdul_spec[1].data
        _y_data=_y_data[..., sel_wl]
        last_row=np.sum(_sel_mod)+first_row
        
        matrix_spec[first_row:last_row,...]=_y_data[_sel_mod,...]
        first_row=last_row
        _hdul_spec.close()
    


    # smooth ref spectrum
    spec_ref_smooth=scind.median_filter(matrix_spec[idx_ref, ...], size=500)
    # error [constant] = normalization factor / snr
    err_spec=np.trapz(matrix_spec[idx_ref, sel_wl_sub],x=wl_sub)/(wl_1-wl_0)/snr
    _pert=np.random.normal(scale=err_spec, size=N_wl_sel)
    spec_ref_pert=(matrix_spec[idx_ref, ...]+_pert)/spec_ref_smooth
                                   
    for i_mod in range(0, N_row):
        # smooth  spectrum
        spec_smooth=scind.median_filter(matrix_spec[i_mod, ...], size=500)
        # error [constant] = normalization factor / snr
        err_spec=np.trapz(matrix_spec[i_mod, sel_wl_sub],x=wl_sub)/(wl_1-wl_0)/snr
        err_spec_norm=err_spec/spec_smooth
        _pert=np.random.normal(scale=err_spec, size=N_wl_sel)
        spec_mod_pert=(matrix_spec[i_mod, ...]+_pert)/spec_smooth
        delta_spec[i_mod]=np.sum(((spec_mod_pert-spec_ref_pert)/err_spec_norm)**2)
        

    # colors
    pert_ur_sub=pert_ur[sel_models]
    pert_gr_sub=pert_gr[sel_models]
    pert_ri_sub=pert_ri[sel_models]
    pert_rz_sub=pert_rz[sel_models]
    # define reference variance
    var_ur=np.var(pert_ur_sub[sel_ref]-pert_ur_sub[idx_ref])
    var_gr=np.var(pert_gr_sub[sel_ref]-pert_gr_sub[idx_ref])
    var_ri=np.var(pert_ri_sub[sel_ref]-pert_ri_sub[idx_ref])
    var_rz=np.var(pert_rz_sub[sel_ref]-pert_rz_sub[idx_ref])

    delta_col=(pert_ur_sub-pert_ur_sub[idx_ref])**2/var_ur+\
        (pert_gr_sub-pert_gr_sub[idx_ref])**2/var_gr+\
        (pert_ri_sub-pert_ri_sub[idx_ref])**2/var_ri+\
        (pert_rz_sub-pert_rz_sub[idx_ref])**2/var_rz
                    
    
    delta2=delta_spec+delta_col
    
    std_ref=np.std(delta2[sel_ref])
    mean_ref=np.mean(delta2[sel_ref])
    perc_84_ref=np.percentile(delta2[sel_ref], 84)
    
    #print(np.mean(delta2[sel_ref]),np.mean(delta_col[sel_ref]),
    #      np.mean(delta_spec[sel_ref]))
    
    bins=50
    bin_edges=np.histogram(dage_n_sub, bins=bins)[1]
    bin_centre=np.zeros(bins)
    perc_delta_84=np.zeros(bins)
    perc_delta_16=np.zeros(bins)
    std=np.zeros(bins)
    
    for i in range(0,bins):
        bin_centre[i]=(bin_edges[i]+bin_edges[i+1])/2
        _sel=((dage_n_sub<bin_edges[i+1]) & (dage_n_sub>bin_edges[i]))
        perc_delta_84[i]=np.percentile(delta2[_sel], 84)
        perc_delta_16[i]=np.percentile(delta2[_sel], 16)
        
    inter_delta=lambda x: np.interp(x, bin_centre, perc_delta_16-(mean_ref+std_ref))
    dage_n_lim[i_bin]=f_plt.bisection(inter_delta, -1.0, 1.2, 0.01)

#name_col="dage_n_lim_m42"
new_col=Column(dage_n_lim, name=name_col)
t.add_column(new_col)
t.write(work_dir+out_fname,format='fits', overwrite=True)
