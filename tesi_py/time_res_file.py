#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 13 18:38:34 2021

@author: edoardo
"""

import os
import numpy as np
import function_plot as f_plt
import astropy.io.fits as fits
from astropy.table import Table
from astropy.table import Column

work_dir=os.getenv('HOME')+'/Desktop/TESI/models/Sandage_v4.1_Zfix_noburst_cb16MILES_1M/'
z_list=['32', '42', '52', '62', '72']
N_z=np.size(z_list)
n_chunks=50
N_bins=83
snr='500'

t=Table()


for z in range(0, N_z):
    #idx_file=work_dir+'sandage_varZ_v4.1_m'+z_list[z]+'fix_noburst_100k_spec_dcombnull_idx_001.fits'
    par_file=work_dir+'sandage_varZ_v4.1_m'+z_list[z]+'fix_noburst_1M_spec_dcombnull_001_physpar_wagef.fits'
    file_pert=work_dir+'sandage_varZ_v4.1_m'+z_list[z]+'fix_noburst_1M_spec_dcombnull_perterr_SNR'+snr+'_001.fits'
    
    
    #hdul_idx=fits.open(idx_file)
    hdul_par=fits.open(par_file)
    hdul_pert=fits.open(file_pert)

    
    d4000n=hdul_pert[1].data['D4000N_SIGMA200_PERT']
    hdhg=hdul_pert[1].data['HDHG_SIGMA200_PERT']
    Hb=hdul_pert[1].data['Lick_Hb_SIGMA200_PERT']
    mg2fe=hdul_pert[1].data['Mg2Fe_SIGMA200_PERT']
    mgfep=hdul_pert[1].data['MgFe_prime_SIGMA200_PERT']

    
    mag_u=hdul_pert[1].data['ABMAG_U_PERT']
    mag_g=hdul_pert[1].data['ABMAG_G_PERT']
    mag_r=hdul_pert[1].data['ABMAG_R_PERT']
    mag_i=hdul_pert[1].data['ABMAG_I_PERT']
    mag_z=hdul_pert[1].data['ABMAG_Z_PERT']
    
    sigma_u=hdul_pert[1].data['ERR_MAG_U']
    sigma_g=hdul_pert[1].data['ERR_MAG_G']
    sigma_r=hdul_pert[1].data['ERR_MAG_R']
    sigma_i=hdul_pert[1].data['ERR_MAG_I']
    sigma_z=hdul_pert[1].data['ERR_MAG_Z']

    age10=hdul_par[1].data['age10']
    age25=hdul_par[1].data['age25']
    age50=hdul_par[1].data['age50']
    age75=hdul_par[1].data['age75']
    age90=hdul_par[1].data['age90']

    
    
    for i_chunks in range(2, n_chunks+1):
       # _idx_file=work_dir+'sandage_varZ_v4.1_m'+z_list[z]+'fix_noburst_100k_spec_dcombnull_idx_{:03d}.fits'
        _par_file=work_dir+'sandage_varZ_v4.1_m'+z_list[z]+'fix_noburst_1M_spec_dcombnull_{:03d}_physpar_wagef.fits'
        _pert_file=work_dir+'sandage_varZ_v4.1_m'+z_list[z]+'fix_noburst_1M_spec_dcombnull_perterr_SNR'+snr+'_{:03d}.fits'

        
        #_hdul_idx=fits.open(_idx_file.format(i_chunks))
        _hdul_par=fits.open(_par_file.format(i_chunks))
        _hdul_pert=fits.open(_pert_file.format(i_chunks))
   
        _d4000n=_hdul_pert[1].data['D4000N_SIGMA200_PERT']
        _hdhg=_hdul_pert[1].data['HDHG_SIGMA200_PERT']
        _Hb=_hdul_pert[1].data['Lick_Hb_SIGMA200_PERT']
        _mg2fe=_hdul_pert[1].data['Mg2Fe_SIGMA200_PERT']
        _mgfep=_hdul_pert[1].data['MgFe_prime_SIGMA200_PERT']
    
        _mag_u=_hdul_pert[1].data['ABMAG_U_PERT']
        _mag_g=_hdul_pert[1].data['ABMAG_G_PERT']
        _mag_r=_hdul_pert[1].data['ABMAG_R_PERT']
        _mag_i=_hdul_pert[1].data['ABMAG_I_PERT']
        _mag_z=_hdul_pert[1].data['ABMAG_Z_PERT']
        
        _sigma_u=_hdul_pert[1].data['ERR_MAG_U']
        _sigma_g=_hdul_pert[1].data['ERR_MAG_G']
        _sigma_r=_hdul_pert[1].data['ERR_MAG_R']
        _sigma_i=_hdul_pert[1].data['ERR_MAG_I']
        _sigma_z=_hdul_pert[1].data['ERR_MAG_Z']
    
        _age10=_hdul_par[1].data['age10']
        _age25=_hdul_par[1].data['age25']
        _age50=_hdul_par[1].data['age50']
        _age75=_hdul_par[1].data['age75']
        _age90=_hdul_par[1].data['age90']
        
        d4000n=np.append(d4000n, _d4000n)
        hdhg=np.append(hdhg, _hdhg)
        Hb=np.append(Hb, _Hb)
        mg2fe=np.append(mg2fe, _mg2fe)
        mgfep=np.append(mgfep, _mgfep)
        
        mag_u=np.append(mag_u, _mag_u)
        mag_g=np.append(mag_g, _mag_g)
        mag_r=np.append(mag_r, _mag_r)
        mag_i=np.append(mag_i, _mag_i)
        mag_z=np.append(mag_z, _mag_z)
        
        sigma_u=np.append(sigma_u, _sigma_u)
        sigma_g=np.append(sigma_g, _sigma_g)
        sigma_r=np.append(sigma_r, _sigma_r)
        sigma_i=np.append(sigma_i, _sigma_i)
        sigma_z=np.append(sigma_z, _sigma_z)
        
        age10=np.append(age10, _age10)
        age25=np.append(age25, _age25)
        age50=np.append(age50, _age50)
        age75=np.append(age75, _age75)
        age90=np.append(age90, _age90)
        
    sigma6=np.sqrt(sigma_u**2+sigma_r**2)
    sigma7=np.sqrt(sigma_g**2+sigma_r**2)
    sigma8=np.sqrt(sigma_r**2+sigma_i**2)
    sigma9=np.sqrt(sigma_r**2+sigma_z**2)
    
    d1090n50=np.log10((age10-age90)/age50)
    bin_age50=np.histogram(np.log10(age50), bins=N_bins, range=(6.0, 10.15))[1]
    t_res=[0.0]*N_bins
    
    for i in range(0, N_bins):
        idx_ref=((d1090n50<-1.00)&(np.log10(age50)<bin_age50[i+1])& (np.log10(age50)>bin_age50[i]))
        idx_sel=((np.log10(age50)<bin_age50[i+1])& (np.log10(age50)>bin_age50[i]))
        t_res[i]=f_plt.chi_q(d1090n50,d4000n,hdhg,Hb,mg2fe,mgfep,mag_u-mag_r,mag_g-mag_r, mag_r-mag_i, mag_r-mag_z, idx_sel,idx_ref, toll=0.001,  mkplot=False)
        
   
    name_col='Log_d1090n50_min_z'+z_list[z]
    new_col=Column(t_res, name=name_col)
    t.add_column(new_col)

new_file=work_dir+'Time_resol_Zfix1M_SNR'+snr+'_col.fits'

t.write(new_file, format='fits', overwrite=True)
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    