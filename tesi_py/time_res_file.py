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

work_dir=os.getenv('HOME')+'/Desktop/TESI/models/Sandage_v4.1_Zfix_noburst_bc03MILES_100k/'
z_list=['32', '42', '52', '62', '72']
N_z=np.size(z_list)
n_chunks=5
N_bins=30

t=Table()


for z in range(0, N_z):
    idx_file=work_dir+'sandage_varZ_v4.1_m'+z_list[z]+'fix_noburst_100k_spec_dcombnull_idx_001.fits'
    par_file=work_dir+'sandage_varZ_v4.1_m'+z_list[z]+'fix_noburst_100k_spec_dcombnull_001_physpar_wagef.fits'
        
    hdul_idx=fits.open(idx_file)
    hdul_par=fits.open(par_file)
    
    d4000n=hdul_idx[1].data['D4000N'][...,4]
    hdhg=hdul_idx[1].data['HdHg'][...,4]
    Hb=hdul_idx[1].data['Lick_Hb'][...,4]
    mg2fe=hdul_idx[1].data['Mg2Fe'][...,4]
    mgfep=hdul_idx[1].data['MgFe_prime'][...,4]
    
    mag_u=hdul_par[1].data['ABMAG'][...,0]
    mag_g=hdul_par[1].data['ABMAG'][...,1]
    mag_r=hdul_par[1].data['ABMAG'][...,2]
    mag_i=hdul_par[1].data['ABMAG'][...,3]
    mag_z=hdul_par[1].data['ABMAG'][...,4]
    
    age10=hdul_par[1].data['age10']
    age25=hdul_par[1].data['age25']
    age50=hdul_par[1].data['age50']
    age75=hdul_par[1].data['age75']
    age90=hdul_par[1].data['age90']

    
    
    for i_chunks in range(2, n_chunks+1):
        _idx_file=work_dir+'sandage_varZ_v4.1_m'+z_list[z]+'fix_noburst_100k_spec_dcombnull_idx_{:03d}.fits'
        _par_file=work_dir+'sandage_varZ_v4.1_m'+z_list[z]+'fix_noburst_100k_spec_dcombnull_{:03d}_physpar_wagef.fits'
        
        _hdul_idx=fits.open(_idx_file.format(i_chunks))
        _hdul_par=fits.open(_par_file.format(i_chunks))
    
        _d4000n=_hdul_idx[1].data['D4000N'][...,4]
        _hdhg=_hdul_idx[1].data['HdHg'][...,4]
        _Hb=_hdul_idx[1].data['Lick_Hb'][...,4]
        _mg2fe=_hdul_idx[1].data['Mg2Fe'][...,4]
        _mgfep=_hdul_idx[1].data['MgFe_prime'][...,4]
    
        _mag_u=_hdul_par[1].data['ABMAG'][...,0]
        _mag_g=_hdul_par[1].data['ABMAG'][...,1]
        _mag_r=_hdul_par[1].data['ABMAG'][...,2]
        _mag_i=_hdul_par[1].data['ABMAG'][...,3]
        _mag_z=_hdul_par[1].data['ABMAG'][...,4]
    
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
        
        age10=np.append(age10, _age10)
        age25=np.append(age25, _age25)
        age50=np.append(age50, _age50)
        age75=np.append(age75, _age75)
        age90=np.append(age90, _age90)
        
    
    d1090n50=np.log10((age10-age90)/age50)
    bin_age50=np.histogram(np.log10(age50), bins=N_bins, range=(8.65, 10.15))[1]
    t_res=[0.0]*N_bins
    
    for i in range(0, N_bins):
        idx_ref=((d1090n50<-1.00)&(np.log10(age50)<bin_age50[i+1])& (np.log10(age50)>bin_age50[i]))
        idx_sel=((np.log10(age50)<bin_age50[i+1])& (np.log10(age50)>bin_age50[i]))
        t_res[i]=f_plt.chi_q(d1090n50,d4000n,hdhg,Hb,mg2fe,mgfep,mag_u-mag_r,mag_g-mag_r, mag_r-mag_i, mag_r-mag_z, idx_sel,idx_ref, toll=0.0001, xmin=-1.50, mkplot=False)
    
    name_col='Log_d1090n50_min_z'+z_list[z]
    new_col=Column(t_res, name=name_col)
    t.add_column(new_col)


new_file=work_dir+'Time_resol_Zfix.fits'

t.write(new_file, format='fits', overwrite=True)
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    