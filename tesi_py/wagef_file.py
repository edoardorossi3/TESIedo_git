#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 20 12:03:53 2020

@author: edoardo
"""


import astropy.io.fits as fits
#import matplotlib.pyplot as plt
import numpy as np

#hdul_th=fits.open('tesi_par.fits')  #t_f parameters
work_dir='/home/edoardo/Desktop/TESI/models/Sandage_v4.1_Zfix_noburst_bc03MILES_100k/'
z_list=['32', '42', '52', '62', '72']
Nz=np.size(z_list)
N_files=5
N_data=20000   #number of data (rows) in each chunk

for j in range(0, Nz):
    t_f_file=work_dir+'t_f_m'+z_list[j]+'fix_noburst.fits'
    f_name=work_dir+'sandage_varZ_v4.1_m'+z_list[j]+'fix_noburst_100k.fits'   #name of the file from that we get tform, tau, nbursts, ecc...
    
    
    hdul_th=fits.open(t_f_file)
    hdul2=fits.open(f_name)   
    
    log_tobs=hdul2[1].data['LOGTFORM']
    
    t10=hdul_th[1].data['t10']
    age10=np.power(10,log_tobs)-t10
    t90=hdul_th[1].data['t90']
    age90=np.power(10,log_tobs)-t90
    t50=hdul_th[1].data['t50']
    age50=np.power(10, log_tobs)-t50
    t25=hdul_th[1].data['t25']
    age25=np.power(10,log_tobs)-t25
    t75=hdul_th[1].data['t75']
    age75=np.power(10,log_tobs)-t75
    
    
    for i in range (1, N_files+1):
        c1s=fits.Column(name='age10', array=age10[N_data*(i-1):N_data*i], format='E', unit='yr')
        c2s=fits.Column(name='age25', array=age25[N_data*(i-1):N_data*i], format='E',unit='yr')
        c3s=fits.Column(name='age50', array=age50[N_data*(i-1):N_data*i], format='E',unit='yr')
        c4s=fits.Column(name='age75', array=age75[N_data*(i-1):N_data*i], format='E',unit='yr')
        c5s=fits.Column(name='age90', array=age90[N_data*(i-1):N_data*i], format='E',unit='yr')
        
        #sndg_file=work_dir+'sandage_varZ_v4.1eq_spec_dcomb_{:03d}_physpar.fits'
        #hdul=fits.open('par_age_dcomb001.fits')
        sndg_file=work_dir+'sandage_varZ_v4.1_m'+z_list[j]+'fix_noburst_100k_spec_dcombnull_{:03d}_physpar.fits'
        hdul1=fits.open(sndg_file.format(i))
        #hdul1=fits.open(sndg_file)
        
        hdr1=hdul1[1].header
        hdr1.append(('TCOMM18', 'Lookback time at the epoch of 10% of stars formed'))
        hdr1.append(('TCOMM19', 'Lookback time at the epoch of 25% of stars formed'))
        hdr1.append( ('TCOMM20', 'Lookback time at the epoch of 50% of stars formed'))
        hdr1.append( ('TCOMM21', 'Lookback time at the epoch of 75% of stars formed'))
        hdr1.append( ('TCOMM22', 'Lookback time at the epoch of 90% of stars formed'))
        
        col_tot=hdul1[1].columns + c1s+c2s+c3s+c4s+c5s
        
        #new_file='/home/edoardo/Desktop/TESI/models/Sandage_varZ_v4.1eq_bc03MILES_ChFall/sandage_varZ_v4.1eq_spec_dcomb_{:03d}_physpar_wagef.fits'
        new_file=work_dir+'sandage_varZ_v4.1_m'+z_list[j]+'fix_noburst_spec_dcombnull_{:03d}_physpar_wagef.fits'
        hs=fits.BinTableHDU.from_columns(col_tot, header=hdul1[1].header)
        hs.writeto(new_file.format(i), overwrite=True)
        #hs.writeto(new_file, overwrite=True)
        
        #hdul.close()
        hdul1.close()
    
    
    hdul_th.close()
    hdul2.close()
