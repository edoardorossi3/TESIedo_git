#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 12 20:40:46 2021

@author: edoardo
"""

#work_dir: directory in which we want to create the wagef files
#f_name: the file's name from which we get t_obs, tau...
#suffix_physpar: suffix of the physpar file (up to and including dcomb)
#n_chunks: number of physpar files



import numpy as np
import astropy.io.fits as fits
from tqdm import tqdm

work_dir='/home/edoardo/Desktop/TESI/models/Sandage_varZ_v4.2eq_CB16MILES_ChFall/'
f_name='sandage_varZ_v4.2eq.fits'
suffix_physpar='sandage_varZ_v4.2eq_spec_dcomb0p25null'
n_chunks=40

def bisection(f, t_min, t_max, toll):
    if f(t_min)*f(t_max)>=0:
        return np.nan
    while (t_max-t_min)>toll:
        t_m=(t_max+t_min)/2.0
        if f(t_min)*f(t_m)>0:
            t_min=t_m
            t_max=t_max
        elif f(t_max)*f(t_m)>0:
            t_max=t_m
            t_min=t_min
        elif f(t_m)==0:
            return t_m
    return (t_max+t_min)/2.0

def Massb(t, tau, t_obs, A_burst, t_burst, N_burst):
    A=1.0/(1.0-np.exp(-0.5*(t_obs/tau)**2))
    m_sndg=A*(1-np.exp(-0.5*(t/tau)**2))
    for i in range(0, N_burst): 
        if t>=t_burst[i]: m_sndg=m_sndg+A_burst[i]
    return m_sndg


f_name=work_dir+f_name
N_files=n_chunks

hdul=fits.open(f_name)
data=hdul[1].data
N_models=np.size(data)
t10b=[0.0]*N_models
t25b=[0.0]*N_models
t50b=[0.0]*N_models
t75b=[0.0]*N_models
t90b=[0.0]*N_models
log_tobs=hdul[1].data['LOGTFORM']
t_obs=10.0**(log_tobs)

for k in tqdm(range(0, N_models)):
    #t_obs=10.0**(data[k][0]) #t_obs in yr
    tau=10.0**(data[k][1]) #tau in yr
    Nburst=data[k][5]
    A_burst=data[k][7]
    t_burst=[0.0]*6
    age_burst=[0.0]*6
    age_burst=np.where(data[k][6]!=0, 10**(data[k][6]), age_burst) #age_burst in yr
    t_burst=np.where(age_burst==0, t_burst, t_obs[k]-age_burst)
     
                 
    #Mass_tot=Massb(t_obs, tau, t_obs, A_burst, t_burst, Nburst)                
    Mass_tot=1.0+np.sum(A_burst)    
                 
    frac50=lambda x: Massb(x, tau, t_obs[k], A_burst, t_burst, Nburst)/Mass_tot-0.5
    frac25=lambda x: Massb(x, tau, t_obs[k], A_burst, t_burst, Nburst)/Mass_tot-0.25
    frac75=lambda x: Massb(x, tau, t_obs[k], A_burst, t_burst, Nburst)/Mass_tot-0.75
    frac90=lambda x: Massb(x, tau, t_obs[k], A_burst, t_burst, Nburst)/Mass_tot-0.9
    frac10=lambda x: Massb(x, tau, t_obs[k], A_burst, t_burst, Nburst)/Mass_tot-0.1
    
    t50b[k]=bisection(frac50, 0.0, t_obs[k], 0.01)
    t25b[k]=bisection(frac25, 0.0, t_obs[k], 0.01)
    t75b[k]=bisection(frac75, 0.0, t_obs[k], 0.01)
    t90b[k]=bisection(frac90, 0.0, t_obs[k], 0.01)
    t10b[k]=bisection(frac10, 0.0, t_obs[k], 0.01)
    
age10=t_obs-t10b
age90=t_obs-t90b
age50=t_obs-t50b
age25=t_obs-t25b
age75=t_obs-t75b

for i in range (1, N_files+1):
    sndg_file=work_dir+suffix_physpar+'_{:03d}_physpar.fits'
    hdul1=fits.open(sndg_file.format(i))
    
    N_data=np.size(hdul1[1].data)
    
    c1s=fits.Column(name='age10', array=age10[N_data*(i-1):N_data*i], format='E', unit='yr')
    c2s=fits.Column(name='age25', array=age25[N_data*(i-1):N_data*i], format='E',unit='yr')
    c3s=fits.Column(name='age50', array=age50[N_data*(i-1):N_data*i], format='E',unit='yr')
    c4s=fits.Column(name='age75', array=age75[N_data*(i-1):N_data*i], format='E',unit='yr')
    c5s=fits.Column(name='age90', array=age90[N_data*(i-1):N_data*i], format='E',unit='yr')
    
    
    
    hdr1=hdul1[1].header
    hdr1.append(('TCOMM18', 'Lookback time at the epoch of 10% of stars formed'))
    hdr1.append(('TCOMM19', 'Lookback time at the epoch of 25% of stars formed'))
    hdr1.append( ('TCOMM20', 'Lookback time at the epoch of 50% of stars formed'))
    hdr1.append( ('TCOMM21', 'Lookback time at the epoch of 75% of stars formed'))
    hdr1.append( ('TCOMM22', 'Lookback time at the epoch of 90% of stars formed'))
    
    col_tot=hdul1[1].columns + c1s+c2s+c3s+c4s+c5s
    
    new_file=work_dir+suffix_physpar+'_{:03d}_physpar_wagef.fits'
    hs=fits.BinTableHDU.from_columns(col_tot, header=hdul1[1].header)
    hs.writeto(new_file.format(i), overwrite=True)
    
    hdul1.close()


hdul.close()
        