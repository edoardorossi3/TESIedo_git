#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec  1 19:59:10 2020

@author: edoardo
"""

import numpy as np
import astropy.io.fits as fits

def bisection(f, t_min, t_max, toll):
    if f(t_min)*f(t_max)>=0:
        return None
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
    #m_sndg=m_sndg/(1-np.exp(-0.5*(t_obs/tau)**2)+np.sum(A_burst))
    return m_sndg

work_dir='/home/edoardo/Desktop/TESI/models/Sandage_varZ_v4.2eq_CB16MILES_ChFall/'

f_name=work_dir+'sandage_varZ_v4.2eq.fits'  #name of the file from that we get tform, tau, nbursts, ecc...
hdul=fits.open(f_name)
data=hdul[1].data
N_models=np.size(data)
#AGE=[0.0]*N_models
#age50=[0.0]*N_models
#tau_n=[0.0]*N_models
t10b=[0.0]*N_models
t25b=[0.0]*N_models
t50b=[0.0]*N_models
t75b=[0.0]*N_models
t90b=[0.0]*N_models

for k in range(0, N_models):
    t_obs=10.0**(data[k][0]) #t_obs in yr
    tau=10.0**(data[k][1]) #tau in yr
    Nburst=data[k][5]
    A_burst=data[k][7]
    t_burst=[0.0]*6
    age_burst=[0.0]*6
    age_burst=np.where(data[k][6]!=0, 10**(data[k][6]), age_burst) #age_burst in yr
    t_burst=np.where(age_burst==0, t_burst, t_obs-age_burst)
     
                 
    #Mass_tot=Massb(t_obs, tau, t_obs, A_burst, t_burst, Nburst)                
    Mass_tot=1.0+np.sum(A_burst)    
                 
    frac50=lambda x: Massb(x, tau, t_obs, A_burst, t_burst, Nburst)/Mass_tot-0.5
    frac25=lambda x: Massb(x, tau, t_obs, A_burst, t_burst, Nburst)/Mass_tot-0.25
    frac75=lambda x: Massb(x, tau, t_obs, A_burst, t_burst, Nburst)/Mass_tot-0.75
    frac90=lambda x: Massb(x, tau, t_obs, A_burst, t_burst, Nburst)/Mass_tot-0.9
    frac10=lambda x: Massb(x, tau, t_obs, A_burst, t_burst, Nburst)/Mass_tot-0.1
    
    t50b[k]=bisection(frac50, 0.0, t_obs, 0.01)
    t25b[k]=bisection(frac25, 0.0, t_obs, 0.01)
    t75b[k]=bisection(frac75, 0.0, t_obs, 0.01)
    t90b[k]=bisection(frac90, 0.0, t_obs, 0.01)
    t10b[k]=bisection(frac10, 0.0, t_obs, 0.01)
    
    #t_sndg=(-t_obs*np.exp(-0.5*(t_obs/tau)**2)+np.sqrt(np.pi/2.0)*tau*math.erf(t_obs/(np.sqrt(2.0)*tau)))/(1.0-np.exp(-0.5*(t_obs/tau)**2))
    #t_mw=(t_sndg+np.sum(A_burst*t_burst))/Mass_tot
    
    
    #AGE[k]=t_obs-t_mw
    #age50[k]=t_obs-t50b
    
#c1=fits.Column(name='age', array=AGE, format='E')
#c2=fits.Column(name='age50', array=age50, format='E')
#f=fits.BinTableHDU.from_columns([c1,c2])
#f.writeto('age_50_analitic.fits')

c3=fits.Column(name='t10', array=t10b, format='E')    
c4=fits.Column(name='t25', array=t25b, format='E')
c5=fits.Column(name='t50', array=t50b, format='E')
c6=fits.Column(name='t75', array=t75b, format='E')
c7=fits.Column(name='t90', array=t90b, format='E')
h=fits.BinTableHDU.from_columns([c3,c4,c5,c6,c7])
h.writeto(work_dir+'t_f.fits', overwrite=True)
    
                    
                     
                     
                     
                     
                     