#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 25 11:40:26 2020

@author: edoardo
"""

import astropy.io.fits as fits


hdul=fits.open('t_f.fits')
hdul1=fits.open('age_50_taun.fits')


N_models=500000

AGE=[0.0]*N_models
age50=[0.0]*N_models
tau_n=[0.0]*N_models
t25b=[0.0]*N_models
t50b=[0.0]*N_models
t75b=[0.0]*N_models

for k in range(0,N_models):
    AGE[k]=hdul1[1].data[k][0]
    age50[k]=hdul1[1].data[k][1]
    tau_n[k]=hdul1[1].data[k][2]
    t25b[k]=hdul[1].data[k][0]
    t50b[k]=hdul[1].data[k][1]
    t75b[k]=hdul[1].data[k][2]

c1=fits.Column(name='age', array=AGE, format='E')
c2=fits.Column(name='age50', array=age50, format='E')
c3=fits.Column(name='tau/t_obs', array=tau_n, format='E')
#f=fits.BinTableHDU.from_columns([c1,c2,c3])
#f.writeto('age_50_taun.fits')

c4=fits.Column(name='t25', array=t25b, format='E')
c5=fits.Column(name='t50', array=t50b, format='E')
c6=fits.Column(name='t75', array=t75b, format='E')
h=fits.BinTableHDU.from_columns([c1,c2,c3,c4,c5,c6])
h.writeto('tesi.fits')
