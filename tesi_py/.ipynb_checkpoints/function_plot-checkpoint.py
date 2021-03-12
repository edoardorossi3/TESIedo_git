#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  9 14:48:49 2021

@author: edoardo
"""
import scipy.stats as stats
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm

def rms_1684(z):
    perc_16=np.percentile(z,16)
    perc_84=np.percentile(z,84)
    r_1684=(perc_84-perc_16)/2.0
    return r_1684

def density_map(x,y,par,statistic,title="title"):
    
    median_par=stats.binned_statistic_2d(x, y, par, bins=50, statistic=statistic)
    
    y_g,x_g=np.meshgrid(median_par.y_edge,median_par.x_edge)
    fig, axs=plt.subplots(figsize=(5,5))
    im=axs.pcolormesh(x_g,y_g,(median_par.statistic),cmap=cm.rainbow)
    fig.colorbar(im, ax=axs)
    
