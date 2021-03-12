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

#density_map_5p is a thesis format...

def density_map_5p(x,y,par,mock_par,mock_err, par_name='', x_label='', y_label='', vmin=[], vmax=[], nx=3, ny=2, figsize=(5,5)):
    import function_plot    
    
    median_par=stats.binned_statistic_2d(x, y, par, bins=50, statistic='median')
    y_g, x_g=np.meshgrid(median_par.y_edge, median_par.x_edge)
    bias_par=stats.binned_statistic_2d(x, y, mock_par-par, bins=50, statistic='median')
    rms_1684_bias=stats.binned_statistic_2d(x, y, mock_par-par, bins=50, statistic=function_plot.rms_1684)
    median_bayes_err=stats.binned_statistic_2d(x,y,mock_err, bins=50, statistic='median')
    err_norm=stats.binned_statistic_2d(x,y,(mock_par-par)/(mock_err), bins=50, statistic=function_plot.rms_1684) 
    
    
    fig1,axs1=plt.subplots(ny,nx,figsize=figsize)
    
    _im=axs1[0,0].pcolormesh(x_g, y_g, median_par.statistic, cmap=cm.gist_rainbow, vmin=vmin[0], vmax=vmax[0])
    fig1.colorbar(_im, ax=axs1[0,0])
    axs1[0,0].set_title('log'+par_name)
    
    _im=axs1[0,2].pcolormesh(x_g,y_g, err_norm.statistic, cmap=cm.hot, vmin=vmin[2], vmax=vmax[2])
    axs1[0,2].set_facecolor('grey')
    #_im=axs1[0,2].pcolormesh(x_g,y_g,rms_1684_bias.statistic/median_bayes_err.statistic , cmap=cm.rainbow, vmin=vmin[2], vmax=vmax[2], facecolor='grey')
    fig1.colorbar(_im, ax=axs1[0,2])
    axs1[0,2].set_title(par_name+'rms1684_(out-in)/errbayes')

    _im=axs1[0,1].pcolormesh(x_g, y_g, bias_par.statistic, cmap=cm.Spectral, vmin=vmin[1], vmax=vmax[1])
    fig1.colorbar(_im, ax=axs1[0,1])
    axs1[0,1].set_title(par_name+'_bias_out-in')

    axs1[1,0].scatter(par, mock_par, s=0.01)
    axs1[1,0].plot(par, par, color='red')

    _im=axs1[1,1].pcolormesh(x_g, y_g, rms_1684_bias.statistic, cmap=cm.hot, vmin=vmin[3], vmax=vmax[3])
    fig1.colorbar(_im, ax=axs1[1,1])
    axs1[1,1].set_facecolor('grey')
    axs1[1,1].set_title(par_name+'_rms1684_out-in')
    #print('[1,1]:', np.nanmax(rms_1684_bias.statistic), np.nanmin(rms_1684_bias.statistic))
    #c=np.nanargmax(rms_1684_bias.statistic, axis=0)
    
    _im=axs1[1,2].pcolormesh(x_g, y_g, median_bayes_err.statistic, cmap=cm.hot, vmin=vmin[4], vmax=vmax[4])
    fig1.colorbar(_im, ax=axs1[1,2])
    axs1[1,2].set_facecolor('grey')
    axs1[1,2].set_title(par_name+'_err_bayes')
    #print('[1,2]:', np.nanmax(median_bayes_err.statistic), np.nanmin(median_bayes_err.statistic))
    #a=np.nanargmin(median_bayes_err.statistic)
    #idx_x=a%np.size(x_g)
    #idx_y=int(a/np.size(x_g))
    #print('hdhg_min:', np.reshape(y_g, -1)[a])
    #print('d4000n_min:', np.reshape(x_g, -1)[a])
    #print('total 0:', np.nansum(median_bayes_err.statistic==0.0))
    axs1[1,2].set_xlabel(x_label)
    axs1[1,0].set_xlabel(par_name)
    axs1[0,0].set_ylabel(y_label)
    axs1[1,1].set_xlabel(x_label)
    axs1[1,0].set_ylabel(par_name+'_mock')
    
    return fig1


def diff_density_map(x,y,par1,par2,statistic,name1='',name2='',xlabel='',ylabel='',figsize=(10,25),vmin=None,vmax=None):
    stat_diffpar=stats.binned_statistic_2d(x,y,par1-par2,bins=50,statistic=statistic)
    
    y_g,x_g=np.meshgrid(stat_diffpar.y_edge, stat_diffpar.x_edge)
    
    fig, ax=plt.subplots(figsize=figsize)
    im=ax.pcolormesh(x_g, y_g, stat_diffpar.statistic,cmap=cm.Spectral,vmin=vmin, vmax=vmax)
    fig.colorbar(im, ax=ax)
    ax.set_facecolor('grey')
    ax.set_title(name1+'-'+name2)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    
    return fig









