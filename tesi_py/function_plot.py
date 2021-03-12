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

def density_map(x,y,par,statistic,n_size_par,title=[''], x_label='', y_label='',vmin=[], vmax=[], nx=1, ny=1, figsize=(5,5)):
    
    n_stat=np.size(statistic)
    
    
    if (nx==1 & ny==1):
        
        stat_par=stats.binned_statistic_2d(x, y, par[0:n_size_par], bins=50, statistic=statistic[0])
        
        y_g,x_g=np.meshgrid(stat_par.y_edge,stat_par.x_edge)
        fig, axs=plt.subplots(figsize=figsize)
        im=axs.pcolormesh(x_g,y_g,(stat_par.statistic),cmap=cm.rainbow, vmin=vmin[0], vmax=vmax[0])
        fig.colorbar(im, ax=axs)
        axs.set_title(title[0])
        axs.set_xlabel(x_label)
        axs.set_ylabel(y_label)
        
    if (nx*ny != n_stat):
        print('something is wrong! Take a look on the number of the density maps!')
        
    if(nx!=0 & ny!=0 & nx*ny==n_stat):
        fig, axs=plt.subplots(nx, ny, figsize=figsize)
        
        for i_x in range(0, nx):
            for i_y in range(0, ny):    
                for i_stat in range(0,n_stat):
                    _stat_par=stats.binned_statistic_2d(x,y,par[(i_stat*n_size_par):(n_size_par+i_stat*n_size_par)])
                    _y_g, _x_g=np.meshgrid(_stat_par.y_edge, _stat_par.x_edge)
                    _im= axs[i_x,i_y].pcolormesh(_x_g,_y_g,(_stat_par.statistic),cmap=cm.rainbow, vmin=vmin[i_stat], vmax=vmax[i_stat])
                    fig.colorbar(_im, ax=axs[i_x,i_y])
                    axs[i_x, i_y].set_title(title[i_stat])
                axs[nx-1, i_y].set_ylabel(y_label)
            axs[i_x, ny-1].set_xlabel(x_label)

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