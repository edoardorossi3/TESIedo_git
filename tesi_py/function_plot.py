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

def perc_84(z):
    perc84=np.percentile(z,84)
    return perc84

def perc_16(z):
    perc16=np.percentile(z,16)
    return perc16

#density_map_5p is a thesis format...

def density_map_5p(x,y,par,mock_par,mock_err, par_name='', x_label='', y_label='', vmin=[], vmax=[], nx=3, ny=3, figsize=(5,10)):
    import function_plot 
    
    #nan_par=np.isnan(par)
    #nonan_par= ~nan_par
    #par=par[nonan_par]
    #nan_mock_par=np.isnan(mock_par)
    #nonan_mock_par= ~nan_mock_par
    #mock_par=mock_par[nonan_mock_par]
    
    idx_fin=np.isfinite(par*mock_par)
    par=par[idx_fin]
    mock_par=mock_par[idx_fin]
    mock_err=mock_err[idx_fin]
    x_fin=x[idx_fin]
    y_fin=y[idx_fin]
    idx_nofin= ~idx_fin
    print('total deleted (no finite values):', np.sum(idx_nofin))
    
    median_par=stats.binned_statistic_2d(x_fin, y_fin, par, bins=50, statistic='median')
    median_mock=stats.binned_statistic_2d(x_fin, y_fin, mock_par, bins=50, statistic='median')
    y_g, x_g=np.meshgrid(median_par.y_edge, median_par.x_edge)
    bias_par=stats.binned_statistic_2d(x_fin, y_fin, (mock_par-par), bins=50, statistic='median')
    rms_1684_bias=stats.binned_statistic_2d(x_fin, y_fin, (mock_par-par), bins=50, statistic=function_plot.rms_1684)
    median_bayes_err=stats.binned_statistic_2d(x_fin,y_fin,mock_err, bins=50, statistic='median')
    err_norm=stats.binned_statistic_2d(x_fin,y_fin,(mock_par-par)/(mock_err), bins=50, statistic=function_plot.rms_1684)
    
    median_in=stats.binned_statistic(par, par, statistic='median', bins=50)
    median_out=stats.binned_statistic(par, mock_par, statistic='median', bins=50)
    
    perc_84_in=stats.binned_statistic((par),(par), statistic=function_plot.perc_84, bins=50)
    perc_84_out=stats.binned_statistic((par),(mock_par), statistic=function_plot.perc_84, bins=50)
    perc_16_in=stats.binned_statistic((par),(par), statistic=function_plot.perc_16, bins=50)
    perc_16_out=stats.binned_statistic((par),(mock_par), statistic=function_plot.perc_16, bins=50)

    
    
    fig1,axs1=plt.subplots(ny,nx,figsize=figsize)
    
    _im=axs1[0,0].pcolormesh(x_g, y_g, median_par.statistic, cmap=cm.gist_rainbow, vmin=vmin[0], vmax=vmax[0])
    fig1.colorbar(_im, ax=axs1[0,0])
    axs1[0,0].set_title(par_name)
    
    _im=axs1[0,1].pcolormesh(x_g, y_g, median_mock.statistic, cmap=cm.gist_rainbow, vmin=vmin[0], vmax=vmax[0])
    fig1.colorbar(_im, ax=axs1[0,1])
    axs1[0,1].set_title(par_name+'_mock')
    
    #axs1[0,1].axis('off')
    axs1[0,2].axis('off')
    
    istpar=np.histogram((par), bins=50)
    istmock=np.histogram((mock_par), bins=50)
    frac_par=istpar[0]/np.size(par)
    frac_mock=istmock[0]/np.size(mock_par)
    frac_par=np.append(frac_par[0],frac_par)
    frac_mock=np.append(frac_mock[0], frac_mock)
    axs1[1,0].step((istpar[1]),frac_par,label='par')
    axs1[1,0].step((istmock[1]),frac_mock,label='mock')
    axs1[1,0].legend(loc='upper right')
    med_par=np.median(par)
    med_mock=np.median(mock_par)
    axs1[1,0].plot([med_par,med_par], axs1[1,0].get_ylim())
    axs1[1,0].plot([med_mock,med_mock], axs1[1,0].get_ylim())
  
    
    _im=axs1[1,2].pcolormesh(x_g,y_g, err_norm.statistic, cmap=cm.hot, vmin=vmin[2], vmax=vmax[2])
    axs1[1,2].set_facecolor('#d8dcd6')
    #_im=axs1[0,2].pcolormesh(x_g,y_g,rms_1684_bias.statistic/median_bayes_err.statistic , cmap=cm.rainbow, vmin=vmin[2], vmax=vmax[2], facecolor='grey')
    fig1.colorbar(_im, ax=axs1[1,2])
    axs1[1,2].set_title(par_name+'rms1684_(out-in)/errbayes')

    _im=axs1[1,1].pcolormesh(x_g, y_g, bias_par.statistic, cmap=cm.Spectral, vmin=vmin[1], vmax=vmax[1])
    fig1.colorbar(_im, ax=axs1[1,1])
    axs1[1,1].set_title(par_name+'_bias_out-in')
    axs1[1,1].set_facecolor('#d8dcd6')


    axs1[2,0].scatter((par), (mock_par), s=0.1)
    axs1[2,0].plot((par), (par), color='red')
    axs1[2,0].plot(perc_84_in.statistic, perc_84_out.statistic, color='orange')
    axs1[2,0].plot(perc_16_in.statistic, perc_16_out.statistic, color='orange')
    axs1[2,0].plot(median_in.statistic, median_out.statistic, color='green')


    _im=axs1[2,1].pcolormesh(x_g, y_g, rms_1684_bias.statistic, cmap=cm.hot, vmin=vmin[3], vmax=vmax[3])
    fig1.colorbar(_im, ax=axs1[2,1])
    axs1[2,1].set_facecolor('#d8dcd6')
    axs1[2,1].set_title(par_name+'_rms1684_out-in')
    #print('[1,1]:', np.nanmax(rms_1684_bias.statistic), np.nanmin(rms_1684_bias.statistic))
    #c=np.nanargmax(rms_1684_bias.statistic, axis=0)
    
    _im=axs1[2,2].pcolormesh(x_g, y_g, median_bayes_err.statistic, cmap=cm.hot, vmin=vmin[4], vmax=vmax[4])
    fig1.colorbar(_im, ax=axs1[2,2])
    axs1[2,2].set_facecolor('#d8dcd6')
    axs1[2,2].set_title(par_name+'_err_bayes')
    #print('[1,2]:', np.nanmax(median_bayes_err.statistic), np.nanmin(median_bayes_err.statistic))
    #a=np.nanargmin(median_bayes_err.statistic)
    #idx_x=a%np.size(x_g)
    #idx_y=int(a/np.size(x_g))
    #print('hdhg_min:', np.reshape(y_g, -1)[a])
    #print('d4000n_min:', np.reshape(x_g, -1)[a])
    #print('total 0:', np.nansum(median_bayes_err.statistic==0.0))
    axs1[2,2].set_xlabel(x_label)
    axs1[2,0].set_xlabel(par_name)
    axs1[0,0].set_ylabel(y_label)
    axs1[1,0].set_ylabel('fraction of model')
    axs1[2,1].set_xlabel(x_label)
    axs1[2,0].set_ylabel(par_name+'_mock')
    
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

def prior_comp(par_binned, par, mock_par, n_lim=3, figsize=(25,15), limits=[None, None, None]):
    fig, axs=plt.subplots(2,2,figsize=figsize)
    
    median=np.median(par)
    istpar_tot=np.histogram(par,bins=50)
    frac_partot=istpar_tot[0]/np.size(par)
    frac_partot=np.append(frac_partot[0],frac_partot)

    
    _i=np.argwhere(par_binned<limits[0])
    _idx=_i.reshape(np.shape(_i)[0])
    istpar=np.histogram(par[_idx], bins=50)
    istmock=np.histogram(mock_par[_idx], bins=50)
    frac_par=istpar[0]/np.size(par[_idx])
    frac_mock=istmock[0]/np.size(mock_par[_idx])
    frac_par=np.append(frac_par[0],frac_par)
    frac_mock=np.append(frac_mock[0], frac_mock)
    axs[0,0].step(istpar_tot[1],frac_partot, label='total')
    axs[0,0].step(istpar[1], frac_par, label='par')
    axs[0,0].step(istmock[1], frac_mock, label='mock')
    axs[0,0].legend(loc='upper right')
    axs[0,0].plot([median, median], axs[0,0].get_ylim())
    
    _i=np.argwhere(np.logical_and((par_binned<limits[1]),(par_binned>limits[0])))
    _idx=_i.reshape(np.shape(_i)[0])
    istpar=np.histogram(par[_idx], bins=50)
    istmock=np.histogram(mock_par[_idx], bins=50)
    frac_par=istpar[0]/np.size(par[_idx])
    frac_mock=istmock[0]/np.size(mock_par[_idx])
    frac_par=np.append(frac_par[0],frac_par)
    frac_mock=np.append(frac_mock[0], frac_mock)
    axs[0,1].step(istpar_tot[1],frac_partot,label='total')
    axs[0,1].step(istpar[1], frac_par,label='par')
    axs[0,1].step(istmock[1], frac_mock,label='mock')
    axs[0,1].legend(loc='upper right')
    axs[0,1].plot([median, median], axs[0,1].get_ylim())

    
    _i=np.argwhere( np.logical_and((par_binned<limits[2]),(par_binned>limits[1])) )
    _idx=_i.reshape(np.shape(_i)[0])
    istpar=np.histogram(par[_idx], bins=50)
    istmock=np.histogram(mock_par[_idx], bins=50)
    frac_par=istpar[0]/np.size(par[_idx])
    frac_mock=istmock[0]/np.size(mock_par[_idx])
    frac_par=np.append(frac_par[0],frac_par)
    frac_mock=np.append(frac_mock[0], frac_mock)
    axs[1,0].step(istpar_tot[1],frac_partot,label='total')
    axs[1,0].step(istpar[1], frac_par,label='par')
    axs[1,0].step(istmock[1], frac_mock,label='mock')
    axs[1,0].legend(loc='upper right')
    axs[1,0].plot([median, median], axs[1,0].get_ylim())

    
    _i=np.argwhere(par_binned>limits[2])
    _idx=_i.reshape(np.shape(_i)[0])
    istpar=np.histogram(par[_idx], bins=50)
    istmock=np.histogram(mock_par[_idx], bins=50)
    frac_par=istpar[0]/np.size(par[_idx])
    frac_mock=istmock[0]/np.size(mock_par[_idx])
    frac_par=np.append(frac_par[0],frac_par)
    frac_mock=np.append(frac_mock[0], frac_mock)
    axs[1,1].step(istpar_tot[1],frac_partot,label='total')
    axs[1,1].step(istpar[1], frac_par,label='par')
    axs[1,1].step(istmock[1], frac_mock,label='mock')
    axs[1,1].legend(loc='upper right')
    axs[1,1].plot([median, median], axs[1,1].get_ylim())


    return fig

def scatter_comp(par_binned, par, mock_par, limits=[None, None, None], figsize=(25,15)):
    fig, axs=plt.subplots(2,2,figsize=figsize)
    
    _i=np.argwhere(par_binned<limits[0])
    _idx=_i.reshape(np.shape(_i)[0])
    axs[0,0].scatter(par[_idx], mock_par[_idx], s=0.1)
    axs[0,0].plot([5.8,10.5], [5.8,10.5])


    
    _i=np.argwhere(np.logical_and((par_binned<limits[1]),(par_binned>limits[0])))
    _idx=_i.reshape(np.shape(_i)[0])
    axs[0,1].scatter(par[_idx], mock_par[_idx], s=0.1)
    axs[0,1].plot([5.8,10.5], [5.8,10.5])
   
    _i=np.argwhere( np.logical_and((par_binned<limits[2]),(par_binned>limits[1])) )
    _idx=_i.reshape(np.shape(_i)[0])
    axs[1,0].scatter(par[_idx], mock_par[_idx], s=0.1)
    axs[1,0].plot([5.8,10.5], [5.8,10.5])

    _i=np.argwhere(par_binned>limits[2])
    _idx=_i.reshape(np.shape(_i)[0])
    axs[1,1].scatter(par[_idx], mock_par[_idx], s=0.8)
    axs[1,1].plot([5.8,10.5], [5.8,10.5])
    
    axs[0,0].set_ylabel('out')
    axs[1,0].set_ylabel('out')
    axs[1,0].set_xlabel('in')
    axs[1,1].set_xlabel('in')
    return fig
