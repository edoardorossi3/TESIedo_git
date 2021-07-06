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
import os
import astropy.io.fits as fits

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

def density_map_5p(x,y,par,mock_par,mock_err,x_err, y_err, par_name='', x_label='', y_label='', vmin=[], vmax=[], nx=3, ny=3, figsize=(5,10), s=0.01, N_bins=50, statistic=False):
    import function_plot 
    from matplotlib.transforms import Bbox
    from astropy.visualization import simple_norm
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
    
    err_x, err_y=function_plot.ellipse(np.mean(x_err), np.mean(y_err), x_c=np.percentile(x,84), y_c=np.percentile(y, 84))
    print('total deleted (no finite values):', np.sum(idx_nofin))
    
    density_in=stats.binned_statistic_2d(x_fin, y_fin, par,bins=N_bins, statistic='count')
    
    mean_par=stats.binned_statistic_2d(x_fin, y_fin, par, bins=N_bins, statistic='mean', expand_binnumbers=True)
    std_par=stats.binned_statistic_2d(x_fin, y_fin, par, bins=N_bins, statistic='std')
    std_par=np.where(std_par.statistic==0, np.nan, std_par.statistic)
    
    p16_par=stats.binned_statistic_2d(x_fin, y_fin, par, bins=N_bins, statistic=function_plot.perc_16)
    p84_par=stats.binned_statistic_2d(x_fin, y_fin, par, bins=N_bins, statistic=function_plot.perc_84)

    mean_mock=stats.binned_statistic_2d(x_fin, y_fin, mock_par, bins=N_bins, statistic='mean')
    std_mock=stats.binned_statistic_2d(x_fin, y_fin, mock_par, bins=N_bins, statistic='std')
    std_mock=np.where(std_mock.statistic==0, np.nan, std_mock.statistic)

    p16_mock=stats.binned_statistic_2d(x_fin, y_fin, mock_par, bins=N_bins, statistic=function_plot.perc_16)
    p84_mock=stats.binned_statistic_2d(x_fin, y_fin, mock_par, bins=N_bins, statistic=function_plot.perc_84)

    y_g, x_g=np.meshgrid(mean_par.y_edge, mean_par.x_edge)
    bias_par=stats.binned_statistic_2d(x_fin, y_fin, (mock_par-par), bins=N_bins, statistic='mean')
    mean2diff_bias=(stats.binned_statistic_2d(x_fin, y_fin, (mock_par-par)**2, bins=N_bins, statistic='mean'))
    
    
    median_bayes_err=stats.binned_statistic_2d(x_fin,y_fin,mock_err, bins=N_bins, statistic='mean')
    err_norm=stats.binned_statistic_2d(x_fin,y_fin,(mock_par-par)/(mock_err), bins=N_bins, statistic=function_plot.rms_1684)
    
    running_med_in=stats.binned_statistic(par, par, bins=N_bins, statistic='mean') 
    #running_p16_in=stats.binned_statistic(par, par, bins=20, statistic=function_plot.perc_16) 
    #running_p84_in=stats.binned_statistic(par, par, bins=20, statistic=function_plot.perc_84) 

    running_med_out=stats.binned_statistic(par, mock_par, bins=N_bins, statistic='mean')
    #running_p16_out=stats.binned_statistic(mock_par, mock_par, bins=20, statistic=function_plot.perc_16)
    #running_p84_out=stats.binned_statistic(mock_par, mock_par, bins=20, statistic=function_plot.perc_84)
    
    #perc_84_in=stats.binned_statistic((par),(par), statistic=function_plot.perc_84, bins=50)
    #perc_84_out=stats.binned_statistic((par),(mock_par), statistic=function_plot.perc_84, bins=50)
    #perc_16_in=stats.binned_statistic((par),(par), statistic=function_plot.perc_16, bins=50)
    #perc_16_out=stats.binned_statistic((par),(mock_par), statistic=function_plot.perc_16, bins=50)

    mean_in=np.reshape(mean_par.statistic, -1)
    mean_out=np.reshape(mean_mock.statistic, -1)
    
    perc_84_in=np.reshape(p84_par.statistic, -1)
    perc_84_out=np.reshape(p84_mock.statistic, -1)
    perc_16_in=np.reshape(p16_par.statistic, -1)
    perc_16_out=np.reshape(p16_mock.statistic, -1)
    
    
    fig1,axs1=plt.subplots(ny,nx,figsize=figsize)
    
    for ix in range(0,3):
        for iy in range(0,3):
            axs1[ix, iy].tick_params(labelsize=14)

    #plt.tight_layout(pad=2, h_pad=2, w_pad=2)
    
    _im=axs1[0,0].pcolormesh(x_g, y_g, mean_par.statistic, cmap=cm.gnuplot2, vmin=vmin[0], vmax=vmax[0])
    cbar=fig1.colorbar(_im, ax=axs1[0,0])
    cbar.set_label(r'$<$'+par_name+r'$_{in}>$', fontsize=20)
    cbar.ax.tick_params(labelsize=18)
    #axs1[0,0].set_title('mean in')
    axs1[0,0].set_facecolor('#d8dcd6')
    axs1[0,0].plot(err_x, err_y, color='red')
    
    norm=simple_norm(std_par, 'asinh',min_cut=vmin[1], max_cut=vmax[1])
    _im=axs1[0,1].pcolormesh(x_g, y_g, std_par, cmap=cm.nipy_spectral, norm=norm)
    cbar=fig1.colorbar(_im, ax=axs1[0,1])
    cbar.set_label(r'$r.m.s._{in}$', fontsize=20)
    cbar.ax.tick_params(labelsize=18)
    #axs1[0,1].set_title('std in')
    axs1[0,1].set_facecolor('#d8dcd6')
    axs1[0,1].plot(err_x, err_y, color='red')

    tmp=density_in.statistic
    tmp=np.where(tmp==0, np.nan, tmp)
    _im=axs1[0,2].pcolormesh(x_g, y_g, tmp, cmap=cm.copper_r)
    cbar=fig1.colorbar(_im, ax=axs1[0,2])
    cbar.set_label(r'$models_{in}$', fontsize=20)
    cbar.ax.tick_params(labelsize=18)
    axs1[0,2].set_facecolor('#d8dcd6')
    axs1[0,2].plot(err_x, err_y, color='red')
    
    _im=axs1[1,0].pcolormesh(x_g, y_g, mean_mock.statistic, cmap=cm.gnuplot2, vmin=vmin[0], vmax=vmax[0])
    cbar=fig1.colorbar(_im, ax=axs1[1,0])
    cbar.set_label(r'$<$'+par_name+r'$_{out}>$', fontsize=20)
    cbar.ax.tick_params(labelsize=18)
    axs1[1,0].set_facecolor('#d8dcd6')
    axs1[1,0].plot(err_x, err_y, color='red')
        
    norm=simple_norm(std_mock, 'asinh', min_cut=vmin[1], max_cut=vmax[1])
    _im=axs1[1,1].pcolormesh(x_g, y_g, std_mock, cmap=cm.nipy_spectral, norm=norm)
    cbar=fig1.colorbar(_im, ax=axs1[1,1])
    cbar.set_label(r'$r.m.s._{out}$', fontsize=20)
    cbar.ax.tick_params(labelsize=18)
    axs1[1,1].set_facecolor('#d8dcd6')
    axs1[1,1].plot(err_x, err_y, color='red')
    
    
    #axs1[1,2].errorbar(mean_in, mean_out,marker='o', xerr=(perc_84_in-perc_16_in)/2.0, yerr=(perc_84_out-perc_16_out)/2.0, linestyle='None', elinewidth=0.5, ms=s)
    

   
    for ix in range(0, N_bins):
        for iy in range(0, N_bins):
            _idx_sel=np.where((x>mean_par.x_edge[ix])&(x<mean_par.x_edge[ix+1])&(y>mean_par.y_edge[iy])&(y<mean_par.y_edge[iy+1]))[0]
            _par=par[_idx_sel]
            _mock_par=mock_par[_idx_sel]
            if (np.size(_idx_sel)>5):
                cov=np.cov(_par, _mock_par)
                pearson=stats.pearsonr(_par, _mock_par)[0]
                #print(_idx_sel)
                #print(ix, iy,cov)
                eig_val, eig_vec=np.linalg.eig(cov)
                v0=eig_vec[:,0]
                v1=eig_vec[:,1]
                phi=np.arctan(v0[1]/v0[0])
                x_ell, y_ell=function_plot.ellipse(np.sqrt(eig_val[0]),np.sqrt(eig_val[1]), -phi, x_c=np.mean(_par), y_c=np.mean(_mock_par))
                axs1[1,2].scatter(np.mean(_par), np.mean(_mock_par),marker="X", color='red', s=0.1)
                #axs1[1,2].plot([np.mean(_par), np.mean(_par)+np.sqrt(eig_val[0])*np.cos(phi)], [np.mean(_mock_par), np.mean(_mock_par)+np.sqrt(eig_val[0])*np.sin(phi)],linewidth=0.05, color='red')
                #axs1[1,2].plot([np.mean(_par), np.mean(_par)-np.sqrt(eig_val[1])*np.sin(phi)], [np.mean(_mock_par), np.mean(_mock_par)+np.sqrt(eig_val[1])*np.cos(phi)], linewidth=0.05,color='red')
                color='red'
                widht=0.04
                #if ((np.sqrt(eig_val[1])/np.sqrt(eig_val[0])<0.6)&(phi<np.pi/3)&(phi>np.pi/6)):
                if (pearson>0.7):
                    color='green'
                    widht=0.5
                
                axs1[1,2].plot(x_ell, y_ell, color=color,linewidth=widht)

    
    axs1[1,2].axis('equal')
    axs1[1,2].plot([np.min(par), np.max(par)], [np.min(par), np.max(par)], color='black')
    pos=axs1[1,2].get_position()
    pos_new=Bbox([[pos.xmin, pos.ymin], [pos.xmax-0.045, pos.ymax]])
    axs1[1,2].set_position(pos_new)
    #axs1[1,2].set_facecolor('#d8dcd6')

    _im=axs1[2,0].pcolormesh(x_g, y_g, bias_par.statistic, cmap=cm.RdYlGn, vmin=vmin[2], vmax=vmax[2])
    cbar=fig1.colorbar(_im, ax=axs1[2,0])
    cbar.set_label(r'$<out-in>$', fontsize=20)
    cbar.ax.tick_params(labelsize=18)
    #axs1[2,0].set_title('mean out-in')
    axs1[2,0].set_facecolor('#d8dcd6')
    axs1[2,0].plot(err_x, err_y, color='red')

    norm=simple_norm(np.sqrt(mean2diff_bias.statistic), 'asinh', min_cut=vmin[3], max_cut=vmax[3])
    _im=axs1[2,1].pcolormesh(x_g, y_g, np.sqrt(mean2diff_bias.statistic), cmap=cm.nipy_spectral, norm=norm)
    cbar=fig1.colorbar(_im, ax=axs1[2,1])
    cbar.set_label(r'$\sqrt{<(out-in)^2>}$', fontsize=20)
    cbar.ax.tick_params(labelsize=18)
    axs1[2,1].set_facecolor('#d8dcd6')
    #axs1[2,1].set_title(r'$\sqrt{mean_{out-in}^2+std_{out-in}}$')
    axs1[2,1].plot(err_x, err_y, color='red')
   
    norm=simple_norm(median_bayes_err.statistic, 'asinh', min_cut=vmin[4], max_cut=vmax[4])
    _im=axs1[2,2].pcolormesh(x_g, y_g, median_bayes_err.statistic, cmap=cm.nipy_spectral, norm=norm)
    cbar=fig1.colorbar(_im, ax=axs1[2,2])
    cbar.set_label(r'$err_{bayes}$', fontsize=20)
    cbar.ax.tick_params(labelsize=18)
    axs1[2,2].set_facecolor('#d8dcd6')
    axs1[2,2].plot(err_x, err_y, color='red')
    
    #axs1[0,1].axis('off')
    #axs1[0,2].axis('off')
    
    axs1[2,2].set_xlabel(x_label, size=18)
    #axs1[2,2].set_ylabel(y_label, size=18)

    axs1[2,0].set_xlabel(x_label, size=18)
    axs1[1,2].set_xlabel(par_name+r'$_{in}$', size=18)

    axs1[0,0].set_ylabel(y_label, size=18)
    axs1[0,0].set_xlabel(x_label, size=18)
    
    #axs1[0,1].set_ylabel(y_label, size=18)
    axs1[0,1].set_xlabel(x_label, size=18)
    
    #axs1[1,1].set_ylabel(y_label, size=18)
    axs1[1,1].set_xlabel(x_label, size=18)
    
    #axs1[0,2].set_ylabel(y_label, size=18)
    axs1[0,2].set_xlabel(x_label, size=18)
    
    axs1[1,0].set_ylabel(y_label, size=18)
    axs1[1,0].set_xlabel(x_label, size=18)

    axs1[2,1].set_xlabel(x_label, size=18)
    #axs1[2,1].set_ylabel(y_label, size=18)

    axs1[2,0].set_ylabel(y_label, size=18)
    axs1[1,2].set_ylabel(par_name+r'$_{out}$', size=18)
    
    if (statistic):
        return std_mock, bias_par.statistic, np.sqrt(mean2diff_bias.statistic), median_bayes_err.statistic, x_g, y_g
    else:
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
    
    return fig, ax, stat_diffpar.statistic


def density_map(x,y,par,statistic,name='',xlabel='',ylabel='',figsize=(10,25),vmin=None,vmax=None, bins=50):
    stat_diffpar=stats.binned_statistic_2d(x,y,par,bins=bins,statistic=statistic)
    
    y_g,x_g=np.meshgrid(stat_diffpar.y_edge, stat_diffpar.x_edge)
    
    fig, ax=plt.subplots(figsize=figsize)
    im=ax.pcolormesh(x_g, y_g, stat_diffpar.statistic,cmap=cm.Spectral,vmin=vmin, vmax=vmax)
    fig.colorbar(im, ax=ax)
    ax.set_facecolor('grey')
    ax.set_title(name)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    
    print(np.sum(stat_diffpar.statistic))
    
    return fig




def prior_comp(par_binned, par, mock_par, n_lim=3, figsize=(15,10), limits=[None, None, None],name_binned_par=''):
    
    idx_fin=np.isfinite(par*mock_par)
    #par=par[idx_fin]
    #mock_par=mock_par[idx_fin]
    idx_nofin= ~idx_fin
    par[idx_nofin]=6.0
    mock_par[idx_nofin]=6.0
    print('total deleted (no finite values):', np.sum(idx_nofin))
    
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
    axs[0,0].set_title(name_binned_par+'<'+str(limits[0]))  
    
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
    axs[0,1].set_title(str(limits[0])+'<'+name_binned_par+'<'+str(limits[1]))

    
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
    axs[1,0].set_title(str(limits[1])+'<'+name_binned_par+'<'+str(limits[2]))

    
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
    axs[1,1].set_title(name_binned_par+'>'+str(limits[2]))


    return fig

def scatter_comp(par_binned, par, mock_par, limits=[None, None, None], figsize=(15,10), name_binned_par=''):
    
    idx_fin=np.isfinite(par*mock_par)
    #par=par[idx_fin]
    #mock_par=mock_par[idx_fin]
    idx_nofin= ~idx_fin
    par[idx_nofin]=6.0
    mock_par[idx_nofin]=6.0
    print('total deleted (no finite values):', np.sum(idx_nofin))
    
    fig, axs=plt.subplots(2,2,figsize=figsize)
    
    _i=np.argwhere(par_binned<limits[0])
    _idx=_i.reshape(np.shape(_i)[0])
    axs[0,0].scatter(par[_idx], mock_par[_idx], s=0.1)
    axs[0,0].plot([np.min(par[_idx]),np.max(par[_idx])], [np.min(par[_idx]),np.max(par[_idx])], color='red')
    axs[0,0].set_title(name_binned_par+'<'+str(limits[0]))


    
    _i=np.argwhere(np.logical_and((par_binned<limits[1]),(par_binned>limits[0])))
    _idx=_i.reshape(np.shape(_i)[0])
    axs[0,1].scatter(par[_idx], mock_par[_idx], s=0.1)
    axs[0,1].plot([np.min(par[_idx]),np.max(par[_idx])], [np.min(par[_idx]),np.max(par[_idx])], color='red')
    axs[0,1].set_title(str(limits[0])+'<'+name_binned_par+'<'+str(limits[1]))

   
    _i=np.argwhere( np.logical_and((par_binned<limits[2]),(par_binned>limits[1])) )
    _idx=_i.reshape(np.shape(_i)[0])
    axs[1,0].scatter(par[_idx], mock_par[_idx], s=0.1)
    axs[1,0].plot([np.min(par[_idx]),np.max(par[_idx])], [np.min(par[_idx]),np.max(par[_idx])], color='red')
    axs[1,0].set_title(str(limits[1])+'<'+name_binned_par+'<'+str(limits[2]))


    _i=np.argwhere(par_binned>limits[2])
    _idx=_i.reshape(np.shape(_i)[0])
    axs[1,1].scatter(par[_idx], mock_par[_idx], s=0.8)
    axs[1,1].plot([np.min(par[_idx]),np.max(par[_idx])], [np.min(par[_idx]),np.max(par[_idx])], color='red')
    axs[1,1].set_title(name_binned_par+'>'+str(limits[2]))

    
    axs[0,0].set_ylabel('out')
    axs[1,0].set_ylabel('out')
    axs[1,0].set_xlabel('in')
    axs[1,1].set_xlabel('in')
    return fig

def scatter_norm(par_binned, par, mock_par, limits=[None, None, None], figsize=(15, 10), name_binned_par='', name_par=''):
    
    idx_fin=np.isfinite(par*mock_par)
    #par=par[idx_fin]
    #mock_par=mock_par[idx_fin]
    idx_nofin= ~idx_fin
    par[idx_nofin]=6.0
    mock_par[idx_nofin]=6.0
    print('total deleted (no finite values):', np.sum(idx_nofin))
    
    fig,axs=plt.subplots(1,2,figsize=figsize)
    
    median_tot=np.median(par)
    median_arr=[0.0]*4
    mock_median_arr=[0.0]*4
    
    _i=np.argwhere(par_binned<limits[0])
    _idx=_i.reshape(np.shape(_i)[0])
    median_arr[0]=np.median(par[_idx])
    mock_median_arr[0]=np.median(mock_par[_idx])
    _medpar_arr=[median_arr[0]]*np.size(par[_idx])
    _medmock_arr=[mock_median_arr[0]]*np.size(mock_par[_idx])
    axs[0].scatter(par[_idx], mock_par[_idx], color='red', s=0.6)
    axs[1].scatter(par[_idx]-_medpar_arr,mock_par[_idx]-_medmock_arr, color='blue', s=0.6)

    
    _i=np.argwhere(np.logical_and((par_binned<limits[1]),(par_binned>limits[0])))
    _idx=_i.reshape(np.shape(_i)[0])
    median_arr[1]=np.median(par[_idx])
    mock_median_arr[1]=np.median(mock_par[_idx])
    _medpar_arr=[median_arr[1]]*np.size(par[_idx])
    _medmock_arr=[mock_median_arr[1]]*np.size(mock_par[_idx])
    axs[0].scatter(par[_idx], mock_par[_idx], color='orange',s=0.6)
    axs[1].scatter(par[_idx]-_medpar_arr,mock_par[_idx]-_medmock_arr, color='blue', s=0.6)

    
    _i=np.argwhere( np.logical_and((par_binned<limits[2]),(par_binned>limits[1])) )
    _idx=_i.reshape(np.shape(_i)[0])
    median_arr[2]=np.median(par[_idx])
    mock_median_arr[2]=np.median(mock_par[_idx])
    _medpar_arr=[median_arr[2]]*np.size(par[_idx])
    _medmock_arr=[mock_median_arr[2]]*np.size(mock_par[_idx])
    axs[0].scatter(par[_idx], mock_par[_idx], color='violet',s=0.6)
    axs[1].scatter(par[_idx]-_medpar_arr,mock_par[_idx]-_medmock_arr, color='blue', s=0.6)


    _i=np.argwhere(par_binned>limits[2])
    _idx=_i.reshape(np.shape(_i)[0])
    median_arr[3]=np.median(par[_idx])
    mock_median_arr[3]=np.median(mock_par[_idx]) 
    _medpar_arr=[median_arr[3]]*np.size(par[_idx])
    _medmock_arr=[mock_median_arr[3]]*np.size(mock_par[_idx])
    axs[0].scatter(par[_idx], mock_par[_idx], color='blue',s=0.6)
    axs[1].scatter(par[_idx]-_medpar_arr,mock_par[_idx]-_medmock_arr, color='blue', s=0.6)


    #axs.scatter(par, mock_par, s=0.8)
    axs[0].scatter(median_arr[0], mock_median_arr[0], color='red',s=100, marker="X", label=name_binned_par+'<'+str(limits[0]))
    axs[0].scatter(median_arr[1], mock_median_arr[1], color='orange',s=100, marker="X", label=str(limits[0])+'<'+name_binned_par+'<'+str(limits[1]))    
    axs[0].scatter(median_arr[2], mock_median_arr[2], color='violet',s=100, marker="X", label=str(limits[1])+'<'+name_binned_par+'<'+str(limits[2]))
    axs[0].scatter(median_arr[3], mock_median_arr[3], color='blue',s=100, marker="X", label=name_binned_par+'>'+str(limits[2]))
    axs[0].plot([median_tot, median_tot], axs[0].get_ylim(), color='red')
    axs[1].plot([0.0, 0.0], axs[1].get_ylim(), color='red')
    axs[1].plot(axs[1].get_xlim(), [0.0, 0.0], color='red')
    axs[1].plot(axs[1].get_xlim(), axs[1].get_xlim(), color='red')

    axs[0].set_xlabel(name_par+'_in')
    axs[0].set_ylabel(name_par+'_out')
    axs[0].legend(loc='upper left')
    axs[1].set_xlabel(name_par+'_in-median')
    axs[1].set_ylabel(name_par+'_out-median')
    
def idx_resol(par, idx_1, idx_2, idx_3, idx_4, idx_5,par_name='', idx_name=['','','','',''], figsize=(10,5),s=1):
    
    fig, axs=plt.subplots(2, 3, figsize=figsize)
    
    axs[0,0].scatter(par, idx_1, s=s)
    axs[0,1].scatter(par, idx_2, s=s)
    axs[0,2].scatter(par, idx_3, s=s)
    axs[1,0].scatter(par, idx_4, s=s)
    axs[1,1].scatter(par, idx_5, s=s)
    
    axs[0,0].set_xlabel(par_name)
    axs[0,1].set_xlabel(par_name)
    axs[0,2].set_xlabel(par_name)
    axs[1,0].set_xlabel(par_name)
    axs[1,1].set_xlabel(par_name)
    
    axs[0,0].set_ylabel(idx_name[0])
    axs[0,1].set_ylabel(idx_name[1])
    axs[0,2].set_ylabel(idx_name[2])
    axs[1,0].set_ylabel(idx_name[3])
    axs[1,1].set_ylabel(idx_name[4])
    
    axs[1,2].axis('off')
    return fig
    
def idx_resol_stat(par,idx_1, idx_2, idx_3, idx_4,  idx_5,par_name='', idx_name=['','','','',''],statistic='median',bins=50, figsize=(10,5),s=1):
    
    fig, axs=plt.subplots(2, 3, figsize=figsize)
    axs[1,2].axis('off')
    
    stat_1=stats.binned_statistic(par, idx_1, statistic=statistic, bins=bins)
    stat_2=stats.binned_statistic(par, idx_2, statistic=statistic, bins=bins)
    stat_3=stats.binned_statistic(par, idx_3, statistic=statistic, bins=bins)
    stat_4=stats.binned_statistic(par, idx_4, statistic=statistic, bins=bins)
    stat_5=stats.binned_statistic(par, idx_5, statistic=statistic, bins=bins)
    
    
    axs[0,0].plot(stat_1.bin_edges[:-1],stat_1.statistic)
    axs[0,1].plot(stat_2.bin_edges[:-1],stat_2.statistic)
    axs[0,2].plot(stat_3.bin_edges[:-1],stat_3.statistic)
    axs[1,0].plot(stat_4.bin_edges[:-1],stat_4.statistic)
    axs[1,1].plot(stat_5.bin_edges[:-1],stat_5.statistic)
    
    axs[0,0].set_xlabel(par_name)
    axs[0,1].set_xlabel(par_name)
    axs[0,2].set_xlabel(par_name)
    axs[1,0].set_xlabel(par_name)
    axs[1,1].set_xlabel(par_name)
    
    axs[0,0].set_ylabel(idx_name[0])
    axs[0,1].set_ylabel(idx_name[1])
    axs[0,2].set_ylabel(idx_name[2])
    axs[1,0].set_ylabel(idx_name[3])
    axs[1,1].set_ylabel(idx_name[4])
    return fig
    
    
def idx_resol_stat4(par,i2,i3,i4,idx_1, idx_2, idx_3, idx_4,  idx_5,x_name='',par_name=['','',''], idx_name=['','','','',''],statistic='median',bins=50, figsize=(10,5),s=1):
    
    fig, axs=plt.subplots(2, 3, figsize=figsize)
    axs[1,2].axis('off')
    
    par2=par[i2]
    par3=par[i3]
    par4=par[i4]
    
    
    
    idx2_1=idx_1[i2]
    idx2_2=idx_2[i2]
    idx2_3=idx_3[i2]
    idx2_4=idx_4[i2]
    idx2_5=idx_5[i2]
    
    idx3_1=idx_1[i3]
    idx3_2=idx_2[i3]
    idx3_3=idx_3[i3]
    idx3_4=idx_4[i3]
    idx3_5=idx_5[i3]
    
    idx4_1=idx_1[i4]
    idx4_2=idx_2[i4]
    idx4_3=idx_3[i4]
    idx4_4=idx_4[i4]
    idx4_5=idx_5[i4]
    
    
    
    stat2_1=stats.binned_statistic(par2, idx2_1, statistic=statistic, bins=bins)
    stat2_2=stats.binned_statistic(par2, idx2_2, statistic=statistic, bins=bins)
    stat2_3=stats.binned_statistic(par2, idx2_3, statistic=statistic, bins=bins)
    stat2_4=stats.binned_statistic(par2, idx2_4, statistic=statistic, bins=bins)
    stat2_5=stats.binned_statistic(par2, idx2_5, statistic=statistic, bins=bins)
    
    stat3_1=stats.binned_statistic(par3, idx3_1, statistic=statistic, bins=bins)
    stat3_2=stats.binned_statistic(par3, idx3_2, statistic=statistic, bins=bins)
    stat3_3=stats.binned_statistic(par3, idx3_3, statistic=statistic, bins=bins)
    stat3_4=stats.binned_statistic(par3, idx3_4, statistic=statistic, bins=bins)
    stat3_5=stats.binned_statistic(par3, idx3_5, statistic=statistic, bins=bins)
    
    stat4_1=stats.binned_statistic(par4, idx4_1, statistic=statistic, bins=bins)
    stat4_2=stats.binned_statistic(par4, idx4_2, statistic=statistic, bins=bins)
    stat4_3=stats.binned_statistic(par4, idx4_3, statistic=statistic, bins=bins)
    stat4_4=stats.binned_statistic(par4, idx4_4, statistic=statistic, bins=bins)
    stat4_5=stats.binned_statistic(par4, idx4_5, statistic=statistic, bins=bins)
    
    
    
    axs[0,0].plot(stat2_1.bin_edges[:-1],stat2_1.statistic, label=par_name[0])
    axs[0,1].plot(stat2_2.bin_edges[:-1],stat2_2.statistic, label=par_name[0])
    axs[0,2].plot(stat2_3.bin_edges[:-1],stat2_3.statistic, label=par_name[0])
    axs[1,0].plot(stat2_4.bin_edges[:-1],stat2_4.statistic, label=par_name[0])
    axs[1,1].plot(stat2_5.bin_edges[:-1],stat2_5.statistic, label=par_name[0])
    
    axs[0,0].plot(stat3_1.bin_edges[:-1],stat3_1.statistic, label=par_name[1])
    axs[0,1].plot(stat3_2.bin_edges[:-1],stat3_2.statistic, label=par_name[1])
    axs[0,2].plot(stat3_3.bin_edges[:-1],stat3_3.statistic, label=par_name[1])
    axs[1,0].plot(stat3_4.bin_edges[:-1],stat3_4.statistic, label=par_name[1])
    axs[1,1].plot(stat3_5.bin_edges[:-1],stat3_5.statistic, label=par_name[1])
    
    axs[0,0].plot(stat4_1.bin_edges[:-1],stat4_1.statistic, label=par_name[2])
    axs[0,1].plot(stat4_2.bin_edges[:-1],stat4_2.statistic, label=par_name[2])
    axs[0,2].plot(stat4_3.bin_edges[:-1],stat4_3.statistic, label=par_name[2])
    axs[1,0].plot(stat4_4.bin_edges[:-1],stat4_4.statistic, label=par_name[2])
    axs[1,1].plot(stat4_5.bin_edges[:-1],stat4_5.statistic, label=par_name[2])
    
    
    axs[0,0].set_xlabel(x_name)
    axs[0,1].set_xlabel(x_name)
    axs[0,2].set_xlabel(x_name)
    axs[1,0].set_xlabel(x_name)
    axs[1,1].set_xlabel(x_name)
    
    axs[0,0].set_ylabel(idx_name[0])
    axs[0,1].set_ylabel(idx_name[1])
    axs[0,2].set_ylabel(idx_name[2])
    axs[1,0].set_ylabel(idx_name[3])
    axs[1,1].set_ylabel(idx_name[4])
    
    axs[0,0].legend(loc='upper right')
    axs[0,1].legend(loc='upper left')
    axs[0,2].legend(loc='upper right')
    axs[1,0].legend(loc='upper right')
    axs[1,1].legend(loc='upper left')
    
    return fig


def idx_resol_stat4col(par,i2,i3,i4,idx_1, idx_2, idx_3, idx_4, x_name='d1090n50',par_name=['age50<9.0','age50_9.0-9.5','age50>9.5'], idx_name=['u-r','g-r','r-i','r-z'],statistic='median',bins=50, figsize=(10,5),s=1):
    
    fig, axs=plt.subplots(2, 3, figsize=figsize)
    axs[1,2].axis('off')
    axs[1,1].axis('off')
    
    par2=par[i2]
    par3=par[i3]
    par4=par[i4]
    
    
    
    idx2_1=idx_1[i2]
    idx2_2=idx_2[i2]
    idx2_3=idx_3[i2]
    idx2_4=idx_4[i2]
    
    idx3_1=idx_1[i3]
    idx3_2=idx_2[i3]
    idx3_3=idx_3[i3]
    idx3_4=idx_4[i3]
    
    idx4_1=idx_1[i4]
    idx4_2=idx_2[i4]
    idx4_3=idx_3[i4]
    idx4_4=idx_4[i4]
    
    
    
    stat2_1=stats.binned_statistic(par2, idx2_1, statistic=statistic, bins=bins)
    stat2_2=stats.binned_statistic(par2, idx2_2, statistic=statistic, bins=bins)
    stat2_3=stats.binned_statistic(par2, idx2_3, statistic=statistic, bins=bins)
    stat2_4=stats.binned_statistic(par2, idx2_4, statistic=statistic, bins=bins)
    
    stat3_1=stats.binned_statistic(par3, idx3_1, statistic=statistic, bins=bins)
    stat3_2=stats.binned_statistic(par3, idx3_2, statistic=statistic, bins=bins)
    stat3_3=stats.binned_statistic(par3, idx3_3, statistic=statistic, bins=bins)
    stat3_4=stats.binned_statistic(par3, idx3_4, statistic=statistic, bins=bins)
    
    stat4_1=stats.binned_statistic(par4, idx4_1, statistic=statistic, bins=bins)
    stat4_2=stats.binned_statistic(par4, idx4_2, statistic=statistic, bins=bins)
    stat4_3=stats.binned_statistic(par4, idx4_3, statistic=statistic, bins=bins)
    stat4_4=stats.binned_statistic(par4, idx4_4, statistic=statistic, bins=bins)
    
    
    
    axs[0,0].plot(stat2_1.bin_edges[:-1],stat2_1.statistic, label=par_name[0])
    axs[0,1].plot(stat2_2.bin_edges[:-1],stat2_2.statistic, label=par_name[0])
    axs[0,2].plot(stat2_3.bin_edges[:-1],stat2_3.statistic, label=par_name[0])
    axs[1,0].plot(stat2_4.bin_edges[:-1],stat2_4.statistic, label=par_name[0])
    
    axs[0,0].plot(stat3_1.bin_edges[:-1],stat3_1.statistic, label=par_name[1])
    axs[0,1].plot(stat3_2.bin_edges[:-1],stat3_2.statistic, label=par_name[1])
    axs[0,2].plot(stat3_3.bin_edges[:-1],stat3_3.statistic, label=par_name[1])
    axs[1,0].plot(stat3_4.bin_edges[:-1],stat3_4.statistic, label=par_name[1])
    
    axs[0,0].plot(stat4_1.bin_edges[:-1],stat4_1.statistic, label=par_name[2])
    axs[0,1].plot(stat4_2.bin_edges[:-1],stat4_2.statistic, label=par_name[2])
    axs[0,2].plot(stat4_3.bin_edges[:-1],stat4_3.statistic, label=par_name[2])
    axs[1,0].plot(stat4_4.bin_edges[:-1],stat4_4.statistic, label=par_name[2])
    
    
    axs[0,0].set_xlabel(x_name)
    axs[0,1].set_xlabel(x_name)
    axs[0,2].set_xlabel(x_name)
    axs[1,0].set_xlabel(x_name)
    
    axs[0,0].set_ylabel(idx_name[0])
    axs[0,1].set_ylabel(idx_name[1])
    axs[0,2].set_ylabel(idx_name[2])
    axs[1,0].set_ylabel(idx_name[3])
    
    axs[0,0].legend(loc='upper right')
    axs[0,1].legend(loc='upper left')
    axs[0,2].legend(loc='upper right')
    axs[1,0].legend(loc='upper right')
    
    return fig
    
    
def chi_q(par,idx1, idx2, idx3, idx4, idx5, mag1, mag2, mag3, mag4, isel, iref,sigma1=None,sigma2=None,sigma3=None,sigma4=None,sigma5=None,sigma6=None,sigma7=None,sigma8=None,sigma9=None, val_ref=-1.00,figsize=(15,10), title='', ylim=[None,None], xlim=[None, None], bins=50, xmin=-1.3, xmax=1.0, toll=0.1, mkplot=True, sigma_obs=False):
    import function_plot as f_plt
    
    
    par_sel=par[isel]
    
    idx1_sel=idx1[isel]
    idx2_sel=idx2[isel]
    idx3_sel=idx3[isel]
    idx4_sel=idx4[isel]
    idx5_sel=idx5[isel]
    
    idx1_ref=np.mean(idx1[iref])
    idx2_ref=np.mean(idx2[iref])
    idx3_ref=np.mean(idx3[iref])
    idx4_ref=np.mean(idx4[iref])
    idx5_ref=np.mean(idx5[iref])
    
    mag1_sel=mag1[isel]
    mag2_sel=mag2[isel]
    mag3_sel=mag3[isel]
    mag4_sel=mag4[isel]
    
    
    mag1_ref=np.mean(mag1[iref])
    mag2_ref=np.mean(mag2[iref])
    mag3_ref=np.mean(mag3[iref])
    mag4_ref=np.mean(mag4[iref])
    
    if sigma_obs:
        sigma_idx1=np.mean(sigma1[iref])
        sigma_idx2=np.mean(sigma2[iref])
        sigma_idx3=np.mean(sigma3[iref])
        sigma_idx4=np.mean(sigma4[iref])
        sigma_idx5=np.mean(sigma5[iref])
        sigma_mag1=np.mean(sigma6[iref])
        sigma_mag2=np.mean(sigma7[iref])
        sigma_mag3=np.mean(sigma8[iref])
        sigma_mag4=np.mean(sigma9[iref])
    else:
        sigma_idx1=np.sqrt(np.mean((idx1[iref]-idx1_ref)**2))
        sigma_idx2=np.sqrt(np.mean((idx2[iref]-idx2_ref)**2))
        sigma_idx3=np.sqrt(np.mean((idx3[iref]-idx3_ref)**2))
        sigma_idx4=np.sqrt(np.mean((idx4[iref]-idx4_ref)**2))
        sigma_idx5=np.sqrt(np.mean((idx5[iref]-idx5_ref)**2))
        
        #sigma_idx1=(f_plt.perc_84(idx1[iref])-f_plt.perc_16(idx1[iref]))/2.0
        #sigma_idx2=(f_plt.perc_84(idx2[iref])-f_plt.perc_16(idx2[iref]))/2.0
        #sigma_idx3=(f_plt.perc_84(idx3[iref])-f_plt.perc_16(idx3[iref]))/2.0
        #sigma_idx4=(f_plt.perc_84(idx4[iref])-f_plt.perc_16(idx4[iref]))/2.0
        #sigma_idx5=(f_plt.perc_84(idx5[iref])-f_plt.perc_16(idx5[iref]))/2.0
        
        
        sigma_mag1=np.sqrt(np.mean((mag1[iref]-mag1_ref)**2))
        sigma_mag2=np.sqrt(np.mean((mag2[iref]-mag2_ref)**2))
        sigma_mag3=np.sqrt(np.mean((mag3[iref]-mag3_ref)**2))
        sigma_mag4=np.sqrt(np.mean((mag4[iref]-mag4_ref)**2))
        
        
        
        #sigma_mag1=(f_plt.perc_84(mag1[iref])-f_plt.perc_16(mag1[iref]))/2.0
        #sigma_mag2=(f_plt.perc_84(mag2[iref])-f_plt.perc_16(mag2[iref]))/2.0
        #sigma_mag3=(f_plt.perc_84(mag3[iref])-f_plt.perc_16(mag3[iref]))/2.0
        #sigma_mag4=(f_plt.perc_84(mag4[iref])-f_plt.perc_16(mag4[iref]))/2.0
    
    
    chi_q_idx=((idx1_sel-idx1_ref)/sigma_idx1)**2+((idx2_sel-idx2_ref)/sigma_idx2)**2+((idx3_sel-idx3_ref)/sigma_idx3)**2+((idx4_sel-idx4_ref)/sigma_idx4)**2+((idx5_sel-idx5_ref)/sigma_idx5)**2
    chi_q_mag=((mag1_sel-mag1_ref)/sigma_mag1)**2+((mag2_sel-mag2_ref)/sigma_mag2)**2+((mag3_sel-mag3_ref)/sigma_mag3)**2+((mag4_sel-mag4_ref)/sigma_mag4)**2
    
    chi_q_balmer=((idx1_sel-idx1_ref)/sigma_idx1)**2+((idx2_sel-idx2_ref)/sigma_idx2)**2+((idx3_sel-idx3_ref)/sigma_idx3)**2
    chi_q_mgfe=((idx4_sel-idx4_ref)/sigma_idx4)**2+((idx5_sel-idx5_ref)/sigma_idx5)**2

    chi_q=chi_q_idx+chi_q_mag
    
    bin_edges=np.histogram(par_sel, bins=bins)[1]
    bin_size=bin_edges[1]-bin_edges[0]
    
    for i_bin in range(1, bins, 2):
        bin_edges[i_bin]=bin_edges[i_bin]+bin_size*0.5
    
    #median_chi_idx=stats.binned_statistic(par_sel, chi_q_idx, statistic='median', bins=bins)
    median_chi_idx=f_plt.running_median(par_sel,chi_q_idx, bin_edges)
    p16_chi_idx=f_plt.running_perc(par_sel, chi_q_idx, bin_edges,16)
    p84_chi_idx=f_plt.running_perc(par_sel, chi_q_idx, bin_edges,84)

    p16_chi_balmer=f_plt.running_perc(par_sel, chi_q_balmer, bin_edges,16)
    p16_chi_mgfe=f_plt.running_perc(par_sel, chi_q_mgfe, bin_edges,16)

    #median_chi_mag=stats.binned_statistic(par_sel, chi_q_mag, statistic='median', bins=bins)
    median_chi_mag=f_plt.running_median(par_sel,chi_q_mag, bin_edges)
    p16_chi_mag=f_plt.running_perc(par_sel, chi_q_mag, bin_edges,16)
    p84_chi_mag=f_plt.running_perc(par_sel, chi_q_mag, bin_edges,84)
    
    #median_chi=stats.binned_statistic(par_sel, chi_q, statistic='median', bins=bins)
    median_chi=f_plt.running_median(par_sel,chi_q, bin_edges)
    p16_chi=f_plt.running_perc(par_sel, chi_q, bin_edges,16)
    p84_chi=f_plt.running_perc(par_sel, chi_q, bin_edges,84)
    
    idx_ref_chi=np.where(par_sel<val_ref)[0]
    median_chi_ref=np.median(chi_q[idx_ref_chi])
    #p84_chi_ref=f_plt.perc_84(chi_q[idx_ref_chi])
    #p16_chi_ref=f_plt.perc_16(chi_q[idx_ref_chi])
    
    median_chi_ref_idx=np.median(chi_q_idx[idx_ref_chi])
    #p84_chi_ref_idx=f_plt.perc_84(chi_q_idx[idx_ref_chi])
    #p16_chi_ref_idx=f_plt.perc_16(chi_q_idx[idx_ref_chi])
    
    median_chi_ref_mag=np.median(chi_q_mag[idx_ref_chi])
    #p84_chi_ref_mag=f_plt.perc_84(chi_q_mag[idx_ref_chi])
    #p16_chi_ref_mag=f_plt.perc_16(chi_q_mag[idx_ref_chi])
    
    #upper_lim=stats.chi2.ppf(0.84, 9)
    #upper_lim_idx=stats.chi2.ppf(0.84, 5)
    #upper_lim_mag=stats.chi2.ppf(0.84,4)
    
    upper_lim=np.percentile(chi_q[idx_ref_chi],84)
    upper_lim_idx=np.percentile(chi_q_idx[idx_ref_chi],84)
    upper_lim_mag=np.percentile(chi_q_mag[idx_ref_chi],84)
    
    upper_lim_balmer=np.percentile(chi_q_balmer[idx_ref_chi],84)
    upper_lim_mgfe=np.percentile(chi_q_mgfe[idx_ref_chi],84)

    x=[0.0]*bins
    
    for i in range(0,bins):
        x[i]=(bin_edges[i+1]+bin_edges[i])/2.0
    
    if mkplot:
        fig, axs=plt.subplots(1,3,figsize=figsize)
        axs[0].scatter(par_sel,(chi_q),s=1)
        axs[0].plot(x, (median_chi), color='red')
        axs[0].plot(x, p16_chi, color='orange')
        axs[0].plot(x, p84_chi, color='orange')
        
        axs[0].plot([np.min(par_sel), np.max(par_sel)], [median_chi_ref, median_chi_ref], color='red')
        axs[0].plot([np.min(par_sel), np.max(par_sel)], [upper_lim, upper_lim],color='orange')
        #axs[0].plot([np.min(par_sel), np.max(par_sel)],[p16_chi_ref,p16_chi_ref], color='#ff028d')
        axs[0].set_yscale("log")
        
        axs[0].set_ylim(ylim)
        
        axs[0].set_xlabel('d1090n50')
        axs[0].set_ylabel('chi_q')
        
        axs[1].scatter(par_sel,(chi_q_idx),s=1)
        axs[1].plot(x, (median_chi_idx), color='red')
        axs[1].plot(x, p16_chi_idx, color='orange')
        axs[1].plot(x, p84_chi_idx, color='orange')
        
        axs[1].plot([np.min(par_sel), np.max(par_sel)], [median_chi_ref_idx, median_chi_ref_idx], color='red')
        axs[1].plot([np.min(par_sel), np.max(par_sel)], [upper_lim_idx, upper_lim_idx],color='orange')
       # axs[1].plot([np.min(par_sel), np.max(par_sel)],[p16_chi_ref_idx,p16_chi_ref_idx], color='#ff028d')
        
    
        axs[1].set_xlabel('d1090n50')
        axs[1].set_ylabel('chi_q_idx')
        
        axs[1].set_yscale("log")
        axs[1].set_ylim(ylim)
    
        
        axs[2].scatter(par_sel,(chi_q_mag),s=1)
        axs[2].plot(x, (median_chi_mag), color='red')
        axs[2].plot(x, p16_chi_mag, color='orange')
        axs[2].plot(x, p84_chi_mag, color='orange')
        
        axs[2].plot([np.min(par_sel), np.max(par_sel)], [median_chi_ref_mag, median_chi_ref_mag], color='red')
        axs[2].plot([np.min(par_sel), np.max(par_sel)], [upper_lim_mag, upper_lim_mag],color='orange')
        #axs[2].plot([np.min(par_sel), np.max(par_sel)],[p16_chi_ref_mag,p16_chi_ref_mag], color='#ff028d')
        
        axs[2].set_yscale("log")
        axs[2].set_ylim(ylim)
        
        axs[2].set_xlabel('d1090n50')
        axs[2].set_ylabel('chi_q_col')
        
        axs[0].set_facecolor('#d8dcd6')
        axs[1].set_facecolor('#d8dcd6')
        axs[2].set_facecolor('#d8dcd6')
    
        axs[1].set_title(title)
    
    #inter_max=np.interp(xmax,x, median_chi.statistic)
    #inter_min=np.interp(xmin, x, median_chi.statistic)
    inter_chi=lambda t: np.interp(t, x, p16_chi-upper_lim)
    inter_chi_idx=lambda t: np.interp(t, x, p16_chi_idx-upper_lim_idx)
    inter_chi_mag=lambda t: np.interp(t, x, p16_chi_mag-upper_lim_mag)
    inter_chi_balmer=lambda t: np.interp(t, x, p16_chi_balmer-upper_lim_balmer)
    inter_chi_mgfe=lambda t: np.interp(t, x, p16_chi_mgfe-upper_lim_mgfe)
    
    xmax_tot=-1.0
    while((inter_chi(xmin)*inter_chi(xmax_tot)>0) and (xmax_tot<1.0)):
        xmax_tot=xmax_tot+0.05
        
    xmax_idx=-1.0
    while((inter_chi_idx(xmin)*inter_chi_idx(xmax_idx)>0) and (xmax_idx<1.0)):
        xmax_idx=xmax_idx+0.05
        
    xmax_mag=-1.0
    while((inter_chi_mag(xmin)*inter_chi_mag(xmax_mag)>0) and (xmax_mag<1.0)):
        xmax_mag=xmax_mag+0.05
    
    x_m=f_plt.bisection(inter_chi, xmin, xmax_tot, toll)
    x_m_idx=f_plt.bisection(inter_chi_idx, xmin, xmax, toll)
    x_m_mag=f_plt.bisection(inter_chi_mag, xmin, xmax, toll)
    x_m_balmer=f_plt.bisection(inter_chi_balmer, xmin, xmax, toll)
    x_m_mgfe=f_plt.bisection(inter_chi_mgfe, xmin, xmax, toll)

    if mkplot:
        axs[0].scatter(x_m, inter_chi(x_m)+upper_lim, s=20, color='#ff028d')
    
    print('d1090n50 limit tot:', x_m)
    print('d1090n50 limit idx:', x_m_idx)
    print('d1090n50 limit col:', x_m_mag)

    print('sigma_idx1:', sigma_idx1)
    print('sigma_idx2:', sigma_idx2)
    print('sigma_idx3:', sigma_idx3)
    print('sigma_idx4:', sigma_idx4)
    print('sigma_idx5:', sigma_idx5)
    print('sigma_col1:', sigma_mag1)
    print('sigma_col2:', sigma_mag2)
    print('sigma_col3:', sigma_mag3)
    print('sigma_col4:', sigma_mag4)
    
    if mkplot:
        return fig, x_m
    else:
        return chi_q




def chi_q_comp_idx(par,idx1, idx2, idx3, idx4, idx5,isel, iref,val_ref=-1.00, name_par='d1090n50', name_idx=['D4000n','HdHg','Hb','Mg2Fe','MgFep'],figsize=(15,10), title='', ylim=[None,None], xlim=[None, None], bins=50):
    
    import function_plot as f_plt
    fig, axs=plt.subplots(2,3,figsize=figsize)    
    axs[1,2].axis('off')
    
    par_sel=par[isel]
    
    idx1_sel=idx1[isel]
    idx2_sel=idx2[isel]
    idx3_sel=idx3[isel]
    idx4_sel=idx4[isel]
    idx5_sel=idx5[isel]
    
    idx1_ref=np.mean(idx1[iref])
    idx2_ref=np.mean(idx2[iref])
    idx3_ref=np.mean(idx3[iref])
    idx4_ref=np.mean(idx4[iref])
    idx5_ref=np.mean(idx5[iref])

   
    #sigma_idx1=(f_plt.perc_84(idx1[iref])-f_plt.perc_16(idx1[iref]))/2.0
    #sigma_idx2=(f_plt.perc_84(idx2[iref])-f_plt.perc_16(idx2[iref]))/2.0
    #sigma_idx3=(f_plt.perc_84(idx3[iref])-f_plt.perc_16(idx3[iref]))/2.0
    #sigma_idx4=(f_plt.perc_84(idx4[iref])-f_plt.perc_16(idx4[iref]))/2.0
    #sigma_idx5=(f_plt.perc_84(idx5[iref])-f_plt.perc_16(idx5[iref]))/2.0
    
    sigma_idx1=np.sqrt(np.mean((idx1[iref]-idx1_ref)**2))
    sigma_idx2=np.sqrt(np.mean((idx2[iref]-idx2_ref)**2))
    sigma_idx3=np.sqrt(np.mean((idx3[iref]-idx3_ref)**2))
    sigma_idx4=np.sqrt(np.mean((idx4[iref]-idx4_ref)**2))
    sigma_idx5=np.sqrt(np.mean((idx5[iref]-idx5_ref)**2))

    chi_q1=((idx1_sel-idx1_ref)/sigma_idx1)**2
    chi_q2=((idx2_sel-idx2_ref)/sigma_idx2)**2
    chi_q3=((idx3_sel-idx3_ref)/sigma_idx3)**2
    chi_q4=((idx4_sel-idx4_ref)/sigma_idx4)**2
    chi_q5=((idx5_sel-idx5_ref)/sigma_idx5)**2
    
    #axs[1,2].scatter(par_sel, idx1_sel, s=1)
    print('sigma idx1:', sigma_idx1)
    #print('sigma idx1 p:', sigma_idx1_p)
    print('sigma idx2:', sigma_idx2)
    #print('sigma idx2 p:', sigma_idx2_p)
    print('sigma idx3:', sigma_idx3)
    #print('sigma idx3 p:', sigma_idx3_p)
    print('sigma idx4:', sigma_idx4)
    #print('sigma idx4 p:', sigma_idx4_p)
    print('sigma idx5:', sigma_idx5)
    #print('sigma idx5 p:', sigma_idx5_p)
    idx_ref_chi=np.where(par_sel<val_ref)[0]
    
    median_chi_ref1=np.median(chi_q1[idx_ref_chi])
    p84_chi_ref1=f_plt.perc_84(chi_q1[idx_ref_chi])
    #p16_chi_ref1=f_plt.perc_16(chi_q1[idx_ref_chi])

    median_chi_ref2=np.median(chi_q2[idx_ref_chi])
    p84_chi_ref2=f_plt.perc_84(chi_q2[idx_ref_chi])
    #p16_chi_ref2=f_plt.perc_16(chi_q2[idx_ref_chi])
    
    median_chi_ref3=np.median(chi_q3[idx_ref_chi])
    p84_chi_ref3=f_plt.perc_84(chi_q3[idx_ref_chi])
    #p16_chi_ref3=f_plt.perc_16(chi_q3[idx_ref_chi])
    
    median_chi_ref4=np.median(chi_q4[idx_ref_chi])
    p84_chi_ref4=f_plt.perc_84(chi_q4[idx_ref_chi])
    #p16_chi_ref4=f_plt.perc_16(chi_q4[idx_ref_chi])
    
    median_chi_ref5=np.median(chi_q5[idx_ref_chi])
    p84_chi_ref5=f_plt.perc_84(chi_q5[idx_ref_chi])
    #p16_chi_ref5=f_plt.perc_16(chi_q5[idx_ref_chi])

    #median_chi_idx1=stats.binned_statistic(par_sel, chi_q1, statistic='median', bins=bins)
    #median_chi_idx2=stats.binned_statistic(par_sel, chi_q2, statistic='median', bins=bins)
    #median_chi_idx3=stats.binned_statistic(par_sel, chi_q3, statistic='median', bins=bins)
    #median_chi_idx4=stats.binned_statistic(par_sel, chi_q4, statistic='median', bins=bins)
    #median_chi_idx5=stats.binned_statistic(par_sel, chi_q5, statistic='median', bins=bins)
    
    bin_edges=np.histogram(par_sel, bins=bins)[1]
    
    median_chi_idx1=f_plt.running_median(par_sel,chi_q1, bin_edges)
    median_chi_idx2=f_plt.running_median(par_sel,chi_q2, bin_edges)
    median_chi_idx3=f_plt.running_median(par_sel,chi_q3, bin_edges)
    median_chi_idx4=f_plt.running_median(par_sel,chi_q4, bin_edges)
    median_chi_idx5=f_plt.running_median(par_sel,chi_q5, bin_edges)
    
    x=[0.0]*bins
    
    for i in range(0,bins):
        x[i]=(bin_edges[i+1]+bin_edges[i])/2.0

    axs[0,0].scatter(par_sel, chi_q1, s=1)
    axs[0,0].plot(x, median_chi_idx1, color='red')
    axs[0,0].plot([np.min(par_sel), np.max(par_sel)], [median_chi_ref1, median_chi_ref1], color='orange')
    axs[0,0].plot([np.min(par_sel), np.max(par_sel)], [p84_chi_ref1, p84_chi_ref1],color='#ff028d')
    #axs[0,0].plot([np.min(par_sel), np.max(par_sel)],[p16_chi_ref1,p16_chi_ref1], color='#ff028d')
    axs[0,0].set_facecolor('#d8dcd6')
    axs[0,0].set_xlabel(name_par)
    axs[0,0].set_ylabel('chi_q '+name_idx[0])
    axs[0,0].set_ylim(ylim)
    axs[0,0].set_yscale("log")

    axs[0,1].scatter(par_sel, chi_q2, s=1)
    axs[0,1].plot(x, median_chi_idx2, color='red')
    axs[0,1].plot([np.min(par_sel), np.max(par_sel)], [median_chi_ref2, median_chi_ref2], color='orange')
    axs[0,1].plot([np.min(par_sel), np.max(par_sel)], [p84_chi_ref2, p84_chi_ref2],color='#ff028d')
    #axs[0,1].plot([np.min(par_sel), np.max(par_sel)],[p16_chi_ref2,p16_chi_ref2], color='#ff028d')
    axs[0,1].set_facecolor('#d8dcd6')
    axs[0,1].set_xlabel(name_par)
    axs[0,1].set_ylabel('chi_q '+name_idx[1])
    axs[0,1].set_ylim(ylim)
    axs[0,1].set_title(title)
    axs[0,1].set_yscale("log")

    axs[0,2].scatter(par_sel, chi_q3, s=1)
    axs[0,2].plot(x, median_chi_idx3, color='red')
    axs[0,2].plot([np.min(par_sel), np.max(par_sel)], [median_chi_ref3, median_chi_ref3], color='orange')
    axs[0,2].plot([np.min(par_sel), np.max(par_sel)], [p84_chi_ref3, p84_chi_ref3],color='#ff028d')
    #axs[0,2].plot([np.min(par_sel), np.max(par_sel)],[p16_chi_ref3,p16_chi_ref3], color='#ff028d')
    axs[0,2].set_facecolor('#d8dcd6')
    axs[0,2].set_xlabel(name_par)
    axs[0,2].set_ylabel('chi_q '+name_idx[2])
    axs[0,2].set_ylim(ylim)
    axs[0,2].set_yscale("log")
    
    axs[1,0].scatter(par_sel, chi_q4, s=1)
    axs[1,0].plot(x, median_chi_idx4, color='red')
    axs[1,0].plot([np.min(par_sel), np.max(par_sel)], [median_chi_ref4, median_chi_ref4], color='orange')
    axs[1,0].plot([np.min(par_sel), np.max(par_sel)], [p84_chi_ref4, p84_chi_ref4],color='#ff028d')
   # axs[1,0].plot([np.min(par_sel), np.max(par_sel)],[p16_chi_ref4,p16_chi_ref4], color='#ff028d')
    axs[1,0].set_facecolor('#d8dcd6')
    axs[1,0].set_xlabel(name_par)
    axs[1,0].set_ylabel('chi_q '+name_idx[3])
    axs[1,0].set_ylim(ylim)
    axs[1,0].set_yscale("log")
    
    axs[1,1].scatter(par_sel, chi_q5, s=1)
    axs[1,1].plot(x, median_chi_idx5, color='red')
    axs[1,1].plot([np.min(par_sel), np.max(par_sel)], [median_chi_ref5, median_chi_ref5], color='orange')
    axs[1,1].plot([np.min(par_sel), np.max(par_sel)], [p84_chi_ref5, p84_chi_ref5],color='#ff028d')
    #axs[1,1].plot([np.min(par_sel), np.max(par_sel)],[p16_chi_ref5,p16_chi_ref5], color='#ff028d')
    axs[1,1].set_facecolor('#d8dcd6')
    axs[1,1].set_xlabel(name_par)
    axs[1,1].set_ylabel('chi_q '+name_idx[4])
    axs[1,1].set_ylim(ylim)
    axs[1,1].set_yscale("log")
    
    return fig

def chi_q_comp_col(par,idx1, idx2, idx3, idx4,isel, iref,val_ref=-1.00, name_par='d1090n50', name_idx=['u-r','g-r','r-i','r-z'],figsize=(15,10), title='', ylim=[None,None], xlim=[None, None], bins=50):
    
    import function_plot as f_plt
    fig, axs=plt.subplots(2,3,figsize=figsize)    
    axs[1,2].axis('off')
    axs[1,1].axis('off')
    
    par_sel=par[isel]
    
    idx1_sel=idx1[isel]
    idx2_sel=idx2[isel]
    idx3_sel=idx3[isel]
    idx4_sel=idx4[isel]
        
    idx1_ref=np.mean(idx1[iref])
    idx2_ref=np.mean(idx2[iref])
    idx3_ref=np.mean(idx3[iref])
    idx4_ref=np.mean(idx4[iref])

   
    sigma_idx1=np.sqrt(np.mean((idx1[iref]-idx1_ref)**2))
    sigma_idx2=np.sqrt(np.mean((idx2[iref]-idx2_ref)**2))
    sigma_idx3=np.sqrt(np.mean((idx3[iref]-idx3_ref)**2))
    sigma_idx4=np.sqrt(np.mean((idx4[iref]-idx4_ref)**2))
    
    #sigma_idx1=(f_plt.perc_84(idx1[iref])-f_plt.perc_16(idx1[iref]))/2.0
    #sigma_idx2=(f_plt.perc_84(idx2[iref])-f_plt.perc_16(idx2[iref]))/2.0
    #sigma_idx3=(f_plt.perc_84(idx3[iref])-f_plt.perc_16(idx3[iref]))/2.0
    #sigma_idx4=(f_plt.perc_84(idx4[iref])-f_plt.perc_16(idx4[iref]))/2.0

    chi_q1=((idx1_sel-idx1_ref)/sigma_idx1)**2
    chi_q2=((idx2_sel-idx2_ref)/sigma_idx2)**2
    chi_q3=((idx3_sel-idx3_ref)/sigma_idx3)**2
    chi_q4=((idx4_sel-idx4_ref)/sigma_idx4)**2
    
    idx_ref_chi=np.where(par_sel<val_ref)[0]
    
    median_chi_ref1=np.median(chi_q1[idx_ref_chi])
    p84_chi_ref1=f_plt.perc_84(chi_q1[idx_ref_chi])
    #p16_chi_ref1=f_plt.perc_16(chi_q1[idx_ref_chi])

    median_chi_ref2=np.median(chi_q2[idx_ref_chi])
    p84_chi_ref2=f_plt.perc_84(chi_q2[idx_ref_chi])
    #p16_chi_ref2=f_plt.perc_16(chi_q2[idx_ref_chi])
    
    median_chi_ref3=np.median(chi_q3[idx_ref_chi])
    p84_chi_ref3=f_plt.perc_84(chi_q3[idx_ref_chi])
    #p16_chi_ref3=f_plt.perc_16(chi_q3[idx_ref_chi])
    
    median_chi_ref4=np.median(chi_q4[idx_ref_chi])
    p84_chi_ref4=f_plt.perc_84(chi_q4[idx_ref_chi])
    #p16_chi_ref4=f_plt.perc_16(chi_q4[idx_ref_chi])
    

   # median_chi_idx1=stats.binned_statistic(par_sel, chi_q1, statistic='median', bins=bins)
    #median_chi_idx2=stats.binned_statistic(par_sel, chi_q2, statistic='median', bins=bins)
    #median_chi_idx3=stats.binned_statistic(par_sel, chi_q3, statistic='median', bins=bins)
    #median_chi_idx4=stats.binned_statistic(par_sel, chi_q4, statistic='median', bins=bins)
    
    bin_edges=np.histogram(par_sel, bins=bins)[1]
    
    median_chi_idx1=f_plt.running_median(par_sel,chi_q1, bin_edges)
    median_chi_idx2=f_plt.running_median(par_sel,chi_q2, bin_edges)
    median_chi_idx3=f_plt.running_median(par_sel,chi_q3, bin_edges)
    median_chi_idx4=f_plt.running_median(par_sel,chi_q4, bin_edges)
    
    x=[0.0]*bins
    
    for i in range(0,bins):
        x[i]=(bin_edges[i+1]+bin_edges[i])/2.0

    axs[0,0].scatter(par_sel, chi_q1, s=1)
    axs[0,0].plot(x, median_chi_idx1, color='red')
    axs[0,0].plot([np.min(par_sel), np.max(par_sel)], [median_chi_ref1, median_chi_ref1], color='orange')
    axs[0,0].plot([np.min(par_sel), np.max(par_sel)], [p84_chi_ref1, p84_chi_ref1],color='#ff028d')
    #axs[0,0].plot([np.min(par_sel), np.max(par_sel)],[p16_chi_ref1,p16_chi_ref1], color='#ff028d')
    axs[0,0].set_facecolor('#d8dcd6')
    axs[0,0].set_xlabel(name_par)
    axs[0,0].set_ylabel('chi_q '+name_idx[0])
    axs[0,0].set_ylim(ylim)
    axs[0,0].set_yscale("log")

    axs[0,1].scatter(par_sel, chi_q2, s=1)
    axs[0,1].plot(x, median_chi_idx2, color='red')
    axs[0,1].plot([np.min(par_sel), np.max(par_sel)], [median_chi_ref2, median_chi_ref2], color='orange')
    axs[0,1].plot([np.min(par_sel), np.max(par_sel)], [p84_chi_ref2, p84_chi_ref2],color='#ff028d')
    #axs[0,1].plot([np.min(par_sel), np.max(par_sel)],[p16_chi_ref2,p16_chi_ref2], color='#ff028d')
    axs[0,1].set_facecolor('#d8dcd6')
    axs[0,1].set_xlabel(name_par)
    axs[0,1].set_ylabel('chi_q '+name_idx[1])
    axs[0,1].set_ylim(ylim)
    axs[0,1].set_title(title)
    axs[0,1].set_yscale("log")

    axs[0,2].scatter(par_sel, chi_q3, s=1)
    axs[0,2].plot(x, median_chi_idx3, color='red')
    axs[0,2].plot([np.min(par_sel), np.max(par_sel)], [median_chi_ref3, median_chi_ref3], color='orange')
    axs[0,2].plot([np.min(par_sel), np.max(par_sel)], [p84_chi_ref3, p84_chi_ref3],color='#ff028d')
    #axs[0,2].plot([np.min(par_sel), np.max(par_sel)],[p16_chi_ref3,p16_chi_ref3], color='#ff028d')
    axs[0,2].set_facecolor('#d8dcd6')
    axs[0,2].set_xlabel(name_par)
    axs[0,2].set_ylabel('chi_q '+name_idx[2])
    axs[0,2].set_ylim(ylim)
    axs[0,2].set_yscale("log")
    
    axs[1,0].scatter(par_sel, chi_q4, s=1)
    axs[1,0].plot(x, median_chi_idx4, color='red')
    axs[1,0].plot([np.min(par_sel), np.max(par_sel)], [median_chi_ref4, median_chi_ref4], color='orange')
    axs[1,0].plot([np.min(par_sel), np.max(par_sel)], [p84_chi_ref4, p84_chi_ref4],color='#ff028d')
    #axs[1,0].plot([np.min(par_sel), np.max(par_sel)],[p16_chi_ref4,p16_chi_ref4], color='#ff028d')
    axs[1,0].set_facecolor('#d8dcd6')
    axs[1,0].set_xlabel(name_par)
    axs[1,0].set_ylabel('chi_q '+name_idx[3])
    axs[1,0].set_ylim(ylim)
    axs[1,0].set_yscale("log")
    
    
    return fig
 

def running_median(par,values,bin_edges):
    
    N_edges=np.size(bin_edges)
    
    median=np.array([0.0]*(N_edges-1))
    
    
    for i_bin in range(0, N_edges-1):
        _i_sel=((par>bin_edges[i_bin])&(par<(bin_edges[i_bin+1])))
        median[i_bin]=np.median(values[_i_sel])
        
    return median
        
def running_perc(par,values,bin_edges,perc):
    
    N_edges=np.size(bin_edges)
    
    percentile=np.array([0.0]*(N_edges-1))
    
    
    for i_bin in range(0, N_edges-1):
        _i_sel=((par>bin_edges[i_bin])&(par<(bin_edges[i_bin+1])))
        percentile[i_bin]=np.percentile(values[_i_sel], perc)
        
    return percentile
                
        
def scatter_hist(par1, par2, bins=50, figsize=(10,10), s=1, name_par=''):
    from mpl_toolkits.axes_grid1 import make_axes_locatable

    
    fig, axs=plt.subplots(figsize=figsize)
    axs.scatter((par1), (par2), s=s)
    divider=make_axes_locatable(axs)
    histx=divider.append_axes('top', 1.5, pad=0.1, sharex=axs)
    histy=divider.append_axes('right', 1.5, pad=0.1, sharey=axs)
    histx.xaxis.set_tick_params(labelbottom=False)
    histy.yaxis.set_tick_params(labelleft=False)
    
    histx.hist((par1), bins=bins)
    histy.hist((par2), bins=bins, orientation='horizontal')
    axs.set_xlabel(name_par+r'$_{in}$', size=16)
    axs.set_ylabel(name_par+r'$_{out}$', size=16)
    
    return fig

def ellipse(a,b,phi=0,x_c=0.0,y_c=0.0):
    import numpy as np
    
    theta=np.arange(0, 2*np.pi, 0.01)
    
    xpos=a*np.cos(theta)
    ypos=b*np.sin(theta)
    
    newx=xpos*np.cos(phi)+ypos*np.sin(phi)+x_c
    newy=-xpos*np.sin(phi)+ypos*np.cos(phi)+y_c
    
    return newx, newy
    
    
    
    
    
    
           
        