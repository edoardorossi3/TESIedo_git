pro idx_specrebin
 
 ;prog_dir='~/Desktop/TESIedo_git/'
 fits_dir='~/Desktop/TESI/sandage_dcomb/'
 name_file='sandage_varZ_v4.1eq_spec_dcomb090n_001.fits'
 fits_file=fits_dir+name_file
 
 idx_model=0
 SNR=20.0
 dwl=1.0
 wl_maxmin=[3700.,5600.] ;angstrom unit
 SPEC=mrdfits(fits_file, 1)
 wl=mrdfits(fits_file, 0)
 
 flx=flux_integral(wl, SPEC[*,idx_model], wl_maxmin,err_in=[])
 mean_flx=flx/(wl_maxmin[1]-wl_maxmin[0])
 err=mean_flx/SNR
 
 n_wlreb=fix((wl_maxmin[1]-wl_maxmin[0])/dwl)+1
 wl_rebin=findgen(n_wlreb)*dwl+wl_maxmin[0]
 idx_lim=where(wl gt wl_maxmin[0] and wl lt wl_maxmin[1], nlim_wl)
 
 
 SPEC_interp=interpol(SPEC[wl[idx_lim],idx_model], wl[idx_lim], wl_rebin)
 
 mean_flx_arr=fltarr(n_elements(wl[idx_lim]))+mean_flx
 pp1=plot( wl[idx_lim],spec[wl[idx_lim], 0])
 ;pp=plot(wl_rebin,spec_interp, color='green')
 pp2=plot(wl[idx_lim], mean_flx_arr,color='red', /overplot)
 
 
 stop
 
end 