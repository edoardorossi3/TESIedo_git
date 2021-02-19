pro idx_specrebin
 
 common idxinfo
 read_idxinfo
 ;prog_dir='~/Desktop/TESIedo_git/'
 fits_dir='~/Desktop/TESI/models/sandage_dcomb/'
 prefix_file='sandage_varZ_v4.1eq_spec_'
 suffix_spec='dcomb_'
 suffix_idx='dcomb_idx_'
 suffix_idx_werr='dcomb_idx_s_werr'
 file_ext='.fits'
 sndg_num=indgen(2)
 tot_file=n_elements(sndg_num)
 
 SNR20=20.0
 dwl=1.0
 data_row=create_struct('IDX', 0L) ; !values.f_nan is a float

 ;for i=1, tot_file-1 do begin
  i=1
  
  ;data_row=create_struct('IDX', !values.f_nan, 'D4000n_0', !values.f_nan, 'D4000n_200', !values.f_nan,'D4000n20',!values.f_nan, 'err20_D4000n', !values.f_nan, $
  ;                       'HdHg_0', !values.f_nan, 'HdHg_200', !values.f_nan, 'HdHg20', !values.f_nan, 'err20_HdHg', !values.f_nan,  $
  ;                       'Lick_Hb_0', !values.f_nan, 'Lick_Hb_200', !values.f_nan, 'Lick_Hb20', !values.f_nan, 'err20_Lick_Hb', !values.f_nan,  $
  ;                       'Mg2Fe_0', !values.f_nan, 'Mg2Fe_200', !values.f_nan, 'Mg2Fe20', !values.f_nan, 'err20_Mg2Fe', !values.f_nan,  $
  ;                       'MgFep_0', !values.f_nan, 'MgFep_200', !values.f_nan, 'MgFep20', !values.f_nan, 'err20_MgFep', !values.f_nan)
  ;data_table=replicate(data_row, 12500)

  spec_file=fits_dir+prefix_file+suffix_spec+string(sndg_num[i], format='(I3.3)')+file_ext
  idx_file=fits_dir+prefix_file+suffix_idx+string(sndg_num[i], format='(I3.3)')+file_ext
  idx_werr_file=fits_dir+prefix_file+suffix_idx_werr+string(sndg_num[i], format='(I3.3)')+file_ext
  
  idx_table=mrdfits(idx_file,1)
  SPEC=mrdfits(spec_file, 1)
  wl=mrdfits(spec_file, 0)
  
  linenames=indexpars.index
  idx_names=tag_names(idx_table)
  Nidx_names=n_elements(idx_names)
  for j=1, Nidx_names-1 do begin
    data_row=create_struct(data_row, idx_names[j]+'_0', !values.f_nan, idx_names[j]+'_200', !values.f_nan, idx_names[j]+'20', !values.f_nan, 'err20_'+idx_names[j], !values.f_nan)
 
  endfor
 ;stop
  ;idx_model=0
 
  n_model=(size(SPEC))[2]
  idx_lim=where(wl gt 3700. and wl lt 8860., nlim_wl)
  wl_maxmin=[wl[idx_lim[0]], wl[idx_lim[-1]]]
  
  colnames=tag_names(data_row)
  ;stop
  data_table=replicate(data_row, n_model)
  data_table.idx=idx_table.idx
  
  
  ;stop
  

  for idx_model=0, n_model-1 do begin
    
    for j=1, Nidx_names-1 do begin
      col_num0=(where(strcmp(colnames,idx_names[j]+'_0',/fold_case)))[0]
      col_num200=(where(strcmp(colnames,idx_names[j]+'_200',/fold_case)))[0]
      tmp=idx_table[idx_model].(j)
      data_table[idx_model].(col_num0)=tmp[0]
      data_table[idx_model].(col_num200)=tmp[4]

          
    endfor
  
 ;mwrfits, data_table, idx_werr_file, /create
  ;stop
  
  
   ;model_num=12500*i+idx_model
   flx=flux_integral(wl, SPEC[*,idx_model], wl_maxmin)
   mean_flx=flx/(wl_maxmin[1]-wl_maxmin[0])
   err20=mean_flx/SNR20
 
   n_wlreb=fix((wl_maxmin[1]-wl_maxmin[0])/dwl)+1
   wl_rebin=findgen(n_wlreb)*dwl+wl_maxmin[0]
   idx_lim=where(wl gt wl_maxmin[0] and wl lt wl_maxmin[1], nlim_wl)
 
   SPEC_interp=interpol(SPEC[*,idx_model], wl, wl_rebin, /spline)

 
 
   err20_arr=fltarr(n_elements(wl_rebin))+err20
  ; stop
   
   
   for k=1, Nidx_names-1 do begin
     if (where(strcmp(linenames, idx_names[k], /fold_case)) ne -1) then begin
        col_num20=(where(strcmp(idx_names[k]+'20', colnames, /fold_case)))
        col_numerr20=(where(strcmp('err20_'+idx_names[k], colnames, /fold_case)))
        ;tmp=(idx_table[idx_model].(k))[0]
        tmp_arr=line_index(wl_rebin, SPEC_interp, err20_arr, indexpars[where(strcmp(indexpars.index, idx_names[k], /fold_case))])
        data_table[idx_model].(col_num20)=tmp_arr[0]
        data_table[idx_model].(col_numerr20)=tmp_arr[1]
     endif
     
     if (where(strcmp(idx_names[k], 'd4000n', /fold_case)) ne -1) then begin
        col_num20=(where(strcmp(idx_names[k]+'20', colnames, /fold_case)))
        col_numerr20=(where(strcmp('err20_'+idx_names[k], colnames, /fold_case)))
        tmp_arr=d4000_index(wl_rebin, SPEC_interp, err20_arr, narrow=1)
        data_table[idx_model].(col_num20)=tmp_arr[0]
        data_table[idx_model].(col_numerr20)=tmp_arr[1]
     endif
     
     if (where(strcmp(idx_names[k], 'd4000', /fold_case)) ne -1) then begin
       col_num20=(where(strcmp(idx_names[k]+'20', colnames, /fold_case)))
       col_numerr20=(where(strcmp('err20_'+idx_names[k], colnames, /fold_case)))
       tmp_arr=d4000_index(wl_rebin, SPEC_interp, err20_arr, narrow=0)
       data_table[idx_model].(col_num20)=tmp_arr[0]
       data_table[idx_model].(col_numerr20)=tmp_arr[1]
     endif 
   
     ;stop
    ; if ((where(strcmp(linenames, idx_names[k], /fold_case)) eq -1) and (where(strcmp(idx_names[k], 'd4000n', /fold_case)) eq -1) and $
     ;     (where(strcmp(idx_names[k], 'd4000', /fold_case)) eq -1)) then begin
     if  ((where(strcmp(idx_names[k],'MgFe_prime',  /fold_case)) ne -1) or (where(strcmp(idx_names[k], 'mg2fe', /fold_case)) ne -1) or $
            (where(strcmp(idx_names[k], 'HdHg', /fold_case)) ne -1)) then begin
       col_num20=(where(strcmp(idx_names[k]+'20', colnames, /fold_case)))
       col_numerr20=(where(strcmp('err20_'+idx_names[k], colnames, /fold_case)))
       tmp_arr=comp_index(wl_rebin, SPEC_interp, err20_arr, idx_names[k])
       data_table[idx_model].(col_num20)=tmp_arr[0]
       data_table[idx_model].(col_numerr20)=tmp_arr[1]
     endif
     
   endfor
  stop
   
  stop 
  endfor
 
 mwrfits, data_table, idx_werr_file, /create
 stop
 ;endfor 
 
 ;stop
 
end

pro check_interp
  
  common idxinfo
  read_idxinfo
  
  fits_dir='~/Desktop/TESI/models/sandage_dcomb/'
  prefix_file='sandage_varZ_v4.1eq_spec_'
  suffix_spec='dcomb090n_'
  suffix_idx='dcomb090n_idx_'
  suffix_idx_werr='dcomb090n_idx_werr'
  file_ext='.fits'
  sndg_num=1
  tot_file=n_elements(sndg_num)

  SNR=20.0
  dwl=1.0
  idx_model=0
  spec_file=fits_dir+prefix_file+suffix_spec+string(sndg_num, format='(I3.3)')+file_ext
  idx_file=fits_dir+prefix_file+suffix_idx+string(sndg_num, format='(I3.3)')+file_ext
  idx_werr_file=fits_dir+prefix_file+suffix_idx_werr+string(sndg_num, format='(I3.3)')+file_ext
  
  idx_table=mrdfits(idx_file, 1)
  d4000n_0=(idx_table[idx_model].d4000n)[0]
  hdhg_0=(idx_table[idx_model].hdhg)[0]
  mg2fe_0=(idx_table[idx_model].mg2fe)[0]
  mgfep_0=(idx_table[idx_model].mg2fe_prime)[0]
  SPEC=mrdfits(spec_file, 1)
  wl=mrdfits(spec_file, 0)
  ;n_model=(size(SPEC))[2]
  
  idx_lim=where(wl gt 3700. and wl lt 5600., nlim_wl)
  wl_maxmin=[wl[idx_lim[0]], wl[idx_lim[-1]]]
  n_wlreb=fix((wl_maxmin[1]-wl_maxmin[0])/dwl)+1
  wl_rebin=findgen(n_wlreb)*dwl+wl_maxmin[0]

  SPEC_interp_l=interpol(SPEC[*,idx_model], wl, wl_rebin)
  SPEC_interp_s=interpol(SPEC[*,idx_model], wl, wl_rebin, /spline)
  flx=flux_integral(wl, SPEC[*,idx_model], wl_maxmin)
  mean_flx=flx/(wl_maxmin[1]-wl_maxmin[0])
  err=mean_flx/SNR

  
   err_arr=fltarr(n_elements(wl_rebin))+err
   
   d4000n_arr_l=d4000_index(wl_rebin, SPEC_interp_l, err_arr, narrow=1)
   HdHg_arr_l=comp_index(wl_rebin, SPEC_interp_l, err_arr, 'HdHg')
   mg2fe_arr_l=comp_index(wl_rebin, SPEC_interp_l, err_arr, 'Mg2fe')
   mgfep_arr_l=comp_index(wl_rebin, SPEC_interp_l, err_arr, 'MgFep')
   lick_hb_arr_l=line_index(wl_rebin, SPEC_interp_l, err_arr, indexpars[where(strcmp(indexpars.index,'Lick_Hb'))])
   
   d4000n_arr_s=d4000_index(wl_rebin, SPEC_interp_s, err_arr, narrow=1)
   HdHg_arr_s=comp_index(wl_rebin, SPEC_interp_s, err_arr, 'HdHg')
   mg2fe_arr_s=comp_index(wl_rebin, SPEC_interp_s, err_arr, 'Mg2fe')
   mgfep_arr_s=comp_index(wl_rebin, SPEC_interp_s, err_arr, 'MgFep')
   lick_hb_arr_s=line_index(wl_rebin, SPEC_interp_s, err_arr, indexpars[where(strcmp(indexpars.index,'Lick_Hb'))])

   d4000n_l=d4000n_arr_l[0]
   err_d4000n_l=d4000n_arr_l[1]
   hdhg_l=HdHg_arr_l[0]
   err_hdhg_l=HdHg_arr_l[1]
   mg2fe_l=mg2fe_arr_l[0]
   err_mg2fe_l=mg2fe_arr_l[1]
   mgfep_l=mgfep_arr_l[0]
   err_mgfep_l=mgfep_arr_l[1]
   lick_hb_l=lick_hb_arr_l[0]
   err_lick_hb_l=lick_hb_arr_l[1]
   
   d4000n_s=d4000n_arr_s[0]
   err_d4000n_s=d4000n_arr_s[1]
   hdhg_s=HdHg_arr_s[0]
   err_hdhg_s=HdHg_arr_s[1]
   mg2fe_s=mg2fe_arr_s[0]
   err_mg2fe_s=mg2fe_arr_s[1]
   mgfep_s=mgfep_arr_s[0]
   err_mgfep_s=mgfep_arr_s[1]
   lick_hb_s=lick_hb_arr_s[0]
   err_lick_hb_s=lick_hb_arr_s[1]
   
   stop
   end
 