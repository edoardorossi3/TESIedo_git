pro idx_specrebin
 
 common idxinfo
 read_idxinfo
 ;prog_dir='~/Desktop/TESIedo_git/'
 fits_dir='~/Desktop/TESI/models/sandage_dcomb/'
 prefix_file='sandage_varZ_v4.1eq_spec_'
 suffix_spec='dcomb090n_'
 suffix_idx='dcomb_idx_'
 suffix_idx_werr='dcomb_idx_werr'
 end_file='.fits'
 sndg_num=indgen(41)
 tot_file=n_elements(sndg_num)
 
 SNR=20.0
 dwl=1.0
 
 for i=1, tot_file-1 do begin
  
  data_row=create_struct('IDX', !values.f_nan, 'D4000n_0', !values.f_nan, 'D4000n_200', !values.f_nan,'D4000n',!values.f_nan, 'err_D4000n', !values.f_nan,$
                         'HdHg_0', !values.f_nan, 'HdHg_200', !values.f_nan, 'HdHg', !values.f_nan, 'err_HdHg', !values.f_nan,$
                         'Lick_Hb_0', !values.f_nan, 'Lick_Hb_200', !values.f_nan, 'Lick_Hb', !values.f_nan, 'err_Lick_Hb', !values.f_nan, $
                         'Mg2Fe_0', !values.f_nan, 'Mg2Fe_200', !values.f_nan, 'Mg2Fe', !values.f_nan, 'err_Mg2Fe', !values.f_nan, $
                         'MgFep_0', !values.f_nan, 'MgFep_200', !values.f_nan, 'MgFep', !values.f_nan, 'err_MgFep', !values.f_nan)
  data_table=replicate(data_row, 12500)
  
  spec_file=fits_dir+prefix_file+suffix_spec+string(sndg_num[i], format='(I3.3)')+end_file
  idx_file=fits_dir+prefix_file+suffix_idx+string(sndg_num[i], format='(I3.3)')+end_file
  idx_werr_file=fits_dir+prefix_file+suffix_idx_werr+string(sndg_num[i], format='(I3.3)')+end_file
 
  ;idx_model=0
 
  ;wl_maxmin=[3700.,5600.] ;angstrom unit
  SPEC=mrdfits(spec_file, 1)
  wl=mrdfits(spec_file, 0)
  n_model=(size(SPEC))[2]
  idx_lim=where(wl gt 3700. and wl lt 5600., nlim_wl)
  wl_maxmin=[wl[idx_lim[0]], wl[idx_lim[-1]]]
  
  colnames=tag_names(data_row)
  
  idx_table=mrdfits(idx_file,1)
  data_table.idx=idx_table.idx
  d4000table=idx_table.d4000n
  d4000n_0=d4000table[0,*]
  d4000n_200=d4000table[4,*]
  hdhgtable=idx_table.hdhg
  hdhg_0=hdhgtable[0,*]
  hdhg_200=hdhgtable[4,*]
  hbtable=idx_table.lick_hb
  lick_hb_0=hbtable[0,*]
  lick_hb_200=hbtable[4,*]
  mg2fetable=idx_table.mg2fe
  mg2fe_0=mg2fetable[0,*]
  mg2fe_200=mg2fetable[4,*]
  mgfeptable=idx_table.mgfe_prime
  mgfep_0=mgfeptable[0,*]
  mgfep_200=mgfeptable[4,*]
  
  
  ;col_idx=(where(strcmp(colnames,'idx',/fold_case)))[0]
  col_d4000_0=(where(strcmp(colnames,'D4000n_0',/fold_case)))[0]
  col_d4000_200=(where(strcmp(colnames,'D4000n_200',/fold_case)))[0]
  col_d4000=(where(strcmp(colnames,'D4000n',/fold_case)))[0]
  col_d4000_err=(where(strcmp(colnames,'err_D4000n',/fold_case)))[0]

  col_hdhg_0=(where(strcmp(colnames,'HdHg_0',/fold_case)))[0]
  col_hdhg_200=(where(strcmp(colnames,'HdHg_200',/fold_case)))[0]
  col_HdHg=(where(strcmp(colnames,'HdHg',/fold_case)))[0]
  col_hdhg_err=(where(strcmp(colnames,'err_HdHg',/fold_case)))[0]
  col_lick_hb_0=(where(strcmp(colnames,'Lick_hb_0',/fold_case)))[0]
  col_lick_hb_200=(where(strcmp(colnames,'Lick_hb_200',/fold_case)))[0]
  col_lick_hb=(where(strcmp(colnames,'lick_hb',/fold_case)))[0]
  col_lick_hb_err=(where(strcmp(colnames,'err_lick_hb',/fold_case)))[0]
  col_mg2fe_0=(where(strcmp(colnames,'Mg2fe_0',/fold_case)))[0]
  col_mg2fe_200=(where(strcmp(colnames,'Mg2fe_200',/fold_case)))[0]
  col_mg2fe=(where(strcmp(colnames,'Mg2Fe',/fold_case)))[0]
  col_mg2fe_err=(where(strcmp(colnames,'err_Mg2Fe',/fold_case)))[0]
  col_mgfep_0=(where(strcmp(colnames,'mgfep_0',/fold_case)))[0]
  col_mgfep_200=(where(strcmp(colnames,'mgfep_200',/fold_case)))[0]
  col_mgfep=(where(strcmp(colnames,'mgfep',/fold_case)))[0]
  col_mgfep_err=(where(strcmp(colnames,'err_mgfep',/fold_case)))[0]

  data_table[*].(col_d4000_0)=d4000n_0[*]
  data_table[*].(col_d4000_200)=d4000n_200[*]
  data_table[*].(col_hdhg_0)=hdhg_0[*]
  data_table[*].(col_hdhg_200)=hdhg_200[*]
  data_table[*].(col_lick_hb_0)=lick_hb_0[*]
  data_table[*].(col_lick_hb_200)=lick_hb_200[*]
  data_table[*].(col_mg2fe_0)=mg2fe_0[*]
  data_table[*].(col_mg2fe_200)=mg2fe_200[*]
  data_table[*].(col_mgfep_0)=mgfep_0[*]
  data_table[*].(col_mgfep_200)=mgfep_200[*]


  ;stop
  for idx_model=0, n_model-1 do begin
  
   ;model_num=12500*i+idx_model
   flx=flux_integral(wl, SPEC[*,idx_model], wl_maxmin)
   mean_flx=flx/(wl_maxmin[1]-wl_maxmin[0])
   err=mean_flx/SNR
 
   n_wlreb=fix((wl_maxmin[1]-wl_maxmin[0])/dwl)+1
   wl_rebin=findgen(n_wlreb)*dwl+wl_maxmin[0]
   idx_lim=where(wl gt wl_maxmin[0] and wl lt wl_maxmin[1], nlim_wl)
 
   SPEC_interp=interpol(SPEC[*,idx_model], wl, wl_rebin)
 
 
 
   err_arr=fltarr(n_elements(wl))+err
   ;mean_flx_arr=fltarr(n_elements(wl[idx_lim]))+mean_flx
   ;pp1=plot( wl_rebin,SPEC_interp)
   ;pp2=plot(wl[idx_lim], mean_flx_arr,color='red', /overplot)
   ;pp3=plot(wl[idx_lim], err_arr, color='blue', /overplot)
   d4000n_arr=d4000_index(wl, SPEC[*,idx_model], err_arr, narrow=1)
   HdHg_arr=comp_index(wl, SPEC[*,idx_model], err_arr, 'HdHg')
   mg2fe_arr=comp_index(wl, SPEC[*,idx_model], err_arr, 'Mg2fe')
   mgfep_arr=comp_index(wl, SPEC[*,idx_model], err_arr, 'MgFep')
   lick_hb_arr=line_index(wl, SPEC[*,idx_model], err_arr, indexpars[where(strcmp(indexpars.index,'Lick_Hb'))])
   
   d4000n=d4000n_arr[0]
   err_d4000n=d4000n_arr[1]
   hdhg=HdHg_arr[0]
   err_hdhg=HdHg_arr[1]
   mg2fe=mg2fe_arr[0]
   err_mg2fe=mg2fe_arr[1]
   mgfep=mgfep_arr[0]
   err_mgfep=mgfep_arr[1]
   lick_hb=lick_hb_arr[0]
   err_lick_hb=lick_hb_arr[1]
   
   data_table[idx_model].(col_d4000)=d4000n
   data_table[idx_model].(col_d4000_err)=err_d4000n
   data_table[idx_model].(col_HdHg)=hdhg
   data_table[idx_model].(col_hdhg_err)=err_hdhg
   data_table[idx_model].(col_mg2fe)=mg2fe
   data_table[idx_model].(col_mg2fe_err)=err_mg2fe
   data_table[idx_model].(col_mgfep)=mgfep
   data_table[idx_model].(col_mgfep_err)=err_mgfep
   data_table[idx_model].(col_lick_hb)=lick_hb
   data_table[idx_model].(col_lick_hb_err)=err_lick_hb
  ;stop 
  endfor
 stop
 endfor 
 ;stop
 
end 