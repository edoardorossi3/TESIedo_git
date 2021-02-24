pro index_errors,SNR
  common idxinfo
  read_idxinfo

  if n_elements(SNR) eq 0  then SNR=20. ; if SNR not given, assume default SNR=20.

  ;prog_dir='~/Desktop/TESIedo_git/'
  models_dir=getenv('SEDLIBRARIES_DIR')+'/Sandage_varZ_v4.1eq_bc03MILES_ChFall/' ;'~/Desktop/TESI/models/sandage_dcomb/'
  file_prefix='sandage_varZ_v4.1eq_spec_'
  spec_suffix='dcomb_'
  idx_suffix='dcomb_idx_'
  idx_suffix_werr='dcomb_idxerr_'
  file_ext='.fits'
  ;; sndg_num=indgen(2) ; vuoi un array intero che parta da 1, non lo generi facendone uno in piu' e partendo poi dal secondo elemento
  ;tot_file=n_elements(sndg_num)
  ; in ogni caso non serve lavorare con un array se tanto fai un loop, ti basta sapere il numero di chunks
  n_chuncks=1 ;40

  optical_range=[3500.,9000.];[2500.,25000.] ;; wl range for the optical spectrum
  SNR_range=[4000.,5350.] ;; wl range to define <SNR>, roughly corresponding to SDSS-g passband
  dwl_rsmp=1.0
  ; define the resampled wl array only once
  n_wlrsmp=fix((optical_range[1]-optical_range[0])/dwl_rsmp)+1
  wl_rsmp=findgen(n_wlrsmp)*dwl_rsmp+optical_range[0]

  data_row=create_struct('IDX', 0L) ; !values.f_nan is a float

  for i_chunk=1, n_chuncks do begin  ; dai al contatore un nome che ti ricordi su che cosa stai facendo il loop.
    ;Se chiami in modo anonimo i j k l m... i contatori in una serie di nested loops, non ci si capisce piu' nulla

    spec_file=models_dir+file_prefix+spec_suffix+string(i_chunk, format='(I03)')+file_ext  ;; preferibile I03 piuttosto che I3.3  (se e' intero non ha decimali)
    idx_file=models_dir+file_prefix+idx_suffix+string(i_chunk, format='(I03)')+file_ext
    idx_werr_file=models_dir+file_prefix+idx_suffix_werr+string(i_chunk, format='(I03)')+file_ext

    wl=mrdfits(spec_file, 0)
    SPEC=mrdfits(spec_file, 1) & n_models=(size(SPEC))[2] ; keep information close in space...
    idx_table=mrdfits(idx_file,1)
    ;; check for matching sizes
    if n_elements(idx_table) ne n_models then message,'Index table and SPECTRA array have different sizes'

    ; define optical index array to limit the size of the spectral array to work with
    optical=where(wl ge optical_range[0] and wl le optical_range[1], n_optical)

    ; structure data_row must be created only once - should be done either outside the loop or only the first time or if it's not defined yet
    if n_tags(data_row) eq 1 then begin ;; if data row has only one tag it means it only contains the IDX tag and must be initialized
      linenames=indexpars.index
      idx_names=(tag_names(idx_table))[1:-1] ; skip first tag: IDX = not an index name
      Nidx_names=n_elements(idx_names)
      for i_index=0, Nidx_names-1 do begin ; start from 0, not 1!
        data_row=create_struct(data_row, idx_names[i_index]+'_SIGMA0', !values.f_nan, idx_names[i_index]+'_SIGMA200', !values.f_nan, $
          idx_names[i_index]+'_SIGMA0_RSMP', !values.f_nan, idx_names[i_index]+'_ERR_SNR'+string(fix(SNR),format='(I03)'), !values.f_nan)
      endfor
      colnames=tag_names(data_row)
    endif

    ; create data_table structure
    data_table=replicate(data_row, n_models)
    data_table.idx=idx_table.idx

    ;; this loop does not need to be inside the model loop - just assign arrays, not element by element
    for i_index=0, Nidx_names-1 do begin
      col_num0=(where(strcmp(colnames,idx_names[i_index]+'_SIGMA0',/fold_case)))[0]
      col_num200=(where(strcmp(colnames,idx_names[i_index]+'_SIGMA200',/fold_case)))[0]
      _tmp=idx_table.(i_index+1) ; use _ to make it clear is a temporary variable - i_index+1 because 1st column is IDX
      data_table.(col_num0)=reform(_tmp[0,*]) ; reform is used to reduce dimensionality when size is 1 - eg here an array [1,12500] becomes [12500]
      data_table.(col_num200)=reform(_tmp[4,*])
    endfor

    ;stop


    for idx_model=0, n_models-1 do begin
      ;;  common data
      ; idx_lim=where(wl gt 3700. and wl lt 8860., nlim_wl)  ;; idx_lim stands for??
      ; wl_maxmin=[wl[idx_lim[0]], wl[idx_lim[-1]]]
      ;flx=flux_integral(wl, SPEC[*,idx_model], wl_maxmin) ;; why this complication to redefine wl ranges?
      flx=flux_integral(wl[optical], SPEC[optical,idx_model], SNR_range)
      mean_flx=flx/(SNR_range[1]-SNR_range[0])
      mean_err=mean_flx/SNR
      SPEC_rsmp=interpol(SPEC[optical,idx_model], wl[optical], wl_rsmp, /spline)
      err_spec=replicate(mean_err,n_wlrsmp)

      for i_index=0, Nidx_names-1 do begin
        ;; the following two lines are independent of the type of index, take out of the if-then clause
        col_num_rsmp=(where(strcmp(idx_names[i_index]+'_SIGMA0_RSMP', colnames, /fold_case)))
        col_num_err=(where(strcmp(idx_names[i_index]+'_ERR_SNR'+string(fix(SNR),format='(I03)'), colnames, /fold_case)))
        _index=[!values.f_nan,!values.f_nan]
        ;if (where(strcmp(linenames, idx_names[k], /fold_case)) ne -1) then begin ;;; where returns an array and you test a condition by comparing with a scalar
        ; look for a matching simple line index in linenames
        if total(strcmp(linenames, idx_names[i_index], /fold_case)) eq 1 then begin
          _index=line_index(wl_rsmp, SPEC_rsmp, err_spec, idx_names[i_index])
        endif else begin ; use nested else to avoid wasting time in checks
          if strcmp(idx_names[i_index], 'd4000n', /fold_case) then begin ; use where only if you need to find the locations for TRUE condition within an array
            _index=d4000_index(wl_rsmp, SPEC_rsmp, err_spec, /narrow)
          endif else  begin
            if strcmp(idx_names[i_index], 'd4000', /fold_case) then begin
              _index=d4000_index(wl_rsmp, SPEC_rsmp, err_spec, narrow=0)
            endif else begin
              if total(strcmp(cmp_index_list, idx_names[i_index], /fold_case)) eq 1 then begin
                _index=comp_index(wl_rsmp, SPEC_rsmp, err_spec, idx_names[i_index])
              endif
            endelse
          endelse
        endelse
        ;; the following two lines are independent of the type of index, take out of the if-then clause
        data_table[idx_model].(col_num_rsmp)=_index[0]
        data_table[idx_model].(col_num_err)=_index[1]
      endfor
    endfor
    ;stop
    
    ; TODO: take header from idx_table file and modify to update column info and add extra info (eg, spline interpolation, SNR etc)
    mwrfits, data_table, idx_werr_file, /create
    hdr=HEADFITS(idx_werr_file, exten=1)
    sxaddpar, hdr, 'COMMENT', '_ERR_SNR'+string(fix(SNR),format='(I03)')+': SNR='+string(fix(SNR),format='(F5.1)')+' per '+string(dwl_rsmp,format='(F3.1)')+' AA',after='TFIELDS'
    sxaddpar, hdr, 'COMMENT', '_SIGMA0_RSMP: 0 velocity dispersion, resmapled at '+string(dwl_rsmp,format='(F3.1)')+' AA',after='TFIELDS'
    sxaddpar, hdr, 'COMMENT', '_SIGMA200: 200 km/s velocity dispersion from SEDLIBRARY',after='TFIELDS'
    sxaddpar, hdr, 'COMMENT', '_SIGMA0: 0 velocity dispersion from SEDLIBRARY',after='TFIELDS'
    sxaddpar, hdr, 'RSMP_TYPE', 'SPLINE', 'Type of resampling interpolant' ,after='TFIELDS'

    sxaddpar, hdr, 'SNR', SNR, 'SNR per '+string(dwl_rsmp,format='(F3.1)')+' AA for error computation' ,after='TFIELDS'
    sxaddpar, hdr, 'DLAMBDA_RSMP', dwl_rsmp, 'Resampled pixel size in AA for error computation',after='TFIELDS'
    sxaddpar, hdr, 'SPEC_ORIG', file_prefix+spec_suffix+string(i_chunk, format='(I03)')+file_ext, 'Orig. SEDLIBRARY spec',after='TFIELDS'
    sxaddpar, hdr, 'IDX_ORIG', file_prefix+idx_suffix+string(i_chunk, format='(I03)')+file_ext, 'Orig. SEDLIBRARY indxs',after='TFIELDS'
    modfits, idx_werr_file, 0, hdr, exten=1
    
  ;stop
  endfor
  
  
  stop

end

pro check_interp

  common idxinfo
  read_idxinfo

  models_dir='~/Desktop/TESI/models/sandage_dcomb/'
  file_prefix='sandage_varZ_v4.1eq_spec_'
  spec_suffix='dcomb090n_'
  idx_suffix='dcomb090n_idx_'
  idx_suffix_werr='dcomb090n_idx_werr'
  file_ext='.fits'
  sndg_num=1
  tot_file=n_elements(sndg_num)

  SNR=20.0
  dwl_reb=1.0
  idx_model=0
  spec_file=models_dir+file_prefix+spec_suffix+string(sndg_num, format='(I3.3)')+file_ext
  idx_file=models_dir+file_prefix+idx_suffix+string(sndg_num, format='(I3.3)')+file_ext
  idx_werr_file=models_dir+file_prefix+idx_suffix_werr+string(sndg_num, format='(I3.3)')+file_ext

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
  n_wlreb=fix((wl_maxmin[1]-wl_maxmin[0])/dwl_reb)+1
  wl_rebin=findgen(n_wlreb)*dwl_reb+wl_maxmin[0]

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
