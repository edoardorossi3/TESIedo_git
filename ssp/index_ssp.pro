pro index_ssp

  common idxinfo
  read_idxinfo
  SSPdir='/home/edoardo/Desktop/TESI/models/ssp/'
  SSP_prefix='bc2003_hr_xmiless_m'
  SSP_suffix='_chab_ssp_ised.fits'
  SSP_idxsuffix='_chab_ssp_idx.fits'
  ssp_num=['22', '32', '42', '52', '62', '72']
  compidx_names=['HdHg','Mg2Fe','MgFep']
  
  lineidx_names=indexpars.index
  data_row=create_struct('age', !values.f_nan)
  for i=0, n_elements(lineidx_names)-1 do begin
    data_row=create_struct(data_row,lineidx_names[i], !values.f_nan)
  endfor
  data_row=create_struct(data_row, 'D4000n', !values.f_nan, 'D4000', !values.f_nan, 'HdHg', !values.f_nan, 'Mg2Fe', !values.f_nan, 'MgFep', !values.f_nan)
  
  colnames=tag_names(data_row)
  ;; read SSP
  ;SSPfname='/Users/SteMaGE/tools/SEDlibrary/SSP/bc2003_hr_xmiless_m62_chab_ssp_ised.fits'


for i_ssp=0, n_elements(ssp_num)-1 do begin
  SSPfname=SSPdir+SSP_prefix+ssp_num[i_ssp]+SSP_suffix
  wl=mrdfits(SSPfname,0)
  uvoptnir=where(wl gt 3000. and wl lt 25000.,n_uvoptnir)
  SEDgrid=mrdfits(SSPfname,1)
  phys_table=mrdfits(SSPfname, 2)
  ;stop
  n_age=(size(SEDgrid))[2]
  data_table=replicate(data_row, n_age)
  data_table.age=phys_table.age
  for age_idx=0, n_age-1 do begin
    for j_idx=0, n_elements(lineidx_names)-1 do begin
      _tmp=line_index(wl[uvoptnir], SEDgrid[uvoptnir, age_idx], [], indexpars[j_idx])
      col_num=(where(strcmp(colnames,lineidx_names[j_idx],/fold_case)))[0]
      data_table[age_idx].(col_num)= _tmp[0]
    endfor
    _tmp=d4000_index(wl[uvoptnir], SEDgrid[uvoptnir, age_idx], [], /narrow)
    col_num=(where(strcmp(colnames,'D4000n',/fold_case)))[0]
    data_table[age_idx].(col_num)= _tmp[0]
    _tmp=d4000_index(wl[uvoptnir], SEDgrid[uvoptnir, age_idx], [])
    col_num=(where(strcmp(colnames,'D4000',/fold_case)))[0]
    data_table[age_idx].(col_num)= _tmp[0]
    for j_idx=0, n_elements(compidx_names)-1 do begin
      _tmp=comp_index(wl[uvoptnir], SEDgrid[uvoptnir, age_idx], [], compidx_names[j_idx])
      col_num=(where(strcmp(colnames,compidx_names[j_idx],/fold_case)))[0]
      data_table[age_idx].(col_num)= _tmp[0]
    endfor
    
    ;stop
  endfor
mwrfits, data_table, SSPdir+SSP_prefix+ssp_num[i_ssp]+SSP_idxsuffix, /create
;  stop

endfor
  
;  D4000=fltarr(221)
;  HgHd=fltarr(221)
;  err=[]
;
;  for age_idx=0, 220 do begin
;    ;age_idx=10
;    SNR=20.
;    narrow=1
;    D4000[age_idx]=(D4000_index(wl[uvoptnir],SEDgrid[uvoptnir,age_idx],err,narrow=narrow))(0)
;    HgHd[age_idx]=(comp_index(wl[uvoptnir],SEDgrid[uvoptnir,age_idx],err,'HgHd'))(0)
;    ;stop
;  endfor
;  ;stop
;  ;print, D4000
;  age_ssp=mrdfits(SSPfname, 2, columns=1)
;  ;print, age_ssp
;  data_row=create_struct('age_ssp', 0., 'D4000n', 0.)
;  ;mwrfits, D4000, 'ssp_index.fits', 'D4000n'
;
;  stop
end