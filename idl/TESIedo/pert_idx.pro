pro pert_idx

;err_file='sandage_varZ_v4.1eq_spec_dcomb_idxerr_001.fits'
;idx_file='sandage_varZ_v4.1eq_spec_dcomb_idx_001.fits'
;spec_file='sandage_varZ_v4.1eq_spec_dcomb_001.fits'
;perterr_file='sandage_varZ_v4.1eq_spec_dcomb_perterr_001.fits'

prefix_file='sandage_varZ_v4.1eq_spec_dcomb_'
model_dir=getenv('SEDLIBRARIES_DIR')+'/Sandage_varZ_v4.1eq_bc03MILES_ChFall/'

err_mag=0.05
Flux_calib_err=0.03 ;sistematic error on d4000 due to flux calibration
SNR=20.
seed=-10
n_chunks=1  ;numbers of fits file to read

for i_chunk=1, n_chunks do begin
  ;we define strings to read our fits file
  end_err_file='idxerr_'+string(i_chunk, format='(I03)')+'.fits'
  end_idx_file='idx_'+string(i_chunk, format='(I03)')+'.fits'
  end_spec_file=string(i_chunk, format='(I03)')+'.fits'
  end_perterr_file='perterr_'+string(i_chunk, format='(I03)')+'.fits'
  end_physpar_file=string(i_chunk, format='(I03)')+'_physpar'+'.fits'
  
  err_file=model_dir+prefix_file+end_err_file
  idx_file=model_dir+prefix_file+end_idx_file
  spec_file=model_dir+prefix_file+end_spec_file
  perterr_file=model_dir+prefix_file+end_perterr_file
  physpar_file=model_dir+prefix_file+end_physpar_file
  
  ;stop
  ;we read the info from fits files
  wl=mrdfits(spec_file, 0)
  SPEC=mrdfits(spec_file, 1) & n_models=(size(SPEC))[2]
  idx_table=mrdfits(idx_file,1)
  idxwerr_table=mrdfits(err_file, 1)
  phys_table=mrdfits(physpar_file, 1)
  ;idx_names=(tag_names(idx_table))[1:-1] & n_idx=n_elements(idx_names)
  idx_names=['D4000n', 'HDHG', 'LICK_Hb', 'MGFE_PRIME', 'MG2FE'] & n_idx=n_elements(idx_names)
  idxwerr_names=(tag_names(idxwerr_table))[1:-1]
  
  ;We create the struct which will contain idx pert. and the err.
  data_row=create_struct('IDX', 0L)
  for i_idx=0, n_idx-1 do begin
    data_row=create_struct(data_row, idx_names[i_idx]+'_SIGMA200_PERT', !values.f_nan, idx_names[i_idx]+'_ERR_SNR'+string(fix(SNR),format='(I03)'), !values.f_nan)
  endfor
  ;Add ugriz ab mag
  data_row=create_struct(data_row, 'ABMAG_u_pert', !values.f_nan, 'ABMAG_g_pert', !values.f_nan, 'ABMAG_r_pert', !values.f_nan, 'ABMAG_i_pert', !values.f_nan, 'ABMAG_z_pert', !values.f_nan, 'ERR_MAG', err_mag)
  
  
  data_table=replicate(data_row, n_models)
  data_table.idx=idx_table.idx
  
  RND_mag=randomn(seed, n_models)*err_mag

  data_table.abmag_u_pert=phys_table.abmag[0]+RND_mag
  data_table.abmag_g_pert=phys_table.abmag[1]+RND_mag
  data_table.abmag_r_pert=phys_table.abmag[2]+RND_mag
  data_table.abmag_i_pert=phys_table.abmag[3]+RND_mag
  data_table.abmag_z_pert=phys_table.abmag[4]+RND_mag

  ;stop
  data_names=(tag_names(data_table))
  ;Now, we want to calculate the pert. idexes
  
  for i_idx=0, n_idx-1 do begin
    
    n_col200=(where(strcmp(idxwerr_names ,idx_names[i_idx]+'_SIGMA200',/fold_case)))[0]+1
    n_colerr=(where(strcmp(idxwerr_names ,idx_names[i_idx]+'_ERR_SNR'+string(fix(SNR),format='(I03)'),/fold_case)))[0]+1
    n_colpert=(where(strcmp(data_names ,idx_names[i_idx]+'_SIGMA200_PERT',/fold_case)))[0]
    n_colerrout=n_colpert+1
    _tmp200=idxwerr_table[*].(n_col200)
    _tmperr=idxwerr_table[*].(n_colerr)
    if (strcmp(idx_names[i_idx], 'D4000N', /fold_case)) then _tmperr=sqrt(_tmperr^2+Flux_calib_err^2)
    RND=randomn(seed, n_models)*_tmperr
    data_table[*].(n_colpert)=_tmp200+RND
    data_table[*].(n_colerrout)=_tmperr
    
    ;stop
  endfor
  
  ;stop
  
mwrfits,  data_table, perterr_file, /create  
  
;stop
endfor





end