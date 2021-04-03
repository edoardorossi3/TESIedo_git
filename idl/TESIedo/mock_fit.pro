@BaStA_fit.pro
pro mock_fit
  common dirnames
  common filenames
  common dust_cond
  common models_common
  common data_common
  common idxrange
  common filters
  common pdf_common
  
  ;stop
  ;pdf_init,getenv('BASTA_DIR')+'/idl/BaStA/pdf_parfiles/pdf_modpars_stelpop_SanV4.1eqMILES_ChFall_test_agef.par'
  pdf_init, '/home/edoardo/BaStA/idl/BaStA/pdf_parfiles/pdf_modpars_stelpop_SanV4.1eqMILES_ChFall_test_agef.par'

  
  redshift=0.0 ;lavora 0
  vdisp=200. & sigma_instr=sigma_mod ;lavora a 200
  n_chunks=1 ;40
  
  mag_mod_redshift_interpol,redshift
  indx_mod_vdisp_interpol,vdisp
  
  for i_chunk=1, n_chunks do begin
    ;; read mock table
    ;mock_fname='/Users/zibetti/ownCloud/Tesi_ERossi/models/sandage_varZ_v4.1eq_spec_dcomb_perterr_'+string(i_chunk, format='(I03)')+'.fits'
    ;mock_fname=getenv('SEDLIBRARIES_DIR')+'/Sandage_varZ_v4.1eq_bc03MILES_ChFall/mock_ER_001/sandage_varZ_v4.1eq_spec_dcomb_perterr_H_'+string(i_chunk, format='(I03)')+'.fits' ;considero anche la banda H (_perterr_H_)
    mock_fname=getenv('SEDLIBRARIES_DIR')+'/Sandage_varZ_v4.1eq_bc03MILES_ChFall/mock_ER_001/sandage_varZ_v4.1eq_spec_dcomb_perterr_H003_'+string(i_chunk, format='(I03)')+'.fits'
    models_dir=getenv('SEDLIBRARIES_DIR')+'/Sandage_varZ_v4.1eq_bc03MILES_ChFall/mock_ER_001/'
    mock_table=mrdfits(mock_fname,1)
    ;stop
    ;idx_names=['D4000n',  'LICK_Hb', 'HDHG', 'MGFE_PRIME','MG2FE'] & n_idx=n_elements(idx_names)
    ;mag_names=[ 'u', 'g', 'r', 'i', 'z'] & n_mag=n_elements(mag_names)
    n_idx=n_elements(index_to_fit)
    n_mag=n_elements(phbands_to_fit)


    mock_names=tag_names(mock_table)
    n_mock=n_elements(mock_table[*])
    n_par=n_elements(phpars_to_fit)


    tables={}
    row_struct=pdf_phpars.(0)
    for i_par=0, n_par-1 do begin

      tables=create_struct(tables, phpars_to_fit[i_par], replicate(row_struct, n_mock))
    endfor

    ;stop
    for idx_mock=0, n_mock-1 do begin
      for i_idx=0, n_idx-1 do begin
        col_idx=(where(strcmp(index_to_fit[i_idx]+'_sigma200_pert', mock_names, /fold_case)))[0]
        col_err=col_idx+1
        _idx=mock_table[idx_mock].(col_idx)
        _err=mock_table[idx_mock].(col_err)
        index_data_val.(i_idx)=_idx
        index_data_err.(i_idx)=_err
      endfor

      for i_mag=0, n_mag-1 do begin
        col_mag=(where(strcmp('abmag_'+phbands_to_fit[i_mag]+'_pert', mock_names, /fold_case)))[0]
        col_err=col_mag+1
        _mag=mock_table[idx_mock].(col_mag)
        _err=mock_table[idx_mock].(col_err)
        mag_data_val.(i_mag)=_mag
        mag_data_err.(i_mag)=_err
      endfor
      ;stop

      results =pdf_simple(mod_exclude=[idx_mock])
      ;stop
      for i_par=0, n_par-1 do begin
        tables.(i_par)[idx_mock]=results.(i_par)
      endfor


    endfor
    ;stop
    for i_par=0, n_par-1 do begin
      ;filename_ph=models_dir+'mock_file_'+phpars_to_fit[i_par]+'_H'+string(i_chunk, format='(I03)')+'.fits' ;se considero anche la banda H (_H)
      filename_ph=models_dir+'mock_file_'+phpars_to_fit[i_par]+'_H003'+string(i_chunk, format='(I03)')+'.fits'
      ;filename_ph=models_dir+'mock_file_'+phpars_to_fit[i_par]+string(i_chunk, format='(I03)')+'.fits' ;se considero anche la banda H (_H)
      mwrfits, tables.(i_par), filename_ph, /create
    endfor



    stop
    
    
  endfor
  
  end
  
  
  ;; read mock table
  ;mock_fname='/Users/zibetti/ownCloud/Tesi_ERossi/models/sandage_varZ_v4.1eq_spec_dcomb_perterr_001.fits'
  ;mock_fname=getenv('SEDLIBRARIES_DIR')+'/Sandage_varZ_v4.1eq_bc03MILES_ChFall/mock_ER_001/sandage_varZ_v4.1eq_spec_dcomb_perterr_H_001.fits' ;considero anche la banda H (_perterr_H_)
  ;models_dir=getenv('SEDLIBRARIES_DIR')+'/Sandage_varZ_v4.1eq_bc03MILES_ChFall/mock_ER_001/'
  ;mock_table=mrdfits(mock_fname,1)
  ;stop
  ;idx_names=['D4000n',  'LICK_Hb', 'HDHG', 'MGFE_PRIME','MG2FE'] & n_idx=n_elements(idx_names)
  ;mag_names=[ 'u', 'g', 'r', 'i', 'z'] & n_mag=n_elements(mag_names)
  ;n_idx=n_elements(index_to_fit)
  ;n_mag=n_elements(phbands_to_fit)
  
 
  ;mock_names=tag_names(mock_table)
  ;n_mock=n_elements(mock_table[*])
  ;n_par=n_elements(phpars_to_fit)
  
 
  ;tables={}
  ;row_struct=pdf_phpars.(0)
  ;for i_par=0, n_par-1 do begin
    
   ; tables=create_struct(tables, phpars_to_fit[i_par], replicate(row_struct, n_mock))
  ;endfor
  
  ;stop
  ;for idx_mock=0, n_mock-1 do begin
   ; for i_idx=0, n_idx-1 do begin
    ;  col_idx=(where(strcmp(index_to_fit[i_idx]+'_sigma200_pert', mock_names, /fold_case)))[0]
    ;  col_err=col_idx+1
    ;  _idx=mock_table[idx_mock].(col_idx)
    ;  _err=mock_table[idx_mock].(col_err)
     ; index_data_val.(i_idx)=_idx
     ; index_data_err.(i_idx)=_err
    ;endfor
   
  ;  for i_mag=0, n_mag-1 do begin
  ;    col_mag=(where(strcmp('abmag_'+phbands_to_fit[i_mag]+'_pert', mock_names, /fold_case)))[0]
   ;   col_err=col_mag+1
   ;   _mag=mock_table[idx_mock].(col_mag)
    ;  _err=mock_table[idx_mock].(col_err)
    ;  mag_data_val.(i_mag)=_mag
    ;  mag_data_err.(i_mag)=_err
    ;endfor
  ;stop
    
    ;results =pdf_simple(mod_exclude=[idx_mock]) 
    ;stop
    ;for i_par=0, n_par-1 do begin
     ; tables.(i_par)[idx_mock]=results.(i_par)
    ;endfor
    
    
  ;endfor  
  ;stop
  ;for i_par=0, n_par-1 do begin
   ; filename_ph=models_dir+'mock_file_'+phpars_to_fit[i_par]+'_H.fits' ;se considero anche la banda H (_H)
    ;filename_ph=models_dir+'mock_file_'+phpars_to_fit[i_par]+'.fits' ;se considero anche la banda H (_H)
    ;mwrfits, tables.(i_par), filename_ph, /create
  ;endfor
  
  
  
  ;stop

