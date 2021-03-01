@BaStA_fit.pro
pro test_BaStA_fit
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

  
  redshift=0.005
  vdisp=175. & sigma_instr=sigma_mod
  mag_mod_redshift_interpol,redshift
  indx_mod_vdisp_interpol,vdisp
  
  ;; read mock table
  ;mock_fname='/Users/zibetti/ownCloud/Tesi_ERossi/models/sandage_varZ_v4.1eq_spec_dcomb_perterr_001.fits'
  mock_fname='/home/edoardo/Desktop/TESI/models/Sandage_varZ_v4.1eq_bc03MILES_ChFall/sandage_varZ_v4.1eq_spec_dcomb_perterr_001.fits'

  mock_table=mrdfits(mock_fname,1)
  ;stop
  idx_names=['D4000n',  'LICK_Hb', 'HDHG', 'MGFE_PRIME','MG2FE'] & n_idx=n_elements(idx_names)
  mag_names=[ 'u', 'g', 'r', 'i', 'z'] & n_mag=n_elements(mag_names)
  idx_mock=10
  mock_names=tag_names(mock_table)
  ;stop
  n_models=n_elements(mock_table[*])
 
  
    for i_idx=0, n_idx-1 do begin
      col_idx=(where(strcmp(idx_names[i_idx]+'_sigma200_pert', mock_names, /fold_case)))[0]
      col_err=col_idx+1
      _idx=mock_table[idx_mock].(col_idx)
      _err=mock_table[idx_mock].(col_err)
      index_data_val.(i_idx)=_idx
      index_data_err.(i_idx)=_err
    endfor
   
    for i_mag=0, n_mag-1 do begin
      col_mag=(where(strcmp('abmag_'+mag_names[i_mag]+'_pert', mock_names, /fold_case)))[0]
      col_err=col_mag+1
      _mag=mock_table[idx_mock].(col_mag)
      _err=mock_table[idx_mock].(col_err)
      mag_data_val.(i_mag)=_mag
      mag_data_err.(i_mag)=_err
    endfor
  ;stop
    
  results=pdf_simple(mod_exclude=[idx_mock])
    
  stop
end
