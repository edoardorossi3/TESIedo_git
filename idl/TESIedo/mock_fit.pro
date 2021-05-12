@BaStA_fit.pro
pro mock_fit, i_chunk
  common dirnames
  common filenames
  common dust_cond
  common models_common
  common data_common
  common idxrange
  common filters
  common pdf_common
  
  i_chunk=1
  
  ;stop
  ;pdf_init,getenv('BASTA_DIR')+'/idl/BaStA/pdf_parfiles/pdf_modpars_stelpop_SanV4.1eqMILES_ChFall_test_agef.par'
  pdf_init, '/home/edoardo/BaStA/idl/BaStA/pdf_parfiles/pdf_modpars_stelpop_SanV4.2eqMILES_ChFall_test_agef.par'

  
  redshift=0.0 ;lavora 0
  vdisp=200. & sigma_instr=sigma_mod ;lavora a 200
  ;n_chunks=1 ;40
  
  mag_mod_redshift_interpol,redshift
  indx_mod_vdisp_interpol,vdisp
  
 ; for i_chunk=1, n_chunks do begin
    ;; read mock table
    ;mock_fname='/Users/zibetti/ownCloud/Tesi_ERossi/models/sandage_varZ_v4.1eq_spec_dcomb_perterr_'+string(i_chunk, format='(I03)')+'.fits'
    ;mock_fname=getenv('SEDLIBRARIES_DIR')+'/Sandage_varZ_v4.1eq_bc03MILES_ChFall/mock_ER_001/sandage_varZ_v4.1eq_spec_dcomb_perterr_H_'+string(i_chunk, format='(I03)')+'.fits' ;considero anche la banda H (_perterr_H_)
    mock_fname=getenv('SEDLIBRARIES_DIR')+'/Sandage_varZ_v4.2eq_CB16MILES_ChFall/sandage_varZ_v4.2eq_spec_dcomb0p25null_perterr_'+string(i_chunk, format='(I03)')+'.fits'
    models_dir=getenv('SEDLIBRARIES_DIR')+'/Sandage_varZ_v4.2eq_CB16MILES_ChFall/mock_ER_001/'
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
      filename_ph=models_dir+'mock_file_'+phpars_to_fit[i_par]+string(i_chunk, format='(I03)')+'.fits'
      ;filename_ph=models_dir+'mock_file_'+phpars_to_fit[i_par]+string(i_chunk, format='(I03)')+'.fits' ;se considero anche la banda H (_H)
      mwrfits, tables.(i_par), filename_ph, /create
    endfor



    ;stop
    
    
  ;endfor
  
  end
  
  pro parall_mock_fit
  maxthreads=10 ; max number of parallel threads
  obridge=objarr(maxthreads)
  cur_chunk=intarr(maxthreads)
  nchunks=40
  nthreads=0 ; number of active threads
  i_chunk=0
  SNR=20.
  watchinterval=2
  while i_chunk lt nchunks or nthreads gt 0 do begin
    if nthreads lt maxthreads and i_chunk lt nchunks then begin ; still chunks to process and thread free
      ;; look for the first free thread
      for ithread=0,maxthreads-1 do begin
        if not(OBJ_VALID(obridge[ithread])) and i_chunk lt nchunks then begin
          ; launch a new process
          ; get ready for the first chunk to process
          obridge[ithread]=obj_new('IDL_IDLBridge',output=getenv('HOME')+'/idl/log_files/mock_'+string(i_chunk+1, format='(I02)')+'.log') ; write terminal output to the named file
          obridge[ithread]->execute,'@'+pref_get('IDL_STARTUP') ; execute IDL_STARTUP commands
          obridge[ithread]->execute,'.r '+getenv('BASTA_DIR')+'/idl/BaStA/BaStA_index.pro'
          obridge[ithread]->execute,'.r '+getenv('TESI_EDO')+'mock_fit.pro' ;; fill in your idl .pro file to be compiled
          obridge[ithread]->execute,'mock_fit,'+string(i_chunk+1),/nowait  ;; fill in the command line to process a single chunk i_chunk
          cur_chunk[ithread]=i_chunk
          print,'Job '+string(ithread)+' started (chunk='+string(cur_chunk[ithread])+')'
          nthreads=nthreads+1
          i_chunk++
        endif
      endfor
    endif else begin
      ;; all threads are already busy
      wait,watchinterval/2. ;; wait for a while
      for ithread=0,maxthreads-1 do begin
        if OBJ_VALID(obridge[ithread]) then begin
          ;; check if any of the threads has finished
          objstatus=obridge[ithread]->status()
          if objstatus gt 1 then begin
            print,'Job '+string(ithread)+' finished (chunk='+string(cur_chunk[ithread])+')'
            obridge[ithread]->cleanup
            nthreads--
          endif
        endif
      endfor
      wait,watchinterval/2. ;; wait for a while for process to actually terminate
    endelse
  endwhile

end

  
