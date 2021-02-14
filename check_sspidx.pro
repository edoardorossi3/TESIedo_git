pro check_sspidx
  SSPdir='/home/edoardo/Desktop/TESI/models/ssp/'
  milesdir='bc03_miles/'
  SSP_prefix='bc2003_hr_xmiless_m'
  SSP_idxsuffix='_chab_ssp_idx.fits'
  ssp_num=['22', '32', '42', '52', '62', '72']
  
  SSPfname=SSPdir+SSP_prefix+ssp_num[4]+SSP_idxsuffix
  phys_table=mrdfits(SSPfname, 1)
  D4000n=phys_table.d4000n
  HgHd=phys_table.HdHg
  B4_VN=fltarr(221)
  HgHdssp=fltarr(221)
  table=fltarr(220, 23)
  readcol, SSPdir+milesdir+SSP_prefix+ssp_num[4]+'_chab_ssp.7lsindx_sed', v1,v2,v3,v4,v5,v6,v7,v8,v9,hdssp,hgssp,v12,v13,v14,d4000ssp, SKIPLINE=32
  ;stop
  D4000n[0]=0.0
  HgHd[0]=0.0
  
    B4_VN[1:-1]=d4000ssp
    HgHdssp[1:-1]=hgssp+hdssp
  ugualid=where(D4000n eq B4_VN,COUNTd)
  ugualih=where(HgHd eq HgHdssp, COUNTh)
 
  
  if(COUNTd eq n_elements(d4000ssp)) then print, 'Good match for D4000n indexes!' else print, 'Attention! There are differences in D4000n indexes!'
  if(COUNTh eq n_elements(hgssp)) then print, 'Good match for HgHd indexes!' else print, 'Attention! There are differences in HgHd indexes!'
;stop
  
  maxerr_d4000=max(abs(B4_VN-D4000n))
  maxerr_hghd=max(abs(HgHdssp-HgHd)) 
  print, 'max difference d4000:', maxerr_d4000
  print, 'max difference hghd:', maxerr_hghd
  
;stop
end