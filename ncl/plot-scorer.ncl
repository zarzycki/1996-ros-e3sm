begin

minlat=28.5
maxlat=47.0
minlon=-102.0
maxlon=-65.0

pnew = (/ 850., 500. /)       ; desired output levels 
  
plot_CAM=True

era5dir="/Volumes/ZARZYCKI_FLASH/arp5873/ERA5/"
CAMdir="/Volumes/ZARZYCKI_FLASH/arp5873/ens_means/"
  

if (plot_CAM) then
  f = addfile(CAMdir+"/output_control_eam_h2.nc","r")
  timeix=16
  stride=1
  Z3 = f->Z3(timeix,:,::stride,::stride)
  PS = f->PS(timeix,::stride,::stride)
  T = f->T(timeix,:,::stride,::stride)
  Q = f->Q(timeix,:,::stride,::stride)
  U = f->U(timeix,:,::stride,::stride)
  V = f->V(timeix,:,::stride,::stride)
  ;W = mixhum_convert(Q, "q", (/0,0/))
  ;copy_VarCoords(Q,W)

  time = f->time(timeix)
  print(cd_calendar(time,-3))

  hyam = f->hyam                   ; read to memory (optional)
  hybm = f->hybm
  P0mb = 1000.      
  P0 = f->P0

  nlev = dimsizes(hyam)

  PMID  = new ( dimsizes(T), typeof(T) )
  do n=0,nlev-1
    PMID(n,:,:)  = tofloat(P0*hyam(n) + hybm(n)*PS)
  end do

  THETA = pot_temp(PMID, T, -1, False)  ; dim=0

  copy_VarCoords(T,THETA)

  All  = rigrad_bruntv_atm(THETA, U, V, Z3, 3, 0 ) 
  ;RI   = All[0]     ; explicitly extract from list variable  for convenience
  BV   = All[1]     ; must use 'list' syntax  [...]
  ;BUOY = All[2]
  SHR2 = All[3]
  
  WIND = U
  WIND = sqrt(U*U + V*V)
  
  DUDZ = center_finite_diff_n (WIND,Z3,False,0,0)
  DDUDZ  = center_finite_diff_n (DUDZ,Z3,False,0,0)
  DTHETADZ = center_finite_diff_n (THETA,Z3,False,0,0)
  
  ;BV2 = brunt_vaisala_atm(THETA, Z3, 0, 0)  ; opt=0, ndim=1
  
  ;TERM1 = WIND
  ;TERM1 = 9.81 * DTHETADZ / WIND / WIND / THETA 
  
  TERM1 = WIND
  TERM1 = BV * BV / WIND / WIND
  
  TERM2 = WIND
  TERM2 = DDUDZ / WIND

  L2 = WIND
  L2 = TERM1 - TERM2
  
  ;L2 = BV
  ;L2 = BV / (WIND*WIND) - DDUDZ / WIND

else


  zfile = addfile(era5dir+"/e5.oper.an.pl.128_129_z.ll025sc.1996011900_1996011923.nc","r")
  tfile = addfile(era5dir+"/e5.oper.an.pl.128_130_t.ll025sc.1996011900_1996011923.nc","r")
  ufile = addfile(era5dir+"/e5.oper.an.pl.128_131_u.ll025uv.1996011900_1996011923.nc","r")
  vfile = addfile(era5dir+"/e5.oper.an.pl.128_132_v.ll025uv.1996011900_1996011923.nc","r")
  qfile = addfile(era5dir+"/e5.oper.an.pl.128_133_q.ll025sc.1996011900_1996011923.nc","r")

  timeix=13
  stride=8
  
  Z3 = zfile->Z(timeix,:,::stride,::stride)
  Z3 = Z3 / 9.81
  U = ufile->U(timeix,:,::stride,::stride)
  V = vfile->V(timeix,:,::stride,::stride)
  T = tfile->T(timeix,:,::stride,::stride)
  Q = qfile->Q(timeix,:,::stride,::stride)
  
  time = zfile->time(timeix)
  print(cd_calendar(time,-3))
  
  printVarSummary(T)

  PMID = T
  PMID = conform(T,tofloat(T&level),0)
  PMID = PMID * 100.

  THETA = pot_temp(PMID, T, -1, False)  ; dim=0

  copy_VarCoords(T,THETA)

  All  = rigrad_bruntv_atm(THETA, U, V, Z3, 3, 0 ) 
  ;RI   = All[0]     ; explicitly extract from list variable  for convenience
  BV   = All[1]     ; must use 'list' syntax  [...]
  ;BUOY = All[2]
  SHR2 = All[3]
  
  WIND = U
  WIND = sqrt(U*U + V*V)
  
  DUDZ = center_finite_diff_n (WIND,Z3,False,0,0)
  DDUDZ  = center_finite_diff_n (DUDZ,Z3,False,0,0)
  DTHETADZ = center_finite_diff_n (THETA,Z3,False,0,0)
  
  ;BV2 = brunt_vaisala_atm(THETA, Z3, 0, 0)  ; opt=0, ndim=1
  
  ;TERM1 = WIND
  ;TERM1 = 9.81 * DTHETADZ / WIND / WIND / THETA 
  
  TERM1 = WIND
  TERM1 = BV * BV / WIND / WIND
  
  TERM2 = WIND
  TERM2 = DDUDZ / WIND

  L2 = WIND
  L2 = TERM1 - TERM2
    
end if

; convert L2 units
L2 = L2 * 1e7
L2 = where(L2 .lt. -10.0,L2@_FillValue,L2)

; pittsburgh
chooseLat=40.4406
chooseLon=-79.9959

;sterling
;chooseLat=39.0067
;chooseLon=-77.4291

if (plot_CAM .eq. False)
  chooseLon = chooseLon + 360.0
end if

Z3prof=Z3(:,{chooseLat},{chooseLon})
L2prof=L2(:,{chooseLat},{chooseLon})
print(Z3prof+" "+L2prof)

system("/bin/rm -f simple.nc")   ; remove any pre-existing file
ncdf = addfile("simple.nc" ,"c")  ; open output netCDF file

;===================================================================
; create global attributes of the file (optional)
;===================================================================
fAtt               = True            ; assign file attributes
fAtt@title         = "NCL Simple Approach to netCDF Creation"
fAtt@source_file   =  "original-file.nc"
fAtt@Conventions   = "None"
fAtt@creation_date = systemfunc ("date")
fileattdef( ncdf, fAtt )            ; copy file attributes

;===================================================================
; make time an UNLIMITED dimension; recommended  for most applications
;===================================================================
filedimdef(ncdf,"time",-1,True) 

;===================================================================
; output variables directly; NCL will call appropriate functions
; to write the meta data associated with each variable
;===================================================================
ncdf->BV = BV
ncdf->DUDZ = DUDZ
ncdf->SHR2 = SHR2


wks   = gsn_open_wks ("x11","xy")               ; send graphics to PNG file

res                  = True                     ; plot mods desired
res@tiMainString     = "Basic XY plot"          ; add title

res@trYMinF  = 0.0
res@trYMaxF  = 14000.0                  ; max value on y-axis

res@trXMinF  = -20.                ; min value on y-axis
res@trXMaxF  = 60.                  ; max value on y-axis
plot  = gsn_csm_xy (wks,L2prof,Z3prof,res) ; create plot


end
