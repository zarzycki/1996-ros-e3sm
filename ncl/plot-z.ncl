begin

minlat=28.5
maxlat=47.0
minlon=-102.0
maxlon=-65.0

pnew = (/ 850., 500. /)       ; desired output levels 
  
plot_CAM=False


if (plot_CAM) then
  f = addfile("/storage/home/cmz5202/scratch/output_eam_h2.nc","r")
  timeix=16
  stride=4
  Z3 = f->Z3(timeix,:,::stride,::stride)
  PS = f->PS(timeix,::stride,::stride)
  T = f->T(timeix,:,::stride,::stride)
  Q = f->Q(timeix,:,::stride,::stride)
  U = f->U(timeix,:,::stride,::stride)
  V = f->V(timeix,:,::stride,::stride)
  ;W = mixhum_convert(Q, "q", (/0,0/))
  ;copy_VarCoords(Q,W)

  print(cd_calendar(f->time,-3))
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

  RH = relhum (T, Q, PMID)  ; rel humidity
  TD = dewtemp_trh(T,RH)
  delete(RH)

  THETA = pot_temp(PMID, T, -1, False)  ; dim=0
  delete(PMID)

  copy_VarCoords(T,TD)
  copy_VarCoords(T,THETA)
  delete(T)

  All  = rigrad_bruntv_atm(THETA, U, V, Z3, 3, 1 ) 
  ;RI   = All[0]     ; explicitly extract from list variable  for convenience
  BV   = All[1]     ; must use 'list' syntax  [...]
  ;BUOY = All[2]
  SHR2 = All[3]

  printVarSummary(SHR2)
  printVarSummary(BV)
  printVarSummary(U)
  delete(All)

  ;L2 = BV
  ;L2 = BV / SHR2 

  Znew = vinth2p (Z3,hyam,hybm,pnew,PS,1,P0mb,1,True)
  TDnew = vinth2p (TD,hyam,hybm,pnew,PS,1,P0mb,1,True)
  Unew = vinth2p (U,hyam,hybm,pnew,PS,1,P0mb,1,True)
  Vnew = vinth2p (V,hyam,hybm,pnew,PS,1,P0mb,1,True)

  delete(U)
  delete(V)
  delete(Z3)
  delete(TD)

else
  day=19
  era5dir="/storage/home/cmz5202/group/arp5873/ERA5/"
  zfile = addfile(era5dir+"/e5.oper.an.pl.128_129_z.ll025sc.199601"+day+"00_199601"+day+"23.nc","r")
  tfile = addfile(era5dir+"/e5.oper.an.pl.128_130_t.ll025sc.199601"+day+"00_199601"+day+"23.nc","r")
  ufile = addfile(era5dir+"/e5.oper.an.pl.128_131_u.ll025uv.199601"+day+"00_199601"+day+"23.nc","r")
  vfile = addfile(era5dir+"/e5.oper.an.pl.128_132_v.ll025uv.199601"+day+"00_199601"+day+"23.nc","r")
  qfile = addfile(era5dir+"/e5.oper.an.pl.128_133_q.ll025sc.199601"+day+"00_199601"+day+"23.nc","r")

  timeix=12
  stride=2
  
  Z3 = zfile->Z(timeix,{pnew},::stride,::stride)
  Z3 = Z3 / 9.81
  U = ufile->U(timeix,{pnew},::stride,::stride)
  V = vfile->V(timeix,{pnew},::stride,::stride)
  T = tfile->T(timeix,{pnew},::stride,::stride)
  Q = qfile->Q(timeix,{pnew},::stride,::stride)
  
  print(cd_calendar(zfile->time,-3))
  time = zfile->time(timeix)
  print(cd_calendar(time,-3))

  RH = relhum (T, Q, conform(T,pnew*100.,0) )  ; rel humidity
  TD = dewtemp_trh(T,RH)
  copy_VarCoords(T,TD)
  
  print("T:  "+max(T)+" "+min(T))
  print("Q:  "+max(Q)+" "+min(Q))
  print("TD: "+max(TD)+" "+min(TD))
  print("RH: "+max(RH)+" "+min(RH))
  
  Znew = Z3
  TDnew = TD
  Unew = U
  Vnew = V
    
end if





; plotting

wks   = gsn_open_wks ("x11", "coneff")        ; send graphics to PNG file

res                 = True                    ; plot mods desired
res@gsnDraw         = False                   ; don't draw yet
res@gsnFrame        = False                   ; don't advance frame yet
res@gsnAddCyclic = False

res@cnLineLabelsOn      = False        ; Turn on contour line labels

;res@mpMinLatF=minlat
;res@mpMaxLatF=maxlat
;res@mpMinLonF=minlon
;res@mpMaxLonF=maxlon
res@mpProjection = "LambertConformal"
res@mpLimitMode       = "Corners"         ; choose range of map
res@mpLeftCornerLatF  = minlat
res@mpLeftCornerLonF  = minlon
res@mpRightCornerLatF = maxlat
res@mpRightCornerLonF = maxlon
res@mpLambertMeridianF = -96.

res@pmTickMarkDisplayMode = "Always"
res@mpOutlineOn           = True
res@mpOutlineBoundarySets = "GeophysicalAndUSStates"
res@mpDataBaseVersion     = "MediumRes"   ; necessary for mpDataSetName to be effective
res@mpLandFillColor       = "bisque2"
  
res@cnLevelSelectionMode = "ManualLevels"	; manually set the contour levels with the following 3 resources
;res@cnMinLevelValF  = 5160.			; set the minimum contour level
;res@cnMaxLevelValF  = 5760.			; set the maximum contour level
;res@cnLevelSpacingF = 60.			; set the interval between contours
res@cnMinLevelValF  = 1100.			; set the minimum contour level
res@cnMaxLevelValF  = 1600.			; set the maximum contour level
res@cnLevelSpacingF = 30.			; set the interval between contours

plot = new(1,graphic)

plot(0) = gsn_csm_contour_map(wks,Znew({850.},:,:), res)

res2 = res
res2@cnMinLevelValF  = 273.15			; set the minimum contour level
res2@cnMaxLevelValF  = 283.15			; set the maximum contour level
res2@cnLevelSpacingF = 5.0			; set the interval between contours
res2@cnLineThicknessF = 4.0

plot_td = gsn_csm_contour(wks,TDnew({850.},:,:),res2)

overlay(plot(0),plot_td)


panres                     = True
panres@gsnMaximize         = True             ; maximize the plots
;panres@gsnPanelLabelBar    = True             ; turn on the panel label bar
;panres@gsnPaperOrientation = "portrait"       ; make sure the paper orientation is set to portrait
;panres@lbLabelStride       = 2                ; set the label bar label stride to 2
gsn_panel(wks,plot,(/1,1/),panres) 

 
 exit
 

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
ncdf->Z = Znew      
ncdf->TD = TDnew
       
end