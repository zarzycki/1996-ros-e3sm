load "rainratetodbz.ncl"

begin

minlat=30.0
maxlat=50.0
minlon=-95.0
maxlon=-70.0

pnew = (/ 850., 500. /)       ; desired output levels 
  
plot_dbz=True

panel_vertical=False

out_type="pdf"

yyyy=1996
mm=01
days = (/18,19,19,20/)
hours = (/12,00,12,00/)
bools = (/True,False/)

wks   = gsn_open_wks (out_type, "FIG_synoptic_comparison")        ; send graphics to PNG file

nmodels = 2
npanels = dimsizes(days)
plot = new(nmodels*npanels,graphic)

do bb = 0,nmodels-1
do aa = 0,npanels-1

plot_CAM=bools(bb)

dd=days(aa)
hh=hours(aa)

era5dir="/Volumes/ZARZYCKI_FLASH/arp5873/rean_ERA5/"
CAMdir="/Volumes/ZARZYCKI_FLASH/arp5873/ens_means/"
;era5dir="/storage/home/cmz5202/group/arp5873/ERA5/"
;CAMdir="/storage/home/cmz5202/scratch/"

if (plot_CAM) then
  f = addfile(CAMdir+"/output_control_eam_h2.nc","r")
    
  stride=4
  Z3 = f->Z3({cd_inv_calendar(yyyy,mm,dd,hh,0,0,f->time@units,-1)},:,::stride,::stride)
  PS = f->PS({cd_inv_calendar(yyyy,mm,dd,hh,0,0,f->time@units,-1)},::stride,::stride)
  T = f->T({cd_inv_calendar(yyyy,mm,dd,hh,0,0,f->time@units,-1)},:,::stride,::stride)
  Q = f->Q({cd_inv_calendar(yyyy,mm,dd,hh,0,0,f->time@units,-1)},:,::stride,::stride)
  U = f->U({cd_inv_calendar(yyyy,mm,dd,hh,0,0,f->time@units,-1)},:,::stride,::stride)
  V = f->V({cd_inv_calendar(yyyy,mm,dd,hh,0,0,f->time@units,-1)},:,::stride,::stride)

  g = addfile(CAMdir+"/output_control_eam_h0.nc","r")
  PRECT = g->PRECT({cd_inv_calendar(yyyy,mm,dd,hh,0,0,g->time@units,-1)},::stride,::stride)
  if (.not. plot_dbz) then
    PRECT = PRECT*3600000
    PRECT@units="mm/hr"
    PRECT@long_name = "Precipitation rate"
  else
    PRECT = rainratetodbz(PRECT,1)
    PRECT@units="dbZ"
    PRECT@long_name = "Simulated radar reflectivity"
  end if
  PSL = g->PSL({cd_inv_calendar(yyyy,mm,dd,hh,0,0,f->time@units,-1)},::stride,::stride)
  
  ;W = mixhum_convert(Q, "q", (/0,0/))
  ;copy_VarCoords(Q,W)

  hyam = f->hyam                   ; read to memory (optional)
  hybm = f->hybm
  P0mb = 1000.      
  P0 = f->P0
  lat = f->lat(::stride)
  lon = f->lon(::stride)

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
  Tnew = vinth2p (T,hyam,hybm,pnew,PS,1,P0mb,1,True)
  PSLnew = PSL

  delete(U)
  delete(V)
  delete(Z3)
  delete(TD)
  delete(T)
  delete(Q)
  delete(PSL)

else


  zfile = addfile(era5dir+"/e5.oper.an.pl.128_129_z.ll025sc.199601"+dd+"00_199601"+dd+"23.nc","r")
  tfile = addfile(era5dir+"/e5.oper.an.pl.128_130_t.ll025sc.199601"+dd+"00_199601"+dd+"23.nc","r")
  ufile = addfile(era5dir+"/e5.oper.an.pl.128_131_u.ll025uv.199601"+dd+"00_199601"+dd+"23.nc","r")
  vfile = addfile(era5dir+"/e5.oper.an.pl.128_132_v.ll025uv.199601"+dd+"00_199601"+dd+"23.nc","r")
  qfile = addfile(era5dir+"/e5.oper.an.pl.128_133_q.ll025sc.199601"+dd+"00_199601"+dd+"23.nc","r")

  timeix=12
  stride=2
  
  Z3 = zfile->Z({cd_inv_calendar(yyyy,mm,dd,hh,0,0,zfile->time@units,-1)},{pnew},::stride,::stride)
  Z3 = Z3 / 9.81
  U = ufile->U({cd_inv_calendar(yyyy,mm,dd,hh,0,0,ufile->time@units,-1)},{pnew},::stride,::stride)
  V = vfile->V({cd_inv_calendar(yyyy,mm,dd,hh,0,0,vfile->time@units,-1)},{pnew},::stride,::stride)
  T = tfile->T({cd_inv_calendar(yyyy,mm,dd,hh,0,0,tfile->time@units,-1)},{pnew},::stride,::stride)
  Q = qfile->Q({cd_inv_calendar(yyyy,mm,dd,hh,0,0,qfile->time@units,-1)},{pnew},::stride,::stride)
  
  lat = zfile->latitude(::stride)
  lon = zfile->longitude(::stride)
  
  g = addfile(era5dir+"/e5.oper.fc.sfc.meanflux.235_055_mtpr.ll025sc.1996010106_1996020106.nc","r")
  PRECT_3hr = g->MTPR({cd_inv_calendar(yyyy,mm,dd,hh+1,0,0,g->time@units,-1):cd_inv_calendar(yyyy,mm,dd,hh+4,0,0,g->time@units,-1)},::stride,::stride)
  PRECT = dim_avg_n_Wrap(PRECT_3hr,0)
  PRECT = PRECT / 1000.  ; convert from kg/m2/s to m/s to match CAM
  if (.not. plot_dbz) then
    PRECT = PRECT*3600000
    PRECT@units="mm/hr"
    PRECT@long_name = "Precipitation rate"
  else
    PRECT = rainratetodbz(PRECT,1)
    PRECT@units="dbZ"
    PRECT@long_name = "Simulated radar reflectivity"
  end if
  
  hfile = addfile(era5dir+"/e5.oper.an.sfc.128_151_msl.ll025sc.1996010100_1996013123.nc","r")
  PSL = hfile->MSL({cd_inv_calendar(yyyy,mm,dd,hh,0,0,hfile->time@units,-1)},::stride,::stride)

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
  Tnew = T
  PSLnew = PSL
  
  delete(U)
  delete(V)
  delete(Z3)
  delete(TD)
  delete(T)
  delete(Q)
  delete(PSL)
    
end if



;opt_adv = 0  ; use 1 to get gradient components
;cyclic = True
;long_name = "temp advection"
;units     = "K/s"   
;advTcfd   = advect_variable_cfd(Unew,Vnew,Tnew,lat,lon, cyclic, long_name, units, opt_adv)
;advTcfd = -advTcfd






resbase                 = True                    ; plot mods desired
resbase@gsnDraw         = False                   ; don't draw yet
resbase@gsnFrame        = False                   ; don't advance frame yet
resbase@gsnLeftString=""
resbase@gsnCenterString=""
resbase@gsnRightString=""
resbase@cnInfoLabelOn     = False

res = resbase

res@gsnAddCyclic = False

res@lbLabelBarOn        = False           ; turn off individual cb's

res@mpMinLatF=minlat
res@mpMaxLatF=maxlat
res@mpMinLonF=minlon
res@mpMaxLonF=maxlon

;res@mpProjection = "LambertConformal"
;res@mpLimitMode       = "Corners"         ; choose range of map
;res@mpLeftCornerLatF  = minlat
;res@mpLeftCornerLonF  = minlon
;res@mpRightCornerLatF = maxlat
;res@mpRightCornerLonF = maxlon
;res@mpLambertMeridianF = -96.

res@pmTickMarkDisplayMode = "Always"
res@mpOutlineOn           = True
res@mpOutlineBoundarySets = "GeophysicalAndUSStates"
res@mpProvincialLineColor = "grey25"
res@mpNationalLineColor   = "grey25"
res@mpGeophysicalLineColor= "grey25"
res@mpDataBaseVersion     = "MediumRes"   ; necessary for mpDataSetName to be effective
res@mpLandFillColor       = "bisque2"
;res@mpOutlineDrawOrder      = "PreDraw"  ; draw contours first

res@cnFillOn    = True                        ; turn on color fill
res@cnLinesOn   = False                       ; turn off the contour lines
res@cnLineLabelsOn = False                    ; turn the line labels off

if (plot_CAM) then
  datastr="E3SM"
else
  datastr="ERA5"
end if

res@gsnCenterString=datastr+": "+tostring_with_format(yyyy, "%0.4i")+"-"+tostring_with_format(mm, "%0.2i")+"-"+tostring_with_format(days(aa), "%0.2i")+" "+tostring_with_format(hours(aa), "%0.2i")+"Z"

res@cnLevelSelectionMode = "ExplicitLevels"
if (plot_dbz) then
  PRECT = where(PRECT .lt. 2.0,PRECT@_FillValue,PRECT)
  res@cnLevels    = (/10,12.5,15,20,25,30,35,40,45,50,55,60,70/)
  cmap = read_colormap_file("radar")
  cmap_dims = dimsizes(cmap)
  ncolors = cmap_dims(0)
  npairs = cmap_dims(1)
  nnewcolors = 200
  cmap_interp = new((/nnewcolors,npairs/),typeof(cmap))
  do ii = 0,npairs-1
   ;cmap_interp(3:nnewcolors-1,ii) = linint1 (fspan(0,1,ncolors-3), cmap(3:ncolors-1,ii), False, fspan(0,1,nnewcolors-3), 0)
   ;cmap_interp(2:nnewcolors-1,ii) = linint1 (fspan(0,1,ncolors-2), cmap(2:ncolors-1,ii), False, fspan(0,1,nnewcolors-2), 0)
   cmap_interp(:,ii) = linint1 (fspan(0,1,ncolors), cmap(:,ii), False, fspan(0,1,nnewcolors), 0)
  end do
  cmap_interp(0,:) = cmap(0,:)
  cmap_interp(1,:) = cmap(1,:)
else
  PRECT = where(PRECT .lt. 0.08,PRECT@_FillValue,PRECT)
  res@cnLevels    = (/0.1,0.2,0.3,0.4,0.5,0.75,1.0,1.25,1.5,1.75,2.0,2.25,2.5,3.0,3.5,4.0,4.5,5.0,6.0,8.0,10.0,15.0,20.0,25.0,35.0/)
  cmap = read_colormap_file("prcp_1")
  cmap_dims = dimsizes(cmap)
  ncolors = cmap_dims(0)
  npairs = cmap_dims(1)
  nnewcolors = 200
  cmap_interp = new((/nnewcolors,npairs/),typeof(cmap))
  do ii = 0,npairs-1
   ;cmap_interp(3:nnewcolors-1,ii) = linint1 (fspan(0,1,ncolors-3), cmap(3:ncolors-1,ii), False, fspan(0,1,nnewcolors-3), 0)
   ;cmap_interp(2:nnewcolors-1,ii) = linint1 (fspan(0,1,ncolors-2), cmap(2:ncolors-1,ii), False, fspan(0,1,nnewcolors-2), 0)
   cmap_interp(:,ii) = linint1 (fspan(0,1,ncolors), cmap(:,ii), False, fspan(0,1,nnewcolors), 0)
  end do
  cmap_interp(0,:) = cmap(0,:)
  cmap_interp(1,:) = cmap(1,:)
  ;cmap_interp((nnewcolors-1),:) = cmap((ncolors-1),:)
  ;cmap_interp((nnewcolors-2),:) = cmap((ncolors-2),:)
end if
res@cnFillPalette=cmap

if (panel_vertical) then
  plotix = (aa*nmodels)+bb
else
  plotix = (bb*npanels)+aa
end if

print(plotix+"")

plot(plotix) = gsn_csm_contour_map(wks,PRECT(:,:), res)

; plot PSL?
if 0 .eq. 1 then
  res3 = resbase
  res3@cnLineDrawOrder ="PostDraw"
  res3@cnLevelSelectionMode = "ManualLevels"	; manually set the contour levels with the following 3 resources
  res3@cnMinLevelValF  = 95000.
  res3@cnMaxLevelValF  = 104000.
  res3@cnLevelSpacingF = 800.0
  res3@cnLineThicknessF = 2.5
  res3@cnLineColor = "grey15"
  res3@cnLineLabelsOn = False        ; Turn on contour line labels
  res3@cnMonoLineDashPattern =True
  res3@cnLineDashPattern = 15
  plot_psl = gsn_csm_contour(wks,PSLnew,res3)
  overlay(plot(plotix),plot_psl)
end if

res1 = resbase
res1@cnLineDrawOrder ="PostDraw"
res1@cnLevelSelectionMode = "ManualLevels"	; manually set the contour levels with the following 3 resources
;res1@cnMinLevelValF  = 1100.			; set the minimum contour level
;res1@cnMaxLevelValF  = 1600.			; set the maximum contour level
;res1@cnLevelSpacingF = 30.			; set the interval between contours
res1@cnMinLevelValF  = 5160.			; set the minimum contour level
res1@cnMaxLevelValF  = 5760.			; set the maximum contour level
res1@cnLevelSpacingF = 60.			; set the interval between contours
res1@cnLineColor = "black"
res1@cnLineLabelsOn      = False        ; Turn on contour line labels
res1@cnLineThicknessF = 3.0
plot_z = gsn_csm_contour(wks,Znew({500.},:,:), res1)
overlay(plot(plotix),plot_z)

res2 = resbase
res2@cnLineDrawOrder ="PostDraw"
res2@cnLevelSelectionMode = "ManualLevels"	; manually set the contour levels with the following 3 resources
res2@cnMinLevelValF  = 273.
res2@cnMaxLevelValF  = 400.
res2@cnLevelSpacingF = 5.0
res2@cnLineThicknessF = 4.0
res2@cnLineColor = "red"
res2@cnLineLabelsOn      = False        ; Turn on contour line labels
res2@cnMonoLineDashPattern =True
res2@cnLineDashPattern = 15
do ii = 0,19
  TDnew({850.},:,:) = smth9_Wrap(TDnew({850.},:,:), 0.50,  0.25, False)
end do
plot_td = gsn_csm_contour(wks,TDnew({850.},:,:),res2)

overlay(plot(plotix),plot_td)



; buncha cleanup
delete(Unew)
delete(Vnew)
delete(Znew)
delete(TDnew)
delete(Tnew)
delete(PSLnew)
delete(lat)
delete(lon)
delete(PRECT)

end do ; aa loop
end do ; bb loop

panres                     = True
;panres@gsnMaximize         = True             ; maximize the plots
panres@gsnPanelLabelBar    = True             ; turn on the panel label bar
;panres@gsnPaperOrientation = "portrait"       ; make sure the paper orientation is set to portrait
;panres@lbLabelStride       = 2                ; set the label bar label stride to 2

panres@amJust           = "TopLeft"
panres@gsnPanelFigureStringsFontHeightF = 0.014
panres@gsnPanelFigureStrings= (/"a.","b.","c.","d.","e.","f.","g.","h.","i.","j."/)

if (panel_vertical) then
  gsn_panel(wks,plot,(/npanels,nmodels/),panres) 
else
  gsn_panel(wks,plot,(/nmodels,npanels/),panres) 
end if

 
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