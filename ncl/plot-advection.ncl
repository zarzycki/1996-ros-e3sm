begin

minlat=28.5
maxlat=47.0
minlon=-102.0
maxlon=-65.0

pnew = (/ 850. /)       ; desired output levels 
  
plot_CAM=False

if (plot_CAM) then
  f = addfile("/storage/home/cmz5202/scratch/arp-precip/output000_eam_h2.nc","r")
  timeix=16
  stride=4
  PS = f->PS(timeix,::stride,::stride)
  T = f->T(timeix,:,::stride,::stride)
  U = f->U(timeix,:,::stride,::stride)
  V = f->V(timeix,:,::stride,::stride)

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

  Tnew = vinth2p (T,hyam,hybm,pnew,PS,1,P0mb,1,True)
  Unew = vinth2p (U,hyam,hybm,pnew,PS,1,P0mb,1,True)
  Vnew = vinth2p (V,hyam,hybm,pnew,PS,1,P0mb,1,True)

  delete(U)
  delete(V)
  delete(T)
  delete(PS)

else

  day=19
  era5dir="/storage/home/cmz5202/group/arp5873/ERA5/"
  tfile = addfile(era5dir+"/e5.oper.an.pl.128_130_t.ll025sc.199601"+day+"00_199601"+day+"23.nc","r")
  ufile = addfile(era5dir+"/e5.oper.an.pl.128_131_u.ll025uv.199601"+day+"00_199601"+day+"23.nc","r")
  vfile = addfile(era5dir+"/e5.oper.an.pl.128_132_v.ll025uv.199601"+day+"00_199601"+day+"23.nc","r")

  timeix=12
  stride=2
  
  U = ufile->U(timeix,{pnew},::stride,::stride)
  V = vfile->V(timeix,{pnew},::stride,::stride)
  T = tfile->T(timeix,{pnew},::stride,::stride)
  
  print(cd_calendar(tfile->time,-3))
  time = tfile->time(timeix)
  print(cd_calendar(time,-3))

  printVarSummary(T)
  
  Tnew = T(::-1,:)
  Unew = U(::-1,:)
  Vnew = V(::-1,:)
  
  lat = Tnew&latitude
  lon = Tnew&longitude
    
end if

opt_adv = 0  ; use 1 to get gradient components
cyclic = True
long_name = "temp advection"
units     = "K/s"   
advTcfd   = advect_variable_cfd(Unew,Vnew,Tnew,lat,lon, cyclic, long_name, units, opt_adv)
;arr = advect_variable(Unew,Vnew,Tnew,1,long_name,units,opt_adv)
;advTcfd = -arr[0]
;dTdx = arr[1]
;dTdy = arr[2]

advTcfd = -advTcfd
; plotting

wks   = gsn_open_wks ("x11", "coneff")        ; send graphics to PNG file

res                 = True                    ; plot mods desired
res@gsnDraw         = False                   ; don't draw yet
res@gsnFrame        = False                   ; don't advance frame yet
res@gsnAddCyclic = False

res@cnFillMode="RasterFill"
res@cnFillOn              = True
res@cnLinesOn   = False
res@cnRasterSmoothingOn = True

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
  
;res@cnLevelSelectionMode = "ManualLevels"	; manually set the contour levels with the following 3 resources
;res@cnMinLevelValF  = -0.001			; set the minimum contour level
;res@cnMaxLevelValF  =  0.001			; set the maximum contour level
;res@cnLevelSpacingF =  0.0001			; set the interval between contours

plot = new(1,graphic)

plot(0) = gsn_csm_contour_map(wks,Tnew, res)

;res2 = True
;res2@gsnDraw         = False                   ; don't draw yet
;res2@gsnFrame        = False                   ; don't advance frame yet
;res2@cnFillOn              = False
;res2@cnMinLevelValF  = 200			; set the minimum contour level
;res2@cnMaxLevelValF  = 325			; set the maximum contour level
;res2@cnLevelSpacingF = 5.0			; set the interval between contours
;res2@cnLineThicknessF = 1.0
;
;plot_td = gsn_csm_contour(wks,Tnew(:,:),res2)
;
;overlay(plot(0),plot_td)

  vcres                         = True             ; plot mods desired
  vcres@gsnDraw         = False                   ; don't draw yet
  vcres@gsnFrame                = False            ; so we can draw time stamp
  vcres@vcRefAnnoOrthogonalPosF = -1.0             ; move ref vector up
  vcres@vcRefMagnitudeF         = 40.0             ; define vector ref mag
  vcres@vcRefLengthF            = 0.045            ; define length of vec ref
  ;vcres@vcGlyphStyle            = "CurlyVector"    ; turn on curly vectors
  vcres@vcMinDistanceF          = 0.017

  vector = gsn_csm_vector(wks,Unew,Vnew,vcres)

overlay(plot(0),vector)


panres                     = True
panres@gsnMaximize         = True             ; maximize the plots
;panres@gsnPanelLabelBar    = True             ; turn on the panel label bar
;panres@gsnPaperOrientation = "portrait"       ; make sure the paper orientation is set to portrait
;panres@lbLabelStride       = 2                ; set the label bar label stride to 2
gsn_panel(wks,plot,(/1,1/),panres) 
 

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
ncdf->T = Tnew      
ncdf->U = Unew   
ncdf->V = Vnew   
ncdf->advTcfd = advTcfd
end