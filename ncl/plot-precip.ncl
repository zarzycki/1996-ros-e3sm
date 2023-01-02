begin

minlat=28.5
maxlat=47.0
minlon=-102.0
maxlon=-65.0
  
plot_CAM=True

if (plot_CAM) then
  f = addfile("/storage/home/cmz5202/scratch/output000_eam_h1.nc","r")
  timeix=30
  stride=1
  PRECT = f->PRECT_indiv(:,timeix,::stride,::stride)

  PRECT = PRECT*3600000
  PRECT@units="mm/hr"
  PRECT@long_name = "Precipitation rate"
  

  print(cd_calendar(f->time,-3))
  time = f->time(timeix)
  print(cd_calendar(time,-3))
  
else

  print("hh")
  
end if

; plotting

wks   = gsn_open_wks ("x11", "coneff")        ; send graphics to PNG file

nplots=4
plot = new(nplots,graphic)


res                 = True                    ; plot mods desired
res@gsnDraw         = False                   ; don't draw yet
res@gsnFrame        = False                   ; don't advance frame yet
res@gsnAddCyclic = False

res@cnFillMode="RasterFill"
res@cnFillOn              = True
res@cnLinesOn   = False
res@cnRasterSmoothingOn = True

PRECT = where(PRECT .lt. 0.08,PRECT@_FillValue,PRECT)


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

res@cnFillPalette=cmap
res@cnLevelSelectionMode = "ExplicitLevels"                        ; set explicit contour levels
res@cnLevels    = (/0.1,0.2,0.3,0.4,0.5,0.75,1.0,1.25,1.5,1.75,2.0,2.25,2.5,3.0,3.5,4.0,4.5,5.0,6.0,8.0,10.0,15.0,20.0,25.0,35.0/)        ; set levels

res@mpMinLatF=minlat
res@mpMaxLatF=maxlat
res@mpMinLonF=minlon
res@mpMaxLonF=maxlon

res@pmTickMarkDisplayMode = "Always"
res@mpOutlineOn           = True
res@mpOutlineBoundarySets = "GeophysicalAndUSStates"
res@mpDataBaseVersion     = "MediumRes"   ; necessary for mpDataSetName to be effective
res@mpLandFillColor       = "gray50"
res@mpOceanFillColor = "blue"
res@mpInlandWaterFillColor = "blue"

plot(0) = gsn_csm_contour_map(wks,PRECT(2,:,:), res)
plot(1) = gsn_csm_contour_map(wks,PRECT(3,:,:), res)
plot(2) = gsn_csm_contour_map(wks,PRECT(4,:,:), res)
plot(3) = gsn_csm_contour_map(wks,PRECT(9,:,:), res)



panres                     = True
panres@gsnMaximize         = True             ; maximize the plots
;panres@gsnPanelLabelBar    = True             ; turn on the panel label bar
;panres@gsnPaperOrientation = "portrait"       ; make sure the paper orientation is set to portrait
;panres@lbLabelStride       = 2                ; set the label bar label stride to 2
gsn_panel(wks,plot,(/nplots,1/),panres) 


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