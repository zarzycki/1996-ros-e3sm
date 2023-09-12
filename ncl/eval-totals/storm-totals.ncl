begin

BASEDATADIR="/storage/home/cmz5202/group/arp5873_NEW/"

stday=16
enday=21

ds1 = addfile(BASEDATADIR+"/obs_PRISM/PRISM_daily_SRB_prcp_1996.nc","r")
start_date=cd_inv_calendar(1996,01,stday,00,00,00,ds1->time@units, 0)
end_date=cd_inv_calendar(1996,01,enday,00,00,00,ds1->time@units, 0)
precip_ds1=ds1->PPT({start_date:end_date},:,:)
printVarSummary(precip_ds1)

ds2 = addfile(BASEDATADIR+"/obs_L15/livneh_daily_SRB_prcp_1996.nc","r")
start_date=cd_inv_calendar(1996,01,stday,00,00,00,ds2->time@units, 0)
end_date=cd_inv_calendar(1996,01,enday,00,00,00,ds2->time@units, 0)
precip_ds2=ds2->prec({start_date:end_date},:,:)
printVarSummary(precip_ds2)

ds3 = addfile(BASEDATADIR+"/rean_ERA5/e5.oper.fc.sfc.meanflux.235_055_mtpr.ll025sc.1996010106_1996020106.nc","r")
start_date=cd_inv_calendar(1996,01,stday,00,00,00,ds3->time@units, 0)
end_date=cd_inv_calendar(1996,01,enday,00,00,00,ds3->time@units, 0)
precip_ds3=ds3->MTPR({start_date:end_date},::-1,:)
precip_ds3=precip_ds3*86400/24
precip_ds3!1="lat"
precip_ds3!2="lon"
printVarSummary(precip_ds3)

ds4 = addfile(BASEDATADIR+"/E3SM/model_saved/cut_e3sm_h1_for_eval.nc","r")
start_date=cd_inv_calendar(1996,01,stday,00,00,00,ds4->time@units, 0)
end_date=cd_inv_calendar(1996,01,enday,00,00,00,ds4->time@units, 0)
precip_ds4=ds4->PRECT_indiv(:,{start_date:end_date},:,:)
;precip_ds4:=dim_avg_n_Wrap(precip_ds4,0)
;precip_ds4:=dim_max_n_Wrap(precip_ds4,0)
precip_ds4:=precip_ds4(7,:,:,:)
precip_ds4=precip_ds4*86400*1000/8
printVarSummary(precip_ds4)

ds5 = addfile(BASEDATADIR+"/obs_daymet/daymet_v4_daily_SRB_prcp_1996_RLL.nc","r")
start_date=cd_inv_calendar(1996,01,stday,00,00,00,ds5->time@units, 0)
end_date=cd_inv_calendar(1996,01,enday,00,00,00,ds5->time@units, 0)
precip_ds5=ds5->prcp({start_date:end_date},:,:)
printVarSummary(precip_ds5)

ds6 = addfile(BASEDATADIR+"/obs_cpc/precip.1996.nc","r")
start_date=cd_inv_calendar(1996,01,stday,00,00,00,ds6->time@units, 0)
end_date=cd_inv_calendar(1996,01,enday,00,00,00,ds6->time@units, 0)
precip_ds6=ds6->precip({start_date:end_date},::-1,:)
printVarSummary(precip_ds6)

ds7 = addfile(BASEDATADIR+"/rean_JRA/JRA.h1.1996.PRECT.nc","r")
start_date=cd_inv_calendar(1996,01,stday,00,00,00,ds7->time@units, 0)
end_date=cd_inv_calendar(1996,01,enday,00,00,00,ds7->time@units, 0)
precip_ds7=ds7->PRECT({start_date:end_date},:,:)
precip_ds7=precip_ds7 * 86400 * 1000 / 4
printVarSummary(precip_ds7)

precip_ds1&lon = where(precip_ds1&lon .lt. 0, precip_ds1&lon+360., precip_ds1&lon)
precip_ds2&lon = where(precip_ds2&lon .lt. 0, precip_ds2&lon+360., precip_ds2&lon)
precip_ds3&lon = where(precip_ds3&lon .lt. 0, precip_ds3&lon+360., precip_ds3&lon)
precip_ds4&lon = where(precip_ds4&lon .lt. 0, precip_ds4&lon+360., precip_ds4&lon)
precip_ds5&lon = where(precip_ds5&lon .lt. 0, precip_ds5&lon+360., precip_ds5&lon)
precip_ds6&lon = where(precip_ds6&lon .lt. 0, precip_ds6&lon+360., precip_ds6&lon)
precip_ds7&lon = where(precip_ds7&lon .lt. 0, precip_ds7&lon+360., precip_ds7&lon)

total_precip_ds1 = dim_sum_n_Wrap(precip_ds1,0)
total_precip_ds2 = dim_sum_n_Wrap(precip_ds2,0)
total_precip_ds3 = dim_sum_n_Wrap(precip_ds3,0)
total_precip_ds4 = dim_sum_n_Wrap(precip_ds4,0)
total_precip_ds5 = dim_sum_n_Wrap(precip_ds5,0)
total_precip_ds6 = dim_sum_n_Wrap(precip_ds6,0)
total_precip_ds7 = dim_sum_n_Wrap(precip_ds7,0)

minlat=40.
maxlat=43.
minlon=277.
maxlon=289.
resol=0.5
nlatint=toint((maxlat-minlat)*(1./resol)+1)
nlonint=toint((maxlon-minlon)*(1./resol)+1)
masterLat = fspan(minlat,maxlat,nlatint)
masterLon = fspan(minlon,maxlon,nlonint)

masterLat!0="lat"
masterLat&lat = masterLat
masterLat@units = "degrees_north"

masterLon!0="lon"
masterLon&lon = masterLon
masterLon@units = "degrees_east"

ndata = 7
nlat = dimsizes(masterLat)
nlon = dimsizes(masterLon)

masterArr = new((/ndata,nlat,nlon/),"float")
masterArr!0 = "dset"
masterArr!1 = "lat"
masterArr!2 = "lon"
masterArr&lat = masterLat
masterArr&lon = masterLon
masterArr(0,:,:) = linint2(total_precip_ds1&lon,total_precip_ds1&lat,tofloat(total_precip_ds1),False,masterLon,masterLat,0)
masterArr(1,:,:) = linint2(total_precip_ds2&lon,total_precip_ds2&lat,tofloat(total_precip_ds2),False,masterLon,masterLat,0)
masterArr(2,:,:) = linint2(total_precip_ds3&lon,total_precip_ds3&lat,tofloat(total_precip_ds3),False,masterLon,masterLat,0)
masterArr(3,:,:) = linint2(total_precip_ds4&lon,total_precip_ds4&lat,tofloat(total_precip_ds4),False,masterLon,masterLat,0)
masterArr(4,:,:) = linint2(total_precip_ds5&lon,total_precip_ds5&lat,tofloat(total_precip_ds5),False,masterLon,masterLat,0)
masterArr(5,:,:) = linint2(total_precip_ds6&lon,total_precip_ds6&lat,tofloat(total_precip_ds6),False,masterLon,masterLat,0)
masterArr(6,:,:) = linint2(total_precip_ds7&lon,total_precip_ds7&lat,tofloat(total_precip_ds7),False,masterLon,masterLat,0)

labels = (/"PRISM","L15","ERA5","E3SM","Daymet","CPC","JRA"/)


 wks_xy   = gsn_open_wks("x11","xy")                 ; send graphics to PNG file

 res                   = True                      ; plot mods desired
 res@tiMainString      = "Two curve XY plot"       ; add title

 res@xyLineThicknesses = (/  3.0,3.0,3.0,3.0,3.0,3.0,3.0/)          ; make second line thicker
 res@xyLineColors      = (/"blue","red","green","orange","black","purple","goldenrod"/)          ; change line color

; add a legend
 res@pmLegendDisplayMode    = "Always"              ; turn on legend

 res@pmLegendSide           = "Top"                 ; Change location of
 res@pmLegendParallelPosF   = .90                   ; move units right
 res@pmLegendOrthogonalPosF = -0.8                  ; more neg = down

 res@pmLegendWidthF         = 0.12                  ; Change width and
 res@pmLegendHeightF        = 0.25                  ; height of legend.
 res@lgLabelFontHeightF     = .02                   ; change font height
 res@lgPerimOn              = False                 ; no box around
; labels for the legend
 res@xyExplicitLegendLabels = labels


 plot_xy  = gsn_csm_xy (wks_xy,masterLon,dim_avg_n_Wrap(masterArr,1),res) ; create plot



wks = gsn_open_wks("x11","panel")              ; send graphics to PNG file

cmap = read_colormap_file("WhiteBlue")

plot = new(9,graphic)                          ; create a plot array

res                     = True
res@gsnDraw             = False           ; don't draw
res@gsnFrame            = False           ; don't advance frame
res@cnInfoLabelOn       = False           ; turn off cn info label
res@cnFillOn            = True            ; turn on color
res@lbLabelBarOn        = False           ; turn off individual cb's

res@cnFillPalette         = cmap

res@cnLevelSelectionMode =  "ManualLevels"
res@cnMinLevelValF       = 0.
res@cnMaxLevelValF       = 80.
res@cnLevelSpacingF      = 5.

res@mpMinLatF=minlat
res@mpMaxLatF=maxlat
res@mpMinLonF=minlon
res@mpMaxLonF=maxlon

res@gsnAddCyclic = False


res@cnFillMode = "RasterFill"
res@cnLinesOn =  False             ; Turn off contour lines
res@cnLineLabelsOn  = False

do ii = 0,6
  res@gsnCenterString=labels(ii)
  plot(ii) = gsn_csm_contour_map(wks,masterArr(ii,:,:),res)
end do
;************************************************
; create panel
;************************************************
resP                     = True                ; modify the panel plot
resP@gsnPanelMainString = "A plot with a common label bar"
resP@gsnPanelLabelBar    = True                ; add common colorbar
resP@lbLabelFontHeightF  = 0.007               ; make labels smaller

gsn_panel(wks,plot,(/3,3/),resP)               ; now draw as one plot


exit

system("/bin/rm -f simple.nc")   ; remove any pre-existing file
ncdf = addfile("simple.nc" ,"c")  ; open output netCDF file
fAtt               = True            ; assign file attributes
fAtt@Conventions   = "None"
fAtt@creation_date = systemfunc ("date")
fileattdef( ncdf, fAtt )            ; copy file attributes
filedimdef(ncdf,"time",-1,True)
ncdf->PRECIP  = masterArr


end
