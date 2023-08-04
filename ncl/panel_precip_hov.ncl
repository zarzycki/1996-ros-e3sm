load "./utils/rainratetodbz.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/time_axis_labels.ncl"

begin

dbzthresh=0.1
minLat=  40.5
maxLat=  42.
minLon= -79.2
maxLon= -74.8
; 18-21
stday=18
enday=21
;DATADIR="/Users/cmz5202/NetCDF/arp5873_NEW/"

; start npanels counter
npanels=0

; ERA5
f=addfile(DATADIR+"/rean_ERA5/e5.oper.fc.sfc.meanflux.235_055_mtpr.ll025sc.1996010106_1996020106.nc","r")
time := f->time
sttime := cd_inv_calendar(1996,1,stday,00,0,0,time@units, 0)
entime := cd_inv_calendar(1996,1,enday,00,0,0,time@units, 0)
hov_data_era5 = f->MTPR({sttime:entime},{39:43},{275:290})
hov_data_era5 = hov_data_era5 * 60.*60.
hov_data_era5 := rainratetodbz(hov_data_era5,0)
hov_data_era5 = where(hov_data_era5 .lt. dbzthresh,hov_data_era5@_FillValue,hov_data_era5)
npanels = npanels+1

; ERA5
g=addfile(DATADIR+"/ens_means/output_control_eam_h0.nc","r")
time := g->time
sttime := cd_inv_calendar(1996,1,stday,00,0,0,time@units, 0)
entime := cd_inv_calendar(1996,1,enday,00,0,0,time@units, 0)
hov_data_e3sm = g->PRECT({sttime:entime},{39:43},{275-360.:290-360.})
hov_data_e3sm = hov_data_e3sm * 3.6e6
hov_data_e3sm := rainratetodbz(hov_data_e3sm,0)
hov_data_e3sm = where(hov_data_e3sm .lt. dbzthresh,hov_data_e3sm@_FillValue,hov_data_e3sm)
npanels = npanels+1

; JRA
h=addfile(DATADIR+"/rean_JRA/JRA.h1.1996.PRECT.nc","r")
sttime = cd_inv_calendar(1996,1,stday,00,0,0,"hours since 1800-01-01 00:00", 0)
entime = cd_inv_calendar(1996,1,enday,00,0,0,"hours since 1800-01-01 00:00", 0)
hov_data_jra = h->PRECT({sttime:entime},{39:43},{275:290})
hov_data_jra = hov_data_jra* 3.6e6
hov_data_jra := rainratetodbz(hov_data_jra,0)
hov_data_jra = where(hov_data_jra .lt. dbzthresh,hov_data_jra@_FillValue,hov_data_jra)
npanels = npanels+1

; OBS
i=addfile(DATADIR+"/obs_gage/gage_ds_catted.nc","r")
time := i->time
sttime = cd_inv_calendar(1996,1,stday,00,0,0,time@units, 0)
entime = cd_inv_calendar(1996,1,enday,00,0,0,time@units, 0)
hov_data_obs = i->precip({sttime:entime},{39:43},{275-360.:290-360.})
hov_data_obs = hov_data_obs*25.4
hov_data_obs := rainratetodbz(hov_data_obs,0)
hov_data_obs = where(hov_data_obs .lt. dbzthresh,hov_data_obs@_FillValue,hov_data_obs)
npanels = npanels+1

; PRISM
j=addfile(DATADIR+"/obs_PRISM/prism_catted.nc","r")
time := j->time
sttime = cd_inv_calendar(1996,1,stday,00,0,0,time@units, 0)
entime = cd_inv_calendar(1996,1,enday,00,0,0,time@units, 0)
hov_data_prism = j->PPT({sttime:entime},{39:43},{275-360.:290-360.})
;hov_data_prism&time = hov_data_prism&time+1.0
hov_data_prism = hov_data_prism/24.
hov_data_prism := rainratetodbz(hov_data_prism,0)
hov_data_prism = where(hov_data_prism .lt. dbzthresh,hov_data_prism@_FillValue,hov_data_prism)
npanels = npanels+1

; CPC
k=addfile(DATADIR+"/obs_cpc/precip.1996.nc","r")
time := k->time
sttime = cd_inv_calendar(1996,1,stday,00,0,0,time@units, 0)
entime = cd_inv_calendar(1996,1,enday,00,0,0,time@units, 0)
hov_data_cpc = k->precip({sttime:entime},{39:43},{275.:290.})
hov_data_cpc = hov_data_cpc/24.
hov_data_cpc := rainratetodbz(hov_data_cpc,0)
hov_data_cpc = where(hov_data_cpc .lt. dbzthresh,hov_data_cpc@_FillValue,hov_data_cpc)
npanels = npanels+1

; L15
;l=addfile(DATADIR+"/obs_L15/L15_Meteo_199601.nc","r")
;time := l->time
;sttime = cd_inv_calendar(1996,1,stday-1,00,0,0,time@units, 0)
;entime = cd_inv_calendar(1996,1,enday-1,00,0,0,time@units, 0)
;hov_data_l15 = l->Prec({sttime:entime},{39:43},{275-360.:290-360.})
;hov_data_l15&time = hov_data_l15&time+1.0
;hov_data_l15 = hov_data_l15/24.
;hov_data_l15 := rainratetodbz(hov_data_l15,0)
;hov_data_l15 = where(hov_data_l15 .lt. dbzthresh,hov_data_l15@_FillValue,hov_data_l15)
;npanels = npanels+1

;==============================================
; create color plot
;=============================================

wks  = gsn_open_wks ("x11", "FIG_preciphov" )           ; send graphics to PNG file

print("npanels: "+npanels)

plot = new(npanels,graphic)                          ; create a plot array

res                      = True               ; plot mods desired
;res@gsnMaximize = True
;res@gsnPaperOrientation = "landscape"

res@cnFillOn             = True               ; turn on color fill
res@vpWidthF      = 0.2
res@vpHeightF     = 0.7
res@gsnDraw             = False           ; don't draw
res@gsnFrame            = False           ; don't advance frame
res@cnInfoLabelOn       = False           ; turn off cn info label
res@lbLabelBarOn        = False           ; turn off individual cb's

res@cnLinesOn        =  False             ; Turn off contour lines
res@cnFillMode = "RasterFill"

res@gsnLeftString=""
res@gsnRightString=""

res@cnLevelSelectionMode = "ExplicitLevels"
;res@cnLevels    = (/10,12.5,15,20,25,30,35,40,45,50,55,60,70/)
;res@cnLevels    = (/10,12.5,15,20,25,30,35,40,45,50,60,70/)
;res@cnLevels    = (/10,12.5,15,20,25,30,35,40,45,50,60,70/)
res@cnLevels    = (/5,7.5,10,12.5,15,20,25,30,35,40,45,50,60/)
cmap = read_colormap_file("prcp_1")
cmap_dims = dimsizes(cmap)
ncolors = cmap_dims(0)
npairs = cmap_dims(1)
nnewcolors = 200
cmap_interp = new((/nnewcolors,npairs/),typeof(cmap))
do ii = 0,npairs-1
 cmap_interp(:,ii) = linint1 (fspan(0,1,ncolors), cmap(:,ii), False, fspan(0,1,nnewcolors), 0)
end do
cmap_interp(0,:) = cmap(0,:)
cmap_interp(1,:) = cmap(1,:)
res@cnFillPalette=cmap

;res@tmXBMode        = "Manual"
;res@tmXBTickStartF  = -80.
;res@tmXBTickEndF    = -72.
;res@tmXBTickSpacingF= 0.5

 ;res@tmXBMode           = "Explicit"         ; Define own tick mark labels.
 ;theseticks = fspan(-80.,-72.,9)
 ;res@tmXBValues         = theseticks
 ;res@tmXBLabels         = theseticks

res@gsnMajorLonSpacing = 2.0
res@gsnMinorLonSpacing = 1.0
res@gsnCenterStringFontHeightF = 0.022
res@tmXTOn = False
res@tmXBLabelFontHeightF = 0.018
res@tmYLLabelFontHeightF = 0.018

; control tick marks
restick = True
restick@ttmFormat = "%d"
restick@ttmAxis = "YL"
turn_off_internal_ylabels=False

res@gsnCenterString="ERA5"
res_era5 = res
res_era5_tick = restick
res_era5_tick@ttmMajorStride = 24
tmp_era5_time = tofloat(hov_data_era5&time)
copy_VarMeta(hov_data_era5&time,tmp_era5_time)
time_axis_labels( tmp_era5_time,  res_era5, res_era5_tick )
if (turn_off_internal_ylabels) then
  res_era5@tmYLLabelsOn = False
end if
plot(4) = gsn_csm_hov(wks, dim_avg_n_Wrap(hov_data_era5(:,{minLat:maxLat},{360.+minLon:360.+maxLon}),1), res_era5)

res@gsnCenterString="E3SM"
res_e3sm = res
res_e3sm_tick = restick
res_e3sm_tick@ttmMajorStride = 8
tmp_e3sm_time = tofloat(hov_data_e3sm&time)
copy_VarMeta(hov_data_e3sm&time,tmp_e3sm_time)
time_axis_labels( tmp_e3sm_time,  res_e3sm, res_e3sm_tick )
plot(0) = gsn_csm_hov(wks, dim_avg_n_Wrap(hov_data_e3sm(:,{minLat:maxLat},{minLon:maxLon}),1), res_e3sm)

res@gsnCenterString="JRA"
res_jra = res
res_jra_tick = restick
res_jra_tick@ttmMajorStride = 4
tmp_jra_time = tofloat(hov_data_jra&time)
copy_VarMeta(hov_data_jra&time,tmp_jra_time)
time_axis_labels( tmp_jra_time,  res_jra, res_jra_tick )
if (turn_off_internal_ylabels) then
  res_jra@tmYLLabelsOn = False
end if
plot(5) = gsn_csm_hov(wks, dim_avg_n_Wrap(hov_data_jra(:,{minLat:maxLat},{360.+minLon:360.+maxLon}),1), res_jra)

res@gsnCenterString="GAGE"
res_obs = res
res_obs_tick = restick
res_obs_tick@ttmMajorStride = 24
tmp_obs_time = tofloat(hov_data_obs&time)
copy_VarMeta(hov_data_obs&time,tmp_obs_time)
time_axis_labels( tmp_obs_time,  res_obs, res_obs_tick )
if (turn_off_internal_ylabels) then
  res_obs@tmYLLabelsOn = False
end if
plot(1) = gsn_csm_hov(wks, dim_avg_n_Wrap(hov_data_obs(:,{minLat:maxLat},{minLon:maxLon}),1), res_obs)

res@gsnCenterString="PRISM"
res_prism = res
res_prism_tick = restick
res_prism_tick@ttmMajorStride = 1
tmp_prism_time = tofloat(hov_data_prism&time)
copy_VarMeta(hov_data_prism&time,tmp_prism_time)
time_axis_labels( tmp_prism_time,  res_prism, res_prism_tick )
if (turn_off_internal_ylabels) then
  res_prism@tmYLLabelsOn = False
end if
plot(3) = gsn_csm_hov(wks, dim_avg_n_Wrap(hov_data_prism(:,{minLat:maxLat},{minLon:maxLon}),1), res_prism)

res@gsnCenterString="CPC"
res_cpc = res
res_cpc_tick = restick
res_cpc_tick@ttmMajorStride = 1
tmp_cpc_time = tofloat(hov_data_cpc&time)
copy_VarMeta(hov_data_cpc&time,tmp_cpc_time)
time_axis_labels( tmp_cpc_time,  res_cpc, res_cpc_tick )
if (turn_off_internal_ylabels) then
  res_cpc@tmYLLabelsOn = False
end if
plot(2) = gsn_csm_hov(wks, dim_avg_n_Wrap(hov_data_cpc(:,{minLat:maxLat},{360.+minLon:360.+maxLon}),1), res_cpc)

;res@gsnCenterString="L15"
;res_l15 = res
;res_l15_tick = restick
;res_l15_tick@ttmMajorStride = 1
;tmp_l15_time = tofloat(hov_data_l15&time)
;copy_VarMeta(hov_data_l15&time,tmp_l15_time)
;time_axis_labels( tmp_l15_time,  res_l15, res_l15_tick )
;plot(6) = gsn_csm_hov(wks, dim_avg_n_Wrap(hov_data_l15(:,{minLat:maxLat},{minLon:maxLon}),1), res_l15)

;************************************************
; create panel
;************************************************
resP                     = True                ; modify the panel plot
;resP@gsnPanelMainString = "Hovmoller"

letters=(/"a.","b.","c.","d.","e.","f.","g.","h.","i.","j.","k.","l.","m.","n.","o.","p.","q.","r.", "s.", "t.","u.","v.","w.","x.","y.","z.","aa.","bb."/)
resP@gsnPanelFigureStringsFontHeightF = 0.02
resP@gsnPanelFigureStrings=letters
resP@amJust           = "TopLeft"

resP@pmLabelBarWidthF     = 0.8
resP@gsnPanelLabelBar    = True                ; add common colorbar
resP@lbLabelFontHeightF  = 0.01               ; make labels smaller

gsn_panel(wks,plot,(/1,npanels/),resP)               ; now draw as one plot

;system("pdftk FIG_preciphov.pdf cat 1-endeast output FIG_preciphov.pdf")

end
