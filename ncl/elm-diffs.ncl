begin


;vars = (/"FSNO","H2OSNO","H2OSNO_TOP","QSNOMELT","QSNOFRZ","SNOW","SNOWDP","SNOWLIQ","SNOTTOPL","ERRH2OSNO","QSNWCPICE","FGR","FSM","QDRIP","ZBOT","TG","TV","TSA","SWup","SWdown","FIRA","FIRE","LWdown","LWup","TREFMNAV","TREFMXAV","Q2M","FPSN","FSH","Rnet","TLAI","FCTR","QOVER","RH","RH2M","SOILWATER_10CM","ZWT","TSOI_10CM","QSOIL","QVEGE","QVEGT","QTOPSOIL","Qstor","H2OCAN","H2OSFC","QH2OSFC","TWS","HC","SNOdTdzL","SNOW_SINKS","SNOW_SOURCES","SNOLIQFL","TH2OSFC","U10","QFLOOD","QRUNOFF","Qle","Qh","Qair","Tair","RAIN","SNO_T","SNO_Z"/)
TVARLIST = (/"Tair","TSOI_10CM","SNOTTOPL","TG","TV","TSA","TREFMNAV","TREFMXAV","TH2OSFC"/)
ZEROLIST = (/"SNOLIQFL","SNOW_SOURCES","SNOW_SINKS","H2OSFC","H2OCAN","FSNO","H2OSNO","H2OSNO_TOP","SNOWLIQ","SNOWDP","SNOW","RAIN","QSNOMELT","SNOWLIQ"/)
out_type = "x11"

;debug
;vars := (/"TWS"/)

; if vars isn't defined above, use the CLI vars
if (.not. isvar("vars")) then
  vars = (/singlevar/)
end if

do zz = 0,dimsizes(vars)-1

  var = vars(zz)
  
  minlat=39.5
  maxlat=42.0
  minlon=280.0-360.
  maxlon=284.0-360.

  print("*** "+var+"")

  exps = (/"PI", "control", "plus1K", "plus2K", "plus3K", "plus4K"/)

  do ii = 0,dimsizes(exps)-1
    filename="/storage/home/cmz5202/scratch/output_"+exps(ii)+"_elm_h0.nc"
    f = addfile(filename,"r")
    v = f->$var$

    ; this code normalizes to a difference from time t = 0
    ;if (var .eq. "TWS") then
    ;  do kk = 1,dimsizes(f->time)-1
    ;    v(kk,:,:) = v(kk,:,:) - v(0,:,:)
    ;  end do
    ;  v(0,:,:) = 0.0
    ;end if
  
    if (ii .eq. 0) then
      fdims = dimsizes(v)
      ntim = fdims(0)
      nlat = fdims(1)
      nlon = fdims(2)
      nexp = dimsizes(exps)
      arr = new((/nexp,ntim,nlat,nlon/),typeof(v))
      time = f->time
    end if
  
    arr(ii,:,:,:) = v
  
    delete(v)
    delete(f)
  
  end do

  ; get a subdomain spatially
  arr_sub = arr(:,:,{minlat:maxlat},{minlon:maxlon})
  delete(arr)
  
  ; time average
  basin_avg = dim_avg_n_Wrap(arr_sub,(/2,3/))

  ; get min/max automatically for plotting
  minvar = min(arr_sub)
  maxvar = max(arr_sub)
  vardiff = maxvar-minvar
  
  minbasin = min(basin_avg)
  maxbasin = max(basin_avg)
  basdiff = maxbasin-minbasin

  if (1 .eq. 1) then
  
    wks1   = gsn_open_wks (out_type, "./img/elm_timeseries_"+var)        ; send graphics to PNG file
  
    res                   = True                      ; plot mods desired
    res@tiMainString      = "Two curve XY plot"       ; add title

    res@xyMonoDashPattern = True
    res@xyLineThicknessF = 5.0
    res@xyLineColors      = (/"blue","green","yellow","orange","red","black"/)          ; change line color

    if (any(TVARLIST .eq. var)) then
      res@gsnYRefLine = 273.15
      res@gsnYRefLineDashPattern = 1
      res@gsnYRefLineThicknessF = 3.0
    end if

    plot  = gsn_csm_xy (wks1,time,basin_avg,res) ; create plot
  
    delete(res)
    delete(plot)
  
  end if
  
  
  if (1 .eq. 1) then

    wks   = gsn_open_wks (out_type, "./img/elm_panels_"+var)        ; send graphics to PNG file

    plot = new(nexp,graphic)

    res                 = True                    ; plot mods desired
    res@gsnDraw         = False                   ; don't draw yet
    res@gsnFrame        = False                   ; don't advance frame yet
    res@gsnAddCyclic = False
    res@lbLabelBarOn = False
    res@cnLineLabelsOn  = False

    res@cnFillMode="RasterFill"
    res@cnFillOn              = True
    res@cnLinesOn   = False
    ;res@cnRasterSmoothingOn = True

    cmap = read_colormap_file("MPL_gnuplot2")
    res@cnFillPalette=cmap
  
    res@cnLevelSelectionMode = "ManualLevels"	; manually set the contour levels with the following 3 resources
    minmaxscaling=0.05
    nlevels=20
    res@cnMinLevelValF  = minvar+minmaxscaling*vardiff
    res@cnMaxLevelValF  = maxvar-minmaxscaling*vardiff
    res@cnLevelSpacingF = vardiff/nlevels

    res@mpMinLatF=minlat
    res@mpMaxLatF=maxlat
    res@mpMinLonF=minlon
    res@mpMaxLonF=maxlon

    res@pmTickMarkDisplayMode = "Always"
    res@mpOutlineOn           = True
    res@mpOutlineBoundarySets = "GeophysicalAndUSStates"
    res@mpDataBaseVersion     = "MediumRes"   ; necessary for mpDataSetName to be effective

    res@gsnLeftString=""

    res2 = True
    res2@gsnDraw         = False                   ; don't draw yet
    res2@gsnFrame        = False                   ; don't advance frame yet
    res2@cnFillOn        = False
    res2@cnLevelSelectionMode = "ExplicitLevels"   ; set explicit contour levels
    
    res2@cnLineThicknessF = 5.0
    res2@gsnLeftString=""
    res2@gsnRightString=""
    res2@gsnCenterString=""
    res2@cnLineLabelsOn = False
    res2@cnInfoLabelOn     = False                   ; turn off info label
    res2@cnLineColor = "green"
    
    do jj = 0,dimsizes(time)-1
      thisTime = jj
      print(""+cd_calendar(time(thisTime),-3))
      do ii = 0,nexp-1
        res@gsnCenterString = exps(ii)
        plot(ii) = gsn_csm_contour_map(wks,arr_sub(ii,thisTime,:,:), res)
        if (any(TVARLIST .eq. var)) then
          res2@cnLevels    = (/ 273.15 /)   ; set levels
          plotover = gsn_csm_contour(wks,arr_sub(ii,thisTime,:,:),res2)
          overlay(plot(ii),plotover)
          delete(plotover)
        end if
        if (any(ZEROLIST .eq. var)) then
          res2@cnLevels    = (/ 0.00001 /)   ; set levels
          plotover = gsn_csm_contour(wks,arr_sub(ii,thisTime,:,:),res2)
          overlay(plot(ii),plotover)
          delete(plotover)
        end if
      end do
    
      panres                     = True
      ;panres@gsnMaximize         = True             ; maximize the plots
      panres@gsnPanelLabelBar    = True             ; turn on the panel label bar
      panres@gsnPanelMainString = tostring(arr_sub@long_name)+"  "+cd_calendar(time(thisTime),-3)
      ;panres@gsnPaperOrientation = "portrait"       ; make sure the paper orientation is set to portrait
      ;panres@lbLabelStride       = 2                ; set the label bar label stride to 2
      gsn_panel(wks,plot,(/3,2/),panres) 
    end do
    
    delete(panres)
    delete(res)
    delete(res2)
    delete(plot)
    
  end if

end do



end