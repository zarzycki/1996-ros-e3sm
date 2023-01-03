begin


vars = (/"FSNO","H2OSNO","H2OSNO_TOP","QSNOMELT","QSNOFRZ","SNOW","SNOWDP","SNOWLIQ","SNOTTOPL","ERRH2OSNO","QSNWCPICE","FGR","FSM","QDRIP","ZBOT","TG","TV","TSA","SWup","SWdown","FIRA","FIRE","LWdown","LWup","TREFMNAV","TREFMXAV","Q2M","FPSN","FSH","Rnet","TLAI","FCTR","QOVER","RH","RH2M","SOILWATER_10CM","ZWT","TSOI_10CM","QSOIL","QVEGE","QVEGT","QTOPSOIL","Qstor","H2OCAN","H2OSFC","QH2OSFC","TWS","HC","SNOdTdzL","SNOW_SINKS","SNOW_SOURCES","SNOLIQFL","TH2OSFC","U10","QFLOOD","QRUNOFF","Qle","Qh","Qair","Tair","RAIN","SNO_T","SNO_Z"/)

do zz = 0,dimsizes(vars)-1

  var = vars(zz)
  
  print("*** "+var+"")

  ;var = "TSOI_10CM"
  minlat=39.5
  maxlat=42.0
  minlon=280.0-360.
  maxlon=284.0-360.


  exps = (/"PI", "control", "plus1K", "plus2K", "plus3K", "plus4K"/)

  do ii = 0,dimsizes(exps)-1
    filename="/storage/home/cmz5202/scratch/output_"+exps(ii)+"_elm_h0.nc"
    f = addfile(filename,"r")
    v = f->$var$
  
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

  thisTime = 18
  print(cd_calendar(time(thisTime),-3))

  arr_sub = arr(:,22,{minlat:maxlat},{minlon:maxlon})
  printVarSummary(arr_sub)
  delete(arr)

  minvar = min(arr_sub)
  maxvar = max(arr_sub)
  vardiff = maxvar-minvar

  wks   = gsn_open_wks ("png", "./img/elm_"+var)        ; send graphics to PNG file

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
  res@cnRasterSmoothingOn = True

  cmap = read_colormap_file("MPL_gnuplot")
  res@cnFillPalette=cmap
  
  res@cnLevelSelectionMode = "ManualLevels"	; manually set the contour levels with the following 3 resources
  res@cnMinLevelValF  = minvar+0.1*vardiff
  res@cnMaxLevelValF  = maxvar-0.1*vardiff
  res@cnLevelSpacingF = vardiff/20

  res@mpMinLatF=minlat
  res@mpMaxLatF=maxlat
  res@mpMinLonF=minlon
  res@mpMaxLonF=maxlon

  res@pmTickMarkDisplayMode = "Always"
  res@mpOutlineOn           = True
  res@mpOutlineBoundarySets = "GeophysicalAndUSStates"
  res@mpDataBaseVersion     = "MediumRes"   ; necessary for mpDataSetName to be effective

  do ii = 0,nexp-1
    res@gsnCenterString = exps(ii)
    plot(ii) = gsn_csm_contour_map(wks,arr_sub(ii,:,:), res)
  end do

  panres                     = True
  panres@gsnMaximize         = True             ; maximize the plots
  panres@gsnPanelLabelBar    = True             ; turn on the panel label bar
  panres@gsnPanelMainString = var
  ;panres@gsnPaperOrientation = "portrait"       ; make sure the paper orientation is set to portrait
  ;panres@lbLabelStride       = 2                ; set the label bar label stride to 2
  gsn_panel(wks,plot,(/3,2/),panres) 

end do



end