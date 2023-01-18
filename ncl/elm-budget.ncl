load "./utils/shapefile_utils.ncl"




begin
  out_type = "pdf"
  
  ;minlat=39.5
  ;maxlat=42.0
  ;minlon=280.0-360.
  ;maxlon=284.0-360.
  minlat=35.5
  maxlat=47.0
  minlon=275.0-360.
  maxlon=290.0-360.

  exps = (/"PI", "control", "plus1K", "plus2K", "plus3K", "plus4K"/)
  ;exps = (/"control"/)
  filedir="/Volumes/ZARZYCKI_FLASH/arp5873/ens_means/"
  
  do ii = 0,dimsizes(exps)-1
    filename=filedir+"/output_"+exps(ii)+"_elm_h0.nc"
    f = addfile(filename,"r")

    time = f->time
    print(cd_calendar(time,-3))
    
    qrain = f->RAIN(:,:,:)
    qsnow = f->SNOW(:,:,:)
    QSOIL = f->QSOIL(:,:,:)
    QVEGE = f->QVEGE(:,:,:)
    QVEGT = f->QVEGT(:,:,:)
    EvEg = QSOIL
    EvEg = QSOIL + QVEGE + QVEGT
    
    qover = f->QOVER(:,:,:)
    qh2osfc = f->QH2OSFC(:,:,:)
    qrunoff = f->QRUNOFF(:,:,:)
    qrunoff = qrunoff - qover - qh2osfc
    qsnwcpice = f->QSNWCPICE(:,:,:)
    qrunoff = qrunoff + qsnwcpice

    h2osfc = f->H2OSFC(:,:,:)
    h2osno = f->H2OSNO(:,:,:)  
    
    BUDGET = qrain
    BUDGET = qrain + qsnow - EvEg - qrunoff - qover - qh2osfc
    
    BUDGET_int = BUDGET
    BUDGET_int = 0.0
    RAIN_int = BUDGET_int
    SNOW_int = BUDGET_int
    EvEg_int = BUDGET_int
    qrunoff_int = BUDGET_int
    qh2osfc_int = BUDGET_int
    qover_int = BUDGET_int
    h2osfc_int = BUDGET_int
    h2osno_int = BUDGET_int

    do kk = 1,dimsizes(time)-1
      BUDGET_int(kk,:,:) = BUDGET_int((kk-1),:,:) + BUDGET((kk-1),:,:)*21600.
      RAIN_int(kk,:,:) = RAIN_int((kk-1),:,:) + qrain((kk-1),:,:)*21600.
      SNOW_int(kk,:,:) = SNOW_int((kk-1),:,:) + qsnow((kk-1),:,:)*21600.
      EvEg_int(kk,:,:) = EvEg_int((kk-1),:,:) + EvEg((kk-1),:,:)*21600.
      qrunoff_int(kk,:,:) = qrunoff_int((kk-1),:,:) + qrunoff((kk-1),:,:)*21600.
      qh2osfc_int(kk,:,:) = qh2osfc_int((kk-1),:,:) + qh2osfc((kk-1),:,:)*21600.
      qover_int(kk,:,:) = qover_int((kk-1),:,:) + qover((kk-1),:,:)*21600.
      
      ; these don't need integration
      h2osno_int(kk,:,:) = h2osno(kk,:,:) - h2osno(0,:,:)
      h2osfc_int(kk,:,:) = h2osfc(kk,:,:) - h2osfc(0,:,:)
    end do
    
    BUDGET_int@long_name = "Storage change"
    BUDGET_int@units     = "mm"

    if (ii .eq. 0) then
      fdims = dimsizes(BUDGET)
      ntim = fdims(0)
      nlat = fdims(1)
      nlon = fdims(2)
      nexp = dimsizes(exps)
      arr = new((/10,nexp,ntim,nlat,nlon/),typeof(BUDGET))
      arr!2="time"
      arr&time=f->time
      arr!3="lat"
      arr&lat=f->lat
      arr!4="lon"
      arr&lon=f->lon
    end if
  
    arr(3,ii,:,:,:) = RAIN_int
    arr(4,ii,:,:,:) = SNOW_int
    arr(5,ii,:,:,:) = -EvEg_int
    arr(7,ii,:,:,:) = -qh2osfc_int
    arr(6,ii,:,:,:) = -qover_int
    arr(8,ii,:,:,:) = -qrunoff_int
    arr(0,ii,:,:,:) = h2osno_int
    arr(1,ii,:,:,:) = h2osfc_int    
    arr(2,ii,:,:,:) = (BUDGET_int - h2osfc_int - h2osno_int)
  end do

  ; get a subdomain spatially
  ;arr_sub = arr(:,:,:,{minlat:maxlat},{minlon:maxlon})

  shpfilename="/Users/cmz5202/Software/1996-ros-e3sm/shapes/srb.shp"
  shp = addfile(shpfilename,"r")
  print_shapefile_info(shpfilename)
  plot_shapefile(shpfilename)
  
  if (min(arr&lon) .gt. 0.0) then
    print("Need to flip lon")
    lon = arr&lon
    newlon = lon - 360.0
    arr&lon=newlon
    delete(newlon)
    delete(lon)
  end if

  Opt=True
  Opt@debug=True
  Opt@return_mask=True
  ; get 2-D mask of 1s and 0s
  MASK_2D = shapefile_mask_data(arr(0,0,0,:,:),shpfilename,Opt)
  MASK = conform (arr,MASK_2D,(/3,4/))
  
  arr_sub = arr
  arr_sub = where(MASK .gt. 0.5,arr_sub,arr_sub@_FillValue)

 ;system("/bin/rm -f simple.nc")   ; remove any pre-existing file
 ;ncdf = addfile("simple.nc" ,"c")  ; open output netCDF file
 ;ncdf->MASK_2D  = MASK_2D           
 ;ncdf->arr_sub  = arr_sub 
 ;ncdf->arr  = arr      

       
  ; time average
  basin_avg = dim_avg_n_Wrap(arr_sub,(/3,4/))



  
    wks   = gsn_open_wks (out_type, "./img/elm_timeseries_BUDGET")        ; send graphics to PNG file
  
    plot = new(9,graphic)
  
    res                   = True                      ; plot mods desired
        res@gsnDraw         = False                   ; don't draw yet
    res@gsnFrame        = False                   ; don't advance frame yet
    
    
    
    ;res@tiMainString      = "Two curve XY plot"       ; add title

    res@xyMonoDashPattern = True
    res@xyLineThicknessF = 4.0
    
    ; abby's colors
  cmap = (/ (/000, 183, 182/), \
            (/026, 000, 248/), \
            (/000, 117, 022/), \
            (/255, 000, 028/), \
            (/190, 000, 179/), \
            (/182, 185, 041/) /)/255.
             
    res@xyLineColors      = cmap

  res@trYMinF  = -100.                   ; min value on y-axis
  res@trYMaxF  =  100.                   ; max value on y-axis

  res@tmXMajorGrid                = True          ; implement x grid 
  res@tmXMajorGridThicknessF      = 0.5           ; 2.0 is default
  res@tmXMajorGridLineDashPattern = 2             ; select short dash lines
  res@tmYMajorGrid                = True          ; implement y grid
  res@tmYMajorGridThicknessF      = 0.5           ; 2.0 is default
  res@tmYMajorGridLineDashPattern = 2             ; select short dash lines
  res@tmXMajorGridLineColor = "gray50"
  res@tmYMajorGridLineColor = "gray50"
  res@tiXAxisString     = "" 
  res@tmGridDrawOrder        = "Draw" 
  
  res@tmXBMode          = "Explicit"              ; explicit labels
  res@tmXBValues        = (/0.5,1.5,2.5,3.5,4.5,5.5,6.5/)
  res@tmXBLabels        = (/"1/16","1/17","1/18","1/19","1/20","1/21","1/22"/)
  
  ; turn off top tick marks
  res@tmXTOn           = False
  
    res@tiYAxisString     = "~F8~D~F14~Q~B~snow~N~"
    plot(0)  = gsn_csm_xy (wks,time,basin_avg(0,:,:),res) ; create plot
    res@tiYAxisString     = "~F8~D~F14~Q~B~sfc~N~"  
    plot(1)  = gsn_csm_xy (wks,time,basin_avg(1,:,:),res) ; create plot
    res@tiYAxisString     = "~F34~;~F8~D~F14~Q~B~soil~N~"  
    plot(2)  = gsn_csm_xy (wks,time,basin_avg(2,:,:),res) ; create plot
    
    ; middle row
    res@tiYAxisString     = "~F14~q~B~rain~N~~F8~D~F14~t"  
    plot(3)  = gsn_csm_xy (wks,time,basin_avg(3,:,:),res) ; create plot
    res@tiYAxisString     = "~F14~q~B~snow~N~~F8~D~F14~t"  
    plot(4)  = gsn_csm_xy (wks,time,basin_avg(4,:,:),res) ; create plot
    res@tiYAxisString     = "~F14~(E~B~v~N~+E~B~g~N~)"
    plot(5)  = gsn_csm_xy (wks,time,basin_avg(5,:,:),res) ; create plot
    
    ;bottom row
    res@tiXAxisString     = "Date (00Z, 1996)"
    res@tiYAxisString     = "~F14~q~B~over~N~~F8~D~F14~t"  
    plot(6)  = gsn_csm_xy (wks,time,basin_avg(6,:,:),res) ; create plot
    res@tiYAxisString     = "~F14~q~B~h2osfc~N~~F8~D~F14~t"
    plot(7)  = gsn_csm_xy (wks,time,basin_avg(7,:,:),res) ; create plot
    
    
    res@xyExplicitLegendLabels = (/"Pre_Ind", "Control", "Plus_1K", "Plus_2K", "Plus_3K", "Plus_4K"/)
    res@lgItemOrder  = (/ 5,4,3,2,1,0 /)   ; Reorder the legends
    res@pmLegendDisplayMode    = "Always"   ; Display a legend.
    res@lgPerimOn              = True      ; No legend perimeter.
    res@pmLegendOrthogonalPosF = -0.61            ; Move legend inside plot
    res@pmLegendParallelPosF   = 0.76            ; Move legend to right
    res@pmLegendWidthF         = 0.24             ; Change width and height
    res@pmLegendHeightF        = 0.24
    res@lgBoxMinorExtentF      = 0.4       ; Shorten the legend lines
    res@lgPerimFill            = "SolidFill"      ; Fill legend box w/white
    res@lgPerimFillColor       = "White"          ; so it masks XY curves
    res@tiYAxisString     = "~F14~q~B~other~N~~F8~D~F14~t"
    plot(8)  = gsn_csm_xy (wks,time,basin_avg(8,:,:),res) ; create plot
    
    ; uncomment to debug budget
    ; move from LHS to RHS and sum to make sure all terms = 0
    ;basin_avg(6,:,:) = -basin_avg(6,:,:)
    ;basin_avg(7,:,:) = -basin_avg(7,:,:)
    ;basin_avg(8,:,:) = -basin_avg(8,:,:)
    ;check_sum = dim_sum_n_Wrap(basin_avg,0)
    ;plot(8)  = gsn_csm_xy (wks,time,check_sum,res) ; create plot
    

  
      resP                     = True
      letters=(/"a.","b.","c.","d.","e.","f.","g.","h.","i.","j.","k.","l.","m.","n.","o.","p.","q.","r.", "s.", "t.","u.","v.","w.","x.","y.","z.","aa.","bb."/)
      resP@gsnPanelFigureStringsFontHeightF = 0.02
      resP@gsnPanelFigureStrings=letters
      resP@amJust           = "TopLeft"
      ;panres@gsnMaximize         = True             ; maximize the plots
      ;panres@gsnPanelLabelBar    = True             ; turn on the panel label bar
      ;panres@gsnPanelMainString = tostring(arr_sub@long_name)+"  "+cd_calendar(time(thisTime),-3)
      ;panres@gsnPaperOrientation = "portrait"       ; make sure the paper orientation is set to portrait
      ;panres@lbLabelStride       = 2                ; set the label bar label stride to 2
      gsn_panel(wks,plot,(/3,3/),resP) 
      

  
end

