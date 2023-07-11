load "/Users/cmz5202/Software/1996-ros-e3sm/ncl/rainratetodbz.ncl"

begin

yyyy_str="1996"
mm_str="01"
;dd_str="18"
debug=False

print(dd_str)

mydate=yyyy_str+mm_str+dd_str

;DATADIR="/Users/cmz5202/NetCDF/arp5873_NEW/gage/"
gauge_data = DATADIR+"/raw_data/gage.hrly.prcp."+mydate

print("Loading "+gauge_data)

gauge_lines = asciiread(gauge_data,-1,"string")

years=toint(str_squeeze(str_get_cols(gauge_lines,0,3)))
months=toint(str_squeeze(str_get_cols(gauge_lines,5,6)))
days=toint(str_squeeze(str_get_cols(gauge_lines,8,9)))
hours=toint(str_squeeze(str_get_cols(gauge_lines,11,12)))
minutes=toint(str_squeeze(str_get_cols(gauge_lines,14,15)))
siteids= str_squeeze(str_get_cols(gauge_lines,17,24))
precips=tofloat(str_squeeze(str_get_cols(gauge_lines,29,36)))

ngaugeobs=dimsizes(years)

;     1 -  4      Year     -- Time
;      6 -  7      Month    -- of
;      9 - 10      Day      -- the
;     12 - 13      Hour     -- Observation
;     15 - 16      Minute   -- (All digits, all times are GMT)
;     18 - 25      Site ID
;     27 - 28      Element Code  (PP=Precipitation)
;     30 - 37      Element Value (Real number with 2 decimal places)
;     39 - 39      Revision Code (0=No, 1=Yes)
;     42 - 45      Duration Code (2001=1 24-hour day,  5004=Time period
;                                 begining at 7 am local time and ending
;                                 at the time of observation)
;     47 - 47      Type Code     (F=Forcast, P=Processed, R=Reading)
;     49 - 49      Source Code   (G=GOES, M=Meteor Burst, P=Phone, Z=None)
;     51 - 51      Extreme Code  (N=Minimum of day, X=Maximum of Day, Z=None)
;     53 - 60      Data Source   (ID of NWS office sending report, if known)
;     62 - 62      Changed       (0=No, 1=Yes) Has the original value been
;                                 changed? Only 1 update is allowed.
;     64 - 64      Quality Flag  (Z=None, B=Buddy Check, C=Climatology Check,
;                                 R=Bad Zero Check, E=Estimated value,
;                                 V=Value verified, M=Manual overide of B/C/R)
if debug then
  print(years+" "+months+" "+days+" "+hours+" "+minutes+" "+siteids+" "+precips)
end if

filename = "processed-stations.csv"

lines = asciiread(filename,-1,"string")

delim = ","
station_id1  = str_squeeze(str_get_field(lines,1,delim))
station_id2  = str_squeeze(str_get_field(lines,2,delim))
station_lat  = tofloat(str_get_field(lines,3,delim))
station_lon  = tofloat(str_get_field(lines,4,delim))

;print(station_lat+" "+station_lon)
;print(ind(ismissing(station_lat)))

low   = -0.00001
high  =  0.00001
unf   = random_uniform(low, high, dimsizes(station_lon))
station_lat = station_lat + unf
station_lon = station_lon - unf


station_lat!0="ncol"
station_lon!0="ncol"

if debug then
  print("Name is '" + station_id1 + "', group is " + station_id2 + ", size is " + station_lat)
end if










nstations = dimsizes(station_id1)
ntimes = 24

precip = new((/ntimes,nstations/),"float")
precip!0 = "time"
precip!1 = "ncol"
precip = 0.0

times = fspan(0.0,23.0,24)
time_units = "hours since "+yyyy_str+"-"+mm_str+"-"+dd_str+" 00:00"
times!0 = "time"
times&time = times
times@units=time_units

hits = new((/ntimes,nstations/),"integer")
hits = 0

do ii = 0,ngaugeobs-1
  if ( mod(ii,500).eq.0) then
    print ("ii = "+ii+"/"+(ngaugeobs-1))
  end if
  thissta=siteids(ii)
  ix = str_match_ind(station_id1,thissta)
  if (dimsizes(ix) .gt. 1) then
    print("uh oh!")
    continue
  end if
  if (ismissing(ix)) then
    ix = str_match_ind(station_id2,thissta)
    if (ismissing(ix)) then
      continue
    end if
  end if

  thisTime = cd_inv_calendar(years(ii),months(ii),days(ii),hours(ii),minutes(ii),0,time_units, 0)
  timeix = ind_nearest_coord(thisTime,times,0)

  precip(timeix,ix) = (/precips(ii)/)
  hits(timeix,ix) = 1

end do

resol=0.25
minlat=29.0
maxlat=53.0
minlon=-100.0
maxlon=-70.0
lats = fspan(minlat,maxlat,toint((maxlat-minlat)/resol) + 1)
lons = fspan(minlon,maxlon,toint((maxlon-minlon)/resol) + 1)
lats!0 = "lat"
lons!0 = "lon"
lons@units="degrees_east"
lats@units="degrees_north"
lats&lat = lats
lons&lon = lons
;zgrd  = natgrid(station_lon, station_lat, precip, lons,lats)    ;  longitude x latitude

zgrd = new((/ntimes,dimsizes(lons),dimsizes(lats)/),"float")
zgrd!0 = "time"
zgrd!1 = "lon"
zgrd!2 = "lat"
zgrd&time = times
zgrd&lon = lons
zgrd&lat = lats

printVarSummary(zgrd)

;opt  = True           ; method=1 and domain=1 are set (domain was set to 0 in NCL V6.3.0 and earlier)
;opt@distmx = 50000.      ; kilometers
;zgrd := zgrd(time|:,lat|:,lon|:)           ; x_order(lat,lon,time)
;dssetp("exp", 2.5)

;nnsetp("igr", 0)
;nnsetp("non", 1)
;nnsetp("bI", 1.1)
;nnsetp("bJ", 7.0)

do jj = 0,ntimes-1
  print(jj+"")
  ixs := ind(hits(jj,:) .eq. 1)
  ;natgrid dsgrid2
  ;zgrd(jj,:,:) = (/natgrid(station_lon(ixs), station_lat(ixs), precip(jj,ixs), lons, lats)/)
  zgrd(jj,:,:) = (/dsgrid2(station_lon(ixs), station_lat(ixs), precip(jj,ixs), lons, lats)/)

  ;Opt                = True
  ;Opt@ForceOverwrite = True
  ;Opt@PrintTimings   = True
  ;Opt@Debug = True
  ;unstructured_to_ESMF("tmp.nc",station_lat(ixs),station_lon(ixs),Opt)

  print("done")
end do

x_reorder = zgrd(time|:,lat|:,lon|:)           ; x_order(lat,lon,time)

cut_sub=True
if cut_sub then
  x_reorder := x_reorder(:,{39.:44.},{-80.:-74.})
end if

system("/bin/rm -f "+DATADIR+"/ds_simple_"+mydate+".nc")   ; remove any pre-existing file
ncdf = addfile(DATADIR+"/ds_simple_"+mydate+".nc" ,"c")  ; open output netCDF file

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
ncdf->precip = x_reorder


end