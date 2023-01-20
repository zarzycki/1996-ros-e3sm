load "../utils/shapefile_utils.ncl"

begin

minlat=39.0
maxlat=43.5
minlon=-79.0
maxlon=-74.0
;data="JRA"

BASEPATH=getenv("BASEPATH")
SOFTPATH=getenv("SOFTPATH")

if (data .eq. "obs") then
  f=addfile(BASEPATH+"/obs_AZsnow/4km_SWE_Depth_WY1996_v01.nc","r")
  time = f->time
  startdate = cd_inv_calendar(1996,01,15,0,0,0,time@units, 0)
  enddate = cd_inv_calendar(1996,01,23,0,0,0,time@units, 0)
  SWE = short2flt( f->SWE({startdate:enddate},{minlat:maxlat},{minlon:maxlon}) )
else if (data .eq. "obs_snowdepth") then
  f=addfile(BASEPATH+"/obs_AZsnow/4km_SWE_Depth_WY1996_v01.nc","r")
  time = f->time
  startdate = cd_inv_calendar(1996,01,10,0,0,0,time@units, 0)
  enddate = cd_inv_calendar(1996,01,23,0,0,0,time@units, 0)
  SWE = short2flt( f->DEPTH({startdate:enddate},{minlat:maxlat},{minlon:maxlon}) )
else if (data .eq. "NLDAS") then
  f=addfile(BASEPATH+"/obs_NLDAS/NLDAS_2_VIC_199601.nc","r")
  time = f->time
  startdate = cd_inv_calendar(1996,01,15,0,0,0,time@units, 0)
  enddate = cd_inv_calendar(1996,01,23,0,0,0,time@units, 0)
  SWE = f->WEASD_GDS0_SFC({startdate:enddate},{minlat:maxlat},{minlon:maxlon})
  SWE!0 = "time"
  SWE!1 = "lat"
  SWE!2 = "lon"
else if (data .eq. "L15") then
  f=addfile("../L15_Fluxes_199601.nc","r")
  time = f->time
  startdate = cd_inv_calendar(1996,01,15,0,0,0,time@units, 0)
  enddate = cd_inv_calendar(1996,01,23,0,0,0,time@units, 0)
  SWE = f->SWE({startdate:enddate},{minlat:maxlat},{minlon:maxlon})
else if (data .eq. "JRA") then
  f=addfile(BASEPATH+"/rean_JRA/JRA.h1.1996.SWE.nc","r")
  time = f->time
  startdate = cd_inv_calendar(1996,01,15,0,0,0,time@units, 0)
  enddate = cd_inv_calendar(1996,01,23,0,0,0,time@units, 0)
  minlon=minlon+360.0
  maxlon=maxlon+360.0
  SWE = f->SWE({startdate:enddate},{minlat:maxlat},{minlon:maxlon})
else if (data .eq. "ERA5") then
  f=addfile(BASEPATH+"/rean_ERA5/e5.oper.an.sfc.128_141_sd.ll025sc.1996010100_1996013123.nc","r")
  time = f->time
  startdate = cd_inv_calendar(1996,01,15,0,0,0,time@units, 0)
  enddate = cd_inv_calendar(1996,01,23,0,0,0,time@units, 0)
  minlon=minlon+360.0
  maxlon=maxlon+360.0
  SWE = f->SD({startdate:enddate},{minlat:maxlat},{minlon:maxlon})
  SWE!0 = "time"
  SWE!1 = "lat"
  SWE!2 = "lon"
  SWE = SWE*1000. ; convert m to mm
else
  exit
end if
end if
end if
end if
end if
end if

shpfilename=SOFTPATH+"/shapes/srb.shp"
shp = addfile(shpfilename,"r")
print_shapefile_info(shpfilename)
plot_shapefile(shpfilename)

if (min(SWE&lon) .gt. 0.0) then
  print("Need to flip lon")
  lon = SWE&lon
  newlon = lon - 360.0
  SWE&lon=newlon
  delete(newlon)
  delete(lon)
end if

Opt=True
Opt@debug=True
Opt@return_mask=True
; get 2-D mask of 1s and 0s
MASK_2D = shapefile_mask_data(SWE(0,:,:),shpfilename,Opt)
MASK = conform (SWE,MASK_2D,(/1,2/))

SWE_sub = SWE
SWE_sub = where(MASK .gt. 0.5,SWE_sub,SWE_sub@_FillValue)

SWE_avg = dim_avg_n_Wrap(SWE_sub,(/1,2/))

print(SWE_avg+" "+SWE_avg&time+" "+cd_calendar(SWE_avg&time,-3))

newtime = SWE_avg&time
if (data .eq. "ERA5" .or. data .eq. "JRA") then
  newtime = newtime
else
  newtime = newtime + 0.5
end if
SWE_avg&time = newtime

system("/bin/rm -f "+data+"_swe.nc")   ; remove any pre-existing file
ncdf = addfile(data+"_swe.nc" ,"c")  ; open output netCDF file

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
ncdf->SWE_avg = SWE_avg

end
