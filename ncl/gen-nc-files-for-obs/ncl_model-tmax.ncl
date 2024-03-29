load "../utils/shapefile_utils.ncl"

begin

minlat=39.0
maxlat=43.5
minlon=-79.0
maxlon=-74.0

;minmax="min"

BASEPATH=getenv("BASEPATH")
SOFTPATH=getenv("SOFTPATH")
  
f=addfile(BASEPATH+"/ens_means/output_control_eam_h0.nc","r")

time = f->time

startdate = cd_inv_calendar(1996,01,15,0,0,0,time@units, 0)
enddate = cd_inv_calendar(1996,01,23,0,0,0,time@units, 0)

var = f->TREFHT({startdate:enddate},{minlat:maxlat},{minlon:maxlon})

newtime = var&time
newtime = newtime - 0.25
var&time = newtime

opt = True
opt@nval_crit = 4    ; require at least 4 values per day for a daily mean to be calculated.
SWE = calculate_daily_values (var, minmax, 0, opt)

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

print(SWE_avg+" "+cd_calendar(SWE_avg&time,-3))

system("/bin/rm -f mod_t"+minmax+".nc")   ; remove any pre-existing file
ncdf = addfile("mod_t"+minmax+".nc" ,"c")  ; open output netCDF file

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
