begin

DATADIR="/Users/cmz5202/NetCDF/arp5873_NEW/obs_PRISM/"

fili = systemfunc("ls "+DATADIR+"/PRISM_PPT_*.nc | sort -n")
print(fili)
f    = addfiles(fili,"r")

PPT = f[:]->PPT
printVarSummary(PPT)

time = fspan(1.-1,31.-1,31)
time := todouble(time)
time!0="time"
time&time = time
time@units = "days since 1996-01-01 00:00"

delete(PPT&time)
PPT&time=time


outfile=DATADIR+"prism_catted.nc"
system("/bin/rm -vf "+outfile)   ; remove any pre-existing file
ncdf = addfile(outfile ,"c")  ; open output netCDF file
print("Writing: "+outfile)

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
ncdf->PPT=PPT


end