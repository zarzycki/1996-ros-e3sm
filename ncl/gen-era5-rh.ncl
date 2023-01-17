begin

f = addfile("/Volumes/Untitled/arp5873/ERA5/e5.oper.an.sfc.128_167_2t.ll025sc.1996010100_1996013123.nc","r")
g = addfile("/Volumes/Untitled/arp5873/ERA5/e5.oper.an.sfc.128_168_2d.ll025sc.1996010100_1996013123.nc","r")

T = f->VAR_2T(:,:,:)
TD = g->VAR_2D(:,:,:)

RH = relhum_ttd(T,TD,0)

copy_VarCoords(T,RH)

outname = "/Volumes/Untitled/arp5873/ERA5/e5.oper.an.sfc.128_999_rh.ll025sc.1996010100_1996013123.nc"
system("/bin/rm -f "+outname)   ; remove any pre-existing file
ncdf = addfile(outname ,"c")  ; open output netCDF file
filedimdef(ncdf,"time",-1,True) 
ncdf->RH  = RH

end