begin

dirs = (/"1996011512","1996011600","1996011612","1996011700","1996011712"/)
ensembles = (/"004","005","006"/)
;experiment = "PI"
;modelcomponent = "eam"
;hname = "h2"
print_holder_var=False
grid = "ne0conus30x8"   ; ne30 or ne0conus30x8
restrict_domain=False
;minlat=36.0
;maxlat=46.0
;minlon=269.0-360.
;maxlon=294.0-360.
minlat=27.5
maxlat=55.0
minlon=-105.0
maxlon=-65.0

print(modelcomponent+" "+hname)

if (modelcomponent .eq. "eam" .and. hname .eq. "h0")
  vars = (/"U10","PRECT","PRECSL","PRECSC","U850","V850","PSL","T850","U500","V500","T500","T300","T200","Z300","Z200","FLUT","SST","U200","V200","OMEGA500","OMEGA850","TMQ","Z500","TUQ","TVQ","TUH","TVH","CAPE","TMCLDLIQ","TMCLDICE","TMNUMLIQ","TMNUMICE","TMRAINQM","TMSNOWQM","TMNUMRAI","TMNUMSNO","WSPDSRFMX","WSPDSRFAV","TREFHT","TREFHTMN","TREFHTMX","QREFHT","U10","RHREFHT","SHFLX","LHFLX","QFLX","SNOWHLND","SNOWHICE","PCLDTOP","PCLDBOT","PBLH","TGCLDCWP","TGCLDLWP","TGCLDIWP","CLDTOT","CLDLOW","CLDMED","CLDHGH"/)
  arr_2D = new(dimsizes(vars),"logical")
  arr_2D = True
end if
if (modelcomponent .eq. "eam" .and. hname .eq. "h1")
  vars = (/"PRECT","PRECC","PRECL","PRECSL","PRECSC"/)
  arr_2D = (/True,True,True,True,True/)
end if
if (modelcomponent .eq. "eam" .and. hname .eq. "h2")
  vars = (/"U","V","T","Q","Z3","PS"/)
  arr_2D = (/False,False,False,False,False,True/)
end if
if (modelcomponent .eq. "eam" .and. hname .eq. "h3")
  vars = (/"PS","TBOT","UBOT","VBOT","QBOT","FLDS","FSDS","PRECT","QREFHT","RHREFHT","TREFHT","TREFHTMN","TREFHTMX","U10","WSPDSRFMX","WSPDSRFAV"/)
  arr_2D = new(dimsizes(vars),"logical")
  arr_2D = True
end if
if (modelcomponent .eq. "elm" .and. hname .eq. "h0")
  vars = (/"FSNO","H2OSNO","H2OSNO_TOP","QSNOMELT","QSNOFRZ","SNOW","SNOWDP","SNOWLIQ","SNOTTOPL","ERRH2OSNO","QSNWCPICE","FGR","FSM","QDRIP","ZBOT","TG","TV","TSA","SWup","SWdown","FIRA","FIRE","LWdown","LWup","TREFMNAV","TREFMXAV","Q2M","FPSN","FSH","Rnet","TLAI","FCTR","QOVER","RH","RH2M","SOILWATER_10CM","ZWT","TSOI_10CM","QSOIL","QVEGE","QVEGT","QTOPSOIL","Qstor","H2OCAN","H2OSFC","QH2OSFC","TWS","HC","SNOdTdzL","SNOW_SINKS","SNOW_SOURCES","SNOLIQFL","TH2OSFC","U10","QFLOOD","QRUNOFF","Qle","Qh","Qair","Tair","RAIN","SNO_T","SNO_Z"/)
  arr_2D = new(dimsizes(vars),"logical")
  arr_2D = True
  ; cheat and set last two vars to 3D since they are snow T and Z of layers
  nvarstmp = dimsizes(vars)
  arr_2D(nvarstmp-1)=False
  arr_2D(nvarstmp-2)=False
  delete(nvarstmp)
end if
if (modelcomponent .eq. "mosart" .and. hname .eq. "h0")
  vars = (/"RIVER_DISCHARGE_OVER_LAND_LIQ","RIVER_DISCHARGE_TO_OCEAN_LIQ"/)
  arr_2D = new(dimsizes(vars),"logical")
  arr_2D = True
end if

if (.not. isvar("vars")) then
  print("uh oh!")
  exit
end if

ndirs = dimsizes(dirs)
nens = dimsizes(ensembles)
nvars = dimsizes(vars)



do zz = 0,nvars-1

varname = vars(zz)
is_2D = arr_2D(zz)
print("is_2D?: "+is_2D)

do kk = 0,nens-1
do ii = 0,ndirs-1

   thisMember = (kk*ndirs) + ii
   print(""+thisMember)

   thisdir = dirs(ii)
   thisensemble = ensembles(kk)

   filedir = DATADIR+"/E3SM/RoS-F2010C5-"+grid+"-"+thisensemble+"-"+experiment+"/"+thisdir+"/*."+modelcomponent+"."+hname+"*.nc_regrid.nc"
   print("loading "+filedir)

   fils = systemfunc ("ls "+filedir)
   ;print(fils)

   f    = addfiles (fils, "r")
   ListSetType (f, "cat")
   if (restrict_domain)
     if (is_2D) then
       var = f[:]->$varname$(:,{minlat:maxlat},{minlon:maxlon})
     else
       var = f[:]->$varname$(:,:,{minlat:maxlat},{minlon:maxlon})
     end if
   else
     if (is_2D) then
       var = f[:]->$varname$(:,:,:)
     else
       var = f[:]->$varname$(:,:,:,:)
     end if
   end if

   if (ii .eq. 0 .and. kk .eq. 0)
     lat = var&lat
     lon = var&lon
     time = var&time
     nlat = dimsizes(lat)
     nlon = dimsizes(lon)
     ntim = dimsizes(time)

     if (.not. is_2D) then
       lev = var&lev
       nlev = dimsizes(lev)
     end if

     if (.not. is_2D) then
       holderVar = new((/ndirs*nens,ntim,nlev,nlat,nlon/),typeof(var))
       holderVar!0 = "ens"
       holderVar!1 = "time"
       holderVar!2 = "lev"
       holderVar!3 = "lat"
       holderVar!4 = "lon"
       holderVar&time = time
       holderVar&lev = lev
       holderVar&lat = lat
       holderVar&lon = lon
     else
       holderVar = new((/ndirs*nens,ntim,nlat,nlon/),typeof(var))
       holderVar!0 = "ens"
       holderVar!1 = "time"
       holderVar!2 = "lat"
       holderVar!3 = "lon"
       holderVar&time = time
       holderVar&lat = lat
       holderVar&lon = lon
     end if

     copy_VarAtts(var, holderVar)

     count_adds = new(ntim,"integer")
     count_adds!0 = "time"
     count_adds&time=time
     ; set to zero
     count_adds = 0

   end if

   thistime = var&time
   thistime_newunits = cd_convert(thistime,time@units)

   do jj = 0,ntim-1
     test = thistime_newunits - time(jj)
     if (min(test) .lt. 0.0001) then
       jx = closest_val(time(jj),thistime_newunits)
       ;print(jj+" "+jx)
       if (.not. is_2D) then
         holderVar(thisMember,jj,:,:,:) = (/var(jx,:,:,:)/)
       else
         holderVar(thisMember,jj,:,:) = (/var(jx,:,:)/)
       end if
       count_adds(jj) = count_adds(jj) + 1
     else
       ;print("nope!")
     end if
   end do

   delete(fils)
   if (kk .ne. (nens-1) .or. ii .ne. (ndirs-1))
     delete(f)
   end if
   delete(var)
   delete(thistime)
   delete(thistime_newunits)
   delete(filedir)
end do
end do

; calculate ensemble average
ens_avg = dim_avg_n_Wrap(holderVar,0)

; define name of individual members
varname_member = varname+"_indiv"
count_name = "count_"+varname

setfileoption("nc","Format","LargeFile")

outputdir=DATADIR+"/ens_means/"
system("mkdir -p "+outputdir)
outputname="output_"+experiment+"_"+modelcomponent+"_"+hname+".nc"
filename=outputdir+"/"+outputname
if (zz .eq. 0) then
  system("/bin/rm -f "+filename)   ; remove any pre-existing file
  ncdf = addfile(filename ,"c")  ; open output netCDF file
  fAtt = True            ; assign file attributes
  fAtt@creation_date = systemfunc ("date")
  fileattdef( ncdf, fAtt )            ; copy file attributes
  ;filedimdef(ncdf,"time",-1,True)
else
  ncdf = addfile(filename ,"w")  ; open output netCDF file
end if

print("***** WRITING VARIABLE:  "+varname)
print("***** FILENAME: "+filename)
print("*****")

if (print_holder_var) then
  ncdf->$varname_member$  = holderVar
end if
ncdf->$varname$  = ens_avg
if (zz .eq. 0) then
  ncdf->count =  count_adds
end if

; if this is our first file and we have a 3D var, add these to the ensemble file
if (zz .eq. 0 .and. any( arr_2D .eq. False) .and. modelcomponent .eq. "eam" ) then
  ncdf->hyai = f[0]->hyai
  ncdf->hybi = f[0]->hybi
  ncdf->hyam = f[0]->hyam
  ncdf->hybm = f[0]->hybm
  ncdf->P0 = f[0]->P0
end if

delete(lat)
delete(lon)
delete(time)
delete(ens_avg)
delete(holderVar)
delete(count_adds)
delete(f)
delete(ncdf)

list_vars()

end do ; zz

exit

end
