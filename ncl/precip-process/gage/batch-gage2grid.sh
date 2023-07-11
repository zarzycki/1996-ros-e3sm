#!/bin/bash -l

DATADIR="/Users/cmz5202/NetCDF/arp5873_NEW/obs_gage/"

# $(seq -w 1 31)
for i in $(seq -w 17 23); do
  ncl gage-to-grid.ncl 'dd_str="'$i'"' 'DATADIR="'$DATADIR'"'
done

cd $DATADIR
ncrcat ds_simple_*.nc gage_ds_catted.nc
ncks -4 -L 1 gage_ds_catted.nc gage_ds_catted.nc
rm -v ds_simple_*.nc

