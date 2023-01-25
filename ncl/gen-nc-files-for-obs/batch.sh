#!/bin/bash

#export BASEPATH=/Volumes/ZARZYCKI_FLASH/arp5873/
#export SOFTPATH=/Users/cmz5202/Software/1996-ros-e3sm/
#export BASEPATH=/storage/home/cmz5202/group/arp5873_NEW/
#export SOFTPATH=/storage/home/cmz5202/scratch/1996-ros-e3sm/
export BASEPATH=/glade/u/home/zarzycki/scratch/arp5873_NEW/
export SOFTPATH=/glade/u/home/zarzycki/work/sw/1996-ros-e3sm/

ncl ncl_model-precip.ncl
ncl ncl_model-tmax.ncl 'minmax="min"'
ncl ncl_model-tmax.ncl 'minmax="max"'
ncl ncl_obs-precip.ncl
ncl ncl_obs-tmax.ncl 'var="tmax"'
ncl ncl_obs-tmax.ncl 'var="tmin"'
ncl ncl_obs-swe.ncl 'data="obs"'
ncl ncl_obs-swe.ncl 'data="obs_snowdepth"'
ncl ncl_obs-swe.ncl 'data="NLDAS"'
ncl ncl_obs-swe.ncl 'data="JRA"'
ncl ncl_obs-swe.ncl 'data="ERA5"'
