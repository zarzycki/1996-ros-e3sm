#!/bin/bash

#PBS -l nodes=1:ppn=20
#PBS -l walltime=8:00:00
#PBS -A open
#PBS -j oe
#PBS -N process

set -e

### Perform analysis
#BASEPATH=/Volumes/ZARZYCKI_FLASH/arp5873/
#SOFTPATH=/Users/cmz5202/Software/1996-ros-e3sm/
BASEPATH=/storage/home/cmz5202/group/arp5873_NEW/
SOFTPATH=/storage/home/cmz5202/scratch/1996-ros-e3sm/

### mamba create -n pettett -c conda-forge python=3.9.7 metpy=1.0 netCDF4=1.5.8 numpy=1.20.1 xarray=0.16.2 xesmf=0.5.2 rasterio=1.2.10 cartopy=0.20.1 matplotlib=3.3.4 pandas=1.2.2 scipy=1.7.1 pint=0.16.1 gdal=3.3.3 poppler=21.09.0 geopandas=0.10.2
### conda activate pettett
### mamba env create -f pettett2.yml

### Navigate to my script directly
cd $SOFTPATH

### Activate conda env
conda activate pettett

### Create files
python3 create_file.py $BASEPATH
python3 create_file_2.py $BASEPATH

### Do analysis / create figs
python3 Station_obs.py $BASEPATH
python3 Sim_Comparison.py $BASEPATH
python3 Spatial_Comp.py $BASEPATH $SOFTPATH
python3 discharge.py $BASEPATH $SOFTPATH
python3 Delta.py
