#!/bin/bash -l

#PBS -l nodes=1:ppn=20
#PBS -l walltime=8:00:00
#PBS -A open
#PBS -j oe
#PBS -N process

### Navigate to my script directly
cd /storage/home/cmz5202/1996-ros-e3sm/
### mamba create -n pettett -c conda-forge python=3.9.7 metpy=1.0 netCDF4=1.5.8 numpy=1.20.1 xarray=0.16.2 xesmf=0.5.2 rasterio=1.2.10 cartopy=0.20.1 matplotlib=3.3.4 pandas=1.2.2 scipy=1.7.1 pint=0.16.1 gdal=3.3.3 poppler=21.09.0 geopandas=0.10.2
### conda activate pettett

### Create files
#python create_file.py
#python create_file_2.py

### Perform analysis
python Station_obs.py
python Sim_Comparison.py
python Spatial_Comp.py
python discharge.py

python Delta.py