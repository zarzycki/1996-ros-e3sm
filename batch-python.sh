#!/bin/bash -l

#PBS -l nodes=1:ppn=20
#PBS -l walltime=11:58:00
#PBS -A open
#PBS -j oe
#PBS -N arp-py-process

##PBS -N arp-1996-analysis
##PBS -A P93300642
##PBS -l select=1:ncpus=6:mem=200GB
##PBS -l walltime=8:00:00
##PBS -q casper
##PBS -j oe

source ~/.bashrc

retry_script() {
    local script_name=$1
    local success=0
    set +e
    for i in {1..3}; do
        python3 $script_name
        if [[ $? -eq 0 ]]; then
            echo "----- $script_name executed successfully!"
            success=1
            break
        else
            echo "----- $script_name failed. Attempt $i/3."
        fi
        sleep 5
    done
    set -e  # Re-enable exit-on-error
    if [[ $success -eq 1 ]]; then
        echo "----- Proceeding with further commands..."
        return 0
    else
        echo "----- $script_name failed after 5 attempts."
        return 1
    fi
}

set -e

PROCESS_DATA="false"   # true or false
BASEPATH=/storage/home/cmz5202/group/arp5873_NEW/
SOFTPATH=/storage/home/cmz5202/work/sw/1996-ros-e3sm

### mamba create -n pettett -c conda-forge python=3.9.7 metpy=1.0 netCDF4=1.5.8 numpy=1.20.1 xarray=0.16.2 xesmf=0.5.2 rasterio=1.2.10 cartopy=0.20.1 matplotlib=3.3.4 pandas=1.2.2 scipy=1.7.1 pint=0.16.1 gdal=3.3.3 poppler=21.09.0 geopandas=0.10.2
### conda activate pettett
### mamba env create -f pettett2.yml

### Navigate to my script directly
cd $SOFTPATH

### Activate conda env
conda activate pettett

### Create files
if [ "$PROCESS_DATA" = "true" ]; then
  echo "Processing data"
  python3 create_file.py $BASEPATH
  python3 create_file_2.py $BASEPATH
fi

### Do analysis / create figs
echo "Running Python analysis"
retry_script "Station_obs.py $BASEPATH"
retry_script "Sim_Comparison.py $BASEPATH"
retry_script "Spatial_Comp.py $BASEPATH $SOFTPATH"
retry_script "discharge.py $BASEPATH $SOFTPATH"
retry_script "Delta.py"
