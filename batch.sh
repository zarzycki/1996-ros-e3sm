#!/bin/bash -l

#PBS -l nodes=1:ppn=20
#PBS -l walltime=8:00:00
#PBS -A open
#PBS -j oe
#PBS -N abby.process

### Navigate to my script directly
cd /storage/home/cmz5202/scratch/abby_test/pyfiles/

### Create files
python3 create_file.py
python3 create_file_2.py

### Perform analysis
python3 Station_obs.py
python3 Sim_Comparison.py
python3 Spatial_Comp.py
python3 discharge.py
python3 Delta.py
