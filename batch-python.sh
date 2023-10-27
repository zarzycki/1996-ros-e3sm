#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --ntasks=12
#SBATCH --mem=100GB
#SBATCH --time=12:00:00
#SBATCH --partition=open

##PBS -N arp-1996-analysis
##PBS -A P93300642
##PBS -l select=1:ncpus=6:mem=200GB
##PBS -l walltime=8:00:00
##PBS -q casper
##PBS -j oe

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

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

PROCESS_DATA="false"   # true or false
BASEPATH=/storage/group/cmz5202/default/arp5873/
SOFTPATH=/scratch/cmz5202/1996-ros-e3sm

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

set -e

date >> timing.txt

source ~/.bashrc

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

date >> timing.txt
echo "-------" >> timing.txt