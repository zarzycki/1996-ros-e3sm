#!/bin/bash -l

##### ICDS ROAR COLLAB
#SBATCH --nodes=1
#SBATCH --ntasks=6
#SBATCH --mem=18GB
#SBATCH --time=11:00:00
#SBATCH --partition=open

##### Cheyenne
##PBS -N arp-ens-avg
##PBS -A P93300642
##PBS -l select=1:ncpus=6:mem=100GB
##PBS -l walltime=10:00:00
##PBS -q casper
##PBS -j oe

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

NUMCORES=6                     # cores for GNU parallel
PROCESS_DATA="false"           # true or false
BASEPATH=/storage/group/cmz5202/default/arp5873/
SOFTPATH=/scratch/cmz5202/1996-ros-e3sm

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Timing file
date >> timing.txt

# This code gets the conda environment activated on the compute nodes
# May have be changed for different systems
source ~/.bashrc
conda activate pettett

if [ "$PROCESS_DATA" = "true" ]; then

  echo "Generating ens_mean data"

  source `which env_parallel.bash`

  TIMESTAMP=`date +%s%N`
  COMMANDFILE=commands.${TIMESTAMP}.txt

  cd ${SOFTPATH}/ncl

  harr=( "h0" "h1" "h2" "h3" )
  models=( "mosart" "eam" "elm" )
  exps=( "PI" "control" "plus1K" "plus2K" "plus3K" "plus4K" )

  for kk in "${exps[@]}" ; do
    for jj in "${models[@]}" ; do
      for ii in "${harr[@]}" ; do
        NCLCOMMAND="ncl -n create-ens-avg.ncl
            'DATADIR = \"'${BASEPATH}'\"'
            'hname = \"'${ii}'\"'
            'experiment = \"'${kk}'\"'
            'modelcomponent = \"'${jj}'\"' "
        echo ${NCLCOMMAND} >> ${COMMANDFILE}
      done
    done
  done

  parallel --jobs ${NUMCORES} --workdir $PWD < ${COMMANDFILE}

  rm -v ${COMMANDFILE}

fi

#### Now perform NCL analysis

cd ${SOFTPATH}/ncl
echo "Running NCL analysis code"
ncl plot-z.ncl 'DATADIR = "'${BASEPATH}'"'
ncl panel_precip_hov.ncl 'DATADIR = "'${BASEPATH}'"'
ncl elm-budget.ncl 'DATADIR = "'${BASEPATH}'"' 'SOFTDIR = "'${SOFTPATH}'"'
ncl elm-diffs.ncl 'singlevar="TSOI_10CM"' 'DATADIR = "'${BASEPATH}'"' 'outdir = "./"'
ncl elm-diffs.ncl 'singlevar="FSNO"' 'DATADIR = "'${BASEPATH}'"' 'outdir = "./"'

date >> timing.txt
echo "-------" >> timing.txt