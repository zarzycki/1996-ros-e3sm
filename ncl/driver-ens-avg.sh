#!/bin/bash -l

##### ICDS ROAR
##PBS -l nodes=1:ppn=6:himem
##PBS -l pmem=18gb
##PBS -l walltime=11:00:00
##PBS -A open
##PBS -j oe
##PBS -N abby_ncl

##### Cheyenne
#PBS -N arp-ens-avg
#PBS -A P93300642
#PBS -l select=1:ncpus=6:mem=100GB
#PBS -l walltime=10:00:00
#PBS -q casper
#PBS -j oe

source /glade/u/home/zarzycki/.bashrc

NUMCORES=6
PROCESS_DATA="false"   # true or false
BASEPATH=/glade/u/home/zarzycki/scratch/arp5873_NEW/
SOFTPATH=/glade/u/home/zarzycki/work/sw/1996-ros-e3sm/

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

  for kk in "${exps[@]}"
  do
    for jj in "${models[@]}"
    do
      for ii in "${harr[@]}"
      do
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