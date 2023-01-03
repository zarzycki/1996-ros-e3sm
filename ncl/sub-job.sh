#!/bin/bash -l

#PBS -l nodes=1:ppn=20
#PBS -l walltime=12:00:00
#PBS -A open
#PBS -j oe
#PBS -N abby_ncl

module load parallel

NUMCORES=5
TIMESTAMP=`date +%s%N`
COMMANDFILE=commands.${TIMESTAMP}.txt

cd /storage/home/cmz5202/1996-ros-e3sm/ncl
     
harr=( "h0" "h1" "h2" "h3" )
models=( "mosart" "eam" "elm" )
exps=( "PI" "control" "plus1K" "plus2K" "plus3K" "plus4K" )
for ii in "${harr[@]}"
do
  for jj in "${models[@]}"
  do
    for kk in "${exps[@]}"
    do
      NCLCOMMAND="ncl -n create-ens-avg.ncl 
          'hname = \"'${ii}'\"' 
          'experiment = \"'${kk}'\"' 
          'modelcomponent = \"'${jj}'\"' "
      echo ${NCLCOMMAND} >> ${COMMANDFILE}
    done
  done 
done

parallel --jobs ${NUMCORES} --workdir $PWD < ${COMMANDFILE}

rm -v ${COMMANDFILE}