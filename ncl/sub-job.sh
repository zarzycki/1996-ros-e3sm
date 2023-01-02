#!/bin/bash -l

#PBS -l nodes=1:ppn=20
#PBS -l walltime=2:00:00
#PBS -A open
#PBS -j oe
#PBS -N abby_ncl

module load parallel

# use numcores = 32 for 1deg runs, 16 for high-res runs
NUMCORES=5
TIMESTAMP=`date +%s%N`
COMMANDFILE=commands.${TIMESTAMP}.txt

cd /storage/home/cmz5202/arp-ncl
     
harr=( "h0" "h1" "h2" "h3" )
models=( "mosart" "eam" "elm" )
for ii in "${harr[@]}"
do
  for jj in "${models[@]}"
  do
    NCLCOMMAND="ncl -n create-ens-avg.ncl 
        'hname = \"'${ii}'\"' 
        'modelcomponent = \"'${jj}'\"' "
    echo ${NCLCOMMAND} >> ${COMMANDFILE}
  done 
done

parallel --jobs ${NUMCORES} --workdir $PWD < ${COMMANDFILE}

rm -v ${COMMANDFILE}