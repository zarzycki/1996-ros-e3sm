#!/bin/bash -l

#PBS -l nodes=1:ppn=20
#PBS -l walltime=1:00:00
#PBS -A open
#PBS -j oe
#PBS -N abby_elm

module load parallel

NUMCORES=5
TIMESTAMP=`date +%s%N`
COMMANDFILE=commands.${TIMESTAMP}.txt

cd /storage/home/cmz5202/1996-ros-e3sm/ncl
     
vars=( "FSNO" "H2OSNO" "H2OSNO_TOP" "QSNOMELT" "QSNOFRZ" "SNOW" "SNOWDP" "SNOWLIQ" "SNOTTOPL" "ERRH2OSNO" "QSNWCPICE" "FGR" "FSM" "QDRIP" "ZBOT" "TG" "TV" "TSA" "SWup" "SWdown" "FIRA" "FIRE" "LWdown" "LWup" "TREFMNAV" "TREFMXAV" "Q2M" "FPSN" "FSH" "Rnet" "TLAI" "FCTR" "QOVER" "RH" "RH2M" "SOILWATER_10CM" "ZWT" "TSOI_10CM" "QSOIL" "QVEGE" "QVEGT" "QTOPSOIL" "Qstor" "H2OCAN" "H2OSFC" "QH2OSFC" "TWS" "HC" "SNOdTdzL" "SNOW_SINKS" "SNOW_SOURCES" "SNOLIQFL" "TH2OSFC" "U10" "QFLOOD" "QRUNOFF" "Qle" "Qh" "Qair" "Tair" "RAIN" "SNO_T" "SNO_Z" )

for ii in "${vars[@]}"
do
  NCLCOMMAND="ncl -n elm-diffs.ncl 
      'singlevar = \"'${ii}'\"' "
  echo ${NCLCOMMAND} >> ${COMMANDFILE}
done

parallel --jobs ${NUMCORES} --workdir $PWD < ${COMMANDFILE}

rm -v ${COMMANDFILE}