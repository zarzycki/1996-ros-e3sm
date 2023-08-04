# Pettett and Zarzycki, JHM

---

### Extract prerequisites

This directory points to the data downloaded from DataCommons.

```
/glade/u/home/zarzycki/scratch/arp5873_NEW/
```

This points to the source code pulled from Zenodo (or cloned from Github)

```
/glade/u/home/zarzycki/work/sw/1996-ros-e3sm/
```

### Install conda environment from YAML

```
cd ${SOFTPATH}
conda env create -f pettett.yml
## mamba env create -f pettett.yml     # alternatively
```

NOTE: This set of packages seems to not work with arm64 (e.g., M1 Mac). This is likely because the versions required pre-date the widespread release of M1. Other architectures like x86_64 are needed.

### Edit all .sh scripts in code directories

```
### Cheyenne
BASEPATH=/glade/u/home/zarzycki/scratch/arp5873_NEW/
SOFTPATH=/glade/u/home/zarzycki/work/sw/1996-ros-e3sm/

### ICDS
BASEPATH=/storage/home/cmz5202/group/arp5873_NEW/
SOFTPATH=/storage/home/cmz5202/work/sw/1996-ros-e3sm
```

🔴 Also, set/export those directories in the current shell environment.

🔴 This also requires setting any requisite batch job preambles at the header of each file (e.g., for sbatch, qsub, etc.)

🔴 To load the relevant conda environment from above, one may need to source a file like `.bashrc` depending on how their login shell and/or compute nodes are configured.

### Unzip relevant shapefiles

```
cd ${SOFTPATH}/shapes/
unzip srb.zip 
unzip ne_10m_rivers_lake_centerlines.zip 
```

### Generate processed datasets and perform associated analyses

```
## NCL gen + analyses
cd ${SOFTPATH}/ncl
qsubcasper driver-ens-avg.sh

## Python gen + analyses
cd ${SOFTPATH}
qsubcasper batch-python.sh
```

NOTE: Set `PROCESS_DATA` in each script to "false" if `proc_arp` and `ens_means` folders are already generated in ${BASEPATH}. Generating the intermediate files takes a few hours on Cheyenne.

NOTE: `qsubcasper driver-elm.sh` in the NCL directory generates non-required ELM difference plots.