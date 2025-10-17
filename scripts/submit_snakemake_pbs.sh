#!/bin/bash

##################
# Inspiration: https://github.com/SchlossLab/snakemake_cluster_tutorial
####  PBS settings

####  Job name
#PBS -N submit_snakemake

####  Resources

#PBS -l nodes=1:ppn=1,mem=16000mb
#PBS -l walltime=72:00:00

####  Account and return

#PBS -M michiel.perneel@ugent.be

#PBS -o /data/gent/vo/001/gvo00125/vsc43619/spring_campaign_2023/

#### Load Snakemake conda environment
source activate snakemake_7
source /data/gent/vo/001/gvo00125/vsc43619/spring_campaign_2023/config.sh

##################

##  Change to the directory from which you submit the job, if running
##  from within a job
if [ -d "$PBS_O_WORKDIR" ] ; then
    cd $PBS_O_WORKDIR
fi

# Initiating snakemake and running workflow in cluster mode
snakemake --profile hpc_config/ --latency-wait 60 --use-conda --rerun-incomplete \
    --conda-cleanup-pkgs --conda-prefix ${conda_env} --conda-frontend conda \
    --max-jobs-per-second 5