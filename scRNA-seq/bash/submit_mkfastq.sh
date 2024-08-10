#!/bin/bash
# Set up path info
export PROJECTDIR=/path/to/your/git/repo
export SEQDIR=${PROJECTDIR}/analysis/data/raw_data/sequencing_data # STORE YOUR RUN RAW DATA HERE

cd ${SEQDIR}

# Make log directory
mkdir ${PROJECTDIR}/logs
export LOGDIR=${PROJECTDIR}/logs

# Names of sequencing runs
export RUN1=211123_NB502110_0378_AH2VJVBGXJ

# Sample sheet
export SAMPLESHEET1=${PROJECTDIR}/bash/samplesheet.csv

## Process BCL sequence files to fastq files
module load cellranger/5.0.0
cellranger mkfastq --run ${SEQDIR}/${RUN1} --csv=${SAMPLESHEET1} --localcores=8 --localmem=60
