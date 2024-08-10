#!/bin/bash
# Set up path info
export PROJECTDIR=/path/to/your/git/repo
export REFDIR=${PROJECTDIR}/analysis/data/raw_data/reference
export SEQDIR=${PROJECTDIR}/analysis/data/raw_data/sequencing_data
export REFDIR=${REFDIR}/mm10-PR8-MA10/mm10-PR8-MA10

# Paths for FASTQ files
export RUN1=${SEQDIR}/H2VJVBGXJ/outs/fastq_path # STORE FASTQ HERE

cd ${SEQDIR}

# CellRanger count must be run individually for each 10X lane
module load cellranger/5.0.0

# BALF1
cellranger count --id BALF1_10X \
	--sample=BALF1 \
  --fastqs=${RUN1} \
  --transcriptome=${REFDIR} \
  --expect-cells=5000 \
  --nosecondary \
  --maxjobs=200

# BALF2
cellranger count --id BALF2_10X \
	--sample=BALF2 \
  --fastqs=${RUN1} \
  --transcriptome=${REFDIR} \
  --expect-cells=5000 \
  --nosecondary \
  --maxjobs=200
