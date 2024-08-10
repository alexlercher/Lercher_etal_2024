#!/bin/bash
# Set up path info
export PROJECTDIR=/path/to/your/git/repo
export REFDIR=${PROJECTDIR}/analysis/data/raw_data/reference # STORE YOUR GENBANK MT952602.1 GFF3 DOWNLOAD HERE

module load cufflinks/2.2.0

gffread ${REFDIR}/SARSCoV2-MA10-MT952602-1.gff3 -T -o ${REFDIR}/SARSCoV2-MA10-MT952602-1.gtf
