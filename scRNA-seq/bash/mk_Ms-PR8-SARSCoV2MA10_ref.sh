#!/bin/bash
# Set up paths
export PROJECTDIR=/path/to/your/git/repo
export REFDIR=${PROJECTDIR}/analysis/data/raw_data/reference # STORE YOUR IAV AND MA10 TRANSCRIPTOMES HERE
export CRDIR=/path/to/your/cellranger/5.0.0/refdata-cellranger/refdata-gex-mm10-2020-A
export LOGDIR=${PROJECTDIR}/logs

# Mouse reference paths
## Original source: Cellranger v5.0.0 refdata
export MM10GTF=${CRDIR}/genes/genes.gtf
export MM10FASTA=${CRDIR}/fasta/genome.fa
export MM10OUT=${REFDIR}/mm10

# SARS-CoV-2/MA10 reference paths
## Original source: NCBI GenBank MT952602.1
export MA10GTF=${REFDIR}/SARSCoV2-MA10-MT952602-1.gtf
export MA10FASTA=${REFDIR}/SARSCoV2-MA10-MT952602-1.fasta
export MA10OUT=${REFDIR}/

# IAV A/Puerto Rico/8/1934 H1N1 reference paths
## Original source: NCBI Reference Sequence NC_002023
export PR8GTF=${REFDIR}/IAV_PR8.gtf
export PR8FASTA=${REFDIR}/IAV_PR8.fasta
export PR8OUT=${REFDIR}

# Make output directory
OUTDIR=${REFDIR}/mm10-PR8-MA10
if [[ ! -d $OUTDIR ]]; then
  mkdir $OUTDIR
fi

cd $OUTDIR

# Make joint mm10/MA10/PR8 fasta and gtf files
touch ${OUTDIR}/mm10-PR8-MA10.fa
cat $MM10FASTA > ${OUTDIR}/mm10-PR8-MA10.fa
cat $MA10FASTA >> ${OUTDIR}/mm10-PR8-MA10.fa
cat $PR8FASTA >> ${OUTDIR}/mm10-PR8-MA10.fa

touch ${OUTDIR}/mm10-PR8-MA10.gtf
cat $MM10GTF > ${OUTDIR}/mm10-PR8-MA10.gtf
cat $MA10GTF >> ${OUTDIR}/mm10-PR8-MA10.gtf
cat $PR8GTF >> ${OUTDIR}/mm10-PR8-MA10.gtf

# Filter GTF, then make joint reference
## From Cellranger/Software/Pipelines/Advanced/References
module load cellranger/5.0.0
cellranger mkgtf ${OUTDIR}/mm10-PR8-MA10.gtf ${OUTDIR}/mm10-PR8-MA10-filtered.gtf \
  --attribute=gene_biotype:protein_coding \
	--attribute=gene_biotype:lincRNA \
	--attribute=gene_biotype:antisense; \
cellranger mkref --ref-version=2021-11-30 \
  --genome=mm10-PR8-MA10 \
	--fasta=${OUTDIR}/mm10-PR8-MA10.fa \
	--genes=${OUTDIR}/mm10-PR8-MA10-filtered.gtf \
	--nthreads=8
