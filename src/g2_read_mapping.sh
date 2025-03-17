#!/bin/bash

#SBATCH -a 1-90
#SBATCH --mem 32G
#SBATCH -o /dev/null
#SBATCH -e /dev/null

samples=($(awk -F '\t' 'NR>1 {print $1}' data/sra_accession.tsv))
sample=${samples[$SGE_TASK_ID-1]}
reference=~/ref/grcg7b/GCF_016699485.2.fa

shopt -s expand_aliases
alias samtools="apptainer exec /usr/local/biotools/s/samtools:1.18--h50ea8bc_1 samtools"
alias bwa="apptainer exec /usr/local/biotools/b/bwa:0.7.17--h5bf99c6_8 bwa"

mkdir -p bam/${sample}

bwa mem -t 4 -M ${reference} \
  fastaq/cleaned/${sample}_clean_1.fq.gz fastaq/cleaned/${sample}_clean_2.fq.gz | \
  samtools view -b -F 4 > bam/${sample}/${sample}.bam
samtools collate -Ou bam/${sample}/${sample}.bam | \
  samtools fixmate - - -mu | \
  samtools sort - -u | \
  samtools markdup - -r bam/${sample}/${sample}.cfsm.bam
samtools index bam/${sample}/${sample}.cfsm.bam

samtools flagstat -O tsv bam/${sample}/${sample}.bam > bam/${sample}/${sample}.stat
samtools flagstat -O tsv bam/${sample}/${sample}.cfsm.bam > bam/${sample}/${sample}.cfsm.stat
samtools coverage bam/${sample}/${sample}.cfsm.bam > bam/${sample}/${sample}.cov

rm bam/${sample}/${sample}.bam
