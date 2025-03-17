#!/bin/bash

#SBATCH --mem 32G
#SBATCH -o /dev/null
#SBATCH -e /dev/null


samples=($(ls bam | sort -V))
reference=~/ref/grcg7b/GCF_016699485.2.fa

shopt -s expand_aliases
alias bcftools="apptainer exec /usr/local/biotools/b/bcftools:1.18--h8b25389_0 bcftools"

for sample in ${samples[@]}; do
  bcftools mpileup -f ${reference} bam/${sample}/${sample}.cfsm.bam \
    --max-depth 1000 -a FORMAT/AD --no-BAQ --regions NC_053523.1 | \
    bcftools call -c -Oz -o bam/${sample}/${sample}.mt.vcf.gz
  bcftools index bam/${sample}/${sample}.mt.vcf.gz
done
