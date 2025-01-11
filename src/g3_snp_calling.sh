#$ -S /bin/bash
#$ -cwd
#$ -t 1-89
#$ -l medium
#$ -l s_vmem=64G
#$ -l mem_req=64G
#$ -o /dev/null
#$ -e /dev/null


samples=($(ls bam | sort -V))
sample=${samples[$SGE_TASK_ID-1]}
reference=~/ref/grcg7b/GCF_016699485.2.fa

shopt -s expand_aliases
alias bcftools="apptainer exec /usr/local/biotools/b/bcftools:1.18--h8b25389_0 bcftools"

bcftools mpileup -f ${reference} bam/${sample}/${sample}.cfsm.bam \
  --max-depth 1000 -a FORMAT/AD,FORMAT/DP --no-BAQ | \
  bcftools call -vm -Oz -o bam/${sample}/${sample}.vcf.gz
bcftools index bam/${sample}/${sample}.vcf.gz
bcftools stats bam/${sample}/${sample}.vcf.gz > bam/${sample}/${sample}.vcf.stat

mode_depth=$(cat bam/${sample}/${sample}.vcf.stat | awk -F '\t' '$1 == "DP"' | sort -t $'\t' -k7,7nr | awk -F '\t' 'NR==1 {print $3}')
max_depth=$((${mode_depth} * 2))
min_depth=$((${mode_depth} / 2))

bcftools filter bam/${sample}/${sample}.vcf.gz \
  -i "QUAL>30 && INFO/DP>=${min_depth} && INFO/DP<=${max_depth}" \
  -Oz -o bam/${sample}/${sample}.q.vcf.gz
bcftools index bam/${sample}/${sample}.q.vcf.gz
bcftools stats bam/${sample}/${sample}.q.vcf.gz > bam/${sample}/${sample}.q.vcf.stat
