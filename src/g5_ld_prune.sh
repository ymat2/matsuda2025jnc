#$ -S /bin/bash
#$ -cwd
#$ -o /dev/null
#$ -e /dev/null

shopt -s expand_aliases
alias bcftools="apptainer exec /usr/local/biotools/b/bcftools:1.18--h8b25389_0 bcftools"
alias plink2="apptainer exec /usr/local/biotools/p/plink2:2.00a5--h4ac6f70_0 plink2"

workdir=~/popstr/structure
vcf=~/popstr/vcf/popstr.vcf.gz
snpvcf=~/popstr/vcf/popstr.snp.vcf.gz
prefix=popstr.snp

cd ${workdir}

bcftools view -v snps -m2 -M2 -Oz ${vcf} > ${snpvcf}
bcftools index ${snpvcf}

plink2 --vcf ${snpvcf} \
  --allow-extra-chr \
  --double-id \
  --set-missing-var-ids @:# \
  --maf 0.01 \
  --indep-pairwise 50 10 0.2 \
  --out ${prefix}

plink2 --vcf ${snpvcf} \
  --allow-extra-chr \
  --double-id \
  --set-missing-var-ids @:# \
  --extract ${prefix}.prune.in \
  --make-bed \
  --out ${prefix}

rm ${prefix}.prune.*

plink2 --bfile ${prefix} --pca --allow-extra-chr --double-id --out popstr.pca.before
