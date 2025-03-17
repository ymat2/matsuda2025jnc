#!/bin/bash

#SBATCH -o /dev/null
#SBATCH -e /dev/null

shopt -s expand_aliases
alias bcftools="apptainer exec /usr/local/biotools/b/bcftools:1.18--h8b25389_0 bcftools"
alias plink2="apptainer exec /usr/local/biotools/p/plink2:2.00a5--h4ac6f70_0 plink2"

workdir=~/popstr/structure
snpvcf=~/popstr/vcf/popstr.snp.vcf.gz
prefix=popstr.snp

[ ! -e ${workdir} ] && mkdir ${workdir}
cd ${workdir}

### prepare phenotype dataframe by hand for now

echo -e "#FID\tIID\tSEQ" > sequencer.pheno
bcftools query -l ${snpvcf} >> sequencer.pheno
# edit

### Run GWAS using sequencer as phenotype

plink2 --bfile ${prefix} \
  --allow-extra-chr \
  --out sequencer.gwas \
  --pheno sequencer.pheno \
  --pheno-name SEQ \
  --glm allow-no-covars skip-invalid-pheno firth-fallback hide-covar omit-ref \
  --adjust \
  --ci 0.95

cat sequencer.gwas.SEQ.glm.logistic.hybrid.adjusted | awk -F '\t' '{print $2 "\t" $10}' > gwas-result.txt
rm sequencer.gwas.SEQ.glm.logistic.hybrid*

### Compute Fst between sequencer types

plink2 --bfile ${prefix} \
  --allow-extra-chr \
  --out sequencer \
  --pheno sequencer.pheno \
  --fst SEQ method=wc report-variants

cat sequencer.CASE.CONTROL.fst.var | awk -F '\t' '{print $1 "\t" $2 "\t" $5}' > fst-result.txt
rm sequencer.CASE.CONTROL.fst.var

### Extract sites with higher pvalue using --extract

cat gwas-result.txt | awk -F '\t' 'NR>1 { if ($2 > 0.05) { print $1 }}' > sites.txt

plink2 --bfile ${prefix} \
  --allow-extra-chr \
  --double-id \
  --set-missing-var-ids @:# \
  --extract sites.txt \
  --make-bed \
  --out ${prefix}.filt

### PCA after

plink2 --bfile ${prefix}.filt --pca --allow-extra-chr --double-id --out popstr.pca.after

sed -i -e 's/NW_//g' ${prefix}.filt.bim
sed -i -e 's/NC_//g' ${prefix}.filt.bim
sed -i -e 's/\.1//g' ${prefix}.filt.bim
