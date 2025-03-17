#!/bin/bash

#SBATCH -o /dev/null
#SBATCH -e /dev/null

shopt -s expand_aliases
alias bcftools="apptainer exec /usr/local/biotools/b/bcftools:1.18--h8b25389_0 bcftools"

poppy summary -i bam -o out/flag_summary.tsv --type flag
poppy summary -i bam -o out/mapping_summary.tsv --type bam
poppy summary -i bam -o out/vcf_summary.tsv --type vcf

samples_sea=($(ls ~/popstr/bam | sort -V))  # Samples of South-east Asia
samples_jnc=($(ls ~/jpool/bam | sort -V))   # Samples of Japan
samples_abrc=($(ls ~/ABRC/bam | sort -V))   # Samples of ABRC

rm out/use_samples.txt
for s in ${samples_sea[@]}; do echo ~/popstr/bam/$s/$s.q.vcf.gz >> out/use_samples.txt; done
for s in ${samples_jnc[@]}; do echo ~/jpool/bam/$s/$s.q.vcf.gz >> out/use_samples.txt; done

# edit out/use_samples.txt
files=($(cat out/use_samples.txt | grep -v '^#'))

workdir=~/popstr/vcf
vcf=${workdir}/popstr.vcf.gz

[ ! -e ${workdir} ] && mkdir ${workdir}
cd ${workdir}

bcftools merge "${files[@]}" -0 -Oz -o tmp.vcf.gz
bcftools query -l tmp.vcf.gz | awk -F'/' '{print $0 "\t" $2}' > rename.txt
bcftools reheader --samples rename.txt tmp.vcf.gz -o ${vcf}
bcftools index ${vcf}

rm rename.txt tmp.vcf.gz

echo -e "sample\tn_hom\tn_het\tn_indel" > out/snp_count.txt
for s in $(ls bam); do
  n_hom=$(bcftools view --no-header --types snps --genotype hom bam/${s}/${s}.q.vcf.gz | wc -l)
  n_het=$(bcftools view --no-header --types snps --genotype het bam/${s}/${s}.q.vcf.gz | wc -l)
  n_indel=$(bcftools view --no-header --types indels bam/${s}/${s}.q.vcf.gz | wc -l)
  echo -e "${s}\t${n_hom}\t${n_het}\t${n_indel}" >> out/snp_count.txt
done
