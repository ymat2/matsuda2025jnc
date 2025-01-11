#$ -S /bin/bash
#$ -cwd
#$ -V
#$ -l s_vmem=16G
#$ -l mem_req=16G
#$ -o /dev/null
#$ -e /dev/null

shopt -s expand_aliases
alias bcftools="apptainer exec /usr/local/biotools/b/bcftools:1.18--h8b25389_0 bcftools"

files=()
for s in $(ls ~/popstr/bam); do files+=("/home/ymat2/popstr/bam/$s/$s.mt.vcf.gz"); done
for s in $(ls ~/jpool/bam); do files+=("/home/ymat2/jpool/bam/$s/$s.mt.vcf.gz"); done

workdir=~/popstr/vcf
vcf=${workdir}/popstr.mt.vcf.gz

[ ! -e ${workdir} ] && mkdir ${workdir}

bcftools merge "${files[@]}" -0 -Ov |\
  bcftools view --exclude-type indels -Ov |\
  bcftools norm --rm-dup all -Oz -o ${vcf}
bcftools index ${vcf}
