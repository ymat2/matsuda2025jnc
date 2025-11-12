#!/bin/bash

#SBATCH -o /dev/null
#SBATCH -e /dev/null

shopt -s expand_aliases
alias plink2="apptainer exec /usr/local/biotools/p/plink2:2.00a5--h4ac6f70_0 plink2"
alias rapidnj="apptainer exec /usr/local/biotools/r/rapidnj:2.3.2--h4ac6f70_4 rapidnj"

workdir=~/popstr/phylo/rapidnj
prefix=~/popstr/structure/popstr.snp.filt
subvcf=${workdir}/popstr.sub

[ ! -e ${workdir} ] && mkdir -p ${workdir}
cd ${workdir}

plink2 --bfile ${prefix} \
  --allow-extra-chr \
  --double-id \
  --set-missing-var-ids @:# \
  --not-chr NC_052571.1, NC_052572.1, NC_053523.1 \
  --export vcf id-paste=iid\
  --out ${subvcf}

/usr/bin/python3 ~/scripts/vcf2phylip/vcf2phylip.py -i ${subvcf}.vcf --resolve-IUPAC
rm ${subvcf}.vcf*

poppy remove_invariant_site -i ${subvcf}.min4.phy -o ${subvcf}.varsites.fa --format fasta
rm ${subvcf}.min4.phy

rapidnj ${subvcf}.varsites.fa -i fa -o t -b 100 --evolution-model kim -x ${subvcf}.nj.tree
