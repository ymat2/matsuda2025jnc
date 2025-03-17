#!/bin/bash

#SBATCH -a 1-90
#SBATCH -o /dev/null
#SBATCH -e /dev/null

sras=($(awk -F '\t' 'NR>1 {print $1}' data/sra_accession.tsv))
sra=${sras[$SGE_TASK_ID-1]}

workdir=~/popstr/sra
uncleaned=~/popstr/fastaq/uncleaned/

[ ! -e ${workdir} ] && mkdir -p ${workdir}
[ ! -e ${uncleaned} ] && mkdir -p ${uncleaned}

cd ${workdir}
prefetch --max-size 30G ${sra}
fasterq-dump ${sra}/${sra}.sra --outdir ${uncleaned}/${sra}
rm -r ${sra}
