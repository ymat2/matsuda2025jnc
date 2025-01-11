#$ -S /bin/bash
#$ -cwd
#$ -t 1-15:1
#$ -l mem_req=63G
#$ -l s_vmem=63G
#$ -o /dev/null
#$ -e /dev/null


shopt -s expand_aliases
alias admixture="apptainer exec /usr/local/biotools/a/admixture:1.3.0--0 admixture"

workdir=~/popstr/admixture
bed=~/popstr/structure/popstr.snp.filt.bed

[ ! -e ${workdir} ] && mkdir ${workdir}
cd ${workdir}

admixture --cv ${bed} ${SGE_TASK_ID} | tee log${SGE_TASK_ID}.out

grep -h CV log*.out > CV-error.txt
