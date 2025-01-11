#$ -S /bin/bash
#$ -cwd
#$ -t 1-90
#$ -l s_vmem=16G
#$ -l mem_req=16G
#$ -o /dev/null
#$ -e /dev/null

samples=($(awk -F '\t' 'NR>1 {print $1}' data/sra_accession.tsv))
sample=${samples[$SGE_TASK_ID-1]}

shopt -s expand_aliases
alias fastp="apptainer exec /usr/local/biotools/f/fastp:0.23.4--h5f740d0_0 fastp"

workdir=~/popstr/fastaq
[ ! -e ${workdir}/cleaned ] && mkdir ${workdir}/cleaned
[ ! -e ~/popstr/out/qc ] && mkdir -p ~/popstr/out/qc

paired_fastaq=($(ls ${workdir}/uncleaned/${sample} | sort -V))
paired_1=${paired_fastaq[0]}
paired_2=${paired_fastaq[1]}

fastp \
  -i ${workdir}/uncleaned/${sample}/${paired_1} \
  -I ${workdir}/uncleaned/${sample}/${paired_2} \
  -o ${workdir}/cleaned/${sample}_clean_1.fq.gz \
  -O ${workdir}/cleaned/${sample}_clean_2.fq.gz \
  -h ~/popstr/out/qc/${sample}_qc.html \
  -q 30 -u 30

# rm -r ${workdir}/uncleaned/${sample}
