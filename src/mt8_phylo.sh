#$ -S /bin/bash
#$ -cwd
#$ -l medium
#$ -l mem_req=512G
#$ -l s_vmem=512G
#$ -o /dev/null
#$ -e /dev/null


shopt -s expand_aliases
alias rapidnj="apptainer exec /usr/local/biotools/r/rapidnj:2.3.2--h4ac6f70_4 rapidnj"
alias iqtree2="apptainer exec /usr/local/biotools/i/iqtree:2.2.6--h21ec9f0_0 iqtree2"

workdir=~/popstr/phylo/mito
vcf=~/popstr/vcf/popstr.mt.vcf.gz
prefix=popstr.mt
model=MFP

[ ! -e ${workdir} ] && mkdir -p ${workdir}
cd ${workdir}

/usr/bin/python3 ~/bin/vcf2phylip/vcf2phylip.py -i ${vcf} --fasta
/usr/bin/python3 ../../src/sort_mitochondria_seq.py -i ${prefix}.min4.fasta -o ${prefix}.min4.sort.fasta
rm ${prefix}.min4.fasta ${prefix}.min4.phy

rapidnj ${prefix}.min4.sort.fasta -i fa -o t -b 100 --evolution-model kim -x ${prefix}.nj.tree
iqtree2 -s ${prefix}.min4.sort.fasta -m ${model} -B 1000 -nt AUTO --prefix ${prefix} -redo

/usr/bin/python3 ~/bin/MitoToolPy/mitotoolpy-seq.py -s chicken -r whole -i ${prefix}.min4.sort.fasta -o stdout > ${prefix}.mitotool
/usr/bin/python3 ~/bin/MitoToolPy/mitotoolpy-seq.py -s chicken -r dloop -i ${prefix}.min4.sort.fasta -o stdout > ${prefix}.D.mitotool
