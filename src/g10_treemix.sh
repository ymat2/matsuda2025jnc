#$ -S /bin/bash
#$ -cwd
#$ -o /dev/null
#$ -e /dev/null


shopt -s expand_aliases
alias plink1="apptainer exec /usr/local/biotools/p/plink:1.90b4--0 plink"
alias plink2="apptainer exec /usr/local/biotools/p/plink2:2.00a5--h4ac6f70_0 plink2"
alias bcftools="apptainer exec /usr/local/biotools/b/bcftools:1.18--h8b25389_0 bcftools"
alias treemix="apptainer exec /usr/local/biotools/t/treemix:1.13--h125836d_9 treemix"

workdir=~/popstr/treemix
vcf=~/popstr/vcf/popstr.snp.vcf.gz
bfile=~/popstr/structure/popstr.snp.filt
prefix=popstr.sub

[ ! -e ${workdir} ] && mkdir -p ${workdir}
cd ${workdir}


## Make .clust file

cat ../data/sra_accession.tsv | awk 'NR>1 {print $1 "\t" $1 "\t" $2}' > popstr.clust
bcftools query -l ${vcf} | grep -E 'YKD|SM' | awk '{print $1 "\t" $1 "\tShamo"}' >> popstr.clust
bcftools query -l ${vcf} | grep -v -E 'YKD|SM|RR'  | awk '{print $1 "\t" $1 "\tJapan"}' >> popstr.clust
# Need to be edited

awk '{print $3}' popstr.clust | sort | uniq > treemix.list

## Make Treemix input file

plink1 --bfile ${bfile} \
  --freq \
  --missing \
  --within popstr.clust \
  --not-chr NC_052571.1, NC_052572.1, NC_053523.1 \
  --allow-extra-chr \
  --double-id \
  --set-missing-var-ids @:# \
  --out ${prefix}

# gzip ${prefix}.frq.strat

/usr/bin/python3 ../src/plink2treemix.py -i ${prefix}.frq.strat -o ${prefix}.treemix.frq
rm ${prefix}.*miss ${prefix}.frq.strat ${prefix}.nosex
gzip ${prefix}.treemix.frq


for m in {1..8}; do
  for i in {1..5}; do
    echo Running Treemix: M=${m}, iteration=${i} ...
    treemix -i ${prefix}.treemix.frq.gz -m ${m} -o treemix.${i}.${m} -root RJF -bootstrap -k 500 -seed ${i} > treemix_${i}.${m}_log
  done
done
