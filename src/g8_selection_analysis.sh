#$ -S /bin/bash
#$ -cwd
#$ -o /dev/null
#$ -e /dev/null

shopt -s expand_aliases
alias vcftools="apptainer exec /usr/local/biotools/v/vcftools:0.1.16--h9a82719_5 vcftools"
alias bcftools="apptainer exec /usr/local/biotools/b/bcftools:1.18--h8b25389_0 bcftools"

vcf=~/popstr/vcf/popstr.snp.vcf.gz
workdir=~/popstr/selection

[ ! -e ${workdir} ] && mkdir ${workdir}
cd ${workdir}

bcftools query -l ${vcf} | grep -E 'SRR|ERR' > other.txt
bcftools query -l ${vcf} | grep -v -E 'SRR|ERR|SM|YKD|WKNP' > jnc.txt
bcftools query -l ${vcf} | grep -E 'YKD|SM' > shamo.txt
# cat ~/popstr/data/sra_accession.tsv | awk -F '\t' 'NR > 1 && $2 ~ /^Gg_/ {print $1}' > ggallus.txt
# cat ~/popstr/data/sra_accession.tsv | awk -F '\t' 'NR > 1 && $2 !~ /^Gg_/ {print $1}' > others.txt


## Windowed pie

pie_dir=${workdir}/pie
[ ! -e ${pie_dir} ] && mkdir ${pie_dir}

vcftools --gzvcf ${vcf} \
  --keep jnc.txt \
  --keep shamo.txt \
  --window-pi 20000 \
  --window-pi-step 10000 \
  --out ${pie_dir}/japan

vcftools --gzvcf ${vcf} \
  --keep other.txt \
  --window-pi 20000 \
  --window-pi-step 10000 \
  --out ${pie_dir}/other


## Fst

fst_dir=${workdir}/fst
[ ! -e ${fst_dir} ] && mkdir ${fst_dir}

# echo "Start calcurating Fst..."
vcftools --gzvcf ${vcf} --weir-fst-pop jnc.txt --weir-fst-pop other.txt \
  --fst-window-size 20000 --fst-window-step 10000 --out ${fst_dir}/jnc_other

vcftools --gzvcf ${vcf} --weir-fst-pop other.txt --weir-fst-pop shamo.txt \
  --fst-window-size 20000 --fst-window-step 10000 --out ${fst_dir}/other_shamo

vcftools --gzvcf ${vcf} --weir-fst-pop shamo.txt --weir-fst-pop jnc.txt \
  --fst-window-size 20000 --fst-window-step 10000 --out ${fst_dir}/shamo_jnc


## Tajima's D

tajima_dir=${workdir}/tajimaD
[ ! -e ${tajima_dir} ] && mkdir ${tajima_dir}

# echo "Start calcurating Tajima's D statistics..."
vcftools --gzvcf ${vcf} --keep jnc.txt --TajimaD 20000 --out ${tajima_dir}/jnc
vcftools --gzvcf ${vcf} --keep shamo.txt --TajimaD 20000 --out ${tajima_dir}/shamo
vcftools --gzvcf ${vcf} --keep other.txt --TajimaD 20000 --out ${tajima_dir}/other


# echo "Start calcurating XP-EPP ..."

alias selscan="apptainer exec /usr/local/biotools/s/selscan:1.2.0a--h0fdf51a_4 selscan"
alias plink2="apptainer exec /usr/local/biotools/p/plink2:2.00a5--h4ac6f70_0 plink2"

vcf=~/popstr/vcf/varsites.snp.vcf.gz
workdir=~/popstr/xpehh

[ ! -e ${workdir} ] && mkdir ${workdir}
cd ${workdir}

plink2 --vcf ${vcf} \
  --allow-extra-chr \
  --double-id \
  --set-missing-var-ids @:# \
  --export ped \
  --out tmp
cat tmp.map | awk -F'\t' '{print $1 "\t" $2 "\t" NR "\t" $4}' > tmp2.map

bcftools query -l ${vcf} | grep -E 'SRR|ERR' > other.txt
bcftools query -l ${vcf} | grep -v -E 'SRR|ERR|SM|YKD|WKNP|WLHP|WPRP|RIRP|AKNP' > jnc.txt

bcftools view ${vcf} -S jnc.txt -Oz -o jnc.vcf.gz
bcftools index jnc.vcf.gz

bcftools view ${vcf} -S other.txt -Oz -o ref.vcf.gz
bcftools index ref.vcf.gz

selscan --xpehh --vcf jnc.vcf.gz --vcf-ref ref.vcf.gz --map tmp2.map --out jnc.xpehh
