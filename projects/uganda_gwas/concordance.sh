#!/bin/bash

bcftools=/lustre/scratch115/teams/sandhu/software/bcftools/bcftools

calls20=/lustre/scratch114/projects/ug2g/users/tc9/vcfeval/out_bt_concat/20.vcf.gz

ref=/lustre/scratch114/resources/ref/Homo_sapiens/1000Genomes_hs37d5/hs37d5.fa

for v in indels snps,mnps; do

 echo $v

 if [ $v == "indels" ]; then V=snps,mnps; else V=indels; fi

 baseline=/lustre/scratch114/projects/ug2g/users/tc9/vcfeval/out_bt_view/chrom1-22/NIST_RTG_PlatGen_merged_highconfidence_v0.2_Allannotate.$v.vcf.gz

 join -a1 \
  <(\
   $bcftools view -r 20 --trim-alt-alleles $baseline \
   | $bcftools norm -f $ref \
   | $bcftools view -V $V \
   | $bcftools annotate -x FORMAT \
   | grep -v ^# \
   | awk '{split($5,a,","); a[0]=$4; A1=a[substr($10,1,1)]; A2=a[substr($10,3,1)]; printf("%09d %s %s\n", $2, A1, A2)}') \
  <(\
   $bcftools view --trim-alt-alleles $calls20 \
   | $bcftools norm -f $ref \
   | $bcftools view -V $V \
   | $bcftools annotate -x FORMAT \
   | grep -v ^# \
   | awk '{split($5,a,","); a[0]=$4; A1=a[substr($10,1,1)]; A2=a[substr($10,3,1)]; printf("%09d %s %s\n", $2, A1, A2)}') \
 | awk '{if(($2==$4&&$3==$5)||($2==$5&&$3==$4)) {conc++} else {disc++}} END{print conc/(conc+disc)}'

done

