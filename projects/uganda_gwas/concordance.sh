#baseline=/lustre/scratch114/projects/ug2g/users/tc9/vcfeval/out_bt_view/chrom1-22/NIST_RTG_PlatGen_merged_highconfidence_v0.2_Allannotate.indels.vcf.gz
#baselineSNPs=/lustre/scratch114/projects/ug2g/users/tc9/vcfeval/out_bt_view/chrom1-22/NIST_RTG_PlatGen_merged_highconfidence_v0.2_Allannotate.snps,mnps.vcf.gz

bcftools=/lustre/scratch115/teams/sandhu/software/bcftools/bcftools

calls20=/lustre/scratch114/projects/ug2g/users/tc9/vcfeval/out_bt_concat/20.vcf.gz

baseline=ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/NA12878_HG001/NISTv3.3.2/GRCh37/HG001_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_PGandRTGphasetransfer.vcf.gz
baseline=$(basename $baseline)

ref=/lustre/scratch114/resources/ref/Homo_sapiens/1000Genomes_hs37d5/hs37d5.fa

bed=ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/NA12878_HG001/NISTv3.3.2/GRCh37/HG001_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_nosomaticdel.bed
bed20=$(basename $bed .bed).chrom20.bed
if [ ! -s $bed20 ]; then curl -s $bed | grep ^20 > $bed20; fi

for v in indels snps,mnps; do

echo $v

if [ $v == "indels" ]; then V=snps,mnps; else V=indels; fi

# -R $bed20

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

exit

#$bcftools stats -s - \
# $baseline \
# $calls \
#| grep -v ^DP | grep -v ^IDD | grep -v ^QUAL

## The extra awk code is to split MNPs into SNPs
## and avoid printing these,
## if REF/ALT or ALT/ALT is identical to REF/REF.
join -a1 -1 2 -2 2 \
 <(\
  $bcftools view -r 20 --trim-alt-alleles $baselineSNPs \
  | $bcftools norm -f $ref \
  | $bcftools view -v snps,mnps \
  | $bcftools annotate -x FORMAT \
  | grep -v ^# \
  | awk '{
   split($5,a,","); for(i=1;i<=length($4);i++) {
    printf("%s %09d %s ", $1, $2+i-1, substr($4,i,1));
    printf(substr(a[1],i,1));
    for(j=2;j<=length(a);j++) {
     printf(","substr(a[j],i,1))}; printf " "$10"\n"}}' \
  | awk '{
   split($4,a,","); a[0]=$3;
   if(a[substr($5,1,1)-1]!=a[0]||a[substr($5,3,1)-1]!=a[0]) print $0}' \
  ) \
<(\
  $bcftools view --trim-alt-alleles $calls \
  | $bcftools norm -f $ref \
  | $bcftools view -v snps,mnps \
  | $bcftools annotate -x FORMAT \
  | grep -v ^# \
  | awk '{
   split($5,a,","); for(i=1;i<=length($4);i++) {
    printf("%s %09d %s ", $1, $2+i-1, substr($4,i,1));
    printf(substr(a[1],i,1));
    for(j=2;j<=length(a);j++) {
     printf(","substr(a[j],i,1))}; printf " "$10"\n"}}' \
  | awk '{
   split($4,a,","); a[0]=$3;
   if(a[substr($5,1,1)-1]!=a[0]||a[substr($5,3,1)-1]!=a[0]) print $0}' \
  ) \
\
| awk 'NF>5{split($4,a,","); split($8,b,","); if(a[2]==b[1]||a[1]==b[2]) print}'
| awk '{if($3!=$7||$4!=$8) {GTdiff++; print $5,"./."} else {print $5,$9}}' \
| awk '{if((substr($1,1,1)==substr($2,1,1)&&substr($1,3,1)==substr($2,3,1))||(substr($1,1,1)==substr($2,3,1)&&substr($1,3,1)==substr($2,1,1))) {conc++} else {disc++}} END{print conc/(conc+disc)}'

join -a1 -1 2 -2 2 \
 <(\
  $bcftools view -r 20 --trim-alt-alleles $baseline \
  | $bcftools norm -f $ref \
  | $bcftools view -v indels \
  | $bcftools annotate -x FORMAT \
  | grep -v ^# \
  | awk '{printf("%s %09d %s %s %s\n", $1, $2, $4, $5, $10)}') \
 <(\
  $bcftools view --trim-alt-alleles $calls \
  | $bcftools norm -f $ref \
  | $bcftools view -v indels \
  | $bcftools annotate -x FORMAT \
  | grep -v ^# \
  | awk '{printf("%s %09d %s %s %s\n", $1, $2, $4, $5, $10)}') \
| awk '{if($3!=$7||$4!=$8) {GTdiff++; print $5,"./."} else {print $5,$9}}' \
| awk '{if((substr($1,1,1)==substr($2,1,1)&&substr($1,3,1)==substr($2,3,1))||(substr($1,1,1)==substr($2,3,1)&&substr($1,3,1)==substr($2,1,1))) {conc++} else {disc++}} END{print conc/(conc+disc)}'

#| sort | uniq -c | sort -k1n,1
