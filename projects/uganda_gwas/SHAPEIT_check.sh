SHAPEIT=/software/team149/shapeit.v2.r790

chrom=$1

map=/lustre/scratch114/resources/mathgen.stats.ox.ac.uk/impute/1000GP_Phase3/genetic_map_chr${chrom}_combined_b37.txt

ref=/lustre/scratch114/projects/ug2g/users/tc9/IMPUTE

#https://mathgen.stats.ox.ac.uk/genetics_software/shapeit/shapeit.html#reference
#Step2: Alignment of the SNPs between the GWAS dataset and the reference panel
mkdir -p out_SHAPEIT_check
$SHAPEIT \
 -check \
 --input-gen \
  out_bt/$chrom.gen.gz \
  out_bt/$chrom.samples \
 --input-map $map \
 --input-ref \
  $ref/out_merge_ref_panels/1000Gp3_agv_ug2g/$chrom.hap.gz \
  $ref/out_merge_ref_panels/1000Gp3_agv_ug2g/$chrom.legend.gz \
  $ref/1000Gp3_agv_ug2g.SHAPEIT2.sample \
 --output-log out_SHAPEIT_check/$chrom \

#cmd=$cmd" --exclude-snp out_SHAPEIT_check/$chrom.snp.strand.exclude"
#eval $cmd

##Mendel errors
#$SHAPEIT \
# -check \
# --input-gen \
#  out_bt/$chrom.gen.gz \
#  out_bt/$chrom.samples \
# --input-map $map \
# --output-log out_SHAPEIT_check/$chrom \

