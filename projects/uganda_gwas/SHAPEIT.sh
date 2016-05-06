SHAPEIT=/software/team149/shapeit.v2.r790

#https://mathgen.stats.ox.ac.uk/genetics_software/shapeit/shapeit.html#reference

chrom=$1

map=/lustre/scratch114/resources/mathgen.stats.ox.ac.uk/impute/1000GP_Phase3/genetic_map_chr${chrom}_combined_b37.txt
if [ ! -s $map ]; then echo map; echo $map; exit; fi

ref=/lustre/scratch114/projects/ug2g/users/tc9/IMPUTE

gen=out_bt/$chrom
out=out_SHAPEIT/chrom$chrom.pos$pos1-$pos2
mkdir -p $(dirname $out)
if [ -f $out.phased.haps.gz ]; then exit; fi

cmd="$SHAPEIT \
 --input-map $map \
 --output-max $out.phased.haps.gz $out.phased.sample \
 --output-log $out.log \
 --effective-size 17469 \
 --thread 16 \
 --input-gen $gen.gen.gz $gen.samples \
 --input-ref \
  $ref/out_merge_ref_panels/1000Gp3_agv_ug2g/$chrom.hap.gz \
  $ref/out_merge_ref_panels/1000Gp3_agv_ug2g/$chrom.legend.gz \
  $ref/1000Gp3_agv_ug2g.SHAPEIT2.sample \
 --exclude-snp out_SHAPEIT_check/$chrom.snp.strand.exclude \
"

## chrom X
if [ "$chrom" == "X" ]; then
# cmd=$cmd" --input-bed $bfile.bed $bfile.bim $bfile.fam"
 cmd=$cmd" --chrX"
fi

## eval
eval $cmd
