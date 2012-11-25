#!/bin/bash

## dos2unix ~/github/ms23/analysis/QC_SNP_count.sh ; bash ~/github/ms23/analysis/QC_SNP_count.sh

prefixu8=omni2.5-8_20120809_gwa_uganda_gtu
prefixa4=omni2.5-4_20120904_agv_gtu
prefixa8=omni2.5-8_agv_20120910_gtu

strand8=HumanOmni2.5-8v1_A-b37.strand
strand436=HumanOmni2.5M-b36.strand
strand437=HumanOmni2.5M-b37-v2.strand

prefix=$prefixa4
strand=$strand437

bim=$prefix.bim

## N.B.!!!
## bim X=23
## strand X=X

col_bim_id=2
col_bim_chrom=1
col_bim_pos=4
col_strand_id=1
col_strand_chrom=2
col_strand_pos=3

##
## sort
##
if [ ! -s $bim.sorted.id ]; then
cat $bim \
| awk '{print $0,$"'$col_bim_chrom'"":"$"'$col_bim_pos'"}' \
| sort -k$col_bim_id,$col_bim_id \
> $bim.sorted.id
fi

if [ ! -s $bim.sorted.chrompos ]; then
cat $bim \
| awk '{print $0,$"'$col_bim_chrom'"":"$"'$col_bim_pos'"}' \
| sort -k$col_bim_chrom,$col_bim_chrom -k$col_bim_pos,$col_bim_pos \
> $bim.sorted.chrompos
fi

if [ ! -s $strand.sorted.id ]; then
cmd="cat $strand \
| awk '{sub(/X/,23,\$2);sub(/Y/,24,\$2);sub(/XY/,25,\$2);sub(/MT/,26,\$2);\
print \$0,\$\"'\$col_strand_chrom'\"\":\"\$\"'\$col_strand_pos'\"}' \
| sort -k$col_strand_id,$col_strand_id \
> $strand.sorted.id"
echo $cmd
eval $cmd
fi

if [ ! -s $strand.sorted.chrompos ]; then
cmd="cat $strand \
| awk '{sub(/X/,23,\$2);sub(/Y/,24,\$2);sub(/XY/,25,\$2);sub(/MT/,26,\$2);\
print \$0,\$\"'\$col_strand_chrom'\"\":\"\$\"'\$col_strand_pos'\"}' \
| sort -k\$col_strand_chrom,\$col_strand_chrom -k\$col_strand_pos,\$col_strand_pos \
> \$strand.sorted.chrompos"
echo $cmd
eval $cmd
fi

##
## wc
##
echo

echo counts
echo bim
cat $bim | wc -l
echo strand
cat $strand | wc -l

##
## comm
##
echo

echo comm
echo comm by id
cat $bim.sorted.id | awk '{print $"'$col_bim_id'"}' > $bim.sorted.id.id
cat $strand.sorted.id | awk '{print $"'$col_strand_id'"}' > $strand.sorted.id.id
cmd="comm -12 $bim.sorted.id.id $strand.sorted.id.id > comm.id; cat comm.id | wc -l"
echo $cmd
eval $cmd

echo comm by chr pos
cat $bim.sorted.chrompos | awk '{print $7}' > $bim.sorted.chrompos.chrompos
cat $strand.sorted.chrompos | awk '{print $7}' > $strand.sorted.chrompos.chrompos
cmd="comm -12 $bim.sorted.chrompos.chrompos $strand.sorted.chrompos.chrompos > comm.chrompos; cat comm.chrompos | wc -l"
echo $cmd
eval $cmd

echo comm by chr pos id
cat $bim.sorted.id | awk '{print $7":"$"'$col_bim_id'"}' | sort > $bim.sorted.chromposid
cat $strand.sorted.id | awk '{print $7":"$"'$col_strand_id'"}' | sort > $strand.sorted.chromposid
cmd="comm -12 $bim.sorted.chromposid $strand.sorted.chromposid > comm.chromposid; cat comm.chromposid | wc -l"
echo $cmd
eval $cmd

##
##
##
echo 

echo "comm by chrompos but not chromposid"
awk 'BEGIN { FS=":"; OFS=":" } ; {print $1,$2}' comm.chromposid | sort > comm.chromposid.chrompos
sort comm.chrompos -o comm.chrompos
cmd="comm -13 comm.chromposid.chrompos comm.chrompos"
echo $cmd
eval $cmd

echo "comm by id but not chromposid"
awk 'BEGIN { FS=":"; OFS=":" } ; {print $3}' comm.chromposid | sort > comm.chromposid.id
sort comm.id -o comm.id
cmd="comm -13 comm.chromposid.id comm.id | wc -l"
echo $cmd
eval $cmd

##
## exit
##
if [ ! -f $prefix.position.SNPs.sorted ]; then
exit
fi

##
## sort
##
echo

if [ ! -s $bim.sorted.all ]
echo sort
then
## exclude chromosome
cat $bim.sorted.id | awk '{print $"'$col_bim_id'"}' > $bim.sorted.all
cat $bim.sorted.id | awk '{if($1>=1 && $1<=22) print $"'$col_bim_id'"}' > $bim.sorted.auto
cat $bim.sorted.id | awk '{if($1==23) print $"'$col_bim_id'"}' > $bim.sorted.X
cat $bim.sorted.chrompos | awk '{print $1":"$4}' > $bim.sorted.chrpos.all
cat $bim.sorted.chrompos | awk '{if($1>=1 && $1<=22) print $1":"$4}' > $bim.sorted.chrpos.auto
cat $bim.sorted.chrompos | awk '{if($1==23) print $1":"$4}' > $bim.sorted.chrpos.X
## exclusion
sort $prefix.position.SNPs > $prefix.position.SNPs.sorted
sort $prefix.miss.SNPs > $prefix.miss.SNPs.sorted
sort $prefix.duplicates.SNPs > $prefix.duplicates.SNPs.sorted
sort $prefix.duplicates.inclusion.SNPs > $prefix.duplicates.inclusion.SNPs.sorted
## resort by chr:pos instead of by chr and then pos
sort $strand.sorted.chrompos.chrompos > $strand.sorted.chrompos.chrompos.resorted
fi

##
## comm bim and strand (repeat of before...)
##
echo
echo comm bim strand
echo chrpos all strand
comm -12 $strand.sorted.chrompos.chrompos $bim.sorted.chrpos.all | wc -l
echo chrpos auto strand
comm -12 $strand.sorted.chrompos.chrompos $bim.sorted.chrpos.auto | wc -l
echo chrpos X strand
comm -12 $strand.sorted.chrompos.chrompos.resorted $bim.sorted.chrpos.X | wc -l
echo id all strand
comm -12 $strand.sorted.id.id $bim.sorted.id.all | wc -l
echo id auto strand
comm -12 $strand.sorted.id.id $bim.sorted.id.auto | wc -l
echo id X strand
comm -12 $strand.sorted.id.id.resorted $bim.sorted.id.X | wc -l

##
## comm bim and exclusion
##
echo
echo comm bim exclusion
for suffix in all auto X
do
echo
echo $suffix
echo pos
comm -12 $bim.sorted.$suffix $prefix.position.SNPs.sorted | wc -l
echo miss
comm -12 $bim.sorted.$suffix $prefix.miss.SNPs.sorted | wc -l
echo dup
comm -12 $bim.sorted.$suffix $prefix.duplicates.SNPs.sorted | wc -l
comm -12 $bim.sorted.$suffix $prefix.duplicates.inclusion.SNPs.sorted | wc -l
done

##
## comm exclusion1 and exclusion2
##
echo
echo comm exclusion1 exclusion2
echo comm position miss
comm -12 $prefix.position.SNPs.sorted $prefix.miss.SNPs.sorted | wc -l
echo comm position duplicates
comm -12 $prefix.position.SNPs.sorted $prefix.duplicates.SNPs.sorted | wc -l
comm -12 $prefix.position.SNPs.sorted $prefix.duplicates.inclusion.SNPs.sorted | wc -l
echo comm miss duplicates
comm -12 $prefix.miss.SNPs.sorted $prefix.duplicates.SNPs.sorted | wc -l
comm -12 $prefix.miss.SNPs.sorted $prefix.duplicates.inclusion.SNPs.sorted | wc -l

##
## 
##
echo
echo comm exclusion1 exclusion2 exclusion3
echo count intersect auto pos miss
comm -12 $prefix.position.SNPs.sorted $prefix.miss.SNPs.sorted | sort > comm.position.miss.SNPs.sorted
comm -12 $bim.sorted.auto comm.position.miss.SNPs.sorted | wc -l
echo print intersect pos dupexcl
comm -12 $prefix.position.SNPs.sorted $prefix.duplicates.SNPs.sorted
echo print intersect pos dupincl
comm -12 $prefix.position.SNPs.sorted $prefix.duplicates.inclusion.SNPs.sorted
