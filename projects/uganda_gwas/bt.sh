bcftools=/software/hgi/pkglocal/bcftools-1.2/bin/bcftools
tabix=/software/hgi/pkglocal/htslib-1.2.1/bin/tabix

chrom=$1

o=out_bt_concat/$chrom.vcf.gz
o=out_bt/$chrom

if [ $(ls ../pipeline_UG3.3/out_UnifiedGenotyper/$chrom/*.vcf.gz | wc -l) -ne $(ls ../pipeline_UG3.3/out_UnifiedGenotyper/$chrom/*.vcf.gz.tbi | wc -l) ]; then exit; fi

mkdir -p $(dirname $o)
if [ -f $o ]; then exit; fi
touch $o

$bcftools concat \
 $(ls ../pipeline_UG3.3/out_UnifiedGenotyper/$chrom/*.vcf.gz | sort -V) \
 -Ou | $bcftools view -M2 -m2 \
 -Ou | $bcftools convert --gensample $o \

# -Oz -o $o \

#$tabix -p vcf $o

## Add the trio information to the .samples file.
cat "ID_1 ID_2 missing ID_father ID_mother
0 0 0
EGAN00001160764 EGAN00001160764 0 EGAN00001160765 EGAN00001160766
EGAN00001160765 EGAN00001160765 0
EGAN00001160766 EGAN00001160766 0" > $o.samples
