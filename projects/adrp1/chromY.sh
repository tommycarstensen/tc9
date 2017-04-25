#1. Get the major branch SNPs from these two files
#https://raw.githubusercontent.com/23andMe/yhaplo/master/input/representative.SNPs.isogg.2015tree.txt
#https://raw.githubusercontent.com/23andMe/yhaplo/master/input/representative.SNPs.additional.txt
#
#Chromosome coordinates for these sites could be get from this file
#https://raw.githubusercontent.com/23andMe/yhaplo/master/input/isogg.2016.01.04.txt

## sed 's/ \+\t/\t/g' to remove trailing white space from each field
## http://stackoverflow.com/questions/20600982/remove-leading-and-trailing-space-in-field-in-awk

curl -s https://raw.githubusercontent.com/23andMe/yhaplo/master/input/isogg.2016.01.04.txt | sed 's/ \+\t/\t/g' | awk -F$'\t' 'BEGIN{OFS="\t"} {if(NF==6) {print $0} else {print $1,$3,$4,$5,$6,$7}}' > isogg.2016.01.04.trimmed.txt

join -t$'\t' \
 <(curl -s \
  https://raw.githubusercontent.com/23andMe/yhaplo/master/input/representative.SNPs.isogg.2015tree.txt \
  https://raw.githubusercontent.com/23andMe/yhaplo/master/input/representative.SNPs.additional.txt \
  | awk '{print $1,$2}' | grep -v ^[C-D,F-Z][1-3] | cut -d" " -f2 | tr "," "\n" | tr "/" "\n" | sort)  \
 <(cat isogg.2016.01.04.trimmed.txt \
  | awk -F$'\t' 'BEGIN{OFS="\t"}{print $1,$(NF-4),int($(NF-1)),substr($NF,1,1),substr($NF,4,1)}' | sort -t$'\t' -k1,1) \
| awk -F$'\t' 'BEGIN{OFS="\t"} {print $2,$1,$3,$4,$5,1}'

#2. Include any SNPs (if allowed) for haplogroup A, B, E (Africa Y chromosome haplogroups) from file isogg.2016.01.04.txt <https://github.com/23andMe/yhaplo/blob/master/input/isogg.2016.01.04.txt> (excluding indels)
#
#2) Select at least 3 diagnostic SNPs for all the African Y haplogroups (Haplogroup A, B, E and African R1b) to the finest resolution we have for the moment according to ISOGG.

cat isogg.2016.01.04.trimmed.txt | grep -v Investigation | awk -F$'\t' 'BEGIN{OFS="\t"} length($6)==4&&substr($6,2,1)=="-"&&substr($6,3,1)==">"&&(substr($2,1,1)=="A"||substr($2,1,1)=="B"||substr($2,1,1)=="E"){print $2,$1,$5,substr($6,1,1),substr($6,4,1)}' | shuf | awk 'BEGIN{OFS="\t"} {cnt[$1]++; if(cnt[$1]<=5) {print $0,2}}' | sort -k1V,1 -k3n,3

#3) We should also include up to 3 diagnostic SNPs for all of other non-African Y haplogroup to the resolution of the three letter level, to identify male admixture from non-African populations.

cat isogg.2016.01.04.trimmed.txt | grep -v Investigation | grep -v Withdrawn | awk -F$'\t' 'BEGIN{OFS="\t"} length($6)==4&&substr($6,2,1)=="-"&&substr($6,3,1)==">"&&(substr($2,1,1)!="A"&&substr($2,1,1)!="B"&&substr($2,1,1)!="E")&&length($2)<=4{print $2,$1,$5,substr($6,1,1),substr($6,4,1)}' | shuf | awk 'BEGIN{OFS="\t"} {cnt[$1]++; if(cnt[$1]<=3) {print $0,3}}' | sort -k1V,1 -k3n,3

#4. Add all of the African R1b haplogroup SNPs identified from 1000G project below

cat African_R1b.tsv | awk '{print "R1b\t""\t"$1"\t"$2"\t"$3"\t"1}'

echo "R1b1c V88 4862861 C T 1" | tr " " "\t"
