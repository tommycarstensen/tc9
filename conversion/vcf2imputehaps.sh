fgrep -v "#" Uganda_1687_phased_SLRP_nomissing_nounphased.vcf | sed 's/:\S*//g' | sed 's/|/ /g' | sed 's/\t/ /g' | sed 's/\./?/g' | awk '{print $1, $1":"$2, $0}' | awk ' { for (i=1; i<=NF;i++) if( i!=3 && i!=5 && i!=8 && i!=9 && i!=10 && i!=11) printf("%s%s", $i,(i!=NF) ? OFS : ORS)}' > Uganda_1687_phased_SLRP_nomissing_nounphased.haps


fgrep "CHR" Uganda_1687_phased_SLRP.vcf | sed 's/\t/\n/g' | fgrep "APP" | awk -F: '{print $2}' > Uganda_1687_phased_SLRP_nomissing_nounphased.sample

