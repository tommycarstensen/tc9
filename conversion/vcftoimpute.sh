## vcf > haps
fgrep -v "#" Uganda_1687_phased_SLRP_nomissing_nounphased.vcf \
      ## replace colon\S* with nothing
      | sed 's/:\S*//g' \
      ## replace pipe with space
      | sed 's/|/ /g' \
      ## replace tab with space
      | sed 's/\t/ /g' \
      ## replace . with question mark
      | sed 's/\./?/g' \
      ## duplicate column1, column1:column2
      | awk '{print $1, $1":"$2, $0}' \
      ## exclude selected columns
      | awk ' { for (i=1; i<=NF;i++) if( i!=3 && i!=5 && i!=8 && i!=9 && i!=10 && i!=11) printf("%s%s", $i,(i!=NF) ? OFS : ORS)}' \
      > Uganda_1687_phased_SLRP_nomissing_nounphased.haps

## vcf > sample
fgrep "CHR" Uganda_1687_phased_SLRP.vcf \
      | sed 's/\t/\n/g' \
      | fgrep "APP" \
      | awk -F: '{print $2}' \
      > Uganda_1687_phased_SLRP_nomissing_nounphased.sample

