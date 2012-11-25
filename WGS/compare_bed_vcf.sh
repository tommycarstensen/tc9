#!/bin/bash

## T. Carstensen (tc9), M.S. Sandhu (ms23), D. Gurdasani (dg11)
## Wellcome Trust Sanger Institute, 2012

##
## compare genotype BED to WGS VCFs and GENs
##

gtool='/nfs/team149/Software/usr/share/gtool/gtool'
sepjoin='join'
echo $sepjoin

##
## 1) bed to ped and map
##

## 1.1) bed > ped, map
## keep - keep selected samples/individuals
## exclude - exclude SNPs from PED file (ped > ped)
## http://pngu.mgh.harvard.edu/~purcell/plink/dataman.shtml#exclude
in=/lustre/scratch107/projects/agv/phasing_rel/shapeit2/data/omni2.5-8_20120516_gwa_ugand_gtu_autosomes_postsampqc_postsnpqc_flipped
out=omni2.5-8_20120516_gwa_ugand_gtu_autosomes_postsampqc_postsnpqc_flipped
prefixomni=omni2.5-8_20120516_gwa_ugand_gtu_autosomes_postsampqc_postsnpqc_flipped
if [ ! -s backup/$out.ped ]
then
bsub \
     -q normal -M4000000 -R'select[mem>4000] rusage[mem=4000]' -P uganda -o LSF_plink_recode.out -e LSF_plink_recode_normal.err \
plink \
     --noweb \
     --bfile $in \
     --recode \
     --out backup/$out \
     --keep keep_long_double_column.txt \
     --exclude exclude.rsid \

fi

##
## 1.2) exclude SNPs
##
if [ ! -s $out.ped ] && [ ! -s $out.map ]
then

## 1.2) exclude SNPs from PED file (ped > ped)
## http://pngu.mgh.harvard.edu/~purcell/plink/dataman.shtml#exclude
plink \
    --noweb \
    --file backup/$out \
    --recode \
    --out $out \
    --exclude exclude.rsid \

    
## 1.3a) rename individuals in *ped* file
cat \
    $out.ped \
    | cut -c 12-22,34- \
    > $out.ped.tmp; \
	mv $out.ped.tmp $out.ped

## 1.3b) rename markers in *map* file
awk -F '\t' '{print $1, $1":"$4, $3, $4}' \
    $out.map \
    | uniq \
    > $out.map.tmp; \
	mv $out.map.tmp $out.map

fi

#for CHROMOSOME in {1..22} X Y
for CHROMOSOME in {1..22}
do
    echo $CHROMOSOME

    ##
    ## 2) vcf or gen to ped
    ##

    ## 2.1a.1) vcf > ped, map
	if [ "$sepjoin" == "join" ]
	then
	in=out_GATK/$sepjoin/ApplyRecalibration.recalibrated.filtered.vcf
	else
    in=out_GATK/$sepjoin/ApplyRecalibration.recalibrated.filtered.$CHROMOSOME.vcf
	fi
    out=out_VCFtools/$sepjoin/ApplyRecalibration.$CHROMOSOME
	if [ ! -s $out.ped ]
	then
    cmd="bsub -q normal -M2000000 -R'select[mem>2000] rusage[mem=2000]' -P uganda \
	    -o stdout/vcftools.out -e stderr/vcftools.err \
        vcftools --vcf $in --chr $CHROMOSOME --plink --out $out"
	echo $cmd
	eval $cmd
    fi

	##
    ## 2.1b.1) gen2gen; merge IMPUTE2 output files and get rid of duplicate SNPs and INDELs
	##
	if [ ! -s out_IMPUTE2/$sepjoin/IMPUTE2.$CHROMOSOME.gen ]
	then
	echo "IMPUTE2 gen merge and mod (gen > gen)"
	## find all files to be merged (cannot do cat prefix* because of other non-gen files in dir)
	cmd="cat "
	for i in {1..10000}
	do
	if [ ! -s out_IMPUTE2/$sepjoin/chromosome$CHROMOSOME.impute2.gen.${i}_summary ]
	then
	break
	fi
	if [ -s out_IMPUTE2/$sepjoin/chromosome$CHROMOSOME.impute2.gen.$i ]
	then
	cmd=$cmd" out_IMPUTE2/$sepjoin/chromosome$CHROMOSOME.impute2.gen.$i"
	fi
	done
	## only consider SNPs (exclude rows) (and exclude duplicates)
    cmd=$cmd" | awk 'x[\$3]++==0&&length(\$4)==1&&length(\$5)==1'"
	## remove columns --- and rsID
	cmd=$cmd" | cut -f3- -d ' '"
	## add columns CHR:POS CHR:POS
	cmd=$cmd" | awk -v CHROMOSOME=\"$CHROMOSOME\" '{print CHROMOSOME\":\"\$1\" \"CHROMOSOME\":\"\$1\" \"\$0}'"
	## specify file output
    cmd=$cmd" > out_IMPUTE2/$sepjoin/IMPUTE2.$CHROMOSOME.gen"
	echo $cmd
	eval $cmd
	fi

    # ## 2.1b.2) get rid of a single SNP causing trouble...
	# if [ $CHROMOSOME -eq "18" ]
	# then
    # grep -v 13489156 out_IMPUTE2/IMPUTE2.18.gen > IMPUTE2.18.gen2; mv IMPUTE2.18.gen2 out_IMPUTE2/IMPUTE2.18.gen
	# fi
    
	##
	## bash BEAGLE gprobs to gen conversion
    ## 2.1c.1) gprobs > gen (remove header)
	##
	if [ ! -s out_BEAGLE/$sepjoin/BeagleOutput.$CHROMOSOME.gen ]
	then
	echo "gprobs > gen"
	## skip header
	in="out_BEAGLE/$sepjoin/BeagleOutput.${CHROMOSOME}.bgl.ProduceBeagleInput.${CHROMOSOME}.bgl.gprobs"
	out="out_BEAGLE/$sepjoin/BeagleOutput.${CHROMOSOME}.gen"
	s_more="more +2 $in"
	s_tail="tail -n+2 $in"
	s_sed="sed 1d $in"
	## print columns 1 and 2 when using field separator : (i.e. a replacement of field operator)
	s_awk1="awk -F : 'NR!=1 {print \$1, \$2}' $in"
	## print columns 1 and 2 with default field separator (i.e. append columns)
	s_awk2="awk '{print \$1\":\"\$2, \$1\":\"\$2, \$0}'"
	## cut out column 3 and do *not* print this column (i.e. remove "marker" from header and chromosomeID from subsequent rows)
	s_cut="cut -d \" \" -f 3 --complement > $out"
	s_convert="${s_tail} | ${s_awk1} | ${s_awk2} | ${s_cut}"
	echo $s_convert
	eval $s_convert
	fi

	##
	## GTOOL GEN to PED Conversion
	## http://www.well.ox.ac.uk/~cfreeman/software/gwas/gtool.html#GEN to PED Conversion
    ## 2.1b.3 / 2.1c.2) gen and sample > ped and map
	##
	sample="out_BEAGLE/$sepjoin/chromosome1.sample"
    for PREFIX in BEAGLE IMPUTE2
	do
		if [ ! -s out_GTOOL/$sepjoin/$PREFIX.$CHROMOSOME.ped ]
		then
			echo "bsub" $PREFIX $CHROMOSOME
			if [ "$PREFIX" == "BEAGLE" ]
			then
				prefixgen="BeagleOutput"
				else
				prefixgen="IMPUTE2"
			fi

			cmd="bsub \
			-q normal -M4000000 -R'select[mem>4000] rusage[mem=4000]' -P uganda \
			-o stdout/GTOOL.$PREFIX.$CHROMOSOME.out \
			-e stderr/GTOOL.$PREFIX.$CHROMOSOME.err \
			-J\"GTOOL.$PREFIX.$CHROMOSOME\" \
				$gtool -G --g out_$PREFIX/$sepjoin/$prefixgen.$CHROMOSOME.gen \
				--s $sample \
				--ped out_GTOOL/$sepjoin/$PREFIX.$CHROMOSOME.ped \
				--map out_GTOOL/$sepjoin/$PREFIX.$CHROMOSOME.map --threshold 0.9 --snp \
				--log out_GTOOL/$sepjoin/$PREFIX.$CHROMOSOME.log"
			echo $cmd
			eval $cmd
		fi
	done

done

echo "loop2"
#for CHROMOSOME in {1..22} X Y
for CHROMOSOME in {1..22}
do

    ## 2.2a) map > map; add chromosome column to gtool map file
    ## BEAGLE
	if [ -s out_GTOOL/$sepjoin/BEAGLE.$CHROMOSOME.map ] && [ ! -s out_GTOOL/$sepjoin/BEAGLE.$CHROMOSOME.map2 ]
	then
	echo 'BEAGLE map mod (map > map2)' $CHROMOSOME
    awk -F : '{print $1"\t"$2}' out_GTOOL/$sepjoin/BEAGLE.$CHROMOSOME.map | awk -F "\t" '{print $2"\t"$2":"$3"\t"$4"\t"$5}' \
        > BEAGLE.$CHROMOSOME.map2 \
        ; \
        mv BEAGLE.$CHROMOSOME.map2 out_GTOOL/$sepjoin/BEAGLE.$CHROMOSOME.map2
	fi
    ## IMPUTE2
    if [ -s out_GTOOL/$sepjoin/IMPUTE2.$CHROMOSOME.map ] && [ ! -s out_GTOOL/$sepjoin/IMPUTE2.$CHROMOSOME.map2 ]
	then
	echo 'IMPUTE2 map mod (map > map2)' $CHROMOSOME
    awk -F "\t" -v CHROMOSOME="$CHROMOSOME" '{print CHROMOSOME"\t"CHROMOSOME":"$4"\t"$3"\t"$4}' out_GTOOL/$sepjoin/IMPUTE2.$CHROMOSOME.map \
        > IMPUTE2.$CHROMOSOME.map2 \
        ; \
        mv IMPUTE2.$CHROMOSOME.map2 out_GTOOL/$sepjoin/IMPUTE2.$CHROMOSOME.map2
	fi
    ## ApplyRecalibration
    if [ ! -s out_VCFtools/$sepjoin/ApplyRecalibration.$CHROMOSOME.map2 ] && [ -s out_VCFtools/$sepjoin/ApplyRecalibration.$CHROMOSOME.map ]
	then
	echo 'ApplyRecalibration map mod (map > map2)' $CHROMOSOME
    awk -F "\t" '{print $1"\t"$1":"$4"\t"$3"\t"$4}' out_VCFtools/$sepjoin/ApplyRecalibration.$CHROMOSOME.map \
        > ApplyRecalibration.$CHROMOSOME.map2 \
        ; \
        mv ApplyRecalibration.$CHROMOSOME.map2 out_VCFtools/$sepjoin/ApplyRecalibration.$CHROMOSOME.map2
	fi

    ## 2.2b) ped2ped; replace unknowns N with 0
    for PREFIX in BEAGLE IMPUTE2
    do
		## -i in-place replacement
	#    sed -i 's/N/0/g' out_GTOOL/$PREFIX.$CHROMOSOME.ped
		if [ -s out_GTOOL/$sepjoin/$PREFIX.$CHROMOSOME.map ] && [ -s out_GTOOL/$sepjoin/$PREFIX.$CHROMOSOME.ped ] && [ ! -s out_GTOOL/$sepjoin/$PREFIX.$CHROMOSOME.ped2 ]
		then
		echo "sed (ped > ped2) replace N with 0" $PREFIX
	#    sed -i 's/N/0/g' out_GTOOL/$PREFIX.$CHROMOSOME.ped
		sed 's/N/0/g' out_GTOOL/$sepjoin/$PREFIX.$CHROMOSOME.ped > out_GTOOL/$sepjoin/$PREFIX.$CHROMOSOME.ped2
		fi
    done

    #
    # PLINK merge
    #
    # 3) plink merge
    for PREFIX in BEAGLE IMPUTE2 ApplyRecalibration
    do

        if [ "$PREFIX" == "BEAGLE" ] || [ "$PREFIX" == "IMPUTE2" ]
        then
		dn="out_GTOOL"
		extension="ped2"
        else
		dn="out_VCFtools"
		extension="ped"
        fi
        
        prefixplink=plink_merge_${PREFIX}_${CHROMOSOME}
        out=out_plink/$sepjoin/$prefixplink
		in_ped=$dn/$sepjoin/$PREFIX.$CHROMOSOME.$extension
		in_map=$dn/$sepjoin/$PREFIX.$CHROMOSOME.map2
        if [ ! -s $out.diff ] && [ -s $in_ped ]
        then
			echo "plink merge" $PREFIX $CHROMOSOME
            cmd="bsub -q normal -M4000000 -R'select[mem>4000] rusage[mem=4000]' -P uganda \
			-J\"plink$PREFIX$CHROMOSOME\" \
			-o stdout/$prefixplink.out \
			-e stderr/$prefixplink.err \
            plink --file $prefixomni \
                  --merge $in_ped $in_map \
                  --merge-mode 7 \
                  --out $out \
                  --noweb"
			echo $cmd
			eval $cmd

        fi

    done

done

echo "plink merge finished"

# ##
# ## 1) create file of individuals to keep from bed
# ##

# more +3 chromosome$CHROMOSOME.sample | awk '{print $1" "$2}' > keep.sample

#awk '{print $1" "$2}' \
#backup/omni2.5-8_20120516_gwa_ugand_gtu_autosomes_postsampqc_postsnpqc_flipped.ped \
#| grep -f keep_short_single_column.txt > keep_long_double_column.txt

##
## 2) create file of SNPs/rsIDs to exclude from bed
##

# ##
# ## 2a) find duplicates (markers) in map file (*AFTER* renaming otherwise not duplicate lines, make sure lines are sorted)
# ##
# #awk -F '\t' '{print $1":"$4}' backup/omni2.5-8_20120516_gwa_ugand_gtu_autosomes_postsampqc_postsnpqc_flipped.map | uniq -d \
# #    > markers.duplicates1
# awk -F '\t' '{print $4}' backup/omni2.5-8_20120516_gwa_ugand_gtu_autosomes_postsampqc_postsnpqc_flipped.map | uniq -D \
    # > exclude_duplicates.chrpos
# ## convert positions to original markers (use pipe and grep -v instead to avoid file i/o)
# cat backup/omni2.5-8_20120516_gwa_ugand_gtu_autosomes_postsampqc_postsnpqc_flipped.map | grep -f exclude_duplicates.chrpos | awk '{if( NR%2 ){print $2}}' \
    # > exclude_duplicates.rsid

##
## 2b) exclude missing snps (create exclude_missnp.rsid)
##
echo "Exclude missing SNPs"

#awk -F "\t" '{print $2}' out_plink/plink*ApplyRecalibration_merge.missnp >> exclude_1col_chrpos.missnp
echo "parse missnp"
cmd="awk -F \"\\t\" '{print \$2}' out_plink/\$sepjoin/plink_merge_*.missnp | awk -F : '{print \$2}' > exclude_1col_pos.missnp"
echo $cmd
eval $cmd

# risk of excluding identical positions in other chromosomes...
#cat backup/omni2.5-8_20120516_gwa_ugand_gtu_autosomes_postsampqc_postsnpqc_flipped.map | grep -f exclude_1col_chrpos.missnp \
# cat backup/omni2.5-8_20120516_gwa_ugand_gtu_autosomes_postsampqc_postsnpqc_flipped.map | grep -f exclude_1col_pos.missnp \
    # | awk -F "\t" '{print $2}' \
    # > exclude_missnp.rsid
echo "get rsids of missnp pos"
awk 'FNR==NR{a[$0];next}($4 in a){print $2}' \
exclude_1col_pos.missnp \
backup/omni2.5-8_20120516_gwa_ugand_gtu_autosomes_postsampqc_postsnpqc_flipped.map \
> exclude_missnp.rsid
rm exclude_1col_pos.missnp

## merge exclude files and append to existing file to avoid previously missing SNPs not being included
echo "Merge exclusion files"
cat exclude_duplicates.rsid exclude_missnp.rsid >> exclude.rsid
rm exclude_missnp.rsid

## exclude duplicates for aesthetic reasons
sort -u exclude.rsid > exclude.rsid2
mv exclude.rsid2 exclude.rsid