#!/bin/sh



# example script of how to convert bed to ped
## this will also create a map file, right?
http://pngu.mgh.harvard.edu/~purcell/plink/dataman.shtml#recode
bsub -q normal -M4000000 -R'select[mem>4000] rusage[mem=4000]' -P agv -o plink_recode.out -e plink_recode.err \
     plink \
     --bfile /lustre/scratch107/projects/agv/phasing_rel/shapeit2/data/omni2.5-8_20120516_gwa_ugand_gtu_autosomes_postsampqc_postsnpqc_flipped \
     --recode \
     --out omni2.5-8_20120516_gwa_ugand_gtu_autosomes_postsampqc_postsnpqc_flipped \
     --exclude omni2.5-8_20120516_gwa_ugand_gtu_autosomes_postsampqc_postsnpqc_flipped.duplicates \
     --keep keep.sample \

## find duplicates
awk -F '\t' '{print $1":"$4}' omni2.5-8_20120516_gwa_ugand_gtu_autosomes_postsampqc_postsnpqc_flipped.map.rsids | uniq -d > omni2.5-8_20120516_gwa_ugand_gtu_autosomes_postsampqc_postsnpqc_flipped.duplicates
