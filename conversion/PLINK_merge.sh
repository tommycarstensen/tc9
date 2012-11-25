## http://pngu.mgh.harvard.edu/~purcell/plink/dataman.shtml#merge
## merge-mode 7 Report mismatching *non-missing* calls (diff mode -- do not merge)
## check concorande between two ped files

## ped
bsub -q normal -M4000000 -R'select[mem>4000] rusage[mem=4000]' -J"ped" -P uganda -o LSF_plink22_merge_ped.out -e LSF_plink22_merge_ped.err \
plink --file omni2.5-8_20120516_gwa_ugand_gtu_autosomes_postsampqc_postsnpqc_flipped \
      --merge gtool.22.ped gtool.22.map \
      --merge-mode 7 \
      --out plink22_merge_ped \
      --noweb \
      --keep keep.sample

## bed
bsub -q normal -M4000000 -R'select[mem>4000] rusage[mem=4000]' -J"bed" -P uganda -o LSF_plink22_merge_bed.out -e LSF_plink22_merge_bed.err \
plink --bfile /lustre/scratch107/projects/agv/phasing_rel/shapeit2/data/omni2.5-8_20120516_gwa_ugand_gtu_autosomes_postsampqc_postsnpqc_flipped \
      --merge gtool.22.ped gtool.22.map \
      --merge-mode 7 \
      --out plink22_merge_bed \
      --noweb \
      --keep keep.sample