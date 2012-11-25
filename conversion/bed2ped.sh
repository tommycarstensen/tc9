bsub -q normal -M3000000 -R'select[mem>3000] rusage[mem=3000]' -P agv -o plink.out -e plink.err plink --bfile omni2.5-8_20120516_gwa_ugand_gtu --recode --out omni2.5-8_20120516_gwa_ugand_gtu
