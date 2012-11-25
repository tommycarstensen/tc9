#!/software/bin/python

## T. Carstensen (tc9), M.S. Sandhu (ms23), D. Gurdasani (dg11)
## Wellcome Trust Sanger Institute, 2012

## THIS SCRIPT CURRENTLY DOESN'T WORK!!!

import os

def main():

    dbSNP = '/lustre/scratch107/projects/uganda/users/tc9/in_GATK/dbsnp_135.b37.vcf'
    fp_in = '/lustre/scratch107/projects/agv/phasing_rel/shapeit2/data/omni2.5-8_20120516_gwa_ugand_gtu_autosomes_postsampqc_postsnpqc_flipped.bed'
    fn_out = 'omni2.5-8_20120516_gwa_ugand_gtu_autosomes_postsampqc_postsnpqc_flipped.vcf'
    walker = 'VariantsToVCF'
    R = '/lustre/scratch111/resources/vrpipe/ref/Homo_sapiens/1000Genomes/human_g1k_v37.fasta'

    s = 'bsub \
    -e %s.err \
    -M4000000 -R\'select[mem>4000] rusage[mem=4000]\' \
    -J %s \
    ' %(walker,walker,)
    s += ' java -Xmx4g -jar /software/varinf/releases/GATK/GenomeAnalysisTK-1.4-15-gcd43f01/GenomeAnalysisTK.jar \
   -R %s \
   -T VariantsToVCF \
   -o %s \
   --variant:BED %s \
   --dbsnp %s \
   ' %(R,fn_out,fp_in,dbSNP,)
    print s
##    stop
    os.system(s)

    return

if __name__ == '__main__':
    main()
