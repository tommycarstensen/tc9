#!/software/bin/python

## T. Carstensen (tc9), M.S. Sandhu (ms23), D. Gurdasani (dg11)
## Wellcome Trust Sanger Institute, 2012

import os

def main():

    dbSNP = '/lustre/scratch107/projects/uganda/users/tc9/in_GATK/dbsnp_135.b37.vcf'

    for chromosome in range(1,22+1,)+['X','Y',]:

        for prefix in [
            'CombineVariants',
            'ApplyRecalibration.recalibrated.filtered',
            ]:
            fp_out = 'out_GATK/VariantEval_%s%s' %(
                prefix,chromosome,
                )
            if os.path.isfile(fp_out):
                continue
            fp_in = 'out_GATK/sep/%s.%s.vcf' %(prefix,chromosome,)
            if not os.path.isfile(fp_in):
                print fp_in
                stop
                continue
            
            s = 'bsub \
            -o tmp.out -e tmp.err \
            -M16000000 -R\'select[mem>16000] rusage[mem=16000]\' \
            java -Xmx1g \
            -jar /software/varinf/releases/GATK/GenomeAnalysisTK-1.4-15-gcd43f01/GenomeAnalysisTK.jar \
            -T VariantEval \
            --eval %s \
            -R /lustre/scratch111/resources/vrpipe/ref/Homo_sapiens/1000Genomes/human_g1k_v37.fasta \
            --dbsnp %s \
            -o %s \
            ' %(fp_in, dbSNP, fp_out,)
            print s
    ##        stop
            os.system(s)
##            stop

    return

if __name__ == '__main__':
    main()
