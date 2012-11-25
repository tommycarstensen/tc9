#!/software/bin/python

## T. Carstensen (tc9), M.S. Sandhu (ms23), D. Gurdasani (dg11)
## Wellcome Trust Sanger Institute, 2012

## http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_walkers_varianteval_VariantEval.html

import os, time

def main():

    dbSNP = '/lustre/scratch107/projects/uganda/users/tc9/in_GATK/dbsnp_135.b37.vcf'
##    dbSNP = 'omni2.5-8_20120516_gwa_ugand_gtu_autosomes_postsampqc_postsnpqc_flipped.vcf'
    walker = 'VariantEval'
##    comp = omnivcf

##    import subprocess
####    'cat tmp.out'
####    fd = subprocess.Popen(['ls','-l',])
##    fd = subprocess.Popen(['cat','tmp.out',],stdout=subprocess.PIPE,)
##    while True:
##        line = fd.stdout.readline()
##        if line == '':
##            break
##        print line
##    stop

    if not os.path.isfile(dbSNP):
        print 'dbSNP missing', dbSNP
        return

    for fp_in in [
        'out_GATK/join/CombineVariants.vcf',
        'out_GATK/join/ApplyRecalibration.recalibrated.filtered.vcf',
        'SelectVariants_discordance1.vcf',
        'SelectVariants_discordance2.vcf',
        'SelectVariants_concordance.vcf',
##        'mp15_vqsr.vcf',
        'out_mp15/vqsr/vqsr_combinevariants.vcf',
        ]:
        if '/' in fp_in:
            index1 = fp_in.index('/')+1
        else:
            index1 = 0
        index2 = fp_in.index('.')
        fp_out = 'out_GATK/VariantEval_%s.txt' %(
##        fp_out = 'out_GATK/VariantEval_%s_omni.txt' %(
            fp_in[index1:index2],
            )
        if os.path.isfile(fp_out):
            print 'output exists', fp_out
            continue

##        fp_in = 'out_GATK/%s' %(fp_in)

        if not os.path.isfile(fp_in):
            print 'input does not exist', fp_in
            continue
        ## input file recently generated and possibly still being written to
        if time.time()-os.path.getmtime(fp_in) < 10*60:
            print 'input is a new file', fp_in
            continue
        
        s = 'bsub \
        -o %s.out -e %s.err \
        -M4000000 -R\'select[mem>4000] rusage[mem=4000]\' \
        -J %s \
        java -Xmx4g \
        -jar /software/varinf/releases/GATK/GenomeAnalysisTK-1.4-15-gcd43f01/GenomeAnalysisTK.jar \
        -T VariantEval \
        --eval %s \
        -R /lustre/scratch111/resources/vrpipe/ref/Homo_sapiens/1000Genomes/human_g1k_v37.fasta \
        --dbsnp %s \
        ' %(walker,walker,walker,fp_in, dbSNP,)
##        --comp %s \ %(comp)
        s += '\
        -o %s \
        ' %(fp_out,)
##        print s
##        stop
        os.system(s)
##        stop

    return

if __name__ == '__main__':
    main()
