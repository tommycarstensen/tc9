#!/software/bin/python

## T. Carstensen (tc9), M.S. Sandhu (ms23), D. Gurdasani (dg11)
## Wellcome Trust Sanger Institute, 2012

import os, sys

def main():

    ## split chips for each population

    l_populations = [u'Muganda', u'KIKUYU', u'Mandinka', u'ZULU', u'IBO', u'KALENJIN', u'Fula', u'Murundi', u'Munyarwanda', u'AMHARA', u'Sotho', u'Wolloff', u'SOMALI', u'Jola', u'OROMO', u'GA-ADANGBE']

    d_io = {
        'omni2.5-8_agv_20120910_gtu':'omni2.5-8_agv_20120910_gtu',
        'omni2.5-4_20120904_agv_gtu_aaa_sexupdated':'omni2.5-4_20120904_agv_gtu_aaa',
        }

    for bfile in d_io.keys():
##    for bfile in [
##        'omni2.5-8_agv_20120910_gtu',
##        'omni2.5-4_20120904_agv_gtu_aaa',
##        ]:
        for population in l_populations:

            out = '%s_%s' %(d_io[bfile],population,)
            if os.path.isfile('%s.bed' %(out)):
                continue

            samples = 'samples/%s.samples' %(population)
            keep = '%s_%s.keep' %(bfile,population,)

            ## I would prefer a guaranteed one to one relatioship...
            ## but genotyping IDs missing in the Excel sheet...
            cmd = "grep -f %s %s.fam | awk '{print $1,$2}' > %s" %(samples,bfile,keep,)
            os.system(cmd)
            count = int(os.popen('cat %s | wc -l' %(keep)).read())
            ## no samples in bfile
            if count == 0:
                os.remove(keep)
                print 'skip', bfile, population
                continue

            cmd = 'plink --bfile %s --keep %s --noweb --nonfounders --make-bed --out %s' %(
                bfile,keep,out,
                )

            bsub = "bsub -M4000000 -R'select[mem>4000] rusage[mem=4000]' "
            bsub += '-P agv '
            bsub += '-q normal '
            bsub += "-J'%s' " %(population)
            bsub += '%s' %(cmd)

            os.system(bsub)

            continue

        continue

    return


if __name__ == '__main__':
    main()
