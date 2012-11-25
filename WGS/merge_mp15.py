#!/software/bin/python

## T. Carstensen (tc9), M.S. Sandhu (ms23), D. Gurdasani (dg11)
## Wellcome Trust Sanger Institute, 2012

import os

def main():

    dn_in = '/lustre/scratch111/projects/uganda/vcf_intermediate-vqsr-2/pooled'

    fn_out = 'mp15_vqsr.vcf'

    s = 'bsub \
    -e %s.err \
    -M4000000 -R\'select[mem>4000] rusage[mem=4000]\' \
    ' %('tmp',)
    s += 'zcat '
    for chromosome in range(1,22+1,)+['X','Y',]:
        fn_in = '%s.vqsr.filt.vcf.gz' %(chromosome)
        fp_in = os.path.join(dn_in,fn_in,)
        s += ' %s ' %(fp_in)
    s += ' > %s' %(fn_out)
    print s
    stop
    os.system(s)

    return

if __name__ == '__main__':
    main()
