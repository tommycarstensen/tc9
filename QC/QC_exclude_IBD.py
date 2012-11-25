#!/software/bin/python

## T. Carstensen (tc9), M.S. Sandhu (ms23), D. Gurdasani (dg11)
## Wellcome Trust Sanger Institute, 2012

import sys, os
sys.path.append('/nfs/users/nfs_t/tc9/lib/python2.7/site-packages')
import xlrd

def main(fp_genome,IBD_min=.05,):

    fp_out = '%s.filtered' %(fp_genome)
    if not os.path.isfile(fp_out):
        cmd = "awk '{if ($10>0%.2f) print}' %s > %s" %(IBD_min,fp_genome,fp_out,)
        print cmd
        os.system(cmd)

    fd = open(fp_out,'r')
    lines = fd.readlines()
    fd.close()
    d_count = {1:[],}
    ## populate dictionary
    for line in lines:
        l = line.split()
        FID1 = l[0]
        d[FID1] = 0
    for line in lines:
        l = line.split()
        FID1 = l[0]
        FID2 = l[2]
        d[FID1] += 1
        d[FID2] += 1
        if d[FID1] == 1:
            d_count[1] += [FID1]
        else:
            d_count[d[FID1]-1] += [FID1]
        if d[FID2] == 1:
            d_count[1] += [FID2]
        stop

    return

if __name__ == '__main__':
    main('omni2.5-8_20120809_gwa_uganda_gtu.genome')
    main()
