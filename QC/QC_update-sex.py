#!/software/bin/python

## T. Carstensen (tc9), M.S. Sandhu (ms23), D. Gurdasani (dg11)
## Wellcome Trust Sanger Institute, 2012

import sys, os
sys.path.append('/nfs/users/nfs_t/tc9/lib/python2.7/site-packages')
import xlrd

def main():

    ##
    ## update-sex.txt
    ##
    workbook = '/lustre/scratch107/projects/agv/data/1000genomes_all2141_150812.xlsx'
    wb = xlrd.open_workbook(workbook)
    sh = wb.sheet_by_index(0)
    l_header = sh.row_values(0)
    col_sex = l_header.index(u'Gender')
    col_iid = l_header.index(u'Broad_ID')
    fd = open('update-sex.txt','w')
    for rownum in range(1,sh.nrows):
        row = sh.row_values(rownum)
        sex = row[col_sex]
        iid = row[col_iid]

        ## convert float to string
        if type(iid) == float:
            iid = str(int(iid))
            pass

        fid = iid
        if sex == 'male':
            sex = 1
        elif sex == 'female':
            sex = 2
        else:
            print sex
            stop
        fd.write('%s\t%s\t%s\n' %(fid,iid,sex,))
    fd.close()

    ##
    ## --update-sex
    ## http://pngu.mgh.harvard.edu/~purcell/plink/dataman.shtml#updatefam
    ##
    bfile = 'omni2.5-4_20120904_agv_gtu_aaa'
    cmd = 'plink --bfile %s \\\n' %(bfile)
    cmd += '--noweb --nonfounders --allow-no-sex --make-bed \\\n'
    cmd += '--out %s_sexupdated \\\n' %(bfile)
    cmd += '--update-sex update-sex.txt'
    
    bsub = "bsub -M4000000 -R'select[mem>4000] rusage[mem=4000]' "
    bsub += '-P agv '
    bsub += '-q normal '
    bsub += "-J'%s' " %('update-sex')
    bsub += '%s' %(cmd)

    os.system(bsub)

    return

if __name__ == '__main__':
    main()
