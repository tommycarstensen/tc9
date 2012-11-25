#!/software/bin/python

## T. Carstensen (tc9), M.S. Sandhu (ms23), D. Gurdasani (dg11)
## Wellcome Trust Sanger Institute, 2012

import sys, os
sys.path.append('/nfs/users/nfs_t/tc9/lib/python2.7/site-packages')
import xlrd

def main():

    if not os.path.isdir('samples'):
        os.mkdir('samples')

    l_xlsx = [
        '/lustre/scratch107/projects/agv/data/AGV_Omni-quad_720samples-to-recall_final_150812.xlsx',
        '/lustre/scratch107/projects/agv/data/AGV_Omni-8_1465samples_to-recall_final_150812.xlsx',
        '/lustre/scratch107/projects/agv/data/1000genomes_all2141_150812.xlsx',
        'GPC_8072_trb_IDlinks_TC_021012.xlsx',
##        '1000G_pops_acronyms.xlsx',
        ]

    d_group2sample = {}
##        lines_keep = []
##        lines_remove = []
    for workbook in l_xlsx:
        print 'workbook', workbook
        if not os.path.isfile(workbook):
            continue
        wb = xlrd.open_workbook(workbook)
        sh = wb.sheet_by_index(0)
        l_header = sh.row_values(0)
        if workbook == 'GPC_8072_trb_IDlinks_TC_021012.xlsx':
            col_pop = l_header.index(u'trb')
            col_sample = l_header.index(u'sangerID')
        else:
            try:
                col_pop = l_header.index(u'Population') ## 1000g
                col_sample = l_header.index(u'Coriell_Sample_ID') ## 1000g
            except:
                col_pop = l_header.index(u'Ethnicity') ## agv
                col_sample = l_header.index(u'Sanger Sample Name') ## agv
        for rownum in range(1,sh.nrows):
            row = sh.row_values(rownum)
            sample = row[col_sample]
            group = row[col_pop].replace(' ','')
            ## convert float to int
            if type(sample) == float:
                sample = int(sample)
            ## convert to str
            sample = str(sample)
            ## N.B. ERROR!!!
            if sample == 'HG00096':
                group = 'GBR'
            if str(group) == '':
                if workbook == 'GPC_8072_trb_IDlinks_TC_021012.xlsx':
                    group = 'Unknown'
                else:
                    print row
                    stop
            if not group in d_group2sample.keys():
                d_group2sample[group] = []
            if not sample in d_group2sample[group]:
                d_group2sample[group] += [sample]

    for group in d_group2sample.keys():
        fp_out = 'samples/%s.samples' %(group.replace(' ',''))
        fd = open(fp_out,'w')
        fd.write('\n'.join(d_group2sample[group]))
        fd.close()
        print fp_out
##                if workbook == '1000genomes_all2141_150812.xlsx':
##                    lines_remove += ['%s %s\n' %(sample,sample,)]
##
####        fd = open('keep.txt','w')
####        fd.writelines(lines_keep)
####        fd.close()
##
####        fd = open('remove.txt','w')
####        fd.writelines(lines_remove)
####        fd.close()
##
##        fd = open('d_sample2group.txt','w')
##        fd.write(str(d_sample2group))
##        fd.close()
##
####        print 'reading and evaluating d_sample2group.txt'
####        fd = open('d_sample2group.txt','r')
####        s = fd.read()
####        fd.close()
####        d_sample2group = eval(s)

    return

if __name__ == '__main__':
    main()
