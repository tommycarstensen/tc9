## Tommy Carstensen, Wellcome Trust Sanger Institute, June 2012
## GNU General Public License

import sys
sys.path.append('C:\Users\Tommy Carstensen\Dropbox\python\misc')
import gnuplot
sys.path.append('C:\Users\Tommy Carstensen\Dropbox\python\math')
import statistics

def main():

    d_pedigree = read_pedigree()

    d_gwas_id = read_households()

    compare_pedigree_and_IBD(d_pedigree, d_gwas_id,)

    return


def read_households():

    '''read raw data and parse information about households'''

    print 'reading households'

    fd = open('in/reldata_allgpc_90212.tsv','r')
    lines = fd.readlines()
    fd.close()

    l_headers = lines[0].strip().split('\t')

    d_vhno = {}
##    d_idno = {}
    d_gwas_id = {}
    for i_line in range(1, len(lines)):
        line = lines[i_line]
        l = line.strip().split('\t')
        ## no known relationship and not a member of a household
        if len(l) == 2:
            continue
        idno = int(l[l_headers.index('idno')])
        gwas_id = l[l_headers.index('gwas_id')][1:-1]
        vhno = int(l[l_headers.index('vhno')])
        rel1 = l[l_headers.index('rel1')]
##        d_idno[idno] = {'gwas_id':gwas_id,'vhno':vhno,}
        d_gwas_id[gwas_id] = {
            'idno':idno, 'vhno':vhno, 'rel1':rel1,
            }

    return d_gwas_id


def read_pedigree():

    ## read output from pedigree.py
    fd = open('out/pedigree.tsv','r')
    lines = fd.readlines()
    fd.close()

    ## create a dictionary with predicted PI_HAT values
    d_pedigree = {}
    for i_line in range(len(lines)):

        line = lines[i_line]
        
        l = line.strip().split('\t')

        IID = l[1].strip()
        if IID == '':
            continue

        d_pedigree = pedigree_append_parents(
            d_pedigree,IID,l,
            )

        d_pedigree = pedigree_append_siblings(
            d_pedigree,IID,l,
            )

    return d_pedigree


def pedigree_append_siblings(d_pedigree,IID,l):

    for i_sibling in range(8,len(l),3):
        IID_sibling = l[i_sibling+1].strip()
        ambiguity_sibling = int(l[i_sibling+2])
        if not IID_sibling in d_pedigree.keys():
            d_pedigree[IID_sibling] = {}
        d_pedigree[IID_sibling][IID] = {
            'PI_HAT':0.5,'ambiguity':ambiguity_sibling,
            'relationship':'sibling',
            }
        d_pedigree[IID][IID_sibling] = {
            'PI_HAT':0.5,'ambiguity':ambiguity_sibling,
            'relationship':'sibling',
            }

    return d_pedigree


def pedigree_append_parents(d_pedigree,IID,l,):

    IID_father = l[3].strip()
    ambiguity_father = int(l[4])
    IID_mother = l[6].strip()
    ambiguity_mother = int(l[7])
    for ID in [IID,IID_father,IID_mother,]:
        if not ID in d_pedigree.keys():
            d_pedigree[ID] = {}

    d_pedigree[IID][IID_father] = {
        'PI_HAT':0.5,'ambiguity':ambiguity_father,'relationship':'child',
        }
    d_pedigree[IID][IID_mother] = {
        'PI_HAT':0.5,'ambiguity':ambiguity_mother,'relationship':'child',
        }
    d_pedigree[IID_father][IID] = {
        'PI_HAT':0.5,'ambiguity':ambiguity_father,'relationship':'child',
        }
    d_pedigree[IID_mother][IID] = {
        'PI_HAT':0.5,'ambiguity':ambiguity_mother,'relationship':'child',
        }

    return d_pedigree


def compare_pedigree_and_IBD(d_pedigree,d_gwas_id,):

    ## OMFG this function needs to be rewritten and split into sub routines
    ## it's too bloody long...

    print 'comparing genealogy and genetics'

    d_stats = init_compare_pedigree_and_IBD()
    l_skip = init_skip_gwas_id()
    l_header = None
    l_skip_vhno = init_skip_vhno()

    fd = open('in/exclusion.txt','r')
    s = fd.read()
    fd.close()
    l_exclusion = eval(s)

    d_stats = loop_lines(
        d_stats,
        l_skip,l_skip_vhno,l_exclusion,
        l_header,
        d_pedigree,
        d_gwas_id,
        )

    for relative in d_stats.keys():
        for ambiguity in d_stats[relative].keys():
            print relative, ambiguity, d_stats[relative][ambiguity]

    return


def init_skip_vhno():

    fd = open('out/error_self_relation.txt','r')
    lines = fd.readlines()
    fd.close()
    l_skip_vhno = [int(line.split()[0]) for line in lines[1:]]

    return l_skip_vhno


def loop_lines(
    d_stats,
    l_skip,l_skip_vhno,l_exclusion,
    l_header,
    d_pedigree,
    d_gwas_id,
    ):

    fd = open('in/compare_ibs_pairs_LDpruned.genome','r')
    for line in fd:
        l = line.strip().split()
        if not l_header:
            l_header = l
            continue
        IID1 = l[l_header.index('IID1')]
        IID2 = l[l_header.index('IID2')]
        PI_HAT_calculated = PI_HAT = float(l[l_header.index('PI_HAT')])

        ## temporary print statement
        if len(
            set(['197833_B01_APP5119039', '197827_G04_APP5118956', '197827_B02_APP5118931', '197845_E06_APP5211012', '197834_B09_APP5118209', '197833_A05_APP5119077', '197831_B01_APP5118258', '197833_C05_APP5119079', '197834_C09_APP5118210', '197828_E12_APP5119273', '197828_E03_APP5119185', '197827_H10_APP5119017', '197829_D11_APP5119381', '197845_A10_APP5211046', '197833_A06_APP5119085', '197832_G10_APP5118634', '197831_C01_APP5118259', '197834_A09_APP5118208', '197827_C02_APP5118932', '197835_C11_APP5118504'])
            &
            set([IID1,IID2,])
            ) == 2:
            print IID1,IID2,PI_HAT_calculated

        ## are individuals genealogically related?
        bool_related = True
        if not IID1 in d_pedigree.keys():
            bool_related = False
        if bool_related == True:
            if not IID2 in d_pedigree[IID1].keys():
                bool_related = False

        ## check that individuals that are
        ## not genealogically related
        ## are not genetically related
        if bool_related == False:
            check_genealogy_relation_0(
                PI_HAT_calculated, IID1, IID2,
                l_skip, l_skip_vhno, l_exclusion,
                d_gwas_id,
                )
            continue
        ## check that individuals that are
        ## genealogically related
        ## are also genetically related
        else:
            d_stats = check_genealogy_relation_1(
                PI_HAT_calculated, IID1, IID2,
                d_pedigree, d_stats,
                )
    fd.close()

    return d_stats


def init_compare_pedigree_and_IBD():

    d_columns = {
        'out/pedigree_relationship_incorrectly_identified.txt':[
            'PI_HAT_predicted','PI_HAT_calculated',
            'relationship','ambiguity',
            'IID1','IID2',
            ],
        'out/pedigree_relationship_not_identified.txt':[
            'PI_HAT_calculated',
            'IID1','IID2',
            'same_household',
            'rel1_IID1','rel1_IID2',
            ],
        }
    d_columns['out/pedigree_relationship_incorrectly_identified_0.txt'] = d_columns['out/pedigree_relationship_incorrectly_identified.txt']
    d_columns['out/pedigree_relationship_incorrectly_identified_1.txt'] = d_columns['out/pedigree_relationship_incorrectly_identified.txt']
    del d_columns['out/pedigree_relationship_incorrectly_identified.txt']
    for fn in d_columns.keys():
        l_columns = d_columns[fn]
        fd = open(fn,'w')
        fd.write('##%s\n' %('\t'.join(l_columns)))
        fd.close()

    d_stats = {
        'child':{
            'ambiguity1':{'outlier1':0,'outlier0':0,},
            'ambiguity0':{'outlier1':0,'outlier0':0,},
            },
        'sibling':{
            'ambiguity1':{'outlier1':0,'outlier0':0,},
            'ambiguity0':{'outlier1':0,'outlier0':0,},
            },
        }

    return d_stats


def check_genealogy_relation_1(
    PI_HAT_calculated, IID1, IID2, d_pedigree,
    d_stats,
    ):

    '''check'''

    PI_HAT_predicted = d_pedigree[IID1][IID2]['PI_HAT']
    ambiguity = d_pedigree[IID1][IID2]['ambiguity']
    relationship = d_pedigree[IID1][IID2]['relationship']
    if abs(PI_HAT_calculated-PI_HAT_predicted) > 0.1:
##            print PI_HAT_predicted, PI_HAT_calculated, IID1, IID2
        fp = 'out/pedigree_relationship_incorrectly_identified_%i.txt' %(ambiguity)
        fd2 = open(fp,'a')
        fd2.write(
            '%.2f %.2f %7s %1i %21s %21s\n' %(
                PI_HAT_predicted, PI_HAT_calculated,
                relationship, ambiguity,
                IID1, IID2,
                )
            )
        fd2.close()
        d_stats[relationship]['ambiguity%i' %(ambiguity)]['outlier1'] += 1
    else:
        d_stats[relationship]['ambiguity%i' %(ambiguity)]['outlier0'] += 1

    return d_stats


def check_genealogy_relation_0(
    PI_HAT_calculated, IID1, IID2,
    l_skip, l_skip_vhno, l_exclusion,
    d_gwas_id,
    ):

    '''check genetic relation, when no genealogical relation'''

    if PI_HAT_calculated > 0.40:

        ## skip individuals not assigned to any households
        if IID1 in l_skip or IID2 in l_skip:
            pass
        else:
            ## only append to list if from same household
            vhno1 = d_gwas_id[IID1]['vhno']
            vhno2 = d_gwas_id[IID2]['vhno']
            if vhno1 != vhno2:
                return
            ## was this household skipped when doing the pedigree analysis
            ## and should it thus be skipped now?
            if vhno1 in l_skip_vhno or vhno2 in l_skip_vhno:
                return
            ## only append to list if relationship given
            if d_gwas_id[IID1]['rel1'] == '' or d_gwas_id[IID2]['rel1'] == '':
                return
            ## only append to list if individual was not excluded
            idno1 = d_gwas_id[IID1]['idno']
            idno2 = d_gwas_id[IID2]['idno']
            if idno1 in l_exclusion or idno2 in l_exclusion:
                return
            ## only append to list if relationship data is available
            if vhno1 == vhno2:
                fd2 = open('out/pedigree_relationship_not_identified.txt','a')
                fd2.write(
                    '%.2f %21s %21s %5s %1s %1s\n' %(
                        PI_HAT_calculated, IID1, IID2, vhno1==vhno2,
                        d_gwas_id[IID1]['rel1'], d_gwas_id[IID2]['rel1'],
                        )
                    )
                fd2.close()

    return


def init_skip_gwas_id():

    '''gwas_ids not assigned to any households'''

    s_skip = '''
197835_A01_APP5118386
197828_D10_APP5119251
197816_D04_APP5117762
197845_B06_APP5118936
197816_D09_APP5117815
197845_G11_APP5211060
197831_H09_APP5118348
197843_H07_APP5118719
197835_D11_APP5118505
'''

    l_skip = s_skip.split()

    return l_skip


def junk():

    ##                ## rare cases (1 occurence)
    ##                '4_3', ## parent of child (20037013 / 20037057) ## redundant
    ##                '4_4', ## parent of parent (17679082)
    ##                '5_5', ## sibling of sibling (7052144)
    ##                '5_8', ## sibling of grandchild (20125188)
    ##                '8_5', ## grandchild of sibling
    ##                '5_4', ## sibling of parent (12919077) ## 43 year difference???
    ##                ## rare cases (less than 10 occurences)
    ##                '3_4', ## child of parent (20037057 / 20037013)
    ##                '3_8', ## grandchild of child (17112122)
    ##                '2_8', ## spouse of grandchild (11761051 / 11761040)
    ##                '2_4', ## spouse of parent (20037126 / 20037013)
    ##                ## ask Kenneth/Deepti
    ##                '8_3', ## grandchild of child
    ##                '8_4', ## grandchild of parent (17270101 / 17270145)

    return

if __name__ == '__main__':
    main()
