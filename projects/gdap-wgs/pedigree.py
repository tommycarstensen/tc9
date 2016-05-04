## Tommy Carstensen
## Genome Research Ltd
## May-June 2012, Dec2015-May2016
## GNU General Public License

import os
import argparse
import gzip

skipped_households = (
##    '1100358', '200050', '1700018', '700287', '1700434',
##    '2000584', '', '600021', '100090', '1000509', '600403',
##    '300859', '1900803', '1700474', '1900003', '400072',
##    '2400004', '700284', '100141', '2300483', '2300952',
##    '100851', '100153', '100852', '2500094', '100028',
##    '600955', '100071', '300131', '1400092', '600059', '1700137',
##    '1800453', '2100485', '200011', '2500481', '1800486', '1600555',
##    '300132', '800015', '500117', '500006', '1500462', '2400516',
    )

def main():

    args = parse_args()

##    d_idno[idno] = {
##        'relatives':{
##            ## own generation
##            'spouses':[],
##            'siblings':[],
##            ## next and future generations
##            'children':[],
##            'grandchildren':[],
##            ## previous and past generations
##            'parents':[],
##            'grandparents':[],
##            },
##        'sex':sex,
##        'agec':agec,
##        'num':num,
##        'gwas_id':gwas_id,
##        'rel1':{rel1:[[num1,0,]],},
##        'num1':num1,
##        'vhno':vhno,
##        }

    with open('duplicate_num.txt', 'w') as f:
        print('idno', 'agec', 'sex', 'num', 'relh', 'num1', 'rel1', 'vhno', 'sangerID', file=f, sep='\t')

    d_idno = {}
    d_vhno = {}
    d_sangerID2idno = {}
    d_idno2sangerID = {}
    with open(args.tsv) as f:

        for line in f:
            header = line.rstrip().split()
            break

        for line in f:
            l = line.rstrip('\n').split('\t')
            idno = l[header.index('idno')]
            sangerID = l[header.index('sangerID')]
            vhno = l[header.index('vhno')]

            d = {}
            for k, v in zip(header, l):
                d[k] = v

            if vhno in skipped_households:
                continue

            if d['num'] == '' and d['relh'] == '' and d['num1'] == '' and d['rel1'] == '':
                continue

            ## Can't skip samples without sangerID,
            ## because that will break pedigrees.
##            if sangerID == '':
##                continue

            if idno in ('19287029', '8009050', '1059058',):
                continue
            
            d_idno[idno] = {
                'relatives':{
                    ## own generation
                    'spouses':set(),
                    'siblings':set(),
                    ## next and future generations
                    'children':set(),
                    'grandchildren':set(),
                    ## previous and past generations
                    'parents':set(),
                    'grandparents':set(),
                    },
                }
            assert len(header) == len(l)
            for k, v in zip(header, l):
                d_idno[idno][k] = v

            d_idno = manual_correction(d_idno, idno)

            ## Add to d_vhno *after* manual correction.
            vhno = d_idno[idno]['vhno']
            num = d_idno[idno]['num']
            try:
                if num in d_vhno[vhno].keys():
                    with open('duplicate_num.txt', 'a') as f:
                        print(d['idno'], d['agec'], d['sex'], d['num'], d['relh'], d['num1'], d['rel1'], d['vhno'], d['sangerID'], file=f, sep='\t')
                        d = d_idno[d_vhno[vhno][num]]
                        print(d['idno'], d['agec'], d['sex'], d['num'], d['relh'], d['num1'], d['rel1'], d['vhno'], d['sangerID'], file=f, sep='\t')
                        continue
                d_vhno[vhno][num] = idno
            except KeyError:
                d_vhno[vhno] = {num: idno}

            d_sangerID2idno[d_idno[idno]['sangerID']] = idno
            d_idno2sangerID[idno] = d_idno[idno]['sangerID']

###gwas_id: chipID or blank
##idno: individual ID
##vno: village number
##hno: household number? unique in entire dataset or unique per vno?
##agec: age
##sex: 1 male, 2 female
##num: number within household
###ridno: respondent ID number if respondent is not answering for themselves
##relh: relationship to head of household
##num1: household member number
##rel1: relationship in household
## 1) head of household
## 2) spouse
## 3) child
## 4) parent
## 5) brother/sister
## 6) other
## 7) no relationship
## 8) grandchild
## 9) grandparent?
###vhno: village household number? probably unique concatenation of vno and hno.
##sangerID

## num1: Person related to household members on enumeration list
## rel1: What relationship does the participant have in this household in relation to num1?
## num: Number of individual within the household
## relh: What relationship does the participant have to the head of this household?

    del vhno

    with open('relh_eq_1_rel1_eq_1_num_ne_num1.txt', 'w') as f:
        print('idno', 'agec', 'sex', 'num', 'relh', 'num1', 'rel1', 'vhno', 'sangerID', file=f, sep='\t')
    with open('num1_not_in_household.txt', 'w') as f:
        print('idno', 'agec', 'sex', 'num', 'relh', 'num1', 'rel1', 'vhno', 'sangerID', file=f, sep='\t')

    with open('rel1_eq_1_relh_ne_1.txt', 'w') as f:
        print('idno', 'agec', 'sex', 'num', 'relh', 'num1', 'rel1', 'vhno', 'sangerID', file=f, sep='\t')
    with open('relh_eq_1_rel1_ne_1.txt', 'w') as f:
        print('idno', 'agec', 'sex', 'num', 'relh', 'num1', 'rel1', 'vhno', 'sangerID', file=f, sep='\t')

    with open('rel_to_self.txt', 'w') as f:
        print('idno', 'agec', 'sex', 'num', 'relh', 'num1', 'rel1', 'vhno', 'sangerID', file=f, sep='\t')
    with open('child_older_than_parent.txt', 'w') as f:
        print('idno', 'agec', 'sex', 'num', 'relh', 'num1', 'rel1', 'vhno', 'sangerID', file=f, sep='\t')
    with open('grandchild_older_than_grandparent.txt', 'w') as f:
        print('idno', 'agec', 'sex', 'num', 'relh', 'num1', 'rel1', 'vhno', 'sangerID', file=f, sep='\t')
    with open('parent_younger_than_child.txt', 'w') as f:
        print('idno', 'agec', 'sex', 'num', 'relh', 'num1', 'rel1', 'vhno', 'sangerID', file=f, sep='\t')
    with open('spouse_less_than_15.txt', 'w') as f:
        print('idno', 'agec', 'sex', 'num', 'relh', 'num1', 'rel1', 'vhno', 'sangerID', file=f, sep='\t')

    with open('num_num1_or_rel1_missing.txt', 'w') as f:
        print('idno', 'agec', 'sex', 'num', 'relh', 'num1', 'rel1', 'vhno', 'sangerID', file=f, sep='\t')

    with open('multiple_heads.txt', 'w') as f:
        print('idno', 'agec', 'sex', 'num', 'relh', 'num1', 'rel1', 'vhno', 'sangerID', file=f, sep='\t')

    with open('head_less_than_15.txt', 'w') as f:
        print('idno', 'agec', 'sex', 'num', 'relh', 'num1', 'rel1', 'vhno', 'sangerID', file=f, sep='\t')

##    i = 0
##    for idno, d in d_idno.items():
    for idno, d in sorted(d_idno.items()):
##        i += 1
##        print('i', i)

        ## ERROR2
        if d['rel1'] in ('3', '8', '2', '5', '4') and d['num'] == d['num1']:
            with open('rel_to_self.txt', 'a') as f:
                print(d['idno'], d['agec'], d['sex'], d['num'], d['relh'], d['num1'], d['rel1'], d['vhno'], d['sangerID'], file=f, sep='\t')
            continue

        ## ERROR5
        if d['num1'] == '' or d['num'] == '' or d['rel1'] == '':
            with open('num_num1_or_rel1_missing.txt', 'a') as f:
                print(d['idno'], d['agec'], d['sex'], d['num'], d['relh'], d['num1'], d['rel1'], d['vhno'], d['sangerID'], file=f, sep='\t')
            continue

        if d['sangerID'] != '' and d['num1'] == d['num'] and d['relh'] != '1' and d['rel1'] != '1' and d['rel1'] not in ('6', '7'):
            print('vhno', d['vhno'], idno, d['relh'], d['rel1'], d['num1'], d['num'], d)
            stop

##        ## Head of household error.
##        ## ERROR3B
##        if d['rel1'] != '1' and d['relh'] == '1':
##            with open('relh_eq_1_rel1_ne_1.txt', 'a') as f:
##                print(d['idno'], d['agec'], d['sex'], d['num'], d['relh'], d['num1'], d['rel1'], d['vhno'], d['sangerID'], file=f, sep='\t')
##            continue

        ## Head of household.
        if d['relh'] == '1' and d['rel1'] == '1':
            ## ERROR7
            if float(d['agec']) < 15:
                with open('head_less_than_15.txt', 'a') as f:
                    print(d['idno'], d['agec'], d['sex'], d['num'], d['relh'], d['num1'], d['rel1'], d['vhno'], d['sangerID'], file=f, sep='\t')
            ## ERROR3D
            if d['num'] != d['num1']:
                with open('relh_eq_1_rel1_eq_1_num_ne_num1.txt', 'a') as f:
                    print(d['idno'], d['agec'], d['sex'], d['num'], d['relh'], d['num1'], d['rel1'], d['vhno'], d['sangerID'], file=f, sep='\t')
            else:
                cnt = 0
                for idno_i in d_vhno[d['vhno']].values():
                    if d_idno[idno_i]['rel1'] == '1' and d_idno[idno_i]['relh'] == '1':
                        cnt += 1
                ## ERROR3A
                if cnt > 1:
                    with open('multiple_heads.txt', 'a') as f:
                        for idno_i in d_vhno[d['vhno']].values():
                            d = d_idno[idno_i]
                            if d_idno[idno_i]['rel1'] == '1' and d_idno[idno_i]['relh'] == '1':
                                print(d['idno'], d['agec'], d['sex'], d['num'], d['relh'], d['num1'], d['rel1'], d['vhno'], d['sangerID'], file=f, sep='\t')
            continue


        ## idno of the relative.
        try:
            idno_num1 = d_vhno[d['vhno']][d['num1']]
        except KeyError:
            ## ERROR4
            with open('num1_not_in_household.txt', 'a') as f:
                print(d['idno'], d['agec'], d['sex'], d['num'], d['relh'], d['num1'], d['rel1'], d['vhno'], d['sangerID'], file=f, sep='\t')
            continue

        ## 2 Spouse.
        if d['rel1'] == '2':
            d['relatives']['spouses'].add(idno_num1)
            d['relatives'][idno_num1] = 2
            d_idno[idno_num1]['relatives']['spouses'].add(idno)
            d_idno[idno_num1]['relatives'][idno] = 2
            if float(d_idno[idno]['agec']) < 15 or float(d_idno[idno_num1]['agec']) < 15:
                with open('spouse_less_than_15.txt', 'a') as f:
                    print(d['idno'], d['agec'], d['sex'], d['num'], d['relh'], d['num1'], d['rel1'], d['vhno'], d['sangerID'], file=f, sep='\t')
                    ## Dict of spouse.
                    d = d_idno[d_vhno[d['vhno']][d['num1']]]
                    print(d['idno'], d['agec'], d['sex'], d['num'], d['relh'], d['num1'], d['rel1'], d['vhno'], d['sangerID'], file=f, sep='\t')
                    continue
            ## children of the spouse could be children of the other spouse
##            for idno_child in d_idno[idno]['relatives']['children']:
            add_relation(d_idno, idno, idno_num1, 3, 3, 4)
            add_relation(d_idno, idno_num1, idno, 3, 3, 4)
##            for idno_child, _ in d_idno[idno]['relatives'].items():
##                if idno_child == '11093042':
##                    stop66
##                if not _ == 3:
##                    continue
##                if idno_child == '11093042':
##                    stop77
##                if not idno_child in d_idno[idno_num1]['relatives'].keys():
##                    d_idno[idno_num1]['relatives'][idno_child] = 3
##                    d_idno[idno_child]['relatives'][idno_num1] = 4
##                    if idno == '15014094':  # tmp!!!
##                        print('---', idno, idno_num1, idno_child)
            ## grandchildren of the spouse could be grandchildren of the other spouse
            add_relation(d_idno, idno, idno_num1, 8, 8, 9)
            add_relation(d_idno, idno_num1, idno, 8, 8, 9)
####            for idno_grandchild in d_idno[idno]['relatives']['grandchildren']:
##            for idno_grandchild, _ in d_idno[idno]['relatives'].items():
##                if not _ == 8:
##                    continue
##                if not idno_grandchild in d_idno[idno_num1]['relatives'].keys():
##                    d_idno[idno_num1]['relatives'][idno_grandchild] = 8
##                    d_idno[idno_grandchild]['relatives'][idno_num1] = 9
##                    if idno == '15014094':  # tmp!!!
##                        print('###', idno, idno_num1, idno_grandchild)
##            if idno == '15014094':
##                print(d_idno[idno]['relatives'])
##                print(d_idno[idno_num1]['relatives'])
##                stop1

##            if idno == '15014094':
##                print(idno, d_idno[idno]['relatives'])
##                try:
##                    print(11093042, d_idno['11093042']['relatives']['11093122'])
##                except KeyError:
##                    print(11093042, d_idno['11093042']['relatives'])
##                try:
##                    print(11093122, d_idno['11093122']['relatives']['11093042'])
##                except KeyError:
##                    print(11093122, d_idno['11093122']['relatives'])
##                stop2a
##            if idno == '11093042':
##                print(d_idno[idno]['relatives'])
##                print(d_idno[idno_num1]['relatives'])
##                print(d_idno['11093122']['relatives'])
##                stop2b

            continue

        ## 3 Child
        if d['rel1'] == '3':
            ## Child must be younger than parent.
            try:
                assert float(d_idno[idno]['agec']) < float(d_idno[idno_num1]['agec'])
            except AssertionError:
                ## ERROR1A1
                with open('child_older_than_parent.txt', 'a') as f:
                    print(d['idno'], d['agec'], d['sex'], d['num'], d['relh'], d['num1'], d['rel1'], d['vhno'], d['sangerID'], file=f, sep='\t')
                    ## Dict of parent.
                    d = d_idno[d_vhno[d['vhno']][d['num1']]]
                    print(d['idno'], d['agec'], d['sex'], d['num'], d['relh'], d['num1'], d['rel1'], d['vhno'], d['sangerID'], file=f, sep='\t')
                continue
            d['relatives']['parents'].add(idno_num1)
            d_idno[idno_num1]['relatives']['children'].add(idno)
            d['relatives'][idno_num1] = 4
            d_idno[idno_num1]['relatives'][idno] = 3
            ## grandparent of the child could be grandparent of siblings
            ## xxx
            ## spouse of the parent could be parent of the child
            add_relation(d_idno, idno_num1, idno, 2, 4, 3)
            ## child of the parent could be sibling of child
            add_relation(d_idno, idno_num1, idno, 3, 5, 5)
            ## child of the parent (num1) could be parent of grandchild of parent
            add_relation(d_idno, idno_num1, idno, 8, 4, 3)
            continue

        ##
        ## 4 Parent
        ##
        if d['rel1'] == '4':
            ## Parent must be younger than child.
            try:
                assert float(d_idno[idno]['agec']) > float(d_idno[idno_num1]['agec'])
            except AssertionError:
                ## ERROR1A2
                with open('parent_younger_than_child.txt', 'a') as f:
                    print(d['idno'], d['agec'], d['sex'], d['num'], d['relh'], d['num1'], d['rel1'], d['vhno'], d['sangerID'], file=f, sep='\t')
                    ## Dict of child.
                    d = d_idno[d_vhno[d['vhno']][d['num1']]]
                    print(d['idno'], d['agec'], d['sex'], d['num'], d['relh'], d['num1'], d['rel1'], d['vhno'], d['sangerID'], file=f, sep='\t')
                continue
            d['relatives']['children'].add(idno_num1)
            d_idno[idno_num1]['relatives']['parents'].add(idno)
            d['relatives'][idno_num1] = 3
            d_idno[idno_num1]['relatives'][idno] = 4
            ## parent (idno) of the child could be parent (idno) of sibling of child
            add_relation(d_idno, idno_num1, idno, 5, 3, 4)
            continue

        ## 5 Siblings
        if d['rel1'] == '5':
            d['relatives']['siblings'].add(idno_num1)
            d_idno[idno_num1]['relatives']['siblings'].add(idno)
            d['relatives'][idno_num1] = 5
            d_idno[idno_num1]['relatives'][idno] = 5
            ## sibling2 of the sibling1 could be child of parent of sibling1
            add_relation(d_idno, idno_num1, idno, 4, 4, 3)
            ## sibling of the sibling could be sibling of sibling
            add_relation(d_idno, idno_num1, idno, 5, 5, 5)
            add_relation(d_idno, idno, idno_num1, 5, 5, 5)
            continue

        ## 6 Other
        if d['rel1'] == '6':
            d['relatives'][idno_num1] = 6
            d_idno[idno_num1]['relatives'][idno] = 6
            continue

        ## 7 No relationship
        if d['rel1'] == '7':
            continue

        ##
        ## 8 Grandchild
        ##
        if d['rel1'] == '8':
            ## Grandchild must be younger than grandparent.
            try:
                assert float(d_idno[idno]['agec']) < float(d_idno[idno_num1]['agec'])
            except AssertionError:
                ## ERROR1B
                with open('grandchild_older_than_grandparent.txt', 'a') as f:
                    print(d['idno'], d['agec'], d['sex'], d['num'], d['relh'], d['num1'], d['rel1'], d['vhno'], d['sangerID'], file=f, sep='\t')
                    d = d_idno[d_vhno[d['vhno']][d['num1']]]
                    print(d['idno'], d['agec'], d['sex'], d['num'], d['relh'], d['num1'], d['rel1'], d['vhno'], d['sangerID'], file=f, sep='\t')
                continue
            d['relatives']['grandparents'].add(idno_num1)
            d_idno[idno_num1]['relatives']['grandchildren'].add(idno)
            ## grandparent
            d['relatives'][idno_num1] = 9
            ## grandchild
            d_idno[idno_num1]['relatives'][idno] = 8

            ## spouses of the grandparent could be a grandparent of the granchild
##            for idno_spouse in d_idno[idno_num1]['relatives']['spouses']:
            for idno_spouse, _ in d_idno[idno_num1]['relatives'].items():
                if not _ == 2:
                    continue
##                d['relatives']['grandparents'].add(idno_spouse)
##                d_idno[idno_spouse]['relatives']['grandchildren'].add(idno)
                if not idno_spouse in d['relatives'].keys():
                    d['relatives'][idno_spouse] = 9
                    d_idno[idno_spouse]['relatives'][idno] = 8
            ## children of the grandparent could be parents of the grandchild
##            for idno_child in d_idno[idno_num1]['relatives']['children']:
            for idno_child, _ in d_idno[idno_num1]['relatives'].items():
                if not _ == 3:
                    continue
##                d['relatives']['parents'].add(idno_child)
##                d_idno[idno_child]['relatives']['children'].add(idno)
                if not idno_child in d['relatives'].keys():
                    d['relatives'][idno_child] = 4
                    d_idno[idno_child]['relatives'][idno] = 3
            ## parents of the grandchild could be children of the grandparent
##            for idno_parent in d_idno[idno]['relatives']['parents']:
            for idno_parent, _ in d_idno[idno]['relatives'].items():
                if not _ == 4:
                    continue
##                d_idno[idno_parent]['relatives']['parents'].add(idno_num1)
##                d_idno[idno_num1]['relatives']['children'].add(idno_parent)
                if not idno_parent in d_idno[idno_num1]['relatives'].keys():
                    d_idno[idno_num1]['relatives'][idno_parent] = 4
                    d_idno[idno_parent]['relatives'][idno_num1] = 3
            ## the grandchild could be a sibling of other grandchildren
            for idno_grandchild, _ in d_idno[idno_num1]['relatives'].items():
                if not _ == 8:
                    continue
                if not idno_grandchild in d_idno[idno]['relatives'].keys():
                    d_idno[idno]['relatives'][idno_grandchild] = 5
                    d_idno[idno_grandchild]['relatives'][idno] = 5

##            if idno == '11093122':
##                print(idno, d_idno[idno]['relatives'])
####                print(idno_num1, d_idno[idno_num1]['relatives']['11093042'])
##                print(11093042, d_idno['11093042']['relatives'])
##                stop2a
##            if idno == '11093042':
##                print(d_idno[idno]['relatives'])
##                print(d_idno[idno_num1]['relatives'])
##                print(d_idno['11093122']['relatives'])
##                stop2b

            continue

        ## Head of household error.
        ## ERROR3C
        if d['rel1'] == '1' and d['relh'] != '1':
            with open('rel1_eq_1_relh_ne_1.txt', 'a') as f:
                print(d['idno'], d['agec'], d['sex'], d['num'], d['relh'], d['num1'], d['rel1'], d['vhno'], d['sangerID'], file=f, sep='\t')
            continue

##        ## Head of household error.
##        if d['relh'] == '1' and d['rel1'] != '1':
##            with open('relh_eq_1_rel1_ne_1.txt', 'a') as f:
##                print(d['idno'], d['agec'], d['sex'], d['num'], d['relh'], d['num1'], d['rel1'], d['vhno'], d['sangerID'], file=f, sep='\t')
##            continue

        print(idno)
        print('num1', d['num1'])
        print('***rel1', d['rel1'])
        print('num', d['num'])
        print('relh', d['relh'])
        print('sex', d['sex'])
        print('vhno', d['vhno'])
        print('ridno', d['ridno'])
        print(d.keys())
        print('idno', d_vhno[d['vhno']][d['num']])
        print('xxxx')
        stop

    ## Assertions
    for idno, d in d_idno.items():
        ## No more than 4 grand parents.
        assert len(d['relatives']['grandparents']) <= 4
        ## No more than 2 parents.
        assert len(d['relatives']['parents']) <= 2
        ## No more than 1 husband for females.
        if d['sex'] == 2:
            assert len(d['relatives']['spouses']) <= 1
        ## Siblings must share at least one parent,
        ## if one of them have both parents specified.
        if d['relatives']['siblings'] and len(d['relatives']['parents']) == 2:
            for idno_sibling in d['relatives']['siblings']:
                if len(d_idno[idno_sibling]['relatives']['parents']) == 2:
                    assert len(
                        d_idno[idno_sibling]['relatives']['parents'] & d_idno[idno]['relatives']['parents']) > 0

    d_seq2array = {}
    with open('../../../../metadata/samples.tsv') as f:
        for line in f:
            l = line.split('\t')
            seqID = l[0]
            arrayID = l[1]
            d_seq2array[seqID] = arrayID

    ## Suddenly Javier wants "non-genotyped males with WGS for both parents".
    ## WTF!!! Wish he could make up his mind!!! FFS!!!
    d_array2seq = dict((v,k) for k,v in d_seq2array.items())
    for idno, d in d_idno.items():
        ## Males only.
        if d['sex'] != '1':
            continue
        set_parents = d['relatives']['parents'] | set((k for k in d['relatives'].keys() if d['relatives'][k] == 4))
        if len(set_parents) != 2:
            continue
        try:
            ID_child = d_idno2sangerID[idno]
        except KeyError:
            continue
        if ID_child == '':
            continue
        try:
            ID_child_EGAN = d_array2seq[ID_child]
            bool_child_sequenced = 1
        except KeyError:
            ID_child_EGAN = ''
            bool_child_sequenced = 0

        ID_father = None
        ID_mother = None
        for idno_parent in set_parents:
            try:
                ID_parent = d_array2seq[d_idno2sangerID[idno_parent]]
            except KeyError:
                break
            if d_idno[idno_parent]['sex'] == '1':
                ID_father = ID_parent
            else:
                ID_mother = ID_parent

        if ID_father and ID_mother:
            print(
                ID_child,
                bool_child_sequenced,
                '{}/{}/{}'.format(idno, ID_child, ID_child_EGAN),
                '{}/{}'.format(ID_father, d_seq2array[ID_father]),
                '{}/{}'.format(ID_mother, d_seq2array[ID_mother]),
                sep='\t')
    stop_2016may03

    d_ped = {}
    with gzip.open(args.genome, 'rt') as f, open('unexplained_rel.txt', 'w') as f2:
        for line in f:
            break
        for line in f:
            l = line.split()

            ## Parse the IBD from the split line.
            PI_HAT = float(l[9])
##            if PI_HAT < 0.125+0.025:

            ## Make sure PI_HAT is close to 0.5 (parent/child).
            if PI_HAT < 0.4:
                continue

            ## Parse IDs from the split line.
            IID1 = l[1]
            IID2 = l[3]

            try:
                ID1 = d_seq2array[IID1]
            except:
                ID1 = IID1[IID1.rindex('_')+1:]

            try:
                ID2 = d_seq2array[IID2]
            except:
                ID2 = IID2[IID2.rindex('_')+1:]

            print(IID1, IID2, ID1, ID2, PI_HAT)

            ## Get idnos.
            try:
                idno1 = d_sangerID2idno[ID1]
            except KeyError:
                continue
            try:
                idno2 = d_sangerID2idno[ID2]
            except KeyError:
                continue

            ## Skip if different households.
            if d_idno[idno1]['vhno'] != d_idno[idno2]['vhno']:
                continue
            ## Skip if from a household to be skipped.
            if d_idno[idno1]['vhno'] in skipped_households:
                continue

            print(idno1, idno2)
            print('vhno', d_idno[idno1]['vhno'])
##            print(idno1, d_idno[idno1]['relatives'])
##            print(idno2, d_idno[idno2]['relatives'])
            try:
                print(idno1, d_idno[idno2]['relatives'][idno1])
                print(idno2, d_idno[idno1]['relatives'][idno2])
            except KeyError:
                print(
                    d_idno[idno1]['vhno'],
                    idno1, idno2,
                    d_idno[idno2]['rel1'], d_idno[idno1]['rel1'],
                    ID1, ID2, IID1, IID2, PI_HAT,
                    file=f2, sep='\t')
                continue
            rel1 = d_idno[idno2]['relatives'][idno1]
            rel2 = d_idno[idno1]['relatives'][idno2]
            if rel1 == 3:
                assert rel2 == 4
                try:
                    d_ped[idno1].add(idno2)
                except KeyError:
                    d_ped[idno1] = set([idno2])
            if rel2 == 3:
                assert rel1 == 4
                try:
                    d_ped[idno2].add(idno1)
                except KeyError:
                    d_ped[idno2] = set([idno1])

    ## without grandparents
    with open('pedigree.txt', 'w') as f:
        for i, (idno_child, idno_parents) in enumerate(d_ped.items()):
            try:
                assert len(idno_parents) <= 2
            except AssertionError:
                print(d_idno[idno_child]['vhno'], idno_child, idno_parents)
##            if len(idno_parents) < 2:
##                continue
            sangerID_mother = 'NA'+str(i)
            sangerID_father = 'NA'+str(i)
            for idno_parent in idno_parents:
                if d_idno[idno_parent]['sex'] == '1':
                    sangerID_father = d_idno2sangerID[idno_parent]
                else:
                    assert d_idno[idno_parent]['sex'] == '2'
                    sangerID_mother = d_idno2sangerID[idno_parent]
            print(
                d_idno2sangerID[idno_child],
                sangerID_mother, sangerID_father,
                file=f, sep='\t')

    ## with grandparents
    with open('pedigree.txt', 'w') as f:
        for i, (idno_child, idno_parents) in enumerate(d_ped.items()):
##        for i, idno_child in enumerate(d_ped.keys()):
##            ## Make a new mutable object instead of mutating the reference.
##            idno_parents = set(d_ped[idno_child])
            try:
                assert len(idno_parents) <= 2
            except AssertionError:
                print(d_idno[idno_child]['vhno'], idno_child, idno_parents)
##            if len(idno_parents) < 2:
##                idno_parents.add('NA')
####                continue
            sangerID_mother = 'NA.F'+str(i)
            sangerID_father = 'NA.M'+str(i)
            sangerID_grandparents = []
            for idno_parent in idno_parents:
                if d_idno[idno_parent]['sex'] == '1':
                    sangerID_father = d_idno2sangerID[idno_parent]
                    idno_father = idno_parent
                else:
                    assert d_idno[idno_parent]['sex'] == '2'
                    sangerID_mother = d_idno2sangerID[idno_parent]
                    idno_mother = idno_parent
                try:
                    ## Make a new mutable object.
                    idno_grandparents = set(d_ped[idno_parent])
                except KeyError:
                    idno_grandparents = set(['NA'+str(i), 'NA'+str(i)])
                if len(idno_grandparents) == 1:
                    idno_grandparents.add('NA')
                assert len(idno_grandparents) == 2
                for idno_grandparent in idno_grandparents:
                    try:
                        sangerID_grandparent = d_idno2sangerID[idno_grandparent]
                    except KeyError:
                        sangerID_grandparent = 'NA'
                    sangerID_grandparents.append(sangerID_grandparent)
            if len(idno_parents) == 1:
                sangerID_grandparents.append('NA')
                sangerID_grandparents.append('NA')
            print(
                d_idno2sangerID[idno_child],
                sangerID_mother, sangerID_father,
                '\t'.join(sangerID_grandparents),
                file=f, sep='\t')

    return


def add_relation(d_idno, idno1, idno2, rel1a, rel1b, rel1c):

##    d_relno2reltxt = {
##        4: 'parents',
##        3: 'children',
##        2: 'spouses',
##        5: 'siblings',
##        8: 'grandchildren',
##        9: 'grandparents',
##        }

    if idno1 == '14033067':
        stop1
    for idno_rel, _ in d_idno[idno1]['relatives'].items():
        ## Skip if self relation; e.g. sibling of self.
        if idno_rel == idno2:
            continue
        ## Skip if not correct type of relation.
        if not _ == rel1a:
            continue
        ## Append relation.
        if not idno_rel in d_idno[idno2]['relatives'].keys():
            d_idno[idno2]['relatives'][idno_rel] = rel1b
            d_idno[idno_rel]['relatives'][idno2] = rel1c
##            d_idno[idno2]['relatives'][d_relno2reltxt[rel1b]].add(idno_rel)
##            d_idno[idno_rel]['relatives'][d_relno2reltxt[rel1c]].add(idno2)

    return


def manual_correction(d_idno, idno):

    if idno == '12085135':
        d_idno[idno]['num1'] = '1'
    if idno == '21719075':
        d_idno[idno]['relh'] = '1'
    if idno in (
        '10044016', '10044118', '10044129', '10044141', '10044152'):
        d_idno[idno]['num1'] = '4'
    if idno == '1091010':
        d_idno[idno]['relh'] = '8'
    if idno == '17137025':
        d_idno[idno]['relh'] = '2'
        d_idno[idno]['rel1'] = '2'
    if idno in (
        '17046053', '17806015', '17806059'):
        d_idno[idno]['num1'] = '1'
    if idno == '25015012':
        d_idno[idno]['num1'] = '3'
    if idno == '19025074':
        d_idno[idno]['num'] = '1'
    if idno == '16047011':
        d_idno[idno]['num1'] = '1'
    if idno == '14120151':
        d_idno[idno]['num'] = '2'
    if idno == '1131073':
        d_idno[idno]['relh'] = '1'
    if idno == '1131095':
        d_idno[idno]['relh'] = '2'
    if idno == '19003158':
        d_idno[idno]['num1'] = '3'
    if idno == '6297062':
        d_idno[idno]['num1'] = '3'
    if idno == '6297073':
        d_idno[idno]['num1'] = '3'
    if idno == '17868011':
        d_idno[idno]['relh'] = '1'
    if idno == '3131019':
        d_idno[idno]['relh'] = '1'
        d_idno[idno]['num1'] = '1'
    if idno == '6107026':
        d_idno[idno]['rel1'] = '2'
    if idno == '20660019':
        d_idno[idno]['relh'] = '1'
    if idno == '12050033':
        d_idno[idno]['relh'] = '1'
    if idno == '9096115':
        d_idno[idno]['relh'] = '1'
    if idno == '10133106':
        d_idno[idno]['relh'] = '1'
    if idno == '9148088':
        d_idno[idno]['relh'] = '1'
    if idno == '17018054':
        d_idno[idno]['rel1'] = '8'
    if idno == '17021025':
        d_idno[idno]['rel1'] = '2'
    if idno == '8138117':
        d_idno[idno]['rel1'] = '8'
    if idno == '8138026':
        d_idno[idno]['rel1'] = '2'
    if idno == '11138011':
        d_idno[idno]['num'] = '1'
        d_idno[idno]['num1'] = '1'
        d_idno[idno]['rel1'] = '1'
    if idno == '11138022':
        d_idno[idno]['num1'] = '1'
    if idno == '24602031':
        d_idno[idno]['rel1'] = '8'
    if idno == '23411016':
        d_idno[idno]['rel1'] = '2'
    if idno == '10487013':
        d_idno[idno]['num1'] = '2'
    if idno == '1127054':
        d_idno[idno]['relh'] = '8'
    if idno == '2022086':
        d_idno[idno]['num1'] = '1'
    if idno == '11093020':
        d_idno[idno]['num1'] = '9'
    if idno == '1021026':
        d_idno[idno]['num1'] = '1'
    if idno == '11068279':
        d_idno[idno]['num1'] = '1'
    if idno == '1127054':
        d_idno[idno]['relh'] = '8'
        d_idno[idno]['num1'] = '1'
    if idno == '12048021':
        d_idno[idno]['num1'] = '1'
    if idno == '12054024':
        d_idno[idno]['rel1'] = '1'
        d_idno[idno]['relh'] = '1'
    if idno == '12097125':
        d_idno[idno]['num1'] = '1'
    if idno == '12551019':
        d_idno[idno]['num1'] = '1'
    if idno == '14802032':
        d_idno[idno]['num1'] = '1'
    if idno == '16022016':
        d_idno[idno]['rel1'] = '1'
        d_idno[idno]['relh'] = '1'
    if idno == '17510011':
        d_idno[idno]['num1'] = '2'
    if idno == '18078179':
        d_idno[idno]['num1'] = '1'
    if idno == '18608041':
        d_idno[idno]['num1'] = '1'
    if idno == '20030056':
        d_idno[idno]['num1'] = '1'
    if idno == '21038095':
        d_idno[idno]['num1'] = '3'
    if idno == '23132018':
        d_idno[idno]['num1'] = '1'
    if idno == '25090108':
        d_idno[idno]['rel1'] = '1'
        d_idno[idno]['relh'] = '1'
    if idno == '3007026':
        d_idno[idno]['num'] = '2'
    if idno == '3027048':
        d_idno[idno]['num1'] = '2'
    if idno == '9086041':
        d_idno[idno]['num1'] = '2'
    if idno == '9196010':
        d_idno[idno]['relh'] = '1'
        d_idno[idno]['rel1'] = '1'
    if idno == '14111126':
        d_idno[idno]['relh'] = '3'
        d_idno[idno]['rel1'] = '3'
    if idno == '14704111':
        d_idno[idno]['num1'] = '3'
        d_idno[idno]['rel1'] = '5'
    if idno == '4112024':
        d_idno[idno]['relh'] = '3'
        d_idno[idno]['rel1'] = '3'
    if idno == '6138057':
        d_idno[idno]['relh'] = '3'
        d_idno[idno]['rel1'] = '3'
        d_idno[idno]['num1'] = '8'
    if idno == '2057020':
        d_idno[idno]['num1'] = '1'
    if idno == '14952012':
        d_idno[idno]['num1'] = '4'
        d_idno[idno]['rel1'] = '5'
    if idno == '11049078':
        d_idno[idno]['rel1'] = '5'
        d_idno[idno]['num1'] = '1'
    if idno == '6026060':
        d_idno[idno]['rel1'] = '3'
    if idno == '5036100':
        d_idno[idno]['num1'] = '3'
        d_idno[idno]['rel1'] = '3'
    if idno == '21162055':
        d_idno[idno]['relh'] = '3'
        d_idno[idno]['rel1'] = '3'
    if idno == '21162077':
        d_idno[idno]['relh'] = '3'
        d_idno[idno]['rel1'] = '3'
    if idno == '8015053':
        d_idno[idno]['relh'] = '3'
        d_idno[idno]['rel1'] = '3'
    if idno == '6150031':
        d_idno[idno]['num1'] = '1'
        d_idno[idno]['rel1'] = '5'
    if idno == '15211044':
        d_idno[idno]['relh'] = '3'
        d_idno[idno]['rel1'] = '3'
    if idno == '20512024':
        d_idno[idno]['relh'] = '2'
        d_idno[idno]['rel1'] = '2'
    if idno == '22031255':
        d_idno[idno]['relh'] = '3'
        d_idno[idno]['rel1'] = '3'

    return d_idno


def parse_args():

    parser = argparse.ArgumentParser()

    parser.add_argument('--tsv', help='phenotypic relationships')

    parser.add_argument('--genome', help='genotypic relationships')

    args = parser.parse_args()

    return args


def calculate_stats_for_deepti(d_idno,d_vhno,):

    count_both = 0
    count_father = 0
    count_mother = 0
    for idno in d_idno.keys():

        ## not a child of anyone
        if not 3 in d_idno[idno]['rel1'].keys():
            continue

        vhno = d_idno[idno]['vhno']

        father = False
        mother = False

        for parent in d_idno[idno]['rel1'][3]:
            ## household deleted during curation? due to self relation...
            if not vhno in d_vhno.keys():
                continue
            if not parent[0] in d_vhno[vhno]['num'].keys():
                continue
            idno_parent = d_vhno[vhno]['num'][parent[0]]['idno']

##            for parent in d_idno[idno]['relatives']['parents']:
##                idno_parent = parent[0]

            sex_parent = d_idno[idno_parent]['sex']
            if sex_parent == 1:
                father = True
            else:
                mother = True

        if father == True and mother == True:
            count_both += 1
        if father == True:
            count_father += 1
        if mother == True:
            count_mother += 1
            
##            num = d_idno[idno]['num']
##            num1 = d_idno[idno]['num1']
##            vhno = d_idno[idno]['vhno']
##            ## household deleted during curation? due to self relation...
##            if not vhno in d_vhno.keys():
##                continue
##
##            ## relative left household
##            if not num1 in d_vhno[vhno]['num'].keys():
##                continue
##            else:
##                idno_num1 = d_vhno[vhno]['num'][num1]['idno']
##
##            ## add to count
##            if d_idno[idno_num1]['sex'] == 1:
##                count_fathers += 1
##            else:
##                count_mothers += 1

##            ## no children
##            if not 4 in d_idno[idno]['rel1'].keys():
##                continue
##            for num1,amb in d_idno[idno]['rel1'][4]:
##                if not num1 in d_vhno[vhno]['num'].keys():
##                    continue
##                idno_num1 = d_vhno[vhno]['num'][num1]['idno']
##                ## add to count
##                if d_idno[idno_num1]['sex'] == 1:
##                    count_fathers += 1
##                else:
##                    count_mothers += 1

    print('father', count_father, 'mother', count_mother, 'both', count_both)

    return


def check_age2(d_idno,):

    '''check age differences after building ambiguous relationships'''

    for idno in d_idno.keys():
        agec = d_idno[idno]['agec']
        for relation in d_idno[idno]['relatives'].keys():
            for relative in d_idno[idno]['relatives'][relation]:
                idno_relative = relative[0]
                agec_relative = d_idno[idno_relative]['agec']
                ambiguity = relative[1]
                if relation == 'spouses':
                    if abs(agec-agec_relative) > (87-22):
                        print()
                        print(idno, idno_relative)
                        print(agec, agec_relative)
                        print('spouse')
                        stop_spouse
                elif relation == 'siblings':
                    if abs(agec_relative-agec) > (68-16):
                        print()
                        print(agec, agec_relative)
                        print(idno, idno_relative)
                        print('sibling')
                        stop_sibling
                elif relation == 'grandparents':
                    if agec > agec_relative or agec_relative-agec > (106-15):
                        print()
                        print(agec_relative, agec)
                        print(idno_relative, idno)
                        stop_grandparent
                elif relation == 'grandchildren':
                    if agec < agec_relative or agec-agec_relative > (106-15):
                        print()
                        print(agec, agec_relative)
                        print(idno_relative, idno)
                        stop_grandchild
                elif relation == 'parents':
                    if (
                        agec > agec_relative
                        or
                        (ambiguity == 0 and agec_relative-agec > (81-19))
                        or
                        (ambiguity == 1 and agec_relative-agec > (88-15))
                        ):
                        print()
                        print(ambiguity)
                        print(agec_relative, agec)
                        print(idno_relative, idno)
                        print('parent')
                        stop_parent
                        if agec > agec_relative:
                            stop_parent
                elif relation == 'children':
                    if (
                        agec < agec_relative
                        or
                        (ambiguity == 0 and agec-agec_relative > (81-19))
                        or
                        (ambiguity == 1 and agec-agec_relative > (88-15))
                        ):
                        print()
                        print(agec_relative, agec)
                        print(idno_relative, idno)
                        print('child')
                        stop_children
                        if agec < agec_relative:
                            stop_children
                else:
                    print(relation, agec, agec_relative)
                    stop_new_relation

    return


def read_gene_comparison_fast(treshold,):

    '''this method is fast (26 seconds / 2 seconds), but it's dependent on the columns of compare_ibs_pairs_LDpruned.genome always having the same width'''

    d_PI_HAT = {}

    fn = 'in/compare_ibs_pairs_LDpruned.genome'
    if not os.path.isfile(fn):
        return d_PI_HAT
    
    fd = open(fn)
    for line in fd:
        break
    for line in fd:

        ## parse data
        IID1 = line[2:23]
        IID2 = line[48:69]
        PI_HAT_calculated = PI_HAT = float(line[127:133])
##            if IID1 == '197828_B05_APP5119199' and IID2 == '197828_F09_APP5119245':
##                print PI_HAT_calculated
##            if IID2 == '197828_B05_APP5119199' and IID1 == '197828_F09_APP5119245':
##                print PI_HAT_calculated
        if PI_HAT_calculated < treshold:
            continue
        d_PI_HAT = self.append_to_d_PI_HAT(d_PI_HAT,IID1,IID2,PI_HAT,)

    fd.close()

    return d_PI_HAT


def append_to_d_PI_HAT(d_PI_HAT,IID1,IID2,PI_HAT,):

    if not IID1 in d_PI_HAT.keys():
        d_PI_HAT[IID1] = {}
    d_PI_HAT[IID1][IID2] = PI_HAT
    ## do reverse (which is slow...)
    if not IID2 in d_PI_HAT.keys():
        d_PI_HAT[IID2] = {}
    d_PI_HAT[IID2][IID1] = PI_HAT

    return d_PI_HAT


def read_gene_comparison_slow(treshold,):

    '''this method is slow (32 seconds / 7 seconds), but it's not dependent on the columns of compare_ibs_pairs_LDpruned.genome always having the same width'''

    fd = open('in/compare_ibs_pairs_LDpruned.genome')
    l_header = None
    d_PI_HAT = {}
    for line in fd:

        ## split line
        l = line.strip().split()

        ## parse header
        if not l_header:
            l_header = l
            continue

        ## parse data
        IID1 = l[l_header.index('IID1')]
        IID2 = l[l_header.index('IID2')]
        PI_HAT_calculated = PI_HAT = float(l[l_header.index('PI_HAT')])
        if PI_HAT_calculated < treshold:
            continue
        d_PI_HAT = self.append_to_d_PI_HAT(d_PI_HAT,IID1,IID2,PI_HAT,)

    fd.close()

    return d_PI_HAT


def find_errors_and_exclude_from_dataset2(d_idno,):

    return d_idno



def line2parameters(line,l_headers,):
    
    l = line.strip().split('\t')
    gwas_id = l[l_headers.index('gwas_id')][1:-1]
    idno = int(l[l_headers.index('idno')])
    d = {
        'idno':idno,
        'gwas_id':gwas_id,
        }
    ## no known relationship and not a member of a household
    if len(l) == 2:
        vhno = num = num1 = rel1 = sex = agec = ''
    else:
        vhno = int(l[l_headers.index('vhno')])
        num = int(l[l_headers.index('num')])
        num1 = l[l_headers.index('num1')]
        if num1 != '':
            num1 = int(num1)
        rel1 = l[l_headers.index('rel1')]
        if rel1 != '':
            rel1 = int(rel1)
        sex = int(l[l_headers.index('sex')])
        agec = int(l[l_headers.index('agec')])
    d['idno'] = idno
    d['vhno'] = vhno
    d['num'] = num
    d['num1'] = num1
    d['rel1'] = rel1
    d['sex'] = sex
    d['agec'] = agec
    d['gwas_id'] = gwas_id

    return d


def append_line_to_dictionary_of_individuals(
    self,
    d_idno, idno,
    num, rel1, num1,
    sex, agec,
    gwas_id,
    vhno,
    ):

    d_idno[idno] = {
        'relatives':{
            ## own generation
            'spouses':[],
            'siblings':[],
            ## next and future generations
            'children':[],
            'grandchildren':[],
            ## previous and past generations
            'parents':[],
            'grandparents':[],
            },
        'sex':sex,
        'agec':agec,
        'num':num,
        'gwas_id':gwas_id,
        'rel1':{rel1:[[num1,0,]],},
        'num1':num1,
        'vhno':vhno,
        }

    return d_idno


def append_line_to_dictionary_of_households(
    self,
    d_vhno,vhno,
    idno,
    num, rel1, num1,
    sex, agec,
    ):

    ## append vhno to dictionary
    if not vhno in d_vhno.keys():
        d_vhno[vhno] = {'num':{},'rel1':{},}
    ## append num to dictionary
    d_vhno[vhno]['num'][num] = {
        'idno':idno,'num1':num1,'rel1':rel1,
        ## add sex to know if parent is female/mother or male/father
        'sex':sex,
        ## add age to check for age discrepancies
        'agec':agec,
        }
    ## append rel1 to dictionary
    if not rel1 in d_vhno[vhno]['rel1'].keys():
        d_vhno[vhno]['rel1'][rel1] = []
    d_vhno[vhno]['rel1'][rel1] += [num]

    return d_vhno


def loop_lines(lines, l_headers,):

    d_idno = {}
    d_vhno = {}
    for i_line in range(1, len(lines)):

        line = lines[i_line]

        d = self.line2parameters(line,l_headers,)

        idno = d['idno']
        vhno = d['vhno']
        num = d['num']
        num1 = d['num1']
        rel1 = d['rel1']
        sex = d['sex']
        agec = d['agec']
        gwas_id = d['gwas_id']

        ##
        ## manual deletion
        ##
        fd = open('in/exclusion.txt','r')
        s = fd.read()
        fd.close()
        l_exclusion = eval(s)
        if idno in l_exclusion:
            continue

        ##
        ## manual intervention
        ##
        ## head should always be related to self
        if rel1 == 1 and num != num1:
            num1 = num

        d_idno = self.append_line_to_dictionary_of_individuals(
            d_idno, idno,
            num, rel1, num1,
            sex, agec,
            gwas_id,
            vhno,
            )

        d_vhno = self.append_line_to_dictionary_of_households(
            d_vhno,vhno,
            idno,
            num, rel1, num1,
            sex, agec,
            )

    return d_vhno, d_idno


def redflag_household_multiple_heads(d_vhno,vhno,):

    ## remove households with multiple heads of household
    if len(d_vhno[vhno]['rel1'][1]) > 1:
        if self.verbose == True:
##                print 'temporarily removing vhno %i from dataset \
##because of multiple heads of household' %(vhno)
            fd = open('out/redflag_multiple_heads.txt','a')
            fd.write('%s\n' %(vhno))
            fd.close()

##            del d_vhno[vhno]
##            bool_continue = True
##            num = None

        bool_continue = False
        num = d_vhno[vhno]['rel1'][1][0]

    else:
        num = d_vhno[vhno]['rel1'][1][0]
        bool_continue = False

    return d_vhno, bool_continue, num


def check_age(d_vhno, vhno,):

    for num in d_vhno[vhno]['num'].keys():
        ## age of individual
        agec = d_vhno[vhno]['num'][num]['agec']
        ## number of related individual
        num1 = d_vhno[vhno]['num'][num]['num1']
        ## relationship type
        rel1 = d_vhno[vhno]['num'][num]['rel1']
        if not num1 in d_vhno[vhno]['num'].keys():
            continue
        ## age of related individual
        agec_num1 = d_vhno[vhno]['num'][num1]['agec']
        if rel1 == 1:
            if abs(agec-agec_num1) > 20: ## should be zero, but sometimes multiple heads...
                print(vhno, num, num1)
                print(agec, agec_num1)
                stopage1
        elif rel1 == 2:
            ## no more than xxx years between spouses
            if abs(agec-agec_num1) > abs(22-87):
                print(vhno, num, num1)
                print(agec, agec_num1)
                stopage2
        elif rel1 == 3:
            ## no more than xxx years between child and parent
            if agec-agec_num1 > 73-15:
                print(vhno, num, num1)
                print(agec, agec_num1)
                stopage3
        elif rel1 == 4:
            ## no more than xxx years between parent and child
            if agec-agec_num1 > 85-30:
                print(vhno, num, num1)
                print(agec, agec_num1)
                stopage4
        elif rel1 == 5:
            ## no more than xxx years between siblings
            if abs(agec-agec_num1) > abs(68-16):
                print(vhno, num, num1)
                print(agec, agec_num1)
                stopage5
        elif rel1 in [6,7,]:
            pass
        elif rel1 == 8:
            ## no more than xxx years between grandchild and parent
            if agec-agec_num1 > 93-13:
                print(vhno, num, num1)
                print(agec, agec_num1)
                stopage8
        elif rel1 in [21,22,]:
            pass
        elif rel1 == '':
            pass
        else:
            print(rel1, agec, agec_num1)
            stop

    return


def print_removal(vhno,num,rel1,num1,):

    if rel1 == 3:
        print('temporarily removing vhno %i from dataset \
because individual %i is child of self' %(vhno,num,))
    elif rel1 == 2:
        print('temporarily removing vhno %i from dataset \
because individual %i is spouse of self' %(vhno,num,))
    elif rel1 == 5:
        print('temporarily removing vhno %i from dataset \
because individual %i is sibling of self' %(vhno,num,))
    elif rel1 == 8:
        print('temporarily removing vhno %i from dataset \
because individual %i is grandchild of self' %(vhno,num,))
    else:
        print('relationship to self', vhno, num, rel1, num1)
        stop

    return


def error_self_relation(d_vhno,vhno,):

    bool_continue = False
    ## loop over individuals of household
    for num in d_vhno[vhno]['num'].keys():
        rel1 = d_vhno[vhno]['num'][num]['rel1']
        num1 = d_vhno[vhno]['num'][num]['num1']
        ## individual related to self and not head of household?
        if rel1 != 1 and num == num1:
            ## spouse,child,sibling or grandchild of self?
            if not rel1 in [2,3,5,8,]:
                continue
            bool_continue = True
            del d_vhno[vhno]
            if self.verbose == True:
##                    self.print_removal(vhno,num,rel1,num1,)
                fd = open('out/error_self_relation.txt','a')
                fd.write('%s %s %s %s\n' %(vhno,num,rel1,num1,))
                fd.close()
            break

    return d_vhno, bool_continue


def find_errors_and_exclude_from_dataset1(d_vhno,d_idno,):

    '''this data curation step is extremly important,
as the pedigree table contains errors, and although some of them are obvious,
it would be guessing and assuming at the end of the day to fix them.'''

    if self.verbose == True:
        d_columns = {
            'out/error_self_relation.txt':['vhno','num','rel1','num1',],
            'out/redflag_multiple_heads.txt':['vhno',],
            'out/redflag_head_not_related_to_self.txt':['vhno',],
            }
        for fn in d_columns.keys():
            l_columns = d_columns[fn]
            fd = open(fn,'w')
            fd.write('##%s\n' %('\t'.join(l_columns)))
            fd.close()

    ## how many households prior to cleanup?
    count_households1 = len(d_vhno.keys())

    ## loop over households
    for vhno in d_vhno.keys():

        ## check households with head of household
        if 1 in d_vhno[vhno]['rel1'].keys():

            ## remove household if multiple heads of household
            (
                d_vhno, bool_continue, num,
                ) = self.redflag_household_multiple_heads(
                    d_vhno,vhno,
                    )
            if bool_continue == True:
                continue
            else:
                ## remove household if head is not related to self
                if d_vhno[vhno]['num'][num]['num1'] != num:
                    if self.verbose == True:
##                        print 'temporarily removing vhno %i from dataset \
##because head of household is not related to self' %(vhno)
                        fd = open('out/redflag_head_not_related_to_self.txt','a')
                        fd.write('%s\n' %(vhno,))
                        fd.close()
##                    del d_vhno[vhno]
##                    continue

        self.check_age(d_vhno, vhno,)

        ## remove households in which a person is related to self
        d_vhno, bool_continue = self.error_self_relation(d_vhno,vhno,)
        if bool_continue == True:
            continue

    ## how many households after cleanup?
    count_households2 = len(d_vhno.keys())

    if self.verbose == True:
        print('removed %i of %i households (%.1f%%) because of data errors' %(
            count_households1-count_households2,
            count_households1,
            100*(count_households1-count_households2)/float(count_households1)
            ))

    return d_vhno, d_idno


def sub_head(d_vhno,vhno,num,d_idno,idno,):

    num1 = d_vhno[vhno]['num'][num]['num1']

    ## head related to self
    if num1 == num:
        bool_append = False
        pass
    ## head not related to self (error!)
    else:
        bool_append = False
        pass

    return d_idno, bool_append


def sub_spouse(d_vhno,vhno,num,d_idno,idno,):

    '''spouse can be male/female
and spouse can be the parent of child of spouse'''

    bool_return, bool_append, num1, idno_num1 = self.initiate_sub_relative(
        d_vhno,vhno,num,
        )
    if bool_return == True: return d_idno, bool_append

    ##
    ## 1) append non-ambiguous spouse to *spouse*
    ##
    d_idno, bool_append = self.append_non_ambiguous_relative(
        d_idno,
        idno, idno_num1,
        num, num1,
        'spouses',2,
        'spouses',2,
        )

    ##
    ## 2) append ambiguous child (child of non-ambiguous spouse)
    ## to current *spouse* (if any) and vice versa
    ##
    ## does spouse have a child?
    if not 3 in d_idno[idno_num1]['rel1'].keys():
        bool_append_sub = False
    else:
        d_idno, bool_append_sub = self.append_ambiguous_relative(
            d_vhno,vhno,d_idno,
            num,idno,
            num1,idno_num1,
            ## what is relationship of ambiguous relative to current relative
            'children', 3,
            ## what is relationship of current relative to ambiguous relative
            'parents', 4,
            ## what is relationship of non-ambiguous relative (spouse)
            ## to ambiguous relative (child)
            ## non-ambiguous relative (spouse) is *parent* of ambiguous relative
            ## i.e. which relatives of the non-ambiguous relative are we looping over?
            4,
            )
        if bool_append_sub == True: bool_append = True

    return d_idno, bool_append


def initiate_sub_relative(d_vhno,vhno,num,):

    bool_append = False
    bool_return = False

    ## num of relative
    num1 = d_vhno[vhno]['num'][num]['num1']

    ## relative no longer in household?
    if not num1 in d_vhno[vhno]['num'].keys():
        bool_return = True
        return bool_return, bool_append, num1, None

    ## idno of relative
    idno_num1 = d_vhno[vhno]['num'][num1]['idno']

    return bool_return, bool_append, num1, idno_num1


def sub_grandchild(d_vhno,vhno,num,d_idno,idno,):

    '''PI_HAT for a grandchild is .25'''

    bool_return, bool_append, num1, idno_num1 = self.initiate_sub_relative(
        d_vhno,vhno,num,
        )
    if bool_return == True: return d_idno, bool_append

    ##
    ## 1) append non-ambiguous grandparent to *grandchild*
    ##
    d_idno, bool_append_sub = self.append_non_ambiguous_relative(
        d_idno,
        idno, idno_num1, num, num1,
        'grandchildren',8,
        'grandparents',9, ## 9 is designated to grandparents by Tommy
        )
    if bool_append_sub == True: bool_append = True

    ##
    ## 2) append ambiguous parent (child of grandparent)
    ## to current *grandchild* (if any) and vice versa
    ##
    ## does grandparent have a child? i.e. is the grandparent a parent of someone?
    if 4 in d_idno[idno_num1]['rel1'].keys():
        d_idno, bool_append_sub = self.append_ambiguous_relative(
            d_vhno,vhno,d_idno,
            num,idno,
            num1,idno_num1,
            ## what is relationship of ambiguous relative to current relative
            'parents', 4,
            ## what is relationship of current relative to ambiguous relative
            'children', 3,
            ## what is relationship of non-ambiguous relative (grandparent)
            ## to ambiguous relative (child)
            ## i.e. which relatives of the non-ambiguous relative are we looping over?
            4,
            )
        if bool_append_sub == True: bool_append = True

    ##
    ## 3) append ambiguous sibling (grandchild of grandparent)
    ## to current *grandchild* (if any) and vice versa
    ##
    d_idno, bool_append_sub = self.append_ambiguous_relative(
        d_vhno,vhno,d_idno,
        num,idno,
        num1,idno_num1,
        ## what is relationship of ambiguous relative to current relative
        'siblings', 5,
        ## what is relationship of current relative to ambiguous relative
        'siblings', 5,
        ## what is relationship of non-ambiguous relative (grandparent)
        ## to ambiguous relative (sibling)
        ## i.e. which relatives of the non-ambiguous relative are we looping over?
        9,
        )
    if bool_append_sub == True: bool_append = True

    return d_idno, bool_append


def sub_child(d_vhno,vhno,num,d_idno,idno,):

    '''PI_HAT for a child is .5'''

    bool_return, bool_append, num1, idno_num1 = self.initiate_sub_relative(
        d_vhno,vhno,num,
        )
    if bool_return == True: return d_idno, bool_append

    ##
    ## 1) append non-ambiguous parent to *child*
    ##
    d_idno, bool_append_sub = self.append_non_ambiguous_relative(
        d_idno,
        idno, idno_num1, num, num1,
        'children',3,
        'parents',4,
        )
    if bool_append_sub == True: bool_append = True

    ##
    ## 2) append ambiguous parent (spouse of non-ambiguous parent)
    ## to current *child* (if any) and vice versa
    ##
    ## does parent have a spouse?
    if 2 in d_idno[idno_num1]['rel1'].keys():
        d_idno, bool_append_sub = self.append_ambiguous_relative(
            d_vhno,vhno,d_idno,
            num,idno,
            num1,idno_num1,
            ## what is relationship of ambiguous relative to current relative
            'parents', 4,
            ## what is relationship of current relative to ambiguous relative
            'children', 3,
            ## what is relationship of non-ambiguous relative (parent)
            ## to ambiguous relative (spouse)
            2,
            )
        if bool_append_sub == True: bool_append = True

    ##
    ## 3a) append non-ambiguous sibling (child of non-ambiguous parent)
    ## to current *child*
    ##
    for child in d_idno[idno_num1]['relatives']['children']:
        idno_sibling = child[0]
        ## do not append to self
        if idno == idno_sibling:
            continue
        ## do not append ambiguous child as non-ambiguous sibling
        if child[1] == 1:
            continue
        ## number within household
        num_sibling = d_idno[idno_sibling]['num']
##            ## shared parent? and thus non-ambiguous siblings
##            if idno_num1 != d_vhno[vhno]['num'][d_vhno[vhno]['num'][num_sibling]['num1']]['idno']:
##                continue
        d_idno, bool_append_sub = self.append_non_ambiguous_relative(
            d_idno,
            idno, idno_sibling, num, num_sibling,
            'siblings',5,
            'siblings',5,
            )
        if bool_append_sub == True: bool_append = True

    ##
    ## 3b) append ambiguous sibling (child of non-ambiguous parent)
    ## to current *child* (if any) and vice versa
    ##
    ## parent has a child (is a parent)!
    if not 4 in d_idno[idno_num1]['rel1'].keys():
        print(idno_num1, d_idno[idno_num1]['rel1'].keys())
        stop_not_expected
    d_idno, bool_append_sub = self.append_ambiguous_relative(
        d_vhno,vhno,d_idno,
        num,idno,
        num1,idno_num1,
        ## what is relationship of ambiguous relative to current relative
        'siblings', 5,
        ## what is relationship of current relative to ambiguous relative
        'siblings', 5,
        ## what is relationship of non-ambiguous relative (parent)
        ## to ambiguous relative (child)
        ## non-amibguous relative is a *parent* of the ambiguous relative
        ## current relative is a sibling of the ambiguous relative
        ## i.e. which relatives of the non-ambiguous relative are we looping over?
        4,
        )
    if bool_append_sub == True: bool_append = True

    return d_idno, bool_append


def sub_parent(d_vhno,vhno,num,d_idno,idno,):

    '''PI_HAT for a parent is .5'''

    bool_return, bool_append, num1, idno_num1 = self.initiate_sub_relative(
        d_vhno,vhno,num,
        )
    if bool_return == True: return d_idno, bool_append

    ##
    ## 1) append non-ambiguous child to *parent*
    ##
    d_idno, bool_append_sub = self.append_non_ambiguous_relative(
        d_idno,
        idno, idno_num1, num, num1,
        'parents',4,
        'children',3,
        )
    if bool_append_sub == True: bool_append = True

    return d_idno, bool_append


def append_non_ambiguous_relative(
    self,
    d_idno,
    idno, idno_num1, num, num1,
    current2nonambiguous, rel1,
    nonambiguous2current, rel1_reverse,
    ):

    bool_append = False

    ##        
    ## add non-ambiguous relative for current individual
    ##

    ## has this already been appended?
    if (
        [idno_num1,0,] in d_idno[idno]['relatives'][nonambiguous2current]
        ):
        pass
    else:
        d_idno[idno]['relatives'][nonambiguous2current] += [[idno_num1,0,]]
        bool_append = True

    if not rel1 in d_idno[idno]['rel1'].keys():
        d_idno[idno]['rel1'][rel1] = []
    if not [num1,0,] in d_idno[idno]['rel1'][rel1]:
        d_idno[idno]['rel1'][rel1] += [[num1,0,]]
        bool_append = True ## not really necessary...

    ##
    ## add non-ambiguous individual for relative of current individual
    ##

    ## has this already been appended?
    if (
        [idno,0,] in d_idno[idno_num1]['relatives'][current2nonambiguous]
        ):
        pass
    else:
        d_idno[idno_num1]['relatives'][current2nonambiguous] += [[idno,0,]]

        relative = [num,0,]
        if not rel1_reverse in d_idno[idno_num1]['rel1'].keys():
            d_idno[idno_num1]['rel1'][rel1_reverse] = []
        if not relative in d_idno[idno_num1]['rel1'][rel1_reverse]:
            d_idno[idno_num1]['rel1'][rel1_reverse] += [[num,0,]]
        else: ## for example spouse of spouse (idno 16020032)
            if bool_append == False:
                stop
            pass

        bool_append = True

    return d_idno, bool_append


def append_to_dict(
    self, d_idno,
    idno_a,idno_b,
    str_relation,
    rel1,
    num1,
    ):

    '''a sub function of append_ambiguous_relative'''

    bool_append = False

    ## append idno
    if [idno_a,0,] in d_idno[idno_b]['relatives'][str_relation]:
        bool_append = False
        pass
    elif [idno_a,1,] in d_idno[idno_b]['relatives'][str_relation]:
        bool_append = False
    else:
        d_idno[idno_b]['relatives'][str_relation] += [[idno_a,1,]]
        bool_append = True

    ## append num
    if not rel1 in d_idno[idno_b]['rel1'].keys():
        d_idno[idno_b]['rel1'][rel1] = []
    if [num1,0,] in d_idno[idno_b]['rel1'][rel1]:
        pass
    elif [num1,1,] in d_idno[idno_b]['rel1'][rel1]:
        bool_append = False
    else:
        d_idno[idno_b]['rel1'][rel1] += [[num1,1,]]
        bool_append = True

    return d_idno, bool_append


def append_ambiguous_relative(
    self,
    d_vhno,vhno,d_idno,
    num,idno,
    num1,idno_num1,
##        parent2,
##        str_current2nonambiguous, #int_current2,
##        str_nonambiguous2current, #int_nonambiguous2,
    str_ambiguous2current, int_ambiguous2current,
    str_current2ambiguous, int_current2ambiguous,
    int_nonambiguous2ambiguous,
    ):

    if [int_ambiguous2current,int_current2ambiguous,int_nonambiguous2ambiguous] not in [
        [5,5,9,], ## grandchild (loop over 9=grandchildren to find sibling/sibling)
        [4,3,4,], ## grandchild (loop over 4=children to find parent/child)
        [4,3,2,], ## child (loop over spouse to find parent/child)
        [5,5,4,], ## child (child1 of parent2, parent2 of child3)
        [4,3,3,], ## sibling (parent1, child2=sibling2, sibling3=child3)
        [5,5,5,], ## sibling (sibling1, sibling2=sibling2, sibling3=sibling3)
        [3,4,4,], ## spouse (child1 of parent2, parent2 of child1, parent2 of child3)
        ]:
        print([int_ambiguous2current,int_current2ambiguous,int_nonambiguous2ambiguous])
        stop

    bool_append = False

    idno_current = idno
    idno_nonambiguous = idno_num1
    num_current = num
    num_nonambiguous = num1

    ## nonambiguous relative of current relative
    ## does not have a relative
    ## that could be a nonambiguous relative
    ## of the current relative
    if not int_nonambiguous2ambiguous in d_idno[idno_nonambiguous]['rel1'].keys():
        return d_idno, bool_append
    ## loop over ambiguous relatives
    for relative_ambiguous in d_idno[idno_nonambiguous]['rel1'][int_nonambiguous2ambiguous]:
        num_ambiguous = relative_ambiguous[0]
        ## ambiguous relative of current relative still in household?
        if not num_ambiguous in d_vhno[vhno]['num'].keys():
            continue
        idno_ambiguous = d_vhno[vhno]['num'][num_ambiguous]['idno']

        ## do not append current relative to self
        if idno_current == idno_ambiguous:
            continue

##            ## do not append parent if a non-abiguous parent of that sex is already given
##            if (
##                int_ambiguous2current == 3
##                and
##                int_current2ambiguous == 4
##                and
##                int_nonambiguous2ambiguous == 4
##                ):
##                bool_continue = self.avoid_multiple_parents_and_false_predictions(
##                    d_idno,
##                    idno,idno_num1,idno_ambiguous,
##                    )
##                if bool_continue == True:
##                    continue

        ## do not append ambiguous child to parent and vice versa
        ## if child is older than parent
        if (
            (
                ## current relative is child of ambiguous relative
                int_current2ambiguous == 3
                and
                ## ambiguous relative is parent of current relative
                int_ambiguous2current == 4
                and
                ## child is older than parent?
                d_idno[idno]['agec'] > d_idno[idno_ambiguous]['agec']
                )
            or
            (
                ## current relative is parent of ambiguous relative
                int_current2ambiguous == 4
                and
                ## ambiguous relative is child of current relative
                int_ambiguous2current == 3
                and
                ## child is older than parent?
                d_idno[idno]['agec'] < d_idno[idno_ambiguous]['agec']
                )
            ):
                continue

        ##
        ## append ambiguous relative to current individual
        ##
        d_idno, bool_append_sub = self.append_to_dict(
            d_idno,
            idno_ambiguous, idno_current,
            str_ambiguous2current,
            int_current2ambiguous,
            num_ambiguous,
            )
        if bool_append_sub == True: bool_append = True

        ##
        ## append current individual to ambiguous relative
        ##
        d_idno, bool_append_sub = self.append_to_dict(
            d_idno,
            idno_current, idno_ambiguous,
            str_current2ambiguous,
            int_ambiguous2current,
            num_current,
            )
        if bool_append_sub == True: bool_append = True

    return d_idno, bool_append


def avoid_multiple_parents_and_false_predictions(
    self,
    d_idno,
    idno, idno_num1, idno_ambiguous,
    ):

    idno_parent = idno_num1
    sex_parent = d_idno[idno_parent]['sex']

    bool_continue = False

    for parent in d_idno[idno_ambiguous]['relatives']['parents']:
        ## parent already in list?
        if idno_parent == parent[0]:
            continue
        ## already a parent of same sex?
        if sex_parent == d_idno[parent[0]]['sex']:
            bool_continue = True
            break

    return bool_continue


def sub_sibling(d_vhno,vhno,num,d_idno,idno,):

    '''PI_HAT for a sibling is .5 and .25 for a half sibling'''
    ## siblings are only related to head, spouse, parent, sibling and grandchild
    ## siblings are not related to children

    bool_return, bool_append, num1, idno_num1 = self.initiate_sub_relative(
        d_vhno,vhno,num,
        )
    if bool_return == True: return d_idno, bool_append

    ##
    ## 1) append non-ambiguous sibling to *sibling*
    ##
    d_idno, bool_append_sub = self.append_non_ambiguous_relative(
        d_idno,
        idno, idno_num1, num, num1,
        'siblings',5,
        'siblings',5,
        )
    if bool_append_sub == True: bool_append = True

    ##
    ## 2) append ambiguous parent (parent of non-ambiguous sibling)
    ## to current *sibling* (if any) and vice versa
    ##
    d_idno, bool_append_sub = self.append_ambiguous_relative(
        d_vhno,vhno,d_idno,
        num,idno,
        num1,idno_num1,
        ## what is relationship of ambiguous relative to current relative
        'parents', 4,
        ## what is relationship of current relative to ambiguous relative
        'children', 3,
        ## what is relationship of non-ambiguous relative (sibling)
        ## to ambiguous relative (child)
        ## i.e. which relatives of the non-ambiguous relative are we looping over?
        3, ## nonambiguous sibling (of current sibling) is *child* of parent
        )
    if bool_append_sub == True: bool_append = True

    ##
    ## 3) append ambiguous sibling (sibling of non-ambiguous sibling)
    ## to current *sibling* (if any) and vice versa
    ##
    d_idno, bool_append_sub = self.append_ambiguous_relative(
        d_vhno,vhno,d_idno,
        num,idno,
        num1,idno_num1,
        ## what is relationship of ambiguous relative to current relative
        'siblings', 5,
        ## what is relationship of current relative to ambiguous relative
        'siblings', 5,
        ## what is relationship of non-ambiguous relative (sibling)
        ## to ambiguous relative (sibling)
        ## i.e. which relatives of the non-ambiguous relative are we looping over?
        5, ## nonambiguous sibling (of current sibling) is *sibling* of sibling
        )
    if bool_append_sub == True: bool_append = True

    return d_idno, bool_append


def loop_households(d_vhno, d_idno, d_PI_HAT,):

##        d_stats = {}

    ## loop over households
    for vhno in d_vhno.keys():

        bool_append = True

        ## keep looping over household,
        ## until no new relationships are appended
        while bool_append == True:
##            if True:

            bool_append = False

            ## loop over individuals in household
            for num in d_vhno[vhno]['num'].keys():
                idno = d_vhno[vhno]['num'][num]['idno']
                rel1 = d_vhno[vhno]['num'][num]['rel1']
                num1 = d_vhno[vhno]['num'][num]['num1']

                if idno == 2057031:
                    print(rel1)
                    stop

##                    if num1 not in d_vhno[vhno]['num'].keys():
##                        continue
##                    k = '%i_%i' %(rel1,d_vhno[vhno]['num'][num1]['rel1'],)
##                    if not k in d_stats.keys():
##                        d_stats[k] = 0
##                    d_stats[k] += 1

                ## head of household (related to self)
                if rel1 == 1:
                    d_idno, bool_append_sub = self.sub_head(d_vhno,vhno,num,d_idno,idno,)
                    if bool_append_sub == True: bool_append = True
                ## spouse (of head, spouse, child, sibling)
                elif rel1 == 2:
                    d_idno, bool_append_sub = self.sub_spouse(
                        d_vhno,vhno,num,d_idno,idno,
                        )
                    if bool_append_sub == True: bool_append = True
                ## child (of head, spouse, child, sibling)
                elif rel1 == 3:
                    d_idno, bool_append_sub = self.sub_child(
                        d_vhno,vhno,num,d_idno,idno,
                        )
                    if bool_append_sub == True: bool_append = True
                ## parent (of head, spouse)
                elif rel1 == 4:
                    d_idno, bool_append_sub = self.sub_parent(
                        d_vhno,vhno,num,d_idno,idno,
                        )
                    if bool_append_sub == True: bool_append = True
                ## sibling (of head, spouse) (never sibling of child)
                elif rel1 == 5:
                    d_idno, bool_append_sub = self.sub_sibling(
                        d_vhno,vhno,num,d_idno,idno,
                        )
                    if bool_append_sub == True: bool_append = True
                ## other
                elif rel1 == 6:
                    pass
                ## no relationship
                elif rel1 == 7:
                    pass
                ## grandchild (of head, spouse)
                elif rel1 == 8:
                    d_idno, bool_append_sub = self.sub_grandchild(
                        d_vhno,vhno,num,d_idno,idno,
                        )
                    if bool_append_sub == True: bool_append = True
                ## ???
                elif rel1 in [21,22,]:
                    pass
                elif rel1 == '':
                    pass
                else:
                    print('vhno', vhno)
                    print('idno', idno)
                    print('rel1', rel1)
                    stop_unknown_rel1

                ## end of loop over household members

            ## end of while loop

        ## end of loop over households

    ## end of function

##        print d_stats
##        stop

    return d_idno


def write_ids(d_idno,d_PI_HAT,):

    ## initiate lines
    lines = []

    ## loop over idnos
    for idno in d_idno.keys():

        gwas_id = d_idno[idno]['gwas_id']
        if gwas_id == '':
            gwas_id = 'N/A'

        ## initiate line
        line = '%8i\t' %(idno)

        line += '%21s\t' %(gwas_id)

        ##
        ## append parents to line
        ##
        idno_father = 'N/A'
        idno_mother = 'N/A'
        ab_father = 1
        ab_mother = 1
        gwas_id_father = 'N/A'
        gwas_id_mother = 'N/A'

        ## more than two parents? e.g. grandchild of head with multiple children
        ## identiy correct parent by PI_HAT if possible
        ## do not rule out false positive if no true positive
        if len(d_idno[idno]['relatives']['parents']) > 2:
            if self.verbose == True: ## tmp
                print('multiple parents', idno, d_idno[idno]['relatives']['parents'])
            gwas_id = d_idno[idno]['gwas_id']
            if gwas_id != '':
                l_parents_PI_HAT_treshold = []
                for parent in d_idno[idno]['relatives']['parents']:
                    idno_parent = parent[0]
                    gwas_id_parent = d_idno[idno_parent]['gwas_id']
                    if gwas_id_parent == '':
                        continue
                    PI_HAT = d_PI_HAT[gwas_id][gwas_id_parent]
                    if PI_HAT > 0.375: ## 0.25 if nephew/niece/uncle/aunt, 0.50 if parent/child
                        l_parents_PI_HAT_treshold += [idno_parent]
                ## delete false positives if (and only if) true positive is found
                if len(l_parents_PI_HAT_treshold) > 0:
                    for idno_parent in l_parents_PI_HAT_treshold:
                        sex_parent = d_idno[idno_parent]['sex']
                        for parent in list(d_idno[idno]['relatives']['parents']):
                            ## self
                            if idno_parent == parent[0]:
                                continue
                            ## only delete parents of same sex
                            ## i.e. keep false positives of opposite sex if no true positive
                            if sex_parent == d_idno[parent[0]]['sex']:
                                d_idno[idno]['relatives']['parents'].remove(parent)
                if len(l_parents_PI_HAT_treshold) > 2:
                    stop
            (
                idno_father, ab_father, gwas_id_father,
                idno_mother, ab_mother, gwas_id_mother,
                ) = self.get_mom_and_dad_from_list_of_parents(
                    d_idno, idno,
                    idno_father, ab_father, gwas_id_father,
                    idno_mother, ab_mother, gwas_id_mother,
                    )
        else:
            (
                idno_father, ab_father, gwas_id_father,
                idno_mother, ab_mother, gwas_id_mother,
                ) = self.get_mom_and_dad_from_list_of_parents(
                    d_idno, idno,
                    idno_father, ab_father, gwas_id_father,
                    idno_mother, ab_mother, gwas_id_mother,
                    )

        line += '%8s\t%21s\t%1i\t%8s\t%21s\t%1i\t' %(
            idno_father,
            gwas_id_father,
            ab_father,
            idno_mother,
            gwas_id_mother,
            ab_mother,
            )

        ## append siblings to line
        for sibling in d_idno[idno]['relatives']['siblings']:
            idno_sibling = sibling[0]
            gwas_id_sibling = d_idno[idno_sibling]['gwas_id']
            line += '%8s\t%21s\t%1i\t' %(
                idno_sibling,
                gwas_id_sibling,
                sibling[1],
                )

##            if len(d_idno[idno]['relatives']['siblings']) == 1 and idno_father != 'N/A' and idno_mother != 'N/A':
##                print idno, 'one sibling'
##            if len(d_idno[idno]['relatives']['siblings']) > 1 and idno_father != 'N/A' and idno_mother != 'N/A':
##                print idno, 'multiple siblings'

        ## add newline
        line = line[:-1]+'\n'

        ## append line to lines
        lines += [line]

    ## write lines
    fd = open('out/pedigree.tsv','w')
    fd.writelines(lines)
    fd.close()

    return


def get_mom_and_dad_from_list_of_parents(
    self, d_idno, idno,
    idno_father, ab_father, gwas_id_father, idno_mother, ab_mother, gwas_id_mother,
    ):

    for parent in d_idno[idno]['relatives']['parents']:
        idno_parent = parent[0]
        sex = d_idno[idno_parent]['sex']
        gwas_id_parent = d_idno[idno_parent]['gwas_id']
        if sex == 1:
            idno_father = idno_parent
            ab_father = parent[1]
            gwas_id_father = gwas_id_parent
            if gwas_id_father == '':
                gwas_id_father = 'N/A'
        elif sex == 2:
            idno_mother = idno_parent
            ab_mother = parent[1]
            gwas_id_mother = gwas_id_parent
            if gwas_id_mother == '':
                gwas_id_mother = 'N/A'

    return idno_father, ab_father, gwas_id_father, idno_mother, ab_mother, gwas_id_mother


if __name__ == '__main__':
    main()
