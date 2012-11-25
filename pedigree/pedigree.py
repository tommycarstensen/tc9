## Tommy Carstensen, Wellcome Trust Sanger Institute, May-June 2012
## GNU General Public License

## This script will generate the following output as agreed on Friday 18May2012
## FatherID Ab MotherID Ab SiblingID1 Ab1 SiblingID2 Ab2 ...
## ambiguous = 1, not ambiguous = 0

## built-ins
import os

class pedigree():


    def main(self,):

        lines, l_headers = self.read_lines_and_header()

        ## loop lines prior to looping households and individuals
        ## in order to have all information about a household and an individual
        ## prior to looping over them
        ## i.e. do unilateral relationships
        d_vhno, d_idno = self.loop_lines(lines, l_headers,)

        self.calculate_stats_for_deepti(d_idno,d_vhno,)

        ## loop genetic data to
        ## rule out false positives
        ## e.g. decide between two possible parents of same sex
        d_PI_HAT = self.read_gene_comparison_fast(self.treshold)
##        d_PI_HAT = self.read_gene_comparison_slow(self.treshold)

        ## exclude errors *prior* to building all relationship
        d_vhno, d_idno = self.find_errors_and_exclude_from_dataset1(d_vhno,d_idno,)

        ## identify non-ambiguous and ambiguous relationships
        ## i.e. do bilateral and extended relationships
        d_idno = self.loop_households(d_vhno, d_idno, d_PI_HAT,)

        ## check errors *after* building relationships
        d_idno = self.find_errors_and_exclude_from_dataset2(d_idno,)

        self.calculate_stats_for_deepti(d_idno,d_vhno,)

        ## write results to a table
        ## N.B. Individuals with multiple parents are skipped
        ## N.B. None of these skipped individuals have gwas_ids
        self.write_ids(d_idno,d_PI_HAT,)

        ## check that parents are not younger than children etc.
        ## *after* building relationships
        self.check_age2(d_idno,)

        return


    def calculate_stats_for_deepti(self,d_idno,d_vhno,):

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

        print 'father', count_father, 'mother', count_mother, 'both', count_both

        return


    def check_age2(self,d_idno,):

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
                            print
                            print idno, idno_relative
                            print agec, agec_relative
                            print 'spouse'
                            stop_spouse
                    elif relation == 'siblings':
                        if abs(agec_relative-agec) > (68-16):
                            print
                            print agec, agec_relative
                            print idno, idno_relative
                            print 'sibling'
                            stop_sibling
                    elif relation == 'grandparents':
                        if agec > agec_relative or agec_relative-agec > (106-15):
                            print
                            print agec_relative, agec
                            print idno_relative, idno
                            stop_grandparent
                    elif relation == 'grandchildren':
                        if agec < agec_relative or agec-agec_relative > (106-15):
                            print
                            print agec, agec_relative
                            print idno_relative, idno
                            stop_grandchild
                    elif relation == 'parents':
                        if (
                            agec > agec_relative
                            or
                            (ambiguity == 0 and agec_relative-agec > (81-19))
                            or
                            (ambiguity == 1 and agec_relative-agec > (88-15))
                            ):
                            print
                            print ambiguity
                            print agec_relative, agec
                            print idno_relative, idno
                            print 'parent'
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
                            print
                            print agec_relative, agec
                            print idno_relative, idno
                            print 'child'
                            stop_children
                            if agec < agec_relative:
                                stop_children
                    else:
                        print relation, agec, agec_relative
                        stop_new_relation

        return


    def read_gene_comparison_fast(self,treshold,):

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


    def append_to_d_PI_HAT(self,d_PI_HAT,IID1,IID2,PI_HAT,):
    
        if not IID1 in d_PI_HAT.keys():
            d_PI_HAT[IID1] = {}
        d_PI_HAT[IID1][IID2] = PI_HAT
        ## do reverse (which is slow...)
        if not IID2 in d_PI_HAT.keys():
            d_PI_HAT[IID2] = {}
        d_PI_HAT[IID2][IID1] = PI_HAT

        return d_PI_HAT


    def read_gene_comparison_slow(self,treshold,):

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


    def find_errors_and_exclude_from_dataset2(self,d_idno,):

        return d_idno


    def read_lines_and_header(self,):

        fd = open('in/reldata_allgpc_90212.tsv','r')
        lines = fd.readlines()
        fd.close()

        l_headers = lines[0].strip().split('\t')

        return lines, l_headers


    def line2parameters(self,line,l_headers,):
        
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


    def loop_lines(self, lines, l_headers,):

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


    def redflag_household_multiple_heads(self,d_vhno,vhno,):

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


    def check_age(self, d_vhno, vhno,):

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
                    print vhno, num, num1
                    print agec, agec_num1
                    stopage1
            elif rel1 == 2:
                ## no more than xxx years between spouses
                if abs(agec-agec_num1) > abs(22-87):
                    print vhno, num, num1
                    print agec, agec_num1
                    stopage2
            elif rel1 == 3:
                ## no more than xxx years between child and parent
                if agec-agec_num1 > 73-15:
                    print vhno, num, num1
                    print agec, agec_num1
                    stopage3
            elif rel1 == 4:
                ## no more than xxx years between parent and child
                if agec-agec_num1 > 85-30:
                    print vhno, num, num1
                    print agec, agec_num1
                    stopage4
            elif rel1 == 5:
                ## no more than xxx years between siblings
                if abs(agec-agec_num1) > abs(68-16):
                    print vhno, num, num1
                    print agec, agec_num1
                    stopage5
            elif rel1 in [6,7,]:
                pass
            elif rel1 == 8:
                ## no more than xxx years between grandchild and parent
                if agec-agec_num1 > 93-13:
                    print vhno, num, num1
                    print agec, agec_num1
                    stopage8
            elif rel1 in [21,22,]:
                pass
            elif rel1 == '':
                pass
            else:
                print rel1, agec, agec_num1
                stop

        return


    def print_removal(self,vhno,num,rel1,num1,):

        if rel1 == 3:
            print 'temporarily removing vhno %i from dataset \
because individual %i is child of self' %(vhno,num,)
        elif rel1 == 2:
            print 'temporarily removing vhno %i from dataset \
because individual %i is spouse of self' %(vhno,num,)
        elif rel1 == 5:
            print 'temporarily removing vhno %i from dataset \
because individual %i is sibling of self' %(vhno,num,)
        elif rel1 == 8:
            print 'temporarily removing vhno %i from dataset \
because individual %i is grandchild of self' %(vhno,num,)
        else:
            print 'relationship to self', vhno, num, rel1, num1
            stop

        return


    def error_self_relation(self,d_vhno,vhno,):

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


    def find_errors_and_exclude_from_dataset1(self,d_vhno,d_idno,):

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
            print 'removed %i of %i households (%.1f%%) because of data errors' %(
                count_households1-count_households2,
                count_households1,
                100*(count_households1-count_households2)/float(count_households1)
                )

        return d_vhno, d_idno


    def sub_head(self,d_vhno,vhno,num,d_idno,idno,):

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


    def sub_spouse(self,d_vhno,vhno,num,d_idno,idno,):

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


    def initiate_sub_relative(self,d_vhno,vhno,num,):

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


    def sub_grandchild(self,d_vhno,vhno,num,d_idno,idno,):

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


    def sub_child(self,d_vhno,vhno,num,d_idno,idno,):

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
            print idno_num1, d_idno[idno_num1]['rel1'].keys()
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


    def sub_parent(self,d_vhno,vhno,num,d_idno,idno,):

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
            print [int_ambiguous2current,int_current2ambiguous,int_nonambiguous2ambiguous]
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


    def sub_sibling(self,d_vhno,vhno,num,d_idno,idno,):

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


    def loop_households(self, d_vhno, d_idno, d_PI_HAT,):

##        d_stats = {}

        ## loop over households
        for vhno in d_vhno.keys():

##            if not vhno in [20382]: ## tmp
##                continue

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
                        print rel1
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
                        print 'vhno', vhno
                        print 'idno', idno
                        print 'rel1', rel1
                        stop_unknown_rel1

                    ## end of loop over household members

                ## end of while loop

            ## end of loop over households

        ## end of function

##        print d_stats
##        stop

        return d_idno


    def write_ids(self,d_idno,d_PI_HAT,):

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
                    print 'multiple parents', idno, d_idno[idno]['relatives']['parents']
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


    def __init__(self,):

        self.verbose = False
        ## PI_HAT treshold
        self.treshold = 0.05

        return


if __name__ == '__main__':
    instance_pedigree = pedigree()
    instance_pedigree.main()
