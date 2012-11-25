#!/software/bin/python

## T. Carstensen (tc9), M.S. Sandhu (ms23), D. Gurdasani (dg11)
## Wellcome Trust Sanger Institute, 2012

import os, sys

def main():

    ## if two samples are related, then only one of them is removed
    ## but which one is removed is random
    ## i.e. the first occurence among sorted samples is removed

    bfile, d_IID2IIDs = parse_genome()

    d_count2IIDs = count_relations(d_IID2IIDs)

    l_exclude = generate_exclusion(d_IID2IIDs,d_count2IIDs,bfile=bfile,)

    write_exclusion_list(l_exclude)

    return


def write_exclusion_list(l_exclude,):

    pi_hat_max = float(sys.argv[sys.argv.index('--pi_hat_max')+1])
    bfile = sys.argv[sys.argv.index('--bfile')+1]

    fd = open('%s.genome%.2f.samples' %(bfile,pi_hat_max,),'w')
    fd.write('\n'.join(['%s\t%s' %(IIDa,IIDa,) for IIDa in l_exclude]))
    fd.close()

    return


def generate_exclusion(d_IID2IIDs,d_count2IIDs,verbose=False,bfile=None,):

    l_exclude = []
    while d_count2IIDs != {}:
        count_max = max(d_count2IIDs.keys())
        if count_max == 0:
            stop
            break
        l_IIDas = list(d_count2IIDs[count_max])
        ## sort so first occurence(s) are deleted first
        l_IIDas.sort()
        for IIDa in l_IIDas:
##            print IIDa, count_max
            l_IIDbs = list(d_IID2IIDs[IIDa])
            if verbose == True:
                print count_max, IIDa, l_IIDbs
                pass
            for IIDb in l_IIDbs:
                if verbose == True:
                    print count_max, IIDa,IIDb
                    pass
                count = len(d_IID2IIDs[IIDb])
                d_IID2IIDs[IIDa].remove(IIDb)
                d_IID2IIDs[IIDb].remove(IIDa)
                d_count2IIDs[count].remove(IIDb)
                if count > 1:
                    if not count-1 in d_count2IIDs.keys():
                        d_count2IIDs[count-1] = []
                        pass
                    d_count2IIDs[count-1] += [IIDb]
                    pass
                if d_count2IIDs[count] == []:
                    del d_count2IIDs[count]
                    pass
                continue
            d_count2IIDs[count_max].remove(IIDa)
            del d_IID2IIDs[IIDa]
            l_exclude += [IIDa]
            break
        if d_count2IIDs[count_max] == []:
            del d_count2IIDs[count_max]
            pass
        continue

    l_exclude.sort()

    return l_exclude


def count_relations(d_IID2IIDs):

    d_count2IIDs = {}
    for IID,l_IIDs in d_IID2IIDs.items():
        count = len(l_IIDs)
        if not count in d_count2IIDs.keys():
            d_count2IIDs[count] = []
            pass
        d_count2IIDs[count] += [IID]

    return d_count2IIDs


def parse_genome():

    bfile = sys.argv[sys.argv.index('--bfile')+1]
    pi_hat_max = float(sys.argv[sys.argv.index('--pi_hat_max')+1])

    fp_genome = '%s.genome' %(bfile)
    fd = open(fp_genome,'r')
    d_IID2IIDs = {}
    ## skip header
    for line in fd: break
    ## loop lines
    for line in fd:
        l = line.split()
        IID1 = l[1]
        IID2 = l[3]
        try:
            PI_HAT = float(l[9])
        except:
            print line
            stop
        if PI_HAT < pi_hat_max:
            continue
        for IIDa,IIDb in [
            [IID1,IID2,],
            [IID2,IID1,],
            ]:
            if not IIDa in d_IID2IIDs.keys():
                d_IID2IIDs[IIDa] = []
                pass
            d_IID2IIDs[IIDa] += [IIDb]
    fd.close()

    return bfile, d_IID2IIDs


def unit_test():

    d_IID2IIDs = {
        'a':['b','c','d','e','f',],
        'b':['a','d','e','g','h',],
        'c':['a','d','i',],
        'd':['a','b','c','h',],
        'e':['a','b','j',],
        'f':['a',],
        'g':['b',],
        'h':['b','d',],
        'i':['c',],
        'j':['e',],
        'k':['l',],
        'l':['k',],
        'm':['n','o',],
        'n':['m','o',],
        'o':['m','n',],
        }
    d_count2IIDs = count_relations(d_IID2IIDs)
    l_exclude = generate_exclusion(d_IID2IIDs,d_count2IIDs,verbose=True,)
    print l_exclude == ['a', 'b', 'c', 'd', 'e', 'k', 'm', 'n']

    d_IID2IIDs = {
        'a':['b',],
        'b':['a',],
        }
    d_count2IIDs = count_relations(d_IID2IIDs)
    l_exclude = generate_exclusion(d_IID2IIDs,d_count2IIDs,verbose=True,)
    print l_exclude == ['a',]

    d_IID2IIDs = {
        'a':['b','c',],
        'b':['a','c',],
        'c':['a','b',],
        }
    d_count2IIDs = count_relations(d_IID2IIDs)
    l_exclude = generate_exclusion(d_IID2IIDs,d_count2IIDs,verbose=True,)
    print l_exclude == ['a','b',]

    return


if __name__ == '__main__':
    main()
##    unit_test()
