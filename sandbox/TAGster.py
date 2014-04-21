#!python3

## Tommy Carstensen, Wellcome Trust Sanger Institute, April 2014

## the script currently does not allow multiple chromosomes to be present in a single haps file

import argparse
import collections
import os
import sys


def main():

    d_args = argparser()

    d_tell = {}
    file_matrix = d_args['out']+'.LDmatrix'
    with open(file_matrix,'w') as f_matrix:
        cnt = collections.Counter()
        for file in d_args['in']:
            print(file)
            sys.stdout.flush()
            d_LD = {}
            with open(file) as f:
                for line in f:
                    l = line.rstrip().split()
                    ID1 = l[0]
                    ID2 = l[1]
                    pos1 = int(l[2])
                    pos2 = int(l[3])
                    if abs(pos2-pos1) > d_args['window']:
                        continue
                    LD = float(l[4])
                    MAF1 = float(l[5])
                    MAF2 = float(l[6])
                    if MAF1 < d_args['min_MAF'] or MAF2 < d_args['min_MAF']:
                        continue
                    if LD < d_args['min_LD']:
                        continue
                    cnt[ID1] += 1
                    cnt[ID2] += 1
                    for IDa,IDb in ((ID1,ID2),(ID2,ID1)):
                        try:
                            d_LD[IDa].append(IDb)
                        except KeyError:
                            d_LD[IDa] = [IDb]
            for ID,l_IDs in d_LD.items():
                tell = f_matrix.tell()
                f_matrix.write('%s\n' %(' '.join(l_IDs)))
                try:
                    d_tell[ID].append(tell)
                except KeyError:
                    d_tell[ID] = [tell]
            
##            break ## tmp!!!

        del d_LD

    print('len(setS)', len(d_tell.keys()), len(cnt.items()))
    setT = set()
    setQ = set()
    if d_args['preselected']:
        with open(d_args['preselected']) as f:
            for line in f:
                ID = line.rstrip()
                setT.add(ID)
                setQ.append(allinLD)
                setP.remove(ID)
                del cnt[ID]
                for ID2 in d_LD[ID]:
                    d_LD[ID2].remove(ID)
                    cnt[ID2] -= 1
                del d_LD[ID]
##                setT = set(f.read().rstrip().split('\n'))

    d_cnt2IDs = {}
    for ID,count in cnt.items():
        try:
            d_cnt2IDs[count].append(ID)
        except KeyError:
            d_cnt2IDs[count] = [ID]
##    del cnt

    import time
    T1 = time.time()
    d_time = {i:0 for i in range(10)}
    most_common_key = max(d_cnt2IDs.keys())
    with open(file_matrix) as f:
        while len(setT) < d_args['max_tagSNP']:
            T2 = time.time()
            step = 0
            t1 = time.time()
            ## Find the most common element,
            ## which is the rate limiting step in this loop.
            ##
            ## Algorithm could be ~10 times faster
            ## by only finding the 10 most common elements every 10th turn,
            ## which would not be 100% correct
            ##
            ## Alternatively algorithm could be made faster by
            ## using count as key in a dict and sets of IDs as value.
            ## The time complexity of set.add() and set.remove() is O(1)
            t2 = time.time(); d_time[step] += t2-t1; step += 1; t1 = t2
            while True:
                print('a',most_common_key)
                try:
                    ID = d_cnt2IDs[most_common_key].pop()
                    break
                except IndexError:
                    del d_cnt2IDs[most_common_key]
                    most_common_key = max(d_cnt2IDs.keys())
                    print('b',most_common_key)
            ## 10 times faster but not 100% correct solution:
##            if not most_common:
##                most_common = list(reversed(cnt.most_common(10)))
##                if not most_common:
##                    break
##            ID = most_common.pop()[0]
            t2 = time.time(); d_time[step] += t2-t1; step += 1; t1 = t2
            ## Break loop if cnt is exhausted.
            print(len(setT),ID,most_common_key,round((T2-T1)/(len(setT)+1),3),len(d_tell[ID]))
            print(len(setT),ID,most_common_key,len(d_tell[ID]))
            t2 = time.time(); d_time[step] += t2-t1; step += 1; t1 = t2
            setT |= set([ID])
            del cnt[ID]
            t2 = time.time(); d_time[step] += t2-t1; step += 1; t1 = t2
            ## File seek is the second slowest step in this loop.
            for seek in d_tell[ID]:
                f.seek(seek)
                line = f.readline()
                for ID2 in line.rstrip().split():
                    ## tag SNP
                    if ID2 in setT:
                        continue
                    ## tagged
                    count = cnt[ID2]
                    d_cnt2IDs[count].remove(ID2)
                    try:
                        d_cnt2IDs[count-1].append(ID2)
                    except KeyError:
                        d_cnt2IDs[count-1] = [ID2]
                    cnt[ID2] -= 1
                    if cnt[ID2] < 0:
                        print(cnt[ID2])
                        print(ID,ID2)
                        stoptmp
                    setQ.add(ID2) ## tmp!!!
            t2 = time.time(); d_time[step] += t2-t1; step += 1; t1 = t2
            del d_tell[ID]
            t2 = time.time(); d_time[step] += t2-t1; step += 1; t1 = t2
            if len(setT) > 1000: break ## tmp!!!
            sys.stdout.flush()

    print('d_time',d_time)

    with open(d_args['out']+'.tagSNPs','w') as f:
        for ID in setT:
            f.write('%s\n' %(ID))
    with open(d_args['out']+'.tagged','w') as f:
        for ID in setQ:
            f.write('%s\n' %(ID))

    ## Clean up large temporary files.
    os.remove(file_matrix)

    print('len(setT)',len(setT))
    print('len(setQ)',len(setQ))
    print('len(setS)', len(d_tell.keys()), len(cnt.items()))

    return


def argparser():

    parser = argparse.ArgumentParser()

    parser.add_argument('--in', required=True, nargs='+')
    parser.add_argument('--out', required=True)
    parser.add_argument('--min_LD', type=float, default=.8)
    parser.add_argument('--min_MAF', type=float, default=.01)
    parser.add_argument('--window', type=int, default=250000)
    parser.add_argument('--bool_verbose', '--verbose', action='store_true',)
    parser.add_argument('--preselected')
    parser.add_argument('--max_tagSNP', type=int, default=1000000)

    d_args = vars(parser.parse_args())

    return d_args


if __name__ == '__main__':
    main()
