import os, sys
sys.path.append('/nfs/users/nfs_t/tc9/github/tc9/misc')
import gnuplot


bim4 = 'omni2.5-4_20120904_agv_gtu.bim'
bim8 = 'omni2.5-8_agv_20120910_gtu.bim'

strand8 = 'HumanOmni2.5-8v1_A-b37.strand'
strand437 = 'HumanOmni2.5M-b37-v2.strand'
strand436 = 'HumanOmni2.5M-b37-v2.strand'


gnuplot.venn3(
    i1=1,i2=1,i3=1,i4=1,i5=1,i6=1,i7=1,
    text1='%s %i' %(strand436,2449906,),
    text2='%s %i' %(strand437,2449626,),
    text3='%s %i' %(strand8,2379514,),
    sum1 = 280, sum2 = 0, sum3 = 2376802,
    sum4 = 
    )

    
d_count = {}
for fn in (bim4,bim8,strand436,strand437,strand8):
    cmd = 'cat %s | wc -l' %(fn)
    i = int(os.popen(cmd).read())
    d_count[fn] = i

## show that strand files are subsets of bim files
for fn1,fn2 in [
    [strand436,bim4,],
    [strand437,bim4,],
    [strand8,bim8,],
    ]:
    cmd = 'comm -12 %s.sorted.id.id %s.sorted.id.id | wc -l' %(fn1,fn2,)
    i = int(os.popen(cmd).read())
    print fn1,fn2,i
    d_count[fn1+fn2] = i

    gnuplot.venn2(
        i1=d_count[fn1],i2=d_count[fn2],i3=i,
        text1=bim4,text2=strand436,
        )
    stop
