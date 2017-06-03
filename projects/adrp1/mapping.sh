bin=/software/gapi/pkg/star/2.5.2b/bin/STAR

samtools=/software/hgi/pkglocal/samtools-1.3/bin/samtools
tabix=/software/hgi/pkglocal/htslib-1.3/bin/tabix

build=$1
if [ $build -eq 37 ]; then
affix=37.75
else
affix=38.83
fi

pass=$2
if [ $pass -ne 1 -a $pass -ne 2 ]; then echo pass $pass; exit; fi

sampleID=$3
outFileNamePrefix=out_STAR/$affix/pass$pass/$sampleID/

if [ -d $outFileNamePrefix ]; then exit; fi
mkdir -p $outFileNamePrefix
chmod 777 $outFileNamePrefix

files1=$4
files2=$5
runThreadN=$6

cmd="$bin \
 --runMode alignReads \
 --runThreadN $runThreadN \
 --genomeDir genome_index_$affix \
 --readFilesIn $files1 $file2 \
 --readFilesCommand zcat \
 --outFileNamePrefix $outFileNamePrefix \
 --outSAMtype BAM SortedByCoordinate \
 --quantMode GeneCounts \
"

echo $cmd

if [ $pass -eq 2 ]; then
## limitSjdbInsertNsj
## maximum number of junction to be inserted to the genome on the fly
## at the mapping stage, including those from annotations and those detected in
## the 1st step of the 2-pass run
cmd=$cmd" --limitSjdbInsertNsj 3000000 "
echo $cmd
cmd=$cmd" --sjdbFileChrStartEnd out_STAR/$affix/pass1/*/SJ.out.tab "
fi

eval $cmd

$samtools index $outFileNamePrefix/Aligned.sortedByCoord.out.bam

