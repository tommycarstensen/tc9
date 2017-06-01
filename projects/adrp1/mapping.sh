bin=/software/gapi/pkg/star/2.5.2b/bin/STAR

samtools=/software/hgi/pkglocal/samtools-1.3/bin/samtools
tabix=/software/hgi/pkglocal/htslib-1.3/bin/tabix

# --genomeDir genome_index_37.75 \
# --genomeFastaFiles Homo_sapiens.GRCh37.75.dna.primary_assembly.fa \
# --sjdbGTFfile Homo_sapiens.GRCh37.75.gtf \

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

cmd="$bin \
 --runMode alignReads \
 --runThreadN 12 \
 --genomeDir genome_index_$affix \
 --readFilesIn $files1 $file2 \
 --readFilesCommand zcat \
 --outFileNamePrefix $outFileNamePrefix \
 --outSAMtype BAM Unsorted SortedByCoordinate \
 --quantMode GeneCounts \
"

if [ $pass -eq 2 ]; then
cmd=$cmd&" --sjdbFileChrStartEnd out_STAR/$affix/pass1/*/SJ.out.tab \ "
fi

echo $cmd
eval $cmd

# --outSAMtype BAM Unsorted SortedByCoordinate \

# --readFilesIn fasta/*.fasta.gz \

# --readFilesIn \
#  $(ls split_fasta/*/1.fasta.gz | sort -V | tr "\n" ",") \
#  $(ls split_fasta/*/2.fasta.gz | sort -V | tr "\n" ",") \

# --readFilesIn bam/*#[1-9]*.bam \

#$samtools sort $outFileNamePrefix/Aligned.out.sorted.bam
$samtools index $outFileNamePrefix/Aligned.sortedByCoord.out.bam

