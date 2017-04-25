bin=/software/gapi/pkg/star/2.5.2b/bin/STAR

samtools=/software/hgi/pkglocal/samtools-1.3/bin/samtools
tabix=/software/hgi/pkglocal/htslib-1.3/bin/tabix

#build=38
#bsub -q long -G adrp1 -R 'select[mem>37900] rusage[mem=37900]' -M37900 -n24 -R'span[hosts=1]' -k "$(pwd)/checkpoint$build method=blcrkill 600" cr_run bash mapping.sh $build 
#bsub -q long -G ug2g -R 'select[mem>37900] rusage[mem=37900]' -M37900 -n24 -R'span[hosts=1]' bash mapping.sh $build

# --genomeDir genome_index_37.75 \
# --genomeFastaFiles Homo_sapiens.GRCh37.75.dna.primary_assembly.fa \
# --sjdbGTFfile Homo_sapiens.GRCh37.75.gtf \

# just request 12 cores...

build=$1
if [ $build -eq 38 ]; then
affix=38.83
else
affix=37.75
fi

pass=$2
if [ $pass -ne 1 -a $pass -ne 2 ]; then echo pass $pass; exit; fi

for sampleID in $(cat ../../../metadata/RNA_imeta.txt | cut -f3 | sort -u); do

laneIDs=$(cat ../../../metadata/RNA_imeta.txt | awk -v sampleID=$sampleID '$3==sampleID{print $1}')

laneID1=$(echo $laneIDs | cut -d" " -f1)
laneID2=$(echo $laneIDs | cut -d" " -f2)

outFileNamePrefix=out_STAR/$affix/pass$pass/$sampleID/

if [ -d $outFileNamePrefix ]; then continue; fi
mkdir -p $outFileNamePrefix
chmod 777 $outFileNamePrefix

cmd="$bin \
 --runMode alignReads \
 --runThreadN 8 \
 --genomeDir genome_index_$affix \
 --readFilesIn \
  split_fasta/$laneID1/1.fasta.gz,split_fasta/$laneID2/1.fasta.gz \
  split_fasta/$laneID1/2.fasta.gz,split_fasta/$laneID2/2.fasta.gz \
 --readFilesCommand zcat \
 --outFileNamePrefix $outFileNamePrefix \
 --outSAMtype BAM Unsorted SortedByCoordinate \
 --quantMode GeneCounts \
"

if [ $pass -eq 2 ]; then
cmd=$cmd&" --sjdbFileChrStartEnd out_STAR/$affix/pass1/[HGN]*/SJ.out.tab \ "
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

done

