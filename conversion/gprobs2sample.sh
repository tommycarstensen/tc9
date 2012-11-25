$CHROMOSOME = 22; \
## parse header of BEAGLE gprobs file
head -1 out_BEAGLE/sep/BeagleOutput.$CHROMOSOME.bgl.ProduceBeagleInput.$CHROMOSOME.bgl.gprobs | \
## replace blank space with newline globally
sed 's/ /\n/g' | \
## exclude line with text marker
fgrep -v "marker" | \
## exclude line with text allele
fgrep -v "allele" | \
## remove duplicate lines (individual IDs)
uniq | \
## duplicate column 1 and append a 3rd NA column
awk '{print $1, $1, "NA"}' > \
## send to sample file
chromosome$CHROMOSOME_body.sample;

## append header
$CHROMOSOME = 22; \
echo "ID_1 ID_2 missing\n0 0 0" > header.txt; cat header.txt chromosome$CHROMOSOME_body.sample > chromosome$CHROMOSOME.sample; rm header.txt chromosome$CHROMOSOME_body.sample
