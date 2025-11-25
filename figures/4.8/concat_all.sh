#!/bin/bash
in="$1"
out="all_results.tsv"
outfile="$in.$out"



# Process first file separately to get header
first_file=$(find "$in" -maxdepth 1 -type f | head -1)
if [ -n "$first_file" ]; then
    head -2 "$first_file" > "$in.$out"
    echo "" >> "$in.$out"
fi

# Process remaining files (skip first file, only take line 2)
find "$in" -maxdepth 1 -type f ! -path "$first_file" -exec awk 'FNR==2 {print; nextfile}' {} + >> "$in.$out"

GTA=$(grep "gta" "$in.$out" | wc -l)
STAR=$(grep "star" "$in.$out" | wc -l)
MINI=$(grep "minimap" "$in.$out" | wc -l)
BWA=$(grep "bwa" "$in.$out" | wc -l)
BOWTIE=$(grep "bowtie2" "$in.$out" | wc -l)
HISAT=$(grep "hisat" "$in.$out" | wc -l)


echo "gta: $GTA"
echo "start: $STAR"
echo "minimap: $MINI"
echo "hisat: $HISAT"
echo "bwa: $BWA"
echo "bowtie2: $BOWTIE"
