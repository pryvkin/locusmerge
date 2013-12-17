#!/bin/bash

if [ $# -lt 2 ]; then
    echo "USAGE: $0 bam outprefix" >&2
    exit 1
fi

bam=$1
outpre=$2

chrlen="$(mktemp)"

# take chr lengths from BAM hdr (needed by bedtools genomecov)
echo "Reading BAM header..." >&2
samtools view -H $bam | \
  awk '/^@SQ/ {
    sub(/^SN:/, "", $2)
    sub(/^LN:/, "", $3)
    print $2"\t"$3 }' > $chrlen

# compute genome coverage and output as a bedgraph
echo "Computing coverage on '+' strand..." >&2
bedtools genomecov -ibam $bam -g $chrlen -split -bg -strand "+" > ${outpre}.pos.bg
echo "Computing coverage on '-' strand..." >&2
bedtools genomecov -ibam $bam -g $chrlen -split -bg -strand "-" > ${outpre}.neg.bg

rm -f $chrlen

