#!/bin/bash

if [ $# -lt 2 ]; then
    echo "USAGE: $0 in_sortedbyname.bam loci.bed" >&2
    exit 1
fi

maxhits=5

./bam_filter_by_nh $1 $maxhits /dev/stdout | \
  bedtools intersect -s -wo -abam - -b $2 | \
  awk '{print $4"\t"$16"\t"($3-$2)}' | \
  ./locusmerge_pairs - | sort | uniq -c | awk '{print $2"\t"$3"\t"$1}' | \
  ./locusmerge_similarity -
