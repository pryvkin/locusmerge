#!/bin/bash

if [ $# -lt 4 ]; then
    echo "USAGE: $0 clusters locus_bed bam out_bed" >&2
    exit 1
fi

clusters=$1
loci=$2
bam=$3
outbed=$4

join <(sort -k4,4 $clusters) <(sort -k4,4 $loci) -1 4 -2 4 -t "$(echo -e "\t")" | \
  awk 'BEGIN{OFS="\t"} {print $11,$12,$13,$1,$14,$15}' | \
  bedtools coverage -split -counts -s -abam $bam -b - | \
  awk 'BEGIN{OFS="\t"} {$5 = $7; print}' | \
  cut -f1-6 | sort -k4,4 > $outbed

