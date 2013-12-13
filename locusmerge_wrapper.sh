#!/bin/bash

if [ $# -lt 4 ]; then
  echo "USAGE: byname.bam loci.bed out.sim tmpdir" >&2
  exit 1
fi

export TMPDIR=$4

bash ~/code/locusmerge/locusmerge.sh $1 $2 > $3
