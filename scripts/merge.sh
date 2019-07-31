#!/bin/bash

DIR=~/scratch/PS
READS=$1
OUT=$2
SeqP=/usr/local/seqprep

if [ -e tmp ]; then
    rm tmp
fi

ls -v $READS/*P* | sed -n 'p;n' > R1
ls -v $READS/*P* | sed -n 'n;p' > R2

paste R1 R2  > tmp

while read for rev; do
    NAME=$(basename -- "$for")
    NAME="${NAME%.fq.gz}"

    $SeqP/SeqPrep -f $for -r $rev -1 $OUT/${NAME}_1U.fq.gz \
        -2 $OUT/${NAME}_2U.fq.gz -s $OUT/${NAME}_merged.fq.gz
done < tmp

rm tmp R1 R2

