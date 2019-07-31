#!/bin/bash

DIR=~/scratch/PS
RAW=$DIR/reads
OUT=$DIR/trimmed

ADAPTER=$DIR/trimmed/adapter.fa

if [ -e tmp ]; then
    rm tmp
fi

ls -v $RAW/*.gz | sed -n 'p;n' > R1
ls -v $RAW/*.gz | sed -n 'n;p' > R2

paste R1 R2  > tmp

while read for rev; do
    #get names of files for naming!
    echo "Trimming $for and $rev"
    NAME=$(basename -- "$for")
    NAME="${NAME%_R1_001.fastq.gz}"

#trimming, TruSeq adapters
java -jar /usr/local/trimmomatic/trimmomatic-0.36.jar PE -phred33 -threads 4 -trimlog $OUT/${NAME}_trim.log \
    $for $rev -baseout $OUT/$NAME.fq.gz ILLUMINACLIP:/usr/local/trimmomatic/adapters/TruSeq2-PE.fa:2:30:10:5:TRUE \
    LEADING:3 TRAILING:3 SLIDINGWINDOW:3:30 MINLEN:36
    echo "done"
done < tmp

rm tmp R1 R2
