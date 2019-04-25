#!/bin/bash

#Usage: mapping.sh sampleDir resultsDir genomeDir
# input arguments: sample directory, name of sample, and directory in which the reads are stored
sampleDir=$1

resultsDir=$2

genomeDir=$3

if [ ! -d "$sampleDir" ]; then
    echo "Samples directory does not exist";
    exit
fi

if [ -e tmp ]; then
    rm tmp
fi

#ls -v $sampleDir/*.gz | sed -n '1~4p' > R1
#ls -v $sampleDir/*.gz | sed -n '2~4p' > R2
#ls -v $sampleDir/*.gz | sed -n '3~4p' > R3
#ls -v $sampleDir/*.gz | sed -n '4~4p' > R4

#paste R1 R2 R3 R4 > tmp


# 0. Trim quality
# java -jar /usr/local/trimmomatic/trimmomatic-0.36.jar SE -threads 7 -trimlog trim.log ~/Eutrema/fieldSamples/$2.fastq $resultsDir/1pass/$ID/$2_qual30.fastq LEADING:3 TRAILING:3 SLIDINGWINDOW:3:30 MINLEN:36

# 1. 1pass
#generate reference genome
#STAR --runMode genomeGenerate --genomeDir $genomeDir --sjdbGTFfile $genomeDir/*.gtf --genomeFastaFiles $genomeDir/*.fa --runThreadN 5
    
for file in `ls -1 $1/*.gz`; do
    #get names of files for naming!
    ID=$(basename -- "$file")
    ID="${ID%_allMe*}"
    echo $ID
    # make directories for analysis
    mkdir -p $resultsDir/1pass/$ID/ \
             $resultsDir/2pass/$ID/ \
             $resultsDir/2pass/$ID/genome/ \
             $resultsDir/2pass/$ID/outfiles/
    echo "1pass"
    STARlong --genomeDir $genomeDir \
             --readFilesIn $file \
             --readFilesCommand zcat \
             --runThreadN 10 \
             --limitBAMsortRAM 20000000000 \
             --outSAMtype BAM SortedByCoordinate \
             --outFileNamePrefix $resultsDir/1pass/$ID/${ID}_ 2>&1 $resultsDir/2pass/$ID/outfiles/0-1pass.out
    echo "done"
    echo "2pass"
    # 2. 2pass
    echo "Genome gen"
    STAR --runMode genomeGenerate \
         --genomeDir $resultsDir/2pass/$ID/genome/\
         --sjdbGTFfile $genomeDir/*.gtf \
         --genomeFastaFiles $genomeDir/*.fa \
         --sjdbFileChrStartEnd $resultsDir/1pass/$ID/${ID}_SJ.out.tab \
         --sjdbOverhang 75 \
         --runThreadN 5 2>&1 $resultsDir/2pass/$ID/outfiles/1-2pass.out

    echo "2pass map"
    STARlong --genomeDir $resultsDir/2pass/$ID/genome/ \
             --readFilesIn $file \
             --readFilesCommand zcat \
             --runThreadN 10 \
             --limitBAMsortRAM 20000000000 \
             --outSAMtype BAM SortedByCoordinate \
             --outFileNamePrefix $resultsDir/2pass/$ID/${ID}_ 2>&1 $resultsDir/2pass/$ID/outfiles/2-2passalign.out
    echo "done"
done 

