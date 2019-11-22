#!/bin/bash
#run Q0RTs
#
# Usage: QoRTs.sh <2pass_loc> <out_dir>

BAMDIR=$1
QoDIR=/home/lucy/bin/hartleys-QoRTs-099881f
GTF=/home/lucy/scratch/v1.0/*gtf

for I in `ls -1 $1/*/*sorted*.bam`; do
    NAME=$(basename -- "$I" )
    NAME="${NAME%_Aligned*}"
    mkdir $2/$NAME

    java -jar $QoDIR/QoRTs.jar QC --prefilterImproperPairs --maxReadLength 301 --maxPhredScore 64 $I $GTF $2/$NAME
done
