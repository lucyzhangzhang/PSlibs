#!/bin/bash
#Usage: catMerge.sh <in> <out>
if [ -e tmp ]; then
    rm tmp
fi

ls -v $1/ps* | sed -n '1~6p' > R1
ls -v $1/ps* | sed -n '2~6p' > R2
ls -v $1/ps* | sed -n '3~6p' > R3
ls -v $1/ps* | sed -n '4~6p' > R4
ls -v $1/ps* | sed -n '5~6p' > R5
ls -v $1/ps* | sed -n '6~6p' > R6
#ls -v $1/ps* | sed -n '7~10p' > R7
#ls -v $1/ps* | sed -n '8~10p' > R8
#ls -v $1/ps* | sed -n '9~10p' > R9
#ls -v $1/ps* | sed -n '10~10p' > R10

paste R1 R2 R3 R4 R5 R6  > tmp

while read a b c d e f; do
    NAME=$( basename -- "$a" )
    NAME="${NAME%_L00*}"

    echo $NAME
    cat $a $b $c $d $e $f > $2/${NAME}_allMerged.fq.gz
done < tmp
rm tmp R*

