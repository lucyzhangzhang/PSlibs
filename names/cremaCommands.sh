#!/bin/bash

#output generated from /home/lucy/R/Eutrema/PS/scripts/PS_pairwise.R
NAME_DIR=/home/lucy/R/Eutrema/PS/names
CREMA=/home/lucy/bin/crema
PROC=3
cd $NAME_DIR
ls -1 */*/*/names > names


#has everything in it
conda activate renv

#generate multi fasta

#go on the cluster
#ssh lucy@114
while read line; do
    xargs samtools faidx /2/scratch/lucy/misc/2transcripts.fa < $line > $line.fa
done < names

#exit the cluster
#exit

#run cpat
conda deactivate

conda activate cpat

parallel -j $PROC cpat.py -g {}.fa -o {}.cpat -x $CREMA/cpat_models/ath_hexamer.txt -d $CREMA/cpat_models/ath_logit.RData :::: names

conda deactivate

#run diamond
parallel -j $PROC diamond blastx -d swissprot.dmnd -q {}.fa -o {}.diamond \
    -e 0.001 -k 5 --matrix BLOSUM62 --gapopen 11 --gapextend 1 --more-sensitive \
    -f 6 qseqid pident length qframe qstart qend sstart send evalue bitscore :::: names

#run prediction
conda activate renv

parallel -j $PROC python3 $CREMA/bin/predict.py -f {}.fa -c {}.cpat -d {}.diamond :::: names
