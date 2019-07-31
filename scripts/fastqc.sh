#!/bin/sh

#Fastqc
fastqc -o ~/scratch/PS/fastqc -f fastq -d ~/scratch/PS/fastqc -t 4 `ls *.gz | tr '\n' ' ' | sed 's/.$//'`
