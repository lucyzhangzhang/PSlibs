#!/bin/bash

xargs samtools faidx /2/scratch/lucy/misc/2transcripts.fa < $1 > $1.fa
