# Transcriptomic Analysis of Pi-S treatment libraries

Here is the documentation of my methods for processing the Pi-S treatment libraries, my scripts and locations of the data

## Read locations

```
# Raw reads
cd /home/lucy/R/Eutrema/PS/PS/reads

# Trimmed reads
cd /home/lucy/R/Eutrema/PS/PS/trimmed

# Read counts
cd /home/lucy/R/Eutrema/PS/QoRTs
```

## Attempt to find novel transcripts
```
cd /home/lucy/R/Eutrema/PS/StringTIe
```
There were many novel transcripts almost all of which overlap with already annotated transcripts and differ by a could of base-pairs. In fact, some of the annotated `DROUGHT` transcripts also overlap with the `Thhalv`s, eventually one should look through the annotation file and remove the overlapping `DROUGHT` transcripts.

## How does `results()` work in DESeq2

Why don't I get the same DEGs as A vs. B when I switch the reference level to B vs. A?

`results()` command extracts the log2FoldChange comparison values from the DESeqDataSet object generated from `DESeq()`.

## Match between two patterns (`Perl`)
```
perl -pe 
```

## Custom `crema` instructions

Before running `crema`, you have to [filter out ribosomal reads!](./scripts/ribosomeFilter.sh).

(Or don't filter, Elizabeth says that some ribosomal reads can also be lncRNA depending on how they were processed)

(In the spirit of keeping everything, might not even filter FPKMs less than 1 because lncRNA that are related to stress seems to have very low expression in general)

```
#!/bin/bash
```

Output generated from /home/lucy/R/Eutrema/PS/scripts/PS_pairwise.R
```
NAME_DIR=/home/lucy/R/Eutrema/PS/names
CREMA=/home/lucy/bin/crema
PROC=3
cd $NAME_DIR
ls -1 */*/*/names > names
```

Activate Conda environment that has everything already downloaded. Check crema GitHub for software requirements.

```
conda activate renv
```

Use batch processing to generate multi fasta

```
#go on the cluster
#ssh lucy@114
while read line; do
    xargs samtools faidx /2/scratch/lucy/misc/2transcripts.fa < $line > $line.fa
done < names

#exit the cluster
#exit
```

Run CPAT, make sure it's in an environment with Python 2 installed

```
conda deactivate

conda activate cpat

parallel -j $PROC cpat.py -g {}.fa -o {}.cpat \
    -x $CREMA/cpat_models/ath_hexamer.txt \
    -d $CREMA/cpat_models/ath_logit.RData :::: names

conda deactivate
```

Run Diamond BLASTX

```
parallel -j $PROC diamond blastx -d swissprot.dmnd -q {}.fa -o {}.diamond \
    -e 0.001 -k 5 --matrix BLOSUM62 --gapopen 11 --gapextend 1 --more-sensitive \
    -f 6 qseqid pident length qframe qstart qend sstart send evalue bitscore :::: names
```

Run the lncRNA prediction tool

```
conda activate renv
parallel -j $PROC python3 $CREMA/bin/predict.py -f {}.fa -c {}.cpat -d {}.diamond :::: names
```
## Weighted gene co-expression network analysis

With R package WGCNA (no longer supported), however still a good diagnostic analysis to explore data. 

See: [K-means clustering paper](https://bmcsystbiol.biomedcentral.com/articles/10.1186/s12918-017-0420-6#article-info)

The above package is also deprecated and their recommended alternative seems to be optimized for animal tissue types.

[Script](./scripts/PS_wgcna.R)

## Looking for QTLs

See: [*Arabidopsis* QTL mapping paper](https://bmcplantbiol.biomedcentral.com/articles/10.1186/s12870-019-1996-3)

See if these QTLs (in supp file 4 of the paper), particularly the PHO and PUE group are found in *Eutrema* Yukon and if their expression levels remark anything interesting.

* Use Phytozome/Phytomine Python API to get homolog genes and sequences (really powerful tool)
* [R script](./scripts/PS_QTL.R) calls the [Python script](./scripts/phytozomeQTL.py) which queries Phytozome through the Phytomine API
* For genes with no match in *Arabidopsis*, extract the transcript sequence and run BLAST in *Eutrema* transcriptome
* Try to annotate those transcripts

Three common genes mapped in these two QTL groups

```
intersect(PHO, PUE)

  [1] "AT1G24020" "AT1G47970" "AT1G47980"
```
#### QTL groupings

The number of matches may exceed mismatches because one *Arabidopsis* gene may match to multiple *Eutrema* genes.

QTL Group| *Arabidopsis* |*Eutrema* Matched|*Eutrema* BLAST| No match 
---------|--------------:|----------------:|--------------:|---------:
 PHO     | 25            | 18              | 3             | 5        
 PUE     | 43            | 15              | 5             | 4     
 SUL     | 26            | 22              | 3             | 3         
 IP6     | 18            | 17              | 1             | 1        

**PUE** - Phosphate use efficiency
**PHO** - Leaf phosphate content
**SUL** - Leaf sulfate content
**IP6** - Leaf phytate content
