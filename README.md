# Transcriptomic Analysis of Pi-S treatment libraries

Here is the documentation of my methods for processing the Pi-S treatment libraries, my scripts and locations of the data

## Read locations

```
# Raw reads
cd /2/scratch/lucy/PS/reads

# Trimmed reads
cd /2/scratch/lucy/PS/trimmed
```

## Attempt to find novel transcripts
```
cd /2/scratch/lucy/PS/StringTIe
```
There were many novel transcripts almost all of which overlap with already annotated transcripts and differ by a could of base-pairs. In fact, some of the annotated `DROUGHT` transcripts also overlap with the `Thhalv`s, eventually one should look through the annotation file and remove the overlapping `DROUGHT` transcripts.
