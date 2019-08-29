#!/bin/bash
# HOW TO USE CREMA

CREMA=/home/lucy/bin/crema

# generate protein database to remove proteins from analysis
diamond makedb --in uniprot_sprot.fasta -d swissprot.dmnd

# run diamond blast and output into specific columns for predict.py
# you HAVE to specify which columns to keep in the ouput format
diamond blastx -d swissprot.dmnd -q your_transcript_fasta_file.fa -o diamond_output.txt \
    -e 0.001 -k 5 --matrix BLOSUM62 --gapopen 11 --gapextend 1 --more-sensitive \
    -f 6 qseqid pident length qframe qstart qend sstart send evalue bitscore

# or parallel processing
parallel -j 3 diamond blastx -d swissprot.dmnd -q {} -o {}.diamond \
    -e 0.001 -k 5 --matrix BLOSUM62 --gapopen 11 --gapextend 1 --more-sensitive \
    -f 6 qseqid pident length qframe qstart qend sstart send evalue bitscore :::: names.fa

# extract fasta files of the relevant transcripts from the transcript multifasta
# transcript multifasta is generated with gffread
xargs samtools fadix ~/scratch/misc/2transcripts.fa < ~/R/Eutrema/PS/names/HS.dn.names > ~/R/Eutrema/PS/names/HS.dn.fa

# activate python2 conda environment
# this environment should have CPAT=1.2.1 isntalled
conda activate python2

# cpat.py to calcuate transcript features
cpat.py -g ~/R/Eutrema/PS/names/HS.up.fa -o ~/R/Eutrema/PS/names/HS.up.cpat.out \
    -x $CREMA/cpat_models/ath_hexamer.txt -d $CREMA/cpat_models/ath_logit.RData

conda deactivate

# activate r environment
# this environment should also have biopython and scikit-learn-0.20.0 installed
conda activate renv

# run prediction tool using the full path, calls R
python3 $CREMA/bin/predict.py -f your_transcript_fasta_file.fa -c cpat_output.txt -d diamond_output.txt

conda deactivate

# END


# In the case of batch processing:
ls */*/*/names > names
ls */*/*/names.fa > names.fa

# Batch processing of transcript files 
# xargs samtools faidx /2/scratch/lucy/misc/2transcripts.fa < $1 > $1.fa
parallel -j 3 samtoolExtract.sh :::: names.fa

# or parallel processing
parallel -j 3 diamond blastx -d swissprot.dmnd -q {} -o {}.diamond \
    -e 0.001 -k 5 --matrix BLOSUM62 --gapopen 11 --gapextend 1 --more-sensitive \
    -f 6 qseqid pident length qframe qstart qend sstart send evalue bitscore :::: names.fa

# predict.py outputs everything to the current directory
