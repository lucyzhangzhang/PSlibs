#!/bin/sh

# blastn -query -o -outfmt 6 qseqid title sblastnames -value 0.001 -num_threads 10 -db /scratch/blastdb/nt

# usage ribosomeFilter.sh -f multifata -b blast.out

echo "Running from $HOSTNAME";

help_msg="Removes all the ribosomal sequences based on blast output,\n
\n
Usage:\n
    $(basename "$0") [-h] -f <multifasta.fa> [-o output.fasta] [-b blast out]\n
    \n
    where:\n
    -h	shows help message.\n
    -f	input fasta to be filtered.\n
    -o	output fasta filename (default: <name>-ribofiltered.fa).\n
    -b	-outfmt 6 of blastn.\n
    -t	Number of threads to use (default: 1)."

# defining options
while getopts ":h:f:o:b:t:" option; do
    case "$option" in
        h)
            echo -e $help_msg >&2;
            exit;
            ;;
        f)
            INFILE=$OPTARG >&2;
            echo Input fasta: $INFILE;
            if [[ ! -e $INFILE ]]; then
                echo "Input file doesn't exist!";
                exit;
            fi
            filename=$(basename -- "$INFILE");
            filename="${filename%.*}";
            
            OUTFILE="$filename-ribofiltered.fa";
            ;;
        o)
            OUTFILE=${OPTARG:-$OUTFILE}; >&2;
            ;;
        b)
            BLAST=${OPTARG:-blastout} >&2;
            ;;
        t)
            PTHREAD=${OPTARG:-1} >&2;
            ;;
        \?)
            echo "ERROR: Invalid option. -$OPTARG" >&2;
            echo -e $help_msg;
            exit 1;
            ;;
        :)
            echo "ERROR: Option -$OPTARG requires an argument." >&2;
            echo -e $help_msg;
            exit 1;
            ;;
        *)
            echo "ERROR: No such option: -$OPTARG" >&2;
            exit 1;
            ;;
    esac
done


if ((OPTIND == 1))
then
    echo "ERROR: No options specified";
    exit 1;
fi


echo "Output fasta: $OUTFILE";

# blastn -query $INFILE -o blastout -outfmt "6 stitle sblastnames scomnames" -evalue 0.001 -num_threads 5 -db /scratch/blastdbt

InName=( $(awk '{print $1}' $BLAST | uniq) );
#awk '{print $1}' $BLAST | uniq | head -n 100  > tmpname

let "counter=0"

echo "Blast file: $BLAST";

#grepper() {
#    if [ $1 == DROUGHT* ]; then
#        counter=$((counter+1))
#    fi
#    if [ ! -z `grep -qi "$1.*(ribosom|protein).*" $BLAST` ]
#    then
#    else
#        echo $1 >> tmpFilterName;
#    fi
#}

for name in "${InName[@]}"; do
    if grep -Eqi "$name.*(protein|ribosom).*" $BLAST; then
        ((counter++))
    else
        echo $name >> tmpFilterName;
    fi
done

#export -f grepper

#parallel --bar -j $PTHREAD -k grepper ::: "${InName[@]}";
#parallel --bar -j $PTHREAD -k grepper :::: tmpname

xargs samtools faidx /2/scratch/lucy/misc/2transcripts.fa < tmpFilterName > ribosomPostFilter.fa;

if [ -e tmpFilterName ]; then
    rm tmpFilterName;
fi

echo "Complete!";
echo "$counter genes were filtered!";
echo "(╯°□°)╯︵ ┻━┻";
