#!/bin/bash
GC=
RC=
NC=
if tty -s < /dev/fd/1 2> /dev/null; then
    GC='\e[0;32m'
    RC='\e[0;31m'
    NC='\e[0m'
fi

function usage {
 echo "Usage: get_annotation_stats.sh [options]"
 echo "Options:"
 echo " -a FILE MANDATORY:annotation file in GFF format default: none"
 echo " -g FILE MANDATORY:genome fasta file default: none"
}

function error_exit {
    dddd=$(date)
    echo -e "${RC}[$dddd]${NC} $1" >&2
    exit "${2:-1}"
}

#parsing arguments
while [[ $# > 0 ]]
do
    key="$1"

    case $key in
        -a|--annotation)
            GFF="$2"
            if [ ! -s $GFF ];then error_exit "genome file $GFF is empty or does not exist!";fi
            shift
            ;;
        -g|--genome)
            GENOME="$2"
            if [ ! -s $GENOME ];then error_exit "genome file $GENOME is empty or does not exist!";fi
            shift
            ;;
    esac
    shift
done

echo -n "Number of genes: ";awk -F'\t' '{if($3=="gene")print $0}' $GFF |wc -l
echo -n "Number of protein coding genes: ";awk -F'\t' '{if($3=="gene")print $0}' $GFF | grep protein_coding | grep -v  'pseudo=true' | wc -l
echo -n "Number of processed pseudo genes: ";awk -F'\t' '{if($3=="gene")print $0}' $GFF | grep 'pseudo=true' |wc -l
echo -n "Number of protein coding transcripts: ";awk -F'\t' '{if($3=="mRNA")print $0}' $GFF |grep -v 'pseudo=true' |wc -l
echo -n "Number of long non-coding RNAs: "; awk -F '\t' '{if($3=="lnc_RNA") print}'  $GFF |wc -l
echo -n "Number of processed pseudo gene transcripts: ";awk -F'\t' '{if($3=="mRNA")print $0}' $GFF | grep 'pseudo=true' |wc -l
echo -n "Number of distinct proteins: "; gffread -y /dev/stdout -g $GENOME $GFF |ufasta one | grep -v '^>' |sort -S 10% |uniq |wc -l
echo -n "Number of functional protein coding transcripts: ";awk '{if($3=="mRNA")print $0}' $GFF |grep Similar |wc -l
