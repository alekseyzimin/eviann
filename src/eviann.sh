#!/bin/bash
#this pipeline generates genome annotation using hisat2, Stringtie2 and maker
GENOME="genome.fa"
PROTEINFILE="proteins.fa"
GENOMEFILE="genome.fa"
RNASEQ_PAIRED="paired"
RNASEQ_UNPAIRED="unpaired"
ALT_EST="altest"
BATCH_SIZE=1000000
MAX_INTRON=100000
UNIPROT="uniprot.fa"
MYPATH="`dirname \"$0\"`"
MYPATH="`( cd \"$MYPATH\" && pwd )`"
SNAP=0
PID=$$
export PATH=$MYPATH:$PATH;
set -o pipefail
NUM_THREADS=1
GC=
RC=
NC=
if tty -s < /dev/fd/1 2> /dev/null; then
    GC='\e[0;32m'
    RC='\e[0;31m'
    NC='\e[0m'
fi

trap abort 1 2 15
function abort {
log "Aborted"
rm -rf /dev/shm/tmp_$PID
kill -9 0
exit 1
}

function usage {
 echo "Usage: eugene.sh [options]"
 echo "Options:"
 echo "-t <number of threads, default:1>"
 echo "-g <MANDATORY:genome fasta file with full path>"
 echo "-p <file containing list of filenames paired Illumina reads from RNAseq experiment, one pair of /path/filename per line; if files are fasta, add "fasta" as the third field on the line>"
 echo "-u <file containing list of filenames of unpaired Illumina reads from RNAseq experiment, one /path/filename per line; if files are fasta, add "fasta" as the third field on the line>"
 echo "-e <fasta file with transcripts from related species>"
 echo "-r <MANDATORY:fasta file of protein sequences to be used with the transcripts for annotation>"
 echo "-s <MANDATORY:fasta file of uniprot proteins>"
 echo "-m <max intron size, default: 100000>"
 echo "-d <add de novo gene finding pass with SNAP>"
 echo "-v <verbose flag>"
 echo "One or more of the -p -u or -e must be supplied."
 echo "De novo gene finding pass will find more exons at the expense of many false positives."
}

log () {
    dddd=$(date)
    echo -e "${GC}[$dddd]${NC} $@"
}


function error_exit {
    dddd=$(date)
    echo -e "${RC}[$dddd]${NC} $1" >&2
    exit "${2:-1}"
}
if [ $# -lt 1 ];then
  usage
  error_exit ""
fi

#parsing arguments
while [[ $# > 0 ]]
do
    key="$1"

    case $key in
        -t|--threads)
            NUM_THREADS="$2"
            shift
            ;;
        -g|--genome)
            GENOMEFILE="$2"
            shift
            ;;
        -p|--paired)
            RNASEQ_PAIRED="$2"
            shift
            ;;
        -e|--est)
            ALT_EST="$2"
            shift
            ;;
        -r|--proteins)
            PROTEINFILE="$2"
            shift
            ;;
        -s|--swissprot)
            UNIPROT="$2"
            shift
            ;;
        -u|--unpaired)
            RNASEQ_UNPAIRED="$2"
            shift
            ;;
        -m|--max-intron)
            MAX_INTRON="$2"
            shift
            ;;
        -v|--verbose)
            set -x
            ;;
        -d|--denovo)
            SNAP=1
            ;;
        -h|--help|-u|--usage)
            usage
            exit 0
            ;;
        *)
            echo "Unknown option $1"
            exit 1        # unknown option
            ;;
    esac
    shift
done

GENOME=`basename $GENOMEFILE`
PROTEIN=`basename $PROTEINFILE`
#checking is dependencies are installed
for prog in $(echo "ufasta hisat2 stringtie2 gffread blastp tblastn makeblastdb");do
  which $prog
  if [ $? -gt 0 ];then error_exit "$prog not found the the PATH";fi
done

#checking inputs
mkdir -p tttt && cd tttt
if [ ! -s $RNASEQ_PAIRED ] && [ ! -s $RNASEQ_UNPAIRED ]  && [ ! -s $ALT_EST ];then
  cd .. && error_exit "Must specify at least one non-empty file with filenames of RNAseq reads with -p or -u or a file with ESTs from the same or closely related species with -e.  Paths for ALL files must be ABSOLUTE."
fi
if [ ! -s $UNIPROT ];then
  cd  .. && error_exit "File with uniprot sequences is missing or specified improperly, please supply it with -s </path_to/uniprot_file.fa> with an ABSOLUTE Path"
fi
if [ ! -s $PROTEINFILE ];then
  cd .. && error_exit "File with proteing sequences for related species is missing or specified improperly, please supply it with -r </path_to/proteins_file.fa> with an ABSOLUTE Path"
fi
cd .. && rm -rf tttt

#path to genome does not have to be absolute
if [ ! -s $GENOMEFILE ];then
  cd .. && error_exit "File with genome sequence is missing or specified improperly, please supply it with -g </path_to/genome_file.fa>"
fi

#first we align
if [ ! -e align-build.success ];then
  log "building HISAT2 index"
  hisat2-build $GENOMEFILE $GENOME.hst 1>hisat2-build.out 2>&1 && touch align-build.success && rm -f align.success
fi

if [ ! -e align.success ];then
  log "aligning RNAseq reads"
  echo "#!/bin/bash" >hisat2.sh
  if [ -s $ALT_EST ];then
    echo "if [ ! -s tissue0.bam ];then hisat2 $GENOME.hst -f --dta -p $NUM_THREADS -U $ALT_EST 2>tissue0.err | samtools view -bhS /dev/stdin > tissue0.bam.tmp && mv tissue0.bam.tmp tissue0.bam;fi" >> hisat2.sh
  fi
  if [ -s $RNASEQ_PAIRED ];then
    awk 'BEGIN{n=1}{
    if($NF == "fasta"){
      print "if [ ! -s tissue"n".bam ];then hisat2 '$GENOME'.hst -f --dta -p '$NUM_THREADS' -1 "$1" -2 "$2" 2>tissue"n".err | samtools view -bhS /dev/stdin > tissue"n".bam.tmp && mv tissue"n".bam.tmp tissue"n".bam; fi"; n++
    }else{
      print "if [ ! -s tissue"n".bam ];then hisat2 '$GENOME'.hst --dta -p '$NUM_THREADS' -1 "$1" -2 "$2" 2>tissue"n".err | samtools view -bhS /dev/stdin > tissue"n".bam.tmp && mv tissue"n".bam.tmp tissue"n".bam; fi"; n++
    }
    }' $RNASEQ_PAIRED >> hisat2.sh
  fi
  if [ -s $RNASEQ_UNPAIRED ];then
    START=1;
    if [ -s $RNASEQ_PAIRED ];then
      START=`wc -l $RNASEQ_PAIRED | awk '{print int($1)+1}'`
    fi
    awk 'BEGIN{n=int("'$START'");}{
    if($NF == "fasta"){
      print "if [ ! -s tissue"n".bam ];then hisat2 '$GENOME'.hst -f --dta -p '$NUM_THREADS' -U "$1" 2>tissue"n".err | samtools view -bhS /dev/stdin > tissue"n".bam.tmp && mv tissue"n".bam.tmp tissue"n".bam; fi"; n++
    }else{
      print "if [ ! -s tissue"n".bam ];then hisat2 '$GENOME'.hst --dta -p '$NUM_THREADS' -U "$1" 2>tissue"n".err | samtools view -bhS /dev/stdin > tissue"n".bam.tmp && mv tissue"n".bam.tmp tissue"n".bam; fi"; n++
    }}' $RNASEQ_UNPAIRED >> hisat2.sh
  fi
  bash ./hisat2.sh && touch align.success && rm -f sort.success || error_exit "Alignment with HISAT2 failed, please check if reads files exist"
fi

NUM_TISSUES=`ls tissue*.bam| grep -v sorted |wc -l`
if [ -e tissue0.bam ];then
  let NUM_TISSUES=$NUM_TISSUES-1;
fi

if [ ! -e sort.success ];then
  log "sorting alignment files"
  if [ $NUM_TISSUES -gt 0 ];then
    for f in $(seq 1 $NUM_TISSUES);do
      if [ ! -s tissue$f.bam.sorted.bam ] && [ -s tissue$f.bam ];then
        samtools sort -@ $NUM_THREADS -m 1G tissue$f.bam tissue$f.bam.sorted.tmp && mv tissue$f.bam.sorted.tmp.bam tissue$f.bam.sorted.bam
      fi
    done
  fi
  touch sort.success && rm -f stringtie.success
fi

if [ ! -e stringtie.success ] && [ -e sort.success ];then
  if [ $NUM_TISSUES -gt 0 ];then
    log "assembling transcripts with Stringtie"
    for f in $(seq 1 $NUM_TISSUES);do
      if [ ! -s tissue$f.bam.sorted.bam.gtf ];then
        stringtie2 -p $NUM_THREADS tissue$f.bam.sorted.bam -o tissue$f.bam.sorted.bam.gtf.tmp && mv tissue$f.bam.sorted.bam.gtf.tmp tissue$f.bam.sorted.bam.gtf
      fi
    done
    OUTCOUNT=`ls tissue*.bam.sorted.bam.gtf|wc -l`
    if [ $OUTCOUNT -eq $NUM_TISSUES ];then
      log "merging transcripts"
      stringtie2 --merge tissue*.bam.sorted.bam.gtf  -o $GENOME.gtf.tmp && mv $GENOME.gtf.tmp $GENOME.gtf
    else
      error_exit "one or more Stringtie jobs failed"
    fi
  fi
  touch stringtie.success && rm -f split.success
fi

if [ ! -e protein_align.success ];then
  protein2genome.sh -t $NUM_THREADS -a $GENOMEFILE -p $PROTEINFILE -m $MAX_INTRON
  if [ -s $GENOME.$PROTEIN.palign.gff ];then
    touch protein_align.success
  fi
fi

if [ ! -e merge.success ];then 
  log "Merging protein and transcript alignments"
  gffcompare -D -o $GENOME.palign.uniq $GENOME.$PROTEIN.palign.gff && \
  gffcompare -o $GENOME.protref -r $GENOME.palign.uniq.combined.gtf $GENOME.gtf && \
  gffread -F $GENOME.palign.uniq.combined.gtf | \
  combine_gene_protein_gff.pl <(gffread -F $GENOME.protref.annotated.gtf ) > $GENOME.gff.tmp && \
  mv $GENOME.gff.tmp $GENOME.gff && \
  gffread -g $GENOMEFILE -w $GENOME.transcripts.fasta -y $GENOME.proteins.fasta $GENOME.gff && \
  touch merge.success && rm -f functional.success pseudo_detect.success 
fi

if [ ! -e functional.success ];then
  log "performing functional annotation" && \
  makeblastdb -in $UNIPROT -input_type fasta -dbtype prot -out uniprot && \
  blastp -db uniprot -query $GENOME.proteins.fasta -out  $GENOME.maker2uni.blastp -evalue 0.000001 -outfmt 6 -num_alignments 1 -seg yes -soft_masking true -lcase_masking -max_hsps 1 -num_threads $NUM_THREADS && \
  my_maker_functional_gff $UNIPROT $GENOME.maker2uni.blastp $GENOME.gff > $GENOME.functional_note.gff.tmp && mv $GENOME.functional_note.gff.tmp $GENOME.functional_note.gff && \
  my_maker_functional_fasta $UNIPROT $GENOME.maker2uni.blastp $GENOME.proteins.fasta > $GENOME.functional_note.proteins.fasta.tmp  && mv $GENOME.functional_note.proteins.fasta.tmp $GENOME.functional_note.proteins.fasta && \
  my_maker_functional_fasta $UNIPROT $GENOME.maker2uni.blastp $GENOME.transcripts.fasta > $GENOME.functional_note.transcripts.fasta.tmp  && mv $GENOME.functional_note.transcripts.fasta.tmp $GENOME.functional_note.transcripts.fasta && \
  touch functional.success && rm -rf pseudo_detect.success
fi

if [ -e functional.success ] && [ ! -e pseudo_detect.success ];then
  log "detecting and annotating pseudogenes"
  ufasta extract -v -f <(awk '{if($3=="gene" || $3=="exon") print $0" "$3}' $GENOME.gff |uniq -c -f 9  | awk '{if($1==1 && $4=="exon"){split($10,a,":");split(a[1],b,"="); print b[2]}}' ) $GENOME.proteins.fasta > $GENOME.proteins.mex.fasta.tmp && mv $GENOME.proteins.mex.fasta.tmp $GENOME.proteins.mex.fasta && \
  ufasta extract -f <(awk '{if($3=="gene" || $3=="exon") print $0" "$3}' $GENOME.gff |uniq -c -f 9  | awk '{if($1==1 && $4=="exon"){split($10,a,":");split(a[1],b,"="); print b[2]}}' ) $GENOME.proteins.fasta > $GENOME.proteins.sex.fasta.tmp && mv $GENOME.proteins.sex.fasta.tmp $GENOME.proteins.sex.fasta && \
  makeblastdb -dbtype prot  -input_type fasta -in  $GENOME.proteins.mex.fasta -out $GENOME.proteins.mex && \
  blastp -db $GENOME.proteins.mex -query $GENOME.proteins.sex.fasta -out  $GENOME.sex2mex.blastp -evalue 0.000001 -outfmt "6 qseqid qlen length pident bitscore" -num_alignments 1 -seg yes -soft_masking true -lcase_masking -max_hsps 1 -num_threads $NUM_THREADS && \
  perl -ane '{if($F[3]>90 && $F[2]/($F[1]+1)>0.90){$pseudo{$F[0]}=1;}}END{open(FILE,"'$GENOME'.functional_note.gff");while($line=<FILE>){chomp($line);@f=split(/\s+/,$line);print $line; ($id,$junk)=split(/;/,$f[8]);if($f[2] eq "gene" && defined($pseudo{substr($id,3)."-mRNA-1"})){ print "pseudo=true;\n";}else{print "\n"}}}' $GENOME.sex2mex.blastp > $GENOME.functional_note.pseudo_label.gff.tmp && mv $GENOME.functional_note.pseudo_label.gff.tmp $GENOME.functional_note.pseudo_label.gff && touch pseudo_detect.success
fi

if [ -e functional.success ] && [ -e pseudo_detect.success ];then
  log "Output annotation is in $GENOME.functional_note.pseudo_label.gff $GENOME.functional_note.proteins.fasta $GENOME.functional_note.transcripts.fasta"
fi


