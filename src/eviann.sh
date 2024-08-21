#!/bin/bash
#this pipeline generates genome annotation using hisat2, Stringtie2 and maker
PROTEINFILE="$PWD/uniprot_sprot.nonred.85.fasta"
GENOMEFILE="na"
RNASEQ="na"
ALT_EST="na"
export BATCH_SIZE=1000000
export MAX_INTRON=250000
export MIN_TPM=0.25
export DEBUG=0
UNIPROT="uniprot_sprot.nonred.85.fasta"
MYPATH="`dirname \"$0\"`"
MYPATH="`( cd \"$MYPATH\" && pwd )`"
PID=$$
export PATH=$MYPATH:$MYPATH/SNAP:$PATH;
set -o pipefail
NUM_THREADS=1
FUNCTIONAL=0
USE_SNAP=1
LIFTOVER=0
GC=
RC=
NC=
if tty -s < /dev/fd/1 2> /dev/null; then
    GC='\e[0;32m'
    RC='\e[0;31m'
    NC='\e[0m'
fi

trap abort 1 2 3 9 15
function abort {
log "Aborted"
kill -9 0
exit 1
}

function usage {
 echo "Usage: eviann.sh [options]"
 echo "Options:"
 echo " -t INT     number of threads, default: 1"
 echo " -g FILE    MANDATORY:genome fasta file default: none"
 echo " -r FILE    file containing list of filenames of reads from transcriptome sequencing experiments, default: none
 
  FORMAT OF THIS FILE:
  Each line in the file must refer to sequencing data from a single experiment.
  Please combine runs so that one file/pair/triplet of files contains a single sample.  
  The lines are in the following format:
 
 /path/filename /path/filename /path/filename tag
  or
 /path/filename /path/filename tag
  or
 /path/filename tag

  Fields are space-separated, no leading space. \"tag\" indicates type of data referred to in the preceding fields.  Possible values are:
 
  fastq -- indicates the data is Illumina RNA-seq in fastq format, expects one or a pair of /path/filename.fastq before the tag
  fasta -- indicates the data is Illumina RNA-seq in fasta format, expects one or a pair of /path/filename.fasta before the tag
  bam -- indicates the data is aligned Illumina RNA-seq reads, expects one /path/filename.bam before the tag
  bam_isoseq -- indicates the data is aligned PacBio Iso-seq reads, expects one /path/filename.bam before the tag
  isoseq -- indicates the data is PacBio Iso-seq reads in fasta or fastq format, expects one /path/filename.(fasta or fastq) before the tag
  mix -- indicates the data is from the sample sequenced with both Illumina RNA-seq provided in fastq format and long reads (Iso-seq or Oxford Nanopore) in fasta/fastq format, expects three /path/filename before the tag
  bam_mix -- indicates the data is from the same sample sequenced with both Illumina RNA-seq provided in bam format and long reads (Iso-seq or Oxford Nanopore) in bam format, expects two /path/filename.bam before the tag
 
  Absense of a tag assumes fastq tag and expects one or a pair of /path/filename.fastq on the line.
 "
 echo " -e FILE    fasta file with assembled transcripts from related species, default: none"
 echo " -p FILE    fasta file with protein sequences from (preferrably multiple) related species, uniprot proteins are used of this file is not provided, default: none"
 echo " -m INT     max intron size, default: 250000"
 echo " -l         liftover mode, optimizes internal parameters for annotation liftover; also useful when supplying proteins from a single species, default: not set"
 echo " -f         perform functional annotation, default: not set"
 echo " --debug    keep intermediate output files, default: not set"
 echo " --verbose  verbose run, default: not set"
 echo " --version  report version and exit, default: not set"
 echo ""
 echo " -r or -e MUST be supplied."
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
        -r|--rnaseq)
            RNASEQ="$2"
            shift
            ;;
        -e|--est)
            ALT_EST="$2"
            shift
            ;;
        -p|--proteins)
            PROTEINFILE="$2"
            shift
            ;;
        -s|--swissprot)
            UNIPROT="$2"
            shift
            ;;
        -l|--liftover)
            LIFTOVER=1
            log "Liftover mode ON"
            ;;
        -f|--functional)
            FUNCTIONAL=1
            log "Will perform functional annotation"
            ;;
        -m|--max-intron)
            MAX_INTRON="$2"
            shift
            ;;
        --verbose)
            set -x
            ;;
        --version)
            echo "version 2.0.0"
            exit 0
            ;;
        --debug)
            DEBUG=1
            log "DEBUG mode on, will not clean up intermediate files"
            ;;
        -h|--help|-u|--usage)
            usage
            exit 255
            ;;
        *)
            echo "Unknown option $1"
            usage
            exit 1        # unknown option
            ;;
    esac
    shift
done

#get absolute paths
GENOMEFILE=`realpath $GENOMEFILE`
PROTEINFILE=`realpath $PROTEINFILE`
RNASEQ=`realpath $RNASEQ`
RNASEQ_UNPAIRED=`realpath RNASEQ_UNPAIRED`
ALT_EST=`realpath $ALT_EST`
UNIPROT=`realpath $UNIPROT`

#get filenames to use as prefixes
GENOME=`basename $GENOMEFILE`
PROTEIN=`basename $PROTEINFILE`

#checking is dependencies are installed
log "Checking dependencies"
for prog in $(echo "ufasta stringtie gffread gffcompare blastp tblastn makeblastdb miniprot snap TransDecoder.Predict TransDecoder.LongOrfs");do
  echo -n "Checking for $prog in $MYPATH ... " && \
  which $prog || error_exit "$prog not found in $MYPATH, please make sure installation of EviAnn ran correctly!";
done
for prog in $(echo "minimap2 hisat2 hisat2-build");do
  echo -n "Checking for $prog on the PATH... " && \
  which $prog || log "WARNING! $prog not found the the PATH, it may or may not be needed, but we ask that it is installed!";
done
echo "Checking if TransDecoder is properly installed and works"
$MYPATH/TransDecoder.Predict --version 1>/dev/null 2>&1  || error_exit "TransDecoder seems to be missing some Perl dependencies. Please run $MYPATH/TransDecoder.Predict to see what is missing."
$MYPATH/TransDecoder.LongOrfs --version 1>/dev/null 2>&1  || error_exit "TransDecoder seems to be missing some Perl dependencies. Please run $MYPATH/TransDecoder.LongOrfs to see what is missing."
log "All dependencies checks passed"

#unpack uniprot
if [ ! -s $UNIPROT ];then
  log "Unpacking Uniprot database" && \
  gunzip -c $MYPATH/uniprot_sprot.nonred.85.fasta.gz > uniprot_sprot.nonred.85.fasta.tmp && \
  mv uniprot_sprot.nonred.85.fasta.tmp uniprot_sprot.nonred.85.fasta 
fi

#checking inputs
if [ ! -s $RNASEQ ] && [ ! -s $ALT_EST ];then
  error_exit "Must specify at least one non-empty file with RNA sequencing data with -r or a file with ESTs from the same or closely related species with -e"
fi
if [ ! -s $UNIPROT ];then
  error_exit "File with uniprot sequences is missing or specified improperly, please supply it with -s </path_to/uniprot_file.fa>"
fi
if [ ! -s $PROTEINFILE ];then
  echo "WARNING: proteins from related species are not specified, or file $PROTEINFILE is missing. Using Uniprot proteins as fallback option" && \
  export PROTEINFILE=`realpath $UNIPROT` && \
  export PROTEIN=`basename $UNIPROT`
fi
if [ ! -s $GENOMEFILE ];then
  error_exit "File with genome sequence is missing or specified improperly, please supply it with -g </path_to/genome_file.fa>"
fi

NUM_PROTEINS=`grep '>' $PROTEINFILE |wc -l`
if [ $NUM_PROTEINS -lt 1 ];then
  error_exit "Invalid format for the proteins file $PROTEINFILE, must be in fasta format"
fi

if [ ! -e transcripts_assemble.success ];then
  if [ -s $ALT_EST ];then
    log "Processing transcripts from related species"
    if [ ! -s tissue0.bam.sorted.bam.gtf ];then 
      minimap2 -a -u f -x splice:hq -t $NUM_THREADS -G $MAX_INTRON $GENOMEFILE $ALT_EST 2>tissue0.err | \
      $MYPATH/samtools view -bhS /dev/stdin > tissue0.bam.tmp && \
      mv tissue0.bam.tmp tissue0.bam  
      samtools view -h tissue0.bam | \
      tee >(grep ^@ > tissue0.header) | \
      grep -v ^@ | \
      sort -S 50% -nrk5,5 |\
      perl -ane '{if($h{$F[0]}<2 || $F[4] >0){$h{$F[0]}++;print}}' > tissue0.filter && \
      cat tissue0.header tissue0.filter | samtools view -bhS /dev/stdin | \
      samtools sort -@ $NUM_THREADS -m 1G /dev/stdin tissue0.bam.sorted.tmp && mv tissue0.bam.sorted.tmp.bam tissue0.bam.sorted.bak && \
      samtools view -h tissue0.bam.sorted.bak | fix_splice_junctions.pl $ALT_EST | \
      perl -e '{while($line=<STDIN>){if($line =~ /^@/){print $line}else{chomp($line);@f=split(/\t/,$line);for(my $i=1;$i<=6;$i++){print $f[0],".$i\t",join("\t",@f[1..$#f]),"\n";}}}}' |samtools view -bhS /dev/stdin > tissue0.bam.sorted.tmp.bam && \
      mv tissue0.bam.sorted.tmp.bam tissue0.bam.sorted.bam && \
      rm -f tissue0.bam.sorted.bak && \
      stringtie -g 1 -m 100 -t -j 1 -l REFSTRG -f 0.01 -c 1 -p $NUM_THREADS tissue0.bam.sorted.bam -o tissue0.bam.sorted.bam.gtf.tmp && \
      mv tissue0.bam.sorted.bam.gtf.tmp tissue0.bam.sorted.bam.gtf
    fi
  fi
  echo "#!/bin/bash" > hisat_stringtie.sh
  if [ -s $RNASEQ ];then
    CHECK_RNASEQ=`grep -P "^>|^@" $RNASEQ | wc -l`
    if [ $CHECK_RNASEQ -gt 0 ];then
      error_exit "Wrong format of the RNA sequencing description file $RNASEQ, please check your inputs and option switches!"
    fi
    log "Aligning and building transcripts from RNAseq reads"
    awk 'BEGIN{n=1}{
      if(NF == 4){
        if($NF == "mix"){
          print "if [ ! -s tissue"n".bam.sorted.bam.gtf ];then\nif [ ! -s "$1" ];then echo \"Input data file "$1" does not exist or empty!\";exit 1;fi\nFIRSTCHAR=`zcat -f "$1" | head -n 1 | cut -b 1`\nif [ $FIRSTCHAR = \">\" ];then\n hisat2 '$GENOME'.hst --dta -p '$NUM_THREADS' --min-intronlen 20 --max-intronlen '$MAX_INTRON' -1 "$1" -2 "$2" 2>tissue"n".err | '$MYPATH'/samtools view -bhS /dev/stdin > tissue"n".bam.tmp && mv tissue"n".bam.tmp tissue"n".bam && minimap2 -a -u f -x splice -t '$NUM_THREADS' -G '$MAX_INTRON' '$GENOMEFILE' "$3" 2>tissue"n".err | '$MYPATH'/samtools view -bhS /dev/stdin > tissue"n"_lr.bam.tmp && mv tissue"n"_lr.bam.tmp tissue"n"_lr.bam && samtools sort -@ '$NUM_THREADS' -m 1G tissue"n".bam tissue"n".bam.sorted.tmp && mv tissue"n".bam.sorted.tmp.bam tissue"n".bam.sorted.bam && samtools sort -@ '$NUM_THREADS' -m 1G tissue"n"_lr.bam tissue"n"_lr.bam.sorted.tmp && mv tissue"n"_lr.bam.sorted.tmp.bam tissue"n"_lr.bam.sorted.bam && rm tissue"n".bam tissue"n"_lr.bam && ./run_stringtie_mix.sh tissue"n".bam.sorted.bam tissue"n"_lr.bam.sorted.bam && rm tissue"n"_lr.bam.sorted.bam || exit 1\nelse echo \"WARNING! Invalid fasta format files "$1" or "$2", ignoring them\"\nfi\nfi";
          n++;
        }
      }else if(NF == 3){
        if($NF == "fasta"){
          print "if [ ! -s tissue"n".bam.sorted.bam.gtf ];then\nif [ ! -s "$1" ];then echo \"Input data file "$1" does not exist or empty!\";exit 1;fi\nFIRSTCHAR=`zcat -f "$1" | head -n 1 | cut -b 1`\nif [ $FIRSTCHAR = \">\" ];then\n hisat2 '$GENOME'.hst -f --dta -p '$NUM_THREADS' --min-intronlen 20 --max-intronlen '$MAX_INTRON' -1 "$1" -2 "$2" 2>tissue"n".err | '$MYPATH'/samtools view -bhS /dev/stdin > tissue"n".bam.tmp && mv tissue"n".bam.tmp tissue"n".bam && samtools sort -@ '$NUM_THREADS' -m 1G tissue"n".bam tissue"n".bam.sorted.tmp && mv tissue"n".bam.sorted.tmp.bam tissue"n".bam.sorted.bam && rm tissue"n".bam && ./run_stringtie.sh tissue"n".bam.sorted.bam || exit 1\nelse echo \"WARNING! Invalid fasta format files "$1" or "$2", ignoring them\"\nfi\nfi"; 
          n++;
        }else if($NF == "fastq"){
          print "if [ ! -s tissue"n".bam.sorted.bam.gtf ];then\nif [ ! -s "$1" ];then echo \"Input data file "$1" does not exist or empty!\";exit 1;fi\nFIRSTCHAR=`zcat -f "$1" | head -n 1 | cut -b 1`\nif [ $FIRSTCHAR = \"@\" ];then\n hisat2 '$GENOME'.hst --dta -p '$NUM_THREADS' --min-intronlen 20 --max-intronlen '$MAX_INTRON' -1 "$1" -2 "$2" 2>tissue"n".err | '$MYPATH'/samtools view -bhS /dev/stdin > tissue"n".bam.tmp && mv tissue"n".bam.tmp tissue"n".bam && samtools sort -@ '$NUM_THREADS' -m 1G tissue"n".bam tissue"n".bam.sorted.tmp && mv tissue"n".bam.sorted.tmp.bam tissue"n".bam.sorted.bam && rm tissue"n".bam && ./run_stringtie.sh tissue"n".bam.sorted.bam || exit 1\nelse echo \"WARNING! Invalid fastq format files "$1" or "$2", ignoring them\"\nfi\nfi";
          n++;
        }else if($NF == "bam_mix"){
          print "if [ ! -s tissue"n".bam.sorted.bam.gtf ];then\nif [ ! -s "$1" ];then echo \"Input data file "$1" does not exist or empty!\";exit 1;fi\ncp "$1" tissue"n".bam.tmp && mv tissue"n".bam.tmp tissue"n".bam && cp "$2" tissue"n"_lr.bam.tmp && mv tissue"n"_lr.bam.tmp tissue"n"_lr.bam && samtools sort -@ '$NUM_THREADS' -m 1G tissue"n".bam tissue"n".bam.sorted.tmp && mv tissue"n".bam.sorted.tmp.bam tissue"n".bam.sorted.bam && samtools sort -@ '$NUM_THREADS' -m 1G tissue"n"_lr.bam tissue"n"_lr.bam.sorted.tmp && mv tissue"n"_lr.bam.sorted.tmp.bam tissue"n"_lr.bam.sorted.bam && rm tissue"n".bam tissue"n"_lr.bam && ./run_stringtie_mix.sh tissue"n".bam.sorted.bam tissue"n"_lr.bam.sorted.bam && rm tissue"n"_lr.bam.sorted.bam || exit 1\nfi";
          n++;
        }
      }else if(NF == 2){
        if($NF == "bam"){
          print "if [ ! -s tissue"n".bam.sorted.bam.gtf ];then\nif [ ! -s "$1" ];then echo \"Input data file "$1" does not exist or empty!\";exit 1;fi\n cp "$1" tissue"n".bam.tmp && mv tissue"n".bam.tmp tissue"n".bam && samtools sort -@ '$NUM_THREADS' -m 1G tissue"n".bam tissue"n".bam.sorted.tmp && mv tissue"n".bam.sorted.tmp.bam tissue"n".bam.sorted.bam && rm tissue"n".bam && ./run_stringtie.sh tissue"n".bam.sorted.bam || exit 1\nfi"; 
          n++;
        }else if($NF == "bam_isoseq"){
          print "if [ ! -s tissue"n".bam.sorted.bam.gtf ];then\nif [ ! -s "$1" ];then echo \"Input data file "$1" does not exist or empty!\";exit 1;fi\n cp "$1" tissue"n".bam.tmp && mv tissue"n".bam.tmp tissue"n".bam && samtools sort -@ '$NUM_THREADS' -m 1G tissue"n".bam tissue"n".bam.sorted.tmp && mv tissue"n".bam.sorted.tmp.bam tissue"n".bam.sorted.bam && rm tissue"n".bam && ./run_stringtie_lr.sh tissue"n".bam.sorted.bam || exit 1\nfi"; 
          n++;
        }else if($NF == "isoseq"){
          print "if [ ! -s tissue"n".bam.sorted.bam.gtf ];then\nif [ ! -s "$1" ];then echo \"Input data file "$1" does not exist or empty!\";exit 1;fi\nFIRSTCHAR=`zcat -f "$1" | head -n 1 | cut -b 1`\nif [ $FIRSTCHAR = \">\" ] || [ $FIRSTCHAR = \"@\" ];then\n minimap2 -a -u f -x splice:hq -t '$NUM_THREADS' -G '$MAX_INTRON' '$GENOMEFILE' "$1" 2>tissue"n".err | '$MYPATH'/samtools view -bhS /dev/stdin > tissue"n".bam.tmp && mv tissue"n".bam.tmp tissue"n".bam && samtools sort -@ '$NUM_THREADS' -m 1G tissue"n".bam tissue"n".bam.sorted.tmp && mv tissue"n".bam.sorted.tmp.bam tissue"n".bam.sorted.bam && rm tissue"n".bam && ./run_stringtie_lr.sh tissue"n".bam.sorted.bam || exit 1\nelse echo \"WARNING! Invalid fasta format file "$1", ignoring it\"\nfi\nfi";
          n++;
        }else if($NF == "fasta"){
          print "if [ ! -s tissue"n".bam.sorted.bam.gtf ];then\nif [ ! -s "$1" ];then echo \"Input data file "$1" does not exist or empty!\";exit 1;fi\nFIRSTCHAR=`zcat -f "$1" | head -n 1 | cut -b 1`\nif [ $FIRSTCHAR = \">\" ];then\n hisat2 '$GENOME'.hst -f --dta -p '$NUM_THREADS' --min-intronlen 20 --max-intronlen '$MAX_INTRON' -U "$1" 2>tissue"n".err | '$MYPATH'/samtools view -bhS /dev/stdin > tissue"n".bam.tmp && mv tissue"n".bam.tmp tissue"n".bam && samtools sort -@ '$NUM_THREADS' -m 1G tissue"n".bam tissue"n".bam.sorted.tmp && mv tissue"n".bam.sorted.tmp.bam tissue"n".bam.sorted.bam && rm tissue"n".bam && ./run_stringtie.sh tissue"n".bam.sorted.bam || exit 1\nelse echo \"WARNING! Invalid fasta format file "$1", ignoring it\"\nfi\nfi";
          n++;
        }else if($NF == "fastq"){
          print "if [ ! -s tissue"n".bam.sorted.bam.gtf ];then\nif [ ! -s "$1" ];then echo \"Input data file "$1" does not exist or empty!\";exit 1;fi\nFIRSTCHAR=`zcat -f "$1" | head -n 1 | cut -b 1`\nif [ $FIRSTCHAR = \"@\" ];then\n hisat2 '$GENOME'.hst --dta -p '$NUM_THREADS' --min-intronlen 20 --max-intronlen '$MAX_INTRON' -U "$1" 2>tissue"n".err | '$MYPATH'/samtools view -bhS /dev/stdin > tissue"n".bam.tmp && mv tissue"n".bam.tmp tissue"n".bam && samtools sort -@ '$NUM_THREADS' -m 1G tissue"n".bam tissue"n".bam.sorted.tmp && mv tissue"n".bam.sorted.tmp.bam tissue"n".bam.sorted.bam && rm tissue"n".bam && ./run_stringtie.sh tissue"n".bam.sorted.bam || exit 1\nelse echo \"WARNING! Invalid fastq format file "$1", ignoring it\"\nfi\nfi";
          n++;
        }else{
          print "if [ ! -s tissue"n".bam.sorted.bam.gtf ];then\nif [ ! -s "$1" ];then echo \"Input data file "$1" does not exist or empty!\";exit 1;fi\nFIRSTCHAR=`zcat -f "$1" | head -n 1 | cut -b 1`\nif [ $FIRSTCHAR = \"@\" ];then\n hisat2 '$GENOME'.hst --dta -p '$NUM_THREADS' --min-intronlen 20 --max-intronlen '$MAX_INTRON' -1 "$1" -2 "$2" 2>tissue"n".err | '$MYPATH'/samtools view -bhS /dev/stdin > tissue"n".bam.tmp && mv tissue"n".bam.tmp tissue"n".bam && samtools sort -@ '$NUM_THREADS' -m 1G tissue"n".bam tissue"n".bam.sorted.tmp && mv tissue"n".bam.sorted.tmp.bam tissue"n".bam.sorted.bam && rm tissue"n".bam && ./run_stringtie.sh tissue"n".bam.sorted.bam || exit 1\nelse echo \"WARNING! Invalid fastq format files "$1" or "$2", ignoring them\"\nfi\nfi"
          n++;
        }
      }else if(NF == 1){
        print "if [ ! -s tissue"n".bam.sorted.bam.gtf ];then\nif [ ! -s "$1" ];then echo \"Input data file "$1" does not exist or empty!\";exit 1;fi\nFIRSTCHAR=`zcat -f "$1" | head -n 1 | cut -b 1`\nif [ $FIRSTCHAR = \"@\" ];then\n hisat2 '$GENOME'.hst -f --dta -p '$NUM_THREADS' --min-intronlen 20 --max-intronlen '$MAX_INTRON' -U "$1" 2>tissue"n".err | '$MYPATH'/samtools view -bhS /dev/stdin > tissue"n".bam.tmp && mv tissue"n".bam.tmp tissue"n".bam && samtools sort -@ '$NUM_THREADS' -m 1G tissue"n".bam tissue"n".bam.sorted.tmp && mv tissue"n".bam.sorted.tmp.bam tissue"n".bam.sorted.bam && rm tissue"n".bam && ./run_stringtie.sh tissue"n".bam.sorted.bam || exit 1\nelse echo \"WARNING! Invalid fastq format file "$1", ignoring it\"\nfi\nfi";
        n++;
      }
    }' $RNASEQ >> hisat_stringtie.sh
  fi
  echo "#!/bin/bash
  $MYPATH/stringtie -p $NUM_THREADS \$1 -o \$1.gtf.tmp && \\
    awk -F '\t' 'BEGIN{flag=0}{if(\$3==\"transcript\"){n=split(\$9,a,\";\");for(i=1;i<=n;i++){if(a[i] ~ /TPM/){ m=split(a[i],b,\"\\\"\");tpm=b[m-1];}else if(a[i] ~ /FPKM/){ m=split(a[i],b,\"\\\"\");fpkm=b[m-1];}}if(fpkm > $MIN_TPM || tpm > $MIN_TPM ) flag=1; else flag=0;}if(flag){print \$0}}' \$1.gtf.tmp > \$1.gtf.filtered.tmp && \\
    mv \$1.gtf.filtered.tmp \$1.gtf && \\
    rm -f \$1.gtf.tmp " > run_stringtie.sh && \
    chmod 0755 run_stringtie.sh && \
  echo "#!/bin/bash
  $MYPATH/stringtie -p $NUM_THREADS \$1 -L -o \$1.gtf.tmp && \\
    awk -F '\t' 'BEGIN{flag=0}{if(\$3==\"transcript\"){n=split(\$9,a,\";\");for(i=1;i<=n;i++){if(a[i] ~ /TPM/){ m=split(a[i],b,\"\\\"\");tpm=b[m-1];}else if(a[i] ~ /FPKM/){ m=split(a[i],b,\"\\\"\");fpkm=b[m-1];}}if(fpkm > $MIN_TPM || tpm > $MIN_TPM ) flag=1; else flag=0;}if(flag){print \$0}}' \$1.gtf.tmp > \$1.gtf.filtered.tmp && \\
    mv \$1.gtf.filtered.tmp \$1.gtf  && \\
    rm -f \$1.gtf.tmp " > run_stringtie_lr.sh && \
    chmod 0755 run_stringtie_lr.sh && \
  echo "#!/bin/bash
  $MYPATH/stringtie -p $NUM_THREADS \$1 \$2 --mix -o \$1.gtf.tmp && \\
    awk -F '\t' 'BEGIN{flag=0}{if(\$3==\"transcript\"){n=split(\$9,a,\";\");for(i=1;i<=n;i++){if(a[i] ~ /TPM/){ m=split(a[i],b,\"\\\"\");tpm=b[m-1];}else if(a[i] ~ /FPKM/){ m=split(a[i],b,\"\\\"\");fpkm=b[m-1];}}if(fpkm > $MIN_TPM || tpm > $MIN_TPM ) flag=1; else flag=0;}if(flag){print \$0}}' \$1.gtf.tmp > \$1.gtf.filtered.tmp && \\
    mv \$1.gtf.filtered.tmp \$1.gtf  && \\
    rm -f \$1.gtf.tmp " > run_stringtie_mix.sh && \
  chmod 0755 run_stringtie.sh && \
  chmod 0755 run_stringtie_lr.sh && \
  chmod 0755 run_stringtie_mix.sh && \
  HISAT=`grep hisat2 hisat_stringtie.sh| wc -l` && \
  if [ $HISAT -gt 0 ] && [ ! -e align-build.success ];then
    log "Building HISAT2 index" && \
    hisat2-build $GENOMEFILE $GENOME.hst 1>/dev/null 2>&1 && \
    touch align-build.success || error_exit "Building HISAT2 index failed, check your inputs"
  fi
  bash ./hisat_stringtie.sh && \
  touch transcripts_assemble.success && \
  rm -f transcripts_merge.success || error_exit "Alignment with HISAT2 or transcript assembly with StringTie failed, please check if reads files exist and formatted correctly"
fi

NUM_TISSUES=`ls tissue*.bam.sorted.bam|wc -l`
if [ -e tissue0.bam.sorted.bam ];then
  let NUM_TISSUES=$NUM_TISSUES-1;
fi

if [ -e transcripts_assemble.success ] && [ ! -e  transcripts_merge.success ];then
  OUTCOUNT=`ls tissue*.bam.sorted.bam.gtf|wc -l`
  if [ $OUTCOUNT -eq 1 ];then
    # a single tissue, we add 1 for the number of samples and the TPMs
    cat tissue*.bam.sorted.bam.gtf | grep -v '^#' | perl -F'\t' -ane '{
      if($F[2] eq "transcript"){
        if($F[8]=~/transcript_id \"(\S+)\";/){
          $transcript=$1;
          $tpm=1;
          $tpm=$1 if($F[8] =~ /TPM \"(\S+)\";/);
          $new_transcript_id="M$transcript:1:$tpm";
        }
      }
      $F[8]=~s/transcript_id \"$transcript\";/transcript_id \"$new_transcript_id\";/;
      print join("\t",@F);
    }' > $GENOME.gtf.tmp && mv $GENOME.gtf.tmp $GENOME.gtf &&  touch transcripts_merge.success && rm -f merge.success || error_exit "Failed to merge transcripts"
  elif [ $OUTCOUNT -ge $NUM_TISSUES ];then
    log "Merging transcripts"
    gffcompare -STC  tissue*.bam.sorted.bam.gtf  -o $GENOME.tmp -p MSTRG 1>gffcompare.out 2>&1 && \
    awk '{tpm=0;num_samples=0;for(i=4;i<=NF;i++){if($i ~ /^q/){num_samples++;split($i,a,"|");if(a[5]>tpm){tpm=a[5]}}}print $1" "tpm" "num_samples}'  $GENOME.tmp.tracking > $GENOME.max_tpm.samples.tmp && \
    mv $GENOME.max_tpm.samples.tmp $GENOME.max_tpm.samples.txt && \
    perl -F'\t' -ane '
    BEGIN{
      $thresh=int('$NUM_TISSUES')/50;
      open(FILE,"'$GENOME'.max_tpm.samples.txt");
      while($line=<FILE>){
        chomp($line);
        @f=split(/\s+/,$line);
        $tpm{$f[0]}=$f[1];
        $samples{$f[0]}=$f[2];
      }
    }{
      if($F[2] eq "transcript"){
        $flag=1;
        if($F[8]=~/transcript_id \"(\S+)\";/){
          $transcript=$1;
          $new_transcript_id="$transcript:$samples{$transcript}:$tpm{$transcript}";
          $flag=0 if($samples{$transcript}<$thresh);
        }else{
          die("no transcript id found on line ".join("\t",@F));
        }
      }
      $F[8]=~s/transcript_id \"$transcript\";/transcript_id \"$new_transcript_id\";/;
      print join("\t",@F) if($flag);
    }' $GENOME.tmp.combined.gtf  > $GENOME.tmp2.combined.gtf && \
    mv $GENOME.tmp2.combined.gtf $GENOME.gtf && \
    rm -f $GENOME.tmp.{combined.gtf,tracking,loci} $GENOME.tmp && touch transcripts_merge.success && rm -f merge.success || error_exit "Failed to merge transcripts"
  else
    error_exit "one or more Stringtie jobs failed to run properly"
  fi
fi

if [ ! -e protein2genome.deduplicate.success ];then
  log "Deduplicating input proteins"
  $MYPATH/ufasta one $PROTEINFILE | \
  awk '{if($0 ~ /^>/){header=$1}else{print header,$1}}' |\
  sort  -S 10% -k2,2 |\
  uniq -f 1 |\
  awk '{print $1"\n"$2}' | \
  tr ':' '_' > $PROTEIN.uniq.tmp && \
  mv $PROTEIN.uniq.tmp $PROTEIN.uniq && \
  touch protein2genome.deduplicate.success && \
  rm -f protein2genome.align.success || error_exit "Failed in deduplicating proteins"
fi
PROTEIN=$PROTEIN.uniq

if [ ! -e protein2genome.align.success ];then
  log "Aligning proteins to the genome with miniprot"
  miniprot -p 0.95 -N 20 -k 5 -t $NUM_THREADS -G $MAX_INTRON --gff $GENOMEFILE $PROTEIN 2>miniprot.err | \
    $MYPATH/convert_miniprot_gff.pl > $GENOME.$PROTEIN.palign.gff.tmp && \
    mv $GENOME.$PROTEIN.palign.gff.tmp $GENOME.$PROTEIN.palign.gff && \
    rm -f merge.success && \
    touch protein2genome.align.success || error_exit "Alignment of proteins to the genome with miniprot failed, please check miniprot.err"
fi

if [ -e transcripts_merge.success ] && [ -e protein2genome.align.success ] && [ ! -e merge.success ];then
  log "Deriving gene models from protein and transcript alignments"
  rm -f merge.unused.success merge.u.success snap.success && \
#we fix suspect intons in the protein alignment files.  If an intron has never been seen before, switch it to the closest one that has been seen
  gffread -F  <( fix_suspect_introns.pl $GENOME.gtf < $GENOME.$PROTEIN.palign.gff ) > $GENOME.palign.fixed.gff.tmp && \
  mv $GENOME.palign.fixed.gff.tmp $GENOME.palign.fixed.gff && \
#here we use the "fixed" protein alignments as reference and compare our transcripts. This annotates each transcript with a protein match and a match code
  gffcompare -T -o $GENOME.protref -r $GENOME.palign.fixed.gff $GENOME.gtf && \
  rm -f $GENOME.protref.{loci,tracking,stats} $GENOME.protref && \
  cat $GENOME.palign.fixed.gff | filter_by_class_code.pl $GENOME.protref.annotated.gtf | \
  filter_by_local_abundance.pl > $GENOME.transcripts_to_keep.txt.tmp && \
  mv $GENOME.transcripts_to_keep.txt.tmp $GENOME.transcripts_to_keep.txt && \
  mv $GENOME.protref.annotated.gtf $GENOME.protref.annotated.gtf.bak && \
  perl -F'\t' -ane 'BEGIN{open(FILE,"'$GENOME'.transcripts_to_keep.txt");while($line=<FILE>){chomp($line);$h{$line}=1}}{if($F[8]=~/transcript_id \"(\S+)\";/){print if(defined($h{$1}));}}' $GENOME.protref.annotated.gtf.bak > $GENOME.protref.annotated.gtf.tmp && \
  mv $GENOME.protref.annotated.gtf.tmp $GENOME.protref.annotated.gtf && rm -f $GENOME.protref.annotated.gtf.bak && \
  perl -F'\t' -ane 'BEGIN{open(FILE,"'$GENOME'.transcripts_to_keep.txt");while($line=<FILE>){chomp($line);$h{$line}=1}}{if($F[8]=~/transcript_id \"(\S+)\";/){print if(defined($h{$1}));}}' $GENOME.gtf > $GENOME.abundanceFiltered.gtf.tmp && \
  mv $GENOME.abundanceFiltered.gtf.tmp $GENOME.abundanceFiltered.gtf && \
#here we combine the transcripts and protein matches
#unused proteins gff file contains all protein alignments that did not match the transcripts; we will use them later
#this produces files $GENOME.{k,u}.gff.tmp  and $GENOME.unused_proteins.gff.tmp
  cat $GENOME.palign.fixed.gff |  combine_gene_protein_gff.pl $GENOME <( gffread -F $GENOME.protref.annotated.gtf ) $GENOMEFILE 1>combine.out 2>&1 && \
  mv $GENOME.k.gff.tmp $GENOME.k.gff && \
  if [ $DEBUG -lt 1 ];then
    rm -rf $GENOME.protref.annotated.gtf $GENOME.transcripts_to_keep.txt
  fi 
  mv $GENOME.u.gff.tmp $GENOME.u.gff && \
  mv $GENOME.unused_proteins.gff.tmp $GENOME.unused_proteins.gff && \
  if [ -s $GENOME.unused_proteins.gff ];then
    log "Filtering unused protein only loci" && \
#we use preliminary "k" file to train SNAP
    if [ $USE_SNAP -gt 0 ];then
      log "Training SNAP and predicting protein coding genes to use in filtering aligned proteins" && \
      mkdir -p SNAP && \
      (cd SNAP && \
        rm -rf * && \
        awk '{print $1}' $GENOMEFILE > $GENOME.onefield.fa && \
        $MYPATH/SNAP/gff3_to_zff.pl $GENOME.onefield.fa <(gffread -F ../$GENOME.k.gff |grep -P '\-mRNA\-1$|\-mRNA-1;' | perl -F'\t' -ane '{if($F[2] eq "mRNA"){$flag=0;if($F[8] =~ /Class==;/){$flag=1}}print if($flag)}') > $GENOME.zff && \
        $MYPATH/SNAP/fathom -categorize 1000 $GENOME.zff $GENOME.onefield.fa && \
        $MYPATH/SNAP/fathom -export 1000 -plus uni.* && \
        $MYPATH/SNAP/forge export.ann export.dna && \
        $MYPATH/SNAP/hmm-assembler.pl $GENOME . > ../$GENOME.hmm.tmp && \
        mv ../$GENOME.hmm.tmp ../$GENOME.hmm && \
        $MYPATH/ufasta split -i $GENOME.onefield.fa $GENOME.onefield.1.fa $GENOME.onefield.2.fa $GENOME.onefield.3.fa $GENOME.onefield.4.fa $GENOME.onefield.5.fa $GENOME.onefield.6.fa $GENOME.onefield.7.fa $GENOME.8.onefield.fa && \
        echo '#!/bin/bash
        if [ -s '$GENOME'.onefield.$1.fa ];then
        '$MYPATH'/SNAP/snap -quiet -gff ../'$GENOME'.hmm '$GENOME'.onefield.$1.fa > '$GENOME'.onefield.$1.gff
        else
        echo "" > '$GENOME'.onefield.$1.gff
        fi' > run_snap.sh && \
        chmod 0755 run_snap.sh && \
        seq 1 8 | xargs -P 8 -I {} ./run_snap.sh {} && \
        cat $GENOME.onefield.{1,2,3,4,5,6,7,8}.gff | \
        perl -F'\t' -ane 'BEGIN{$id="";}{
            if($#F == 8){
              if(not($id eq $F[8])){
                if(not($id eq "")){
                  print "$chrom\tSNAP\tmRNA\t$start\t$end\t.\t$dir\t.\tID=$id";
                  print join("",@exons);
                }
                $id=$F[8];
                $chrom=$F[0];
                $dir=$F[6];
                $start=$F[3];
                $end=$F[4];
                @exons=();
              }
              $F[2]="exon";
              $F[8]="Parent=$F[8]";
              if($dir eq "+"){
                $end=$F[4];
              }else{
                $start=$F[3];
              }
              push(@exons,join("\t",@F));
              $F[2]="cds";
              push(@exons,join("\t",@F));
            }
          }END{
            if(not($id eq "")){
              print "$chrom\tSNAP\tmRNA\t$start\t$end\t.\t$dir\t.\tID=$id";
              print join("",@exons);
            }
          }' > $GENOME.snap.gff.tmp && \
        mv $GENOME.snap.gff.tmp ../$GENOME.snap.gff) && \
      gffcompare -T -r $GENOME.snap.gff $GENOME.unused_proteins.gff -o $GENOME.snapcompare && \
      grep -P 'class_code "(=|c)"' $GENOME.snapcompare.annotated.gtf |\
      grep -v "contained_in" |\
      perl -F'\t' -ane '{if($F[8]=~/^transcript_id\s\"(\S+)\";\sgene_id/){print $1,"\n"}}' > $GENOME.snap_match.txt.tmp && \
      mv $GENOME.snap_match.txt.tmp $GENOME.snap_match.txt && \
      touch snap.success && \
      if [ $DEBUG -lt 1 ];then
        rm -rf SNAP $GENOME.snapcompare $GENOME.snapcompare.{annotated.gtf,loci,tracking}
      fi
    fi
    gffread -V -y $GENOME.unused.faa -g $GENOMEFILE $GENOME.unused_proteins.gff && \
    ufasta one $GENOME.unused.faa |\
    awk '{if($1 ~ /^>/){name=substr($1,2)}else{split(name,a,":");print name" "$1}}' |sort -k2,2 -S 10% |uniq -c -f 1 |awk '{print $2" "$1}' > $GENOME.protein_count.txt.tmp && \
    mv $GENOME.protein_count.txt.tmp $GENOME.protein_count.txt && \
    gffread --cluster-only <(awk '{if($3=="cds" || $3=="transcript") print $0}' $GENOME.unused_proteins.gff) | \
    filter_unused_proteins.pl $GENOMEFILE $GENOME.unused_proteins.gff $GENOME.snap_match.txt $GENOME.protein_count.txt $LIFTOVER > $GENOME.best_unused_proteins.gff.tmp && \
    mv $GENOME.best_unused_proteins.gff.tmp $GENOME.best_unused_proteins.gff && touch merge.unused.success
    if [ $DEBUG -lt 1 ];then
      rm -rf $GENOME.protein_count.txt $GENOME.unused.faa $GENOME.unused_proteins.gff
    fi
  else
    echo "" > $GENOME.best_unused_proteins.gff && touch merge.unused.success
  fi 
  #these are u's -- no match to a protein
  if [ -s $GENOME.u.gff ];then
    log "Looking for ORFs in transcripts with no protein matches"
    gffread -g $GENOMEFILE -w $GENOME.lncRNA.fa $GENOME.u.gff && \
    rm -rf $GENOME.lncRNA.fa.transdecoder* && \
    TransDecoder.LongOrfs -S -t $GENOME.lncRNA.fa 1>transdecoder.LongOrfs.out 2>&1 && \
    makeblastdb -in $UNIPROT -input_type fasta -dbtype prot -out uniprot 1>makeblastdb.out 2>&1 && \
    blastp -query $GENOME.lncRNA.fa.transdecoder_dir/longest_orfs.pep -db uniprot  -max_target_seqs 1 -outfmt 6  -evalue 0.000001 -num_threads $NUM_THREADS 2>blastp1.out > $GENOME.lncRNA.blastp.tmp && \
    mv $GENOME.lncRNA.blastp.tmp $GENOME.lncRNA.u.blastp && \
    TransDecoder.Predict -t $GENOME.lncRNA.fa --single_best_only --retain_blastp_hits $GENOME.lncRNA.u.blastp 1>transdecoder.Predict.out 2>&1
    if [ -s $GENOME.lncRNA.fa.transdecoder.gff3 ];then
      add_cds_to_gff.pl <(awk -F '\t' 'BEGIN{flag=0}{if($3=="gene"){if($9~/ORF type:complete/){flag=1}else{flag=0}}if(flag){print}}' $GENOME.lncRNA.fa.transdecoder.gff3) <  $GENOME.u.gff | \
      perl -F'\t' -ane 'BEGIN{
          open(FILE,"'$GENOME'.lncRNA.u.blastp");
          while($line=<FILE>){
            $line=~/^(.+)\.p\d+/;
            $f=$1;
            $f=~s/_lncRNA//;
            $h{$f}=1;
          }
        }{
          unless($F[2] eq "exon" || $F[2] eq "gene" || $F[2] =~ /_UTR/ || $F[8] =~/_lncRNA/){
            @f=split(/;|:/,$F[8]); 
            if(defined($h{substr($f[0],3)})){
              if($F[2] eq "mRNA"){
                print join("\t",@F[0..1]),"\tgene\t",join("\t",@F[3..8]);
              }else{
                print join("\t",@F);
              }
            }
          }
        }' > $GENOME.u.cds.gff.tmp && \
      mv $GENOME.u.cds.gff.tmp $GENOME.u.cds.gff && touch merge.u.success
    else
      echo "#gff" > $GENOME.u.cds.gff && touch merge.u.success
    fi
    rm -rf $GENOME.lncRNA.fa $GENOME.lncRNA.u.blastp pipeliner.*.cmds $GENOME.lncRNA.fa.transdecoder_dir  $GENOME.lncRNA.fa.transdecoder_dir.__checkpoints $GENOME.lncRNA.fa.transdecoder_dir.__checkpoints_longorfs transdecoder.LongOrfs.out $GENOME.lncRNA.fa.transdecoder.{cds,pep} blastp1.out transdecoder.Predict.out 
  else
    echo "#gff" > $GENOME.u.cds.gff 
  fi
  if [ $DEBUG -lt 1 ];then
    rm -rf $GENOME.u.gff $GENOME.lncRNA.fa.transdecoder.{gff3,bed}
  fi
  log "Working on final merge"
  gffcompare -T -r $GENOME.snap.gff $GENOME.abundanceFiltered.gtf -o $GENOME.protref.snap && \
  if [ -s $GENOME.best_unused_proteins.gff ];then
#now we have additional proteins produced by transdecoder, let's use them all
    gffcompare -STC $GENOME.best_unused_proteins.gff $GENOME.abundanceFiltered.gtf -o $GENOME.all && \
    gffread $GENOME.palign.fixed.gff $GENOME.u.cds.gff \
      <(cat $GENOME.protref.snap.annotated.gtf |\
      perl -F'\t' -ane '{
        if($F[8]=~/cmp_ref "(\S+)"; class_code "(=|k)"/){
          $h{$1}=1;
        }
      }END{
        $flag=0;
        open(FILE,"'$GENOME'.snap.gff");
        while($line=<FILE>){
          chomp($line);
          @F=split(/\t/,$line);
          if($F[2] eq "mRNA"){
            $flag=defined($h{substr($F[8],3)}) ? 1:0;
            $F[2]="gene";
          }
          print join("\t",@F),"\n" if($flag);
        }
      }') >  $GENOME.palign.all.gff && \
    gffcompare -T -o $GENOME.protref.all -r $GENOME.palign.all.gff $GENOME.all.combined.gtf && \
    log "Checking for and repairing broken ORFs" && \
    cat $GENOME.palign.all.gff | filter_by_class_code.pl $GENOME.protref.all.annotated.gtf | gffread -F > $GENOME.protref.all.annotated.class.gff.tmp && \
    mv $GENOME.protref.all.annotated.class.gff.tmp $GENOME.protref.all.annotated.class.gff && \
    cat $GENOME.palign.all.gff |  check_cds.pl $GENOME $GENOME.protref.all.annotated.class.gff $GENOMEFILE 1>check_cds.out 2>&1 && \
    mv $GENOME.good_cds.fa.tmp $GENOME.good_cds.fa && \
    mv $GENOME.broken_cds.fa.tmp $GENOME.broken_cds.fa && \
    mv $GENOME.broken_ref.txt.tmp $GENOME.broken_ref.txt && \
    cat $PROTEINFILE <(gffread -y /dev/stdout -g $GENOMEFILE $GENOME.u.cds.gff) | ufasta extract -f $GENOME.broken_ref.txt > $GENOME.broken_ref.faa.tmp &&\
    mv $GENOME.broken_ref.faa.tmp $GENOME.broken_ref.faa && \
    rm -rf $GENOME.broken_cds.fa.transdecoder* && \
    TransDecoder.LongOrfs -S -t $GENOME.broken_cds.fa 1>transdecoder.LongOrfs.out 2>&1 && \
    makeblastdb -in $GENOME.broken_ref.faa -input_type fasta -dbtype prot -out broken_ref 1>makeblastdb.out 2>&1 && \
    blastp -query $GENOME.broken_cds.fa.transdecoder_dir/longest_orfs.pep -db broken_ref  -max_target_seqs 1 -outfmt 6  -evalue 0.000001 -num_threads $NUM_THREADS 2>blastp2.out > $GENOME.broken_cds.blastp.tmp && \
    mv $GENOME.broken_cds.blastp.tmp $GENOME.broken_cds.blastp && \
    TransDecoder.Predict -t $GENOME.broken_cds.fa --single_best_only --retain_blastp_hits $GENOME.broken_cds.blastp 1>transdecoder.Predict.out 2>&1
    if [ -s $GENOME.broken_cds.fa.transdecoder.bed ];then
      awk -F '\t' '{if(NF>8 && !($4 ~/ORF_type:internal/)) print $1" "$7" "$8}' $GENOME.broken_cds.fa.transdecoder.bed  > $GENOME.fixed_cds.txt.tmp && \
      mv $GENOME.fixed_cds.txt.tmp $GENOME.fixed_cds.txt
    fi
    rm -rf transdecoder.Predict.out $GENOME.broken_cds.fa pipeliner.*.cmds $GENOME.broken_cds.fa.transdecoder_dir  $GENOME.broken_cds.transdecoder_dir.__checkpoints $GENOME.broken_cds.fa.transdecoder_dir.__checkpoints_longorfs transdecoder.LongOrfs.out $GENOME.broken_cds.fa.transdecoder.{cds,pep,gff3} && \
    cat $GENOME.palign.all.gff |  combine_gene_protein_gff.pl $GENOME $GENOME.protref.all.annotated.class.gff $GENOMEFILE $GENOME.fixed_cds.txt $GENOME.hmm <(perl -F'\t' -ane '{if($F[2] eq "transcript"){print "$1 $3\n" if($F[8] =~ /transcript_id "(.+)"; gene_id "(.+)"; oId "(.+)"; tss_id "(.+)"; num_samples "(.+)";$/);}}'  $GENOME.all.combined.gtf) 1>combine.out 2>&1 && \
    gffread -F --keep-exon-attrs --keep-genes $GENOME.k.gff.tmp $GENOME.u.gff.tmp | awk '{if($0 ~ /^# gffread/){print "# EviAnn automated annotation"}else{print $0}}' > $GENOME.gff.tmp && \
    mv $GENOME.gff.tmp $GENOME.gff  && \
    rm -f $GENOME.{k,u,unused_proteins}.gff.tmp && \
    touch merge.success && rm -f pseudo_detect.success functional.success || error_exit "Merging transcript and protein evidence failed."
  else
    gffread $GENOME.palign.fixed.gff $GENOME.u.cds.gff \
      <(cat $GENOME.protref.snap.annotated.gtf |\
      perl -F'\t' -ane '{
        if($F[8]=~/cmp_ref "(\S+)"; class_code "(=|k)"/){
          $h{$1}=1;
        }
      }END{
        $flag=0;
        open(FILE,"'$GENOME'.snap.gff");
        while($line=<FILE>){
          chomp($line);
          @F=split(/\t/,$line);
          if($F[2] eq "mRNA"){
            $flag=defined($h{substr($F[8],3)}) ? 1:0;
            $F[2]="gene";
          }
          print join("\t",@F),"\n" if($flag);
        }
      }') >  $GENOME.palign.all.gff && \
    gffcompare -T -o $GENOME.protref.all -r $GENOME.palign.all.gff $GENOME.abundanceFiltered.gtf && \
    log "Checking for and repairing broken ORFs" && \
    cat $GENOME.palign.all.gff | filter_by_class_code.pl $GENOME.protref.all.annotated.gtf | gffread -F > $GENOME.protref.all.annotated.class.gff.tmp && \
    mv $GENOME.protref.all.annotated.class.gff.tmp $GENOME.protref.all.annotated.class.gff && \
    cat $GENOME.palign.all.gff |  check_cds.pl $GENOME $GENOME.protref.all.annotated.class.gff $GENOMEFILE 1>check_cds.out 2>&1 && \
    mv $GENOME.good_cds.fa.tmp $GENOME.good_cds.fa && \
    mv $GENOME.broken_cds.fa.tmp $GENOME.broken_cds.fa && \
    mv $GENOME.broken_ref.txt.tmp $GENOME.broken_ref.txt && \
    cat $PROTEINFILE <(gffread -y /dev/stdout -g $GENOMEFILE $GENOME.u.cds.gff) |ufasta extract -f $GENOME.broken_ref.txt > $GENOME.broken_ref.faa.tmp &&\
    mv $GENOME.broken_ref.faa.tmp $GENOME.broken_ref.faa && \
    rm -rf $GENOME.broken_cds.fa.transdecoder* && \
    TransDecoder.LongOrfs -t $GENOME.broken_cds.fa 1>transdecoder.LongOrfs.out 2>&1 && \
    makeblastdb -in $GENOME.broken_ref.faa -input_type fasta -dbtype prot -out broken_ref 1>makeblastdb.out 2>&1 && \
    blastp -query $GENOME.broken_cds.fa.transdecoder_dir/longest_orfs.pep -db broken_ref  -max_target_seqs 1 -outfmt 6  -evalue 0.000001 -num_threads $NUM_THREADS 2>blastp3.out > $GENOME.broken_cds.blastp.tmp && \
    mv $GENOME.broken_cds.blastp.tmp $GENOME.broken_cds.blastp && \
    TransDecoder.Predict -t $GENOME.broken_cds.fa --single_best_only --retain_blastp_hits $GENOME.broken_cds.blastp 1>transdecoder.Predict.out 2>&1
    if [ -s $GENOME.broken_cds.fa.transdecoder.bed ];then
      awk -F '\t' '{if(NF>8 && !($4 ~/ORF_type:internal/)) print $1" "$7" "$8}' $GENOME.broken_cds.fa.transdecoder.bed  > $GENOME.fixed_cds.txt.tmp && \
      mv $GENOME.fixed_cds.txt.tmp $GENOME.fixed_cds.txt
    fi
    rm -rf transdecoder.Predict.out $GENOME.broken_cds.fa pipeliner.*.cmds $GENOME.broken_cds.fa.transdecoder_dir  $GENOME.broken_cds.transdecoder_dir.__checkpoints $GENOME.broken_cds.fa.transdecoder_dir.__checkpoints_longorfs transdecoder.LongOrfs.out $GENOME.broken_cds.fa.transdecoder.{cds,pep,gff3} && \
    cat $GENOME.palign.all.gff |  combine_gene_protein_gff.pl $GENOME $GENOME.protref.all.annotated.class.gff $GENOMEFILE $GENOME.fixed_cds.txt $GENOME.hmm <(perl -F'\t' -ane '{if($F[2] eq "transcript"){print "$1 $3\n" if($F[8] =~ /transcript_id "(.+)"; gene_id "(.+)"; oId "(.+)"; tss_id "(.+)"; num_samples "(.+)";$/);}}'  $GENOME.all.combined.gtf) 1>combine.out 2>&1 && \
    gffread -F --keep-exon-attrs --keep-genes $GENOME.k.gff.tmp $GENOME.u.gff.tmp | \
    awk '{if($0 ~ /^# gffread/){print "# EviAnn automated annotation"}else{print $0}}' > $GENOME.gff.tmp && \
    mv $GENOME.gff.tmp $GENOME.gff && \
    rm -f $GENOME.{k,u,unused_proteins}.gff.tmp  && \
    touch merge.success && rm -f pseudo_detect.success functional.success || error_exit "Merging transcript and protein evidence failed."
  fi
  rm -rf broken_ref.{pjs,ptf,pto,pot,pdb,psq,phr,pin} makeblastdb.out blastp2.out && \
  if [ $DEBUG -lt 1 ];then
    rm -f $GENOME.all.{loci,stats,tracking} 
    rm -f $GENOME.protref.all.{loci,stats,tracking}
    rm -f $GENOME.protref.snap.{loci,stats,tracking}
    rm -rf $GENOME.all.{combined,redundant}.gtf $GENOME.all $GENOME.palign.all.gff $GENOME.protref.all.annotated.class.gff $GENOME.protref.all.annotated.gtf $GENOME.protref.all $GENOME.good_cds.fa $GENOME.broken_cds.fa $GENOME.broken_ref.{txt,faa} $GENOME.broken_cds.{blastp,fa.transdecoder.bed} $GENOME.fixed_cds.txt 
  fi
  gffread -S -g $GENOMEFILE -w $GENOME.transcripts.fasta $GENOME.gff && \
  if [ ! -e merge.unused.success ];then
    log "Merging proteins without transcript match failed, results will be incomplete!"
  fi
  if [ ! -e merge.u.success ];then
    log "Merging transcripts without protein match failed, results will be incomplete!"
  fi
  if [ ! -e snap.success ];then
    log "Ab initio gene finding with snap failed, results may be suboptimal!"
  fi
  rm -f merge.{unused,u}.success snap.success
fi

if [ -e merge.success ] && [ ! -e pseudo_detect.success ];then
  log "Detecting and annotating processed pseudogenes" && \
  gffread -S -g $GENOMEFILE -y $GENOME.proteins.fasta $GENOME.gff && \
  ufasta extract -f <(awk -F'\t' '{if($3=="exon"){print substr($9,8);}}'  $GENOME.gff|uniq -d) $GENOME.proteins.fasta > $GENOME.proteins.mex.fasta.tmp && \
  mv $GENOME.proteins.mex.fasta.tmp $GENOME.proteins.mex.fasta && \
  ufasta extract -f <(awk -F'\t' '{if($3=="exon"){print substr($9,8);}}'  $GENOME.gff|uniq -c | perl -ane '{($gene,$junk)=split(/-/,$F[1]);$max_count{$gene}=$F[0] if($max_count{$gene}<$F[0]);$transcripts{$gene}.="$F[1]\n";}END{foreach $k(keys %max_count){if($max_count{$k}==1){print $transcripts{$k}}}}') $GENOME.proteins.fasta > $GENOME.proteins.sex.fasta.tmp && \
  mv $GENOME.proteins.sex.fasta.tmp $GENOME.proteins.sex.fasta && \
  if [ -s $GENOME.proteins.mex.fasta ] && [ -s $GENOME.proteins.sex.fasta ];then
    makeblastdb -dbtype prot  -input_type fasta -in  $GENOME.proteins.mex.fasta -out $GENOME.proteins.mex 1>makeblastdb.sex2mex.out 2>&1 && \
    blastp -db $GENOME.proteins.mex -query $GENOME.proteins.sex.fasta -out  $GENOME.sex2mex.blastp.tmp -evalue 0.000001 -outfmt "6 qseqid qlen length pident bitscore sseqid" -num_alignments 5 -seg yes -soft_masking true -lcase_masking -max_hsps 1 -num_threads $NUM_THREADS 1>blastp5.out 2>&1 && \
    mv $GENOME.sex2mex.blastp.tmp $GENOME.sex2mex.blastp && \
    perl -ane '{
      @f1=split(/-/,$F[0]);
      @f2=split(/-/,$F[5]);
      if($F[3]>90 && $F[2]/($F[1]+1)>0.9 && not($f1[0] eq $f2[0])){
        $pseudo{$f1[0]}=1;
      }
    }END{
      open(FILE,"'$GENOME'.gff");
      while($line=<FILE>){
        chomp($line);
        @f=split(/\t/,$line);
        ($id,$junk)=split(/;/,$f[8]);
        ($loc_id,$junk)=split(/-/,$id);
        if($f[2] eq "mRNA" || $f[2] eq "gene"){
          if(defined($pseudo{substr($loc_id,3)})){
            print "$line;pseudo=true\n";
            $printcds=0;
          }else{
            print "$line\n";
            $printcds=1;
          }
        }elsif($f[2] eq "CDS"){
          print "$line\n" if($printcds);
        }else{
          print "$line\n";
        } 
      }
    }' $GENOME.sex2mex.blastp > $GENOME.pseudo_label.gff.tmp && \
    gffread -S -g $GENOMEFILE -y $GENOME.proteins.fasta $GENOME.pseudo_label.gff && \
    mv $GENOME.pseudo_label.gff.tmp $GENOME.pseudo_label.gff && touch pseudo_detect.success || error_exit "Detection of pseudogenes failed, you can use annotation in $GENOME.gff without pseudo-gene labels"
  else
    cp $GENOME.gff $GENOME.pseudo_label.gff.tmp && mv $GENOME.pseudo_label.gff.tmp $GENOME.pseudo_label.gff && touch pseudo_detect.success || error_exit "Detection of pseudogenes failed, you can use annotation in $GENOME.gff without pseudo-gene labels"
  fi
  if [ $DEBUG -lt 1 ];then
    rm -rf $GENOME.proteins.mex.p?? $GENOME.proteins.{s,m}ex.fasta  makeblastdb.sex2mex.out blastp5.out $GENOME.sex2mex.blastp 
  fi
fi

if [ -e merge.success ] && [ -e pseudo_detect.success ];then
  if [ $FUNCTIONAL -gt 0 ];then
    if [ ! -e functional.success ];then
      log "Performing functional annotation" && \
      makeblastdb -in $UNIPROT -input_type fasta -dbtype prot -out uniprot 1>makeblastdb.out 2>&1 && \
      blastp -db uniprot -query $GENOME.proteins.fasta -out  $GENOME.maker2uni.blastp.tmp -evalue 0.000001 -outfmt 6 -num_alignments 1 -seg yes -soft_masking true -lcase_masking -max_hsps 1 -num_threads $NUM_THREADS 1>blastp4.out 2>&1 && \
      mv $GENOME.maker2uni.blastp.tmp $GENOME.maker2uni.blastp && \
      my_maker_functional_gff $UNIPROT $GENOME.maker2uni.blastp $GENOME.pseudo_label.gff > $GENOME.functional_note.pseudo_label.gff.tmp && mv $GENOME.functional_note.pseudo_label.gff.tmp $GENOME.functional_note.pseudo_label.gff && \
      my_maker_functional_fasta $UNIPROT $GENOME.maker2uni.blastp $GENOME.proteins.fasta > $GENOME.functional_note.proteins.fasta.tmp  && mv $GENOME.functional_note.proteins.fasta.tmp $GENOME.functional_note.proteins.fasta && \
      my_maker_functional_fasta $UNIPROT $GENOME.maker2uni.blastp $GENOME.transcripts.fasta > $GENOME.functional_note.transcripts.fasta.tmp  && mv $GENOME.functional_note.transcripts.fasta.tmp $GENOME.functional_note.transcripts.fasta && \
      rm -rf blastp4.out && \
      touch functional.success  || error_exit "Functional annotation failed"
    fi
    log "Output functionally annotated GFF is in $GENOME.functional_note.pseudo_label.gff, proteins are in  $GENOME.functional_note.proteins.fasta, transcripts are in $GENOME.functional_note.transcripts.fasta"
    echo -n "Number of functional protein coding transcripts: ";awk '{if($3=="mRNA")print $0}' $GENOME.functional_note.pseudo_label.gff |grep Similar |wc -l
  else
    log "Output annotation GFF is in $GENOME.pseudo_label.gff, proteins are in  $GENOME.proteins.fasta, transcripts are in $GENOME.transcripts.fasta"
  fi
  echo -n "Number of genes: ";awk -F'\t' '{if($3=="gene")print $0}' $GENOME.pseudo_label.gff |wc -l
  echo -n "Number of processed pseudo gene transcripts: ";awk -F'\t' '{if($3=="mRNA")print $0}' $GENOME.pseudo_label.gff| grep 'pseudo=true' |wc -l
  echo -n "Number of processed pseudo genes: ";awk -F'\t' '{if($3=="gene")print $0}' $GENOME.pseudo_label.gff| grep 'pseudo=true' |wc -l
  echo -n "Number of transcripts: ";awk -F'\t' '{if($3=="mRNA")print $0}' $GENOME.pseudo_label.gff |wc -l
  echo -n "Number of long non-coding RNAs: "; awk -F '\t' '{if($3=="mRNA") print}'  $GENOME.pseudo_label.gff | grep _lncRNA |wc -l
  echo -n "Number of distinct proteins: "; ufasta one $GENOME.proteins.fasta | grep -v '^>' |sort -S 10% |uniq |wc -l
fi

