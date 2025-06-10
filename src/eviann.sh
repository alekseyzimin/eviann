#!/bin/bash
#this pipeline generates genome annotation using hisat2, Stringtie2 and maker
PROTEINFILE="$PWD/uniprot_sprot.fasta"
PLOIDY=2
GENOMEFILE="na"
CDSFILE="na"
RNASEQ="na"
ALT_EST="na"
export BATCH_SIZE=1000000
export MAX_INTRON=1
export MIN_TPM=0.25
export DEBUG=0
export PARTIAL=0
export LNCRNATPM=1.0
EXTRA_GFF="na"
UNIPROT="$PWD/uniprot_sprot.fasta"
MYPATH="`dirname \"$0\"`"
MYPATH="`( cd \"$MYPATH\" && pwd )`"
PID=$$
export PATH=$MYPATH:$PATH;
set -o pipefail
NUM_THREADS=1
FUNCTIONAL=0
LIFTOVER=0
JUNCTION_THRESHOLD=4
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
 echo " -t INT           number of threads, default: 1"
 echo " -g FILE          MANDATORY:genome fasta file default: none"
 echo " -r FILE          file containing list of filenames of reads from transcriptome sequencing experiments, default: none
 
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
 echo " -e FILE               fasta file with assembled transcripts from related species, default: none"
 echo " -p FILE               fasta file with protein sequences from (preferrably multiple) related species, uniprot proteins are used of this file is not provided, default: none"
 echo " -s FILE               fasta file with UniProt-SwissProt proteins to use in functional annotation or if proteins from close relatives are not available. "
 echo "                         EviAnn will download the most recent version of this database automatically."
 echo "                         To use a different version, supply it with this switch. The database is available at:"
 echo "                         https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz"
 echo " -m INT                max intron size, default: auto-determined as sqrt(genome size in kb)*1000; this setting will override automatically estimated value"
 echo " --partial             include transcripts with partial (mising start or stop codon) CDS in the output"
 echo " -d INT                set ploidy for the genome, this value is used in estimating the maximum intron size, default 2"
 echo " -c FILE               GFF file with CDS sequences for THIS genome to be used in annotations. Each CDS must have gene/transcript/mRNA AND exon AND CDS attributes"
 echo " --lncrnamintpm FLOAT  minimum TPM to include non-coding transcript into the annotation as lncRNA, default: 1.0"
 echo " --liftover            liftover mode, optimizes internal parameters for annotation liftover; also useful when supplying proteins from a single species, default: not set"
 echo " -f|--functional       perform functional annotation, default: not set"
 echo " --extra FILE          extra features to add from an external GFF file.  Feautures MUST have gene records.  Any features that overlap with existing annotations will be ignored"
 echo " --debug               keep intermediate output files, default: not set"
 echo " --verbose             verbose run, default: not set"
 echo " --version             report version and exit."
 echo " --help                display this message and exit."
 echo ""
 echo " IMPORTANT!!! -r or -e MUST be supplied."
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
        -c|--cds)
            CDSFILE="$2"
            shift
            ;;
        -d|--ploidy)
            PLOIDY="$2"
            shift
            ;;
        -s|--swissprot)
            UNIPROT="$2"
            shift
            ;;
        -l|--liftover)
            LIFTOVER=1
            JUNCTION_THRESHOLD=-1000
            log "Liftover mode ON"
            ;;
        -f|--functional)
            FUNCTIONAL=1
            log "Will perform functional annotation"
            ;;
        --partial)
            PARTIAL=1
            log "Will include transcripts with partial (mising start or stop codon) CDS in the output"
            ;;
        --lncrnamintpm)
            LNCRNATPM="$2"
            shift
            ;;
        -m|--max-intron)
            MAX_INTRON="$2"
            shift
            ;;
        --extra)
            EXTRA_GFF="$2"
            shift
            ;;
        --verbose)
            set -x
            ;;
        --version)
            echo "version 2.0.3"
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

#checking inputs
if [ ! -s $RNASEQ ] && [ ! -s $ALT_EST ];then
  error_exit "Must specify at least one non-empty file with RNA sequencing data with -r or a file with ESTs from the same or closely related species with -e"
fi

if [ ! -s $PROTEINFILE ] && [ ! -s $CDSFILE ];then
  echo "WARNING: proteins from related species are not specified, or file $PROTEINFILE is missing. Using Uniprot proteins as fallback option" && \
  if [ ! -s $UNIPROT ];then
    log "Downloading UniProt database" && \
    wget --no-check-certificate https://ftp.uniprot.org/pub/databases/uniprot/knowledgebase/complete/uniprot_sprot.fasta.gz && gunzip uniprot_sprot.fasta.gz || error_exit "Downloading UniProt database failed, please obtain it manually"
  fi
  P=`cd "$(dirname "$UNIPROT")" && pwd`
  F=`basename $UNIPROT`
  export PROTEINFILE=$P/$F
  export PROTEIN=$F
fi

if [ ! -s $GENOMEFILE ];then
  error_exit "File with genome sequence is missing or specified improperly, please supply it with -g </path_to/genome_file.fa>"
fi

#get absolute paths
P=`cd "$(dirname "$GENOMEFILE")" && pwd`
F=`basename $GENOMEFILE`
GENOMEFILE=$P/$F
GENOME=$F
if [ -s $PROTEINFILE ];then
  P=`cd "$(dirname "$PROTEINFILE")" && pwd`
  F=`basename $PROTEINFILE`
  PROTEINFILE=$P/$F
  PROTEIN=$F
fi
if [ -s $RNASEQ ];then
  P=`cd "$(dirname "$RNASEQ")" && pwd`
  F=`basename $RNASEQ`
  RNASEQ=$P/$F
fi
if [ -s $ALT_EST ];then
  P=`cd "$(dirname "$ALT_EST")" && pwd`
  F=`basename $ALT_EST`
  ALT_EST=$P/$F
fi
if [ -s $UNIPROT ];then
  P=`cd "$(dirname "$UNIPROT")" && pwd`
  F=`basename $UNIPROT`
  UNIPROT=$P/$F
fi
if [ -s $CDSFILE ];then
  P=`cd "$(dirname "$CDSFILE")" && pwd`
  F=`basename $CDSFILE`
  CDSFILE=$P/$F
  CDS=$F
fi

#checking if dependencies are installed
log "Checking dependencies"
for prog in $(echo "ufasta stringtie gffread gffcompare miniprot TransDecoder.Predict TransDecoder.LongOrfs");do
  echo -n "Checking for $prog in $MYPATH ... " && \
  which $prog || error_exit "$prog not found in $MYPATH, please make sure installation of EviAnn ran correctly!";
done
for prog in $(echo "minimap2 hisat2 hisat2-build samtools makeblastdb blastp");do
  echo -n "Checking for $prog on the PATH ... " && \
  which $prog || error_exit "ERROR! $prog not found the the PATH!";
done
echo "Checking if TransDecoder is properly installed and works"
if ! $MYPATH/TransDecoder.Predict --version 1>/dev/null;then 
  error_exit "TransDecoder seems to be missing some Perl dependencies."
fi
if ! $MYPATH/TransDecoder.Predict --version 1>/dev/null;then
  error_exit "TransDecoder seems to be missing some Perl dependencies."
fi
log "All dependencies checks passed"

if [ $MAX_INTRON -le 1 ];then
  log "Auto-determining the maximum intron size based on the genome size" && \
  MAX_INTRON=`ufasta n50 -S $GENOMEFILE | perl -ane '$p=int('$PLOIDY'); $p=2 if($p<=0);$m=int(sqrt($F[1]/1000/$p*2)*1000); $m=100000 if($m<100000); print $m;'` && \
  echo "Maximum intron size set to $MAX_INTRON"
fi

if [ ! -e transcripts_assemble.success ];then
  if [ -s $ALT_EST ];then
    log "Processing transcripts from related species"
    if [ ! -s tissue0.bam.sorted.bam.gtf ];then 
      cat <(ufasta sizes -H $GENOMEFILE | awk '{print "@SQ\tSN:"$1"\tLN:"$2}') \
        <(minimap2 -a -u f -x splice:hq -t $NUM_THREADS -G $MAX_INTRON $GENOMEFILE $ALT_EST 2>tissue0.err | grep -v '^@') | \
        samtools view -bhS /dev/stdin > tissue0.bam.tmp && \
      mv tissue0.bam.tmp tissue0.bam  
      samtools view -h tissue0.bam | \
      tee >(grep ^@ > tissue0.header) | \
      grep -v ^@ | \
      sort -S 50% -nrk5,5 |\
      perl -ane '{if($h{$F[0]}<2 || $F[4] >0){$h{$F[0]}++;print}}' > tissue0.filter && \
      cat tissue0.header tissue0.filter | samtools view -bhS /dev/stdin | \
      samtools sort -@ $NUM_THREADS -m 1G /dev/stdin -T sts -O bam| \
      samtools view -h /dev/stdin | fix_splice_junctions.pl $ALT_EST | \
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
    log "Parsing the RNA sequencing data file" && \
    awk 'BEGIN{n=1}{
      if(NF == 4){
        if($NF == "mix"){
          print "if [ ! -s tissue"n".bam.sorted.bam.gtf ];then\nif [ ! -s "$1" ];then echo \"Input data file "$1" does not exist or empty!\";exit 1;fi\nFIRSTCHAR=`zcat -f "$1" | head -n 1 | cut -b 1`\nif [ $FIRSTCHAR = \">\" ] || [ $FIRSTCHAR = \"@\" ];then\n hisat2 -x '$GENOME'.hst -p '$NUM_THREADS' --min-intronlen 20 --max-intronlen '$MAX_INTRON' -1 "$1" -2 "$2" 2>tissue"n".err | samtools view -bhS /dev/stdin > tissue"n".bam.tmp && mv tissue"n".bam.tmp tissue"n".bam && samtools sort -@ '$NUM_THREADS' -m 1G tissue"n".bam -T sts -O bam> tissue"n".bam.sorted.tmp.bam && mv tissue"n".bam.sorted.tmp.bam tissue"n".bam.sorted.bam;\n cat <(ufasta sizes -H '$GENOMEFILE' | awk \047{print \"@SQ\\tSN:\"$1\"\\tLN:\"$2}\047) <(minimap2 -a -u f -x splice -t '$NUM_THREADS' -G '$MAX_INTRON' '$GENOMEFILE' "$3" 2>tissue"n".err| grep -v \047^@\047) | samtools view -bhS /dev/stdin > tissue"n"_lr.bam.tmp && mv tissue"n"_lr.bam.tmp tissue"n"_lr.bam && samtools sort -@ '$NUM_THREADS' -m 1G tissue"n"_lr.bam -T sts -O bam> tissue"n"_lr.bam.sorted.tmp.bam && mv tissue"n"_lr.bam.sorted.tmp.bam tissue"n"_lr.bam.sorted.bam && rm tissue"n".bam tissue"n"_lr.bam;\nif [ -s tissue"n"_lr.bam.sorted.bam ];then\n ./run_stringtie_mix.sh tissue"n".bam.sorted.bam tissue"n"_lr.bam.sorted.bam;\nelse\n ./run_stringtie.sh tissue"n".bam.sorted.bam;\nfi\nfi\nfi";
          n++;
        }
      }else if(NF == 3){
        if($NF == "fasta"){
          print "if [ ! -s tissue"n".bam.sorted.bam.gtf ];then\nif [ ! -s "$1" ];then echo \"Input data file "$1" does not exist or empty!\";exit 1;fi\nFIRSTCHAR=`zcat -f "$1" | head -n 1 | cut -b 1`\nif [ $FIRSTCHAR = \">\" ];then\n hisat2 -x '$GENOME'.hst -f -p '$NUM_THREADS' --min-intronlen 20 --max-intronlen '$MAX_INTRON' -1 "$1" -2 "$2" 2>tissue"n".err | samtools view -bhS /dev/stdin > tissue"n".bam.tmp && mv tissue"n".bam.tmp tissue"n".bam && samtools sort -@ '$NUM_THREADS' -m 1G tissue"n".bam  -T sts -O bam >tissue"n".bam.sorted.tmp.bam && mv tissue"n".bam.sorted.tmp.bam tissue"n".bam.sorted.bam && rm tissue"n".bam && ./run_stringtie.sh tissue"n".bam.sorted.bam || exit 1\nelse echo \"WARNING! Invalid fasta format files "$1" or "$2", ignoring them\"\nfi\nfi"; 
          n++;
        }else if($NF == "fastq"){
          print "if [ ! -s tissue"n".bam.sorted.bam.gtf ];then\nif [ ! -s "$1" ];then echo \"Input data file "$1" does not exist or empty!\";exit 1;fi\nFIRSTCHAR=`zcat -f "$1" | head -n 1 | cut -b 1`\nif [ $FIRSTCHAR = \"@\" ];then\n hisat2 -x '$GENOME'.hst -p '$NUM_THREADS' --min-intronlen 20 --max-intronlen '$MAX_INTRON' -1 "$1" -2 "$2" 2>tissue"n".err | samtools view -bhS /dev/stdin > tissue"n".bam.tmp && mv tissue"n".bam.tmp tissue"n".bam && samtools sort -@ '$NUM_THREADS' -m 1G tissue"n".bam  -T sts -O bam >tissue"n".bam.sorted.tmp.bam && mv tissue"n".bam.sorted.tmp.bam tissue"n".bam.sorted.bam && rm tissue"n".bam && ./run_stringtie.sh tissue"n".bam.sorted.bam || exit 1\nelse echo \"WARNING! Invalid fastq format files "$1" or "$2", ignoring them\"\nfi\nfi";
          n++;
        }else if($NF == "bam_mix"){
          print "if [ ! -s tissue"n".bam.sorted.bam.gtf ];then\nif [ ! -s "$1" ];then echo \"Input data file "$1" does not exist or empty!\";exit 1;fi\ncp "$1" tissue"n".bam.tmp && mv tissue"n".bam.tmp tissue"n".bam && cp "$2" tissue"n"_lr.bam.tmp && mv tissue"n"_lr.bam.tmp tissue"n"_lr.bam && samtools sort -@ '$NUM_THREADS' -m 1G tissue"n".bam  -T sts -O bam>tissue"n".bam.sorted.tmp.bam && mv tissue"n".bam.sorted.tmp.bam tissue"n".bam.sorted.bam && samtools sort -@ '$NUM_THREADS' -m 1G tissue"n"_lr.bam  -T sts -O bam >tissue"n"_lr.bam.sorted.tmp.bam && mv tissue"n"_lr.bam.sorted.tmp.bam tissue"n"_lr.bam.sorted.bam && rm tissue"n".bam tissue"n"_lr.bam && ./run_stringtie_mix.sh tissue"n".bam.sorted.bam tissue"n"_lr.bam.sorted.bam || exit 1\nfi";
          n++;
        }
      }else if(NF == 2){
        if($NF == "bam"){
          print "if [ ! -s tissue"n".bam.sorted.bam.gtf ];then\nif [ ! -s "$1" ];then echo \"Input data file "$1" does not exist or empty!\";exit 1;fi\n cp "$1" tissue"n".bam.tmp && mv tissue"n".bam.tmp tissue"n".bam && samtools sort -@ '$NUM_THREADS' -m 1G tissue"n".bam -T sts -O bam >tissue"n".bam.sorted.tmp.bam && mv tissue"n".bam.sorted.tmp.bam tissue"n".bam.sorted.bam && rm tissue"n".bam && ./run_stringtie.sh tissue"n".bam.sorted.bam || exit 1\nfi"; 
          n++;
        }else if($NF == "bam_isoseq"){
          print "if [ ! -s tissue"n".bam.sorted.bam.gtf ];then\nif [ ! -s "$1" ];then echo \"Input data file "$1" does not exist or empty!\";exit 1;fi\n cp "$1" tissue"n".bam.tmp && mv tissue"n".bam.tmp tissue"n".bam && samtools sort -@ '$NUM_THREADS' -m 1G tissue"n".bam  -T sts -O bam >tissue"n".bam.sorted.tmp.bam && mv tissue"n".bam.sorted.tmp.bam tissue"n".bam.sorted.bam && rm tissue"n".bam && ./run_stringtie_lr.sh tissue"n".bam.sorted.bam || exit 1\nfi"; 
          n++;
        }else if($NF == "isoseq"){
          print "if [ ! -s tissue"n".bam.sorted.bam.gtf ];then\nif [ ! -s "$1" ];then echo \"Input data file "$1" does not exist or empty!\";exit 1;fi\nFIRSTCHAR=`zcat -f "$1" | head -n 1 | cut -b 1`\nif [ $FIRSTCHAR = \">\" ] || [ $FIRSTCHAR = \"@\" ];then\n cat <(ufasta sizes -H '$GENOMEFILE' | awk \047{print \"@SQ\\tSN:\"$1\"\\tLN:\"$2}\047) <(minimap2 -a -u f -x splice:hq -t '$NUM_THREADS' -G '$MAX_INTRON' '$GENOMEFILE' "$1" 2>tissue"n".err | grep -v \047^@\047) | samtools view -bhS /dev/stdin > tissue"n".bam.tmp && mv tissue"n".bam.tmp tissue"n".bam && samtools sort -@ '$NUM_THREADS' -m 1G tissue"n".bam  -T sts -O bam >tissue"n".bam.sorted.tmp.bam && mv tissue"n".bam.sorted.tmp.bam tissue"n".bam.sorted.bam && rm tissue"n".bam && ./run_stringtie_lr.sh tissue"n".bam.sorted.bam || exit 1\nelse echo \"WARNING! Invalid fasta format file "$1", ignoring it\"\nfi\nfi";
          n++;
        }else if($NF == "fasta"){
          print "if [ ! -s tissue"n".bam.sorted.bam.gtf ];then\nif [ ! -s "$1" ];then echo \"Input data file "$1" does not exist or empty!\";exit 1;fi\nFIRSTCHAR=`zcat -f "$1" | head -n 1 | cut -b 1`\nif [ $FIRSTCHAR = \">\" ];then\n hisat2 -x '$GENOME'.hst -f -p '$NUM_THREADS' --min-intronlen 20 --max-intronlen '$MAX_INTRON' -U "$1" 2>tissue"n".err | samtools view -bhS /dev/stdin > tissue"n".bam.tmp && mv tissue"n".bam.tmp tissue"n".bam && samtools sort -@ '$NUM_THREADS' -m 1G tissue"n".bam  -T sts -O bam >tissue"n".bam.sorted.tmp.bam && mv tissue"n".bam.sorted.tmp.bam tissue"n".bam.sorted.bam && rm tissue"n".bam && ./run_stringtie.sh tissue"n".bam.sorted.bam || exit 1\nelse echo \"WARNING! Invalid fasta format file "$1", ignoring it\"\nfi\nfi";
          n++;
        }else if($NF == "fastq"){
          print "if [ ! -s tissue"n".bam.sorted.bam.gtf ];then\nif [ ! -s "$1" ];then echo \"Input data file "$1" does not exist or empty!\";exit 1;fi\nFIRSTCHAR=`zcat -f "$1" | head -n 1 | cut -b 1`\nif [ $FIRSTCHAR = \"@\" ];then\n hisat2 -x '$GENOME'.hst -p '$NUM_THREADS' --min-intronlen 20 --max-intronlen '$MAX_INTRON' -U "$1" 2>tissue"n".err | samtools view -bhS /dev/stdin > tissue"n".bam.tmp && mv tissue"n".bam.tmp tissue"n".bam && samtools sort -@ '$NUM_THREADS' -m 1G tissue"n".bam  -T sts -O bam >tissue"n".bam.sorted.tmp.bam && mv tissue"n".bam.sorted.tmp.bam tissue"n".bam.sorted.bam && rm tissue"n".bam && ./run_stringtie.sh tissue"n".bam.sorted.bam || exit 1\nelse echo \"WARNING! Invalid fastq format file "$1", ignoring it\"\nfi\nfi";
          n++;
        }else{
          print "if [ ! -s tissue"n".bam.sorted.bam.gtf ];then\nif [ ! -s "$1" ];then echo \"Input data file "$1" does not exist or empty!\";exit 1;fi\nFIRSTCHAR=`zcat -f "$1" | head -n 1 | cut -b 1`\nif [ $FIRSTCHAR = \"@\" ];then\n hisat2 -x '$GENOME'.hst -p '$NUM_THREADS' --min-intronlen 20 --max-intronlen '$MAX_INTRON' -1 "$1" -2 "$2" 2>tissue"n".err | samtools view -bhS /dev/stdin > tissue"n".bam.tmp && mv tissue"n".bam.tmp tissue"n".bam && samtools sort -@ '$NUM_THREADS' -m 1G tissue"n".bam  -T sts -O bam >tissue"n".bam.sorted.tmp.bam && mv tissue"n".bam.sorted.tmp.bam tissue"n".bam.sorted.bam && rm tissue"n".bam && ./run_stringtie.sh tissue"n".bam.sorted.bam || exit 1\nelse echo \"WARNING! Invalid fastq format files "$1" or "$2", ignoring them\"\nfi\nfi"
          n++;
        }
      }else if(NF == 1){
        print "if [ ! -s tissue"n".bam.sorted.bam.gtf ];then\nif [ ! -s "$1" ];then echo \"Input data file "$1" does not exist or empty!\";exit 1;fi\nFIRSTCHAR=`zcat -f "$1" | head -n 1 | cut -b 1`\nif [ $FIRSTCHAR = \"@\" ];then\n hisat2 -x '$GENOME'.hst -f -p '$NUM_THREADS' --min-intronlen 20 --max-intronlen '$MAX_INTRON' -U "$1" 2>tissue"n".err | samtools view -bhS /dev/stdin > tissue"n".bam.tmp && mv tissue"n".bam.tmp tissue"n".bam && samtools sort -@ '$NUM_THREADS' -m 1G tissue"n".bam  -T sts -O bam >tissue"n".bam.sorted.tmp.bam && mv tissue"n".bam.sorted.tmp.bam tissue"n".bam.sorted.bam && rm tissue"n".bam && ./run_stringtie.sh tissue"n".bam.sorted.bam || exit 1\nelse echo \"WARNING! Invalid fastq format file "$1", ignoring it\"\nfi\nfi";
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
  log "Aligning and building transcripts from RNAseq reads" && \
  bash ./hisat_stringtie.sh && \
  touch transcripts_assemble.success && \
  rm -f transcripts_merge.success || error_exit "Alignment with HISAT2 or transcript assembly with StringTie failed, please check if reads files exist and formatted correctly"
fi

NUM_TISSUES=`ls tissue*.bam.sorted.bam| grep -v "_lr" |wc -l`
if [ -e tissue0.bam.sorted.bam ];then
  let NUM_TISSUES=$NUM_TISSUES-1;
fi

if [ -e transcripts_assemble.success ] && [ ! -e  transcripts_merge.success ];then
  OUTCOUNT=`ls tissue*.bam.sorted.bam.gtf|wc -l`
  if [ $OUTCOUNT -eq 1 ];then
    # a single tissue, we add 1 for the number of samples and the TPMs
    cat tissue*.bam.sorted.bam.gtf | \
      grep -v '^#' | \
      perl -F'\t' -ane '{
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
      }' > $GENOME.gtf.tmp && mv $GENOME.gtf.tmp $GENOME.merged.gtf && \
      touch transcripts_merge.success && \
      rm -f merge.success || error_exit "Failed to merge transcripts"
  elif [ $OUTCOUNT -ge $NUM_TISSUES ];then
    log "Merging transcripts" && \
    gffcompare -ST tissue*.bam.sorted.bam.gtf  -o $GENOME.tmp -p MSTRG 1>gffcompare.out 2>&1 && \
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
    mv $GENOME.tmp2.combined.gtf $GENOME.merged.gtf && \
    rm -f $GENOME.tmp.{combined.gtf,tracking,loci,redundant.gtf} $GENOME.tmp $GENOME.max_tpm.samples.txt && \
    touch transcripts_merge.success && \
    rm -f merge.success || error_exit "Failed to merge transcripts"
  else
    error_exit "one or more Stringtie jobs failed to run properly"
  fi
fi

if [ ! -e protein2genome.deduplicate.success ];then
  if [ -s $PROTEINFILE ];then 
    log "Deduplicating input proteins"
    ufasta one $PROTEINFILE | \
      awk '{if($0 ~ /^>/){header=$1}else{print header,$1}}' |\
      sort  -S 10% -k2,2 |\
      uniq -f 1 |\
      awk '{print $1"\n"$2}' | \
      tr ':' '_' > $PROTEIN.uniq.tmp && \
    mv $PROTEIN.uniq.tmp $PROTEIN.uniq && \
    touch protein2genome.deduplicate.success && \
    rm -f protein2genome.align.success || error_exit "Failed in deduplicating proteins"
  fi
fi

if [ ! -e protein2genome.align.success ];then
  if [ -s $PROTEINFILE ];then
    log "Aligning proteins to the genome with miniprot"
    #we may need a bigger k for big genomes
    KMERVALUE=`ls -lL $GENOMEFILE | perl -ane '{if($F[4]>1000000000){print "6";}else{print "5"}}'` && \
    miniprot -p 0.95 -N 20 -k $KMERVALUE -t $NUM_THREADS -G $MAX_INTRON --gff $GENOMEFILE $PROTEIN.uniq 2>miniprot.err | \
      convert_miniprot_gff.pl > $GENOME.$PROTEIN.uniq.palign.gff.tmp && \
    mv $GENOME.$PROTEIN.uniq.palign.gff.tmp $GENOME.$PROTEIN.uniq.palign.gff
  else
    touch $GENOME.$PROTEIN.uniq.palign.gff
  fi
  rm -f merge.success && \
  touch protein2genome.align.success || error_exit "Alignment of proteins to the genome with miniprot failed, please check miniprot.err"
fi

if [ -e transcripts_merge.success ] && [ -e protein2genome.align.success ] && [ ! -e merge.success ];then
  log "Deriving gene models from protein and transcript alignments" && \
  if [ ! -s $GENOME.merged.gtf ];then
    error_exit "No transcripts useful for annotation, please check your inputs!"
  fi && \
#we fix suspect introns in the protein alignment files.  If an intron has never been seen before, switch it to the closest one that has been seen
  if [ -s $CDSFILE ] && [ -s $GENOME.$PROTEIN.uniq.palign.gff ];then
    log "Using external CDSs and protein alignments" && \
    perl -F'\t' -ane 'next if($F[0] =~/^#/);$F[6]="+" if(not($F[6] eq "+") && not($F[6] eq "-")); print join("\t",@F);' $CDSFILE |\
      gffread -C -F | \
      perl -F'\t' -ane '{chomp($F[8]);if($F[2] eq "mRNA" || $F[2] eq "transcript"){$pos=($F[4]+$F[3])/2;@f=split(/;/,$F[8]);($junk,$id)=split(/=/,$f[0]);$id.=":$F[0]:$pos"."_EXTERNAL";}elsif(uc($F[2]) eq "CDS"){$F[2]=uc($F[2]); print join("\t",@F[0..7]),"\tParent=$id\n";$F[2]="exon";print join("\t",@F[0..7]),"\tParent=$id\n";}}' |\
      gffread -F | \
      perl -F'\t' -ane '{next if($F[0] =~/^#/);chomp($F[8]);if($F[2] eq "transcript"){$F[2]="gene";$F[8].=";gene$F[8];identity=100.00;similarity=100.00";}print join("\t",@F),"\n";}' > $CDS.CDS.tmp &&
    mv $CDS.CDS.tmp $CDS.CDS && \
    cat $CDS.CDS <(gffread -F  <( fix_suspect_introns.pl $GENOME.merged.gtf < $GENOME.$PROTEIN.uniq.palign.gff )) >  $GENOME.palign.fixed.gff.tmp && \
    mv $GENOME.palign.fixed.gff.tmp $GENOME.palign.fixed.gff && \
    perl -F'\t' -ane 'next if($F[0] =~/^#/);$F[6]="+" if(not($F[6] eq "+") && not($F[6] eq "-")); print join("\t",@F);' $CDSFILE | gffread -y $PROTEIN.cds -g $GENOMEFILE && \
    cat $PROTEIN.uniq $PROTEIN.cds > $PROTEIN.all.tmp && \
    mv $PROTEIN.all.tmp $PROTEIN.all && \
    PROTEINFILE=$PROTEIN.all 
  elif [ -s $GENOME.$PROTEIN.uniq.palign.gff ];then
    log "Using protein alignments" && \
    gffread -F  <( fix_suspect_introns.pl $GENOME.merged.gtf < $GENOME.$PROTEIN.uniq.palign.gff ) > $GENOME.palign.fixed.gff.tmp && \
    mv $GENOME.palign.fixed.gff.tmp $GENOME.palign.fixed.gff 
  elif [ -s $CDSFILE ];then
    log "Using external CDSs only" && \
    perl -F'\t' -ane 'next if($F[0] =~/^#/);$F[6]="+" if(not($F[6] eq "+") && not($F[6] eq "-")); print join("\t",@F);' $CDSFILE |\
      gffread -C -F | \
      perl -F'\t' -ane '{chomp($F[8]);if($F[2] eq "mRNA" || $F[2] eq "transcript"){$pos=($F[4]+$F[3])/2;@f=split(/;/,$F[8]);($junk,$id)=split(/=/,$f[0]);$id.=":$F[0]:$pos"."_EXTERNAL";}elsif(uc($F[2]) eq "CDS"){$F[2]=uc($F[2]); print join("\t",@F[0..7]),"\tParent=$id\n";$F[2]="exon";print join("\t",@F[0..7]),"\tParent=$id\n";}}' |\
      gffread -F | \
      perl -F'\t' -ane '{next if($F[0] =~/^#/);chomp($F[8]);if($F[2] eq "transcript"){$F[2]="gene";$F[8].=";gene$F[8];identity=100.00;similarity=100.00";}print join("\t",@F),"\n";}' > $GENOME.palign.fixed.gff.tmp && \
    mv $GENOME.palign.fixed.gff.tmp $GENOME.palign.fixed.gff && \
    perl -F'\t' -ane 'next if($F[0] =~/^#/);$F[6]="+" if(not($F[6] eq "+") && not($F[6] eq "-")); print join("\t",@F);' $CDSFILE |gffread -y $PROTEIN.cds -g $GENOMEFILE && \
    PROTEINFILE=$PROTEIN.cds 
  fi
#here we use the "fixed" protein alignments as reference and compare our transcripts. This annotates each transcript with a protein match and a match code
  gffcompare -T -o $GENOME.protref -r $GENOME.palign.fixed.gff $GENOME.merged.gtf && \
  #assign_class_code.pl <(trmap $GENOME.palign.fixed.gff $GENOME.merged.gtf -o /dev/stdout) < $GENOME.merged.gtf > $GENOME.protref.annotated.gtf && \
  rm -f $GENOME.protref.{loci,tracking,stats} $GENOME.protref && \
#here we fix missing orientations in the GTF transcripts file
  perl -F'\t' -ane '{
    if($F[2] eq "transcript" && $F[8]=~/^transcript_id "(\S+)"; gene_id/){
      $ori{$1}=$F[6];
    }
  }END{
    open(FILE,"'$GENOME'.merged.gtf");
    while($line=<FILE>){
      chomp($line);
      @f=split(/\t/,$line);
      if($f[6] eq "."){
        $f[8]=~/transcript_id "(\S+)";/;
        $f[6]=$ori{$1} if(defined($ori{$1}));
      }
      print join("\t",@f)."\n" unless($f[6] eq ".");
    }
  }' $GENOME.protref.annotated.gtf > $GENOME.gtf.tmp && \
  mv $GENOME.gtf.tmp $GENOME.gtf && \
#here we combine the transcripts and protein matches
#unused proteins gff file contains all protein alignments that did not match the transcripts; we will use them later
#this produces files $GENOME.{k,u}.gff.tmp  and $GENOME.unused_proteins.gff.tmp
  cat $GENOME.palign.fixed.gff | \
    combine_gene_protein_gff.pl \
      --prefix $GENOME \
      --annotated <( cat $GENOME.palign.fixed.gff | filter_by_class_code.pl $GENOME.protref.annotated.gtf | gffread -F ) \
      --genome $GENOMEFILE \
      --pwms <(echo "") \
      --transdecoder  <(echo "") \
      1>combine.out 2>&1 && \
  mv $GENOME.k.gff.tmp $GENOME.k.gff && \
  extract_utr_transcripts.pl 0 < $GENOME.k.gff > $GENOME.utrs.gff.tmp && \
  mv $GENOME.utrs.gff.tmp $GENOME.utrs.gff && \
  perl -F'\t' -ane '{if($F[2] eq "mRNA"){$protid="$2:$1" if($F[8]=~/EvidenceProteinID=(\S+);EvidenceTranscriptID=(\S+);StartCodon=/);$F[2]="gene";$F[8]="ID=$protid\n";unless(defined($out{$protid})){$gene{$protid}=join("\t",@F);$cds{$protid}="";$beg{$protid}=0;$end{$protid}=0;$ori{$protid}=$F[6];$flag=1;$out{$protid}=1;}else{$flag=0}}elsif($F[2] eq "CDS" && $flag){$beg{$protid}=$F[3] if($beg{$protid}==0); $end{$protid}=$F[4] if($end{$protid}<$F[4]);$F[8]="Parent=$protid\n";$cds{$protid}.=join("\t",@F);$F[2]="exon";$gene{$protid}.=join("\t",@F);}}END{foreach $p(keys %gene){@f=split(/\t/,$gene{$p});$f[3]=$beg{$p};$f[4]=$end{$p}; print join("\t",@f),"$cds{$p}"}}'  $GENOME.k.gff > $GENOME.cds.gff.tmp && \
  mv $GENOME.cds.gff.tmp $GENOME.cds.gff && \
  log "Computing Markov chain matrices at splice junctions" && \
  gffread -F --tlf $GENOME.k.gff |\
    perl -F'\t' -ane 'BEGIN{$n=1}{if($F[8]=~/^ID=(\S+);exonCount=(\S+);exons=(\S+);(.+);EvidenceTranscriptID=(\S+);StartCodon=(.+);Class==;/){$tid=$4;$exons=$3;@f=split(/-/,$exons);for($i=1;$i<$#f;$i++){($c1,$c2)=split(/,/,$f[$i]);$c2--; unless(defined($output{"$F[0]\t$c1\t$c2\t$F[6]"})){print "$F[0]\t$c1\t$c2\tJUNC$n\t1\t$F[6]\n"; $output{"$F[0]\t$c1\t$c2\t$F[6]"}=1;$n++}}}}' | \
    tee >(wc -l > $GENOME.num_introns.txt) | \
    compute_junction_scores_bed.pl $GENOMEFILE 1>$GENOME.pwm.tmp 2>$GENOME.pwm.err && \
  mv $GENOME.pwm.tmp $GENOME.pwm && \
  score_transcripts_with_hmms.pl <(gffread -F $GENOME.gtf) $GENOMEFILE $GENOME.pwm <(cat <(gffread -F --tlf $GENOME.gtf |    perl -F'\t' -ane 'BEGIN{$n=1}{if($F[8]=~/^ID=(\S+);exonCount=(\S+);exons=(\S+);geneID=/){$tid=$1;$exons=$3;@f=split(/-/,$exons);for($i=1;$i<$#f;$i++){($c1,$c2)=split(/,/,$f[$i]);$c2--; unless(defined($output{"$F[0]\t$c1\t$c2\t$F[6]"})){print "$F[0]\t$c1\t$c2\tJUNC\t1\t$F[6]\n"; $output{"$F[0]\t$c1\t$c2\t$F[6]"}=1;$n++}}}}' ) <(gffread -F --tlf $GENOME.palign.fixed.gff |    perl -F'\t' -ane 'BEGIN{$n=1}{if($F[8]=~/^ID=(\S+);exonCount=(\S+);exons=(\S+);CDS=/){$tid=$1;$exons=$3;@f=split(/-/,$exons);for($i=1;$i<$#f;$i++){($c1,$c2)=split(/,/,$f[$i]);$c2--; unless(defined($output{"$F[0]\t$c1\t$c2\t$F[6]"})){print "$F[0]\t$c1\t$c2\tJUNC\t1\t$F[6]\n"; $output{"$F[0]\t$c1\t$c2\t$F[6]"}=1;$n++}}}}')) > $GENOME.transcript_splice_scores.txt.tmp && \
  mv $GENOME.transcript_splice_scores.txt.tmp $GENOME.transcript_splice_scores.txt && \
  perl -F'\t' -ane '{if($F[8] =~ /^transcript_id "(\S+)"; gene_id "(\S+)"; xloc "(\S+)"; cmp_ref "(\S+)"; class_code "(k|=|c)"; tss_id/){print "$1 $4 $5\n"}}' $GENOME.protref.annotated.gtf > $GENOME.reliable_transcripts_proteins.txt && \
  #we now filter the transcripts file using the splice scores, leaving alone the transcripts that do match proteins and rerun combine
  gffcompare -T -r $GENOME.palign.fixed.gff $GENOME.utrs.gff -o $GENOME.readthrough1 && \
  gffcompare -T -r $GENOME.cds.gff $GENOME.utrs.gff -o $GENOME.readthrough2 && \
  cat <(perl -F'\t' -ane '{if($F[2] eq "transcript"){$flag=0;if($F[8]=~/^transcript_id "(\S+)";(.+)class_code "="/){@f=split(/:/,$1);$suffix=substr($f[-1],-2);$tid="$f[0].$suffix:$f[1]:".substr($f[2],0,-3); $F[8]="transcript_id \"$tid\"; gene_id \"$tid\"";$flag=1;print join("\t",@F),"\n"}}elsif($F[2] eq "exon" && $flag){$F[8]="transcript_id \"$tid\"; gene_id \"$tid\"";print join("\t",@F),"\n"}}' $GENOME.readthrough1.annotated.gtf ) \
      <(cat <(detect_readthrough_exons.pl $GENOME.palign.fixed.gff < $GENOME.readthrough1.annotated.gtf) \
            <(detect_readthrough_exons.pl $GENOME.cds.gff < $GENOME.readthrough2.annotated.gtf) | \
        sort -S 5% |\
        uniq | \
        remove_readthrough_exons.pl $GENOME.gtf) | \
    gffread --tlf | \
    perl -F'\t' -ane 'BEGIN{
      open(FILE,"'$GENOME'.num_introns.txt");
      $num_introns=int(<FILE>);
      $index=2;
      $index++ if($num_introns >= 1024);
      $index++ if($num_introns >= 4096);
      open(FILE,"'$GENOME'.transcript_splice_scores.txt");
      while($line=<FILE>){
        chomp($line);
        @f=split(/\s+/,$line);
        $score{$f[0]}=$f[$index];
        $ex_score{$f[0]}=$f[1];
      }
      open(FILE,"'$GENOME'.reliable_transcripts_proteins.txt");
      while($line=<FILE>){
        chomp($line);
        @f=split(/\s+/,$line);
        $reliable{$f[0]}=int('$JUNCTION_THRESHOLD');
      }
    }{
      chomp($F[8]);
      if($F[8] =~ /^ID=(\S+);exonCount=(\S+);exons=(\S+);geneID=(\S+)/){
        $id=$1;
        $exons=$3;
        $geneid=$4;
        ($name,$samples,$tpm)=split(/:/,$id);
        $score{$id}=int('$JUNCTION_THRESHOLD')+1 if(not(defined($score{$id})));
        $score{$id}+=2 if($tpm > 10||$samples>1);
        $score{$id}+=2 if($ex_score{$id}>'$JUNCTION_THRESHOLD');
        $score{$id}+=$reliable{$id};
        if($score{$id}>int('$JUNCTION_THRESHOLD')){
          @f=split(/-/,$exons);
          if($#f>1){
            $transcripts{join("-",@f[1..$#f-1])}.="$id $F[0] $F[6] $f[0] $f[-1] $geneid ";
          }elsif(not(defined($output{"$F[0] $F[6] $exons"}))){
            $output{"$F[0] $F[6] $exons"}=1;
            print;
          }
        }
      }
    }END{
      foreach $c(keys %transcripts){
        @ec=split(/,/,$c);
        $ec=$#ec+1;
        $start=10000000000;
        $end=0;
        @tr=split(/\s/,$transcripts{$c});
        for($i=0;$i<$#tr;$i+=6){
          $start=$tr[$i+3] if($tr[$i+3]<$start);
          $end=$tr[$i+4] if($end<$tr[$i+4]);
        }
        print "$tr[1]\tStringTie\ttranscript\t$start\t$end\t.\t$tr[2]\t.\tID=$tr[0];exonCount=$ec;exons=$start-$c-$end;geneID=$tr[5]\n";
      }
    }' | gffread -T  > $GENOME.spliceFiltered.gtf.tmp && \
  mv $GENOME.spliceFiltered.gtf.tmp $GENOME.spliceFiltered.gtf && \
  if [ $(wc -l $GENOME.spliceFiltered.gtf|awk '{print $1}') -eq 0 ];then error_exit "Transcript file is empty, likely insufficient RNA-seq data, exiting...";fi && \
#we compare and combine filtered proteins and transcripts files
  gffcompare -T -o $GENOME.protref.spliceFiltered -r $GENOME.palign.fixed.gff $GENOME.spliceFiltered.gtf && \
  #assign_class_code.pl <(trmap $GENOME.palign.fixed.gff $GENOME.spliceFiltered.gtf -o /dev/stdout) < $GENOME.spliceFiltered.gtf > $GENOME.protref.spliceFiltered.annotated.gtf && \
  cat $GENOME.palign.fixed.gff | \
    filter_by_class_code.pl $GENOME.protref.spliceFiltered.annotated.gtf | \
    filter_by_local_abundance.pl $GENOME.k.gff 1> $GENOME.transcripts_to_keep.txt.tmp 2>/dev/null && \
  mv $GENOME.transcripts_to_keep.txt.tmp $GENOME.transcripts_to_keep.txt && \
  mv $GENOME.protref.spliceFiltered.annotated.gtf $GENOME.protref.spliceFiltered.annotated.gtf.bak && \
  perl -F'\t' -ane 'BEGIN{open(FILE,"'$GENOME'.transcripts_to_keep.txt");while($line=<FILE>){chomp($line);$h{$line}=1}}{if($F[8]=~/transcript_id \"(\S+)\";/){print if(defined($h{$1}));}}' $GENOME.protref.spliceFiltered.annotated.gtf.bak > $GENOME.protref.spliceFiltered.annotated.gtf.tmp && \
  mv $GENOME.protref.spliceFiltered.annotated.gtf.tmp $GENOME.protref.spliceFiltered.annotated.gtf && \
  rm -f $GENOME.protref.spliceFiltered.annotated.gtf.bak && \
#this run of combine gives us the unused proteins and unused transcripts, we do not care about everything else
  cat $GENOME.palign.fixed.gff | \
    combine_gene_protein_gff.pl \
      --prefix $GENOME \
      --annotated <( gffread -F $GENOME.protref.spliceFiltered.annotated.gtf ) \
      --genome $GENOMEFILE \
      --transdecoder <(echo "") \
      --pwms $GENOME.pwm \
      1>combine.out 2>&1 && \
  mv $GENOME.u.gff.tmp $GENOME.u.gff && \
  mv $GENOME.unused_proteins.gff.tmp $GENOME.unused_proteins.gff && \
  mv $GENOME.k.gff.tmp $GENOME.k.gff && \
#here we detect and remove readthrough trancripts that were not trimmed previously
  perl -F'\t' -ane '{unless($F[2] eq "gene" || $F[0]=~/^#/){if($F[8] =~ /^ID=(\S+);Parent=(\S+);EvidenceProteinID=(\S+);EvidenceTranscriptID=(\S+);StartCodon/){$F[8]="ID=$4";$tid=$4;}else{$F[8]="Parent=$tid"}print join("\t",@F),"\n"}}' $GENOME.k.gff |\
    gffread --cluster-only |\
    detect_readthroughs.pl > $GENOME.readthroughs.txt.tmp && \
  mv $GENOME.readthroughs.txt.tmp $GENOME.readthroughs.txt && \
  echo -n "Found readthrough transcripts: " && wc -l $GENOME.readthroughs.txt && \
  gffread --ids $GENOME.transcripts_to_keep.txt $GENOME.spliceFiltered.gtf | \
    gffread -T --nids $GENOME.readthroughs.txt > $GENOME.abundanceFiltered.spliceFiltered.gtf.tmp && \
  mv $GENOME.abundanceFiltered.spliceFiltered.gtf.tmp $GENOME.abundanceFiltered.spliceFiltered.gtf && \
#here we process proteins that did not match to any transcripts -- we derive CDS-based transcripts from them
  if [ -s $GENOME.unused_proteins.gff ];then
    log "Filtering unused protein only loci" && \
    score_transcripts_with_hmms.pl <(perl -F'\t' -ane '$F[2]="transcript" if($F[2] eq "gene");print join("\t",@F);' $GENOME.unused_proteins.gff) $GENOMEFILE $GENOME.pwm > $GENOME.protein_splice_scores.txt && \
    perl -F'\t' -ane 'BEGIN{
      open(FILE,"'$GENOME'.num_introns.txt");
      $index = 2;
      $num_introns = int(<FILE>);
      $index++ if($num_introns >= 1024);
      $index++ if($num_introns >= 4096);
      open(FILE,"'$GENOME'.protein_splice_scores.txt");
      while($line=<FILE>){
        chomp($line);
        @f=split(/\s+/,$line);
        $score{$f[0]}=$f[$index];
        $ex_score{$f[0]}=$f[1];
      }
    }{
      if($F[2] eq "gene"){
        $id=$1 if($F[8] =~ /^ID=(\S+);geneID/);
        $flag=($id =~/_EXTERNAL$/ || ($score{$id}>'$JUNCTION_THRESHOLD' && $ex_score{$id}>0)) ? 1 : 0;
      }
      print if($flag);
    }' $GENOME.unused_proteins.gff > $GENOME.unused_proteins.spliceFiltered.gff.tmp && \
    mv $GENOME.unused_proteins.spliceFiltered.gff.tmp $GENOME.unused_proteins.spliceFiltered.gff && \
    gffread --cluster-only --tlf $GENOME.unused_proteins.spliceFiltered.gff | \
      filter_unused_proteins.pl \
        $GENOMEFILE \
        $GENOME.unused_proteins.spliceFiltered.gff \
        <(gffread -V -x /dev/stdout -g $GENOMEFILE $GENOME.unused_proteins.spliceFiltered.gff | \
          ufasta one | \
          awk '{if($1 ~ /^>/){name=substr($1,2)}else{if(length($1)%3 == 0) names[$1]=names[$1]" "name}}END{for(a in names){n=split(names[a],b," ");print b[1]" "n}}'\
        ) \
        $GENOME.k.gff \
      > $GENOME.best_unused_proteins.gff.tmp && \
    mv $GENOME.best_unused_proteins.gff.tmp $GENOME.best_unused_proteins.gff && \
    echo -n "Candidate CDS-only transcripts: " && \
    awk -F'\t' '{if($3=="transcript") print;}' $GENOME.best_unused_proteins.gff |wc -l
  fi 
  #these are u's -- no match to a protein, use transdecoder to try to find CDS
  if [ -s $GENOME.u.gff ];then
    log "Looking for ORFs in transcripts with no protein matches"
    gffread -g $GENOMEFILE -w $GENOME.lncRNA.fa $GENOME.u.gff && \
    rm -rf $GENOME.lncRNA.fa.transdecoder* && \
    TransDecoder.LongOrfs -S -t $GENOME.lncRNA.fa 1>transdecoder.LongOrfs.out 2>&1 && \
    TransDecoder.Predict -t $GENOME.lncRNA.fa --single_best_only  1>transdecoder.Predict.out 2>&1
    if [ -s $GENOME.lncRNA.fa.transdecoder.gff3 ];then
      add_cds_to_gff.pl <(awk -F '\t' 'BEGIN{flag=0}{if($3=="gene"){if($9~/ORF type:complete/){flag=1}else{flag=0}}if(flag){print}}' $GENOME.lncRNA.fa.transdecoder.gff3) <  $GENOME.u.gff | \
      gffread -C | \
      perl -F'\t' -ane '{$F[2]="gene" if($F[2] eq "mRNA"); print join("\t",@F);}' > $GENOME.u.cds.gff.tmp && \
      mv $GENOME.u.cds.gff.tmp $GENOME.u.cds.gff 
    fi
    rm -rf $GENOME.lncRNA.fa $GENOME.lncRNA.u.blastp pipeliner.*.cmds $GENOME.lncRNA.fa.transdecoder_dir  $GENOME.lncRNA.fa.transdecoder_dir.__checkpoints $GENOME.lncRNA.fa.transdecoder_dir.__checkpoints_longorfs transdecoder.LongOrfs.out $GENOME.lncRNA.fa.transdecoder.{cds,pep} blastp1.out transdecoder.Predict.out 
  fi
  log "Working on final merge"
#here we combine all transcripts, adding CDSs that did not match any transcript to the transcripts file
  if [ -s $GENOME.best_unused_proteins.gff ];then
    gffread -T $GENOME.best_unused_proteins.gff $GENOME.abundanceFiltered.spliceFiltered.gtf  > $GENOME.all.combined.gtf.tmp && \
    mv $GENOME.all.combined.gtf.tmp $GENOME.all.combined.gtf
  else
    cp $GENOME.abundanceFiltered.spliceFiltered.gtf $GENOME.all.combined.gtf.tmp && \
    mv $GENOME.all.combined.gtf.tmp $GENOME.all.combined.gtf
  fi
#now we have additional proteins produced by transdecoder, let's use them all, along with SNAP proteins that match the transcripts
  if [ -s $GENOME.u.cds.gff ];then
    gffread $GENOME.palign.fixed.gff $GENOME.u.cds.gff >  $GENOME.palign.all.gff
  else
    gffread $GENOME.palign.fixed.gff >  $GENOME.palign.all.gff 
  fi
# the file $GENOME.palign.all.gff contains all CDSs we need to use
  gffcompare -T -o $GENOME.protref.all -r $GENOME.palign.all.gff $GENOME.all.combined.gtf && \
  #assign_class_code.pl <(trmap $GENOME.palign.all.gff $GENOME.all.combined.gtf -o /dev/stdout) < $GENOME.all.combined.gtf > $GENOME.protref.all.annotated.gtf.tmp && mv $GENOME.protref.all.annotated.gtf.tmp $GENOME.protref.all.annotated.gtf && \
  log "Checking for and repairing broken ORFs" && \
  cat $GENOME.palign.all.gff | \
    filter_by_class_code.pl $GENOME.protref.all.annotated.gtf | \
    gffread -F > $GENOME.protref.all.annotated.class.gff.tmp && \
  mv $GENOME.protref.all.annotated.class.gff.tmp $GENOME.protref.all.annotated.class.gff && \
  cat $GENOME.palign.all.gff | \
    check_cds.pl $GENOME $GENOME.protref.all.annotated.class.gff $GENOMEFILE 1>check_cds.out 2>&1 && \
  mv $GENOME.good_cds.fa.tmp $GENOME.good_cds.fa && \
  mv $GENOME.broken_cds.fa.tmp $GENOME.broken_cds.fa && \
  mv $GENOME.broken_ref.txt.tmp $GENOME.broken_ref.txt && \
  if [ -s $GENOME.u.cds.gff ];then
    cat $PROTEINFILE <(gffread -y /dev/stdout -g $GENOMEFILE $GENOME.u.cds.gff) | \
      ufasta extract -f $GENOME.broken_ref.txt > $GENOME.broken_ref.faa.tmp && \
      mv $GENOME.broken_ref.faa.tmp $GENOME.broken_ref.faa
  else
    ufasta extract -f $GENOME.broken_ref.txt $PROTEINFILE > $GENOME.broken_ref.faa.tmp && \
    mv $GENOME.broken_ref.faa.tmp $GENOME.broken_ref.faa
  fi
  rm -rf $GENOME.broken_cds.fa.transdecoder* && \
  TransDecoder.LongOrfs -S -t $GENOME.broken_cds.fa -m 75 1>transdecoder.LongOrfs.out 2>&1 && \
  makeblastdb -in $GENOME.broken_ref.faa -input_type fasta -dbtype prot -out broken_ref 1>makeblastdb.out 2>&1 && \
  blastp -query $GENOME.broken_cds.fa.transdecoder_dir/longest_orfs.pep -db broken_ref  -max_target_seqs 1 -outfmt 6  -evalue 0.000001 -num_threads $NUM_THREADS 2>blastp2.out > $GENOME.broken_cds.blastp.tmp && \
  mv $GENOME.broken_cds.blastp.tmp $GENOME.broken_cds.blastp && \
  TransDecoder.Predict -t $GENOME.broken_cds.fa --single_best_only --retain_blastp_hits $GENOME.broken_cds.blastp 1>transdecoder.Predict.out 2>&1
  if [ -s $GENOME.broken_cds.fa.transdecoder.bed ];then
    perl -F'\t' -ane 'BEGIN{open(FILE,"'$GENOME'.broken_cds.blastp");while($line=<FILE>){@f=split(/\t/,$line);$h{$f[0]}=1;}}{if($F[3]=~/ID=(\S+);GENE/){$id=$1;print "$F[0] $F[6] $F[7]\n" if(defined($h{$id}) && $#F>7 && not($F[3] =~ /ORF_type:internal/));}}' $GENOME.broken_cds.fa.transdecoder.bed > $GENOME.fixed_cds.txt.tmp && \
    mv $GENOME.fixed_cds.txt.tmp $GENOME.fixed_cds.txt
  fi
  rm -rf transdecoder.Predict.out pipeliner.*.cmds $GENOME.broken_cds.fa.transdecoder_dir  $GENOME.broken_cds.transdecoder_dir.__checkpoints $GENOME.broken_cds.fa.transdecoder_dir.__checkpoints_longorfs transdecoder.LongOrfs.out $GENOME.broken_cds.fa.transdecoder.{cds,pep,gff3} && \
  cat $GENOME.palign.all.gff | \
    combine_gene_protein_gff.pl \
      --prefix $GENOME \
      --annotated $GENOME.protref.all.annotated.class.gff \
      --genome $GENOMEFILE \
      --transdecoder $GENOME.fixed_cds.txt \
      --pwms $GENOME.pwm \
      --names <(perl -F'\t' -ane '{if($F[2] eq "transcript"){print "$1 $3\n" if($F[8] =~ /transcript_id "(\S+)"; gene_id "(\S+)"; oId "(\S+)";/);}}'  $GENOME.all.combined.gtf) \
      --final_pass \
      --proteins $PROTEINFILE \
      --output_partial $PARTIAL \
      --lncrnamintpm $LNCRNATPM \
      1>combine.out 2>&1 && \
    mv $GENOME.k.gff.tmp $GENOME.k.gff && \
    mv $GENOME.u.gff.tmp $GENOME.u.gff && \
    rm -f $GENOME.unused_proteins.gff.tmp && \
  touch merge.success && rm -f loci.success pseudo_detect.success functional.success || error_exit "Merging transcript and protein evidence failed."
fi

if [ -e merge.success ] && [ ! -e loci.success ];then
  log "Reassigning loci based on coding sequences" && \
  #here we remove readthrough transcripts and reassign loci; we allow up to 6 nt overlap between CDSs to put them into the same locus
  gffread -F --nids \
      <(perl -F'\t' -ane '{unless($F[2] eq "gene" || $F[0]=~/^#/){if($F[8] =~ /^ID=(\S+);Parent=(\S+);EvidenceProteinID=(\S+);EvidenceTranscriptID=(\S+);StartCodon/){$F[8]="ID=$1";$tid=$1;}else{$F[8]="Parent=$tid"}print join("\t",@F),"\n"}}' $GENOME.k.gff | \
    gffread --cluster-only | \
    detect_readthroughs.pl) $GENOME.k.gff |\
    tee $GENOME.k.nort.gff.tmp |\
    perl -F'\t' -ane '{if($F[2] eq "mRNA"){$protid=$1 if($F[8]=~/^ID=(\S+);EvidenceProteinID/);$F[8]="ID=$protid\n";$gene{$protid}=join("\t",@F);$cds{$protid}="";$beg{$protid}=-1;$end{$protid}=0;$ori{$protid}=$F[6];}elsif($F[2] eq "CDS"){if($F[4]-$F[3]>16){$F[3]+=8;$F[4]-=8;}$beg{$protid}=$F[3] if($beg{$protid}==-1); $end{$protid}=$F[4] if($end{$protid}<$F[4]);$F[8]="Parent=$protid\n";$F[2]="exon";$cds{$protid}.=join("\t",@F);}}END{foreach $p(keys %gene){@f=split(/\t/,$gene{$p});$f[3]=$beg{$p};$f[4]=$end{$p}; print join("\t",@f),"$cds{$p}"}}' | \
    gffread --cluster-only |\
    awk -F'\t' '{if($3=="locus"){split($9,a,";");print substr(a[1],5)" "substr(a[2],13)}}' > $GENOME.locus_transcripts.tmp && \
  mv $GENOME.locus_transcripts.tmp $GENOME.locus_transcripts && \
  mv $GENOME.k.nort.gff.tmp $GENOME.k.nort.gff && \
  gffread -F --keep-exon-attrs --keep-genes --sort-alpha <(reassign_transcripts.pl $GENOME.locus_transcripts < $GENOME.k.nort.gff) $GENOME.u.gff|\
    awk -F '\t' '{if($0 ~ /^# gffread/){print "# EviAnn automated annotation"}else{if($3!=prev){counter=1}if($9 ~ /^Parent=/){print $0";ID="substr($9,8)":"$3":"counter;counter++;}else{print $0}prev=$3;}}' > $GENOME.gff.tmp && \
  mv $GENOME.gff.tmp $GENOME.gff  && \
  touch loci.success && rm -f pseudo_detect.success functional.success || error_exit "Merging transcript and protein evidence failed."
fi

#cleanup
if [ $DEBUG -lt 1 ];then
  rm -f $GENOME.num_introns.txt
  rm -f $GENOME.{k,u,unused_proteins}.gff.tmp
  rm -f broken_ref.{pjs,ptf,pto,pot,pdb,psq,phr,pin} makeblastdb.out blastp2.out
  rm -f $GENOME.unused_proteins.gff $GENOME.u.cds.gff $GENOME.unused_proteins.spliceFiltered.gff
  rm -f $GENOME.abundanceFiltered.spliceFiltered.gtf
  rm -f $GENOME.protref.annotated.gtf $GENOME.protref.spliceFiltered.annotated.gtf $GENOME.reliable_transcripts_proteins.txt $GENOME.{transcript,protein}_splice_scores.txt $GENOME.transcripts_to_keep.txt
  rm -f $GENOME.all.{loci,stats,tracking,combined.gtf,redundant.gtf} $GENOME.all 
  rm -f $GENOME.protref.all.{loci,stats,tracking,annotated.class.gff,annotated.gtf} $GENOME.protref.all
  rm -f $GENOME.protref.spliceFiltered.{loci,tracking,stats} $GENOME.protref.spliceFiltered
  rm -rf $GENOME.palign.all.gff $GENOME.good_cds.fa $GENOME.broken_cds.fa $GENOME.broken_ref.{txt,faa} $GENOME.broken_cds.{blastp,fa.transdecoder.bed} $GENOME.fixed_cds.txt
  rm -f $GENOME.utrs.gff  $GENOME.readthrough{1,2}.* $GENOME.readthrough{1,2} $GENOME.locus_transcripts $GENOME.k.nort.gff $GENOME.cds.gff
fi

if [ -e loci.success ] && [ ! -e pseudo_detect.success ];then
  log "Detecting and annotating processed pseudogenes" && \
  gffread -S -g $GENOMEFILE -y $GENOME.proteins.fasta $GENOME.gff && \
  ufasta extract -f <(awk -F'\t' '{if($3=="exon"){split($9,a,";");print substr(a[1],8);}}'  $GENOME.gff|uniq -d) $GENOME.proteins.fasta > $GENOME.proteins.mex.fasta.tmp && \
  mv $GENOME.proteins.mex.fasta.tmp $GENOME.proteins.mex.fasta && \
  ufasta extract -f <(awk -F'\t' '{if($3=="exon"){split($9,a,";");print substr(a[1],8);}}'  $GENOME.gff|uniq -c | \
    perl -ane '{($gene,$junk)=split(/-/,$F[1]);$max_count{$gene}=$F[0] if($max_count{$gene}<$F[0]);$transcripts{$gene}.="$F[1]\n";}END{foreach $k(keys %max_count){if($max_count{$k}==1){print $transcripts{$k}}}}') $GENOME.proteins.fasta > $GENOME.proteins.sex.fasta.tmp && \
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
            $line=~s/gene_biotype=protein_coding/gene_biotype=processed_pseudogene/;
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
    mv $GENOME.pseudo_label.gff.tmp $GENOME.pseudo_label.gff && \
    rm -f $GENOME.functional_note.pseudo_label.gff add_external.success && \
    touch pseudo_detect.success || error_exit "Detection of pseudogenes failed, you can use annotation in $GENOME.gff without pseudo-gene labels"
  fi
  if [ $DEBUG -lt 1 ];then
    rm -rf $GENOME.proteins.mex.p?? $GENOME.proteins.{s,m}ex.fasta  makeblastdb.sex2mex.out blastp5.out $GENOME.sex2mex.blastp 
  fi
fi

if [ -e loci.success ] && [ -e pseudo_detect.success ];then
#produce transcript and protein sequences
  gffread -S -g $GENOMEFILE -y $GENOME.proteins.fasta -w $GENOME.transcripts.fasta $GENOME.pseudo_label.gff && \
#check if we need to do functional annotation
  if [ $FUNCTIONAL -gt 0 ];then
    if [ ! -e functional.success ];then
      log "Performing functional annotation" && \
      if [ ! -s $UNIPROT ];then 
          log "Downloading UniProt database" && \
          wget --no-check-certificate https://ftp.uniprot.org/pub/databases/uniprot/knowledgebase/complete/uniprot_sprot.fasta.gz && gunzip uniprot_sprot.fasta.gz || error_exit "Downloading UniProt database failed, please obtain it manually" && \
          UNIPROT=uniprot_sprot.fasta
      fi
      makeblastdb -in $UNIPROT -input_type fasta -dbtype prot -out uniprot 1>makeblastdb.out 2>&1 && \
      log "Aligning annotated proteins to UniProt proteins" && \
      blastp -db uniprot -query $GENOME.proteins.fasta -out  $GENOME.maker2uni.blastp.tmp -evalue 0.000001 -outfmt 6 -num_alignments 1 -seg yes -soft_masking true -lcase_masking -max_hsps 1 -num_threads $NUM_THREADS 1>blastp4.out 2>&1 && \
      mv $GENOME.maker2uni.blastp.tmp $GENOME.maker2uni.blastp && \
      my_maker_functional_gff $UNIPROT $GENOME.maker2uni.blastp $GENOME.pseudo_label.gff > $GENOME.functional_note.pseudo_label.gff.tmp && mv $GENOME.functional_note.pseudo_label.gff.tmp $GENOME.functional_note.pseudo_label.gff && \
      my_maker_functional_fasta $UNIPROT $GENOME.maker2uni.blastp $GENOME.proteins.fasta > $GENOME.functional_note.proteins.fasta.tmp  && mv $GENOME.functional_note.proteins.fasta.tmp $GENOME.functional_note.proteins.fasta && \
      my_maker_functional_fasta $UNIPROT $GENOME.maker2uni.blastp $GENOME.transcripts.fasta > $GENOME.functional_note.transcripts.fasta.tmp  && mv $GENOME.functional_note.transcripts.fasta.tmp $GENOME.functional_note.transcripts.fasta && \
      rm -rf blastp4.out add_external.success && \
      touch functional.success  || error_exit "Functional annotation failed"
    fi
    log "Output functionally annotated GFF is in $GENOME.functional_note.pseudo_label.gff, proteins are in  $GENOME.functional_note.proteins.fasta, transcripts are in $GENOME.functional_note.transcripts.fasta"
    echo -n "Number of functional protein coding transcripts: ";awk '{if($3=="mRNA")print $0}' $GENOME.functional_note.pseudo_label.gff |grep Similar |wc -l
  else
    log "Output annotation GFF is in $GENOME.pseudo_label.gff, proteins are in  $GENOME.proteins.fasta, transcripts are in $GENOME.transcripts.fasta"
  fi
  echo -n "Number of genes: ";awk -F'\t' '{if($3=="gene")print $0}' $GENOME.pseudo_label.gff |wc -l
  echo -n "Number of protein coding genes: ";awk -F'\t' '{if($3=="gene")print $0}' $GENOME.pseudo_label.gff| grep protein_coding | grep -v  'pseudo=true' | wc -l
  echo -n "Number of processed pseudo gene transcripts: ";awk -F'\t' '{if($3=="mRNA")print $0}' $GENOME.pseudo_label.gff| grep 'pseudo=true' |wc -l
  echo -n "Number of processed pseudo genes: ";awk -F'\t' '{if($3=="gene")print $0}' $GENOME.pseudo_label.gff| grep 'pseudo=true' |wc -l
  echo -n "Number of transcripts: ";awk -F'\t' '{if($3=="mRNA")print $0}' $GENOME.pseudo_label.gff |wc -l
  echo -n "Number of long non-coding RNAs: "; awk -F '\t' '{if($3=="lnc_RNA") print}'  $GENOME.pseudo_label.gff |wc -l
  echo -n "Number of distinct proteins: "; ufasta one $GENOME.proteins.fasta | grep -v '^>' |sort -S 10% |uniq |wc -l
fi

if [ -s $EXTRA_GFF ] && [ -e merge.success ];then
  if [ ! -e add_external.success ];then
    log "Adding extra features from $EXTRA_GFF"
    if [ -e functional.success ];then
      add_features.pl $GENOME.functional_note.pseudo_label.gff $EXTRA_GFF > $GENOME.functional_note.pseudo_label.extra.gff.tmp && mv $GENOME.functional_note.pseudo_label.extra.gff.tmp $GENOME.functional_note.pseudo_label.extra.gff && \
      touch add_external.success
    else
      add_features.pl $GENOME.pseudo_label.gff $EXTRA_GFF > $GENOME.pseudo_label.extra.gff.tmp && mv $GENOME.pseudo_label.extra.gff.tmp $GENOME.pseudo_label.extra.gff && \
      touch add_external.success
    fi
  fi
  if [ -e functional.success ];then
    log "Annotation GFF with external features added is in $GENOME.functional_note.pseudo_label.extra.gff" && \
    echo -n "Number of genes with external: ";awk -F'\t' '{if($3=="gene")print $0}' $GENOME.functional_note.pseudo_label.extra.gff |wc -l
  else
    log "Annotation GFF with external features added is in $GENOME.pseudo_label.extra.gff" && \
    echo -n "Number of genes with external: ";awk -F'\t' '{if($3=="gene")print $0}' $GENOME.pseudo_label.extra.gff |wc -l
  fi
fi
