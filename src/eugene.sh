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
 echo "--rmlib <fasta file with species-specific repeats produced by RepeatModeler>"
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
        --rmlib)
            RMLIB="$2"
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
for prog in $(echo "ufasta hisat2 stringtie2 maker gffread gff3_merge blastp makeblastdb");do
  which $prog
  if [ $? -gt 0 ];then error_exit "$prog not found the the PATH";fi
done

SNAP_PATH=`which maker`
SNAP_PATH=`dirname $SNAP_PATH`
export SNAP_PATH="$SNAP_PATH/../exe/snap"

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

if [ ! -e split.success ];then 
  log "setting up directories and files"
#now to set up parallel maker we need to split fasta files and transcript files split filenames into batches batch.number
  ufasta sizes -H $GENOMEFILE | perl -ane 'BEGIN{$fn="batch";$index=1;$batch_size=int("'$BATCH_SIZE'");open(FILE,">$fn.$index");$n=0;}{if($n > $batch_size){close(FILE);$index++;open(FILE,">$fn.$index");$n=$F[1];}else{$n+=$F[1];}print FILE $F[0],"\n";}' 
 NUM_BATCHES=`ls batch.* |wc -l`
#now we extract all batches from gtf files and set up directories
  for f in $(seq 1 $NUM_BATCHES);do
    rm -rf $f.dir && mkdir -p $f.dir && \
    touch $f.dir/$f.transcripts.fa && \
    ufasta extract -f batch.$f $GENOMEFILE | ufasta format > $f.dir/$f.fa.tmp && mv $f.dir/$f.fa.tmp $f.dir/$f.fa && \
    perl -ane 'BEGIN{open(FILE,"batch.'$f'");while($line=<FILE>){chomp($line);$h{$line}=1;}}{print if defined($h{$F[0]})}'  $GENOME.$PROTEIN.palign.gff | \
    awk -F '\t' '{if($3=="gene" || $3=="exon") print $0}' | \
    sed s,exonerate:protein2genome:local,masurca, | \
    sed s,gene,transcript, > $f.dir/$f.proteins.gff.tmp && \
    mv $f.dir/$f.proteins.gff.tmp $f.dir/$f.proteins.gff
    if [ -e $GENOME.gtf ];then
      #these are the fasta transcripts by Stringtie
      perl -ane 'BEGIN{open(FILE,"batch.'$f'");while($line=<FILE>){chomp($line);$h{$line}=1;}}{print if defined($h{$F[0]})}' $GENOME.gtf | gffread -g $f.dir/$f.fa -w $f.dir/$f.transcripts.fa.tmp /dev/stdin && mv $f.dir/$f.transcripts.fa.tmp $f.dir/$f.transcripts.fa && \
      #if we have transcripts computed from RNAseq data, then convert them from gtf and merge them with transcripts computed from the protein alignments
      perl -ane 'BEGIN{open(FILE,"batch.'$f'");while($line=<FILE>){chomp($line);$h{$line}=1;}}{print if defined($h{$F[0]})}' $GENOME.gtf | \
      gffread -F >  $f.dir/$f.transcripts.gff.tmp && \
      mv $f.dir/$f.transcripts.gff.tmp $f.dir/$f.transcripts.gff && \
      stringtie2 --merge $f.dir/$f.transcripts.gff $f.dir/$f.proteins.gff | \
      gffread -F > $f.dir/$f.evidence.gff.tmp && mv $f.dir/$f.evidence.gff.tmp $f.dir/$f.evidence.gff
    else
      #no RNAseq transcripts -- simply use aligned proteins as evidence
      (cd $f.dir && ln -sf $f.proteins.gff $f.evidence.gff)
    fi
    #extract proteins fasta - -for some reason this is better than simply supplying aligned proteins gff
    ufasta format $PROTEINFILE > $f.dir/$f.proteins.fa.tmp && mv $f.dir/$f.proteins.fa.tmp $f.dir/$f.proteins.fa
    #ufasta extract -f <(awk -F '=' '{print $NF}' $f.dir/$f.proteins.gff) $PROTEINFILE | ufasta format > $f.dir/$f.proteins.fa.tmp && mv $f.dir/$f.proteins.fa.tmp $f.dir/$f.proteins.fa      
    #finally if there are transcripts from related species, add them as fasta ESTs
    if [ -s tissue0.bam ];then
      mv $f.dir/$f.transcripts.fa $f.dir/$f.transcripts.fa.bak && \
      ufasta extract -f <(samtools view tissue0.bam  | perl -ane 'BEGIN{open(FILE,"batch.'$f'");while($line=<FILE>){chomp($line); $h{$line}=1}}{print $F[0],"\n" if(defined($h{$F[2]}));}' ) $ALT_EST | ufasta format > $f.dir/$f.transcripts.fa.tmp && mv $f.dir/$f.transcripts.fa.tmp $f.dir/$f.transcripts.fa && \
      cat $f.dir/$f.transcripts.fa.bak $f.dir/$f.transcripts.fa > $f.dir/$f.transcripts.fa.tmp && mv $f.dir/$f.transcripts.fa.tmp $f.dir/$f.transcripts.fa
    fi
  done
  touch split.success && rm -f maker1.success
fi

if [ ! -e maker1.success ] && [ -e split.success ];then
  log "running maker"
#first we create a script to run
  maker -CTL
  mkdir -p /dev/shm/tmp_$PID
  echo "#!/bin/bash" > run_maker.sh
  echo "cd \$1 && echo \"running maker in \$PWD\" && maker -cpus 4 -base $GENOME 1>maker.log 2>&1" >> run_maker.sh && \
  chmod 0755 run_maker.sh
  NUM_BATCHES=`ls batch.* |wc -l`
#let's create the appropriate maker ctl files
  for f in $(seq 1 $NUM_BATCHES);do
    sed s,^genome=,genome=$PWD/$f.dir/$f.fa, maker_opts.ctl | \
    sed s,^est2genome=0,est2genome=1, | \
    sed s,^protein2genome=0,protein2genome=1, | \
    sed s,^single_exon=0,single_exon=1, | \
    sed s,^single_length=250,single_length=200, | \
    sed s,^TMP=,TMP=/dev/shm/tmp_$PID, | \
    sed s,^max_dna_len=100000,max_dna_len=1000000, | \
    sed s,^split_hit=10000,split_hit=$MAX_INTRON, | \
    sed s,^cpus=1,cpus=4, | \
    sed s,^min_contig=1,min_contig=1000, > $f.dir/maker_opts.ctl 
    #these only appear if related species ests are given
    if [ -s "$f.dir/$f.transcripts.fa" ];then
      mv $f.dir/maker_opts.ctl $f.dir/maker_opts.ctl.bak && \
      sed s,^est=,est=$PWD/$f.dir/$f.transcripts.fa, $f.dir/maker_opts.ctl.bak > $f.dir/maker_opts.ctl && \
      rm $f.dir/maker_opts.ctl.bak
    fi
    #these are proteins split by batch, they are always there
    if [ -s "$f.dir/$f.proteins.fa" ];then
      mv $f.dir/maker_opts.ctl $f.dir/maker_opts.ctl.bak && \
      sed s,^protein=,protein=$PWD/$f.dir/$f.proteins.fa, $f.dir/maker_opts.ctl.bak > $f.dir/maker_opts.ctl && \
      rm $f.dir/maker_opts.ctl.bak
    fi
    if [ -s "$RMLIB" ];then
      mv $f.dir/maker_opts.ctl $f.dir/maker_opts.ctl.bak && \
      sed s,^rmlib=,rmlib=$RMLIB, $f.dir/maker_opts.ctl.bak > $f.dir/maker_opts.ctl && \
      rm $f.dir/maker_opts.ctl.bak
    fi
    #these are the aligned est evidence, they could be either from de novo or from aligned proteins
    #if [ -s "$f.dir/$f.evidence.gff" ];then
    #  mv $f.dir/maker_opts.ctl $f.dir/maker_opts.ctl.bak && \
    #  sed s,^est_gff=,est_gff=$PWD/$f.dir/$f.evidence.gff, $f.dir/maker_opts.ctl.bak > $f.dir/maker_opts.ctl && \
    #  rm $f.dir/maker_opts.ctl.bak
    #fi
    #so far I did not find any advantage to supplying proteins as aligned gff
    #if [ -s "$f.dir/$f.proteins.gff" ];then
    #  mv $f.dir/maker_opts.ctl $f.dir/maker_opts.ctl.bak && \
    #  sed s,^protein_gff=,protein_gff=$PWD/$f.dir/$f.proteins.gff, $f.dir/maker_opts.ctl.bak > $f.dir/maker_opts.ctl && \
    #  rm $f.dir/maker_opts.ctl.bak
    #fi
    cp maker_exe.ctl $f.dir && cp maker_bopts.ctl $f.dir
  done
#running maker
  ls -d *.dir | xargs -P $NUM_THREADS  -I % ./run_maker.sh % && touch maker1.success && rm -rf /dev/shm/tmp_$PID
fi

if [ $SNAP -eq 1 ];then
log "training SNAP for de novo gene finding"
#SNAP to build HMMs from the maker pass1
  if [ ! -e snap1.success ] && [ -e maker1.success ];then
    log "training gene models pass1"
    echo "#!/bin/bash" > run_snap.sh
    echo "cd \$1 && echo \"running SNAP in \$PWD\" && \\" >> run_snap.sh && \
    echo "f=\`echo \$1 | awk -F '.' '{print \$1}'\` && \\" >> run_snap.sh && \
    echo "rm -rf snap && mkdir -p snap && cd snap && \\" >> run_snap.sh && \
    echo "gff3_merge -d ../$GENOME.maker.output/${GENOME}_master_datastore_index.log -o \$f.gff && \\" >> run_snap.sh && \
    echo "maker2zff -c 1 -e 1 -x 0.25 -l 50 \$f.gff && \\" >> run_snap.sh && \
    echo "$SNAP_PATH/fathom -categorize 1000 genome.ann genome.dna && \\" >> run_snap.sh && \
    echo "$SNAP_PATH/fathom -export 1000 -plus uni.ann uni.dna && \\" >> run_snap.sh && \
    echo "$SNAP_PATH/forge export.ann export.dna && \\" >> run_snap.sh && \
    echo "$SNAP_PATH/hmm-assembler.pl \$f . > ../\$f.hmm.tmp && \\" >> run_snap.sh && \
    echo "mv ../\$f.hmm.tmp ../\$f.hmm \\" >> run_snap.sh && \
    chmod 0755 run_snap.sh && \
    ls -d *.dir | xargs -P $NUM_THREADS  -I % ./run_snap.sh % && touch snap1.success && rm -f maker2.success
  fi

#rerun maker with HMMs built by SNAP
  if [ ! -e maker2.success ] && [ -e snap1.success ];then
    log "running maker with gene models"
#first we create a script to run
    maker -CTL
    mkdir -p /dev/shm/tmp_$PID
    echo "#!/bin/bash" > run_maker2.sh
    echo "cd \$1 && echo \"running maker in \$PWD\" && if [ -e $GENOME.maker.output ];then mv $GENOME.maker.output $GENOME.maker.output_pass1; fi && maker -cpus 4 -base $GENOME 1>maker2.log 2>&1" >> run_maker2.sh && \
    chmod 0755 run_maker2.sh
    NUM_BATCHES=`ls batch.* |wc -l`
#let's create the appropriate maker ctl files
    for f in $(seq 1 $NUM_BATCHES);do
      sed s,^genome=,genome=$PWD/$f.dir/$f.fa, maker_opts.ctl | \
      sed s,^maker_gff=,maker_gff=$PWD/$f.dir/snap/$f.gff, | \
      sed s,^snaphmm=,snaphmm=$PWD/$f.dir/$f.hmm, | \
      sed s,^est_pass=0,est_pass=1, | \
      sed s,^protein_pass=0,protein_pass=1, | \
      sed s,^model_pass=0,model_pass=1, | \
      sed s,^repeat_protein=,"repeat_protein= #", | \
      sed s,^model_org=,"model_org= #", | \
      sed s,^rm_pass=0,rm_pass=1, | \
      sed s,^single_exon=0,single_exon=1, | \
      sed s,^single_length=250,single_length=200, | \
      sed s,^TMP=,TMP=/dev/shm/tmp_$PID, | \
      sed s,^max_dna_len=100000,max_dna_len=1000000, | \
      sed s,^split_hit=10000,split_hit=$MAX_INTRON, | \
      sed s,^cpus=1,cpus=4, | \
      sed s,^min_contig=1,min_contig=1000, > $f.dir/maker_opts.ctl && \
      cp maker_exe.ctl $f.dir && \
      cp maker_bopts.ctl $f.dir
    done
#running maker
    ls -d *.dir | xargs -P $NUM_THREADS  -I % ./run_maker2.sh % && touch maker2.success && rm -f functional.success && rm -rf /dev/shm/tmp_$PID
  fi
fi #if SNAP

if [ ! -e functional.success ];then
  log "concatenating outputs"
  NUM_BATCHES=`ls batch.* |wc -l`
  for f in $(seq 1 $NUM_BATCHES);do
    if [ -e $f.dir/$GENOME.maker.output ];then
      (cd $f.dir/$GENOME.maker.output && gff3_merge -g -d ${GENOME}_master_datastore_index.log -o $f.gff && fasta_merge -d ${GENOME}_master_datastore_index.log -o $f.fasta )
    fi
  done
  gff3_merge -o $GENOME.gff.tmp *.dir/$GENOME.maker.output/*.gff && mv $GENOME.gff.tmp $GENOME.gff && \
  cat *.dir/$GENOME.maker.output/*.all.maker.proteins.fasta > $GENOME.proteins.fasta.tmp  && mv $GENOME.proteins.fasta.tmp $GENOME.proteins.fasta && \
  cat *.dir/$GENOME.maker.output/*.all.maker.transcripts.fasta > $GENOME.transcripts.fasta.tmp && mv $GENOME.transcripts.fasta.tmp $GENOME.transcripts.fasta && \
  log "maker output is in $GENOME.gff $GENOME.proteins.fasta $GENOME.transcripts.fasta" && \
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



