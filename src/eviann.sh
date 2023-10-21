#!/bin/bash
#this pipeline generates genome annotation using hisat2, Stringtie2 and maker
GENOME="genome.fa"
PROTEINFILE="$PWD/uniprot_sprot.nonred.85.fasta"
GENOMEFILE="genome.fa"
RNASEQ_PAIRED="paired"
RNASEQ_UNPAIRED="unpaired"
ALT_EST="altest"
export BATCH_SIZE=1000000
export MAX_INTRON=100000
export MIN_TPM=0.25
export DEBUG=0
UNIPROT="$PWD/uniprot_sprot.nonred.85.fasta"
MYPATH="`dirname \"$0\"`"
MYPATH="`( cd \"$MYPATH\" && pwd )`"
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
 echo "Usage: eviann.sh [options]"
 echo "Options:"
 echo "-t <number of threads, default:1>"
 echo "-g <MANDATORY:genome fasta file with full path>"
 echo "-p <file containing list of filenames of paired Illumina reads from RNAseq experiments, one pair of /path/filename per line; fastq is expected by default, if files are fasta, add \"fasta\" as the third field on the line>"
 echo "-u <file containing list of filenames of unpaired Illumina reads from RNAseq experiments, one /path/filename per line; fastq is expected by default, if files are fasta, add \"fasta\" as the third field on the line>"
 echo "-e <fasta file with transcripts from related species>"
 echo "-r <fasta file of protein sequences from related species>"
 echo "-m <max intron size, default: 100000>"
 echo "--debug <debug flag, if used intermediate output files will be kept>"
 echo "-v <verbose flag>"
 echo "--version version"
 echo ""
 echo "-r AND one or more of the -p -u or -e must be supplied."
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
        --version)
            echo "version 1.0.2"
            exit 0
            ;;
        --debug)
            DEBUG=1
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
for prog in $(echo "ufasta hisat2 minimap2 stringtie gffread blastp tblastn makeblastdb gffcompare TransDecoder.Predict TransDecoder.LongOrfs");do
  echo -n "Checking for $prog on the PATH... " && \
  which $prog || error_exit "$prog not found the the PATH, please install the appropriate package";
done

#unpack uniprot
if [ ! -s $UNIPROT ];then
  log "Unpacking Uniprot database"
  gunzip -c $MYPATH/uniprot_sprot.nonred.85.fasta.gz > uniprot_sprot.nonred.85.fasta.tmp && \
  mv uniprot_sprot.nonred.85.fasta.tmp uniprot_sprot.nonred.85.fasta && \
  makeblastdb -in $UNIPROT -input_type fasta -dbtype prot -out uniprot 1>makeblastdb.out 2>&1
else
  makeblastdb -in $UNIPROT -input_type fasta -dbtype prot -out uniprot 1>makeblastdb.out 2>&1
fi

#checking inputs
mkdir -p tttt && cd tttt
if [ ! -s $RNASEQ_PAIRED ] && [ ! -s $RNASEQ_UNPAIRED ]  && [ ! -s $ALT_EST ];then
  cd .. && rm -rf tttt && error_exit "Must specify at least one non-empty file with filenames of RNAseq reads with -p or -u or a file with ESTs from the same or closely related species with -e.  Paths for ALL files must be ABSOLUTE."
fi
if [ ! -s $UNIPROT ];then
  cd  .. && rm -rf tttt && error_exit "File with uniprot sequences is missing or specified improperly, please supply it with -s </path_to/uniprot_file.fa> with an ABSOLUTE Path"
fi
if [ ! -s $PROTEINFILE ];then
  echo "Warning: proteins from related species are not specified, or file $PROTEINFILE is missing using Uniprot proteins as fallback option" && \
  cd .. && \
  export PROTEINFILE=$PWD/uniprot_sprot.nonred.85.fasta && \
  export PROTEIN=uniprot_sprot.nonred.85.fasta && \
  cd tttt 
fi
cd .. && rm -rf tttt

#path to genome does not have to be absolute
if [ ! -s $GENOMEFILE ];then
  cd .. && error_exit "File with genome sequence is missing or specified improperly, please supply it with -g </path_to/genome_file.fa>"
fi

#first we align
if [ ! -e align-build.success ];then
  log "Building HISAT2 index"
  hisat2-build $GENOMEFILE $GENOME.hst 1>hisat2-build.out 2>&1 && touch align-build.success && rm -f align.success
fi

if [ ! -e align.success ];then
  log "Aligning RNAseq reads"
  echo "#!/bin/bash" >hisat2.sh
  if [ -s $ALT_EST ];then
    echo "if [ ! -s tissue0.bam ];then minimap2 -a -u f -x splice:hq -t $NUM_THREADS -G $MAX_INTRON $GENOMEFILE $ALT_EST 2>tissue0.err | $MYPATH/samtools view -bhS /dev/stdin > tissue0.bam.tmp && mv tissue0.bam.tmp tissue0.bam; fi;" >> hisat2.sh
  fi
  if [ -s $RNASEQ_PAIRED ];then
    awk 'BEGIN{n=1}{
      if(NF==3 && $NF == "fasta"){
        print "if [ ! -s tissue"n".bam ];then hisat2 '$GENOME'.hst -f --dta -p '$NUM_THREADS' --min-intronlen 20 --max-intronlen '$MAX_INTRON' -1 "$1" -2 "$2" 2>tissue"n".err | '$MYPATH'/samtools view -bhS /dev/stdin > tissue"n".bam.tmp && mv tissue"n".bam.tmp tissue"n".bam; fi"; n++
      }else{ 
        if(NF==2 || (NF==3 && $NF == "fastq")){
          print "if [ ! -s tissue"n".bam ];then hisat2 '$GENOME'.hst --dta -p '$NUM_THREADS' --min-intronlen 20 --max-intronlen '$MAX_INTRON' -1 "$1" -2 "$2" 2>tissue"n".err | '$MYPATH'/samtools view -bhS /dev/stdin > tissue"n".bam.tmp && mv tissue"n".bam.tmp tissue"n".bam; fi"; n++
        }
      }
    }' $RNASEQ_PAIRED >> hisat2.sh
  fi
  if [ -s $RNASEQ_UNPAIRED ];then
    START=1;
    if [ -s $RNASEQ_PAIRED ];then
      START=`wc -l $RNASEQ_PAIRED | awk '{print int($1)+1}'`
    fi
    awk 'BEGIN{n=int("'$START'");}{
    if(NF==2 && $NF == "fasta"){
      print "if [ ! -s tissue"n".bam ];then hisat2 '$GENOME'.hst -f --dta -p '$NUM_THREADS' --min-intronlen 20 --max-intronlen '$MAX_INTRON' -U "$1" 2>tissue"n".err | '$MYPATH'/samtools view -bhS /dev/stdin > tissue"n".bam.tmp && mv tissue"n".bam.tmp tissue"n".bam; fi"; n++
    }else{
      if(NF==1 || (NF==2 && $NF == "fastq")){
        print "if [ ! -s tissue"n".bam ];then hisat2 '$GENOME'.hst --dta -p '$NUM_THREADS' --min-intronlen 20 --max-intronlen '$MAX_INTRON' -U "$1" 2>tissue"n".err | '$MYPATH'/samtools view -bhS /dev/stdin > tissue"n".bam.tmp && mv tissue"n".bam.tmp tissue"n".bam; fi"; n++
      }
    }}' $RNASEQ_UNPAIRED >> hisat2.sh
  fi
  bash ./hisat2.sh && touch align.success && rm -f sort.success || error_exit "Alignment with HISAT2 failed, please check if reads files exist"
fi

NUM_TISSUES=`ls tissue*.bam| grep -v sorted |wc -l`
if [ -e tissue0.bam ];then
  let NUM_TISSUES=$NUM_TISSUES-1;
fi

if [ ! -e protein_align.success ];then
  log "Aligning proteins"
  eviprot.sh -t $NUM_THREADS -a $GENOMEFILE -p $PROTEINFILE -m $MAX_INTRON
  if [ -s $GENOME.$PROTEIN.palign.gff ];then
    touch protein_align.success
  fi
fi

if [ ! -e sort.success ];then
  log "Sorting alignment files"
  if [ -s tissue0.bam ];then
    #this is a cludge here, we duplicate each alignment in the sam file n=5 times, adding .i suffix to the transcript name, to force stringtie to drop fewer alignments
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
    rm -f tissue0.bam.sorted.bak
  fi
  if [ $NUM_TISSUES -gt 0 ];then
    for f in $(seq 1 $NUM_TISSUES);do
      if [ ! -s tissue$f.bam.sorted.bam ] && [ -s tissue$f.bam ];then
        samtools sort -@ $NUM_THREADS -m 1G tissue$f.bam tissue$f.bam.sorted.tmp && \
        mv tissue$f.bam.sorted.tmp.bam tissue$f.bam.sorted.bam
      fi
    done
  fi
  touch sort.success && rm -f stringtie.success
fi

if [ ! -e stringtie.success ] && [ -e sort.success ];then
  if [ -s tissue0.bam.sorted.bam ];then
    log "Assembling transcripts from related species with Stringtie"
    stringtie -g 1 -m 100 -t -j 1 -l REFSTRG -f 0.01 -c 1 -p $NUM_THREADS tissue0.bam.sorted.bam -o tissue0.bam.sorted.bam.gtf.tmp && \
    mv tissue0.bam.sorted.bam.gtf.tmp tissue0.bam.sorted.bam.gtf
  fi
  if [ $NUM_TISSUES -gt 0 ];then
    log "Assembling transcripts with Stringtie"
    rm -f stringtie_to_assemble.txt && \
    for f in $(seq 1 $NUM_TISSUES);do
      if [ ! -s tissue$f.bam.sorted.bam.gtf ];then
        echo  tissue$f.bam.sorted.bam >> stringtie_to_assemble.txt
      fi
    done
    if [ -s stringtie_to_assemble.txt ];then
      echo "#!/bin/bash" > run_stringtie.sh && \
      echo "stringtie -p 4 \$1 -o \$1.gtf.tmp && \\">>run_stringtie.sh && \
      echo "awk -F '\t' 'BEGIN{flag=0}{if(\$3==\"transcript\"){n=split(\$9,a,\";\");for(i=1;i<=n;i++){if(a[i] ~ /TPM/){ m=split(a[i],b,\"\\\"\");tpm=b[m-1];}else if(a[i] ~ /FPKM/){ m=split(a[i],b,\"\\\"\");fpkm=b[m-1];}}if(fpkm > '\$MIN_TPM' || tpm > '\$MIN_TPM' ) flag=1; else flag=0;}if(flag){print \$0}}' \$1.gtf.tmp > \$1.gtf.filtered.tmp && \\" >> run_stringtie.sh && \
      echo "mv \$1.gtf.filtered.tmp \$1.gtf  && \\" >> run_stringtie.sh && \
      echo "rm -f \$1.gtf.tmp " >> run_stringtie.sh && \
      chmod 0755 run_stringtie.sh && \
      cat stringtie_to_assemble.txt | xargs -P $(($NUM_THREADS/2+1)) -I {} ./run_stringtie.sh {}
    fi
    rm -f run_stringtie.sh stringtie_to_assemble.txt
  fi
  OUTCOUNT=`ls tissue*.bam.sorted.bam.gtf|wc -l`
  if [ $OUTCOUNT -eq 1 ];then
    # a single tissue, we add 1 for the number of samples and the TPMs
    cat tissue*.bam.sorted.bam.gtf | perl -F'\t' -ane '{
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
    }' > $GENOME.gtf.tmp && mv $GENOME.gtf.tmp $GENOME.gtf
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
    rm -f $GENOME.tmp.combined.gtf
  else
    error_exit "one or more Stringtie jobs failed"
  fi
  touch stringtie.success && rm -f merge.success
fi

if [ ! -e merge.success ];then
  log "Deriving gene models from protein and transcript alignments"
#we first estimate the redundancy of proteins -- this can tell us how many species were used for protein references
  gffcompare -p PCONS -SDT $GENOME.$PROTEIN.palign.gff -o count && \
  NUM_PROT_SPECIES=`perl -F'\t' -ane '{if($F[2] eq "transcript"){if($F[8]=~/.+gene_id "(.+)"; oId .+/){$n++;$h{$1}=1;}}}END{print int($n/scalar(keys(%h))),"\n";}' count.combined.gtf` && \
  rm -f count.combined.gtf  count.{loci,stats,tracking} && \
#now we fix suspect intons in the protein alignment files.  If an intron has never been seen before, switch it to the closest one that has been seen
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
  mv $GENOME.protref.annotated.gtf.tmp $GENOME.protref.annotated.gtf && \
  perl -F'\t' -ane 'BEGIN{open(FILE,"'$GENOME'.transcripts_to_keep.txt");while($line=<FILE>){chomp($line);$h{$line}=1}}{if($F[8]=~/transcript_id \"(\S+)\";/){print if(defined($h{$1}));}}' $GENOME.gtf > $GENOME.abundanceFiltered.gtf.tmp && \
  mv $GENOME.abundanceFiltered.gtf.tmp $GENOME.abundanceFiltered.gtf && \
#here we combine the transcripts and protein matches; we will use only protein CDS's that are contained in the transcripts;some transcripts do not get annotated, we only use "=", "k","j" and "u"
#unused proteins gff file contains all protein alignments that did not match the transcripts; we will use them later
#this produces files $GENOME.{k,j,u}.gff.tmp  and $GENOME.unused_proteins.gff.tmp
  cat $GENOME.palign.fixed.gff |  combine_gene_protein_gff.pl $GENOME <( gffread -F $GENOME.protref.annotated.gtf ) $GENOMEFILE 1>combine.out 2>&1 && \
  mv $GENOME.k.gff.tmp $GENOME.k.gff && \
  mv $GENOME.u.gff.tmp $GENOME.u.gff && \
  mv $GENOME.unused_proteins.gff.tmp $GENOME.unused_proteins.gff && \
  if [ ! -e merge.unused.success ];then
  if [ -s $GENOME.unused_proteins.gff ] && [ ! -e merge.unused.success ];then
    log "Filtering unused protein only loci" && \
    gffread --cluster-only <(awk '{if($3=="cds" || $3=="transcript") print $0}' $GENOME.unused_proteins.gff) > $GENOME.unused_proteins.combined.gff && \
    #here we compute the score for each protein -- the score is the alignemtn similarity listed in palign file
    perl -F'\t' -ane '{
      if($F[2] eq "gene"){
        $similarity{$1}=$4*100+$3 if($F[8]=~/^ID=(\S+);geneID=(\S+);identity=(\S+);similarity=(\S+)/ );
      }
    }END{
      open(FILE,"'$GENOME'.unused_proteins.combined.gff");
      while($line=<FILE>){
        chomp($line);
        @f=split(/\t/,$line);
        if($f[2] eq "transcript"){
          undef($transcript_id);
          undef($gene_id);
          if($f[8]=~/ID=(\S+);locus=(\S+)$/){
            $transcript_id=$1;
            $gene_id=$2;
          }
          if(defined($transcript_id) && defined($gene_id)){
            print "$similarity{$transcript_id} $gene_id $transcript_id\n";
          }
        }
      }
    }' $GENOME.$PROTEIN.palign.gff | \
    sort -nrk1,1 -S 10% |\
    perl -ane 'BEGIN{
      $max_prot_at_locus=1;
      $similarity_threshold=9600;
#if NUM_PROT_SPECIES is 2 or less then it look like we are given a protein homology file for a single species
#then we allow for more extra proteins per locus
      if(int('$NUM_PROT_SPECIES')<2){
        $max_prot_at_locus=2;
        $similarity_threshold=8000;
      }
    }{
      if($F[0]>$similarity_threshold){
        if($h{$F[1]} < $max_prot_at_locus || ($h{$F[1]} < $max_prot_at_locus+1 && ($F[0]>$hs{$F[1]}*.98 && $F[0]>9800))){
          $h{$F[1]}+=1;#this is the number of proteins per locus
          $hs{$F[1]}=$F[0] if(not(defined($hs{$F[1]})));#this is the highest score per locus
          $hn{$F[2]}=1;#we mark the proteins to keep
        }
      }
    }END{
      open(FILE,"'$GENOME'.unused_proteins.gff");
      while($line=<FILE>){
        chomp($line);
        @f=split(/=/,$line);
        print "$line\n" if(defined($hn{$f[-1]}));
      }
    }' > $GENOME.best_unused_proteins.gff
  else
    echo "" > $GENOME.best_unused_proteins.gff
  fi 
  fi && touch merge.unused.success && \
  #these are u's -- no match to a protein
  if [ ! -e merge.u.success ];then
  if [ -s $GENOME.u.gff ];then
    log "Looking for ORFs in transcripts with no protein matches"
    gffread -g $GENOMEFILE -w $GENOME.lncRNA.fa $GENOME.u.gff && \
    rm -rf $GENOME.lncRNA.fa.transdecoder* && \
    TransDecoder.LongOrfs -t $GENOME.lncRNA.fa 1>transdecoder.LongOrfs.out 2>&1 && \
    blastp -query $GENOME.lncRNA.fa.transdecoder_dir/longest_orfs.pep -db uniprot  -max_target_seqs 1 -outfmt 6  -evalue 0.000001 -num_threads $NUM_THREADS 2>blastp1.out > $GENOME.lncRNA.blastp.tmp && \
    mv $GENOME.lncRNA.blastp.tmp $GENOME.lncRNA.u.blastp && \
    TransDecoder.Predict -t $GENOME.lncRNA.fa --single_best_only --retain_blastp_hits $GENOME.lncRNA.u.blastp 1>transdecoder.Predict.out 2>&1
    if [ -s $GENOME.lncRNA.fa.transdecoder.gff3 ];then
      add_cds_to_gff.pl $GENOME.lncRNA.fa.transdecoder.gff3 <  $GENOME.u.gff | \
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
              print join("\t",@F);;
            }
          }
        }
      }' > $GENOME.u.cds.gff.tmp && \
      mv $GENOME.u.cds.gff.tmp $GENOME.u.cds.gff
    else
      echo "#gff" > $GENOME.u.cds.gff
    fi
    rm -rf $GENOME.lncRNA.fa pipeliner.*.cmds $GENOME.lncRNA.fa.transdecoder_dir  $GENOME.lncRNA.fa.transdecoder_dir.__checkpoints $GENOME.lncRNA.fa.transdecoder_dir.__checkpoints_longorfs transdecoder.LongOrfs.out $GENOME.lncRNA.fa.transdecoder.{cds,pep,gff3}
  else
    echo "#gff" > $GENOME.u.cds.gff 
  fi
  fi && touch merge.u.success && \
  log "Working on final merge"
  if [ -s $GENOME.best_unused_proteins.gff ];then
#now we have additional proteins produces by transdecoder, let's use them all
    gffcompare -STC $GENOME.best_unused_proteins.gff $GENOME.abundanceFiltered.gtf -o $GENOME.all && \
    rm -f $GENOME.all.{loci,stats,tracking} && \
    gffread $GENOME.palign.fixed.gff $GENOME.u.cds.gff >  $GENOME.palign.all.gff && \
    gffcompare -T -o $GENOME.protref.all -r $GENOME.palign.all.gff $GENOME.all.combined.gtf && \
    rm -f $GENOME.protref.all.{loci,stats,tracking} && \
    log "Checking for and repairing broken ORFs" && \
    cat $GENOME.palign.all.gff | filter_by_class_code.pl $GENOME.protref.all.annotated.gtf | gffread -F > $GENOME.protref.all.annotated.class.gff.tmp && \
    mv $GENOME.protref.all.annotated.class.gff.tmp $GENOME.protref.all.annotated.class.gff && \
    cat $GENOME.palign.all.gff |  check_cds.pl $GENOME $GENOME.protref.all.annotated.class.gff $GENOMEFILE 1>check_cds.out 2>&1 && \
    mv $GENOME.good_cds.fa.tmp $GENOME.good_cds.fa && \
    mv $GENOME.broken_cds.fa.tmp $GENOME.broken_cds.fa && \
    mv $GENOME.broken_ref.txt.tmp $GENOME.broken_ref.txt && \
    ufasta extract -f $GENOME.broken_ref.txt $PROTEINFILE > $GENOME.broken_ref.faa.tmp &&\
    mv $GENOME.broken_ref.faa.tmp $GENOME.broken_ref.faa && \
    rm -rf $GENOME.broken_cds.fa.transdecoder* && \
    TransDecoder.LongOrfs -t $GENOME.broken_cds.fa 1>transdecoder.LongOrfs.out 2>&1 && \
    makeblastdb -in $GENOME.broken_ref.faa -input_type fasta -dbtype prot -out broken_ref 1>makeblastdb.out 2>&1 && \
    blastp -query $GENOME.broken_cds.fa.transdecoder_dir/longest_orfs.pep -db broken_ref  -max_target_seqs 1 -outfmt 6  -evalue 0.000001 -num_threads $NUM_THREADS 2>blastp2.out > $GENOME.broken_cds.blastp.tmp && \
    mv $GENOME.broken_cds.blastp.tmp $GENOME.broken_cds.blastp && \
    TransDecoder.Predict -t $GENOME.broken_cds.fa --single_best_only --retain_blastp_hits $GENOME.broken_cds.blastp 1>transdecoder.Predict.out 2>&1
    if [ -s $GENOME.broken_cds.fa.transdecoder.bed ];then
      awk -F '\t' '{if(NF>8) print $1" "$7" "$8}' $GENOME.broken_cds.fa.transdecoder.bed  > $GENOME.fixed_cds.txt.tmp && \
      mv $GENOME.fixed_cds.txt.tmp $GENOME.fixed_cds.txt
    fi && \
    rm -rf $GENOME.broken_cds.fa pipeliner.*.cmds $GENOME.broken_cds.fa.transdecoder_dir  $GENOME.broken_cds.transdecoder_dir.__checkpoints $GENOME.broken_cds.fa.transdecoder_dir.__checkpoints_longorfs transdecoder.LongOrfs.out $GENOME.broken_cds.fa.transdecoder.{cds,pep,gff3} && \
    cat $GENOME.palign.all.gff |  combine_gene_protein_gff.pl $GENOME $GENOME.protref.all.annotated.class.gff $GENOMEFILE $GENOME.fixed_cds.txt 1>combine.out 2>&1 && \
    mv $GENOME.k.gff.tmp $GENOME.gff && rm -f $GENOME.{u,unused_proteins}.gff.tmp
  else
    gffread $GENOME.palign.fixed.gff $GENOME.u.cds.gff >  $GENOME.palign.all.gff && \
    gffcompare -T -o $GENOME.protref.all -r $GENOME.palign.all.gff $GENOME.abundanceFiltered.gtf && \
    rm -f $GENOME.protref.all.{loci,stats,tracking} && \
    log "Checking for and repairing broken ORFs" && \
    cat $GENOME.palign.all.gff | filter_by_class_code.pl $GENOME.protref.all.annotated.gtf | gffread -F > $GENOME.protref.all.annotated.class.gff.tmp && \
    mv $GENOME.protref.all.annotated.class.gff.tmp $GENOME.protref.all.annotated.class.gff && \
    cat $GENOME.palign.all.gff |  check_cds.pl $GENOME $GENOME.protref.all.annotated.class.gff $GENOMEFILE 1>check_cds.out 2>&1 && \
    mv $GENOME.good_cds.fa.tmp $GENOME.good_cds.fa && \
    mv $GENOME.broken_cds.fa.tmp $GENOME.broken_cds.fa && \
    mv $GENOME.broken_ref.txt.tmp $GENOME.broken_ref.txt && \
    ufasta extract -f $GENOME.broken_ref.txt $PROTEINFILE > $GENOME.broken_ref.faa.tmp &&\
    mv $GENOME.broken_ref.faa.tmp $GENOME.broken_ref.faa && \
    rm -rf $GENOME.broken_cds.fa.transdecoder* && \
    TransDecoder.LongOrfs -t $GENOME.broken_cds.fa 1>transdecoder.LongOrfs.out 2>&1 && \
    makeblastdb -in $GENOME.broken_ref.faa -input_type fasta -dbtype prot -out broken_ref 1>makeblastdb.out 2>&1 && \
    blastp -query $GENOME.broken_cds.fa.transdecoder_dir/longest_orfs.pep -db broken_ref  -max_target_seqs 1 -outfmt 6  -evalue 0.000001 -num_threads $NUM_THREADS 2>blastp3.out > $GENOME.broken_cds.blastp.tmp && \
    mv $GENOME.broken_cds.blastp.tmp $GENOME.broken_cds.blastp && \
    TransDecoder.Predict -t $GENOME.broken_cds.fa --single_best_only --retain_blastp_hits $GENOME.broken_cds.blastp 1>transdecoder.Predict.out 2>&1
    if [ -s $GENOME.broken_cds.fa.transdecoder ];then
      awk -F '\t' '{if(NF>8) print $1" "$7" "$8}' $GENOME.broken_cds.fa.transdecoder.bed  > $GENOME.fixed_cds.txt.tmp && \
      mv $GENOME.fixed_cds.txt.tmp $GENOME.fixed_cds.txt
    fi && \
    rm -rf $GENOME.broken_cds.fa pipeliner.*.cmds $GENOME.broken_cds.fa.transdecoder_dir  $GENOME.broken_cds.transdecoder_dir.__checkpoints $GENOME.broken_cds.fa.transdecoder_dir.__checkpoints_longorfs transdecoder.LongOrfs.out $GENOME.broken_cds.fa.transdecoder.{cds,pep,gff3} && \
    cat $GENOME.palign.all.gff |  combine_gene_protein_gff.pl $GENOME $GENOME.protref.all.annotated.class.gff $GENOMEFILE $GENOME.fixed_cds.txt 1>combine.out 2>&1 && \
    mv $GENOME.k.gff $GENOME.gff && rm -f $GENOME.{u,unused_proteins}.gff.tmp $GENOME.{u,unused_proteins}.gff
  fi && \
  touch merge.success && rm -f functional.success merge.{unused,j,u}.success
fi

if [ -e merge.success ] && [ ! -e functional.success ];then
  log "Performing functional annotation" && \
  gffread -S -g $GENOMEFILE -w $GENOME.transcripts.fasta -y $GENOME.proteins.fasta $GENOME.gff && \
  blastp -db uniprot -query $GENOME.proteins.fasta -out  $GENOME.maker2uni.blastp.tmp -evalue 0.000001 -outfmt 6 -num_alignments 1 -seg yes -soft_masking true -lcase_masking -max_hsps 1 -num_threads $NUM_THREADS 1>blastp4.out 2>&1 && \
  mv $GENOME.maker2uni.blastp.tmp $GENOME.maker2uni.blastp && \
  my_maker_functional_gff $UNIPROT $GENOME.maker2uni.blastp $GENOME.gff > $GENOME.functional_note.gff.tmp && mv $GENOME.functional_note.gff.tmp $GENOME.functional_note.gff && \
  my_maker_functional_fasta $UNIPROT $GENOME.maker2uni.blastp $GENOME.proteins.fasta > $GENOME.functional_note.proteins.fasta.tmp  && mv $GENOME.functional_note.proteins.fasta.tmp $GENOME.functional_note.proteins.fasta && \
  my_maker_functional_fasta $UNIPROT $GENOME.maker2uni.blastp $GENOME.transcripts.fasta > $GENOME.functional_note.transcripts.fasta.tmp  && mv $GENOME.functional_note.transcripts.fasta.tmp $GENOME.functional_note.transcripts.fasta && \
  rm -rf $GENOME.transcripts.fasta $GENOME.proteins.fasta && \
  touch functional.success && rm -rf pseudo_detect.success || error_exit "Functional annotation failed"
fi 

#if [ ! -e partial_detect.success ] && [ -e functional.success ];then
#  log "Eliminating non-functional partial proteins"
#  ufasta one $GENOME.functional_note.proteins.all.fasta | awk '{if($1~/^>/){if($0 ~ "function unknown"){rn=substr($1,2)}else{rn=""}}else if($1 !~ /^M/ && rn != ""){print rn}}' > $GENOME.non_functional_partial_proteins.txt.tmp && \
#  mv $GENOME.non_functional_partial_proteins.txt.tmp $GENOME.non_functional_partial_proteins.txt && \
#  ufasta extract -v -f $GENOME.non_functional_partial_proteins.txt $GENOME.functional_note.proteins.all.fasta > $GENOME.functional_note.proteins.fasta.tmp && \
#  mv $GENOME.functional_note.proteins.fasta.tmp $GENOME.functional_note.proteins.fasta &&\
#  ufasta extract -v -f $GENOME.non_functional_partial_proteins.txt $GENOME.functional_note.transcripts.all.fasta > $GENOME.functional_note.transcripts.fasta.tmp && \
#  mv $GENOME.functional_note.transcripts.fasta.tmp $GENOME.functional_note.transcripts.fasta &&\
#  perl -ane '{$h{$F[0]}=1}END{open(FILE,"'$GENOME'.functional_note.all.gff");while($line=<FILE>){chomp($line);@gff_fields=split(/\t/,$line);if($gff_fields[8]=~/^ID=(XLOC_\d+-mRNA-\d+).+/){print $line,"\n" if(not(defined($h{$1})));}else{print $line,"\n"}}}' $GENOME.non_functional_partial_proteins.txt | \
#  awk -F '\t' 'BEGIN{prev=""}{if($3=="gene"){geneline=$0}else{if(geneline !=""){print geneline;geneline="";}print $0;}}' > $GENOME.functional_note.gff.tmp && \
#  mv $GENOME.functional_note.gff.tmp $GENOME.functional_note.gff && \
#  touch partial_detect.success && \
#  rm -rf pseudo_detect.success && \
#  if [ $DEBUG -lt 1 ];then
#    rm -rf $GENOME.proteins.fasta $GENOME.transcripts.fasta
#  fi
#fi

if [ -e functional.success ] && [ ! -e pseudo_detect.success ];then
  log "Detecting and annotating processed pseudogenes"
  ufasta extract -v -f <(awk '{if($3=="gene" || $3=="exon") print $0" "$3}' $GENOME.gff |uniq -c -f 9  | awk '{if($1==1 && $4=="exon"){split($10,a,":");split(a[1],b,"="); print b[2]}}' ) $GENOME.functional_note.proteins.fasta > $GENOME.proteins.mex.fasta.tmp && mv $GENOME.proteins.mex.fasta.tmp $GENOME.proteins.mex.fasta && \
  ufasta extract -f <(awk '{if($3=="gene" || $3=="exon") print $0" "$3}' $GENOME.gff |uniq -c -f 9  | awk '{if($1==1 && $4=="exon"){split($10,a,":");split(a[1],b,"="); print b[2]}}' ) $GENOME.functional_note.proteins.fasta > $GENOME.proteins.sex.fasta.tmp && mv $GENOME.proteins.sex.fasta.tmp $GENOME.proteins.sex.fasta && \
  makeblastdb -dbtype prot  -input_type fasta -in  $GENOME.proteins.mex.fasta -out $GENOME.proteins.mex 1>makeblastdb.sex2mex.out 2>&1 && \
  blastp -db $GENOME.proteins.mex -query $GENOME.proteins.sex.fasta -out  $GENOME.sex2mex.blastp.tmp -evalue 0.000001 -outfmt "6 qseqid qlen length pident bitscore" -num_alignments 1 -seg yes -soft_masking true -lcase_masking -max_hsps 1 -num_threads $NUM_THREADS 1>blastp5.out 2>&1 && \
  mv $GENOME.sex2mex.blastp.tmp $GENOME.sex2mex.blastp && \
  perl -ane '{if($F[3]>90 && $F[2]/($F[1]+1)>0.90){$pseudo{$F[0]}=1;}}END{open(FILE,"'$GENOME'.functional_note.gff");while($line=<FILE>){chomp($line);@f=split(/\s+/,$line);print $line; ($id,$junk)=split(/;/,$f[8]);if($f[2] eq "gene" && defined($pseudo{substr($id,3)."-mRNA-1"})){ print "pseudo=true;\n";}else{print "\n"}}}' $GENOME.sex2mex.blastp > $GENOME.functional_note.pseudo_label.gff.tmp && \
  mv $GENOME.functional_note.pseudo_label.gff.tmp $GENOME.functional_note.pseudo_label.gff && \
  rm -f $GENOME.functional_note.gff && \
  if [ $DEBUG -lt 1 ];then
    rm -rf $GENOME.proteins.mex.p?? $GENOME.proteins.{s,m}ex.fasta
  fi && \
  touch pseudo_detect.success || error_exit "Detection of pseudogenes failed"
fi

if [ -e functional.success ] && [ -e pseudo_detect.success ];then
  log "Output annotation GFF is in $GENOME.functional_note.pseudo_label.gff, proteins are in  $GENOME.functional_note.proteins.fasta, transcripts are in $GENOME.functional_note.transcripts.fasta"
  echo "Annotation summary:"
  echo -n "Number of genes: ";awk '{if($3=="gene")print $0}' $GENOME.functional_note.pseudo_label.gff |wc -l
  echo -n "Number of functional genes: "; awk '{if($3=="gene")print $0}' $GENOME.functional_note.pseudo_label.gff| grep Similar |wc -l
  echo -n "Number of processed pseudo genes: ";awk '{if($3=="gene")print $0}' $GENOME.functional_note.pseudo_label.gff| grep 'pseudo=true' |wc -l
  echo -n "Number of transcripts: ";awk '{if($3=="mRNA")print $0}' $GENOME.functional_note.pseudo_label.gff |wc -l
  echo -n "Number of functional protein coding transcripts: ";awk '{if($3=="mRNA")print $0}' $GENOME.functional_note.pseudo_label.gff |grep Similar |wc -l
  echo -n "Number of proteins: "; grep '^>' $GENOME.functional_note.proteins.fasta |wc -l
else
  error_exit "Something went wrong, please check your data"
fi



