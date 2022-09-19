#!/bin/bash
MYPATH="`dirname \"$0\"`"
MYPATH="`( cd \"$MYPATH\" && pwd )`"
export PATH=$MYPATH:$PATH;
set -o pipefail
export NUM_THREADS=16
export MAX_INTRON=200000
BATCH_SIZE=1000000
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
kill -9 0
exit 1
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

function usage {
echo "Usage:"
echo "protein2genome.sh [arguments]"
echo "-a <assembly contigs or scaffolds>"
echo "-p <proteins fasta file>"
echo "-t <number of threads, default: 16>"
echo "-m <max intron size, default: 200000>"
echo ""
echo "External dependencies: must have ufasta, makeblastdb, tblastn and exonerate available on the PATH"
which tblastn
which makeblastdb
which exonerate
which ufasta
}

#parsing arguments
if [[ $# -eq 0 ]];then
usage
exit 1
fi

while [[ $# > 0 ]]
do
    key="$1"
    case $key in
        -t|--threads)
            export NUM_THREADS="$2"
            shift
            ;;
        -a|--assembly)
            export GENOME="$2"
            shift
            ;;
        -p|--proteins)
            export PROTEIN="$2";
            shift
            ;;
        -m|--max_intron)
            export MAX_INTRON="$2";
            shift
            ;;
        -v|--verbose)
            set -x
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

PROTEINN=`basename $PROTEIN`;
GENOMEN=`basename $GENOME`

if [ ! -e protein2genome.protein_align.success ];then
log "Aligning proteins to the genome"
makeblastdb -dbtype nucl -in $GENOME -out $GENOMEN.blastdb 1>/dev/null 2>&1 && \
ufasta sizes -H $PROTEIN | \
perl -ane 'BEGIN{$fn="'$PROTEINN'";$index=1;$batch_size=int("'$BATCH_SIZE'");open(FILE,">$fn.$index.batch");$n=0;}{if($n > $batch_size){close(FILE);$index++;open(FILE,">$fn.$index.batch");$n=$F[1];}else{$n+=$F[1];}print FILE $F[0],"\n";}' && \
NUM_BATCHES=`ls $PROTEINN.*.batch |wc -l` && \
echo "#!/bin/bash" > run_tblastn.sh && \
echo "ufasta extract -f \$1 $PROTEIN > \$1.fa && \\" >>  run_tblastn.sh && \
echo -n "tblastn -db $GENOMEN.blastdb -matrix BLOSUM80 -gapopen 13 -gapextend 2 -max_intron_length $MAX_INTRON -soft_masking true -num_threads " >> run_tblastn.sh && \
echo -n $(($NUM_THREADS/4+1)) >> run_tblastn.sh && \
echo " -outfmt 6 -query \$1.fa  -evalue 1e-8 2>/dev/null | awk '{if(\$3>75) print \$0}' > tblastn.\$1.out && rm -f \$1.fa " >> run_tblastn.sh && \
chmod 0755 run_tblastn.sh && \
ls $PROTEINN.*.batch |xargs -P $NUM_THREADS -I {} ./run_tblastn.sh {} && \
log "Concatenating outputs" && \
cat tblastn.$PROTEINN.*.batch.out | \
sort -k2,2 -k1,1 -k12,12nr -S 10% > $PROTEINN.tblastn.tmp && \
mv $PROTEINN.tblastn.tmp $PROTEINN.tblastn && \
rm -rf  tblastn.$PROTEINN.*.batch.out $PROTEINN.*.batch && \
touch protein2genome.protein_align.success
fi

if [ ! -e protein2genome.filter.success ] && [ -e protein2genome.protein_align.success ];then
log "Filtering protein alignment file"
awk 'BEGIN{protid="";seqid="";mp=0;}{
  coord=($9+$10)/2;
  if($1 != protid || $2 != seqid || ((coord-last_coord)>int("'$MAX_INTRON'") || (last_coord-coord)>int("'$MAX_INTRON'"))){  
    if(mp >0){
      print protid,seqid,mp;
    }
    mp=$3*$4;
    protid=$1;
    seqid=$2;
  }else{
    mp+=$3*$4;
  }
  last_coord=coord;
}END{  
  if(mp >0){
    print protid,seqid,mp;
  }
}' $PROTEINN.tblastn | \
sort -nrk3,3 -S 10% | \
perl -ane '{$h{$F[0]}=$F[1] if(not(defined($h{$F[0]})));}
END{
  $prot="";
  $seq="";
  $padding=5000;

  open(FILE,"'$GENOME'");
  $ctg="";
  $seq="";
  while($line=<FILE>){
    chomp($line);
    if($line=~ /^>/){
      $sequence{$ctg}=$seq if(not($ctg eq ""));
      @f=split(/\s+/,$line);
      $ctg=substr($f[0],1);
      $seq="";
    }else{
      $seq.=$line;
    }
  }
  $sequence{$ctg}=$seq if(not($ctg eq ""));

  open(FILE,"'$PROTEIN'");
  $ctg="";
  $seq="";
  while($line=<FILE>){
    chomp($line);
    if($line=~ /^>/){
      $protsequence{$ctg}=$seq if(not($ctg eq ""));
      @f=split(/\s+/,$line);
      $ctg=substr($f[0],1);
      $seq="";
    }else{
      $seq.=$line;
    }
  }
  $protsequence{$ctg}=$seq if(not($ctg eq ""));

  open(FILE,"'$PROTEINN'.tblastn");
  while($line=<FILE>){
    chomp($line);
    @f=split(/\t/,$line);
    next if(not($h{$f[0]} eq $f[1]));
    if($f[0] eq $prot){
      push(@lines,$line);
    }else{
      if(not($prot eq "")){
        $max_intron=int("'$MAX_INTRON'");
        @lines_filter=();
        @ff=split(/\t/,@lines[0]);
        $min_coord=($ff[8]+$ff[9])/2;
        $max_coord=$min_coord;
        foreach $l (@lines){
          @ff=split(/\t/,$l);
          $coord=($ff[8]+$ff[9])/2;
          $lori="+";
          $lori="-" if($ff[8]>$ff[9]);
          if($lori eq $ori){
            if($coord<=$min_coord){
              if($min_coord-$coord<$max_intron){
                $min_coord=$coord;
                push(@lines_filter,$l);
              }
            }elsif($coord>=$max_coord){
              if($coord-$max_coord<$max_intron){
                $max_coord=$coord;
                push(@lines_filter,$l);
              }
            }else{
              push(@lines_filter,$l);
            }
          }
        }
        @lines_sorted = sort {(split(/\t/,$a))[8] <=> (split(/\t/,$b))[8]} @lines_filter;
        @ff=split(/\t/,$lines_sorted[0]);
        $start = $ff[8] < $ff[9] ? $ff[8] : $ff[9];
        @ff=split(/\t/,$lines_sorted[-1]);
        $end = $ff[8] < $ff[9] ? $ff[9] : $ff[8];
        $start = ($start-$padding)>=0 ? $start-$padding :0;
        open(OUTFILE,">$seq.$prot.taskfile");
        print OUTFILE ">$seq\n",substr($sequence{$seq},$start,$end-$start+$padding),"\n>$prot\n",$protsequence{$prot},"\n#\t$start\t",$end-$start+$padding,"\n" if(defined($sequence{$seq}) && defined($protsequence{$prot}));
        close(OUTFILE);
      }
      @lines=();
      push(@lines,$line);
      $prot=$f[0];
      $seq=$f[1];
      if($f[8]<$f[9]){
        $ori="+";
      }else{
        $ori="-";
      }
    }
  }
}' && touch protein2genome.filter.success 
fi

if [ -e protein2genome.filter.success ] && [ ! -e protein2genome.exonerate_gff.success ];then
log "Running exonerate on the filtered sequences" && \
echo -n '#!/bin/bash
TASKFILE=$1
if [ -s $TASKFILE ] ;then
GENOME=`head -n 1 $TASKFILE | cut -c 2-`
head -n 2 $TASKFILE > /dev/shm/$TASKFILE.fa && \
head -n 4 $TASKFILE |tail -n 2 | tr J I | tr B D | tr Z E > /dev/shm/$TASKFILE.faa && \
PROTLEN=`ufasta sizes /dev/shm/$TASKFILE.faa` && \
tail -n 1 $TASKFILE > $TASKFILE.gff && \
exonerate --model protein2genome -Q protein -T dna -t /dev/shm/$TASKFILE.fa --minintron 10 --maxintron ' > .run_exonerate.sh
echo -n $MAX_INTRON >> .run_exonerate.sh
echo -n ' -q /dev/shm/$TASKFILE.faa --bestn 1 --showtargetgff 2>/dev/null | \
awk '\''BEGIN{flag=0}{if($0 ~ /START OF GFF DUMP/ || $0 ~ /END OF GFF DUMP/){flag++} if(flag==1) print $0}'\'' | \
grep "^$GENOME" >> $TASKFILE.gff && \
ALNLEN=`awk '\''{if($3=="cds") n+=$5-$4}END{print int(n/3)}'\'' $TASKFILE.gff` && \
if [ $ALNLEN -lt $(($PROTLEN-2)) ];then
  tail -n 1 $TASKFILE > $TASKFILE.gff && \
  exonerate --model protein2genome -Q protein --refine full -T dna -t /dev/shm/$TASKFILE.fa --minintron 10 --maxintron ' >> .run_exonerate.sh && \
  echo -n $MAX_INTRON >> .run_exonerate.sh && \
  echo ' -q /dev/shm/$TASKFILE.faa --bestn 1 --showtargetgff 2>/dev/null | \
  awk '\''BEGIN{flag=0}{if($0 ~ /START OF GFF DUMP/ || $0 ~ /END OF GFF DUMP/){flag++} if(flag==1) print $0}'\'' | \
  grep "^$GENOME" >> $TASKFILE.gff
fi 
rm $TASKFILE
fi
rm -rf /dev/shm/$TASKFILE.fa /dev/shm/$TASKFILE.faa ' >> .run_exonerate.sh && \
chmod 0755 .run_exonerate.sh && \
ls |grep .taskfile$ |xargs -P $NUM_THREADS -I {} ./.run_exonerate.sh {} 
rm -f .run_exonerate.sh
touch protein2genome.exonerate_gff.success
fi

if [ -e protein2genome.exonerate_gff.success ] && [ ! -e protein2genome.final_gff_merge.success ];then
log "Merging individual gff files"
ls |grep .taskfile.gff$ | \
xargs cat | \
perl -e '{while($line=<STDIN>){chomp($line);@f=split(/\t/,$line);next if($f[2] eq "similarity"); if($#f==2){$offset=$f[1]}else{print join("\t",@f[0..2]),"\t",$f[3]+$offset,"\t",$f[4]+$offset,"\t",join("\t",@f[5..7]); if($f[2] eq "gene"){@ff=split(/\s+/,$f[-1]);$id=$ff[4];print "\tID=$id;geneID=$id\n";}else{print "\tID=cds-$id;Parent=$id\n";}}}}' > $GENOMEN.$PROTEINN.palign.gff.tmp && mv $GENOMEN.$PROTEINN.palign.gff.tmp $GENOMEN.$PROTEINN.palign.gff && \
ls |grep .taskfile.gff$ | xargs -n 20 rm -f && \
log "Output gff is in $GENOMEN.$PROTEINN.palign.gff"
touch protein2genome.final_gff_merge.success
fi

