#!/bin/bash
MYPATH="`dirname \"$0\"`"
MYPATH="`( cd \"$MYPATH\" && pwd )`"
export PATH=$MYPATH:$PATH;
set -o pipefail
export NUM_THREADS=16
export MAX_INTRON=200000
export MAX_MATCHES=2
export MATCH_RATIO="0.6"
export MYPID=$$
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
rm -rf /dev/shm/tmp$MYPID
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
        -n|--num_matches)
            export MAX_MATCHES="$2";
            shift
            ;;
        -r|--ratio)
            export MATCH_RATIO="$2";
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

#remove duplicates from protein sequences
if [ -s $PROTEIN ];then
  ufasta one $PROTEIN | awk '{if($0 ~ /^>/){header=$1}else{print header,$1}}' |sort  -S 10% -k2,2 |uniq -f 1 |awk '{print $1"\n"$2}' > $PROTEINN.uniq.faa
  PROTEIN=$PROTEINN.uniq.faa
else
  error_exit "Query protein file $PROTEIN is empty or does not exist"
fi


if [ ! -e protein2genome.protein_align.success ];then
log "Aligning proteins to the genome"
mkdir -p protein_align.tmp && \
cp $PROTEIN protein_align.tmp/$PROTEINN && \
cd protein_align.tmp && \
makeblastdb -dbtype nucl -in $GENOME -out $GENOMEN.blastdb 1>/dev/null 2>&1 && \
ufasta sizes -H $PROTEINN | \
perl -ane 'BEGIN{$fn="'$PROTEINN'";$index=1;$batch_size=int("'$BATCH_SIZE'");open(FILE,">$fn.$index.batch");$n=0;}{if($n > $batch_size){close(FILE);$index++;open(FILE,">$fn.$index.batch");$n=$F[1];}else{$n+=$F[1];}print FILE $F[0],"\n";}' && \
NUM_BATCHES=`ls $PROTEINN.*.batch |wc -l` && \
echo "#!/bin/bash" > run_tblastn.sh && \
echo "ufasta extract -f \$1 $PROTEINN > \$1.fa && \\" >>  run_tblastn.sh && \
echo -n "tblastn -db $GENOMEN.blastdb -word_size 5 -threshold 19  -matrix BLOSUM80 -gapopen 13 -gapextend 2 -max_intron_length $MAX_INTRON -soft_masking true -num_threads " >> run_tblastn.sh && \
echo -n $(($NUM_THREADS/4+1)) >> run_tblastn.sh && \
echo " -outfmt 6 -query \$1.fa  -evalue 1e-6 2>/dev/null | awk '{if(\$3>=75) print \$0}' > tblastn.\$1.out && rm -f \$1.fa " >> run_tblastn.sh && \
chmod 0755 run_tblastn.sh && \
ls $PROTEINN.*.batch |xargs -P $NUM_THREADS -I {} ./run_tblastn.sh {} 
rm -f $PROTEINN.*.batch $PROTEINN.*.batch.fa run_tblastn.sh && \
log "Concatenating outputs" && \
cat tblastn.$PROTEINN.*.batch.out | \
sort -S 10% > $PROTEINN.tblastn.tmp && \
mv $PROTEINN.tblastn.tmp ../$PROTEINN.tblastn && \
cd .. && \
rm -rf protein_align.tmp && \
touch protein2genome.protein_align.success && rm -f protein2genome.exonerate_gff.success
fi

if [ ! -e protein2genome.exonerate_gff.success ] && [ -e protein2genome.protein_align.success ];then
log "Filtering protein alignment file"
mkdir -p /dev/shm/tmp$MYPID && \
perl -e '{
  $pathprefix="/dev/shm/tmp'$MYPID'/";
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
  $max_intron=int("'$MAX_INTRON'");
  $max_matches=int("'$MAX_MATCHES'");
  $match_ratio="'$MATCH_RATIO'"+0;
  while($line=<FILE>){
    chomp($line);
    @f=split(/\t/,$line,3);
    next if(not(defined($sequence{$f[1]})) || not(defined($protsequence{$f[0]})));
    if($f[0] eq $prot){
      push(@lines,$line);
    }else{
      if(not($prot eq "")){
        @filenames=();
        @filecontents=();
        #we first sort lines by bitscore in reverse
        @lines_all_sorted = sort {(split(/\t/,$b))[11] <=> (split(/\t/,$a))[11]} @lines;
        #print "ALL LINES:\n",join("\n",@lines_all_sorted),"\n";
        $prev_length=0;
        for(my $i=0;$i<=$#lines_all_sorted;$i++){
          next if($lines_all_sorted[$i] eq "");
          @ff=split(/\t/,$lines_all_sorted[$i]);
          @lines_filter=();
          push(@lines_filter,$lines_all_sorted[$i]);
          #print "NEW center $i $lines_all_sorted[$i]\n";
          $lines_all_sorted[$i]="";
          $start = $ff[8] < $ff[9] ? $ff[8] : $ff[9];
          $end = $ff[8] < $ff[9] ? $ff[9] : $ff[8];
          my $new_seq=$ff[1];
          my $new_ori="+";
          $new_ori="-" if($ff[8]>$ff[9]);
          #here we build a cluster of matches
          $cluster_size=$ff[3]*$ff[2];
          for(my $j=0;$j<=$#lines_all_sorted;$j++){
            next if($lines_all_sorted[$j] eq "");
            @ffc=split(/\t/,$lines_all_sorted[$j]);
            next if(not($ffc[1] eq  $new_seq));
            $lori="+";
            $lori="-" if($ffc[8]>$ffc[9]);
            next if(not($lori eq  $new_ori));
	    $startl = $ffc[8] < $ffc[9] ? $ffc[8] : $ffc[9];
	    $endl = $ffc[8] < $ffc[9] ? $ffc[9] : $ffc[8];
            if($endl > $start-$max_intron && $startl < $end+$max_intron){
              $cluster_size+=$ffc[3]*$ffc[2];
              push(@lines_filter,$lines_all_sorted[$j]);
              $lines_all_sorted[$j]="";
	      $start = $startl if($startl < $start);
              $end = $endl if($endl > $end);
            }
          }
          $start = ($start-$padding)>=0 ? $start-$padding :0;
          $prev_length=$cluster_size if($prev_length==0 || $cluster_size>$prev_length);
          #print "CLUSTER $cluster_size PREV $prev_length $start $end $new_seq\n";
          if($cluster_size >$prev_length*$match_ratio){
            #print "OUTPUT cluster $start $end $new_seq\n";
            push(@filenames,"$new_seq.$prot.$start.taskfile");
            push(@filecontents,">$new_seq\n".substr($sequence{$new_seq},$start,$end-$start+$padding)."\n>$prot.$new_seq.$start\n".$protsequence{$prot}."\n#\t$start\t".($end-$start+$padding)."\n");
          }
        }#for
        for(my $k=0;$k<=$#filenames && $k<$max_matches;$k++){
          open(OUTFILE,">$pathprefix$filenames[$k]");
          print OUTFILE $filecontents[$k];
          close(OUTFILE);
        }
      }#if
      @lines=();
      push(@lines,$line);
      $prot=$f[0];
    }#if
  } #while
  #the last one
      if(not($prot eq "")){
        @filenames=();
        @filecontents=();
        #we first sort lines by bitscore in reverse
        @lines_all_sorted = sort {(split(/\t/,$b))[11] <=> (split(/\t/,$a))[11]} @lines;
        #print "ALL LINES:\n",join("\n",@lines_all_sorted),"\n";
        $prev_length=0;
        for(my $i=0;$i<=$#lines_all_sorted;$i++){
          next if($lines_all_sorted[$i] eq "");
          @ff=split(/\t/,$lines_all_sorted[$i]);
          @lines_filter=();
          push(@lines_filter,$lines_all_sorted[$i]);
          #print "NEW center $i $lines_all_sorted[$i]\n";
          $lines_all_sorted[$i]="";
          $start = $ff[8] < $ff[9] ? $ff[8] : $ff[9];
          $end = $ff[8] < $ff[9] ? $ff[9] : $ff[8];
          my $new_seq=$ff[1];
          my $new_ori="+";
          $new_ori="-" if($ff[8]>$ff[9]);
          #here we build a cluster of matches
          $cluster_size=$ff[3]*$ff[2];
          for(my $j=0;$j<=$#lines_all_sorted;$j++){
            next if($lines_all_sorted[$j] eq "");
            @ffc=split(/\t/,$lines_all_sorted[$j]);
            next if(not($ffc[1] eq  $new_seq));
            $lori="+";
            $lori="-" if($ffc[8]>$ffc[9]);
            next if(not($lori eq  $new_ori));
            $startl = $ffc[8] < $ffc[9] ? $ffc[8] : $ffc[9];
            $endl = $ffc[8] < $ffc[9] ? $ffc[9] : $ffc[8];
            if($endl > $start-$max_intron && $startl < $end+$max_intron){
              $cluster_size+=$ffc[3]*$ffc[2];
              push(@lines_filter,$lines_all_sorted[$j]);
              $lines_all_sorted[$j]="";
              $start = $startl if($startl < $start);
              $end = $endl if($endl > $end);
            }
          }
          $start = ($start-$padding)>=0 ? $start-$padding :0;
          $prev_length=$cluster_size if($prev_length==0 || $cluster_size>$prev_length);
          #print "CLUSTER $cluster_size PREV $prev_length $start $end $new_seq\n";
          if($cluster_size >$prev_length*$match_ratio){
            #print "OUTPUT cluster $start $end $new_seq\n";
            push(@filenames,"$new_seq.$prot.$start.taskfile");
            push(@filecontents,">$new_seq\n".substr($sequence{$new_seq},$start,$end-$start+$padding)."\n>$prot.$new_seq.$start\n".$protsequence{$prot}."\n#\t$start\t".($end-$start+$padding)."\n");
          }
        }#for
        for(my $k=0;$k<=$#filenames && $k<$max_matches;$k++){
          open(OUTFILE,">$pathprefix$filenames[$k]");
          print OUTFILE $filecontents[$k];
          close(OUTFILE);
        }
      }#if
}' && \
log "Running exonerate on the filtered sequences" && \
echo '   A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  Z  X  *
A  5 -2 -2 -2 -1 -1 -1  0 -2 -2 -2 -1 -1 -3 -1  1  0 -3 -2  0 -2 -1 -1 -4
R -2  6 -1 -2 -4  1 -1 -3  0 -3 -3  2 -2 -4 -2 -1 -1 -4 -3 -3 -1  0 -1 -4
N -2 -1  6  1 -3  0 -1 -1  0 -4 -4  0 -3 -4 -3  0  0 -4 -3 -4  5  0 -1 -4
D -2 -2  1  6 -4 -1  1 -2 -2 -4 -5 -1 -4 -4 -2 -1 -1 -6 -4 -4  5  1 -1 -4
C -1 -4 -3 -4  9 -4 -5 -4 -4 -2 -2 -4 -2 -3 -4 -2 -1 -3 -3 -1 -4 -4 -1 -4
Q -1  1  0 -1 -4  6  2 -2  1 -3 -3  1  0 -4 -2  0 -1 -3 -2 -3  0  4 -1 -4
E -1 -1 -1  1 -5  2  6 -3  0 -4 -4  1 -2 -4 -2  0 -1 -4 -3 -3  1  5 -1 -4
G  0 -3 -1 -2 -4 -2 -3  6 -3 -5 -4 -2 -4 -4 -3 -1 -2 -4 -4 -4 -1 -3 -1 -4
H -2  0  0 -2 -4  1  0 -3  8 -4 -3 -1 -2 -2 -3 -1 -2 -3  2 -4 -1  0 -1 -4
I -2 -3 -4 -4 -2 -3 -4 -5 -4  5  1 -3  1 -1 -4 -3 -1 -3 -2  3 -4 -4 -1 -4
L -2 -3 -4 -5 -2 -3 -4 -4 -3  1  4 -3  2  0 -3 -3 -2 -2 -2  1 -4 -3 -1 -4
K -1  2  0 -1 -4  1  1 -2 -1 -3 -3  5 -2 -4 -1 -1 -1 -4 -3 -3 -1  1 -1 -4
M -1 -2 -3 -4 -2  0 -2 -4 -2  1  2 -2  6  0 -3 -2 -1 -2 -2  1 -3 -1 -1 -4
F -3 -4 -4 -4 -3 -4 -4 -4 -2 -1  0 -4  0  6 -4 -3 -2  0  3 -1 -4 -4 -1 -4
P -1 -2 -3 -2 -4 -2 -2 -3 -3 -4 -3 -1 -3 -4  8 -1 -2 -5 -4 -3 -2 -2 -1 -4
S  1 -1  0 -1 -2  0  0 -1 -1 -3 -3 -1 -2 -3 -1  5  1 -4 -2 -2  0  0 -1 -4
T  0 -1  0 -1 -1 -1 -1 -2 -2 -1 -2 -1 -1 -2 -2  1  5 -4 -2  0 -1 -1 -1 -4
W -3 -4 -4 -6 -3 -3 -4 -4 -3 -3 -2 -4 -2  0 -5 -4 -4 11  2 -3 -5 -3 -1 -4
Y -2 -3 -3 -4 -3 -2 -3 -4  2 -2 -2 -3 -2  3 -4 -2 -2  2  7 -2 -3 -3 -1 -4
V  0 -3 -4 -4 -1 -3 -3 -4 -4  3  1 -3  1 -1 -3 -2  0 -3 -2  4 -4 -3 -1 -4
B -2 -1  5  5 -4  0  1 -1 -1 -4 -4 -1 -3 -4 -2  0 -1 -5 -3 -4  5  0 -1 -4
Z -1  0  0  1 -4  4  5 -3  0 -4 -3  1 -1 -4 -2  0 -1 -3 -3 -3  0  5 -1 -4
X -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -4
* -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4  1' > /dev/shm/tmp$MYPID/blosum80.txt && \
echo -n '#!/bin/bash
TASKFILE=$1
TASKFILEN=`basename $1`
if [ -s $TASKFILE ] ;then
GENOME=`head -n 1 $TASKFILE | cut -c 2-`
head -n 2 $TASKFILE > /dev/shm/$TASKFILEN.fa && \
head -n 3 $TASKFILE |tail -n 1 > /dev/shm/$TASKFILEN.faa && \
head -n 4 $TASKFILE |tail -n 1 | tr J I | tr B D | tr Z E >> /dev/shm/$TASKFILEN.faa && \
PROTLEN=`ufasta sizes /dev/shm/$TASKFILEN.faa` && \
tail -n 1 $TASKFILE > $TASKFILE.gff && \
exonerate --model protein2genome  -Q protein -T dna -t /dev/shm/$TASKFILEN.fa -f -100 -p ./blosum80.txt --minintron 21 --maxintron ' > /dev/shm/tmp$MYPID/run_exonerate.sh
echo -n $MAX_INTRON >> /dev/shm/tmp$MYPID/run_exonerate.sh
echo -n ' -q /dev/shm/$TASKFILEN.faa --bestn 1 --showtargetgff 2>/dev/null | tee exonerate.out |\
awk '\''BEGIN{flag=0}{if($0 ~ /START OF GFF DUMP/ || $0 ~ /END OF GFF DUMP/){flag++} if(flag==1) print $0}'\'' | \
grep "^$GENOME" >> $TASKFILE.gff && \
ALNLEN=`awk '\''{if($3=="cds") n+=$5-$4}END{print int(n/3)}'\'' $TASKFILE.gff` && \
if [ $ALNLEN -lt $(($PROTLEN-2)) ] || [ $ALNLEN -gt $(($PROTLEN+2)) ];then
  tail -n 1 $TASKFILE > $TASKFILE.gff && \
  exonerate --model protein2genome -Q protein --refine full -T dna -t /dev/shm/$TASKFILEN.fa -f -100 -p ./blosum80.txt --minintron 21 --maxintron ' >> /dev/shm/tmp$MYPID/run_exonerate.sh && \
  echo -n $MAX_INTRON >> /dev/shm/tmp$MYPID/run_exonerate.sh && \
  echo ' -q /dev/shm/$TASKFILEN.faa --bestn 1 --showtargetgff 2>/dev/null | tee exonerate.out |\
  awk '\''BEGIN{flag=0}{if($0 ~ /START OF GFF DUMP/ || $0 ~ /END OF GFF DUMP/){flag++} if(flag==1) print $0}'\'' | \
  grep "^$GENOME" >> $TASKFILE.gff
fi 
rm $TASKFILE
fi
rm -rf /dev/shm/$TASKFILEN.fa /dev/shm/$TASKFILEN.faa ' >> /dev/shm/tmp$MYPID/run_exonerate.sh && \
chmod 0755 /dev/shm/tmp$MYPID/run_exonerate.sh && \
(cd /dev/shm/tmp$MYPID && ls | grep .taskfile$ |xargs -P $NUM_THREADS -I {} ./run_exonerate.sh {} )
(cd /dev/shm/tmp$MYPID && ls | grep .taskfile.gff$ | xargs cat ) | \
perl -e '{
  while($line=<STDIN>){
    chomp($line);
    @f=split(/\t/,$line);
    next if($f[2] eq "similarity"); 
    if($#f==2){
      $offset=$f[1];
    }else{
      print join("\t",@f[0..2]),"\t",$f[3]+$offset,"\t",$f[4]+$offset,"\t",join("\t",@f[5..7]);
      if($f[2] eq "gene"){
        @ff=split(/\s+/,$f[-1]);
        $id=$ff[4];
        print "\tID=$id;geneID=$id;$ff[9]=$ff[10];$ff[12]=$ff[13]\n";
      }else{
        print "\tID=cds-$id;Parent=$id;$ff[6]=$ff[7];$ff[9]=$ff[10]";
        print ";$ff[12]=$ff[13]" if($ff[12] eq "frameshifts");
        print "\n";
      }
    }
  }
}' > $GENOMEN.$PROTEINN.palign.gff.tmp && mv $GENOMEN.$PROTEINN.palign.gff.tmp $GENOMEN.$PROTEINN.palign.gff && \
rm -rf /dev/shm/tmp$MYPID && \
touch protein2genome.exonerate_gff.success && \
log "Output gff is in $GENOMEN.$PROTEINN.palign.gff"
fi

