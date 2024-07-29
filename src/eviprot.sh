#!/bin/bash
MYPATH="`dirname \"$0\"`"
MYPATH="`( cd \"$MYPATH\" && pwd )`"
export PATH=$MYPATH:$PATH;
set -o pipefail
export NUM_THREADS=16
export MAX_INTRON=200000
export MAX_MATCHES=2
export MAX_MATCHES_LOCUS=50
export MATCH_RATIO="0.6"
export MYPID=$$
BATCH_SIZE=250000
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
rm -rf /dev/shm/tmp$MYPID &
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
echo "eviprot.sh [arguments]"
echo "-a <string: assembly contigs or scaffolds>"
echo "-p <string: proteins fasta file>"
echo "-t <int: number of threads, default: 16>"
echo "-m <int: max intron size, default: 200000>"
echo "-n <int: max matches for a protein, default: 2>"
echo "-l <int: max matches at a locus, default: 50>"
echo "--version report version" 
echo "-h|-u this message"
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
        -l|--num_locus)
            export MAX_MATCHES_LOCUS="$2";
            shift
            ;;
        -r|--ratio)
            export MATCH_RATIO="$2";
            shift
            ;;
        -v|--verbose)
            set -x
            ;;
        --version)
            echo "version 1.0.2"
            exit 0
            ;;
        -h|--help|-u|--usage)
            usage
            exit 0
            ;;
        *)
            echo "Unknown option $1"
            usage
            exit 1        # unknown option
            ;;
    esac
    shift
done

#getting absolute paths
GENOME=`realpath $GENOME`
PROTEIN=`realpath $PROTEIN`

#getting filenames
PROTEINN=`basename $PROTEIN`;
GENOMEN=`basename $GENOME`

#remove duplicates from protein sequences
if [ -s $PROTEIN ];then
  ufasta one $PROTEIN | awk '{if($0 ~ /^>/){header=$1}else{print header,$1}}' |sort  -S 10% -k2,2 |uniq -f 1 |awk '{print $1"\n"$2}' > $PROTEINN.uniq.faa && \
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
echo "if [ ! -e tblastn.\$1.out ]; then" >> run_tblastn.sh && \
echo "ufasta extract -f \$1 $PROTEINN > \$1.fa && \\" >>  run_tblastn.sh && \
echo "$MYPATH/tblastn -db $GENOMEN.blastdb -matrix BLOSUM80 -gapopen 13 -gapextend 2 -task tblastn-fast -subject_besthit -max_intron_length $MAX_INTRON -lcase_masking -soft_masking true -num_threads 6 -outfmt 6 -query \$1.fa  -evalue 1e-6 2>/dev/null | awk '{if(\$3>=75) print \$0}' > tblastn.\$1.out.tmp && mv tblastn.\$1.out.tmp tblastn.\$1.out && rm -f \$1.fa " >> run_tblastn.sh && \
echo "fi" >> run_tblastn.sh && \
chmod 0755 run_tblastn.sh && \
ls $PROTEINN.*.batch |xargs -P $(($NUM_THREADS/4+1)) -I {} ./run_tblastn.sh {} 
NUM_BATCHES_DONE=`ls tblastn.$PROTEINN.*.batch.out |wc -l` && \
if [ $NUM_BATCHES -eq $NUM_BATCHES_DONE ];then
  log "Concatenating outputs" && \
  cat tblastn.$PROTEINN.*.batch.out | \
  sort -S 10% > $PROTEINN.tblastn.tmp && \
  mv $PROTEINN.tblastn.tmp ../$PROTEINN.tblastn
else
  cd .. && error_exit "Alignment of proteins failed, please rerun eviprot"
fi
cd .. && \
rm -rf protein_align.tmp && \
touch protein2genome.protein_align.success && rm -f protein2genome.exonerate_gff.success
fi

if [ ! -e protein2genome.exonerate_gff.success ] && [ -e protein2genome.protein_align.success ];then
log "Filtering protein alignment file"
mkdir -p /dev/shm/tmp$MYPID && \
perl -e '{
  my $pathprefix="/dev/shm/tmp'$MYPID'/";
  my $prot="";
  my $seq="";
  my $padding=5000;
  my @filenames=();
  my @filecontents=();
  my @starts=();
  my @ends=();
  my @scores=();

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
  $max_matches=1 if($max_matches < 1);
  $max_matches_locus=int("'$MAX_MATCHES_LOCUS'");
  $max_matches_locus=1 if($max_matches_locus < 1);

  $match_ratio="'$MATCH_RATIO'"+0;
  while($line=<FILE>){
    chomp($line);
    @f=split(/\t/,$line);
    next if(not(defined($sequence{$f[1]})) || not(defined($protsequence{$f[0]})));
    if($f[0] eq $prot){
      push(@lines,$line);
    }else{
      cluster_alignments(@lines) if(not($prot eq ""));
      @lines=();
      push(@lines,$line);
      $prot=$f[0];
    }#if
  } #while
#the last one
  cluster_alignments(@lines) if(not($prot eq ""));
#final output
#first we determine the best matches for each pair of start and end
  my @cluster_scores=();
  my @cluster_indices=();
  for (my $i=0;$i<=$#scores;$i++){
    $cluster_indices{"$starts[$i] $ends[$i]"}.="$i ";
  }
  foreach my $c(keys %cluster_indices){
    my @cindices=split(/\s+/,$cluster_indices{$c});
    my @cindices_sorted=sort {$scores[$b] <=> $scores[$a]} @cindices;
    #print "DEBUG cluster $cluster_indices{$c}\n";
    for (my $i=0;$i<=$#cindices_sorted && $scores[$cindices_sorted[$i]]>=0.8*$scores[$cindices_sorted[0]] && $i<$max_matches_locus;$i++){
      #print "DEBUG $i $scores[$cindices_sorted[$i]] $filenames[$cindices_sorted[$i]]\n";
      open(OUTFILE,">$pathprefix$filenames[$cindices_sorted[$i]]");
      print OUTFILE $filecontents[$cindices_sorted[$i]];
      close(OUTFILE);
    }
  }

  sub cluster_alignments{
    my @input_lines=@_;
#we first sort lines by bitscore in reverse
    my @lines_all_sorted = sort {(split(/\t/,$b))[11] <=> (split(/\t/,$a))[11]} @input_lines;
    my @lines_all_sorted_coord = sort {(split(/\t/,$a))[8] <=> (split(/\t/,$b))[8]} @input_lines;
#print "ALL LINES:\n",join("\n",@lines_all_sorted),"\n";
    my $prev_length=0;
    my $num_matches=0;
    for(my $i=0;$i<=$#lines_all_sorted;$i++){
      next if($lines_all_sorted[$i] eq "");
      next if($num_matches >= $max_matches);
      my @ff=split(/\t/,$lines_all_sorted[$i]);
      #next if($ff[3] < length($protsequence{$ff[0]})*0.1);
#print "NEW center $i $lines_all_sorted[$i]\n";
      $start = $ff[8] < $ff[9] ? $ff[8] : $ff[9];
      $end = $ff[8] < $ff[9] ? $ff[9] : $ff[8];
      my $new_seq=$ff[1];
      my $new_ori="+";
      $new_ori="-" if($ff[8]>$ff[9]);
#here we build a cluster of matches
      $cluster_size=$ff[11];
      my $center_coord=-1;
#we find the coordinate of the starting line in the sorted lines by coordinate and then go find the boundaries in the protein coordnates
      for(my $j=0;$j<=$#lines_all_sorted_coord;$j++){
        if($lines_all_sorted_coord[$j] eq $lines_all_sorted[$i]){
          $center_coord=$j;
          $j=$#lines_all_sorted_coord+1;
        }
      }
      next if($center_coord == -1);
      #going into the negative direction
      if($center_coord>0){
        for(my $j=$center_coord-1;$j>=0;$j--){
          next if($lines_all_sorted_coord[$j] eq "");
          @ffc=split(/\t/,$lines_all_sorted_coord[$j]);
          next if(not($ffc[1] eq $new_seq));
          my $lori="+";
          $lori="-" if($ffc[8]>$ffc[9]);
          next if(not($lori eq  $new_ori));
          if($new_ori eq "+"){
            if($start-$ffc[8]<$max_intron){
              $start=$ffc[8];
              $cluster_size+=$ffc[11];
              $lines_all_sorted_coord[$j]="";
            }else{
              $j=-1;
            }
          }else{
            if($start-$ffc[9]<$max_intron){
              $start=$ffc[9];
              $cluster_size+=$ffc[11];
              $lines_all_sorted_coord[$j]="";
            }else{
              $j=-1;
            }
          }
        }
      }

      #going into the positive direction
      if($center_coord<$#lines_all_sorted_coord){
        for(my $j=$center_coord+1;$j<=$#lines_all_sorted_coord;$j++){
          next if($lines_all_sorted_coord[$j] eq "");
          @ffc=split(/\t/,$lines_all_sorted_coord[$j]);
          next if(not($ffc[1] eq $new_seq));
          my $lori="+";
          $lori="-" if($ffc[8]>$ffc[9]);
          next if(not($lori eq  $new_ori));
          if($new_ori eq "+"){
            if($ffc[9]-$end<$max_intron){
              $end=$ffc[9];
              $cluster_size+=$ffc[11];
              $lines_all_sorted_coord[$j]="";
            }else{
              $j=$#lines_all_sorted_coord+1;
            }   
          }else{
            if($ffc[8]-$end<$max_intron){
              $end=$ffc[8];
              $lines_all_sorted_coord[$j]="";
              $cluster_size+=$ffc[11];
            }else{
              $j=$#lines_all_sorted_coord+1;
            }
          }
        }
      }
      $start = ($start-$padding)>=0 ? $start-$padding :0;
      $prev_length=$cluster_size if($prev_length==0 || $cluster_size>$prev_length);
#print "CLUSTER $cluster_size PREV $prev_length $start $end $new_seq\n";
      if($cluster_size >$prev_length*$match_ratio){
        $num_matches++;
#print "OUTPUT cluster $start $end $new_seq\n";
        push(@filenames,"$new_seq.$prot.$start.taskfile");
        push(@filecontents,">$new_seq\n".substr($sequence{$new_seq},$start,$end-$start+$padding)."\n>$prot:$new_seq:$start\n".$protsequence{$prot}."\n#\t$start\t".($end-$start+$padding)."\n");
        push(@starts,$start);
        push(@ends,$end);
        push(@scores,$cluster_size)
      }
    }#for
  }#sub
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
* -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4  1' > ./blosum80.txt && \
echo -n '#!/bin/bash
TASKFILE=$1
TASKFILEN=`basename $1`
if [ -s $TASKFILE ] && [ ! -s  exonerate_alignments.tmp/$TASKFILEN.gff ];then
GENOME=`head -n 1 $TASKFILE | cut -c 2-`
head -n 2 $TASKFILE > /dev/shm/$TASKFILEN.fa && \
head -n 3 $TASKFILE |tail -n 1 > /dev/shm/$TASKFILEN.faa && \
head -n 4 $TASKFILE |tail -n 1 | tr J I | tr B D | tr Z E >> /dev/shm/$TASKFILEN.faa && \
PROTLEN=`ufasta sizes /dev/shm/$TASKFILEN.faa` && \
tail -n 1 $TASKFILE  > exonerate_alignments.tmp/$TASKFILEN.gff.tmp && \
exonerate --model protein2genome  -Q protein -T dna --refine full -t /dev/shm/$TASKFILEN.fa -f -100 -p ./blosum80.txt --minintron 21 --maxintron ' > run_exonerate.sh
echo -n $MAX_INTRON >> run_exonerate.sh
echo -n ' -q /dev/shm/$TASKFILEN.faa --bestn 1 --showtargetgff --softmasktarget --seedrepeat 10 2>/dev/null |\
awk '\''BEGIN{flag=0}{if($0 ~ /START OF GFF DUMP/ || $0 ~ /END OF GFF DUMP/){flag++} if(flag==1) print $0}'\'' | \
grep "^$GENOME" >>  exonerate_alignments.tmp/$TASKFILEN.gff.tmp && \
rm $TASKFILE && \
mv  exonerate_alignments.tmp/$TASKFILEN.gff.tmp  exonerate_alignments.tmp/$TASKFILEN.gff || rm -f  exonerate_alignments.tmp/$TASKFILEN.gff.tmp
fi
rm -rf /dev/shm/$TASKFILEN.fa /dev/shm/$TASKFILEN.faa ' >> run_exonerate.sh && \
chmod 0755 run_exonerate.sh && \
mkdir -p exonerate_alignments.tmp && \
ls /dev/shm/tmp$MYPID |awk '{print "/dev/shm/tmp'$MYPID'/"$1}' | grep taskfile$ |xargs -P $NUM_THREADS -I {} ./run_exonerate.sh {} 
ls  exonerate_alignments.tmp | grep taskfile.gff$ | awk '{print "exonerate_alignments.tmp/"$1}'| xargs cat | \
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
        $counter=1;
      }elsif($f[2] eq "exon"){
        @ff=split(/\s+/,$f[-1]);
        print "\tID=$f[2]-$id-$counter;Parent=$id;$ff[6]=$ff[7];$ff[9]=$ff[10]";
        print ";$ff[12]=$ff[13]" if($ff[12] eq "frameshifts");
        print "\n";
        $counter++;
      }else{
        print "\tID=$f[2]-$id-$counter;Parent=$id\n";
        $counter++;
      }
    }
  }
}' > $GENOMEN.$PROTEINN.palign.gff.tmp && mv $GENOMEN.$PROTEINN.palign.gff.tmp $GENOMEN.$PROTEINN.palign.gff && \
rm -rf exonerate_alignments.tmp ./run_exonerate.sh blosum80.txt /dev/shm/tmp$MYPID && \
touch protein2genome.exonerate_gff.success && \
log "Output protein alignments are in $GENOMEN.$PROTEINN.palign.gff"
fi

