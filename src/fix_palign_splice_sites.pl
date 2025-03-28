#!/usr/bin/env perl
#This code adds gene and CDS record for every locus in GFF file
my $geneID="";
my @exons=();
my %transcript_gff;
my $protID="";
my $dir="";
my $scf="";
my $seq="";
my $score_fix_threshold=5;
my $score_threshold=5;
my $fix_range=4;
#these are genetic codes for HMMs
my %code=();
$code{"A"}=0;
$code{"a"}=0;
$code{"C"}=1;
$code{"c"}=1;
$code{"G"}=2;
$code{"g"}=2;
$code{"T"}=3;
$code{"t"}=3;
$code{"N"}=4;
$code{"n"}=4;

#here we load up all transcripts that matched proteins
open(FILE,$ARGV[0]);
while(my $line=<FILE>){
  chomp($line);
  my @gff_fields=split(/\t/,$line);
  my @attributes=split(";",$gff_fields[8]);
  if($gff_fields[2] eq "gene"){
    if(defined($transcript{$geneID})){
      $transcript_gff{$geneID}=[@exons];
    }
    @exons=();
    $locID=substr($attributes[1],7);#this is the gene_id
    $geneID=substr($attributes[0],3);#this is the transcript_id
    $transcript_junction_score{$geneID}=10000;
    $transcript{$geneID}=$line;
  }elsif($gff_fields[2] eq "exon"){
    push(@exons,$line) if(defined($transcript{$geneID}));
  }
}
if(defined($transcript{$geneID})){
  $transcript_gff{$geneID}=[@exons];
}
@exons=();

#we load the genome sequences
open(FILE,$ARGV[1]);
while(my $line=<FILE>){
  chomp($line);
  if($line=~ /^>/){
    if(not($scf eq "")){
      $genome_seqs{$scf}=$seq;
      $seq="";
    }
    my @f=split(/\s+/,$line);
    $scf=substr($f[0],1);
  }else{
    $seq.=$line;
  } 
}   
$genome_seqs{$scf}=$seq if(not($scf eq ""));

#we load SNAP HMMs
if(defined($ARGV[2])){
  open(FILE,$ARGV[2]);
  $line=<FILE>;
  if($line =~ /^zoeHMM/){#check format
    while($line=<FILE>){
      chomp($line);
      if($line=~/^Donor/){
        $line=<FILE>;
        my $i=0;
        while($line=<FILE>){
          last if($line=~/NN TRM/);
          chomp($line);
          $line=~s/^\s+//;
          my @f=split(/\s+/,$line);
          for(my $j=0;$j<=3;$j++){
            $donor_freq[$i][$j]=$f[$j];
          }
          $donor_freq[$i][4]=0;
          $i++;
        }
      }elsif($line=~/^Acceptor/){
        $line=<FILE>;
        my $i=0;
        while($line=<FILE>){
          last if($line=~/NN TRM/);
          chomp($line);
          $line=~s/^\s+//;
          my @f=split(/\s+/,$line);
          for(my $j=0;$j<=3;$j++){
            $acceptor_freq[$i][$j]=$f[$j];
          }
          $acceptor_freq[$i][4]=0;
          $i++;
        }
      }elsif($line=~/^Start/){
        $line=<FILE>;
        my $i=0;
        while($line=<FILE>){
          last if($line=~/NNN TRM/);
          chomp($line);
          $line=~s/^\s+//;
          my @f=split(/\s+/,$line);
          for(my $j=0;$j<=3;$j++){
            $coding_start_freq[$i][$j]=$f[$j];
          }
          $coding_start_freq[$i][4]=0;
          $i++;
        }
      }
    }
  }
}


#we make the transcript sequences for protein coding transcripts and score the transcripts with HMMs
for my $g(keys %transcript_gff){
  my @gff_fields=();
  for(my $j=1;$j<=$#{$transcript_gff{$g}};$j++){
    @gff_fields=split(/\t/,${$transcript_gff{$g}}[$j]);
    die("Genome sequence $gff_fields[0] needed for transcript $g not found!") if(not(defined($genome_seqs{$gff_fields[0]})));
    if($j>0 && defined($ARGV[2])){
      my @gff_fields_prev=split(/\t/,${$transcript_gff{$g}}[$j-1]);
      if($gff_fields[6] eq "+"){
        $donor_seq=uc(substr($genome_seqs{$gff_fields[0]},$gff_fields_prev[4]-3,9));
        $acceptor_seq=uc(substr($genome_seqs{$gff_fields[0]},$gff_fields[3]-28,30));
      }else{
        $donor_seq=uc(substr($genome_seqs{$gff_fields[0]},$gff_fields[3]-7,9));
        $acceptor_seq=uc(substr($genome_seqs{$gff_fields[0]},$gff_fields_prev[4]-3,30));
        $donor_seq=~tr/ACGTNacgtn/TGCANtgcan/;
        $donor_seq=reverse($donor_seq);
        $acceptor_seq=~tr/ACGTNacgtn/TGCANtgcan/;
        $acceptor_seq=reverse($acceptor_seq);
      }
      my $donor_junction_score=score_donor_seq($donor_seq);
      my $acceptor_junction_score=score_acceptor_seq($acceptor_seq);
      my $junction_score=$donor_junction_score+$acceptor_junction_score;
      print "DEBUG $g $gff_fields[6] $donor_seq $acceptor_seq $donor_junction_score $acceptor_junction_score\n";
      if($junction_score<$score_threshold){
        #suggest modification to the junction
        #try to modify donor or acceptor, look for lower score
        #look for alternative GT or AG
        if($donor_junction_score<$acceptor_junction_score){
          print "DEBUG trying to improve donor junction with score $donor_junction_score seq $donor_seq\n";
          if($gff_fields[6] eq "+"){
            for(my $i=-$fix_range;$i<=$fix_range;$i++){
              next if($i==0 || not(uc(substr($genome_seqs{$gff_fields[0]},$gff_fields_prev[4]+$i,2)) eq "GT"));
              $new_donor_seq=uc(substr($genome_seqs{$gff_fields[0]},$gff_fields_prev[4]-3+$i,9));
              $new_donor_junction_score=score_donor_seq($new_donor_seq);
              if($new_donor_junction_score>$donor_junction_score && $acceptor_junction_score+$new_donor_junction_score>=$score_fix_threshold){
                print "DEBUG suggest new donor junction in $g $i $donor_junction_score $new_donor_junction_score $new_donor_seq $donor_seq $gff_fields[6]\n";
                $gff_fields_prev[4]=$gff_fields_prev[4]+$i;
                ${$transcript_gff{$g}}[$j-1]=join("\t",@gff_fields_prev);
                $junction_score=$acceptor_junction_score+$new_donor_junction_score;
                last;
              }
            }
          }else{
            for(my $i=-$fix_range;$i<=$fix_range;$i++){
              next if($i==0 || not(uc(substr($genome_seqs{$gff_fields[0]},$gff_fields[3]+$i-3,2)) eq "AC"));
              $new_donor_seq=uc(substr($genome_seqs{$gff_fields[0]},$gff_fields[3]-7+$i,9));
              $new_donor_seq=~tr/ACGTNacgtn/TGCANtgcan/;
              $new_donor_seq=reverse($new_donor_seq);
              $new_donor_junction_score=score_donor_seq($new_donor_seq);
              if($new_donor_junction_score>$donor_junction_score && $acceptor_junction_score+$new_donor_junction_score>=$score_fix_threshold){
                print "DEBUG suggest new donor junction in $g $i $donor_junction_score $new_donor_junction_score $new_donor_seq $donor_seq $gff_fields[6]\n";
                $gff_fields[3]=$gff_fields[3]+$i;
                ${$transcript_gff{$g}}[$j]=join("\t",@gff_fields);
                $junction_score=$acceptor_junction_score+$new_donor_junction_score;
                last;
              }
            }
          }
        }else{
          if($gff_fields[6] eq "+"){
            for(my $i=-$fix_range;$i<$fix_range;$i++){
              next if($i==0 || not(uc(substr($genome_seqs{$gff_fields[0]},$gff_fields[3]+$i-3,2)) eq "AG"));
              $new_acceptor_seq=uc(substr($genome_seqs{$gff_fields[0]},$gff_fields[3]-28+$i,30));
              $new_acceptor_junction_score=score_acceptor_seq($new_acceptor_seq);
              if($new_acceptor_junction_score>$acceptor_junction_score && $donor_junction_score+$new_acceptor_junction_score>=$score_fix_threshold){
                print "DEBUG suggest new acceptor junction in $g $i $acceptor_junction_score $new_acceptor_junction_score $new_acceptor_seq $acceptor_seq $gff_fields[6]\n";
                $gff_fields[3]=$gff_fields[3]+$i;
                ${$transcript_gff{$g}}[$j]=join("\t",@gff_fields);
                $junction_score=$new_acceptor_junction_score+$donor_junction_score;
                last;
              }
            }
          }else{
            for(my $i=-$fix_range;$i<=$fix_range;$i++){
              next if($i==0 || not(uc(substr($genome_seqs{$gff_fields[0]},$gff_fields_prev[4]+$i,2)) eq "CT"));
              $new_acceptor_seq=uc(substr($genome_seqs{$gff_fields[0]},$gff_fields_prev[4]-3+$i,30));
              $new_acceptor_seq=~tr/ACGTNacgtn/TGCANtgcan/;
              $new_acceptor_seq=reverse($new_acceptor_seq);
              $new_acceptor_junction_score=score_acceptor_seq($new_acceptor_seq);
              if($new_acceptor_junction_score>$acceptor_junction_score && $donor_junction_score+$new_acceptor_junction_score>=$score_fix_threshold){
                print "DEBUG suggest new acceptor junction in $g $i $acceptor_junction_score $new_acceptor_junction_score $new_acceptor_seq $acceptor_seq $gff_fields[6]\n";
                $gff_fields_prev[4]=$gff_fields_prev[4]+$i;
                ${$transcript_gff{$g}}[$j-1]=join("\t",@gff_fields_prev);
                $junction_score=$new_acceptor_junction_score+$donor_junction_score;
                last;
              }
            }
          }
        }
      }
      $transcript_junction_score{$g}=$junction_score if($transcript_junction_score{$g}>$junction_score);
    }
  }
}  

for my $g(keys %transcript_gff){
  print $transcript{$g},"\n";
  for(my $j=0;$j<=$#{$transcript_gff{$g}};$j++){
    print ${$transcript_gff{$g}}[$j],"\n";
  }
  for(my $j=0;$j<=$#{$transcript_gff{$g}};$j++){
    my @ff=split(/\t/,${$transcript_gff{$g}}[$j]);
    print join("\t",@ff[0..1]),"\tcds\t",join("\t",@ff[3..$#ff]),"\n";
  }
}

sub score_donor_seq{
  my $seq=$_[0];
  my $js=0;
  for(my $i=0;$i<9;$i++){
    $js+=$donor_freq[$i][$code{substr($seq,$i,1)}];
  }
  return($js);
}

sub score_acceptor_seq{
  my $seq=$_[0];
  my $js=0;
  for(my $i=0;$i<30;$i++){
    $js+=$acceptor_freq[$i][$code{substr($seq,$i,1)}];
  }
  return($js);
}

