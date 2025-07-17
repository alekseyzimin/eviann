#!/usr/bin/env perl
#this code scores and filters unused protein alignments
my $genome_file=$ARGV[0];
my $unused_proteins_file=$ARGV[1];
my $counts_file=$ARGV[2];
my $k_file=$ARGV[3];
my %similarity;
my %contigs;
my %used_intron_chains=();
my $ext=33;

#read in genome sequence file
open(FILE,$genome_file);
my $seq="";
my $cn="";
while($line=<FILE>){
  chomp($line);
  if($line=~/^>/){
    my @f=split(/\s+/,$line);
    $contigs{$cn}=$seq if(length($seq)>0);
    $seq="";
    $cn=substr($f[0],1);
  }else{
    $seq.=$line;
  }
}
$contigs{$cn}=$seq if(length($seq)>0);

my @unused=();
open(FILE,$unused_proteins_file);
while($line=<FILE>){
  push(@unused,$line);
  my @F=split(/\t/,$line);
  if($F[2] eq "gene"){
    $similarity{$1}=($4+$3)/2 if($F[8]=~/^ID=(\S+);geneID=(\S+);identity=(\S+);similarity=(\S+)/ );
  }
}

open(FILE,$counts_file);
while($line=<FILE>){
  chomp($line);
  my ($name,$c)=split(/\s+/,$line);
  $pcount{$name}=$c;
}

my $num_complete=0;
open(FILE,$k_file);
while($line=<FILE>){
  $num_complete++ if($line=~/Evidence=complete/);
}

my %mito_contigs=(); 
if(defined($ARGV[4])){  
  open(FILE,$ARGV[4]);  
  while(my $line=<FILE>){
    chomp($line);
    $mito_contigs{$line}=1;
  }
}


#clustered file on STDIN
my $gcode=0;
while($line=<STDIN>){
  chomp($line);
  my @f=split(/\t/,$line);
  $gcode=defined($mito_contigs{$f[0]}) ? 1 : 0;
  if($f[2] eq "gene"){
    my $transcript_id;
    my $gene_id;
    my $startcodon;
    my $stopcodon;
    my $intron_chain;
    if($f[8]=~/^ID=(\S+);exonCount=(\S+);exons=(\S+);CDS=(\S+);CDSphase=(\S+);locus=(\S+)/){
      $transcript_id=$1;
      $intron_chain="$f[0]:$f[6]:$3";
      $gene_id=$6;
    }else{
      next;
    }
    $pcount{$transcript_id}=1 if($gcode>0 && not(defined($pcount{$transcript_id})));
    $intron_chains{$transcript_id}=$intron_chain;
    #print "DEBUG $transcript_id $intron_chain\n";
    next if(not(defined($pcount{$transcript_id})));# if this is not defined then there is an in frame stop
    my $start=$f[3];
    my $end=$f[4];
    my $ori=$f[6];
    if($ori eq "+"){
      for(my $i=0;$i<$ext;$i+=3){
        $startcodon=uc(substr($contigs{$f[0]},$start-1-$i,3));
        last if (valid_start($startcodon));
      }
      for(my $i=0;$i<$ext;$i+=3){
        $stopcodon=uc(substr($contigs{$f[0]},$end-3+$i,3));
        last if(valid_stop($stopcodon,$gcode));
      }
    }else{
      for(my $i=0;$i<$ext;$i+=3){
        my $seq=substr($contigs{$f[0]},$end-3+$i,3);
        $seq=~tr/acgtACGT/tgcaTGCA/;
        $startcodon=uc(reverse($seq));
        last if (valid_start($startcodon));
      }
      for(my $i=0;$i<$ext;$i+=3){
        my $seq=substr($contigs{$f[0]},$start-1-$i,3);
        $seq=~tr/acgtACGT/tgcaTGCA/;
        $stopcodon=uc(reverse($seq));
        last if(valid_stop($stopcodon,$gcode));
      }
    }
    my $codon_score=0;
    $codon_score++ if(valid_start($startcodon));
    $codon_score++ if(valid_stop($stopcodon,$gcode));
    if($codon_score>1){
      my $score=100-(100-$similarity{$transcript_id})/$pcount{$transcript_id};#this scoring boosts proteins that have multiple evidence
      push(@scores,"$score $gene_id $transcript_id $codon_score");
      #print "DEBUG $score $gene_id $transcript_id $ori $codon_score\n";
    }
  }
}

my @scores_sorted=sort { (split(/\s+/, $b))[0] <=> (split(/\s+/, $a))[0] } @scores;
my %h=();
my %hs=();
my %hn=();
my $add_thresh=.9999;
#here we figure out what the additional threshold should be based on the ratio of complete to protein-only
#first we only consider complete proteins
for(my $i=0;$i<=$#scores_sorted;$i++){
  my @F=split(/\s+/,$scores_sorted[$i]);
  #print "DEBUG consider $scores_sorted[$i] $hn{$F[1]} $hs{$F[1]}\n";
  if($hn{$F[1]} < 1 || $F[0]>$hs{$F[1]}*$add_thresh){
    $hn{$F[1]}+=1;#this is the number of proteins per locus
    $hs{$F[1]}=$F[0] if(not(defined($hs{$F[1]})));#this is the highest score per locus
    $min_complete_score=$F[0] if($F[0]<$min_complete_score);
    $h{$F[2]}=1;#we mark the proteins to keep
    #print "DEBUG keep $F[2]\n";
  }
}

#here we adjust the threshold for secondary protein alignments based on the ratio of complete to protein_only
$add_thresh-=scalar(keys %h)/$num_complete/150;
#print "DEBUG $add_thresh $num_complete $min_complete_score\n";
my %h=();
my %hs=();
my %hn=();
for(my $i=0;$i<=$#scores_sorted;$i++){
  my @F=split(/\s+/,$scores_sorted[$i]);
  if($hn{$F[1]} < 1 || $F[0]>$hs{$F[1]}*$add_thresh){
    $hn{$F[1]}+=1;#this is the number of proteins per locus
    $hs{$F[1]}=$F[0] if(not(defined($hs{$F[1]})));#this is the highest score per locus
    $h{$F[2]}=1;#we mark the proteins to keep
  }
}

my $flag=0;
foreach my $l(@unused){
  chomp($l);
  my @f=split(/\t/,$l);
  my $id="";
  if($f[2] eq "gene"){
    $f[2]="transcript";
    $id=$1 if($f[8]=~ /ID=(\S+);geneID=(\S+);identity=(\S+);similarity=(\S+)$/);
    #print "#DEBUG $id $intron_chains{$id} $h{$id} $used_intron_chains{$intron_chains{$id}}\n";
    $flag=((defined($h{$id}) || $id=~/_EXTERNAL$/) && not(defined($used_intron_chains{$intron_chains{$id}}))) ? 1 : 0;
    $used_intron_chains{$intron_chains{$id}}=1 if($flag);
  }
  print join("\t",@f),"\n" if($flag);
}

sub valid_start{
  my $codon=uc($_[0]);
  if(length($codon)==3 && $codon eq "ATG"){
    return(1);
  }else{
    return(0);
  }
}

sub valid_stop{
  my $codon=uc($_[0]);
  my $type=$_[1];
  if($type==0){
    if(length($codon)==3 && ($codon eq "TAG" || $codon eq "TAA" || $codon eq "TGA")){
      return(1);
    }else{
      return(0);
    }
  }elsif($type==1){
    if(length($codon)==3 && ($codon eq "AGA" || $codon eq "AGG" || $codon eq "TAA" || $codon eq "TAG" )){
      return(1);
    }else{
      return(0);
    }  
  }else{
    return(0);
  }
}

