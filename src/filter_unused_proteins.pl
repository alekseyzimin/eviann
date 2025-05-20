#!/usr/bin/env perl
#this code scores and filters unused protein alignments
my $genome_file=$ARGV[0];
my $unused_proteins_file=$ARGV[1];
my $counts_file=$ARGV[2];
my $k_file=$ARGV[3];
my %similarity;
my %contigs;

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

#combined file on STDIN
while($line=<STDIN>){
  chomp($line);
  my @f=split(/\t/,$line);
  if($f[2] eq "gene"){
    my $transcript_id;
    my $gene_id;
    my $startcodon;
    my $stopcodon;
    if($f[8]=~/ID=(\S+);locus=(\S+)$/){
      $transcript_id=$1;
      $gene_id=$2;
    }else{
      next;
    }
    next if(not(defined($pcount{$transcript_id})));
    my $start=$f[3];
    my $end=$f[4];
    my $ori=$f[6];
    if($ori eq "+"){
      $startcodon=uc(substr($contigs{$f[0]},$start-1,3));
      $stopcodon1=uc(substr($contigs{$f[0]},$end,3));
      $stopcodon2=uc(substr($contigs{$f[0]},$end-3,3));
    }else{
      my $seq=substr($contigs{$f[0]},$end-3,3);
      $seq=~tr/acgtACGT/tgcaTGCA/;
      $startcodon=uc(reverse($seq));
      $seq=substr($contigs{$f[0]},$start-4,3);
      $seq=~tr/acgtACGT/tgcaTGCA/;
      $stopcodon1=uc(reverse($seq));
      $seq=substr($contigs{$f[0]},$start-1,3);
      $seq=~tr/acgtACGT/tgcaTGCA/;
      $stopcodon2=uc(reverse($seq));
    }
    my $codon_score=0;
    $codon_score++ if($startcodon eq "ATG");
    $codon_score++ if($stopcodon1 eq "TAA" || $stopcodon1 eq "TAG" || $stopcodon1 eq "TGA" || $stopcodon2 eq "TAA" || $stopcodon2 eq "TAG" || $stopcodon2 eq "TGA");
    if($codon_score>1){
      #$similarity{$transcript_id}-=20 if($codon_score==1);
      #$similarity{$transcript_id}=1 if($similarity{$transcript_id}<1);
      my $score=100-(100-$similarity{$transcript_id})/$pcount{$transcript_id};#this scoring boosts proteins that have multiple evidence
      push(@scores,"$score $gene_id $transcript_id $codon_score");
      #print "DEBUG $score $gene_id $transcript_id $ori\n";
    }
  }
}

my @scores_sorted=sort { (split(/\s+/, $b))[0] <=> (split(/\s+/, $a))[0] } @scores;
my %h=();
my %hs=();
my %hn=();
my $min_complete_sens=90;
my $add_thresh=.9999;
#here we figure out what the additional threshold should be based on the ratio of complete to protein-only
#first we only consider complete proteins
for(my $i=0;$i<=$#scores_sorted;$i++){
  my @F=split(/\s+/,$scores_sorted[$i]);
  next if($F[3]<2);
  if($hn{$F[1]} < 1 || $F[0]>$hs{$F[1]}*$add_thresh){
    $hn{$F[1]}+=1;#this is the number of proteins per locus
    $hs{$F[1]}=$F[0] if(not(defined($hs{$F[1]})));#this is the highest score per locus
    $h{$F[2]}=1;#we mark the proteins to keep
  }
}

#here we adjust the threshold for secondary protein alignments based on the ratio of complete to protein_only
$add_thresh-=scalar(keys %h)/$num_complete/150;
print "#$add_thresh $num_complete\n";
my %h=();
my %hs=();
my %hn=();
for(my $i=0;$i<=$#scores_sorted;$i++){
  my @F=split(/\s+/,$scores_sorted[$i]);
  next if($F[3]<2);
  if($hn{$F[1]} < 1 || $F[0]>$hs{$F[1]}*$add_thresh){
    $hn{$F[1]}+=1;#this is the number of proteins per locus
    $hs{$F[1]}=$F[0] if(not(defined($hs{$F[1]})));#this is the highest score per locus
    $h{$F[2]}=1;#we mark the proteins to keep
  }
}
#now we use incomplete proteins only for loci that do not ave any completes, 1 per locus
for(my $i=0;$i<=$#scores_sorted;$i++){
  my @F=split(/\s+/,$scores_sorted[$i]);
  next if($F[3]==2);
  if($hn{$F[1]} < 1 && $F[0]>$min_complete_sens ){
    $hn{$F[1]}+=1;#this is the number of proteins per locus
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
    $flag=defined($h{$id}) ? 1 : 0;
  }
  print join("\t",@f),"\n" if($flag);
}

