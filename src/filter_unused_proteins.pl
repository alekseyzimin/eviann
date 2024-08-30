#!/usr/bin/env perl
#this code scores and filters unused protein alignments
my $genome_file=$ARGV[0];
my $unused_proteins_file=$ARGV[1];
my $counts_file=$ARGV[2];
my $liftover=$ARGV[3];
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

open(FILE,$unused_proteins_file);
while($line=<FILE>){
  my @F=split(/\t/,$line);
  if($F[2] eq "gene"){
    $similarity{$1}=$4+$3/100-1 if($F[8]=~/^ID=(\S+);geneID=(\S+);identity=(\S+);similarity=(\S+)/ );
  }
}

open(FILE,$counts_file);
while($line=<FILE>){
  chomp($line);
  my ($name,$c)=split(/\s+/,$line);
  $pcount{$name}=$c;
}

#combined file on STDIN
while($line=<STDIN>){
  chomp($line);
  my @f=split(/\t/,$line);
  if($f[2] eq "transcript"){
    my $transcript_id;
    my $gene_id;
    my $tcount;
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
      $stopcodon=uc(substr($contigs{$f[0]},$end,3));
    }else{
      my $seq=substr($contigs{$f[0]},$end-3,3);
      $seq=~tr/acgtACGT/tgcaTGCA/;
      $startcodon=uc(reverse($seq));
      $seq=substr($contigs{$f[0]},$start-4,3);
      $seq=~tr/acgtACGT/tgcaTGCA/;
      $stopcodon=uc(reverse($seq));
    }
    $codon_count=0;
    $codon_count++ if($startcodon eq "ATG"); 
    $codon_count++ if($stopcodon eq "TAA" || $stopcodon eq "TAG" || $stopcodon eq "TGA"); 
    if($codon_count>1){
      my $score=100-(100-$similarity{$transcript_id})/$pcount{$transcript_id};#this scoring boosts proteins that have multiple evidence
      push(@scores,"$score $gene_id $transcript_id");
      #print "DEBUG $score $gene_id $transcript_id $ori\n";
    }
  }
}
my @scores_sorted=sort { (split(/\s+/, $b))[0] <=> (split(/\s+/, $a))[0] } @scores;
my $similarity_threshold=0;
if(not(defined($liftover)) || $liftover<1){
  #my @f=split(/\s+/,$scores_sorted[int($#scores_sorted*.95)]);
  my @uniq_scores=();
  my %hs=();
  my %hn=();
  for(my $i=0;$i<=$#scores_sorted;$i++){
    my @F=split(/\s+/,$scores_sorted[$i]);
    if($hn{$F[1]} < 1 || $F[0]>$hs{$F[1]}*.99){
      $hn{$F[1]}+=1;#this is the number of proteins per locus
      $hs{$F[1]}=$F[0] if(not(defined($hs{$F[1]})));#this is the highest score per locus
      push(@uniq_scores,$F[0]);
    }
  }
  $similarity_threshold=$uniq_scores[$#uniq_scores];
}
print "#similarity threshold $similarity_threshold\n";
my %h=();
my %hs=();
my %hn=();
for(my $i=0;$i<=$#scores_sorted;$i++){
  my @F=split(/\s+/,$scores_sorted[$i]);
  if(($hn{$F[1]} < 1 || $F[0]>$hs{$F[1]}*.99) && $F[0]>=$similarity_threshold){
    $hn{$F[1]}+=1;#this is the number of proteins per locus
    $hs{$F[1]}=$F[0] if(not(defined($hs{$F[1]})));#this is the highest score per locus
    $h{$F[2]}=1;#we mark the proteins to keep
  }
}

open(FILE,$unused_proteins_file);
while($line=<FILE>){
  chomp($line);
  my @f=split(/\t/,$line);
  my $id="";
  if($f[2] eq "gene"){
    $f[2]="transcript";
    $id=$1 if($f[8]=~ /ID=(\S+);geneID=(\S+);identity=(\S+);similarity=(\S+)$/);
  }elsif($f[2] eq "exon" || $f[2] eq "cds"){
    $id=$1 if($f[8]=~ /Parent=(\S+)$/);
  }
  print join("\t",@f),"\n" if(defined($h{$id}));
}

