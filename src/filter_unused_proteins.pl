#!/usr/bin/env perl
#this code scores and filters unused protein alignments
my $genome_file=$ARGV[0];
my $unused_proteins_file=$ARGV[1];
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

#combined file on STDIN
my $avg=0;
my $count=1;
while($line=<STDIN>){
  chomp($line);
  my @f=split(/\t/,$line);
  if($f[2] eq "transcript"){
    my $transcript_id;
    my $gene_id;
    my $tcount;
    my $startcodon;
    my $stopcodon;
    if($f[8]=~/ID=(\S+);locus=(\S+);count=(\S+)$/){
      $transcript_id=$1;
      $gene_id=$2;
      $tcount=$3;
    }
    $codon_count=0;
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
    $codon_count++ if($startcodon eq "ATG"); 
    $codon_count++ if($stopcodon eq "TAA" || $stopcodon eq "TAG" || $stopcodon eq "TGA"); 
    if(defined($transcript_id) && defined($gene_id) && $codon_count>0){
      my $score=100-(100-$similarity{$transcript_id})/$tcount;
      push(@scores,"$score $gene_id $transcript_id");
      $avg+=$score;
      $count++;
    }
  }
}
my @scores_sorted=sort { (split(/\s+/, $b))[0] <=> (split(/\s+/, $a))[0] } @scores;
my $similarity_threshold=$avg/$count;
for(my $i=0;$i<=$#scores_sorted;$i++){
  my @F=split(/\s+/,$scores_sorted[$i]);
  if(($h{$F[1]} < 1 || $F[0]>$hs{$F[1]}*.99) && $F[0]>=$similarity_threshold){
    $h{$F[1]}+=1;#this is the number of proteins per locus
    $hs{$F[1]}=$F[0] if(not(defined($hs{$F[1]})));#this is the highest score per locus
    $hn{$F[2]}=1;#we mark the proteins to keep
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
  print join("\t",@f),"\n" if(defined($hn{$id}));
}

