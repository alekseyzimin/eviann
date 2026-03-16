#!/usr/bin/env perl
#we load the genome sequences
my $gff=$ARGV[0];
my %genome_seqs=();
my %masked_genome_seqs=();
my $padding=2000;
while(my $line=<STDIN>){
  chomp($line);
  if($line=~ /^>/){
    if(not($scf eq "")){
      $genome_seqs{$scf}=$seq;
      $seq="";
    }
    $scf=substr((split(/\s+/,$line))[0],1);
  }else{
    $seq.=$line;
  } 
}   
$genome_seqs{$scf}=$seq if(not($scf eq ""));

foreach $g (keys %genome_seqs){
  $masked_genome_seqs{$g}="N"x(length($genome_seqs{$g}));
}

%coord_intervals=();
open(FILE,$gff);
while(my $line=<FILE>){
  chomp($line);
  my @F=split(/\t/,$line);
  my $start=$F[3]-$padding;
  my $end=$F[4]+$padding;
  $start=1 if($start<1);
  $end=length($genome_seqs{$F[0]}) if($end>length($genome_seqs{$F[0]}));
  $coord_intervals{$F[0]}.="$start $end " if($F[2] eq "locus");
}

foreach my $k(keys %coord_intervals){
  my $ctg=$k;
  my @intervals=split(/\s/,$coord_intervals{$k});
  for(my $i=0;$i<$#intervals;$i+=2){
    substr($masked_genome_seqs{$ctg},$intervals[$i]-1,$intervals[$i+1]-$intervals[$i]+1)=substr($genome_seqs{$ctg},$intervals[$i]-1,$intervals[$i+1]-$intervals[$i]+1);
  }
}

foreach $g (keys %masked_genome_seqs){
  print ">$g\n$masked_genome_seqs{$g}\n";
}

