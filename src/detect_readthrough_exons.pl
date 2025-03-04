#!/usr/bin/env perl
$protein_gff=$ARGV[0];

open(FILE,$protein_gff);
while($line=<FILE>){
  chomp($line);
  @F=split(/\t/,$line);
  if($F[2] eq "gene"){
    @f=split(/;/,$F[8]);
    $pbeg{substr($f[0],3)}=$F[3];
    $pend{substr($f[0],3)}=$F[4];
  }
}

my $protid="";
my $transcript="";
while($line=<STDIN>){
  chomp($line);
  @F=split(/\t/,$line);
  if($F[2] eq "transcript"){
    $flag=0;
    if($F[8]=~/^transcript_id "(\S+)"; gene_id (.+); cmp_ref "(\S+)"; class_code "(k|=|j|c|m|n|o|i|y)"/){
      $protid=$3;
      $tid=$1;
      @tr=split(/\./,$1);
      $transcript=join(".",@tr[0..($#tr-1)])." ".$tr[-1];
      $flag=1;
      #print "DEBUG $transcript $F[3] $F[4] $pbeg{$protid} $pend{$protid} $flag\n";
    }
  }elsif($F[2] eq "exon" && $flag){
    #print "DEBUG exon $transcript $F[3] $F[4] $F[6] $pbeg{$protid} $pend{$protid}\n";
    print "$transcript $F[3] $F[4] $F[6] $pbeg{$protid} $pend{$protid}\n" if((($tid =~ /5p$/ && $pend{$protid}>=$F[3])||($tid =~ /3p$/ && $pbeg{$protid}<=$F[4])) && $F[6] eq "+");
    print "$transcript $F[3] $F[4] $F[6] $pbeg{$protid} $pend{$protid}\n" if((($tid =~ /5p$/ && $pbeg{$protid}<=$F[4])||($tid =~ /3p$/ && $pend{$protid}>=$F[3])) && $F[6] eq "-");
  }
}
