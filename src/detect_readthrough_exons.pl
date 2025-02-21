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

while($line=<STDIN>){
  chomp($line);
  @F=split(/\t/,$line);
  if($F[2] eq "transcript"){
    $flag=0;
    if($F[8]=~/^transcript_id "(\S+)"; gene_id (.+); cmp_ref "(\S+)"; class_code "(k|=|j|c|m|n|o)"/){
      $protid=$3;
      $tid=$1;
      @tr=split(/\./,$1);
      $transcript=join(".",@tr[0..($#tr-1)])." ".$tr[-1];
      if($F[6] eq "+"){
        $flag=1 if($tid =~ /5p$/ && $pend{$protid}<$F[4]);
        $flag=1 if($tid =~ /3p$/ && $pbeg{$protid}>$F[3]);
      }else{
        $flag=1 if($tid =~ /5p$/ && $pbeg{$protid}>$F[3]);
        $flag=1 if($tid =~ /3p$/ && $pend{$protid}<$F[4]);
      }
    }
  }elsif($F[2] eq "exon" && $flag){
    print "$transcript $F[3] $F[4] $F[6]\n";
  }
}
