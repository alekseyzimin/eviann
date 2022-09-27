#!/usr/bin/env perl
#exonerate makes mistakes inserting introns of size 11,14,...,59.  This code removes introns of these sizes if they do not match the splice sites by stringtie
my $prev="";
my $prevc=0;
open(FILE,$ARGV[0]);
while($line=<FILE>){
  chomp($line);
  @F=split(/\t/,$line);
  if($F[2] eq "exon" && $prev eq "exon"){
    $introndb{"$F[0] ".($prevc+1)." ".($F[3]-1)}=1;
  }
  $prev=$F[2];
  $prevc=$F[4];
}
my $prev_intron_good=1;
my $first_cds=1;
while($line=<STDIN>){
  chomp($line);
  @f=split(/\t/,$line);
  if($f[2] eq "gene"){
    print join("\t",@{$cds}),"\n" unless($first_cds);
    print join("\t",@{$cds}[0..1]),"\texon\t",join("\t",@{$cds}[3..$#{$cds}]),"\n" unless($first_cds);
    print join("\t",@f),"\n";
    $prev_intron_good=0;
    $first_cds=1;
  }elsif($f[2] eq "cds"){
    if($first_cds){
      $cds=[@f];
      $first_cds=0;
      $prev_intron_good=0;
    }else{
      if($prev_intron_good){
        print join("\t",@{$cds}),"\n";
        print join("\t",@{$cds}[0..1]),"\texon\t",join("\t",@{$cds}[3..$#{$cds}]),"\n";
        $cds=[@f];
      }else{
        if($f[6] eq "+"){
          $cds=[(@{$cds}[0..3],@f[4..$#f])];
        }else{
          $cds=[(@f[0..3],@{$cds}[4..$#{$cds}])]
        }
      }
    }
  }elsif($f[2] eq "intron"){
    $intron_size=$f[4]-$f[3]; 
    $intron_suspect=0; 
    for($i=11;$i<60;$i+=3){
      $intron_suspect=1 if($intron_size==$i);
    }
    $prev_intron_good=1;
    $prev_intron_good=0  if(not(defined($introndb{"$f[0] $f[3] $f[4]"})) && $intron_suspect); 
  }
}
print join("\t",@{$cds}),"\n" unless($first_cds);
print join("\t",@{$cds}[0..1]),"\texon\t",join("\t",@{$cds}[3..$#{$cds}]),"\n" unless($first_cds);
