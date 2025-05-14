#!/usr/bin/env perl
#this code removes readthrough exons or trims them
#MSTRG_00211120:1:1.860589 5p 41614985 41615486 -

$gtf_file=$ARGV[0];
while($line=<STDIN>){
  chomp($line);
  @F=split(/\s/,$line);
  $h{"$F[0] $F[2] $F[3] $F[4]"}=1;
  if(($F[1] eq "5p" && $F[4] eq "+") || ($F[1] eq "3p" && $F[4] eq "-")){
    $cut_exon{"$F[0] $F[2]"}=$F[3];
  }else{
    $cut_exon{"$F[0] $F[3]"}=$F[2];
  }
}

open(FILE,$gtf_file);
while($line=<FILE>){
  @F=split(/\t/,$line);
  $tid=$1 if($F[8]=~/transcript_id "(\S+)";/);
  if($F[2] eq "exon"){
    unless(defined($h{"$tid $F[3] $F[4] $F[6]"})){
      $F[4]=$cut_exon{"$tid $F[4]"} if(defined($cut_exon{"$tid $F[4]"}));
      $F[3]=$cut_exon{"$tid $F[3]"} if(defined($cut_exon{"$tid $F[3]"}));
      print join("\t",@F);
    }
  }else{
    print $line;
  }
}
