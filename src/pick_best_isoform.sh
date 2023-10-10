#!/bin/bash
PREFIX=$1;
ufasta one  $PREFIX.functional_note.proteins.fasta | \
awk '{if($1 ~/^>/){pn=substr($1,2);if($2 ~ /Similar/){s=10000}else{s=0}}else{print pn" "s+length($1)}}' |\
sort -nrk2,2 -S10% |\
perl -ane '{
  @f=split(/-/,$F[0]);
  if(not(defined($h{$f[0]}))){
    $h{$f[0]}=$F[0];
    $hh{$F[0]}=1;
  }
}END{
  open(FILE,"'$PREFIX'.functional_note.pseudo_label.gff");
  while($line=<FILE>){
    chomp($line);
    @f=split(/\t/,$line);
    @ff=split(/;/,$f[8]);
    ($junk,$id)=split(/=/,$ff[0]);
    ($junk,$parent)=split(/=/,$ff[1]);
    print $line,"\n" if(defined($hh{$id}) || defined($hh{$parent}));
  }
}' > $PREFIX.bestisoform.gff
