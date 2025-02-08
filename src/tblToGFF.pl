#!/usr/bin/env perl
my $n=0;
while($l=<STDIN>){
  chomp($l);
  @F=split(/\s+/,$l);
  if($F[0]>=1 && not($F[1]=~/archaea|bacteria|microsporidia/) && $F[12] eq "no"){
    $type="";
    $line=join(" ",@F);
    if($line =~ /spliceosomal RNA/){
      $type="splice_RNA";
    }elsif($line =~ /ribosomal RNA/){
      $type="rRNA";
    }elsif($line =~/Small nucleolar RNA/){
      $type="snoRNA";
    }elsif($line =~/micro RNA|microRNA/){
      $type="miRNA";
    }elsif($line =~/tRNA/){
      $type="tRNA";
    }else{}
    print;
    if(not($type eq "") && $F[16]>60){
      if($F[9]<$F[10]){
        $start=$F[9];
        $end=$F[10];
        $dir="+";
      }else{
        $start=$F[10];
        $end=$F[9];
        $dir="-";
      }
    $id="$type-$n";
    print "$F[3]\tcmsearch\tgene\t$start\t$end\t$F[16]\t$dir\t.\tID=gene-$id;Note=",join("_",@F[26..$#F]),"\n";
    print "$F[3]\tcmsearch\t$type\t$start\t$end\t$F[16]\t$dir\t.\tID=$id;Parent=gene-$id\n";
    print "$F[3]\tcmsearch\texon\t$start\t$end\t$F[16]\t$dir\t.\tID=exon-$id;Parent=$id\n";
    $n++;
    }
  }
}
