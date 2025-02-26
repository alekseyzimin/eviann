#!/usr/bin/env perl
$loci_file=$ARGV[0];
open(FILE,$loci_file);
while($line=<FILE>){
  chomp($line);
  my ($l,$transcripts)=split(/\s+/,$line);
  @t=split(/,/,$transcripts);
  foreach $tr(@t){
    $reassign_locus{$tr}=$l;
    $locus_beg{$l}=200000000000;
    $locus_end{$l}=-1;
  }
}

while($line=<STDIN>){
  print $line if ($line =~/^#/);
  @f=split(/\t/,$line);
  if($f[2] eq "mRNA"){
    if($f[8] =~ /^ID=(\S+);Parent=(\S+);EvidenceProteinID/){
      #this is protein coding transcript 
      $tr=$1;
      $lid=$2;
      if(defined($reassign_locus{$tr}))
        $f[8]=~s/Parent=$lid;EvidenceProteinID/Parent=$reassign_locus{$tr};EvidenceProteinID/;
        $lid=$reassign_locus{$tr};
      }
      $locus{$lid}.=join("\t",$line);
      $locus_beg{$lid}=$f[3] if($locus_beg{$lid}>$f[3]);
      $locus_end{$lid}=$f[4] if($locus_end{$lid}<$f[4]);
      $locis_chr{$lid}=$f[0] unless(defined($locis_chr{$lid}));
      $locis_ori{$lid}=$f[6] unless(defined($locis_ori{$lid}));
      $locis_src{$lid}=$f[1] unless(defined($locis_src{$lid}));
    }else{
      #non-coding transcript
      print $line;
    }
  }elsif($f[8]=~/_lncRNA/){
    print $line;
  }else{
    $locus{$reassign_locus{$tr}}.=$line;
  }
}

foreach $lid(keys %locus){
  print "$locus_chr{$lid}\t$locus_src{$lid}\tgene\t$locus_beg{$lid}\t$locus_end{$lid}\t.\t$locus_ori{$lid}\t.\tID=$lid;geneID=$lid;gene_biotype=protein_coding\n$locus{$lid}";
}

