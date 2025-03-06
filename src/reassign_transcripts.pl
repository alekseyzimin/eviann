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
    $locus_count{$l}=1;
  }
}

while($line=<STDIN>){
  @f=split(/\t/,$line);
  if($line =~/^#/ || $f[8]=~/_lncRNA/){
    print $line;
    next;
  }
  if($f[2] eq "mRNA"){
    if($f[8] =~ /^ID=(\S+);EvidenceProteinID/){
      #this is protein coding transcript 
      $tr=$1;
      ($lid)=split(/-/,$tr,1);
      if(defined($reassign_locus{$tr})){
        my @ff=split(/;/,$f[8]); 
        my @fff=split(/-/,$ff[0]);
        $lid=$reassign_locus{$tr};
        splice(@ff,1,0,"Parent=$lid");
        $ff[-2]="geneID=$lid";
        $tr=$lid."-mRNA-".$locus_count{$lid};
        $ff[0]="ID=$tr";
        $f[8]=join(";",@ff);
        $locus_count{$lid}++;
      }
      $locus{$lid}.=join("\t",@f);
      $locus_beg{$lid}=$f[3] if($locus_beg{$lid}>$f[3]);
      $locus_end{$lid}=$f[4] if($locus_end{$lid}<$f[4]);
      $locus_chr{$lid}=$f[0] unless(defined($locus_chr{$lid}));
      $locus_ori{$lid}=$f[6] unless(defined($locus_ori{$lid}));
      $locus_src{$lid}=$f[1] unless(defined($locus_src{$lid})); 
    }
  }elsif($f[2] eq "gene"){
  }else{
    $f[8]="Parent=$tr\n";
    $locus{$lid}.=join("\t",@f);
  }
}

foreach $lid(keys %locus){
  print "$locus_chr{$lid}\t$locus_src{$lid}\tgene\t$locus_beg{$lid}\t$locus_end{$lid}\t.\t$locus_ori{$lid}\t.\tID=$lid;geneID=$lid;gene_biotype=protein_coding\n$locus{$lid}";
}

