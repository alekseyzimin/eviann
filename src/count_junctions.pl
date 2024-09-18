#!/usr/bin/env perl
#
## EviAnn automated annotation
#NC_004353.4     EviAnn  gene    929     5036    .       +       .       ID=XLOC_000001;geneID=XLOC_000001;type=protein_coding
#NC_004353.4     EviAnn  mRNA    930     5036    .       +       .       ID=XLOC_000001-mRNA-1;Parent=XLOC_000001;EvidenceProteinID=XP_002133780.1;EvidenceTranscriptID=MSTRG_00000479:2:3.681515;StartCodon=atg;StopCodon=taa;Class==;Evidence=complete;
#NC_004353.4     EviAnn  exon    930     1079    .       +       .       ID=XLOC_000001-mRNA-1:exon:1;Parent=XLOC_000001-mRNA-1
#NC_004353.4     EviAnn  exon    1144    1410    .       +       .       ID=XLOC_000001-mRNA-1:exon:2;Parent=XLOC_000001-mRNA-1
#NC_004353.4     EviAnn  exon    1768    3851    .       +       .       ID=XLOC_000001-mRNA-1:exon:3;Parent=XLOC_000001-mRNA-1
#NC_004353.4     EviAnn  exon    4353    4527    .       +       .       ID=XLOC_000001-mRNA-1:exon:4;Parent=XLOC_000001-mRNA-1
#NC_004353.4     EviAnn  exon    4592    4824    .       +       .       ID=XLOC_000001-mRNA-1:exon:5;Parent=XLOC_000001-mRNA-1
#NC_004353.4     EviAnn  exon    4875    5036    .       +       .       ID=XLOC_000001-mRNA-1:exon:6;Parent=XLOC_000001-mRNA-1
#NC_004353.4     EviAnn  cds     930     1079    .       +       .       ID=XLOC_000001-mRNA-1:cds:1;Parent=XLOC_000001-mRNA-1
#NC_004353.4     EviAnn  cds     1144    1410    .       +       .       ID=XLOC_000001-mRNA-1:cds:2;Parent=XLOC_000001-mRNA-1
#NC_004353.4     EviAnn  cds     1768    3851    .       +       .       ID=XLOC_000001-mRNA-1:cds:3;Parent=XLOC_000001-mRNA-1
#NC_004353.4     EviAnn  cds     4353    4527    .       +       .       ID=XLOC_000001-mRNA-1:cds:4;Parent=XLOC_000001-mRNA-1
#NC_004353.4     EviAnn  cds     4592    4824    .       +       .       ID=XLOC_000001-mRNA-1:cds:5;Parent=XLOC_000001-mRNA-1
#NC_004353.4     EviAnn  cds     4875    4995    .       +       .       ID=XLOC_000001-mRNA-1:cds:6;Parent=XLOC_000001-mRNA-1


my @cds=();
my %jcounts=();
while($line=<STDIN>){
  chomp($line);
  push(@gff,$line);
  my @F=split(/\t/,$line);
  if($F[2] eq "gene"){
    if(scalar(@cds)>1){
      for($i=1;$i<=$#cds;$i++){
        my ($p1,$p2)=split(/\s/,$cds[$i-1]);
        my ($c1,$c2)=split(/\s/,$cds[$i]);
        $jcounts{"$seq $p2 $c1"}++;
      }
    }
    @cds=();
    @f=split(/;/,$F[8]);
    $seq=$F[0];
  }elsif($F[2] eq "CDS"){
    push(@cds,"$F[3] $F[4]");
  }
}

foreach $line (@gff){
  my @F=split(/\t/,$line);
  if($F[2] eq "gene"){
    $min_count=1000000000;
    if(scalar(@cds)>1){
      for($i=1;$i<=$#cds;$i++){
        my ($p1,$p2)=split(/\s/,$cds[$i-1]);
        my ($c1,$c2)=split(/\s/,$cds[$i]);
        $min_count=$jcounts{"$seq $p2 $c1"} if($jcounts{"$seq $p2 $c1"}<$min_count);
      } 
    }
    $min_gene_count{$id}=$min_count;
    @cds=();
    @f=split(/;/,$F[8]);
    $id=substr($f[0],3);
    $seq=$F[0];
  }elsif($F[2] eq "CDS"){
    push(@cds,"$F[3] $F[4]");
  }
}

foreach $line (@gff){
  my @F=split(/\t/,$line);
  if($F[2] eq "gene"){
    $F[2]="transcript";
    @f=split(/;/,$F[8]);
    $id=substr($f[0],3);
    $flag=($min_gene_count{$id}>$ARGV[0]) ? 1 : 0;
  }
  print join("\t",@F),";min_j_count=$min_gene_count{$id}\n" if($flag);
}

