#!/usr/bin/env perl 
#3.1_contig_10_(contig_13_RagTag)        EviAnn  exon    6460    6932    .       -       .       ID=XLOC_000151_lncRNA-mRNA-1:exon:5;Parent=XLOC_000151_lncRNA-mRNA-1
#XLOC_000151_lncRNA-mRNA-1       transdecoder    CDS     1       1989
##this code integrates the CDS features detected in transcripts into the gff file
#CDS gff is the first arg, input GFF on STDIN
my $trans_gff=$ARGV[0];
my %cds_start=();
my %cds_end=();
my $transc_len=();
my $transc_l=0;

open(FILE, $trans_gff);
while($line=<FILE>){
  chomp($line);
  my @f=split(/\t/,$line);
  if($f[2] eq "mRNA"){
    $transc_l=$f[4];
  }elsif($f[2] eq "CDS"){
    if(defined($cds_start{$f[0]})){
    #use the longer one
      if($cds_end{$f[0]}-$cds_start{$f[0]} < $f[4]-$f[3]){
        $cds_start{$f[0]}=$f[3];
        $cds_end{$f[0]}=$f[4];
        $transc_len{$f[0]}=$transc_l;
      }
    }else{
      $cds_start{$f[0]}=$f[3];
      $cds_end{$f[0]}=$f[4];
      $transc_len{$f[0]}=$transc_l;
    }
  #print "DEBUG $f[0] $cds_start{$f[0]} $cds_end{$f[0]} $transc_len{$f[0]}\n";
  }
}

my $remember_exons=0;
my @exons=();
close(FILE);
while($line=<STDIN>){
  chomp($line);  
  my @f=split(/\t/,$line);
  if(not($f[2] eq "exon") && $remember_exons){
    #here we output the CDS records
    my @ff=split(/\t/,$exons[0]);
    my @fff=split(/(;|:)/,$ff[8]);
    my $rna_id=substr($fff[0],3);
    if($ff[6] eq "-"){
      my $tstart=$cds_start{$rna_id};
      my $tend=$cds_end{$rna_id};
      $cds_start{$rna_id}=$transc_len{$rna_id}-$tend+1;
      $cds_end{$rna_id}=$transc_len{$rna_id}-$tstart+1;
    }
    #determining cds start and end coordinates in the genome
    my $seen_exon_length=0;
    my $cds_genome_start=-1;
    my $cds_genome_end=-1;
    for my $e(@exons){
      my @fe=split(/\t/,$e);
      $seen_exon_length+=$fe[4]-$fe[3]+1;
      $cds_genome_start=$fe[4]-($seen_exon_length-$cds_start{$rna_id}) if($cds_genome_start == -1 && $cds_genome_end == -1 && $seen_exon_length>$cds_start{$rna_id});
      $cds_genome_end=$fe[4]-($seen_exon_length-$cds_end{$rna_id}) if($cds_genome_end == -1 && $cds_genome_start>-1 && $seen_exon_length>=$cds_end{$rna_id});
    }
    #print "DEBUG outputting CDS for $rna_id $cds_start{$rna_id} $cds_end{$rna_id} $transc_len{$rna_id} $cds_genome_start $cds_genome_end\n";
    for(my $i=0;$i<=$#exons;$i++){
      #print "DEBUG exon $exons[$i]\n";
      my @fe=split(/\t/,$exons[$i]);
      $fe[8]=~s/_lncRNA//g;
      if($cds_genome_start>$fe[4]){
        #this means that the CDS has not started yet, the whole exon is 5' UTR
        print join("\t",@fe[0..1]),"\tfive_prime_UTR\t",join("\t",@fe[3..$#fe]),"\n";
      }elsif($cds_genome_start>=$fe[3] && $cds_genome_start<=$fe[4] && $cds_genome_end>=$fe[3] && $cds_genome_end<=$fe[4]){
        #this means that the CDS starts and ends within $exons[$i]
        print join("\t",@fe[0..1]),"\tfive_prime_UTR\t$fe[3]\t",$cds_genome_start-1,"\t",join("\t",@fe[5..$#fe]),"\n" if($cds_genome_start>$fe[3]);
        print join("\t",@fe[0..1]),"\tcds\t$cds_genome_start\t$cds_genome_end\t",join("\t",@fe[5..$#fe]),"\n";
        print join("\t",@fe[0..1]),"\tthree_prime_UTR\t",$cds_genome_end+1,"\t",join("\t",@fe[4..$#fe]),"\n" if($cds_genome_end<$fe[4]);
      }elsif($cds_genome_start>=$fe[3] && $cds_genome_start<=$fe[4] && $cds_genome_end>$fe[4]){
        #this means that the CDS started within $exons[$i], but ends beyond
        print join("\t",@fe[0..1]),"\tfive_prime_UTR\t$fe[3]\t",$cds_genome_start-1,"\t",join("\t",@fe[5..$#fe]),"\n" if($cds_genome_start>$fe[3]);
        print join("\t",@fe[0..1]),"\tcds\t$cds_genome_start\t",join("\t",@fe[4..$#fe]),"\n";
      }elsif($cds_genome_end>=$fe[3] && $cds_genome_end<=$fe[4] && $cds_genome_start<$fe[3]){
        #this means that the CDS started before but ends inside the current exon
        print join("\t",@fe[0..1]),"\tcds\t$fe[3]\t$cds_genome_end\t",join("\t",@fe[5..$#fe]),"\n";
        print join("\t",@fe[0..1]),"\tthree_prime_UTR\t",$cds_genome_end+1,"\t",join("\t",@fe[4..$#fe]),"\n" if($cds_genome_end<$fe[4]);
      }elsif($cds_genome_end<$fe[3]){
        #this means that the CDS ended and the whole exon is 3' UTR
        print join("\t",@fe[0..1]),"\tthree_prime_UTR\t",join("\t",@fe[3..$#fe]),"\n";
      }else{
        #this means that the whole exon contains CDS
        print join("\t",@fe[0..1]),"\tcds\t",join("\t",@fe[3..$#fe]),"\n";
      }
    }
    $remember_exons=0;
    @exons=();
  }
  if($f[2] eq "mRNA"){
    my @ff=split(/;/,$f[8]);
    my $rna_id=substr($ff[0],3);
    if(defined($cds_start{$rna_id})){
      #read in exons and determine UTR and CDS coordinates
      $remember_exons=1;
      #remove lncRNA from the name
      $f[8]=~s/_lncRNA//g;
      $line=join("\t",@f);
    }
  }elsif($f[2] eq "exon"){
    if( $remember_exons ){
      push(@exons,$line);
      $f[8]=~s/_lncRNA//g;
      $line=join("\t",@f);
    }
  }elsif($f[2] eq "gene"){
    my @ff=split(/;/,$f[8]);
    my @fff=split(/_/,substr($ff[0],3));
    if(defined($cds_start{join("_",@fff[0..2])."-mRNA-1"})||defined($cds_start{join("_",@fff[0..2])."-mRNA-2"})||defined($cds_start{join("_",@fff[0..2])."-mRNA-3"})){
      $f[8]=~s/_lncRNA//g;
      $f[8]=~s/type=lncRNA/type=protein_coding/;
      $line=join("\t",@f);
    }
  }
  print $line,"\n";
}

