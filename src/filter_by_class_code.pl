#!/usr/bin/env perl
#This code adds gene and CDS record for every locus in GFF file
my $geneID="";
my @exons=();
my %loci;
my %protein_cds;
my %protein;
my %transcript_gff;
my %transcript_cds;
my @gene_records_k;
my %gene_record_k;
my @gene_records_j;
my %gene_record_j;
my @gene_records_u;
my %gene_record_u;
my $protID="";
my $dir="";
my $scf="";
my $seq="";
my $ext_length=300;#this controls how much of an extension for the transcripts we allow for protein alignments that go beyond the transcript boundaries

#this is output of gffcompare -D -o protuniq ../GCF_000001735.4_TAIR10.1_genomic.fna.GCF_000001735.4_TAIR10.1_protein.faa.palign.gff
#here we read in the aligned CDS features
while(my $line=<STDIN>){#we just read in the whole file
  chomp($line);
  my @gff_fields=split(/\t/,$line);
  my @attributes=split(";",$gff_fields[8]);
  if($gff_fields[2] eq "gene"){
    if(not($protID eq "")){
      $protein_start{$protID}=$pstart;
      $protein_end{$protID}=$pend;
    }
    @exons=();
    $protID=substr($attributes[0],3);#this is protein name
    $pstart=$gff_fields[3];
    $pend=$gff_fields[4];
    die("error in line $line, protein ID $protID already exists in $protein{$protID}") if(defined($protein{$protID}));
    $protein{$protID}=$line;
  }
}
if(not($protID eq "")){
  $protein_start{$protID}=$pstart;
  $protein_end{$protID}=$pend;
}

#this is output of gffcompare -r protalign.gtf ../GCF_000001735.4_TAIR10.1_genomic.clean.fna.gtf -o protref
#here we load up all transcripts that matched proteins
open(FILE,$ARGV[0]);
my $printflag=0;
while(my $line=<FILE>){
  chomp($line);
  my @gtf_fields=split(/\t/,$line);
  if($gtf_fields[2] eq "transcript"){
    my $tstart=$gtf_fields[3];
    my $tend=$gtf_fields[4];
    my $geneID=$1 if($gtf_fields[8] =~ /transcript_id \"(\S+)\";/);
    my $protID=$1 if($gtf_fields[8] =~ /cmp_ref \"(\S+)\";/);
    my $class_code=$1 if($gtf_fields[8] =~ /class_code \"(\S+)\";/);

    if($class_code =~ /k|=|u|i|y|x/ || ($class_code =~ /n|j|m|o/  && ($protein_start{$protID} > $tstart-$ext_length  && $protein_end{$protID} < $tend+$ext_length))){
      die("Protein $protID is not defined for protein coding transcript $geneID") if(not(defined($protein{$protID})) && not($class_code eq "u"));
      $printflag=1;
    }else{#likely messed up protein?
      $printflag=0;
    }
  }
  print $line,"\n" if($printflag);
}

