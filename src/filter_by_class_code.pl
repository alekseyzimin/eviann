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
my $ext_length=501;

#this is output of gffcompare -D -o protuniq ../GCF_000001735.4_TAIR10.1_genomic.fna.GCF_000001735.4_TAIR10.1_protein.faa.palign.gff
#here we read in the aligned CDS features
while(my $line=<STDIN>){#we just read in the whole file
  chomp($line);
  my @gff_fields=split(/\t/,$line);
  my @attributes=split(";",$gff_fields[8]);
  if($gff_fields[2] eq "gene"){
    if(not($protID eq "")){
      $protein_cds{$protID}=[@exons];
      $protein_start{$protID}=$pstart;
      $protein_end{$protID}=$pend;
    }
    @exons=();
    $protID=substr($attributes[0],3);#this is protein name
    $pstart=$gff_fields[3];
    $pend=$gff_fields[4];
    $pori=$gff_fields[6];
    die("error in line $line, protein ID $protID already exists in $protein{$protID}") if(defined($protein{$protID}));
    $protein{$protID}=$line;
  }elsif($gff_fields[2] eq "CDS"){
    push(@exons,$line);
  }
}
if(not($protID eq "")){
  $protein_cds{$protID}=[@exons];
  $protein_start{$protID}=$pstart;
  $protein_end{$protID}=$pend;
}

#this is output of gffcompare -r protalign.gtf ../GCF_000001735.4_TAIR10.1_genomic.clean.fna.gtf -o protref
#here we load up all transcripts that matched proteins
open(FILE,$ARGV[0]);
my $printflag=0;
while(my $line=<FILE>){
  chomp($line);
  my @gff_fields=split(/\t/,$line);
  my @attributes=split(";",$gff_fields[8]);
  if($gff_fields[2] eq "transcript"){
    $locID=substr($attributes[1],7);#this is the gene_id
    $geneID=substr($attributes[0],3);#this is the transcript_id
    $tstart=$gff_fields[3];
    $tend=$gff_fields[4];
    $tori=$gff_fields[6];
    my $class_code="";
    my $protID="";
    foreach my $attr(@attributes){
      $class_code=substr($attr,11,1) if($attr =~ /^class_code=/);
      $protID=substr($attr,8) if($attr =~ /^cmp_ref=/);
    }
    if($class_code eq "k" || $class_code eq "=" || $class_code eq "u" || (($class_code eq "m" || $class_code eq "j" ||$class_code eq "n" || $class_code eq "c") && ($protein_start{$protID} > $tstart-$ext_length  && $protein_end{$protID} < $tend+$ext_length))){#equal intron chain or contains protein
      die("Protein $protID is not defined for protein coding transcript $geneID") if(not(defined($protein{$protID})));
      $printflag=1;
    }else{#likely messed up protein?
      $printflag=0;
    }
  }
  print $line,"\n" if($printflag);
}

