#!/usr/bin/env perl
#This code extracts all donor/acceptor sequences
#here we load up all transcripts
open(FILE,$ARGV[0]);
while(my $line=<FILE>){
  chomp($line);
  my @gff_fields=split(/\t/,$line);
  my @attributes=split(";",$gff_fields[8]);
  if($gff_fields[2] eq "transcript" || $gff_fields[2] eq "mRNA" ){
    if(defined($transcript{$geneID})){
      $transcript_gff{$geneID}=[@exons];
    }
    @exons=();
    $geneID=substr($attributes[0],3);#this is the transcript_id 
    $transcript{$geneID}=$line;
  }elsif($gff_fields[2] eq "exon"){
    push(@exons,$line) if(defined($transcript{$geneID}));
  }
}
if(defined($transcript{$geneID})){
  $transcript_gff{$geneID}=[@exons];
}
@exons=();

#we load the genome sequences
open(FILE,$ARGV[1]);
print "DEBUG Loading genome sequence\n";
while(my $line=<FILE>){
  chomp($line);
  if($line=~ /^>/){
    if(not($scf eq "")){
      $genome_seqs{$scf}=$seq;
      $seq="";
    }
    my @f=split(/\s+/,$line);
    $scf=substr($f[0],1);
  }else{
    $seq.=$line;
  } 
}   
$genome_seqs{$scf}=$seq if(not($scf eq ""));

#we make the transcript sequences for protein coding transcripts and score the transcripts with HMMs
for my $g(keys %transcript_gff){
  my @gff_fields=();
  for(my $j=1;$j<=$#{$transcript_gff{$g}};$j++){
    @gff_fields=split(/\t/,${$transcript_gff{$g}}[$j]);
    die("Genome sequence $gff_fields[0] needed for transcript $g not found!") if(not(defined($genome_seqs{$gff_fields[0]})));
    my @gff_fields_prev=split(/\t/,${$transcript_gff{$g}}[$j-1]);
    if($gff_fields[6] eq "+"){
      $donor_seq=uc(substr($genome_seqs{$gff_fields[0]},$gff_fields_prev[4],2));
      $acceptor_seq=uc(substr($genome_seqs{$gff_fields[0]},$gff_fields[3]-3,2));
      my $dtype=($donor_seq eq "GT") ? "CANONICAL" : "NONCANONICAL";
      my $atype=($acceptor_seq eq "AG") ? "CANONICAL" : "NONCANONICAL";
      print "DONOR $gff_fields[0] $gff_fields_prev[4] $gff_fields_prev[6] $donor_seq $dtype\n";
      print "ACCEPTOR $gff_fields[0] $gff_fields[3] $gff_fields[6] $acceptor_seq $atype\n";
    }else{
      $donor_seq=uc(substr($genome_seqs{$gff_fields[0]},$gff_fields[3]-3,2));
      $acceptor_seq=uc(substr($genome_seqs{$gff_fields[0]},$gff_fields_prev[4],2));
      $donor_seq=~tr/ACGTNacgtn/TGCANtgcan/;
      $donor_seq=reverse($donor_seq);
      $acceptor_seq=~tr/ACGTNacgtn/TGCANtgcan/;
      $acceptor_seq=reverse($acceptor_seq);
      my $dtype=($donor_seq eq "GT") ? "CANONICAL" : "NONCANONICAL";
      my $atype=($acceptor_seq eq "AG") ? "CANONICAL" : "NONCANONICAL";
      print "ACCEPTOR $gff_fields[0] $gff_fields_prev[4] $gff_fields_prev[6] $acceptor_seq $atype\n";
      print "DONOR $gff_fields[0] $gff_fields[3] $gff_fields[6] $donor_seq $dtype\n";
    }
  }
}  

