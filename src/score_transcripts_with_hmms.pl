#!/usr/bin/env perl
#This code adds gene and CDS record for every locus in GFF file
my $geneID="";
my @exons=();
my %transcript_gff;
my $protID="";
my $dir="";
my $scf="";
my $seq="";
#these are genetic codes for HMMs
my %code=();
$code{"A"}=0;
$code{"a"}=0;
$code{"C"}=1;
$code{"c"}=1;
$code{"G"}=2;
$code{"g"}=2;
$code{"T"}=3;
$code{"t"}=3;
$code{"N"}=4;
$code{"n"}=4;

#here we load up all transcripts that matched proteins
open(FILE,$ARGV[0]);
while(my $line=<FILE>){
  chomp($line);
  my @gff_fields=split(/\t/,$line);
  my @attributes=split(";",$gff_fields[8]);
  if($gff_fields[2] eq "transcript"){
    if(defined($transcript{$geneID})){
      $transcript_gff{$geneID}=[@exons];
    }
    @exons=();
    $locID=substr($attributes[1],7);#this is the gene_id
    $geneID=substr($attributes[0],3);#this is the transcript_id
    $transcript_junction_score{$geneID}=10000;
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

#we load SNAP HMMs
if(defined($ARGV[2])){
  print "DEBUG Loading SNAP HMMs\n";
  open(FILE,$ARGV[2]);
  $line=<FILE>;
  if($line =~ /^zoeHMM/){#check format
    while($line=<FILE>){
      chomp($line);
      if($line=~/^Donor/){
        $line=<FILE>;
        my $i=0;
        while($line=<FILE>){
          last if($line=~/NN TRM/);
          chomp($line);
          $line=~s/^\s+//;
          my @f=split(/\s+/,$line);
          for(my $j=0;$j<=3;$j++){
            $donor_freq[$i][$j]=$f[$j];
          }
          $donor_freq[$i][4]=0;
          $i++;
        }
        #reset frequencies at the splice sites to not penalize non-canonical splicing that much
        #$donor_freq[3][0]=0;
        #$donor_freq[3][1]=0;
        #$donor_freq[3][3]=0;
        #$donor_freq[4][0]=0;
        #$donor_freq[4][1]=0;
        #$donor_freq[4][2]=0;
      }elsif($line=~/^Acceptor/){
        $line=<FILE>;
        my $i=0;
        while($line=<FILE>){
          last if($line=~/NN TRM/);
          chomp($line);
          $line=~s/^\s+//;
          my @f=split(/\s+/,$line);
          for(my $j=0;$j<=3;$j++){
            $acceptor_freq[$i][$j]=$f[$j];
          }
          $acceptor_freq[$i][4]=0;
          $i++;
        }
        #reset frequencies at the splice sites to not penalize non-canonical splicing that much
        #$acceptor_freq[25][1]=0;
        #$acceptor_freq[25][2]=0;
        #$acceptor_freq[25][3]=0;
        #$acceptor_freq[26][0]=0;
        #$acceptor_freq[26][1]=0;
        #$acceptor_freq[26][3]=0;
      }elsif($line=~/^Start/){
        $line=<FILE>;
        my $i=0;
        while($line=<FILE>){
          last if($line=~/NNN TRM/);
          chomp($line);
          $line=~s/^\s+//;
          my @f=split(/\s+/,$line);
          for(my $j=0;$j<=3;$j++){
            $coding_start_freq[$i][$j]=$f[$j];
          }
          $coding_start_freq[$i][4]=0;
          $i++;
        }
      }
    }
  }
}


#we make the transcript sequences for protein coding transcripts and score the transcripts with HMMs
for my $g(keys %transcript_gff){
  my @gff_fields=();
  for(my $j=1;$j<=$#{$transcript_gff{$g}};$j++){
    @gff_fields=split(/\t/,${$transcript_gff{$g}}[$j]);
    die("Genome sequence $gff_fields[0] needed for transcript $g not found!") if(not(defined($genome_seqs{$gff_fields[0]})));
    if($j>0 && defined($ARGV[2])){
      my @gff_fields_prev=split(/\t/,${$transcript_gff{$g}}[$j-1]);
      if($gff_fields[6] eq "+"){
        $donor_seq=uc(substr($genome_seqs{$gff_fields[0]},$gff_fields_prev[4]-3,9));
        $acceptor_seq=uc(substr($genome_seqs{$gff_fields[0]},$gff_fields[3]-28,30));
        #print "DEBUG $g $gff_fields[6] $donor_seq $acceptor_seq\n";
      }else{
        $donor_seq=uc(substr($genome_seqs{$gff_fields[0]},$gff_fields[3]-7,9));
        $acceptor_seq=uc(substr($genome_seqs{$gff_fields[0]},$gff_fields_prev[4]-3,30));
        $donor_seq=~tr/ACGTNacgtn/TGCANtgcan/;
        $donor_seq=reverse($donor_seq);
        $acceptor_seq=~tr/ACGTNacgtn/TGCANtgcan/;
        $acceptor_seq=reverse($acceptor_seq);
        #print "DEBUG $g $gff_fields[6] $donor_seq $acceptor_seq\n";
      }
      my $junction_score=0;
      for(my $i=0;$i<9;$i++){
        $junction_score+=$donor_freq[$i][$code{substr($donor_seq,$i,1)}];
        #print "DEBUG ",substr($donor_seq,$i,1)," ",$donor_freq[$i][$code{substr($donor_seq,$i,1)}],"\n";
      }
      for(my $i=0;$i<30;$i++){
        $junction_score+=$acceptor_freq[$i][$code{substr($acceptor_seq,$i,1)}];
        #print "DEBUG ",substr($acceptor_seq,$i,1)," ",$acceptor_freq[$i][$code{substr($acceptor_seq,$i,1)}],"\n";
      }
      $transcript_junction_score{$g}=$junction_score if($transcript_junction_score{$g}>$junction_score);
    }
  }
  print "$g $transcript_junction_score{$g} $donor_seq $acceptor_seq\n";
}  

