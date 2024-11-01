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
    $transcript_acceptor_score{$geneID}=10000;
    $transcript_donor_score{$geneID}=10000;
    $transcript_hmm_donor_score{$geneID}=10000;
    $transcript_hmm_acceptor_score{$geneID}=10000;
    $transcript_hmm_score{$geneID}=10000;
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
      }elsif($line=~/^SDonor/){
        while($line=<FILE>){
          chomp($line);
          my @f=split(/\s+/,$line);
          last if($line eq "NNNNNN");
          $sdonor{$f[0]}=$f[1];
        }
      }elsif($line=~/^SAcceptor/){
        while($line=<FILE>){
          chomp($line);
          my @f=split(/\s+/,$line);
          last if($line eq "NNNNNN");
          $sacceptor{$f[0]}=$f[1];
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
      my $donor_score=0;
      my $acceptor_score=0;
      for(my $i=0;$i<9;$i++){
        $donor_score+=$donor_freq[$i][$code{substr($donor_seq,$i,1)}];
        #print "DEBUG ",substr($donor_seq,$i,1)," ",$donor_freq[$i][$code{substr($donor_seq,$i,1)}],"\n";
      }
      for(my $i=0;$i<30;$i++){
        $acceptor_score+=$acceptor_freq[$i][$code{substr($acceptor_seq,$i,1)}];
        #print "DEBUG ",substr($acceptor_seq,$i,1)," ",$acceptor_freq[$i][$code{substr($acceptor_seq,$i,1)}],"\n";
      }
      my $junction_score=$donor_score+$acceptor_score;
      my $hmm_donor_score=defined($sdonor{substr($donor_seq,2,1).substr($donor_seq,5,4)})?$sdonor{substr($donor_seq,2,1).substr($donor_seq,5,4)}:-10000;
      my $hmm_acceptor_score=defined($sacceptor{substr($acceptor_seq,21,4).substr($acceptor_seq,27,1)})?$sacceptor{substr($acceptor_seq,21,4).substr($acceptor_seq,27,1)}:-10000;
      my $hmm_score=$hmm_donor_score+$hmm_acceptor_score;
      #print "DEBUG $hmm_score ",$sdonor{substr($donor_seq,0,3).substr($donor_seq,5,4)}," ",$sacceptor{substr($acceptor_seq,21,4).substr($acceptor_seq,27,3)},"\n";
      $transcript_junction_score{$g}=$junction_score if($transcript_junction_score{$g}>$junction_score);
      $transcript_hmm_donor_score{$g}=$hmm_donor_score if($transcript_hmm_donor_score{$g}>$hmm_donor_score);
      $transcript_hmm_score{$g}=$hmm_score if($transcript_hmm_score{$g}>$hmm_score);
      $transcript_hmm_acceptor_score{$g}=$hmm_acceptor_score if($transcript_hmm_acceptor_score{$g}>$hmm_acceptor_score);
      $transcript_donor_score{$g}=$donor_score if($transcript_donor_score{$g}>$donor_score);
      $transcript_acceptor_score{$g}=$acceptor_score if($transcript_acceptor_score{$g}>$donor_score);
    }
  }
  print "$g $transcript_junction_score{$g} $donor_seq $transcript_donor_score{$g} $acceptor_seq $transcript_acceptor_score{$g} $transcript_hmm_score{$g} $transcript_hmm_donor_score{$g} $transcript_hmm_acceptor_score{$g}\n";
}  

