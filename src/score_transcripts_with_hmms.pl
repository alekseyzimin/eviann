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
$code2{"AA"}=0;
$code2{"AC"}=1;
$code2{"AG"}=2;
$code2{"AT"}=3;
$code2{"CA"}=4;
$code2{"CC"}=5;
$code2{"CG"}=6;
$code2{"CT"}=7;
$code2{"GA"}=8;
$code2{"GC"}=9;
$code2{"GG"}=10;
$code2{"GT"}=11;
$code2{"TA"}=12;
$code2{"TC"}=13;
$code2{"TG"}=14;
$code2{"TT"}=15;
my $donor_length=16;
my $acceptor_length=30;

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
      if($line=~/^Donor WMM/){
        my $i=0;
        while($line=<FILE>){
          last if($line=~/NN TRM/);
          chomp($line);
          $line=~s/^\s+//;
          my @f=split(/\s+/,$line);
          for(my $j=0;$j<4;$j++){
            $donor_freq[$i][$j]=$f[$j];
          }
          $donor_freq[$i][4]=0;
          $i++;
        }
      }elsif($line=~/^Donor WAM/){
        my $i=0;
        while($line=<FILE>){
          last if($line=~/NN TRM/);
          chomp($line);
          $line=~s/^\s+//;
          my @f=split(/\s+/,$line);
          for(my $j=0;$j<16;$j++){
            $donor_hmm_freq[$i][$j]=$f[$j];
          }
          $donor_hmm_freq[$i][16]=0;
          $i++;
        }
      }elsif($line=~/^Acceptor WMM/){
        my $i=0;
        while($line=<FILE>){
          last if($line=~/NN TRM/);
          chomp($line);
          $line=~s/^\s+//;
          my @f=split(/\s+/,$line);
          for(my $j=0;$j<4;$j++){
            $acceptor_freq[$i][$j]=$f[$j];
          }
          $acceptor_freq[$i][4]=0;
          $i++;
        }
      }elsif($line=~/^Acceptor WAM/){
        my $i=0;
        while($line=<FILE>){
          last if($line=~/NN TRM/);
          chomp($line);
          $line=~s/^\s+//;
          my @f=split(/\s+/,$line);
          for(my $j=0;$j<16;$j++){
            $acceptor_hmm_freq[$i][$j]=$f[$j];
          }
          $acceptor_hmm_freq[$i][16]=0;
          $i++;
        } 
      }elsif($line=~/^Start/){
        my $i=0;
        while($line=<FILE>){
          last if($line=~/NNN TRM/);
          chomp($line);
          $line=~s/^\s+//;
          my @f=split(/\s+/,$line);
          for(my $j=0;$j<4;$j++){
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
        $donor_seq=uc(substr($genome_seqs{$gff_fields[0]},$gff_fields_prev[4]-3,$donor_length));
        $acceptor_seq=uc(substr($genome_seqs{$gff_fields[0]},$gff_fields[3]-($acceptor_length-2),$acceptor_length));
      }else{
        $donor_seq=uc(substr($genome_seqs{$gff_fields[0]},$gff_fields[3]-($donor_length-2),$donor_length));
        $acceptor_seq=uc(substr($genome_seqs{$gff_fields[0]},$gff_fields_prev[4]-3,$acceptor_length));
        $donor_seq=~tr/ACGTNacgtn/TGCANtgcan/;
        $donor_seq=reverse($donor_seq);
        $acceptor_seq=~tr/ACGTNacgtn/TGCANtgcan/;
        $acceptor_seq=reverse($acceptor_seq);
      }
      my $donor_score=0;
      my $acceptor_score=0;
      my $donor_hmm_score=0;
      my $acceptor_hmm_score=0;
      for(my $i=0;$i<$donor_length;$i++){
        $donor_score+=$donor_freq[$i][$code{substr($donor_seq,$i,1)}];
      }
      for(my $i=0;$i<$acceptor_length;$i++){
        $acceptor_score+=$acceptor_freq[$i][$code{substr($acceptor_seq,$i,1)}];
      }
      for(my $i=0;$i<($donor_length-1);$i++){
        $index=16;
        $index=$code2{substr($donor_seq,$i,2)} if(defined($code2{substr($donor_seq,$i,2)}));
        $donor_hmm_score+=$donor_hmm_freq[$i][$index];
      }
      $donor_hmm_score+=$donor_freq[0][$code{substr($donor_seq,0,1)}];
      for(my $i=0;$i<($acceptor_length-1);$i++){
        $index=16;
        $index=$code2{substr($acceptor_seq,$i,2)} if(defined($code2{substr($acceptor_seq,$i,2)}));
        $acceptor_hmm_score+=$acceptor_hmm_freq[$i][$index];
      }
      $acceptor_hmm_score+=$acceptor_freq[$acceptor_length-1][$code{substr($acceptor_seq,0,1)}];
      #print "DEBUG $donor_seq $donor_score $donor_hmm_score $acceptor_seq $acceptor_score $acceptor_hmm_score\n";
      #this here combined the PWM and WAM scores
      my $junction_score=($donor_hmm_score+$acceptor_hmm_score)*0.46;
      #my $junction_score=($donor_score+$acceptor_score)*0.25+(($donor_hmm_score+$acceptor_hmm_score)*0.5)*0.75;
      my $hmm_donor_score=defined($sdonor{substr($donor_seq,2,1).substr($donor_seq,5,4)})?$sdonor{substr($donor_seq,2,1).substr($donor_seq,5,4)}:-10000;
      my $hmm_acceptor_score=defined($sacceptor{substr($acceptor_seq,21,4).substr($acceptor_seq,27,1)})?$sacceptor{substr($acceptor_seq,21,4).substr($acceptor_seq,27,1)}:-10000;
      my $hmm_score=$hmm_donor_score+$hmm_acceptor_score;
      $transcript_junction_score{$g}=$junction_score if($transcript_junction_score{$g}>$junction_score);
      $transcript_hmm_score{$g}=$hmm_score if($transcript_hmm_score{$g}>$hmm_score);
      $transcript_hmm_donor_score{$g}=$hmm_donor_score if($transcript_hmm_donor_score{$g}>$hmm_donor_score);
      $transcript_hmm_acceptor_score{$g}=$hmm_acceptor_score if($transcript_hmm_acceptor_score{$g}>$hmm_acceptor_score);
      $transcript_donor_score{$g}=$donor_hmm_score if($transcript_donor_score{$g}>$donor_hmm_score);
      $transcript_acceptor_score{$g}=$acceptor_hmm_score if($transcript_acceptor_score{$g}>$acceptor_hmm_score);
    }
  }
  print "$g $transcript_junction_score{$g} $transcript_hmm_score{$g} $transcript_donor_score{$g} $transcript_acceptor_score{$g}\n";
}  

