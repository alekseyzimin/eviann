#!/usr/bin/env perl
#This code adds gene and CDS record for every locus in GFF file
my $geneID="";
my @exons=();
my %transcript_gff;
my $protID="";
my $dir="";
my $scf="";
my $seq="";
#these are genetic codes for Markov chains
my %code=();
my %code2=();
my %code3=();
my @narray=("A","C","G","T");
#initialize code hashes
$n=0;
$n2=0;
$n3=0;
for($i=0;$i<4;$i++){
  $code{$narray[$i]}=$n;
  $n++;
  for($j=0;$j<4;$j++){
    $code2{"$narray[$i]$narray[$j]"}=$n2;
    $n2++;
    for($k=0;$k<4;$k++){
      $code3{"$narray[$i]$narray[$j]$narray[$k]"}=$n3;
      $n3++;
    }
  }
}
my $donor_length;
my $acceptor_length;
my $use_negatives=0;

#here we load up all transcripts that matched proteins
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
    $locID=substr($attributes[1],7);#this is the gene_id
    $geneID=substr($attributes[0],3);#this is the transcript_id
    $transcript_junction0_score{$geneID}=10000;
    $transcript_junction1_score{$geneID}=10000;
    $transcript_junction2_score{$geneID}=10000;
    $transcript_donor2_score{$geneID}=10000;
    $transcript_acceptor2_score{$geneID}=10000;
    $transcript_fix_score{$geneID}=10000;
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
if(-e $ARGV[2]){
  print "DEBUG Loading SNAP HMMs\n";
  open(FILE,$ARGV[2]);
  $line=<FILE>;
  if($line =~ /^zoeHMM/){#check format
    while($line=<FILE>){
      chomp($line);
      if($line=~/^Donor 0HMM/){
        my $i=0;
        while($line=<FILE>){
          last if($line=~/NN TRM/);
          chomp($line);
          $line=~s/^\s+//;
          my @f=split(/\s+/,$line);
          for(my $j=0;$j<4;$j++){
            $donor_freq[$i][$j]=$f[$j];
          }
          $i++;
        }
        $donor_length=$i;
      }elsif($line=~/^Donor 1HMM/){
        my $i=0;
        while($line=<FILE>){
          last if($line=~/NN TRM/);
          chomp($line);
          $line=~s/^\s+//;
          my @f=split(/\s+/,$line);
          for(my $j=0;$j<16;$j++){
            $donor_hmm_freq[$i][$j]=$f[$j];
          }
          $i++;
        }
      }elsif($line=~/^Donor 2HMM/){
        my $i=0;
        while($line=<FILE>){
          last if($line=~/NN TRM/);
          chomp($line); 
          $line=~s/^\s+//;
          my @f=split(/\s+/,$line);
          for(my $j=0;$j<64;$j++){
            $donor_hmm2_freq[$i][$j]=$f[$j];
          }
          $i++;
        }
      }elsif($line=~/^Acceptor 0HMM/){
        my $i=0;
        while($line=<FILE>){
          last if($line=~/NN TRM/);
          chomp($line);
          $line=~s/^\s+//;
          my @f=split(/\s+/,$line);
          for(my $j=0;$j<4;$j++){
            $acceptor_freq[$i][$j]=$f[$j];
          }
          $i++;
        }
        $acceptor_length=$i;
      }elsif($line=~/^Acceptor 1HMM/){
        my $i=0;
        while($line=<FILE>){
          last if($line=~/NN TRM/);
          chomp($line);
          $line=~s/^\s+//;
          my @f=split(/\s+/,$line);
          for(my $j=0;$j<16;$j++){
            $acceptor_hmm_freq[$i][$j]=$f[$j];
          }
          $i++;
        }
      }elsif($line=~/^Acceptor 2HMM/){
        my $i=0;
        while($line=<FILE>){
          last if($line=~/NN TRM/);
          chomp($line);
          $line=~s/^\s+//;
          my @f=split(/\s+/,$line);
          for(my $j=0;$j<64;$j++){
            $acceptor_hmm2_freq[$i][$j]=$f[$j];
          }
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

#we load scores for the negative model
if( -e $ARGV[3]){
  $use_negatives=1;
  print "DEBUG Loading NEG HMMs\n";
  open(FILE,$ARGV[3]);
  $line=<FILE>;
  if($line =~ /^zoeHMM/){#check format
    while($line=<FILE>){
      chomp($line);
      if($line=~/^Donor 0HMM/){
        my $i=0;
        while($line=<FILE>){
          last if($line=~/NN TRM/);
          chomp($line);
          $line=~s/^\s+//;
          my @f=split(/\s+/,$line);
          for(my $j=0;$j<4;$j++){
            $donor_nfreq[$i][$j]=$f[$j];
          }
          $i++;
        }
        $donor_length=$i;
      }elsif($line=~/^Donor 1HMM/){
        my $i=0;
        while($line=<FILE>){
          last if($line=~/NN TRM/);
          chomp($line);
          $line=~s/^\s+//;
          my @f=split(/\s+/,$line);
          for(my $j=0;$j<16;$j++){
            $donor_hmm_nfreq[$i][$j]=$f[$j];
          }
          $i++;
        }
      }elsif($line=~/^Donor 2HMM/){
        my $i=0;
        while($line=<FILE>){
          last if($line=~/NN TRM/);
          chomp($line);
          $line=~s/^\s+//;
          my @f=split(/\s+/,$line);
          for(my $j=0;$j<64;$j++){
            $donor_hmm2_nfreq[$i][$j]=$f[$j];
          }
          $i++;
        }
      }elsif($line=~/^Acceptor 0HMM/){
        my $i=0;
        while($line=<FILE>){
          last if($line=~/NN TRM/);
          chomp($line);
          $line=~s/^\s+//;
          my @f=split(/\s+/,$line);
          for(my $j=0;$j<4;$j++){
            $acceptor_nfreq[$i][$j]=$f[$j];
          }
          $i++;
        }
        $acceptor_length=$i;
      }elsif($line=~/^Acceptor 1HMM/){
        my $i=0;
        while($line=<FILE>){
          last if($line=~/NN TRM/);
          chomp($line);
          $line=~s/^\s+//;
          my @f=split(/\s+/,$line);
          for(my $j=0;$j<16;$j++){
            $acceptor_hmm_nfreq[$i][$j]=$f[$j];
          }
          $i++;
        }
      }elsif($line=~/^Acceptor 2HMM/){
        my $i=0;
        while($line=<FILE>){
          last if($line=~/NN TRM/);
          chomp($line);
          $line=~s/^\s+//;
          my @f=split(/\s+/,$line);
          for(my $j=0;$j<64;$j++){
            $acceptor_hmm2_nfreq[$i][$j]=$f[$j];
          }
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
    #print ${$transcript_gff{$g}}[$j],"\n";
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
      my $rel=2 if($reliable{"$gff_fields_prev[0] $gff_fields_prev[4] $gff_fields[3] $gff_fields[6]"}>=2);
      my $donor_score=0;
      my $acceptor_score=0;
      my $donor_hmm_score=0;
      my $acceptor_hmm_score=0;
      my $donor_hmm2_score=0;
      my $acceptor_hmm2_score=0;
      my $donor_nscore=0;
      my $acceptor_nscore=0;
      my $donor_hmm_nscore=0;
      my $acceptor_hmm_nscore=0;
      my $donor_hmm2_nscore=0;
      my $acceptor_hmm2_nscore=0;
      my $junction0_score=0;
      my $junction1_score=0;
      my $junction2_score=0;
      
      for(my $i=0;$i<$donor_length;$i++){
        $donor_score+=$donor_freq[$i][$code{substr($donor_seq,$i,1)}] if(defined($code{substr($donor_seq,$i,1)}));
      }  
      for(my $i=0;$i<$acceptor_length;$i++){
        $acceptor_score+=$acceptor_freq[$i][$code{substr($acceptor_seq,$i,1)}] if(defined($code{substr($acceptor_seq,$i,1)}));
      }

      for(my $i=0;$i<($donor_length-1);$i++){
        $donor_hmm_score+=$donor_hmm_freq[$i][$code2{substr($donor_seq,$i,2)}] if(defined($code2{substr($donor_seq,$i,2)}));
      }   
      $donor_hmm_score+=$donor_freq[0][$code{substr($donor_seq,0,1)}] if(defined($code{substr($donor_seq,0,1)}));
      for(my $i=0;$i<($acceptor_length-1);$i++){
        $acceptor_hmm_score+=$acceptor_hmm_freq[$i][$code2{substr($acceptor_seq,$i,2)}] if(defined($code2{substr($acceptor_seq,$i,2)}));
      }   
      $acceptor_hmm_score+=$acceptor_freq[0][$code{substr($acceptor_seq,0,1)}] if(defined($code{substr($acceptor_seq,0,1)}));

      for(my $i=0;$i<($donor_length-2);$i++){
        $donor_hmm2_score+=$donor_hmm2_freq[$i][$code3{substr($donor_seq,$i,3)}] if(defined($code3{substr($donor_seq,$i,3)}));
      }   
      $donor_hmm2_score+=$donor_hmm_freq[0][$code2{substr($donor_seq,0,2)}] if(defined($code2{substr($donor_seq,0,2)}));
      for(my $i=0;$i<($acceptor_length-2);$i++){
        $acceptor_hmm2_score+=$acceptor_hmm2_freq[$i][$code3{substr($acceptor_seq,$i,3)}] if(defined($code3{substr($acceptor_seq,$i,3)}));
      }   
      $acceptor_hmm2_score+=$acceptor_hmm_freq[0][$code2{substr($acceptor_seq,0,2)}] if(defined($code2{substr($acceptor_seq,0,2)}));

      #print "DEBUG D $donor_seq $donor_score $donor_hmm_score $donor_hmm2_score A $acceptor_seq $acceptor_score $acceptor_hmm_score $acceptor_hmm2_score\n";
      #here we calibrate the scores so that they are universal
      
      if($use_negatives == 0){
        $junction0_score=($donor_score+$acceptor_score)-4;
        $junction1_score=($donor_hmm_score+$acceptor_hmm_score)*0.45-4;
        $junction2_score=($donor_hmm2_score+$acceptor_hmm2_score)*.32-4;
      }else{
        for(my $i=0;$i<$donor_length;$i++){
          $donor_nscore+=$donor_nfreq[$i][$code{substr($donor_seq,$i,1)}] if(defined($code{substr($donor_seq,$i,1)}));
        }
        for(my $i=0;$i<$acceptor_length;$i++){
          $acceptor_nscore+=$acceptor_nfreq[$i][$code{substr($acceptor_seq,$i,1)}] if(defined($code{substr($acceptor_seq,$i,1)}));	
        }

        for(my $i=0;$i<($donor_length-1);$i++){
          $donor_hmm_nscore+=$donor_hmm_nfreq[$i][$code2{substr($donor_seq,$i,2)}] if(defined($code2{substr($donor_seq,$i,2)}));
        }
        $donor_hmm_nscore+=$donor_nfreq[0][$code{substr($donor_seq,0,1)}] if(defined($code{substr($donor_seq,0,1)}));
        for(my $i=0;$i<($acceptor_length-1);$i++){
          $acceptor_hmm_nscore+=$acceptor_hmm_nfreq[$i][$code2{substr($acceptor_seq,$i,2)}] if(defined($code2{substr($acceptor_seq,$i,2)}));
        }
        $acceptor_hmm_nscore+=$acceptor_nfreq[0][$code{substr($acceptor_seq,0,1)}] if(defined($code{substr($acceptor_seq,0,1)}));

        for(my $i=0;$i<($donor_length-2);$i++){
          $donor_hmm2_nscore+=$donor_hmm2_nfreq[$i][$code3{substr($donor_seq,$i,3)}] if(defined($code3{substr($donor_seq,$i,3)}));
        }
        $donor_hmm2_nscore+=$donor_hmm_nfreq[0][$code2{substr($donor_seq,0,2)}] if(defined($code2{substr($donor_seq,0,2)}));
        for(my $i=0;$i<($acceptor_length-2);$i++){
          $acceptor_hmm2_nscore+=$acceptor_hmm2_nfreq[$i][$code3{substr($acceptor_seq,$i,3)}] if(defined($code3{substr($acceptor_seq,$i,3)}));
        }
        $acceptor_hmm2_nscore+=$acceptor_hmm_nfreq[0][$code2{substr($acceptor_seq,0,2)}] if(defined($code2{substr($acceptor_seq,0,2)}));

        #print "DEBUG D $donor_seq $donor_score $donor_nscore $donor_hmm_score $donor_hmm_nscore $donor_hmm2_score $donor_hmm2_nscore A $acceptor_seq $acceptor_score $acceptor_nscore $acceptor_hmm_score $acceptor_hmm_nscore $acceptor_hmm2_score $acceptor_hmm2_nscore\n";
        $junction0_score=($donor_score-$donor_nscore+$acceptor_score-$acceptor_nscore)+3;
        $junction1_score=($donor_hmm_score-$donor_hmm_nscore+$acceptor_hmm_score-$acceptor_hmm_nscore)*0.45+3;
        $junction2_score=($donor_hmm2_score-$donor_hmm2_nscore+$acceptor_hmm2_score-$acceptor_hmm2_nscore)*.32+3;
      }
      my $fix_donor_score=defined($sdonor{substr($donor_seq,2,7)})?$sdonor{substr($donor_seq,2,7)}:-10000;
      my $fix_acceptor_score=defined($sacceptor{substr($acceptor_seq,$acceptor_length-9,7)})?$sacceptor{substr($acceptor_seq,$acceptor_length-9,7)}:-10000;
      my $fix_score=$fix_donor_score+$fix_acceptor_score;
      $transcript_junction0_score{$g}=$junction0_score if($transcript_junction0_score{$g}>$junction0_score);
      $transcript_junction1_score{$g}=$junction1_score if($transcript_junction1_score{$g}>$junction1_score);
      $transcript_junction2_score{$g}=$junction2_score if($transcript_junction2_score{$g}>$junction2_score);
      $transcript_donor2_score{$g}=$donor_hmm2_score if($transcript_donor2_score{$g}>$donor_hmm2_score);
      $transcript_acceptor2_score{$g}=$acceptor_hmm2_score if($transcript_acceptor2_score{$g}>$acceptor_hmm2_score);
      $transcript_fix_score{$g}=$fix_score if($transcript_fix_score{$g}>$fix_score);
    }
  }
  print "$g $transcript_fix_score{$g} $transcript_junction0_score{$g} $transcript_junction1_score{$g} $transcript_junction2_score{$g} $transcript_donor2_score{$g} $transcript_acceptor2_score{$g}\n";
}  

