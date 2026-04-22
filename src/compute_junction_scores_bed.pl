#!/usr/bin/env perl
#This code adds gene and CDS record for every locus in GFF file
my $geneID="";
my $scf="";
my $seq="";
my $donor_length=16;
my $acceptor_length=30;
my $exon_bases=3;
my $score_floor_value=1e-10;
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

#we load the genome sequences
open(FILE,$ARGV[0]);
#print "DEBUG Loading genome sequence\n";
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

$exon_bases=$ARGV[1] if(defined($ARGV[1]));

my $wd=1;
my $wa=1;
while($line=<STDIN>){
  chomp($line);
#print STDERR "DEBUG: line $line\n";
  my @bed_fields=split(/\t/,$line);
  die("Genome sequence $bed_fields[1] needed for transcript $g not found!") if(not(defined($genome_seqs{$bed_fields[1]})));
  if($bed_fields[3] eq "+"){
    if($bed_fields[0] eq "don"){
      $donor_seq=uc(substr($genome_seqs{$bed_fields[1]},$bed_fields[2]-$exon_bases,$donor_length));
    }elsif($bed_fields[0] eq "acc"){
      $acceptor_seq=uc(substr($genome_seqs{$bed_fields[1]},$bed_fields[2]-($acceptor_length-$exon_bases),$acceptor_length));
    }elsif($bed_fields[0] eq "pair"){
      $pair_seq=uc(substr($genome_seqs{$bed_fields[1]},$bed_fields[2]-$exon_bases,$donor_length))." ".uc(substr($genome_seqs{$bed_fields[1]},$bed_fields[3]-($acceptor_length-$exon_bases),$acceptor_length));
    }
  }else{
    if($bed_fields[0] eq "don"){
      $donor_seq=uc(substr($genome_seqs{$bed_fields[1]},$bed_fields[2]-($donor_length-$exon_bases),$donor_length));
      $donor_seq=~tr/ACGTNacgtn/TGCANtgcan/;
      $donor_seq=reverse($donor_seq);
    }elsif($bed_fields[0] eq "acc"){
      $acceptor_seq=uc(substr($genome_seqs{$bed_fields[1]},$bed_fields[2]-$exon_bases,$acceptor_length));
      $acceptor_seq=~tr/ACGTNacgtn/TGCANtgcan/;
      $acceptor_seq=reverse($acceptor_seq);
    }elsif($bed_fields[0] eq "pair"){
      $pair_seq=uc(substr($genome_seqs{$bed_fields[1]},$bed_fields[3]-$exon_bases,$acceptor_length))." ".uc(substr($genome_seqs{$bed_fields[1]},$bed_fields[2]-($donor_length-$exon_bases),$donor_length));
      $pair_seq=~tr/ACGTNacgtn /TGCANtgcan /;
      $pair_seq=reverse($pair_seq);
    } 
  }
  next if($donor_seq=~/N/ || $acceptor_seq=~/N/);

  if($bed_fields[0] eq "don"){
    print STDERR "DEBUG donor $donor_seq $bed_fields[3] $bed_fields[4] $bed_fields[1] $bed_fields[2]\n";
    $donor7{substr($donor_seq,2,7)}++;
    for(my $i=0;$i<$donor_length;$i++) {$donor_pwm[$i][$code{substr($donor_seq,$i,1)}]++ if(defined($code{substr($donor_seq,$i,1)}));}
    for(my $i=0;$i<($donor_length-1);$i++) {$donor2_pwm[$i][$code2{substr($donor_seq,$i,2)}]++ if(defined($code2{substr($donor_seq,$i,2)}));}
    for(my $i=0;$i<($donor_length-2);$i++) {$donor3_pwm[$i][$code3{substr($donor_seq,$i,3)}]++ if(defined($code3{substr($donor_seq,$i,3)}));}
    $wd++;
  }elsif($bed_fields[0] eq "acc"){
    print STDERR "DEBUG acceptor $acceptor_seq $bed_fields[3] $bed_fields[4] $bed_fields[1] $bed_fields[2]\n";
    $acceptor7{substr($acceptor_seq,$acceptor_length-9,7)}++;
    for(my $i=0;$i<$acceptor_length;$i++) {$acceptor_pwm[$i][$code{substr($acceptor_seq,$i,1)}]++ if(defined($code{substr($acceptor_seq,$i,1)}));}
    for(my $i=0;$i<($acceptor_length-1);$i++) {$acceptor2_pwm[$i][$code2{substr($acceptor_seq,$i,2)}]++ if(defined($code2{substr($acceptor_seq,$i,2)}));}
    for(my $i=0;$i<($acceptor_length-2);$i++) {$acceptor3_pwm[$i][$code3{substr($acceptor_seq,$i,3)}]++ if(defined($code3{substr($acceptor_seq,$i,3)}));}
    $wa++;
  }elsif($bed_fields[0] eq "pair"){
    my $max_len=int(abs($bed_fields[3]-$bed_fields[2])/2);
    print STDERR "DEBUG pair $pair_seq $bed_fields[4] $max_len $bed_fields[1] $bed_fields[2] $bed_fields[3]\n";
  }
}

#OUTPUT PWMs
print "zoeHMM\n";
print "Donor 0HMM\n";
for(my $i=0;$i<$donor_length;$i++){
  for(my $j=0;$j<4;$j++){
    printf("%.3f ", log($donor_pwm[$i][$j]/$wd*4+$score_floor_value));
  }
  print "\n";
}
print "NN TRM\n";
print "Acceptor 0HMM\n";
for(my $i=0;$i<$acceptor_length;$i++){
  for(my $j=0;$j<4;$j++){
    printf("%.3f ",log($acceptor_pwm[$i][$j]/$wa*4+$score_floor_value));
  } 
  print "\n";
} 
print "NN TRM\n";
print "Donor 1HMM\n";
for(my $i=0;$i<($donor_length-1);$i++){
  for(my $j=0;$j<16;$j++){
    printf("%.3f ", log($donor2_pwm[$i][$j]/$wd*16+$score_floor_value));
  }
  print "\n";
}
print "NN TRM\n";
print "Acceptor 1HMM\n";
for(my $i=0;$i<($acceptor_length-1);$i++){
  for(my $j=0;$j<16;$j++){
    printf("%.3f ",log($acceptor2_pwm[$i][$j]/$wa*16+$score_floor_value));
  }
  print "\n";
}
print "NN TRM\n";
print "Donor 2HMM\n";
for(my $i=0;$i<($donor_length-2);$i++){
  for(my $j=0;$j<64;$j++){
    printf("%.3f ", log($donor3_pwm[$i][$j]/$wd*64+$score_floor_value));
  }
  print "\n";
}
print "NN TRM\n";
print "Acceptor 2HMM\n";
for(my $i=0;$i<($acceptor_length-2);$i++){
  for(my $j=0;$j<64;$j++){
    printf("%.3f ",log($acceptor3_pwm[$i][$j]/$wa*64+$score_floor_value));
  }
  print "\n";
}
print "NN TRM\n";
my @keys = sort { $donor7{$b} <=> $donor7{$a} } keys(%donor7);
print "SDonor\n";
for my $k(@keys){
  print "$k ",log($donor7{$k}/$wd*8192+$score_floor_value),"\n";
}
print "NNNNNN\n";
@keys = sort { $acceptor7{$b} <=> $acceptor7{$a} } keys(%acceptor7);
print "SAcceptor\n";
for my $k(@keys){
  print "$k ",log($acceptor7{$k}/$wa*8192+$score_floor_value),"\n";
}
print "NNNNNN\n";
