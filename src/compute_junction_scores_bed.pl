#!/usr/bin/env perl
#This code adds gene and CDS record for every locus in GFF file
my $geneID="";
my @exons=();
my %transcript_gff;
my $protID="";
my $dir="";
my $scf="";
my $seq="";
my $donor_length=15;
my $acceptor_length=30;
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

#initialize PWMs
for(my $i=0;$i<$acceptor_length;$i++){
  for(my $j=0;$j<4;$j++){
    $donor_pwm[$i][$j]=0;
    $acceptor_pwm[$i][$j]=0;
    $start_pwm[$i][$j]=0;
  }
  for(my $j=0;$j<$acceptor_length;$j++){
    $donor_cc_pwm[$i][$j]=0;
    $acceptor_cc_pwm[$i][$j]=0;
  }
}

#NC_003071.7     8427861 8427939 JUNC00018209    1       +
#NC_003076.8     26939570        26939752        JUNC00062166    1       +
my $w=1;
while($line=<STDIN>){
  chomp($line);
#print "DEBUG: processing transcript $g\n";
  my @gff_fields=split(/\t/,$line);
  die("Genome sequence $gff_fields[0] needed for transcript $g not found!") if(not(defined($genome_seqs{$gff_fields[0]})));
  if($gff_fields[5] eq "+"){
    $donor_seq=uc(substr($genome_seqs{$gff_fields[0]},$gff_fields[1]-3,$donor_length));
    $acceptor_seq=uc(substr($genome_seqs{$gff_fields[0]},$gff_fields[2]-($acceptor_length-3),$acceptor_length));
  }else{
    $donor_seq=uc(substr($genome_seqs{$gff_fields[0]},$gff_fields[2]-($donor_length-3),$donor_length));
    $acceptor_seq=uc(substr($genome_seqs{$gff_fields[0]},$gff_fields[1]-3,$acceptor_length));
    $donor_seq=~tr/ACGTNacgtn/TGCANtgcan/;
    $donor_seq=reverse($donor_seq);
    $acceptor_seq=~tr/ACGTNacgtn/TGCANtgcan/;
    $acceptor_seq=reverse($acceptor_seq);
  }
  next if($donor_seq=~/N/ || $acceptor_seq=~/N/);
  $donor7{substr($donor_seq,2,7)}++;
  $acceptor7{substr($acceptor_seq,$acceptor_length-9,7)}++;
  print STDERR "DEBUG donor $donor_seq acceptor $acceptor_seq $gff_fields[5]\n";
  for(my $i=0;$i<$donor_length;$i++) {$donor_pwm[$i][$code{substr($donor_seq,$i,1)}]++ if(defined($code{substr($donor_seq,$i,1)}));}
  for(my $i=0;$i<$acceptor_length;$i++) {$acceptor_pwm[$i][$code{substr($acceptor_seq,$i,1)}]++ if(defined($code{substr($acceptor_seq,$i,1)}));}
  for(my $i=0;$i<($donor_length-1);$i++) {$donor2_pwm[$i][$code2{substr($donor_seq,$i,2)}]++ if(defined($code2{substr($donor_seq,$i,2)}));}
  for(my $i=0;$i<($acceptor_length-1);$i++) {$acceptor2_pwm[$i][$code2{substr($acceptor_seq,$i,2)}]++ if(defined($code2{substr($acceptor_seq,$i,2)}));}
  for(my $i=0;$i<($donor_length-2);$i++) {$donor3_pwm[$i][$code3{substr($donor_seq,$i,3)}]++ if(defined($code3{substr($donor_seq,$i,3)}));}
  for(my $i=0;$i<($acceptor_length-2);$i++) {$acceptor3_pwm[$i][$code3{substr($acceptor_seq,$i,3)}]++ if(defined($code3{substr($acceptor_seq,$i,3)}));}
  $w++;
}

#OUTPUT PWMs
print "zoeHMM\n";
print "Donor 0HMM\n";
for(my $i=0;$i<$donor_length;$i++){
  for(my $j=0;$j<4;$j++){
    printf("%.3f ", log($donor_pwm[$i][$j]/$w*4+$score_floor_value));
  }
  print "\n";
}
print "NN TRM\n";
print "Acceptor 0HMM\n";
for(my $i=0;$i<$acceptor_length;$i++){
  for(my $j=0;$j<4;$j++){
    printf("%.3f ",log($acceptor_pwm[$i][$j]/$w*4+$score_floor_value));
  } 
  print "\n";
} 
print "NN TRM\n";
print "Donor 1HMM\n";
for(my $i=0;$i<($donor_length-1);$i++){
  for(my $j=0;$j<16;$j++){
    printf("%.3f ", log($donor2_pwm[$i][$j]/$w*16+$score_floor_value));
  }
  print "\n";
}
print "NN TRM\n";
print "Acceptor 1HMM\n";
for(my $i=0;$i<($acceptor_length-1);$i++){
  for(my $j=0;$j<16;$j++){
    printf("%.3f ",log($acceptor2_pwm[$i][$j]/$w*16+$score_floor_value));
  }
  print "\n";
}
print "NN TRM\n";
print "Donor 2HMM\n";
for(my $i=0;$i<($donor_length-2);$i++){
  for(my $j=0;$j<64;$j++){
    printf("%.3f ", log($donor3_pwm[$i][$j]/$w*64+$score_floor_value));
  }
  print "\n";
}
print "NN TRM\n";
print "Acceptor 2HMM\n";
for(my $i=0;$i<($acceptor_length-2);$i++){
  for(my $j=0;$j<64;$j++){
    printf("%.3f ",log($acceptor3_pwm[$i][$j]/$w*64+$score_floor_value));
  }
  print "\n";
}
print "NN TRM\n";
my @keys = sort { $donor7{$b} <=> $donor7{$a} } keys(%donor7);
print "SDonor\n";
for my $k(@keys){
  print "$k ",log($donor7{$k}/$w*8192),"\n";
}
print "NNNNNN\n";
@keys = sort { $acceptor7{$b} <=> $acceptor7{$a} } keys(%acceptor7);
print "SAcceptor\n";
for my $k(@keys){
  print "$k ",log($acceptor7{$k}/$w*8192),"\n";
}
print "NNNNNN\n";


