#!/usr/bin/env perl
#This code adds gene and CDS record for every locus in GFF file
my $geneID="";
my @exons=();
my %transcript_gff;
my $protID="";
my $dir="";
my $scf="";
my $seq="";
my $donor_length=16;
my $acceptor_length=30;
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

#here we load up all transcripts that matched proteins
while(my $line=<STDIN>){
  chomp($line);
  my @gff_fields=split(/\t/,$line);
  my @attributes=split(";",$gff_fields[8]);
  if($gff_fields[2] eq "transcript" || $gff_fields[2] eq "mRNA"){
    if(defined($transcript{$geneID})){
      $transcript_gff{$geneID}=[@exons];
      $transcript_cds_start{$geneID}=$CDS_start;
      $transcript_cds_stop{$geneID}=$CDS_stop;
    }
    $CDS_start=-1;
    $CDS_stop=-1;
    $first_exon_start=-1;
    $first_exon_end=-1;
    $last_exon_end=-1;
    $last_exon_start=-1;
    @exons=();
    $locID=substr($attributes[1],7);#this is the gene_id
    $geneID=substr($attributes[0],3);#this is the transcript_id
    $transcript{$geneID}=$line;
  }elsif($gff_fields[2] eq "exon"){
    push(@exons,$line) if(defined($transcript{$geneID}));
    $first_exon_start=$gff_fields[3] if($first_exon_start==-1);
    $first_exon_end=$gff_fields[4] if($first_exon_end==-1);
    $last_exon_start=$gff_fields[3];
    $last_exon_end=$gff_fields[4];
#for now only use CDS that start or end inside the first/last exon
  }elsif(uc($gff_fields[2]) eq "CDS"){
    if($gff_fields[6] eq "+"){
      $CDS_start=$gff_fields[3] if($gff_fields[4]==$first_exon_end && $gff_fields[4]-$gff_fields[3]>=6 && $gff_fields[3]-$first_exon_start>=9);
      $CDS_stop=$gff_fields[4] if($gff_fields[3]==$last_exon_start);
    }else{
      $CDS_stop=$gff_fields[3] if($gff_fields[4]==$first_exon_end);
      $CDS_start=$gff_fields[4] if($gff_fields[3]==$last_exon_start && $gff_fields[4]-$gff_fields[3]>=6 && $last_exon_end-$gff_fields[4]>=9);
    }
  }
}
if(defined($transcript{$geneID})){
  $transcript_gff{$geneID}=[@exons];
  $transcript_cds_start{$geneID}=$CDS_start;
  $transcript_cds_stop{$geneID}=$CDS_stop;
}
@exons=();

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

#we make the transcript sequences for protein coding transcripts and score the transcripts with HMMs
my $w=1;
my $sw=1;
for my $g(keys %transcript_gff){
  #print "DEBUG: processing transcript $g\n";
  my @gff_fields=();
  for(my $j=1;$j<=$#{$transcript_gff{$g}};$j++){
    @gff_fields=split(/\t/,${$transcript_gff{$g}}[$j]);
    die("Genome sequence $gff_fields[0] needed for transcript $g not found!") if(not(defined($genome_seqs{$gff_fields[0]})));
    if($j>0){
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
      $donor6{substr($donor_seq,2,1).substr($donor_seq,5,4)}++;
      $acceptor6{substr($acceptor_seq,$acceptor_length-9,4).substr($acceptor_seq,$acceptor_length-3,1)}++;
      print STDERR "DEBUG donor $donor_seq acceptor $acceptor_seq $gff_fields[6]\n";
      for(my $i=0;$i<$donor_length;$i++) {$donor_pwm[$i][$code{substr($donor_seq,$i,1)}]++;}
      for(my $i=0;$i<$acceptor_length;$i++) {$acceptor_pwm[$i][$code{substr($acceptor_seq,$i,1)}]++;}
      for(my $i=0;$i<($donor_length-1);$i++) {my $index=16; $index=$code2{substr($donor_seq,$i,2)} if(defined($code2{substr($donor_seq,$i,2)}));$donor_cc_pwm[$i][$index]++;}
      for(my $i=0;$i<($acceptor_length-1);$i++) {my $index=16; $index=$code2{substr($acceptor_seq,$i,2)} if(defined($code2{substr($acceptor_seq,$i,2)}));$acceptor_cc_pwm[$i][$index]++;}

      $w++;
    }
  }
  #this is in general incorrect because this ignores splice sites -- need to make transcript sequences
  if($transcript_cds_start{$g}>-1){
    @gff_fields=split(/\t/,$transcript{$g});
    if($gff_fields[6] eq "+"){
      $start_seq=uc(substr($genome_seqs{$gff_fields[0]},$transcript_cds_start{$g}-10,15));
    }else{
      $start_seq=uc(substr($genome_seqs{$gff_fields[0]},$transcript_cds_start{$g}-6,15));
      $start_seq=~tr/ACGTNacgtn/TGCANtgcan/;
      $start_seq=reverse($start_seq);
    }
    for(my $i=0;$i<15;$i++) {$start_pwm[$i][$code{substr($start_seq,$i,1)}]++;}
#print "DEBUG start $start_seq  $gff_fields[6]\n";
    $sw++;
  }
}  

#OUTPUT PWMs
print "zoeHMM\n";
print "Donor WMM\n";
for(my $i=0;$i<$donor_length;$i++){
  for(my $j=0;$j<4;$j++){
    printf("%.3f ", log($donor_pwm[$i][$j]/$w*4+1e-10));
  }
  print "\n";
}
print "NN TRM\n";
print "Acceptor WMM\n";
for(my $i=0;$i<$acceptor_length;$i++){
  for(my $j=0;$j<4;$j++){
    printf("%.3f ",log($acceptor_pwm[$i][$j]/$w*4+1e-10));
  } 
  print "\n";
} 
print "NN TRM\n";
print "Donor WAM\n";
for(my $i=0;$i<($donor_length-1);$i++){
  for(my $j=0;$j<16;$j++){
    printf("%.3f ", log($donor_cc_pwm[$i][$j]/$w*16+1e-10));
  }
  print "\n";
}
print "NN TRM\n";
print "Acceptor WAM\n";
for(my $i=0;$i<($acceptor_length-1);$i++){
  for(my $j=0;$j<16;$j++){
    printf("%.3f ",log($acceptor_cc_pwm[$i][$j]/$w*16+1e-10));
  }
  print "\n";
}
print "NN TRM\n";
print "Start ATG WMM\n";
for(my $i=0;$i<15;$i++){
  for(my $j=0;$j<4;$j++){
    printf("%.3f ",log($start_pwm[$i][$j]/$sw*4+1e-10));
  }
  print "\n";
}
print "NNN TRM\n";
my @keys = sort { $donor6{$b} <=> $donor6{$a} } keys(%donor6);
print "SDonor\n";
for my $k(@keys){
  print "$k ",log($donor6{$k}/$w*10000),"\n";
}
print "NNNNNN\n";
@keys = sort { $acceptor6{$b} <=> $acceptor6{$a} } keys(%acceptor6);
print "SAcceptor\n";
for my $k(@keys){
  print "$k ",log($acceptor6{$k}/$w*10000),"\n";
}
print "NNNNNN\n";


