#!/usr/bin/env perl
#read trmap file
my $trmap=$ARGV[0];
my %code_priority;
$code_priority{"="}=5000000000;
$code_priority{"k"}=4000000000;
$code_priority{"j"}=3000000000;
$code_priority{"n"}=2000000000;
$code_priority{"m"}=1000000000;
my $transcript="";
my @matches=();
open(FILE,$trmap);
#print "DEBUG reading trmap\n";
while($line=<FILE>){
  chomp($line);
  if($line=~/^>/){
    unless($transcript eq ""){
      $best_match=(sort { (split(/\s/,$b))[0]<=>(split(/\s/,$a))[0] } @matches)[0];
      @f=split(/\s/,$best_match);
      $transcript_cds{$transcript}=$f[2];
      $transcript_cds_start{$transcript}=$protein_start{$f[2]};
      $transcript_cds_end{$transcript}=$protein_end{$f[2]};
      $transcript_class{$transcript}=$f[1];
      $transcript_ori{$transcript}=$f[3];
      #print "DEBUG trmap $transcript $f[1] $f[2] $protein_start{$f[2]} $protein_end{$f[2]}\n";
    }
    $transcript=substr((split(/\s+/,$line))[0],1);
    @matches=();
  }else{
    @f=split(/\t/,$line);
    if(defined($code_priority{$f[0]})){
      push(@matches,($code_priority{$f[0]}+$f[4]-$F[3])." $f[0] $f[5] $f[2]");
    }else{
      push(@matches,($f[4]-$F[3])." N $f[5] $f[2]");
    }
  }
}
unless($transcript eq ""){
  $best_match=(sort { (split(/\s/,$b))[0]<=>(split(/\s/,$a))[0] } @matches)[0];
  @f=split(/\s/,$best_match);
  $transcript_cds{$transcript}=$f[2];
  $transcript_cds_start{$transcript}=$protein_start{$f[2]};
  $transcript_cds_end{$transcript}=$protein_end{$f[2]};
  $transcript_class{$transcript}=$f[1];
  $transcript_ori{$transcript}=$f[3];
}

#add transcript class to gtf file
while($line=<STDIN>){
  chomp($line);
  @F=split(/\t/,$line);
  if($F[2] eq "transcript"){
    $F[8].=";" unless($F[8] =~ /;$/);
    $tid=$1 if($F[8] =~ /^transcript_id "(\S+)"/);
    if(defined($transcript_class{$tid})){
      $F[8].=" class_code \"$transcript_class{$tid}\"; cmp_ref \"$transcript_cds{$tid}\";";
      $F[6]=$transcript_ori{$tid} if($F[6] eq ".");
    }else{
      $F[8].=" class_code \"u\";";
    }
  }else{
    $F[6]=$transcript_ori{$tid} if($F[6] eq "." && defined($transcript_ori{$tid}));
  }
  print join("\t",@F),"\n";
}

