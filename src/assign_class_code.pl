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
      print "DEBUG $transcript best match $best_match\n";
      @f=split(/\s/,$best_match);
      $transcript_cds{$transcript}=$f[2];
      $transcript_cds_start{$transcript}=$protein_start{$f[2]};
      $transcript_cds_end{$transcript}=$protein_end{$f[2]};
      $transcript_class{$transcript}=$f[1];
      $transcript_ori{$transcript}=$f[3];
      #print "DEBUG trmap $transcript $f[1] $f[2] $protein_start{$f[2]} $protein_end{$f[2]}\n";
    }
    @F=split(/\s+/,$line);
    $transcript=substr($F[0],1);
    @junctions=();
    @f=split(/\-/,$F[3]);
    $tbeg{$transcript}=$f[0];
    $tend{$transcript}=$f[1];
    for(my $i=1;$i<$#f;$i++){
      push(@junctions,$f[$i]);
    }
    @matches=();
  }else{
    @f=split(/\t/,$line);
    @F=split(/\-/,$f[6]);
    my $pbeg=$F[0];
    my $pend=$F[-1];
    my $overhang_penalty=0;
    $overhang_penalty+=$pbeg-$tbeg{$transcript} if($pbeg<$tbeg{$transcript});
    $overhang_penalty+=$tend{$transcript}-$pend if($pend>$tend{$transcript});
    my $mismatch_penalty=0;
    $mismatch_penalty-=$pbeg-$tbeg{$transcript} if($pbeg>$tbeg{$transcript});
    $mismatch_penalty-=$tend{$transcript}-$pend if($pend<$tend{$transcript});

    my $common_junctions=0;
    for(my $i=1;$i<$#F;$i++){
      for(my $j=0;$j<=$#junctions;$j++){
        $common_junctions++ if($F[$i] eq $junctions[$j]);
      }
    }
    print "DEBUG trmap $transcript $f[5] $f[3] $f[4] $common_junctions $overhang_penalty $mismatch_penalty\n";
    if(defined($code_priority{$f[0]})){
      push(@matches,($code_priority{$f[0]}+$common_junctions*100000+$overhang_penalty*2+$mismatch_penalty)." $f[0] $f[5] $f[2]");
    }else{
      push(@matches,($f[4]-$f[3]+$common_junctions*100000+$overhang_penalty*2+$mismatch_penalty)." N $f[5] $f[2]");
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

