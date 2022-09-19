#!/usr/bin/env perl
#this script uses splice junction site info from  gffread dump to remove extra splice junctions found by minimap
#
my $transcript_file=$ARGV[0];
open(FILE,$transcript_file);
while($line=<FILE>){
  my @splice_junctions=();
  chomp($line);
  @F=split(/\s+/,$line);
  if($F[0] =~ /^>/){
    if($#F==4){
      @f=split(/,/,$F[-1]);
      for($i=1;$i<=$#f;$i++){
        @ff=split(/-/,$f[$i]); 
        push(@splice_junctions,$ff[0]);
      }
      $junctions{substr($F[0],1)}=[@splice_junctions];
    }
  }
}

#we take sam file on STDIN
while($line=<STDIN>){
  if($line=~ /^@/){
    print $line;
    next;
  }
  chomp($line);
  @F=split(/\t/,$line);
  if(not(defined($junctions{$F[0]}))){
    print $line,"\n";
    next;
  }
  $offset=1;
  @f=split(/(\D)/,$F[5]);
  if($F[1]==0){
    for($i=0;$i<$#f;$i++){
      if($f[$i+1] eq "M" || $f[$i+1] eq "S" ||  $f[$i+1] eq "I"){
        $offset+=$f[$i];
      }elsif($f[$i+1] eq "N"){
        #found splice site; check if it exists
        my $found=0;
        foreach $j(@{$junctions{$F[0]}}){
          if($j==$offset){
            $found=1;
            last;
          }
        }
        $f[$i+1]="D" unless($found);
      }
    }
  }else{
    for($i=$#f-1;$i>0;$i--){
      if($f[$i+1] eq "M" || $f[$i+1] eq "S" || $f[$i+1] eq "I"){
        $offset+=$f[$i];
      }elsif($f[$i+1] eq "N"){
        my $found=0;
        foreach $j(@{$junctions{$F[0]}}){
          if($j==$offset){
            $found=1;
            last;
          }
        }
        $f[$i+1]="D" unless($found);
      }
    }
  }
  $F[5]=join("",@f);
  print join("\t",@F),"\n";
}
