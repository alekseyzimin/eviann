#!/usr/bin/env perl
#load genome sequences
open(FILE,$ARGV[0]);
while(my $line=<FILE>){
  chomp($line);
  if($line=~ /^>/){
    if(not($scf eq "")){
      $genome_seqs{$scf}=$seq;
      $seq="";
    }
    $scf=substr((split(/\s+/,$line))[0],1);
  }else{
    $seq.=$line;
  } 
}   
$genome_seqs{$scf}=$seq if(not($scf eq ""));

my $lineno=1;
open(FILE,"samtools view $ARGV[1] |");
while($line=<FILE>){
  @F=split(/\t/,$line,7);
  if($F[5]=~/N/){
    @f=split(/(\D)/,$F[5]);
    $offset=$F[3];
    $ori=".";
    $ori=$1 if($line =~/XS:A:(\S)/);
    for($i=0;$i<$#f;$i+=2){
      if($f[$i+1] eq "M" || $f[$i+1] eq "X" ||  $f[$i+1] eq "D" || $f[$i+1] eq "="){
        $offset+=$f[$i]; 
      }elsif($f[$i+1] eq "N"){
        $junc{"$F[2]\t$offset\t".($offset+$f[$i])}++;
        $junc_ori{"$F[2]\t$offset\t".($offset+$f[$i])}=$ori;
        $junc_line{"$F[2]\t$offset\t".($offset+$f[$i])}.="$lineno ";
	$offset+=$f[$i];
      }
    }
  }
  $lineno++;
}

my %bad_lines=();
foreach $j(keys %junc){
  @f=split(/\t/,$j);
  $orif=0;
  $orir=0;
  $f[1]--;
  $f[2]--;
  next unless(uc(substr($genome_seqs{$f[0]},$f[1]+1,1)) eq "T" || uc(substr($genome_seqs{$f[0]},$f[2]-2,1)) eq "A");
  
  if(uc(substr($genome_seqs{$f[0]},$f[1],1)) eq "C"){
    $orir++;
  }elsif(uc(substr($genome_seqs{$f[0]},$f[1],1)) eq "G"){
    $orif++;
  }
  if(uc(substr($genome_seqs{$f[0]},$f[2]-1,1)) eq "C"){
    $orir++;
  }elsif(uc(substr($genome_seqs{$f[0]},$f[2]-1,1)) eq "G"){
    $orif++;
  }
  next if($orif==$orir);

  $dir=$orif>$orir ? "+" : "-";
  
  if(not($dir eq $junc_ori{$j}) && $junc{$j}<5){
    @f=split(/\s/,$junc_line{$j});
    foreach my $n(@f){
      $bad_lines{$n}=1;
    }
  }
}

$lineno=1;
open(FILE,"samtools view -h $ARGV[1] |");
while($line=<FILE>){
  if($line =~ /^(@|#)/){
    print $line;
    next;
  }else{
    print $line unless($bad_lines{$lineno});
    $lineno++;
  }
}
