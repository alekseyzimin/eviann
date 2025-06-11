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


while($line=<STDIN>){
  @F=split(/\t/,$line,7);
  if($F[5]=~/N/){
    @f=split(/(\D)/,$F[5]);
    $offset=$F[3];
    if($line =~/XS:A:(\S)/){
    for($i=0;$i<$#f;$i+=2){
      if($f[$i+1] eq "M" || $f[$i+1] eq "S" ||  $f[$i+1] eq "I"){
        $offset+=$f[$i]; 
      }elsif($f[$i+1] eq "N" && $f[$i]>10){
        $junc{"$F[2]\t$offset\t".($offset+$f[$i])}++;
        $junc_ori{"$F[2]\t$offset\t".($offset+$f[$i])}=$1;
      }
    }
    }
  }
}
$n=0;
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
  
  print join("\t",@f),"\tJUNC$n\t$junc{$j}\t$dir\t",substr($genome_seqs{$f[0]},$f[1],2)," ",substr($genome_seqs{$f[0]},$f[2]-2,2)," $junc_ori{$j}\n";
  $n++;
}
