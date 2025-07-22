#!/usr/bin/env perl
my $transdecoder_gff=$ARGV[0];
my $input_gff=$ARGV[1];
open(FILE,"gffread --tlf $transdecoder_gff |");
while($line=<FILE>){
  @F=split(/\t/,$line);
  if($F[8]=~/CDS=(\d+):(\d+);CDSphase=/){
    $cds_s{$F[0]}=$1;
    $cds_e{$F[0]}=$2;
  }
}
open(FILE,"gffread --tlf $input_gff |");
while($line=<FILE>){
  chomp($line);
  @F=split(/\t/,$line);
  @attrs=split(/;/,$F[8]);
  $id=substr($attrs[0],3);
  if(defined($cds_s{$id})){
    ($cds_start,$cds_end)=find_cds_start_end(@attrs[2],$F[6],$cds_s{$id},$cds_e{$id});
    print "#FOUND $cds_s{$id} $cds_e{$id}\n";
    print join("\t",@F[0..7])."\t",join(";",@attrs[0..2]),";CDS=$cds_start:$cds_end;CDSphase=0;$attrs[3]\n";
  }else{
    print "$line\n";
  }
}

sub find_cds_start_end{
  my $exon_line=substr($_[0],6);
  my $ori=$_[1];
  my $cds_s=$_[2];
  my $cds_e=$_[3];
  my @exons=split(/,/,$exon_line);
  my $offset=1;
  my $cds_ts=0;
  my $cds_te=0;
  my $tl=0;
  for(my $i=0;$i<=$#exons;$i++){
     my ($s,$e)=split(/-/,$exons[$i]);
     $tl+=$e-$s+1;
  }
  print "#DEBUG transcript length $tl\n";
  if($ori eq "-"){
    $cds_s=$tl-$_[3]+1;
    $cds_e=$tl-$_[2]+1;
  }
  for(my $i=0;$i<=$#exons;$i++){
    my ($s,$e)=split(/-/,$exons[$i]);
    if($cds_s<=$e-$s+$offset){
      $cds_ts=$s+$cds_s-$offset;
      last;
    }else{
      $offset+=$e-$s+1;
    }
  }
  for(my $i=0;$i<=$#exons;$i++){
    my ($s,$e)=split(/-/,$exons[$i]);
    if($offset<=$cds_e && $cds_e<=$e-$s+$offset){
      $cds_te=$s+$cds_e-$offset;
      last;
    }else{
      $offset+=$e-$s+1;
    }
  }
  return($cds_ts,$cds_te);
}

    


