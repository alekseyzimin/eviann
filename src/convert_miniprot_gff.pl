#!/usr/bin/env perl
#this code converts miniprot output to the same stadnard as exonerate, correcting for inclusion or exclusion of the stop codon
@output=();
while($line=<STDIN>){
  chomp($line);
  @F=split(/\t/,$line);
  if($F[2] eq "mRNA"){
    if(scalar(@output) > 1){
      foreach $l(@output){
        print $l;
      }
    }
    @output=();
    $F[2]="gene";
    if($F[8] =~ /^ID=(\S+);Rank=(\S+);Identity=(\d+.\d+);Positive=(\d+.\d+);.*Target=(\S+)\s\d+\s\d+$/){
      $count=1;
      $chrom=$F[0];
      $target=$5;
      $identity=$4*100;
      $similarity=$3*100;
      $chrom=~s/;/_/g;
      $parent="$target:$chrom:$F[3]:$identity:$similarity";
      if(not(defined($already_output{$parent}))){
        $flag=1;
        $already_output{$parent}=1;
      }else{
        $flag=0;
      }
      push(@output,join("\t",@F[0..7])."\tID=$parent;geneID=$parent;identity=$identity;similarity=$similarity\n") if($flag);
    }
  }elsif(uc($F[2]) eq "CDS"){
    if($count>1){
      $F[2]="intron";
      if($F[6] eq "+"){
        $start=$last_base+1;
        $stop=$F[3]-1;
      }else{
        $start=$F[4]+1;
        $stop=$last_base-1;
      }
      push(@output,join("\t",@F[0..1])."\tintron\t$start\t$stop\t".join("\t",@F[5..7])."\tID=intron-$parent-",$count-1,";Parent=$parent\n") if($flag);
    }
    $F[2]="cds";
    push(@output,join("\t",@F[0..7])."\tID=cds-$parent-$count;Parent=$parent\n") if($flag);
    $F[2]="exon";
    push(@output,join("\t",@F[0..7])."\tID=exon-$parent-$count;Parent=$parent\n") if($flag);
    $count++;
    $last_base=($F[6] eq "+")?$F[4]:$F[3];
  }elsif($F[2] eq "stop_codon"){
    @f=split(/\t/,$output[$#output]);
    if($F[6] eq "+"){
      if($F[4]==$f[4]){
        $f[4]-=3;
        $output[$#output]=join("\t",@f);
        @f=split(/\t/,$output[$#output-1]);
        $f[4]-=3;
        $output[$#output-1]=join("\t",@f);
        @f=split(/\t/,$output[0]);
        $f[4]-=3;
        $output[0]=join("\t",@f);
      }
    }else{
      if($F[3]==$f[3]){
        $f[3]+=3;
        $output[$#output]=join("\t",@f);
        @f=split(/\t/,$output[$#output-1]);
        $f[3]+=3;
        $output[$#output-1]=join("\t",@f);
        @f=split(/\t/,$output[0]);
        $f[3]+=3;
        $output[0]=join("\t",@f);
      }
    }
    push(@output,join("\t",@F[0..7])."\tID=stop_codon-$parent-$count;Parent=$parent\n") if($flag);
  }
}
if(scalar(@output) > 1){
  foreach $l(@output){
    print $l;
  }
}

