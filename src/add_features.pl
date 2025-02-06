#!/usr/bin/env perl
$gff_file=$ARGV[0];
$trna_file=$ARGV[1];

#read in gff file with gene records, assume sorted
if(not($gff_file eq "")){
  open(FILE,$gff_file);
  while(my $line=<FILE>){#we just read in the whole file
    chomp($line);
    if($line =~ /^\#/){ print $line,"\n";next;}
    my @gff_fields=split(/\t/,$line);
    if($gff_fields[2] eq "gene"){
      $gene_beg=$gff_fields[3];
      $gene_end=$gff_fields[4];
      $gene_ori=$gff_fields[6];
      $gene_chr=$gff_fields[0];
      push(@gchr,$gene_chr);
      push(@gbeg,$gene_beg);
      push(@gend,$gene_end);
      push(@gori,$gene_ori);
      $dir_factor=$gene_ori eq "+" ? 0 : 0.5;
      $gene_key=$gene_chr." ".($gene_beg+$dir_factor);
      push(@gene_keys,$gene_key);
    }
    $gene_record{$gene_key}.=$line."\n";
  }
}

#read in trnas from external
if(not($trna_file eq "")){
  open(FILE,$trna_file);
  my $contain=0;
  while(my $line=<FILE>){#we just read in the whole file
    chomp($line);
    my @gff_fields=split(/\t/,$line);
    if($gff_fields[2] eq "gene"){
      $contain=0;
      $gene_beg=$gff_fields[3];
      $gene_end=$gff_fields[4];
      $gene_ori=$gff_fields[6];
      $gene_chr=$gff_fields[0];
      for(my $i=0;$i<=$#gchr;$i++){
        if($gene_chr eq $gchr[$i] && $gene_ori eq $gori[$i] && (($gene_beg<=$gend[$i] && $gene_beg>=$gbeg[$i])|| ($gene_end<=$gend[$i] && $gene_end>=$gbeg[$i]) || ($gene_beg<=$gbeg[$i] && $gene_end>=$gend[$i]))){
          $contain=1;
          $i=$#gchr+1;
        }
      }
      if(!$contain){
        $dir_factor=$gene_ori eq "+" ? 0 : 0.5;
        $gene_key="$gene_chr ".($gene_beg+$dir_factor);
        push(@gene_keys,$gene_key);
        $gene_record{$gene_key}.=$line."\n";
      }
    }elsif(!$contain){
      $gene_record{$gene_key}.=$line."\n";
    }
  }
}

#now we added non-intersecting gene records, lets sort and output
my @gene_keys_sorted=sort mysort @gene_keys;
foreach my $g(@gene_keys_sorted){
  print $gene_record{$g};
}

sub mysort{
  my ($chroma,$coorda)=split(/\s+/,$a);
  my ($chromb,$coordb)=split(/\s+/,$b);
  return($chroma cmp $chromb || $coorda <=>$coordb);
}

