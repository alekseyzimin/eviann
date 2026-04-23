#!/usr/bin/env perl
#transcript filtering by score
#
my $prefix=$ARGV[0];
my %score, %reliable;

die("File $prefix.num_introns.txt empty or does not exist") if(not( -s "$prefix.num_introns.txt")); 
open(FILE,"$prefix.num_introns.txt");
my $num_introns=int(<FILE>);
my $index=2;
$index++ if($num_introns >= 1024);
$index++ if($num_introns >= 4096);
open(FILE,"$prefix.transcript_splice_scores.txt");
while($line=<FILE>){
  chomp($line);
  my @f=split(/\s+/,$line);
  $score{$f[0]}=$f[$index];
}

die("File $prefix.num_introns.txt does not exist") if(not( -e "$prefix.reliable_transcripts_proteins.txt"));
open(FILE,"$prefix.reliable_transcripts_proteins.txt");
while($line=<FILE>){
  chomp($line);
  @f=split(/\s+/,$line);
  $reliable{$f[0]}=10;
  $reliable_t{$f[0]}=30 if($f[-1] =~ /=/);
}

my @gff;
while($line=<STDIN>){#tlf format  on STDIN
  chomp($line);
  my @F=split(/\t/,$line);
#first determine the optimal threshold
  push(@gff,join("\t",@F));
  if($F[8] =~ /^ID=(\S+);exonCount=(\S+);exons=(\S+);geneID=(\S+)/){
    push(@scores,$score{$1}) if(defined($reliable_t{$1}) && defined($score{$1}) );
  }
}
      
@scores_sorted=sort {$b <=> $a} @scores;
my $threshold=$scores_sorted[int($#scores_sorted*.999)];
print STDERR "$threshold\n";
      
foreach $t(@gff){
  my @F=split(/\t/,$t);
  if($F[8] =~ /^ID=(\S+);exonCount=(\S+);exons=(\S+);geneID=(\S+)/){
    my $id=$1;
    my $exons=$3;
    my $geneid=$4;
    my ($name,$samples,$tpm)=split(/:/,$id);
    $score{$id}=30 if(not(defined($score{$id})) || $name =~ /^REFSTRG/);
    $score{$id}+=$samples if($samples>1);
    $score{$id}+=$reliable_t{$id};
    if($score{$id}>$threshold){
      @f=split(/-/,$exons);
      if($#f>1){
        $transcripts{join("-",@f[1..$#f-1])}.="$id $F[0] $F[6] $f[0] $f[-1] $geneid ";
      }elsif(not(defined($output{"$F[0] $F[6] $exons"}))){
        $output{"$F[0] $F[6] $exons"}=1;
        print join("\t",@F),"\n";
      }
    }
  }
}

#average UTRs for multi-exon transcripts
foreach $c(keys %transcripts){
  @ec=split(/,/,$c);
  $ec=$#ec+1;
  $start=0;
  $end=0;
  @tr=split(/\s/,$transcripts{$c});
  $n=0;
  for($i=0;$i<$#tr;$i+=6){
    $start+=$tr[$i+3];
    $end+=$tr[$i+4];
    $n++;
  }
  $start=int($start/$n);
  $end=int($end/$n);
  print "$tr[1]\tStringTie\ttranscript\t$start\t$end\t.\t$tr[2]\t.\tID=$tr[0];exonCount=$ec;exons=$start-$c-$end;geneID=$tr[5]\n";
}

