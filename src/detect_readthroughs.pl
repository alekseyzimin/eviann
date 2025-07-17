#!/usr/bin/env perl
#NC_015889.1     EviAnn  mRNA    51      3703    .       +       .       ID=XLOC_000001-mRNA-1;Parent=XLOC_000001;EvidenceProteinID=sp|Q6EMT1.1|NU1M_BOSIN:RecName:_Full_NADH-ubiquinone_oxidoreductase_chain_1__AltName:_Full_NADH_dehydrogenase_subunit_1;EvidenceTranscriptID=MSTRG_00000010:4:4880.048828;StartCodon=ATG;StopCodon=TGA;Class=k;Evidence=complete
#NC_015889.1     EviAnn  exon    51      3703    .       +       .       Parent=XLOC_000001-mRNA-1
#NC_015889.1     EviAnn  CDS     2743    3000    .       +       0       Parent=XLOC_000001-mRNA-1

my $gene="";
my $tid="";
my $gid="";
my %startCDS;
my %endCDS;
while($line=<STDIN>){
  chomp($line);
  @gff_fields=split(/\t/,$line);
  if($gff_fields[2] eq "mRNA"|| $gff_fields[2] eq "gene"){
    if($gff_fields[8] =~ /^ID=(\S+);locus=(\S+)$/){
      $tid=$1;
      $gid=$2;
    }
    $endCDS{$tid}=0;
    $endexon{$tid}=0;
    unless($gid eq $gene){
      check_readthroughs(@transcripts) if(scalar(@transcripts)>4);
      @transcripts=();
      $gene=$gid;
    }
    push(@transcripts,$tid);
  }elsif($gff_fields[2] eq "CDS"){
    $startCDS{$tid}=$gff_fields[3] unless(defined($startCDS{$tid}));
    $endCDS{$tid}=$gff_fields[4] if($endCDS{$tid}<$gff_fields[4]);
  }elsif($gff_fields[2] eq "exon"){
    $startexon{$tid}=$gff_fields[3] unless(defined($startexon{$tid}));
    $endexon{$tid}=$gff_fields[4] if($endexon{$tid}<$gff_fields[4]);
  }
}
check_readthroughs(@transcripts) if(scalar(@transcripts)>4);

sub check_readthroughs{
@tr=@_;
#for($i=0;$i<=$#tr;$i++){
#  print "$i exon $startexon{$tr[$i]} $endexon{$tr[$i]} $tr[$i]\n";
#  print "$i CDS $startCDS{$tr[$i]} $endCDS{$tr[$i]} $tr[$i]\n";
#}
for($skip=0;$skip<=$#tr;$skip++){
  $last_CDS_end=$skip==0?$endCDS{$tr[1]}:$endCDS{$tr[0]};
  $last_exon_end=$skip==0?$endexon{$tr[1]}:$endexon{$tr[0]};
  #print "skipping $skip\n";
  for($i=0;$i<=$#tr;$i++){
    next if($i==$skip);
    if($startexon{$tr[$i]}>$last_exon_end && $startCDS{$tr[$i]}>$last_CDS_end){
      #print "readthrough $tr[$skip]\n";
      print "$tr[$skip]\n";
      $i=$#tr;
    }
    $last_CDS_end=$endCDS{$tr[$i]} if($endCDS{$tr[$i]}>$last_CDS_end);
    $last_exon_end=$endexon{$tr[$i]} if($endexon{$tr[$i]}>$last_exon_end);
  }
}
}
     



