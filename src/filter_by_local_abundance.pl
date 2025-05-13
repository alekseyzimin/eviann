#!/usr/bin/env perl
#this code filters the transcripts by local abundance; for the same location (gene) and the same CDS, it keeps the best expressed ones
#input is the GTF file, where the name of the transcript includes the abundance
#transcript_id "MSTRG_00000160:8"; gene_id "XLOC_000001"; xloc "XLOC_000001"; cmp_ref "NP_051101.1.NC_000932.1.81474"; class_code "k"; tss_id "TSS1";
my %class_factor;
my $GFF_K=$ARGV[0];
$class_factor{"="}=20;
$class_factor{"k"}=22;
$class_factor{"j"}=3;
$class_factor{"q"}=2;

#first we read the k gff file -- these are transcritps that yielded complete proteins
open(FILE,$GFF_K);
while($line=<FILE>){
  $full_cds_transcripts{$1}=1 if($line=~/EvidenceTranscriptID=(\S+:\S+:\S+);Start/);
}

my %transcripts_at_xloc_same_cds=();
while(my $line=<STDIN>){
next if($line=~/^#/);
chomp($line);
@gtf_fields=split(/\t/,$line);
if($gtf_fields[2] eq "transcript"){
  #print "DEBUG $line\n";
  my $id=$1 if($gtf_fields[8] =~ /transcript_id \"(\S+)\";/);
  my $xloc=$1 if($gtf_fields[8] =~ /gene_id \"(\S+)\";/);
  my $protein=$1 if($gtf_fields[8] =~ /cmp_ref \"(\S+)\";/);
  my $class=$1 if($gtf_fields[8] =~ /class_code \"(\S+)\";/);
  #print "DEBUG $id $xloc $protein $class $gtf_fields[8]\n";
  if(not($id =~ /^MSTRG/)){#if not a merged transript -- do not mess with it
    print "$id\n";
  }elsif($class =~ /i|y|u|x/ && defined($id)){#no protein, keep
    $transcripts_at_xloc_same_cds{"$xloc:u"}.="$id:$class ";
  }elsif(defined($id) && defined($xloc) && defined($protein) && defined($class)){#protein defined, examine
    $class="q" if($full_cds_transcripts{$id} && not($class=~/k|j|=/));
    #print "DEBUG $id:$class\n";
    $transcripts_at_xloc_same_cds{"$xloc:$protein"}.="$id:$class ";
  }
}
}
    
for $l(keys %transcripts_at_xloc_same_cds){
  my @transcripts=sort by_abundance split(/\s/,$transcripts_at_xloc_same_cds{$l});
  #print "DEBUG sorted $#transcripts transcripts at $l\n",join("\n",@transcripts),"\n";
  my ($tr,$top_count,$top_tpm,$top_class)=split(/:/,$transcripts[0]);
  my $threshold = weight_function($top_count,$top_tpm,$top_class)>=1 ? weight_function($top_count,$top_tpm,$top_class)**$pow : weight_function($top_count,$top_tpm,$top_class);
  my ($xloc,$protein)=split(/:/,$l);
  my $pow=($protein eq "u") ? 0.75 : 0.5;
  my $threshold=weight_function($top_count,$top_tpm,$top_class)**$pow;
  for(my $i=0;$i<=$#transcripts;$i++){
    my ($tr,$count,$tpm,$class)=split(/:/,$transcripts[$i]);
#print "$tr:$count:$tpm $threshold ",weight_function($count,$tpm,$class),"\n" if(weight_function($count,$tpm,$class) >= $threshold);
    if(weight_function($count,$tpm,$class) >= $threshold){
      print "$tr:$count:$tpm\n";
    }else{
      print STDERR "DISCARD $tr:$count:$tpm $threshold ",weight_function($count,$tpm,$class)," $count,$tpm,$class\n"
    }
  }
}


sub by_abundance{
  my @fa=split(/:/,$a);
  my @fb=split(/:/,$b);
  return(weight_function($fb[1],$fb[2],$fb[3]) <=> weight_function($fa[1],$fa[2],$fa[3]));
}
 
sub weight_function{
  my $class_f=0.5;
  $class_f=$class_factor{$_[2]} if(defined($class_factor{$_[2]}));
  #return($_[0]*($_[1])**0.5*$class_f);
  return($_[0]*log($_[1]+1)*$class_f);
  #return($_[0]*$class_f);
}
