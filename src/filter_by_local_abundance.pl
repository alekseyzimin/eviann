#!/usr/bin/env perl
#this code filters the transcripts by local abundance; for the sale location (gene) and the same CDS, it keeps the best expressed ones
#input is the GTF file, where the name of the transcript includes the abundance
#transcript_id "MSTRG_00000160:8"; gene_id "XLOC_000001"; xloc "XLOC_000001"; cmp_ref "NP_051101.1.NC_000932.1.81474"; class_code "k"; tss_id "TSS1";
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
  if(not($id =~ /^MSTRG/)){#if not a merged transript -- do not mess with it
    print "$id\n";
  }elsif($class eq "u" && defined($id)){#no protein, keep
    $transcripts_at_xloc_same_cds{"$xloc:u"}.="$id:$class ";
  }elsif(defined($id) && defined($xloc) && defined($protein) && ($class eq "j" || $class eq "k" || $class eq "=")){#protein defined, examine
    $transcripts_at_xloc_same_cds{"$xloc:$protein"}.="$id:$class ";
    #print "DEBUG $xloc:$protein $id\n";
  }
}
}
    
for $l(keys %transcripts_at_xloc_same_cds){
  #print "$l\n";
  my @transcripts=sort by_abundance split(/\s/,$transcripts_at_xloc_same_cds{$l});
  my ($tr,$top_count,$top_tpm,$top_class)=split(/:/,$transcripts[0]);
  my $threshold=($top_count*$top_tpm)**.5;
  #print "DEBUG top $tr count $top_count class $class threshold $threshold\n";
  for(my $i=0;$i<=$#transcripts;$i++){
    my ($tr,$count,$tpm,$class)=split(/:/,$transcripts[$i]);
    print "$tr:$count:$tpm\n" if($count*$tpm >= $threshold || ($top_class eq "j" && ($class eq "k" || $class eq "=")));
  }
}


sub by_abundance{
my @fa=split(/:/,$a);
my @fb=split(/:/,$b);
return($fb[1]*$fb[2] <=> $fa[1]*$fa[2]);
}
  
