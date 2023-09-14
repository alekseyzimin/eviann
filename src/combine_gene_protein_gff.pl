#!/usr/bin/env perl
#This code adds gene and CDS record for every locus in GFF file
my $geneID="";
my @exons=();
my %loci;
my %transcript_seqs=();
my %genome_seqs=();
my %protein_cds;
my %protein;
my @outputLOCbeg;
my @outputLOCend;
my @outputLOCchr;
my %transcript_gff;
my %transcript_cds;
my @gene_records_k;
my %gene_record_k;
my @gene_records_j;
my %gene_record_j;
my @gene_records_u;
my %gene_record_u;
my $protID="";
my $dir="";
my $scf="";
my $seq="";
my %used_proteins;
my %suspect_proteins;
my $output_prefix=$ARGV[0];
open(OUTFILE1,">$output_prefix".".k.gff.tmp");
open(OUTFILE3,">$output_prefix".".u.gff.tmp");
open(OUTFILE4,">$output_prefix".".unused_proteins.gff.tmp");

#this is output of gffcompare -D -o protuniq ../GCF_000001735.4_TAIR10.1_genomic.fna.GCF_000001735.4_TAIR10.1_protein.faa.palign.gff
#here we read in the aligned CDS features
while(my $line=<STDIN>){#we just read in the whole file
  chomp($line);
  my @gff_fields=split(/\t/,$line);
  my @attributes=split(";",$gff_fields[8]);
  if($gff_fields[2] eq "gene"){
    if(not($protID eq "")){
      $protein_cds{$protID}=[@exons];
      $protein_start{$protID}=$pstart;
      $protein_end{$protID}=$pend;
    }
    @exons=();
    $protID=substr($attributes[0],3);#this is protein name
    $pstart=$gff_fields[3];
    $pend=$gff_fields[4];
    $pori=$gff_fields[6];
    die("error in line $line, protein ID $protID already exists in $protein{$protID}") if(defined($protein{$protID}));
    $protein{$protID}=$line;
  }elsif($gff_fields[2] eq "CDS"){
    push(@exons,$line);
  }
}
if(not($protID eq "")){
  $protein_cds{$protID}=[@exons];
  $protein_start{$protID}=$pstart;
  $protein_end{$protID}=$pend;
}
@exons=();

#this is output of gffcompare -r protalign.gtf ../GCF_000001735.4_TAIR10.1_genomic.clean.fna.gtf -o protref
#here we load up all transcripts that matched proteins
open(FILE,$ARGV[1]);
while(my $line=<FILE>){
  chomp($line);
  my @gff_fields=split(/\t/,$line);
  my @attributes=split(";",$gff_fields[8]);
  if($gff_fields[2] eq "transcript"){
    if(defined($transcript{$geneID})){
    #here we need to fix the first and the last exons so that the CDS does not stick out from the boundaries of the transcript
      my @gff_fields=split(/\t/,$exons[0]);
      $gff_fields[3]=$protein_start{$transcript_cds{$geneID}} if($gff_fields[3]>$protein_start{$transcript_cds{$geneID}});
      $gff_fields[3]=1 if($gff_fields[3]<1);
      $exons[0]=join("\t",@gff_fields);
      my @gff_fields=split(/\t/,$exons[-1]);
      $gff_fields[4]=$protein_end{$transcript_cds{$geneID}} if($gff_fields[4]<$protein_end{$transcript_cds{$geneID}});
      $exons[-1]=join("\t",@gff_fields);
      $transcript_gff{$geneID}=[@exons];
      $transcript_ori{$geneID}=$tori;
    }elsif(defined($transcript_u{$geneID})){
      unless($#exons==-1){#not interested in single exon
        $transcript_gff_u{$geneID}=[@exons];
      }
    }
    @exons=();
    $locID=substr($attributes[1],7);#this is the gene_id
    $geneID=substr($attributes[0],3);#this is the transcript_id
    $tstart=$gff_fields[3];
    $tend=$gff_fields[4];
    $tori=$gff_fields[6];
    my $class_code="";
    my $protID="";
    foreach my $attr(@attributes){
      $class_code=substr($attr,11,1) if($attr =~ /^class_code=/);
      $protID=substr($attr,8) if($attr =~ /^cmp_ref=/);
    }
    if($class_code eq "k" || $class_code eq "=" || ($class_code eq "j" && $protein_start{$protID}>=$tstart && $protein_end{$protID}<=$tend)){#equal intron chain or contains protein
      $transcript{$geneID}=$line;
      die("Protein $protID is not defined for protein coding transcript $geneID") if(not(defined($protein{$protID})));
      $transcript_cds{$geneID}=$protID;
      $transcript_cds_start{$geneID}=$protein_start{$protID};
      $transcript_cds_end{$geneID}=$protein_end{$protID};
      $transcript_class{$geneID}=$class_code;
      $transcript_origin{$geneID}=$gff_fields[1];
      $transcripts_cds_loci{$locID}.="$geneID ";
    }elsif($class_code eq "u"){#no match to protein or an inconsistent match; we record these and output them without CDS features only if they are the only ones at a locus
      $transcript_u{$geneID}=$line;
      $transcripts_only_loci{$locID}.="$geneID ";
    }else{#likely messed up protein?
      $suspect_proteins{$protID}=1;
    }
  }elsif($gff_fields[2] eq "exon"){
    push(@exons,$line) if(defined($transcript{$geneID}) || defined($transcript_u{$geneID}));
  }
}
if(defined($transcript{$geneID})){
#here we need to fix the first and the last exons so that the CDS does not stick out from the boundaries of the transcript
  my @gff_fields=split(/\t/,$exons[0]);
  $gff_fields[3]=$protein_start{$transcript_cds{$geneID}} if($gff_fields[3]>$protein_start{$transcript_cds{$geneID}});
  $gff_fields[3]=1 if($gff_fields[3]<1);
  $exons[0]=join("\t",@gff_fields);
  my @gff_fields=split(/\t/,$exons[-1]);
  $gff_fields[4]=$protein_end{$transcript_cds{$geneID}} if($gff_fields[4]<$protein_end{$transcript_cds{$geneID}});
  $exons[-1]=join("\t",@gff_fields);
  $transcript_gff{$geneID}=[@exons];
  $transcript_ori{$geneID}=$tori;
}elsif(defined($transcript_u{$geneID})){
  unless($#exons==-1){#not interested in single exon
    $transcript_gff_u{$geneID}=[@exons];
  }
}
@exons=();

#we load the genome sequences
open(FILE,$ARGV[2]);
while(my $line=<FILE>){
  chomp($line);
  if($line=~ /^>/){
    if(not($scf eq "")){
      $genome_seqs{$scf}=$seq;
      $seq="";
    }
    my @f=split(/\s+/,$line);
    $scf=substr($f[0],1);
  }else{
    $seq.=$line;
  } 
}   
$genome_seqs{$scf}=$seq if(not($scf eq ""));

#we make the transcript sequences for protein coding transcripts
for my $g(keys %transcript_gff){
  $transcript_seqs{$g}="";
  my @gff_fields=();
  for(my $j=0;$j<=$#{$transcript_gff{$g}};$j++){
    @gff_fields=split(/\t/,${$transcript_gff{$g}}[$j]);
    die("Genome sequence $gff_fields[0] needed for transcript $g not found!") if(not(defined($genome_seqs{$gff_fields[0]})));
    $transcript_seqs{$g}.=substr($genome_seqs{$gff_fields[0]},$gff_fields[3]-1,$gff_fields[4]-$gff_fields[3]+1);
  }
  if($gff_fields[6] eq "-"){
    $transcript_seqs{$g}=~tr/ACGTNacgtn/TGCANtgcan/;
    $transcript_seqs{$g}=reverse($transcript_seqs{$g});
  }
}  

#here we try to fix the start codons for the aligned CDS feaures -- if the CDS is assigned to the transcript, we then check for a start codon and, if not found, extend the start/stop of the CDS up to the transcript boundary to look for the valid start/stop
for my $g(keys %transcript_cds){
  my @gff_fields_t=split(/\t/,$transcript{$g});
  my $tstart=$gff_fields_t[3];
  my $tend=$gff_fields_t[4];
  print "\nDEBUG protein $transcript_cds{$g} transript $g length ",length($transcript_seqs{$g}),"\n";
  if($transcript_ori{$g} eq "+"){#forward orientation, check for the start codon
    print "DEBUG examining protein $transcript_cds{$g} $protein_start{$transcript_cds{$g}} $protein_end{$transcript_cds{$g}}\n";
#we need to determine the position of the CDS start on the transcript, minding the introns, and CDS length
    my $cds_start_on_transcript=0;
    my $cds_length=0;
    my @gff_fields=();
    my @gff_fields_prev=();
    for(my $j=0;$j<=$#{$transcript_gff{$g}};$j++){
      @gff_fields=split(/\t/,${$transcript_gff{$g}}[$j]);
      if($transcript_cds_start{$g}>=$gff_fields[3] && $transcript_cds_start{$g}<=$gff_fields[4]){
        $cds_start_on_transcript+=$transcript_cds_start{$g}-$gff_fields[3];
        last;
      }elsif($transcript_cds_start{$g}<$gff_fields[3] && $transcript_cds_start{$g}>$gff_fields_prev[4]){#protein ends inside an intron
        $cds_start_on_transcript+=$transcript_cds_start{$g}-$gff_fields[3]-1;
        last;
      }else{
        $cds_start_on_transcript+=($gff_fields[4]-$gff_fields[3]+1);
      }
      @gff_fields_prev=@gff_fields;
    }
#we need to determine the position of the CDS end on the transcript, minding the introns
    for(my $j=0;$j<=$#{$protein_cds{$transcript_cds{$g}}};$j++){
      @gff_fields_p=split(/\t/,${$protein_cds{$transcript_cds{$g}}}[$j]);
      $cds_length+=$gff_fields_p[4]-$gff_fields_p[3]+1;
    }
    my $cds_end_on_transcript=$cds_start_on_transcript+$cds_length;

#now we look at the start codon
    $first_codon=substr($transcript_seqs{$g},$cds_start_on_transcript,3);
    $last_codon=substr($transcript_seqs{$g},$cds_end_on_transcript,3);
    print "DEBUG $first_codon $last_codon start_cds $cds_start_on_transcript end_cds $cds_end_on_transcript protein $transcript_cds{$g} transcript $g cds_length $cds_length transcript length ",length($transcript_seqs{$g})," tstart $tstart pstart $transcript_cds_start{$g} pend $transcript_cds_end{$g} tori $transcript_ori{$g}\n";

    if($cds_start_on_transcript<0 || $cds_end_on_transcript>length($transcript_seqs{$g})){
      $transcript_class{$g}="n";
      next;
    }

    if(not(uc($first_codon) eq "ATG")){
      my $i;
      for($i=$cds_start_on_transcript-3;$i>=0;$i-=3){
        last if(uc(substr($transcript_seqs{$g},$i,3)) eq "ATG");
      }
      if($i>=0){
        print "DEBUG found new start codon upstream at $i\n";
        $cds_start_on_transcript=$i;
      }else{
        print "DEBUG failed to find new start codon upstream\n";
      }
    }
    if(not(uc($last_codon) eq "TAA" || uc($last_codon) eq "TAG" || uc($last_codon) eq "TGA") && $cds_end_on_transcript<length($transcript_seqs{$g})-1){
      my $i;
      for($i=$cds_end_on_transcript+3;$i<length($transcript_seqs{$g});$i+=3){
        last if(uc(substr($transcript_seqs{$g},$i,3)) eq "TAA" || uc(substr($transcript_seqs{$g},$i,3)) eq "TAG" || uc(substr($transcript_seqs{$g},$i,3)) eq "TGA");
      }
      if($i<length($transcript_seqs{$g})){
        print "DEBUG found new stop codon downstream at $i\n";
        $cds_end_on_transcript=$i;
      }else{
        print "DEBUG failed to find new stop codon downstream\n";
      }
    }

#translating back to genome coords
    my $sequence_covered=0;
    for(my $j=0;$j<=$#{$transcript_gff{$g}};$j++){
      @gff_fields=split(/\t/,${$transcript_gff{$g}}[$j]);
      $sequence_covered+=$gff_fields[4]-$gff_fields[3]+1;
      if($sequence_covered>$cds_start_on_transcript){
        $transcript_cds_start{$g}=$gff_fields[4]-($sequence_covered-$cds_start_on_transcript)+1;
        last;
      }
    }
    my $sequence_covered=0;
    for(my $j=0;$j<=$#{$transcript_gff{$g}};$j++){
      @gff_fields=split(/\t/,${$transcript_gff{$g}}[$j]);
      $sequence_covered+=$gff_fields[4]-$gff_fields[3]+1;
      if($sequence_covered>=$cds_end_on_transcript){
        $transcript_cds_end{$g}=$gff_fields[4]-($sequence_covered-$cds_end_on_transcript);
        last;
      }
    }

    $first_codon=substr($transcript_seqs{$g},$cds_start_on_transcript,3);
    $last_codon=substr($transcript_seqs{$g},$cds_end_on_transcript,3);
    $cds_length=$cds_end_on_transcript-$cds_start_on_transcript+1;
    print "DEBUG $first_codon $last_codon start_cds $cds_start_on_transcript end_cds $cds_end_on_transcript protein $transcript_cds{$g} transcript $g cds_length $cds_length transcript length ",length($transcript_seqs{$g})," tstart $tstart pstart $transcript_cds_start{$g} pend $transcript_cds_end{$g} tori $transcript_ori{$g}\n";

  }else{#reverse orientation
    print "DEBUG examining protein $transcript_cds{$g} $protein_start{$transcript_cds{$g}} $protein_end{$transcript_cds{$g}}\n";
#we need to determine the position or the CDS start on the transcript, minding the introns
    my $cds_start_on_transcript=0;
    my $cds_length=0;
    my @gff_fields=();
    my @gff_fields_prev=();
    for(my $j=$#{$transcript_gff{$g}};$j>=0;$j--){
      @gff_fields=split(/\t/,${$transcript_gff{$g}}[$j]);
      if($transcript_cds_end{$g}>=$gff_fields[3] && $transcript_cds_end{$g}<=$gff_fields[4]){
        $cds_start_on_transcript+=$gff_fields[4]-$transcript_cds_end{$g};
        last;
      }elsif($transcript_cds_end{$g}>$gff_fields[4] && $transcript_cds_end{$g}<$gff_fields_prev[3]){
        $cds_start_on_transcript+=$gff_fields[4]-$transcript_cds_end{$g}-1;
        last;
      }else{
        $cds_start_on_transcript+=($gff_fields[4]-$gff_fields[3]+1);
      }
      @gff_fields_prev=@gff_fields;
    }
    #we need to determine the position of the CDS end on the transcript, minding the introns
    for(my $j=0;$j<=$#{$protein_cds{$transcript_cds{$g}}};$j++){
      @gff_fields_p=split(/\t/,${$protein_cds{$transcript_cds{$g}}}[$j]);
      $cds_length+=$gff_fields_p[4]-$gff_fields_p[3]+1;
    }
    my $cds_end_on_transcript=$cds_start_on_transcript+$cds_length;
    $first_codon=substr($transcript_seqs{$g},$cds_start_on_transcript,3);
    $last_codon=substr($transcript_seqs{$g},$cds_end_on_transcript,3);
    print "DEBUG $first_codon $last_codon start_cds $cds_start_on_transcript end_cds $cds_end_on_transcript protein $transcript_cds{$g} transcript $g cds_length $cds_length transcript length ",length($transcript_seqs{$g})," tstart $tstart pstart $transcript_cds_start{$g} pend $transcript_cds_end{$g} tori $transcript_ori{$g}\n";

    if($cds_start_on_transcript<0 || $cds_end_on_transcript>length($transcript_seqs{$g})){
      $transcript_class{$g}="n";
      next;
    }

    if(not(uc($first_codon) eq "ATG")){
      my $i;
      for($i=$cds_start_on_transcript-3;$i>=0;$i-=3){
        last if(uc(substr($transcript_seqs{$g},$i,3)) eq "ATG");
      }
      if($i>=0){
        print "DEBUG found new start codon upstream at $i\n";
        $cds_start_on_transcript=$i;
      }else{
        print "DEBUG failed to find new start codon upstream\n";
      }
    }
    if(not(uc($last_codon) eq "TAA" || uc($last_codon) eq "TAG" || uc($last_codon) eq "TGA") && $cds_end_on_transcript<length($transcript_seqs{$g})-1){
      my $i;
      for($i=$cds_end_on_transcript+3;$i<length($transcript_seqs{$g});$i+=3){
        last if(uc(substr($transcript_seqs{$g},$i,3)) eq "TAA" || uc(substr($transcript_seqs{$g},$i,3)) eq "TAG" || uc(substr($transcript_seqs{$g},$i,3)) eq "TGA");
      } 
      if($i<length($transcript_seqs{$g})){
        print "DEBUG found new stop codon downstream at $i\n";
        $cds_end_on_transcript=$i;
      }else{
        print "DEBUG failed to find new stop codon downstream\n";
      } 
    } 

#translating start and end to genome coords
    my $sequence_covered=0;
    for(my $j=$#{$transcript_gff{$g}};$j>=0;$j--){
      @gff_fields=split(/\t/,${$transcript_gff{$g}}[$j]);
      $sequence_covered+=$gff_fields[4]-$gff_fields[3]+1;
      if($sequence_covered>$cds_start_on_transcript){
        $transcript_cds_end{$g}=$gff_fields[3]+($sequence_covered-$cds_start_on_transcript)-1;
        last;
      }
    }

    my $sequence_covered=0;
    for(my $j=$#{$transcript_gff{$g}};$j>=0;$j--){ 
      @gff_fields=split(/\t/,${$transcript_gff{$g}}[$j]);
      $sequence_covered+=$gff_fields[4]-$gff_fields[3]+1;
      if($sequence_covered>=$cds_end_on_transcript){
        $transcript_cds_start{$g}=$gff_fields[3]+($sequence_covered-$cds_end_on_transcript);
        last;
      }
    }

    $first_codon=substr($transcript_seqs{$g},$cds_start_on_transcript,3);
    $last_codon=substr($transcript_seqs{$g},$cds_end_on_transcript,3);
    $cds_length=$cds_end_on_transcript-$cds_start_on_transcript+1;
    print "DEBUG $first_codon $last_codon start_cds $cds_start_on_transcript end_cds $cds_end_on_transcript protein $transcript_cds{$g} transcript $g cds_length $cds_length transcript length ",length($transcript_seqs{$g})," tstart $tstart pstart $transcript_cds_start{$g} pend $transcript_cds_end{$g} tori $transcript_ori{$g}\n";
  }
}

#process the loci
#print the headers
print OUTFILE1 "##gff-version 3\n# EviAnn automated annotation\n";
for my $locus(keys %transcripts_cds_loci){
  my @output=();
  my $gene_feature="";
  my @transcripts_at_loci=split(/\s+/,$transcripts_cds_loci{$locus});
  my @gff_fields=split(/\t/,$transcript{$transcripts_at_loci[0]});
  my $locus_start=$gff_fields[3];
  my $locus_end=$gff_fields[4];
  my @attributes=split(";",$gff_fields[8]);
  $geneID=substr($attributes[1],7);#this is the XLOC
  my $parent=$geneID."-mRNA-";
  my $transcript_index=0;
  my %output_proteins_for_locus=();
  #we output transcripts by class code, first = then k and then j, and we record which cds we used; if the cds was used for a higher class we skip the transcript
  for my $class ("=","k","j"){
    my $class_success=0;
    for my $t(@transcripts_at_loci){
      $class_success=1 if($class eq "=" || $class eq "k");
      next if(not($transcript_class{$t} eq $class));
      #next if($class eq "j" && $class_success);#not interested in outputting j's if already have = or k here
      my $protID=$transcript_cds{$t};
      #next if(defined($output_proteins_for_locus{$protID}) && $class eq "j");#do not need to output j transcripts matching proteins already output earlier
      $output_proteins_for_locus{$protID}=1;
      $used_proteins{$protID}=1;
      my $note="";
      my @gff_fields_t=split(/\t/,$transcript{$t});
      my @attributes_t=split(";",$gff_fields_t[8]);
      my $start_cds=$transcript_cds_start{$t};
      my $end_cds=$transcript_cds_end{$t};;
      my $transcript_start=$gff_fields_t[3];
      my $transcript_end=$gff_fields_t[4];
      #next if($class eq "j" &&  $end_cds>$transcript_end);#we are not interested in j's that are too short
      my $transcript_cds_start_index=0;
      my $transcript_cds_end_index=$#{$transcript_gff{$t}};
      $locus_start=$gff_fields_t[3] if($gff_fields_t[3]<$locus_start);
      $locus_end=$gff_fields_t[4] if($gff_fields_t[4]>$locus_end);
#here we figure out which exons are UTR
      for(my $i=0;$i<=$#{$transcript_gff{$t}};$i++){
        my @gff_fields=split(/\t/,${$transcript_gff{$t}}[$i]);
        if($gff_fields[4]>=$start_cds && $gff_fields[3]<=$start_cds){
          $transcript_cds_start_index=$i;
          last;
        }
      }
      for(my $i=$#{$transcript_gff{$t}};$i>=0;$i--){
        my @gff_fields=split(/\t/,${$transcript_gff{$t}}[$i]);
        if($gff_fields[3]<=$end_cds && $gff_fields[4]>=$end_cds){
          $transcript_cds_end_index=$i;
          last;
        }
      }
      $locus_start=$transcript_start if($transcript_start<$locus_start);
      $locus_end=$transcript_end if($transcript_end>$locus_end);
#output transcript
      $transcript_index++;
      push(@output,$gff_fields[0]."\tEviAnn\tmRNA\t$transcript_start\t$transcript_end\t".join("\t",@gff_fields_t[5..7])."\tID=$parent$transcript_index;Parent=$geneID;ProteinID=$protID");
#output exons
      my $i=1;
      my $first_j=0;
      my $last_j=$#{$transcript_gff{$t}};
      for (my $j=$first_j;$j<=$last_j;$j++){
        my @gff_fields=split(/\t/,${$transcript_gff{$t}}[$j]);
        if($first_j==$last_j){ #single exon, make sure the boundaries are OK
          my $exon_start=$gff_fields[3];
          my $exon_end=$gff_fields[4];
          $exon_start=$start_cds if($exon_start>$start_cds);
          $exon_end=$end_cds if($exon_end<$end_cds);
          push(@output,$gff_fields[0]."\tEviAnn\texon\t$exon_start\t$exon_end\t".join("\t",@gff_fields[5..7])."\tID=$parent$transcript_index:exon:$i;Parent=$parent$transcript_index");
        }else{
          push(@output,$gff_fields[0]."\tEviAnn\texon\t".join("\t",@gff_fields[3..7])."\tID=$parent$transcript_index:exon:$i;Parent=$parent$transcript_index");
        }
        $i++;
      }
      $i=1;
#output 5'UTR 
      for(my $j=$first_j;$j<=$transcript_cds_start_index;$j++){
        my @gff_fields=split(/\t/,${$transcript_gff{$t}}[$j]);
        if($j==$transcript_cds_start_index){
          push(@output,$gff_fields[0]."\tEviAnn\tfive_prime_UTR\t$gff_fields[3]\t".($start_cds-1)."\t".join("\t",@gff_fields[5..7])."\tID=$parent$transcript_index:5UTR:$i;Parent=$parent$transcript_index") if($start_cds>$gff_fields[3]);
        }else{
          push(@output,$gff_fields[0]."\tEviAnn\tfive_prime_UTR\t".join("\t",@gff_fields[3..7])."\tID=$parent$transcript_index:5UTR:$i;Parent=$parent$transcript_index");
        }
        $i++;
      }
#output cds from exons
      $i=1;
      if($transcript_cds_start_index<$transcript_cds_end_index){
        my @gff_fields=split(/\t/,${$transcript_gff{$t}}[$transcript_cds_start_index]);
        push(@output,$gff_fields[0]."\tEviAnn\tcds\t$start_cds\t".join("\t",@gff_fields[4..7])."\tID=$parent$transcript_index:cds:$i;Parent=$parent$transcript_index$note");
        $i++;
        for(my $j=$transcript_cds_start_index+1;$j<$transcript_cds_end_index;$j++){
          my @gff_fields=split(/\t/,${$transcript_gff{$t}}[$j]);
          push(@output,$gff_fields[0]."\tEviAnn\tcds\t".join("\t",@gff_fields[3..7])."\tID=$parent$transcript_index:cds:$i;Parent=$parent$transcript_index$note");
          $i++;
        }
        my @gff_fields=split(/\t/,${$transcript_gff{$t}}[$transcript_cds_end_index]);
        push(@output,$gff_fields[0]."\tEviAnn\tcds\t$gff_fields[3]\t$end_cds\t".join("\t",@gff_fields[5..7])."\tID=$parent$transcript_index:cds:$i;Parent=$parent$transcript_index$note");
      }else{#single exon
        my @gff_fields=split(/\t/,${$transcript_gff{$t}}[$transcript_cds_start_index]);
        push(@output,$gff_fields[0]."\tEviAnn\tcds\t$start_cds\t$end_cds\t".join("\t",@gff_fields[5..7])."\tID=$parent$transcript_index:cds:$i;Parent=$parent$transcript_index$note");
      }
#output 3'UTR
      $i=1;
      for(my $j=$transcript_cds_end_index;$j<=$last_j;$j++){
        my @gff_fields=split(/\t/,${$transcript_gff{$t}}[$j]);
        if($j==$transcript_cds_end_index){
          push(@output,$gff_fields[0]."\tEviAnn\tthree_prime_UTR\t".($end_cds+1)."\t".join("\t",@gff_fields[4..7])."\tID=$parent$transcript_index:3UTR:$i;Parent=$parent$transcript_index") if($end_cds<$gff_fields[4]);
        }else{
          push(@output,$gff_fields[0]."\tEviAnn\tthree_prime_UTR\t".join("\t",@gff_fields[3..7])."\tID=$parent$transcript_index:3UTR:$i;Parent=$parent$transcript_index");
        }
        $i++
      }
    }#end of transcripts loop
  }#end of class loop
#now we know the locus start ans end, and we can output the gene record
  $gene_record_k{$gff_fields[0]." ".$locus_start}="$gff_fields[0]\tEviAnn\tgene\t$locus_start\t$locus_end\t".join("\t",@gff_fields[5..7])."\tID=$geneID;geneID=$geneID;type=protein_coding\n".join("\n",@output)."\n";
  push(@gene_records_k,$gff_fields[0]." ".$locus_start);
  push(@outputLOCchr,$gff_fields[0]);
  push(@outputLOCbeg,$locus_start);
  push(@outputLOCend,$locus_end);
}

#finally output "intergenic" transcripts
#some of these may be completely messed up
for my $locus(keys %transcripts_only_loci){
  #next if we have seen this locus with a protein match
  next if(defined($transcripts_cds_loci{$locus}));
  my @output=();
  my @transcripts_at_loci=split(/\s+/,$transcripts_only_loci{$locus});
  my @gff_fields=split(/\t/,$transcript_u{$transcripts_at_loci[0]});
  my $locus_start=$gff_fields[3];
  my $locus_end=$gff_fields[4];
  my $geneID=$locus."U_lncRNA";
  my $parent=$geneID."-mRNA-";
  my $transcript_index=0;
  #here we first compute intron junction score as the number of distinct intron junctions over the number of total intron junctions
  my %distinct_intron_junctions=();
  my $total_intron_junctions=1;
  my $junction_score=0;
  if($#transcripts_at_loci>0){
    for my $t(@transcripts_at_loci){
      next if(not(defined($transcript_gff_u{$t})));
      my @gff_fields_t=split(/\t/,$transcript_u{$t});
      my @attributes_t=split(";",$gff_fields_t[8]);
      $locus_start=$gff_fields_t[3] if($gff_fields_t[3]<$locus_start);
      $locus_end=$gff_fields_t[4] if($gff_fields_t[4]>$locus_end);
      for(my $j=1;$j<=$#{$transcript_gff_u{$t}};$j++){
        my @gff_fields_curr=split(/\t/,${$transcript_gff_u{$t}}[$j]);
        my @gff_fields_prev=split(/\t/,${$transcript_gff_u{$t}}[$j-1]);
        $distinct_intron_junctions{"$gff_fields_prev[4] $gff_fields_curr[3]"}=1;
        $total_intron_junctions++;
     }
    }
    $junction_score=scalar(keys %distinct_intron_junctions)/$total_intron_junctions;
  }
  #do not output the locus if there are too many disagreements between the intron junctions
  next if($junction_score>0.66);
  #if we got here we can output the transcript
  for my $t(@transcripts_at_loci){
    next if(not(defined($transcript_gff_u{$t})));
    my @gff_fields_t=split(/\t/,$transcript_u{$t});
    my @attributes_t=split(";",$gff_fields_t[8]);
    $transcript_index++;
    push(@output,"$gff_fields_t[0]\tEviAnn\tmRNA\t".join("\t",@gff_fields_t[3..7])."\tID=$parent$transcript_index;Parent=$geneID;$attributes_t[3]");
    my $i=1;
    for my $x(@{$transcript_gff_u{$t}}){
      my @gff_fields=split(/\t/,$x);
      push(@output,"$gff_fields_t[0]\tEviAnn\t".join("\t",@gff_fields[2..7])."\tID=$parent$transcript_index:exon:$i;Parent=$parent$transcript_index");
      $i++;
    }
  }
  if($transcript_index>0){
    $gene_record_u{$gff_fields[0]." ".$locus_start}="$gff_fields[0]\tEviAnn\tgene\t$locus_start\t$locus_end\t".join("\t",@gff_fields[5..7])."\tID=$geneID;geneID=$geneID;type=lncRNA;junction_score=$junction_score;\n".join("\n",@output)."\n";
    push(@gene_records_u,$gff_fields[0]." ".$locus_start);
  }
}

#output unused proteins
#we will then look at them, pick only one per locus that best matches uniprot and join them in at the second pass
foreach my $p(keys %protein){
  next if(defined($used_proteins{$p}));
  #next if(defined($suspect_proteins{$p}));
  my @gff_fields_p=split(/\t/,$protein{$p});
  my $ptstart=$gff_fields_p[3]-99>0 ? $gff_fields_p[3]-99:1;
  my $ptend=$gff_fields_p[4]+99<=length($genome_seqs{$gff_fields_p[0]}) ? $gff_fields_p[4]+99:length($genome_seqs{$gff_fields_p[0]});

  print OUTFILE4 "$gff_fields_p[0]\tEviAnn\t$gff_fields_p[2]\t",$ptstart,"\t",$ptend,"\t",join("\t",@gff_fields_p[5..$#gff_fields_p]),"\n";
  for(my $j=0;$j<=$#{$protein_cds{$p}};$j++){
    my @gff_fields_c=split(/\t/,${$protein_cds{$p}}[$j]);
    if($#{$protein_cds{$p}}==0){#add "fake" 5' utr and 3' utr
      print OUTFILE4 "$gff_fields_c[0]\tEviAnn\texon\t$ptstart\t$ptend\t",join("\t",@gff_fields_c[5..$#gff_fields_c]),"\n";
    }elsif($j==0){#add "fake" 5' utr
      print OUTFILE4 "$gff_fields_c[0]\tEviAnn\texon\t$ptstart\t",join("\t",@gff_fields_c[4..$#gff_fields_c]),"\n";
    }elsif($j==$#{$protein_cds{$p}}){
      print OUTFILE4 "$gff_fields_c[0]\tEviAnn\texon\t$gff_fields_c[3]\t$ptend\t",join("\t",@gff_fields_c[5..$#gff_fields_c]),"\n";
    }else{
      print OUTFILE4 "$gff_fields_c[0]\tEviAnn\texon\t",join("\t",@gff_fields_c[3..$#gff_fields_c]),"\n";
    }
  }
  for(my $j=0;$j<=$#{$protein_cds{$p}};$j++){
    my @gff_fields_c=split(/\t/,${$protein_cds{$p}}[$j]);
    print OUTFILE4 "$gff_fields_c[0]\tEviAnn\tcds\t",join("\t",@gff_fields_c[3..$#gff_fields_c]),"\n";
  } 
}

#now we sort and output
my @gene_records_sorted=sort mysort @gene_records_k;
my %output=();
foreach $g(@gene_records_sorted){
  if(not(defined($output{$g}))){
    print OUTFILE1 $gene_record_k{$g};
    $output{$g}=1;
  }
}

@gene_records_sorted=sort mysort @gene_records_j;
%output=();
foreach $g(@gene_records_sorted){
  if(not(defined($output{$g}))){
    print OUTFILE2 $gene_record_j{$g};
    $output{$g}=1;
  }
}

@gene_records_sorted=sort mysort @gene_records_u;
%output=();
foreach $g(@gene_records_sorted){
  if(not(defined($output{$g}))){
    print OUTFILE3 $gene_record_u{$g};
    $output{$g}=1;
  }
}

sub mysort{
my ($chroma,$coorda)=split(/\s+/,$a);
my ($chromb,$coordb)=split(/\s+/,$b);
return($chroma cmp $chromb || $coorda <=>$coordb);
}


  
