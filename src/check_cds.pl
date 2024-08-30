#!/usr/bin/env perl
#This code adds gene and CDS record for every locus in GFF file
my $geneID="";
my @exons=();
my %loci;
my %transcript_seqs=();
my %genome_seqs=();
my %protein_cds;
my %protein;
my %transcript_gff;
my %transcript_cds;
my $protID="";
my $dir="";
my $scf="";
my $seq="";
my %used_proteins;
my $output_prefix=$ARGV[0];
open(OUTFILE1,">$output_prefix".".good_cds.fa.tmp");
open(OUTFILE2,">$output_prefix".".broken_cds.fa.tmp");
open(OUTFILE3,">$output_prefix".".broken_ref.txt.tmp");

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
    if($class_code eq "u"){
      $transcript_u{$geneID}=$line;
      $transcripts_only_loci{$locID}.="$geneID ";
    }else{
      $transcript{$geneID}=$line;
      die("Protein $protID is not defined for protein coding transcript $geneID") if(not(defined($protein{$protID})));
      $transcript_cds{$geneID}=$protID;
      $transcript_cds_start{$geneID}=$protein_start{$protID};
      $transcript_cds_start_codon{$geneID}="INVALID";
      $transcript_cds_end{$geneID}=$protein_end{$protID};
      $transcript_cds_end_codon{$geneID}="INVALID";
      $transcript_class{$geneID}=$class_code;
      $transcript_origin{$geneID}=$gff_fields[1];
      $transcripts_cds_loci{$locID}.="$geneID ";
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
    print "DEBUG start_cds $cds_start_on_transcript end_cds $cds_end_on_transcript transcript length ",length($transcript_seqs{$g}),"\n";
    if(($cds_start_on_transcript < 0 || $cds_start_on_transcript > length($transcript_seqs{$g})) && ($cds_end_on_transcript < 0 || $cds_end_on_transcript > length($transcript_seqs{$g}))){#both start and end are messed up -- send to transdecoder
      print OUTFILE2 ">$g\n$transcript_seqs{$g}\n";
      my @pn=split(/:/,$transcript_cds{$g});
      print OUTFILE3 "$pn[0]\n";
      print "DEBUG broken CDS start and stop outside $g\n";
      next;
    }elsif($cds_start_on_transcript < 0 || $cds_start_on_transcript > length($transcript_seqs{$g})){
      $cds_start_on_transcript=$cds_end_on_transcript%3;
      $cds_length=$cds_end_on_transcript-$cds_start_on_transcript;
    }elsif($cds_end_on_transcript < 0 || $cds_end_on_transcript > length($transcript_seqs{$g})){
      $cds_end_on_transcript=length($transcript_seqs{$g})-(length($transcript_seqs{$g})-$cds_start_on_transcript)%3;
      $cds_length=$cds_end_on_transcript-$cds_start_on_transcript;
    }
    print "DEBUG start_cds $cds_start_on_transcript end_cds $cds_end_on_transcript transcript length ",length($transcript_seqs{$g}),"\n";

    if($cds_length %3 >0){
      print "DEBUG CDS length $cds_length not divisible by 3, possible frameshift, adjusting ";
      if(uc(substr($transcript_seqs{$g},$cds_start_on_transcript,3)) eq "ATG"){
        print "end\n";
        $cds_length-=$cds_length%3;
      }elsif(uc(substr($transcript_seqs{$g},$cds_start_on_transcript+$cds_length,3)) eq "TAG" || uc(substr($transcript_seqs{$g},$cds_start_on_transcript+$cds_length,3)) eq "TAA" || uc(substr($transcript_seqs{$g},$cds_start_on_transcript+$cds_length,3)) eq "TGA"){
        $cds_start_on_transcript+=$cds_length%3;
        $cds_length-=$cds_length%3;
        print "beginning\n"
      }else{
        $cds_length-=$cds_length%3;
        print "end\n";
      }
    }
    $cds_end_on_transcript=$cds_start_on_transcript+$cds_length;

    $first_codon=substr($transcript_seqs{$g},$cds_start_on_transcript,3);
    $last_codon=substr($transcript_seqs{$g},$cds_end_on_transcript,3);
    print "DEBUG $first_codon $last_codon start_cds $cds_start_on_transcript end_cds $cds_end_on_transcript protein $transcript_cds{$g} transcript $g cds_length $cds_length transcript length ",length($transcript_seqs{$g})," tstart $tstart pstart $transcript_cds_start{$g} pend $transcript_cds_end{$g} tori $transcript_ori{$g}\n";

#checking for in-frame stop codons
    my $cds_start_on_transcript_original=$cds_start_on_transcript;
    my $cds_end_on_transcript_original=$cds_end_on_transcript;
    $cds_end_on_transcript=check_in_frame_stops($cds_start_on_transcript_original,$cds_end_on_transcript_original,$transcript_seqs{$g});
    print "DEBUG broken CDS in frame stop $cds_start_on_transcript $cds_end_on_transcript $g\n" if(not($cds_end_on_transcript==$cds_end_on_transcript_original));
#fixing start/stop    
    ($cds_start_on_transcript,$cds_end_on_transcript)=fix_start_stop_codon($cds_start_on_transcript,$cds_end_on_transcript,$transcript_seqs{$g});
    $first_codon=substr($transcript_seqs{$g},$cds_start_on_transcript,3);
    $last_codon=substr($transcript_seqs{$g},$cds_end_on_transcript,3);
    $cds_length=$cds_end_on_transcript-$cds_start_on_transcript;
    if((not(uc($last_codon) eq "TAG" || uc($last_codon) eq "TAA" || uc($last_codon) eq "TGA") || not(uc($first_codon) eq "ATG"))||($cds_end_on_transcript-$cds_start_on_transcript)<0.5*($cds_end_on_transcript_original-$cds_start_on_transcript_original)){
      print OUTFILE2 ">$g\n$transcript_seqs{$g}\n";
      my @pn=split(/:/,$transcript_cds{$g});
      print OUTFILE3 "$pn[0]\n";
      print "DEBUG broken CDS no start and stop or too short $g\n";
    }else{
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
      print "DEBUG $first_codon $last_codon start_cds $cds_start_on_transcript end_cds $cds_end_on_transcript protein $transcript_cds{$g} transcript $g cds_length $cds_length transcript length ",length($transcript_seqs{$g})," tstart $tstart pstart $transcript_cds_start{$g} pend $transcript_cds_end{$g} tori $transcript_ori{$g}\n";
      print OUTFILE1 ">$g $transcript_cds_start{$g} $transcript_cds_end{$g}\n$transcript_seqs{$g}\n";
    }
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
    print "DEBUG start_cds $cds_start_on_transcript end_cds $cds_end_on_transcript transcript length ",length($transcript_seqs{$g}),"\n";
    if(($cds_start_on_transcript < 0 || $cds_start_on_transcript > length($transcript_seqs{$g})) && ($cds_end_on_transcript < 0 || $cds_end_on_transcript > length($transcript_seqs{$g}))){#both start and end are messed up
      print OUTFILE2 ">$g\n$transcript_seqs{$g}\n";
      my @pn=split(/:/,$transcript_cds{$g});
      print OUTFILE3 "$pn[0]\n";
      print "DEBUG broken CDS start and stop outside $g\n";
      next;
    }elsif($cds_start_on_transcript < 0 || $cds_start_on_transcript > length($transcript_seqs{$g})){
      $cds_start_on_transcript=$cds_end_on_transcript%3;
      $cds_length=$cds_end_on_transcript-$cds_start_on_transcript;
    }elsif($cds_end_on_transcript < 0 || $cds_end_on_transcript > length($transcript_seqs{$g})){
      $cds_end_on_transcript=length($transcript_seqs{$g})-(length($transcript_seqs{$g})-$cds_start_on_transcript)%3;
      $cds_length=$cds_end_on_transcript-$cds_start_on_transcript;
    }
    print "DEBUG start_cds $cds_start_on_transcript end_cds $cds_end_on_transcript transcript length ",length($transcript_seqs{$g}),"\n";

    if($cds_length %3 >0){
      print "DEBUG CDS length $cds_length not divisible by 3, possible frameshift, adjusting ";
      if(uc(substr($transcript_seqs{$g},$cds_start_on_transcript,3)) eq "ATG"){
        print "end\n";
        $cds_length-=$cds_length%3;
      }elsif(uc(substr($transcript_seqs{$g},$cds_start_on_transcript+$cds_length,3)) eq "TAG" || uc(substr($transcript_seqs{$g},$cds_start_on_transcript+$cds_length,3)) eq "TAA" || uc(substr($transcript_seqs{$g},$cds_start_on_transcript+$cds_length,3)) eq "TGA"){
        $cds_start_on_transcript+=$cds_length%3;
        $cds_length-=$cds_length%3;
        print "beginning\n"
      }else{
        $cds_length-=$cds_length%3;
        print "end\n";
      }
    }
    $cds_end_on_transcript=$cds_start_on_transcript+$cds_length;

    $first_codon=substr($transcript_seqs{$g},$cds_start_on_transcript,3);
    $last_codon=substr($transcript_seqs{$g},$cds_end_on_transcript,3);
    print "DEBUG $first_codon $last_codon start_cds $cds_start_on_transcript end_cds $cds_end_on_transcript protein $transcript_cds{$g} transcript $g cds_length $cds_length transcript length ",length($transcript_seqs{$g})," tstart $tstart pstart $transcript_cds_start{$g} pend $transcript_cds_end{$g} tori $transcript_ori{$g}\n";

#checking for in-frame stop codons
    my $cds_start_on_transcript_original=$cds_start_on_transcript;
    my $cds_end_on_transcript_original=$cds_end_on_transcript;
    $cds_end_on_transcript=check_in_frame_stops($cds_start_on_transcript_original,$cds_end_on_transcript_original,$transcript_seqs{$g});
    print "DEBUG broken CDS in frame stop $cds_start_on_transcript $cds_end_on_transcript $g\n" if(not($cds_end_on_transcript==$cds_end_on_transcript_original));
#fixing start/stop    
    ($cds_start_on_transcript,$cds_end_on_transcript)=fix_start_stop_codon($cds_start_on_transcript,$cds_end_on_transcript,$transcript_seqs{$g});
    $first_codon=substr($transcript_seqs{$g},$cds_start_on_transcript,3);
    $last_codon=substr($transcript_seqs{$g},$cds_end_on_transcript,3);
    $cds_length=$cds_end_on_transcript-$cds_start_on_transcript;
    if((not(uc($last_codon) eq "TAG" || uc($last_codon) eq "TAA" || uc($last_codon) eq "TGA") || not(uc($first_codon) eq "ATG"))||($cds_end_on_transcript-$cds_start_on_transcript)<0.5*($cds_end_on_transcript_original-$cds_start_on_transcript_original)){
      print OUTFILE2 ">$g\n$transcript_seqs{$g}\n";
      my @pn=split(/:/,$transcript_cds{$g});
      print OUTFILE3 "$pn[0]\n";
      print "DEBUG broken CDS no start and stop or too short $g\n";
    }else{
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
      print "DEBUG $first_codon $last_codon start_cds $cds_start_on_transcript end_cds $cds_end_on_transcript protein $transcript_cds{$g} transcript $g cds_length $cds_length transcript length ",length($transcript_seqs{$g})," tstart $tstart pstart $transcript_cds_start{$g} pend $transcript_cds_end{$g} tori $transcript_ori{$g}\n";
      print OUTFILE1 ">$g $transcript_cds_start{$g} $transcript_cds_end{$g}\n$transcript_seqs{$g}\n";
    }
  }
}

sub fix_start_stop_codon{
  my $cds_start_on_transcript=$_[0];
  my $cds_end_on_transcript=$_[1];
  my $transcript_seq=$_[2];
  my $first_codon=substr($transcript_seq,$cds_start_on_transcript,3);
  my $last_codon=substr($transcript_seq,$cds_end_on_transcript,3);

  if(not(uc($first_codon) eq "ATG")){
    my $i;
    my $found=0;
    for($i=$cds_start_on_transcript-3;$i>=0;$i-=3){
      $found=$i if(uc(substr($transcript_seq,$i,3)) eq "ATG");
      #stop if found a stop
      last if(uc(substr($transcript_seq,$i,3)) eq "TAA" || uc(substr($transcript_seq,$i,3)) eq "TAG" || uc(substr($transcript_seq,$i,3)) eq "TGA");
    } 
    if($found>0){
      print "DEBUG found new start codon upstream at $found\n";
      $cds_start_on_transcript=$found;
    }else{ 
      print "DEBUG failed to find new start codon, looking downstream\n";
      for($i=$cds_start_on_transcript+3;$i<$cds_end_on_transcript;$i+=3){
        last if(uc(substr($transcript_seq,$i,3)) eq "ATG");
      }
      if($i<$cds_end_on_transcript){
        print "DEBUG found new start codon downstream at $i\n";
        $cds_start_on_transcript=$i;
      }else{
        print "DEBUG failed to find new start codon\n";
      }
    }
  }
  if(not(uc($last_codon) eq "TAA" || uc($last_codon) eq "TAG" || uc($last_codon) eq "TGA") && $cds_end_on_transcript<length($transcript_seq)-1){
    my $i;
    for($i=$cds_end_on_transcript+3;$i<length($transcript_seq);$i+=3){
      last if(uc(substr($transcript_seq,$i,3)) eq "TAA" || uc(substr($transcript_seq,$i,3)) eq "TAG" || uc(substr($transcript_seq,$i,3)) eq "TGA");
    } 
    if($i<length($transcript_seq)){
      print "DEBUG found new stop codon downstream at $i\n";
      $cds_end_on_transcript=$i;
    }else{
      print "DEBUG failed to find new stop codon downstream\n";
    }
  }
  return($cds_start_on_transcript,$cds_end_on_transcript);
}

sub check_in_frame_stops{
  my $in_frame_stop=-1;  
  my $cds_start_on_transcript=$_[0];
  my $cds_end_on_transcript=$_[1];
  my $transcript_seq=$_[2];
  print "DEBUG checking for in-frame stops $cds_start_on_transcript $cds_end_on_transcript\n";
  for($i=$cds_start_on_transcript;$i<$cds_end_on_transcript;$i+=3){
    if(uc(substr($transcript_seq,$i,3)) eq "TAA" || uc(substr($transcript_seq,$i,3)) eq "TAG" || uc(substr($transcript_seq,$i,3)) eq "TGA"){
      $in_frame_stop=$i;
      last;
    }
  }
  if($in_frame_stop>-1){
    print "DEBUG found in-frame stop at $in_frame_stop\n";
    $cds_end_on_transcript=$in_frame_stop
  }
  return($cds_end_on_transcript);
}
