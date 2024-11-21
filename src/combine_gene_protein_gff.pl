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
my @outputLOCdir;
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
my $ext_length=40;
my $length_fraction=0.75;
my $output_prefix=$ARGV[0];
#these are genetic codes for HMMs
my %code=();
$code{"A"}=0;
$code{"a"}=0;
$code{"C"}=1;
$code{"c"}=1;
$code{"G"}=2;
$code{"g"}=2;
$code{"T"}=3;
$code{"t"}=3;
$code{"N"}=4;
$code{"n"}=4;
open(OUTFILE1,">$output_prefix".".k.gff.tmp");
open(OUTFILE3,">$output_prefix".".u.gff.tmp");
open(OUTFILE4,">$output_prefix".".unused_proteins.gff.tmp");
my $donor_length=16;
my $acceptor_length=30;

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
  }elsif(uc($gff_fields[2]) eq "CDS"){
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
    if($class_code eq "u"){#no match to protein or an inconsistent match; we record these and output them without CDS features only if they are the only ones at a locus
      $transcript_u{$geneID}=$line;
      $transcripts_only_loci{$locID}.="$geneID ";
    }else{
      $transcript{$geneID}=$line;
      die("Protein $protID is not defined for protein coding transcript $geneID") if(not(defined($protein{$protID})));
      $transcript_cds{$geneID}=$protID;
      $transcript_source{$geneID}=$gff_fields[1];
      $transcript_cds_start{$geneID}=$protein_start{$protID};
      $transcript_cds_start_codon{$geneID}="MISSING";
      $transcript_cds_end{$geneID}=$protein_end{$protID};
      $transcript_cds_end_codon{$geneID}="MISSING";
      $transcript_class{$geneID}=$class_code;
      $transcript_origin{$geneID}=$gff_fields[1];
      $transcripts_cds_loci{$locID}.="$geneID ";
    }
    $transcript_class{$geneID}="NA" if($gff_fields[8] =~ /contained_in=/);
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
print "DEBUG Loading genome sequence\n";
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

#we load the adjusted cds start and stop coordinates
open(FILE,$ARGV[3]);
  print "DEBUG Loading CDS start/stop coordinates from TransDecoder\n";
while(my $line=<FILE>){
  chomp($line);
  my ($geneID, $cds_start,$cds_stop)=split(/\s+/,$line);
  $transcript_cds_start_on_transcript{$geneID}=$cds_start;
  $transcript_cds_stop_on_transcript{$geneID}=$cds_stop-3;
}

#we load SNAP HMMs
if(defined($ARGV[4])){
  print "DEBUG Loading PWMs\n";
  open(FILE,$ARGV[4]);
  $line=<FILE>;
  if($line =~ /^zoeHMM/){#check format
    while($line=<FILE>){
      chomp($line);
      if($line=~/^Donor WMM/){
        my $i=0;
        while($line=<FILE>){
          last if($line=~/NN TRM/);
          chomp($line);
          $line=~s/^\s+//;
          my @f=split(/\s+/,$line);
          for(my $j=0;$j<=3;$j++){
            $donor_freq[$i][$j]=$f[$j];
          }
          $donor_freq[$i][4]=0;
          $i++;
        }
      }elsif($line=~/^Acceptor WMM/){
        my $i=0;
        while($line=<FILE>){
          last if($line=~/NN TRM/);
          chomp($line);
          $line=~s/^\s+//;
          my @f=split(/\s+/,$line);
          for(my $j=0;$j<=3;$j++){
            $acceptor_freq[$i][$j]=$f[$j];
          }
          $acceptor_freq[$i][4]=0;
          $i++;
        }
      }elsif($line=~/^Start/){
        my $i=0;
        while($line=<FILE>){
          last if($line=~/NNN TRM/);
          chomp($line);
          $line=~s/^\s+//;
          my @f=split(/\s+/,$line);
          for(my $j=0;$j<=3;$j++){
            $coding_start_freq[$i][$j]=$f[$j];
          }
          $coding_start_freq[$i][4]=0;
          $i++;
        }
      }
    }
  }
}else{
  #no extension
  $ext_length=0;
}

#finally, if available we load the original transcript names
my %original_transcript_name=();
if(defined($ARGV[5])){
  open(FILE,$ARGV[5]);
  while(my $line=<FILE>){
    chomp($line);
    @f=split(/\s+/,$line);
    $original_transcript_name{$f[0]}=$f[1];
  }
}

#we make the transcript sequences for protein coding transcripts and score the transcripts with HMMs
for my $g(keys %transcript_gff){
  $transcript_seqs{$g}="";
  my @gff_fields=();
  for(my $j=0;$j<=$#{$transcript_gff{$g}};$j++){
    @gff_fields=split(/\t/,${$transcript_gff{$g}}[$j]);
    die("Genome sequence $gff_fields[0] needed for transcript $g not found!") if(not(defined($genome_seqs{$gff_fields[0]})));
    if($j==0){
      my $ext="";
      $ext=substr($genome_seqs{$gff_fields[0]},$gff_fields[3]-$ext_length-1,$ext_length) if($gff_fields[3]>=$ext_length+1);
      if($gff_fields[6] eq "+"){
        $transcript_seqs_5pext{$g}=$ext;
      }else{
        $transcript_seqs_3pext{$g}=$ext;
      }
    }
    if($j==$#{$transcript_gff{$g}}){
      my $ext="";
      $ext=substr($genome_seqs{$gff_fields[0]},$gff_fields[4],$ext_length) if($gff_fields[4]<length($genome_seqs{$gff_fields[0]})-$ext_length);
      if($gff_fields[6] eq "+"){
        $transcript_seqs_3pext{$g}=$ext;
      }else{
        $transcript_seqs_5pext{$g}=$ext;
      }
    } 
    $transcript_seqs{$g}.=substr($genome_seqs{$gff_fields[0]},$gff_fields[3]-1,$gff_fields[4]-$gff_fields[3]+1);
  }

  $transcript_seqs_5pext_actual{$g}="";
  $transcript_seqs_3pext_actual{$g}="";
  
  if($gff_fields[6] eq "-"){
    $transcript_seqs{$g}=~tr/ACGTNacgtn/TGCANtgcan/;
    $transcript_seqs{$g}=reverse($transcript_seqs{$g});
    $transcript_seqs_5pext{$g}=~tr/ACGTNacgtn/TGCANtgcan/;
    $transcript_seqs_5pext{$g}=reverse($transcript_seqs_5pext{$g});
    $transcript_seqs_3pext{$g}=~tr/ACGTNacgtn/TGCANtgcan/;
    $transcript_seqs_3pext{$g}=reverse($transcript_seqs_3pext{$g});
  }
}  

#here we try to fix the start codons for the aligned CDS feaures -- if the CDS is assigned to the transcript, we then check for a start codon and, if not found, extend the start/stop of the CDS up to the transcript boundary to look for the valid start/stop
for my $g(keys %transcript_cds){
  next if($transcript_class{$g} eq "NA");
  my @gff_fields_t=split(/\t/,$transcript{$g});
  my $tstart=$gff_fields_t[3];
  my $tend=$gff_fields_t[4];
  print "\nDEBUG protein $transcript_cds{$g} transcript $g $original_transcript_name{$g} length ",length($transcript_seqs{$g}),"\n";
  if($transcript_ori{$g} eq "+"){#forward orientation, check for the start codon
    print "DEBUG examining protein $transcript_cds{$g} $protein_start{$transcript_cds{$g}} $protein_end{$transcript_cds{$g}} from transdecoder $transcript_cds_start_on_transcript{$g} $transcript_cds_stop_on_transcript{$g}\n";
#we need to determine the position of the CDS start on the transcript, minding the introns, and CDS length
    my $cds_start_on_transcript=0;
    my $cds_end_on_transcript=0;
    my $cds_length=0;
    if(defined($transcript_cds_start_on_transcript{$g}) && defined($transcript_cds_stop_on_transcript{$g})){
      $cds_start_on_transcript=$transcript_cds_start_on_transcript{$g};
      $cds_end_on_transcript=$transcript_cds_stop_on_transcript{$g};
      $cds_length=$cds_end_on_transcript-$cds_start_on_transcript;
    }else{
      print "DEBUG computing CDS start and stop position on the transcript\n";
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
      $cds_end_on_transcript=$cds_start_on_transcript+$cds_length;
    }

    print "DEBUG start_cds $cds_start_on_transcript end_cds $cds_end_on_transcript transcript length ",length($transcript_seqs{$g}),"\n";

    if(($cds_start_on_transcript < 0 || $cds_start_on_transcript > length($transcript_seqs{$g})) && ($cds_end_on_transcript < 0 || $cds_end_on_transcript > length($transcript_seqs{$g}))){#both start and end are messed up
      $transcript_class{$g}="NA";
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

    my $cds_start_on_transcript_original=$cds_start_on_transcript;
    my $cds_end_on_transcript_original=$cds_end_on_transcript;
    $first_codon=substr($transcript_seqs{$g},$cds_start_on_transcript,3);
    $last_codon=substr($transcript_seqs{$g},$cds_end_on_transcript,3);
    print "DEBUG $first_codon $last_codon start_cds $cds_start_on_transcript end_cds $cds_end_on_transcript protein $transcript_cds{$g} transcript $g cds_length $cds_length transcript length ",length($transcript_seqs{$g})," tstart $tstart pstart $transcript_cds_start{$g} pend $transcript_cds_end{$g} tori $transcript_ori{$g}\n";

#checking for in-frame stop codons
    ($cds_start_on_transcript,$cds_end_on_transcript)=fix_in_frame_stops_keep_frame($cds_start_on_transcript_original,$cds_end_on_transcript_original,$transcript_seqs{$g});

#fixing start/stop, extending if needed
    ($cds_start_on_transcript,$cds_end_on_transcript)=fix_start_stop_codon($cds_start_on_transcript,$cds_end_on_transcript,$transcript_seqs{$g});
    $first_codon=substr($transcript_seqs{$g},$cds_start_on_transcript,3);
    $last_codon=substr($transcript_seqs{$g},$cds_end_on_transcript,3);
    ($cds_start_on_transcript,$cds_end_on_transcript,$transcript_seqs_5pext_actual{$g},$transcript_seqs_3pext_actual{$g})=fix_start_stop_codon_ext($cds_start_on_transcript,$cds_end_on_transcript,$transcript_seqs{$g},$transcript_seqs_5pext{$g},$transcript_seqs_3pext{$g}) if((length($transcript_seqs_5pext{$g})==$ext_length && length($transcript_seqs_3pext{$g})==$ext_length) && (not(uc($first_codon) eq "ATG") || not(uc($last_codon) eq "TAA" || uc($last_codon) eq "TAG" || uc($last_codon) eq "TGA"))); #try to extend
    #at this point we did all we could to fix the coding region; last thing that remains is to make sure the stop codon is included both into the cds and the transcript
    #from here on we assume that the stop codon is included!!!
    $cds_end_on_transcript+=3;
    #make sure we are not out of the transcript!!!
    $transcript_seqs_3pext_actual{$g}=substr($transcript_seqs_3pext{$g},0,length($transcript_seqs_3pext_actual{$g})+3) if($cds_end_on_transcript > length($transcript_seqs_5pext_actual{$g}.$transcript_seqs{$g}.$transcript_seqs_3pext_actual{$g}));

    if(length($transcript_seqs_5pext_actual{$g})>0 || length($transcript_seqs_3pext_actual{$g})>0){
#updating the transcript
      $transcript_seqs{$g}=$transcript_seqs_5pext_actual{$g}.$transcript_seqs{$g}.$transcript_seqs_3pext_actual{$g};

#updating the transcript records
      @gff_fields=split(/\t/,${$transcript_gff{$g}}[0]);
      my $ori=$gff_fields[6];
      if($ori eq "+"){
        if(length($transcript_seqs_5pext_actual{$g})>0){
          @gff_fields=split(/\t/,${$transcript_gff{$g}}[0]);
          $gff_fields[3]-=length($transcript_seqs_5pext_actual{$g});
          ${$transcript_gff{$g}}[0]=join("\t",@gff_fields);
        }
        if(length($transcript_seqs_3pext_actual{$g})>0){
          @gff_fields=split(/\t/,${$transcript_gff{$g}}[-1]);
          $gff_fields[4]+=length($transcript_seqs_3pext_actual{$g});
          ${$transcript_gff{$g}}[-1]=join("\t",@gff_fields);
        }   
      }else{
        if(length($transcript_seqs_3pext_actual{$g})>0){
          @gff_fields=split(/\t/,${$transcript_gff{$g}}[0]);
          $gff_fields[3]-=length($transcript_seqs_3pext_actual{$g});
          ${$transcript_gff{$g}}[0]=join("\t",@gff_fields);
        }
        if(length($transcript_seqs_5pext_actual{$g})>0){
          @gff_fields=split(/\t/,${$transcript_gff{$g}}[-1]);
          $gff_fields[4]+=length($transcript_seqs_5pext_actual{$g});
          ${$transcript_gff{$g}}[-1]=join("\t",@gff_fields);
        }  
      }
    }

#check to see if we truncated the cds -- maybe it is a special protein?
    if($cds_end_on_transcript-$cds_start_on_transcript < $cds_length*$length_fraction){
      $cds_end_on_transcript=$cds_end_on_transcript_original;
      $cds_start_on_transcript=$cds_start_on_transcript_original;
      print "DEBUG too short can't fix $first_codon $last_codon start_cds $cds_start_on_transcript end_cds $cds_end_on_transcript \n";
      $transcript_class{$g}="NA";
      next;
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
    $last_codon=substr($transcript_seqs{$g},$cds_end_on_transcript-3,3);
    $cds_length=$cds_end_on_transcript-$cds_start_on_transcript;
    $transcript_cds_start_codon{$g}=$first_codon if(uc($first_codon) eq "ATG");
    $transcript_cds_end_codon{$g}=$last_codon if(uc($last_codon) eq "TAA" || uc($last_codon) eq "TAG" || uc($last_codon) eq "TGA");
    $transcript_class{$g}="NA" if($transcript_cds_start_codon{$g} eq "MISSING" || $transcript_cds_end_codon{$g} eq "MISSING");#we eliminate transcripts without a start or a stop
    print "DEBUG $first_codon $last_codon start_cds $cds_start_on_transcript end_cds $cds_end_on_transcript protein $transcript_cds{$g} transcript $g cds_length $cds_length transcript length ",length($transcript_seqs{$g})," tstart $tstart pstart $transcript_cds_start{$g} pend $transcript_cds_end{$g} tori $transcript_ori{$g} class $transcript_class{$g}\n";

  }else{#reverse orientation
    print "DEBUG examining protein $transcript_cds{$g} $protein_start{$transcript_cds{$g}} $protein_end{$transcript_cds{$g}} from transdecoder $transcript_cds_start_on_transcript{$g} $transcript_cds_stop_on_transcript{$g}\n";
#we need to determine the position or the CDS start on the transcript, minding the introns
    my $cds_start_on_transcript=0;
    my $cds_end_on_transcript=0;
    my $cds_length=0;
    if(defined($transcript_cds_start_on_transcript{$g}) && defined($transcript_cds_stop_on_transcript{$g})){#from transdecoder
      $cds_start_on_transcript=$transcript_cds_start_on_transcript{$g};
      $cds_end_on_transcript=$transcript_cds_stop_on_transcript{$g};
      $cds_length=$cds_end_on_transcript-$cds_start_on_transcript;
    }else{
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
      $cds_end_on_transcript=$cds_start_on_transcript+$cds_length;
    }
   
    print "DEBUG start_cds $cds_start_on_transcript end_cds $cds_end_on_transcript transcript length ",length($transcript_seqs{$g}),"\n";

    if(($cds_start_on_transcript < 0 || $cds_start_on_transcript > length($transcript_seqs{$g})) && ($cds_end_on_transcript < 0 || $cds_end_on_transcript > length($transcript_seqs{$g}))){#both start and end are messed up
      $transcript_class{$g}="NA";
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

    my $cds_start_on_transcript_original=$cds_start_on_transcript;
    my $cds_end_on_transcript_original=$cds_end_on_transcript;
    $first_codon=substr($transcript_seqs{$g},$cds_start_on_transcript,3);
    $last_codon=substr($transcript_seqs{$g},$cds_end_on_transcript,3);
    print "DEBUG $first_codon $last_codon start_cds $cds_start_on_transcript end_cds $cds_end_on_transcript protein $transcript_cds{$g} transcript $g cds_length $cds_length transcript length ",length($transcript_seqs{$g})," tstart $tstart pstart $transcript_cds_start{$g} pend $transcript_cds_end{$g} tori $transcript_ori{$g}\n";

#checking for in-frame stop codons
    ($cds_start_on_transcript,$cds_end_on_transcript)=fix_in_frame_stops_keep_frame($cds_start_on_transcript_original,$cds_end_on_transcript_original,$transcript_seqs{$g});

#fixing start/stop, extending if needed
    ($cds_start_on_transcript,$cds_end_on_transcript)=fix_start_stop_codon($cds_start_on_transcript,$cds_end_on_transcript,$transcript_seqs{$g});
    $first_codon=substr($transcript_seqs{$g},$cds_start_on_transcript,3);
    $last_codon=substr($transcript_seqs{$g},$cds_end_on_transcript,3);
    ($cds_start_on_transcript,$cds_end_on_transcript,$transcript_seqs_5pext_actual{$g},$transcript_seqs_3pext_actual{$g})=fix_start_stop_codon_ext($cds_start_on_transcript,$cds_end_on_transcript,$transcript_seqs{$g},$transcript_seqs_5pext{$g},$transcript_seqs_3pext{$g}) if((length($transcript_seqs_5pext{$g})==$ext_length && length($transcript_seqs_3pext{$g})==$ext_length) && (not(uc($first_codon) eq "ATG") || not(uc($last_codon) eq "TAA" || uc($last_codon) eq "TAG" || uc($last_codon) eq "TGA"))); #try to extend
#at this point we did all we could to fix the coding region; last thing that remains is to make sure the stop codon is included both into the cds and the transcript
#from here on we assume that the stop codon is included!!!
    $cds_end_on_transcript+=3;
#make sure we are not out of the transcript!!!
    $transcript_seqs_3pext_actual{$g}=substr($transcript_seqs_3pext{$g},0,length($transcript_seqs_3pext_actual{$g})+3) if($cds_end_on_transcript > length($transcript_seqs_5pext_actual{$g}.$transcript_seqs{$g}.$transcript_seqs_3pext_actual{$g}));

#updating the transcript if extended
    if(length($transcript_seqs_5pext_actual{$g})>0 || length($transcript_seqs_3pext_actual{$g})>0){
      $transcript_seqs{$g}=$transcript_seqs_5pext_actual{$g}.$transcript_seqs{$g}.$transcript_seqs_3pext_actual{$g};
#updating the transcript records
      @gff_fields=split(/\t/,${$transcript_gff{$g}}[0]);
      my $ori=$gff_fields[6];
      if($ori eq "+"){
        if(length($transcript_seqs_5pext_actual{$g})>0){
          @gff_fields=split(/\t/,${$transcript_gff{$g}}[0]);
          $gff_fields[3]-=length($transcript_seqs_5pext_actual{$g});
          ${$transcript_gff{$g}}[0]=join("\t",@gff_fields);
        }
        if(length($transcript_seqs_3pext_actual{$g})>0){
          @gff_fields=split(/\t/,${$transcript_gff{$g}}[-1]);
          $gff_fields[4]+=length($transcript_seqs_3pext_actual{$g});
          ${$transcript_gff{$g}}[-1]=join("\t",@gff_fields);
        }   
      }else{
        if(length($transcript_seqs_3pext_actual{$g})>0){
          @gff_fields=split(/\t/,${$transcript_gff{$g}}[0]);
          $gff_fields[3]-=length($transcript_seqs_3pext_actual{$g});
          ${$transcript_gff{$g}}[0]=join("\t",@gff_fields);
        }
        if(length($transcript_seqs_5pext_actual{$g})>0){
          @gff_fields=split(/\t/,${$transcript_gff{$g}}[-1]);
          $gff_fields[4]+=length($transcript_seqs_5pext_actual{$g});
          ${$transcript_gff{$g}}[-1]=join("\t",@gff_fields);
        } 
      }
    }

    if($cds_end_on_transcript-$cds_start_on_transcript < $cds_length*$length_fraction){
      $cds_end_on_transcript=$cds_end_on_transcript_original;
      $cds_start_on_transcript=$cds_start_on_transcript_original;
      print "DEBUG too short can't fix $first_codon $last_codon start_cds $cds_start_on_transcript end_cds $cds_end_on_transcript \n";
      $transcript_class{$g}="NA";
      next;
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
    $last_codon=substr($transcript_seqs{$g},$cds_end_on_transcript-3,3);
    $cds_length=$cds_end_on_transcript-$cds_start_on_transcript;
    $transcript_cds_start_codon{$g}=$first_codon if(uc($first_codon) eq "ATG");
    $transcript_cds_end_codon{$g}=$last_codon if(uc($last_codon) eq "TAA" || uc($last_codon) eq "TAG" || uc($last_codon) eq "TGA");
    $transcript_class{$g}="NA" if($transcript_cds_start_codon{$g} eq "MISSING" || $transcript_cds_end_codon{$g} eq "MISSING");#we eliminate transcripts without at least a start or a stop
#$transcript_class{$g}="NA" if(($transcript_cds_start_codon{$g} eq "MISSING" || $transcript_cds_end_codon{$g} eq "MISSING") && $transcript_source{$g} eq "EviAnnP");#we eliminate incomplete transcripts derived from protein alignments
      print "DEBUG $first_codon $last_codon start_cds $cds_start_on_transcript end_cds $cds_end_on_transcript protein $transcript_cds{$g} transcript $g cds_length $cds_length transcript length ",length($transcript_seqs{$g})," tstart $tstart pstart $transcript_cds_start{$g} pend $transcript_cds_end{$g} tori $transcript_ori{$g} class $transcript_class{$g}\n";
  }
}

#process the loci
#print the headers
print OUTFILE1 "##gff-version 3\n# EviAnn automated annotation\n";
for my $locus(keys %transcripts_cds_loci){
  my @output=();
  my $output_hc_count=0;
  my $output_p_count=0;
  my $gene_feature="";
  my @transcripts_at_loci=split(/\s+/,$transcripts_cds_loci{$locus});
  my @gff_fields=split(/\t/,$transcript{$transcripts_at_loci[0]});
  my $locus_start=$gff_fields[3];
  my $locus_end=$gff_fields[4];
  my @attributes=split(";",$gff_fields[8]);
  $geneID=substr($attributes[1],7);#this is the XLOC
  my $parent=$geneID."-mRNA-";
  my $transcript_index=0;
  #we output transcripts by class code, first = then k and then j, and we record which cds we used; if the cds was used for a higher class we skip the transcript
  for my $source("StringTie","EviAnnP"){
     for my $class ("=","k","j","m","n"){  
      for my $t(@transcripts_at_loci){
        next if(not($transcript_class{$t} eq $class));
        next if(not($transcript_source{$t} eq $source));
        next if(($class eq "m" || $class eq "n") && ($output_hc_count>0 || $source eq "EviAnnP"));#these are low confidence, we use them as last resort, there should be no such codes for EviAnnP
        $output_hc_count++ if($class eq "=" || $class eq "k" || $class eq "j");
        $output_p_count++ if($source eq "EviAnnP");
        print "DEBUG considering transcript $t class $transcript_class{$t} protein $transcript_cds{$t}\n";
        my $protID=$transcript_cds{$t};
        $used_proteins{$protID}=1;
        my $note="";
        my @gff_fields_t=split(/\t/,$transcript{$t});
        my @attributes_t=split(";",$gff_fields_t[8]);
        my $transcriptID=substr($attributes_t[0],3);#this is the source transcript ID
        $transcriptID=$original_transcript_name{$transcriptID} if(defined($original_transcript_name{$transcriptID}));
        my $start_cds=$transcript_cds_start{$t};
        my $end_cds=$transcript_cds_end{$t};
        my $transcript_start=$gff_fields_t[3];
        my $transcript_end=$gff_fields_t[4];
        $transcript_start=$start_cds if($transcript_start > $start_cds);
        $transcript_end=$end_cds if($transcript_end < $end_cds);
        my $transcript_cds_start_index=0;
        my $transcript_cds_end_index=$#{$transcript_gff{$t}};
        $locus_start=$transcript_start if($transcript_start < $locus_start);
        $locus_end=$transcript_end if($transcript_end > $locus_end);
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
#output transcript
        $transcript_index++;
        print "DEBUG output transcript $t class $transcript_class{$t} protein $transcript_cds{$t}\n";
        my $evidence_type="complete";
        $evidence_type="protein_only" if($source eq "EviAnnP");
        $evidence_type="transcript_only" if($protID =~ /^XLOC_/);
        ($evidenceProtID,$junk)=split(/:/,$protID);
        push(@output,$gff_fields[0]."\tEviAnn\tmRNA\t$transcript_start\t$transcript_end\t".join("\t",@gff_fields_t[5..7])."\tID=$parent$transcript_index;Parent=$geneID;EvidenceProteinID=$evidenceProtID;EvidenceTranscriptID=$transcriptID;StartCodon=$transcript_cds_start_codon{$t};StopCodon=$transcript_cds_end_codon{$t};Class=$class;Evidence=$evidence_type;");
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
#assign UTRs based on strand
        my $first_UTR="five_prime_UTR";
        my $first_UTR_label="5UTR";
        my $second_UTR="three_prime_UTR";
        my $second_UTR_label="3UTR";
        if($gff_fields[6] eq "-"){
          $first_UTR="three_prime_UTR";
          $first_UTR_label="3UTR";
          $second_UTR="five_prime_UTR";
          $second_UTR_label="5UTR";
        }

#output first UTR 
        $i=1;
        for(my $j=$first_j;$j<=$transcript_cds_start_index;$j++){
          my @gff_fields=split(/\t/,${$transcript_gff{$t}}[$j]);
          if($j==$transcript_cds_start_index){
            push(@output,$gff_fields[0]."\tEviAnn\t$first_UTR\t$gff_fields[3]\t".($start_cds-1)."\t".join("\t",@gff_fields[5..7])."\tID=$parent$transcript_index:$first_UTR_label:$i;Parent=$parent$transcript_index") if($start_cds>$gff_fields[3]);
          }else{
            push(@output,$gff_fields[0]."\tEviAnn\t$first_UTR\t".join("\t",@gff_fields[3..7])."\tID=$parent$transcript_index:$first_UTR_label:$i;Parent=$parent$transcript_index");
          }
          $i++;
        }
#output cds from exons
        $i=1;
        if($transcript_cds_start_index<$transcript_cds_end_index){
          my @gff_fields=split(/\t/,${$transcript_gff{$t}}[$transcript_cds_start_index]);
          push(@output,$gff_fields[0]."\tEviAnn\tCDS\t$start_cds\t".join("\t",@gff_fields[4..7])."\tID=$parent$transcript_index:cds:$i;Parent=$parent$transcript_index$note");
          $i++;
          for(my $j=$transcript_cds_start_index+1;$j<$transcript_cds_end_index;$j++){
            my @gff_fields=split(/\t/,${$transcript_gff{$t}}[$j]);
            push(@output,$gff_fields[0]."\tEviAnn\tCDS\t".join("\t",@gff_fields[3..7])."\tID=$parent$transcript_index:cds:$i;Parent=$parent$transcript_index$note");
            $i++;
          }
          my @gff_fields=split(/\t/,${$transcript_gff{$t}}[$transcript_cds_end_index]);
          push(@output,$gff_fields[0]."\tEviAnn\tCDS\t$gff_fields[3]\t$end_cds\t".join("\t",@gff_fields[5..7])."\tID=$parent$transcript_index:cds:$i;Parent=$parent$transcript_index$note");
        }else{#single exon
          my @gff_fields=split(/\t/,${$transcript_gff{$t}}[$transcript_cds_start_index]);
          push(@output,$gff_fields[0]."\tEviAnn\tCDS\t$start_cds\t$end_cds\t".join("\t",@gff_fields[5..7])."\tID=$parent$transcript_index:cds:$i;Parent=$parent$transcript_index$note");
        }
#output second UTR
        $i=1;
        for(my $j=$transcript_cds_end_index;$j<=$last_j;$j++){
          my @gff_fields=split(/\t/,${$transcript_gff{$t}}[$j]);
          if($j==$transcript_cds_end_index){
            push(@output,$gff_fields[0]."\tEviAnn\t$second_UTR\t".($end_cds+1)."\t".join("\t",@gff_fields[4..7])."\tID=$parent$transcript_index:$second_UTR_label:$i;Parent=$parent$transcript_index") if($end_cds<$gff_fields[4]);
          }else{
            push(@output,$gff_fields[0]."\tEviAnn\t$second_UTR\t".join("\t",@gff_fields[3..7])."\tID=$parent$transcript_index:$second_UTR_label:$i;Parent=$parent$transcript_index");
          }
          $i++
        }
        print "DEBUG finished transcript $t class $transcript_class{$t} protein $transcript_cds{$t}\n";
      }#end of transcripts loop
    }#end of souce loop
  }#end of class loop
#now we know the locus start ans end, and we can output the gene record
  if(scalar(@output)>0){
    my $dir_factor=0;
    $dir_factor=0.5 if($gff_fields[6] eq "-");
    $gene_record_k{$gff_fields[0]." ".($locus_start+$dir_factor)}="$gff_fields[0]\tEviAnn\tgene\t$locus_start\t$locus_end\t".join("\t",@gff_fields[5..7])."\tID=$geneID;geneID=$geneID;type=protein_coding\n".join("\n",@output)."\n";
    push(@gene_records_k,$gff_fields[0]." ".($locus_start+$dir_factor));
    push(@outputLOCchr,$gff_fields[0]);
    push(@outputLOCbeg,$locus_start);
    push(@outputLOCend,$locus_end);
    push(@outputLOCdir,$gff_fields[6]);
  }
  #print $gff_fields[0]." ".($locus_start+$dir_factor),"\n",$gene_record_k{$gff_fields[0]." ".($locus_start+$dir_factor)};
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
  #check if we already have something at this locus
  my $skip=0;
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
      #check if any of the transcripts overlap known protein coding loci
      $skip=check_overlap($gff_fields_t[0],$gff_fields_t[6],$locus_start,$locus_end);
      for(my $j=1;$j<=$#{$transcript_gff_u{$t}};$j++){
        my @gff_fields_curr=split(/\t/,${$transcript_gff_u{$t}}[$j]);
        my @gff_fields_prev=split(/\t/,${$transcript_gff_u{$t}}[$j-1]);
        $distinct_intron_junctions{"$gff_fields_prev[4] $gff_fields_curr[3]"}=1;
        $total_intron_junctions++;
      }
    }
    $junction_score=scalar(keys %distinct_intron_junctions)/$total_intron_junctions;
  }
  next if($skip);
  #do not output the locus if there are too many disagreements between the intron junctions
  next if($junction_score>0.66);
  #if we got here we can output the transcripts at the locus
  for my $t(@transcripts_at_loci){
    next if(not(defined($transcript_gff_u{$t})));
    my @gff_fields_t=split(/\t/,$transcript_u{$t});
    next if($gff_fields_t[8] =~ /contained_in/);
    my @attributes_t=split(";",$gff_fields_t[8]);
    my $transcriptID=substr($attributes_t[0],3);#this is the source transcript ID
    $transcriptID=$original_transcript_name{$transcriptID} if(defined($original_transcript_name{$transcriptID}));
    my ($original_name,$num_samples,$tpm)=split(/:/,$transcriptID);
    next if(($num_samples<2 || $tpm<2) && $transcriptID =~ /^MSTRG/ && substr($attributes_t[0],3) =~ /^TCONS/);#require this transcript to be in minimum 2 samples with TPM>=1, unless it is assembled from reference
    next if(scalar(@{$transcript_gff_u{$t}})<2 && substr($attributes_t[0],3) =~ /^TCONS/); #must be multi-exon if final output
    $transcript_index++;
    push(@output,"$gff_fields_t[0]\tEviAnn\tmRNA\t".join("\t",@gff_fields_t[3..7])."\tID=$parent$transcript_index;Parent=$geneID;EvidenceTranscriptID=$transcriptID");
    my $i=1;
    for my $x(@{$transcript_gff_u{$t}}){
      my @gff_fields=split(/\t/,$x);
      push(@output,"$gff_fields_t[0]\tEviAnn\t".join("\t",@gff_fields[2..7])."\tID=$parent$transcript_index:exon:$i;Parent=$parent$transcript_index");
      $i++;
    }
  }
  if($transcript_index>0){
    my $dir_factor=0;
    $dir_factor=0.5 if($gff_fields[6] eq "-");
    $gene_record_u{$gff_fields[0]." ".($locus_start+$dir_factor)}="$gff_fields[0]\tEviAnn\tgene\t$locus_start\t$locus_end\t".join("\t",@gff_fields[5..7])."\tID=$geneID;geneID=$geneID;type=lncRNA;junction_score=$junction_score;\n".join("\n",@output)."\n";
    push(@gene_records_u,$gff_fields[0]." ".($locus_start+$dir_factor));
  }
}

#output unused proteins
#we will then look at them, pick only one per locus that best matches uniprot and join them in at the second pass
my $fake_utr=0;
foreach my $p(keys %protein){
  next if(defined($used_proteins{$p}));
  my @gff_fields_p=split(/\t/,$protein{$p});
  #it is better to allow proteins to detect new isoforms at existing locations
  #next if(check_overlap($gff_fields_p[0],$gff_fields_p[6],$gff_fields_p[3],$gff_fields_p[4]));
  my $ptstart=$gff_fields_p[3]-$fake_utr>0 ? $gff_fields_p[3]-$fake_utr:1;
  my $ptend=$gff_fields_p[4]+$fake_utr<=length($genome_seqs{$gff_fields_p[0]}) ? $gff_fields_p[4]+$fake_utr:length($genome_seqs{$gff_fields_p[0]});

  print OUTFILE4 "$gff_fields_p[0]\tEviAnnP\t$gff_fields_p[2]\t",$ptstart,"\t",$ptend,"\t",join("\t",@gff_fields_p[5..$#gff_fields_p]),"\n";
  for(my $j=0;$j<=$#{$protein_cds{$p}};$j++){
    my @gff_fields_c=split(/\t/,${$protein_cds{$p}}[$j]);
    if($#{$protein_cds{$p}}==0){#add "fake" 5' utr and 3' utr
      print OUTFILE4 "$gff_fields_c[0]\tEviAnnP\texon\t$ptstart\t$ptend\t",join("\t",@gff_fields_c[5..$#gff_fields_c]),"\n";
    }elsif($j==0){#add "fake" 5' utr
      print OUTFILE4 "$gff_fields_c[0]\tEviAnnP\texon\t$ptstart\t",join("\t",@gff_fields_c[4..$#gff_fields_c]),"\n";
    }elsif($j==$#{$protein_cds{$p}}){
      print OUTFILE4 "$gff_fields_c[0]\tEviAnnP\texon\t$gff_fields_c[3]\t$ptend\t",join("\t",@gff_fields_c[5..$#gff_fields_c]),"\n";
    }else{
      print OUTFILE4 "$gff_fields_c[0]\tEviAnnP\texon\t",join("\t",@gff_fields_c[3..$#gff_fields_c]),"\n";
    }
  }
  for(my $j=0;$j<=$#{$protein_cds{$p}};$j++){
    my @gff_fields_c=split(/\t/,${$protein_cds{$p}}[$j]);
    print OUTFILE4 "$gff_fields_c[0]\tEviAnnP\tCDS\t",join("\t",@gff_fields_c[3..$#gff_fields_c]),"\n";
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

@gene_records_sorted=sort mysort @gene_records_u;
%output=();
foreach $g(@gene_records_sorted){
  if(not(defined($output{$g}))){
    print OUTFILE3 $gene_record_u{$g};
    $output{$g}=1;
  }
}

sub check_overlap{
  my $chrom=$_[0];
  my $dir=$_[1];
  my $l_start=$_[2];
  my $l_end=$_[3];
  my $overlap=0;
  for(my $i=0;$i<=$#outputLOCchr;$i++){
    if($chrom eq $outputLOCchr[$i] &&  (($l_start<=$outputLOCend[$i] && $l_start>=$outputLOCbeg[$i])||($l_end<=$outputLOCend[$i] && $l_end>=$outputLOCbeg[$i]))){
      $overlap=1;
      $i=$#outputLOCchr;
    }
  }
  return($overlap);
}

sub mysort{
my ($chroma,$coorda)=split(/\s+/,$a);
my ($chromb,$coordb)=split(/\s+/,$b);
return($chroma cmp $chromb || $coorda <=>$coordb);
}

sub fix_start_stop_codon_ext{
  my $cds_start_on_transcript=$_[0];
  my $cds_end_on_transcript=$_[1];
  my $transcript_seq=$_[2];
  my $transcript_5pext=$_[3];
  my $transcript_3pext=$_[4];
  my $ext_length=length($transcript_5pext);
  #we do not do 5' extension if the first codon is start
  my $first_codon=substr($transcript_seq,$cds_start_on_transcript,3);
  my $last_codon=substr($transcript_seq,$cds_end_on_transcript,3);
  print "DEBUG extending start and stop starting at $cds_start_on_transcript $cds_end_on_transcript\n";
  my ($cds_start_on_transcript_ext,$cds_end_on_transcript_ext)=fix_start_stop_codon($cds_start_on_transcript+$ext_length,$cds_end_on_transcript+$ext_length,$transcript_5pext.$transcript_seq.$transcript_3pext);
  if(uc($first_codon) eq "ATG"){#ignore the extension
    $transcript_5pext="";
    print "DEBUG no 5p extension needed $cds_start_on_transcript_ext\n";
  }elsif($cds_start_on_transcript_ext>$ext_length-1){
    $transcript_5pext="";
    $cds_start_on_transcript=$cds_start_on_transcript_ext-$ext_length;
    print "DEBUG no 5p extension $cds_start_on_transcript_ext\n";
  }else{
    $transcript_5pext=substr($transcript_5pext,$cds_start_on_transcript_ext);
    #check the extension for AG -- acceptor sites, if found, do not extend
    my $found_acceptor=0;
    if(defined($ARGV[4])){
      for(my $j=3;$j<length($transcript_5pext)-2;$j++){
        if(substr($transcript_5pext,$j,2) eq "AG"){
          print "DEBUG found acceptor at $j\n";
          my $index5=$j;
          my $ext_seq=uc($transcript_5pext.$transcript_seq);
          if($index5>=($acceptor_length-5)){
            $score_seq=substr($ext_seq,$index5-($acceptor_length-5),$acceptor_length);
          }elsif($index5<($acceptor_length-5)){
            $score_seq="N"x(($acceptor_length-5)-$index5).substr($ext_seq,0,$index5+5);
          }
#score the extension
          $ext_score=0;
          for(my $i=0;$i<$acceptor_length;$i++){
            $ext_score+=$acceptor_freq[$i][$code{substr($score_seq,$i,1)}];
            print "DEBUG ",substr($score_seq,$i,1)," ",$acceptor_freq[$i][$code{substr($score_seq,$i,1)}],"\n";
          }
          print "DEBUG $index5 $score_seq $ext_score\n";
          if($ext_score > 5){
            $found_acceptor=1;
            $j=length($transcript_5pext);
          }
        }
      }
    }
    if($found_acceptor==0){
      $cds_start_on_transcript=0;
      print "DEBUG extend 5p $cds_start_on_transcript_ext $transcript_5pext ",length($transcript_5pext),"\n";
    }else{
      print "DEBUG reject 5p $cds_start_on_transcript_ext $transcript_5pext ",length($transcript_5pext),"\n";
      $transcript_5pext="";
    }
  }
  if(uc($last_codon) eq "TAA" || uc($last_codon) eq "TAG" || uc($last_codon) eq "TGA"){
    print "DEBUG no 3p extension needed $cds_end_on_transcript_ext\n";
    $transcript_3pext="";
    $cds_end_on_transcript+=length($transcript_5pext);
  }elsif($cds_end_on_transcript_ext<$ext_length+length($transcript_seq)){
    $transcript_3pext="";
    print "DEBUG no 3p extension $cds_end_on_transcript_ext\n";
    $cds_end_on_transcript=$cds_end_on_transcript_ext-$ext_length+length($transcript_5pext);
  }else{
    $transcript_3pext=substr($transcript_3pext,0,$cds_end_on_transcript_ext-length($transcript_seq)-$ext_length+3);
    my $found_donor=0;
    if(defined($ARGV[4])){
      for(my $j=0;$j<length($transcript_3pext)-3;$j++){
        if(substr($transcript_3pext,$j,2) eq "GT"){
          print "DEBUG found acceptor at $j\n";
          my $index3=$j;
          my $ext_seq=uc($transcript_seq.$transcript_3pext);
          if(length($transcript_3pext)-$index3>=($donor_length-3)){
            $score_seq=substr($ext_seq,length($transcript_seq)+$index3-3,$donor_length);
          }else{
            $score_seq=substr($ext_seq,length($transcript_seq)+$index3-3)."N"x(($donor_length-3)-$index3);
          }
          #score the extension
          $ext_score=0;
          for(my $i=0;$i<$donor_length;$i++){
            $ext_score+=$donor_freq[$i][$code{substr($score_seq,$i,1)}];
            print "DEBUG ",substr($score_seq,$i,1)," ",$donor_freq[$i][$code{substr($score_seq,$i,1)}],"\n";
          }
          print "DEBUG $index3 $score_seq $ext_score\n";
          if($ext_score > 4){
            $found_donor=1;
            $j=length($transcript_3pext);
          }
        }
      }
    }
    if($found_donor==0){
      $cds_end_on_transcript=$cds_end_on_transcript_ext-$ext_length+length($transcript_5pext);
      print "DEBUG extend 3p $cds_end_on_transcript_ext $transcript_3pext ",length($transcript_3pext),"\n";
    }else{
      print "DEBUG reject 3p $cds_end_on_transcript_ext $transcript_3pext ",length($transcript_3pext),"\n";
      $transcript_3pext="";
      $cds_end_on_transcript+=length($transcript_5pext);
    }
  }
  return($cds_start_on_transcript,$cds_end_on_transcript,$transcript_5pext,$transcript_3pext);
}

sub fix_start_stop_codon{
  my $cds_start_on_transcript=$_[0];
  my $cds_end_on_transcript=$_[1];
  my $transcript_seq=$_[2];
  my $first_codon=substr($transcript_seq,$cds_start_on_transcript,3);
  my $last_codon=substr($transcript_seq,$cds_end_on_transcript,3);
  print "DEBUG fixing start and stop starting at $cds_start_on_transcript $cds_end_on_transcript $first_codon $last_codon\n";
  my $i;
  my @found=();
  my @fscores=();
  for($i=$cds_start_on_transcript;$i>=0;$i-=3){
    if(uc(substr($transcript_seq,$i,3)) eq "ATG"){
      push(@found,$i);
      my $score=0;
      #if($i-6>=0){#score kozak-like sequence
      #  print "DEBUG kozak ",uc(substr($transcript_seq,$i-6,10)),"\n";
      #  for(my $j=0;$j<=9;$j++){
      #    #print "Letter ",uc(substr($transcript_seq,$i-6+$j,1))," index ",$letter_index{uc(substr($transcript_seq,$i-6+$j,1))}," score ",$freq[$letter_index{uc(substr($transcript_seq,$i-6+$j,1))}][$j],"\n";
      #    $score+=$freq[$letter_index{uc(substr($transcript_seq,$i-6+$j,1))}][$j];
      #  }
      #  $score+=10/($i+10);
      #  print "DEBUG found candidate start at $i with score $score\n";
      #}else{
      #  $score=9;
      #}
      push(@fscores,$score);
    }
#stop if found a stop
    last if(uc(substr($transcript_seq,$i,3)) eq "TAA" || uc(substr($transcript_seq,$i,3)) eq "TAG" || uc(substr($transcript_seq,$i,3)) eq "TGA");
  } 
  if(scalar(@found)>0){
    #my $max_score=-1;
    #for(my $j=0;$j<=$#fscores;$j++){
    #  $max_score=$fscores[$j] if($fscores[$j]>$max_score);
    #}
    #for(my $j=0;$j<=$#found;$j++){
    #  $cds_start_on_transcript=$found[$j] if($fscores[$j] == $max_score);
    #}
    $cds_start_on_transcript=$found[-1];
    print "DEBUG found new start codon upstream at $cds_start_on_transcript\n";
  }else{ 
    print "DEBUG failed to find new start codon, looking downstream\n";
    my $foundU=-1;
    for($i=$cds_start_on_transcript+3;$i<$cds_end_on_transcript;$i+=3){
      if(uc(substr($transcript_seq,$i,3)) eq "ATG"){
        $foundU=$i;
        last;
      }
    }
    if($foundU>-1){
      print "DEBUG found new start codon upstream at $foundU\n";
      $cds_start_on_transcript=$foundU;
    }else{
      print "DEBUG failed to find new start codon\n";
    }
  }
  if(not(uc($last_codon) eq "TAA" || uc($last_codon) eq "TAG" || uc($last_codon) eq "TGA")){
    my $i;
    for($i=$cds_start_on_transcript+3;$i<length($transcript_seq);$i+=3){
      last if(uc(substr($transcript_seq,$i,3)) eq "TAA" || uc(substr($transcript_seq,$i,3)) eq "TAG" || uc(substr($transcript_seq,$i,3)) eq "TGA");
    } 
    if($i<length($transcript_seq)){
      print "DEBUG found new stop codon downstream at $i\n";
      $cds_end_on_transcript=$i;
    }else{
      print "DEBUG failed to find new stop codon downstream\n";
    }
  }else{
    print "DEBUG last codon is stop $last_codon $cds_end_on_transcript\n";
  }
  return($cds_start_on_transcript,$cds_end_on_transcript);
}

sub fix_in_frame_stops_keep_frame{
  my $in_frame_stop=0;
  my $cds_start_on_transcript=$_[0];
  my $cds_end_on_transcript=$_[1];
  my $transcript_seq=$_[2];
  my $frame0_start=$cds_start_on_transcript;
  my $frame0_end=$cds_end_on_transcript;
  print "DEBUG checking for in frame stop starting $cds_start_on_transcript $cds_end_on_transcript\n";
  for($i=$frame0_start;$i<$frame0_end;$i+=3){
    if(uc(substr($transcript_seq,$i,3)) eq "TAA" || uc(substr($transcript_seq,$i,3)) eq "TAG" || uc(substr($transcript_seq,$i,3)) eq "TGA"){
      $in_frame_stop=$i;
      $frame0_end=$i-3;
      last;
    }
  }
  if($in_frame_stop){
    print "DEBUG found in-frame stop at $in_frame_stop not switching frame\n";
    $cds_end_on_transcript=$frame0_end;
  }else{
    print "DEBUG no in frame stop\n";
  }
  return($cds_start_on_transcript,$cds_end_on_transcript);
}

  
