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
my $length_fraction=0.75;
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
    if(defined($transcript{$ID})){
    #here we need to fix the first and the last exons so that the CDS does not stick out from the boundaries of the transcript
      my @gff_fields_t=split(/\t/,$transcript{$ID});
      my @gff_fields=split(/\t/,$exons[0]);
      my @gff_fields_p=split(/\t/,${$protein_cds{$transcript_cds{$ID}}}[0]);
      if($gff_fields[3]>$gff_fields_p[3]){
        print "DEBUG longer 5' for $transcript{$ID}\n";
        print join("\n",@{$protein_cds{$transcript_cds{$ID}}}),"\n";
        print join("\n",@exons),"\n";
        if($gff_fields[3]<=$gff_fields_p[4]){
          $gff_fields[3]=$gff_fields_p[3];
          $gff_fields_t[3]=$gff_fields[3];
          $exons[0]=join("\t",@gff_fields);
        }else{
          print "DEBUG extra exons 5' for $ID\n";#add exons
          my @pexonsb=();
          my @pexonse=();
          my $found_j=0;
          push(@pexonsb,$gff_fields_p[3]);
          push(@pexonse,$gff_fields_p[4]);
          for($i=1;$i<=$#{$protein_cds{$transcript_cds{$ID}}};$i++){
            @gff_fields_p=split(/\t/,${$protein_cds{$transcript_cds{$ID}}}[$i]);
            push(@pexonsb,$gff_fields_p[3]);
            push(@pexonse,$gff_fields_p[4]);
            if($gff_fields_p[4]==$gff_fields[4]){
              $found_j=1;
              last;
            }
          }
          if($found_j){
            #now we remove the firs exon and replace exons
            shift(@exons);
            for($i=$#pexonsb;$i>=0;$i--){
              my $exonline=join("\t",@gff_fields[0..2])."\t$pexonsb[$i]\t$pexonse[$i]\t".join("\t",@gff_fields[5..$#gff_fields]);
              unshift(@exons,$exonline);
            }
            print "DEBUG fixed 5' for $ID\n";
            print join("\n",@exons),"\n";
            $gff_fields_t[3]=$pexonsb[0];
          }
        }
      }
      @gff_fields=split(/\t/,$exons[-1]);
      @gff_fields_p=split(/\t/,${$protein_cds{$transcript_cds{$ID}}}[-1]);
      if($gff_fields[4]<$gff_fields_p[4]){
        print "DEBUG longer 3' for $transcript{$ID}\n";
        print join("\n",@{$protein_cds{$transcript_cds{$ID}}}),"\n";
        print join("\n",@exons),"\n";
        if($gff_fields[4]>=$gff_fields_p[3]){
          $gff_fields[4]=$gff_fields_p[4];
          $gff_fields_t[4]=$gff_fields[4];
          $exons[-1]=join("\t",@gff_fields);
        }else{
          print "DEBUG extra exons 3' for $ID\n";#add exons
          my @pexonsb=();
          my @pexonse=();
          my $found_j=0;
          push(@pexonsb,$gff_fields_p[3]);
          push(@pexonse,$gff_fields_p[4]);
          for($i=$#{$protein_cds{$transcript_cds{$ID}}}-1;$i>=0;$i--){
            @gff_fields_p=split(/\t/,${$protein_cds{$transcript_cds{$ID}}}[$i]);
            unshift(@pexonsb,$gff_fields_p[3]);
            unshift(@pexonse,$gff_fields_p[4]);
            if($gff_fields_p[3]==$gff_fields[3]){
              $found_j=1;
              last;
            }
          }
          if($found_j){
            #now we remove the firs exon and replace exons
            pop(@exons);
            for($i=0;$i<=$#pexonsb;$i++){
              my $exonline=join("\t",@gff_fields[0..2])."\t$pexonsb[$i]\t$pexonse[$i]\t".join("\t",@gff_fields[5..$#gff_fields]);
              push(@exons,$exonline);
            }
            print "DEBUG fixed 3' for $ID\n";
            print join("\n",@exons),"\n";
            $gff_fields_t[4]=$pexonse[-1];
          }
        }
      }
      $transcript{$ID}=join("\t",@gff_fields_t);
      print "DEBUG updated transcript $transcript{$ID}\n";
      $transcript_gff{$ID}=[@exons];
      $transcript_ori{$ID}=$tori;
    }elsif(defined($transcript_u{$ID})){
      if($#exons > 0){#not interested in single exon
        $transcript_gff_u{$ID}=[@exons];
      }
    }
    @exons=();
    $ID=substr($attributes[0],3);#this is the transcript_id
    $tstart=$gff_fields[3];
    $tend=$gff_fields[4];
    $tori=$gff_fields[6];
    my $class_code="";
    my $protID="";
    my $locID="";
    foreach my $attr(@attributes){
      $class_code=substr($attr,11,1) if($attr =~ /^class_code=/);
      $protID=substr($attr,8) if($attr =~ /^cmp_ref=/);
      $locID=substr($attr,7) if($attr =~ /^geneID=/);
      $locID=substr($attr,5) if($attr =~ /^xloc=/);
    }
    if($gff_fields[1] eq "EviAnnP" && $ID =~ /_EXTERNAL$/){
      $protID=$ID;
      $class="=";
    }
    print "DEBUG read transcript $ID class $class_code protein $protID locus $locID\n";
    if($class_code =~ /i|y|u|o|x/){
      $transcript_u{$ID}=$line;
      $transcripts_only_loci{$locID}.="$ID ";
    }elsif(($protID =~ /_EXTERNAL$/ && $class_code =~ /k|=/) || not($protID =~ /_EXTERNAL$/)){
      $transcript{$ID}=$line;
      die("Protein $protID is not defined for protein coding transcript $ID") if(not(defined($protein{$protID})));
      $transcript_cds{$ID}=$protID;
      $transcript_cds_start{$ID}=$protein_start{$protID};
      $transcript_cds_start_codon{$ID}="INVALID";
      $transcript_cds_end{$ID}=$protein_end{$protID};
      $transcript_cds_end_codon{$ID}="INVALID";
      $transcript_class{$ID}=$class_code;
      $transcript_source{$ID}=$gff_fields[1];
      $transcripts_cds_loci{$locID}.="$ID ";
    }
  }elsif($gff_fields[2] eq "exon"){
    push(@exons,$line) if(defined($transcript{$ID}) || defined($transcript_u{$ID}));
  }
}
if(defined($transcript{$ID})){
#here we need to fix the first and the last exons so that the CDS does not stick out from the boundaries of the transcript
  my @gff_fields=split(/\t/,$exons[0]);
  $gff_fields[3]=$protein_start{$transcript_cds{$ID}} if($gff_fields[3]>$protein_start{$transcript_cds{$ID}});
  $exons[0]=join("\t",@gff_fields);
  my @gff_fields=split(/\t/,$exons[-1]);
  $gff_fields[4]=$protein_end{$transcript_cds{$ID}} if($gff_fields[4]<$protein_end{$transcript_cds{$ID}});
  $exons[-1]=join("\t",@gff_fields);
  $transcript_gff{$ID}=[@exons];
  $transcript_ori{$ID}=$tori;
}elsif(defined($transcript_u{$ID})){
  unless($#exons==-1){#not interested in single exon
    $transcript_gff_u{$ID}=[@exons];
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

my %mito_contigs=();
if(defined($ARGV[3])){
  open(FILE,$ARGV[3]);
  while(my $line=<FILE>){
    chomp($line);
    $mito_contigs{$line}=1;
  }
}

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
my $gcode=0;
for my $g(keys %transcript_cds){
  my @gff_fields_t=split(/\t/,$transcript{$g});
  my $tstart=$gff_fields_t[3];
  my $tend=$gff_fields_t[4];
  $gcode= defined($mito_contigs{$gff_fields_t[0]}) ? 1 : 0;
  print "\nDEBUG protein $transcript_cds{$g} transript $g length ",length($transcript_seqs{$g}),"\n";
  #do not mess with external CDSs
  next if( $g =~ /_EXTERNAL$/||($transcript_source{$g} eq "EviAnnP" && $transcript_cds{$g} =~ /_EXTERNAL$/));
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
    my $cds_start_on_transcript_original=$cds_start_on_transcript;
    my $cds_end_on_transcript_original=$cds_end_on_transcript;
    $first_codon=substr($transcript_seqs{$g},$cds_start_on_transcript,3);
    $last_codon=substr($transcript_seqs{$g},$cds_end_on_transcript-3,3);
    print "DEBUG $first_codon $last_codon start_cds $cds_start_on_transcript end_cds $cds_end_on_transcript protein $transcript_cds{$g} transcript $g cds_length $cds_length transcript length ",length($transcript_seqs{$g})," tstart $tstart pstart $transcript_cds_start{$g} pend $transcript_cds_end{$g} tori $transcript_ori{$g}\n";

    $cds_end_on_transcript=check_in_frame_stops($cds_start_on_transcript_original,$cds_end_on_transcript_original,$transcript_seqs{$g});

    if((($cds_start_on_transcript_original < 0 || $cds_start_on_transcript_original > length($transcript_seqs{$g})) && ($cds_end_on_transcript < 0 || $cds_end_on_transcript > length($transcript_seqs{$g}))) || $cds_length %3 >0){#both start and end are messed up -- send to transdecoder
      print OUTFILE2 ">$g\n$transcript_seqs{$g}\n";
      my @pn=split(/:/,$transcript_cds{$g});
      print OUTFILE3 "$pn[0]\n";
      print "DEBUG broken CDS start and stop outside or non divisible or in-frame stop $g\n";
      next;
    }

#fixing start/stop    
    ($cds_start_on_transcript,$cds_end_on_transcript)=fix_start_stop_codon($cds_start_on_transcript,$cds_end_on_transcript,$transcript_seqs{$g});
    $first_codon=substr($transcript_seqs{$g},$cds_start_on_transcript,3);
    $last_codon=substr($transcript_seqs{$g},$cds_end_on_transcript-3,3);
    $cds_length=$cds_end_on_transcript-$cds_start_on_transcript;
    if($cds_end_on_transcript < $cds_start_on_transcript_original || $cds_start_on_transcript > $cds_end_on_transcript_original ||  (not(valid_stop($last_codon,$gcode)) || not(valid_start($first_codon))) || ($cds_end_on_transcript-$cds_start_on_transcript)<$length_fraction*($cds_end_on_transcript_original-$cds_start_on_transcript_original)){
      print OUTFILE2 ">$g\n$transcript_seqs{$g}\n";
      my @pn=split(/:/,$transcript_cds{$g});
      print OUTFILE3 "$pn[0]\n";
      print "DEBUG broken CDS no start and stop or too short $g $first_codon $last_codon $cds_start_on_transcript $cds_end_on_transcript $cds_start_on_transcript_original $cds_end_on_transcript_original\n";
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
      print "DEBUG CDS OK $first_codon $last_codon start_cds $cds_start_on_transcript end_cds $cds_end_on_transcript protein $transcript_cds{$g} transcript $g cds_length $cds_length transcript length ",length($transcript_seqs{$g})," tstart $tstart pstart $transcript_cds_start{$g} pend $transcript_cds_end{$g} tori $transcript_ori{$g}\n";
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
    my $cds_start_on_transcript_original=$cds_start_on_transcript;
    my $cds_end_on_transcript_original=$cds_end_on_transcript;
    $first_codon=substr($transcript_seqs{$g},$cds_start_on_transcript,3);
    $last_codon=substr($transcript_seqs{$g},$cds_end_on_transcript-3,3);
    print "DEBUG $first_codon $last_codon start_cds $cds_start_on_transcript end_cds $cds_end_on_transcript protein $transcript_cds{$g} transcript $g cds_length $cds_length transcript length ",length($transcript_seqs{$g})," tstart $tstart pstart $transcript_cds_start{$g} pend $transcript_cds_end{$g} tori $transcript_ori{$g}\n";

    $cds_end_on_transcript=check_in_frame_stops($cds_start_on_transcript_original,$cds_end_on_transcript_original,$transcript_seqs{$g});
    
    if((($cds_start_on_transcript_original < 0 || $cds_start_on_transcript_original > length($transcript_seqs{$g})) && ($cds_end_on_transcript < 0 || $cds_end_on_transcript > length($transcript_seqs{$g}))) || $cds_length %3 >0){#both start and end are messed up -- send to transdecoder
      print OUTFILE2 ">$g\n$transcript_seqs{$g}\n";
      my @pn=split(/:/,$transcript_cds{$g});
      print OUTFILE3 "$pn[0]\n";
      print "DEBUG broken CDS start and stop outside or non divisible or in-frame stop $g\n";
      next;
    }

#fixing start/stop    
    ($cds_start_on_transcript,$cds_end_on_transcript)=fix_start_stop_codon($cds_start_on_transcript,$cds_end_on_transcript,$transcript_seqs{$g});
    $first_codon=substr($transcript_seqs{$g},$cds_start_on_transcript,3);
    $last_codon=substr($transcript_seqs{$g},$cds_end_on_transcript-3,3);
    $cds_length=$cds_end_on_transcript-$cds_start_on_transcript;
    if($cds_end_on_transcript < $cds_start_on_transcript_original || $cds_start_on_transcript > $cds_end_on_transcript_original  || (not(valid_stop($last_codon,$gcode)) || not(valid_start($first_codon))) || ($cds_end_on_transcript-$cds_start_on_transcript)<$length_fraction*($cds_end_on_transcript_original-$cds_start_on_transcript_original)){
      print OUTFILE2 ">$g\n$transcript_seqs{$g}\n";
      my @pn=split(/:/,$transcript_cds{$g});
      print OUTFILE3 "$pn[0]\n";
      print "DEBUG broken CDS no start and stop or too short $g $first_codon $last_codon $cds_start_on_transcript $cds_end_on_transcript $cds_start_on_transcript_original $cds_end_on_transcript_original\n";
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
      print "DEBUG CDS OK $first_codon $last_codon start_cds $cds_start_on_transcript end_cds $cds_end_on_transcript protein $transcript_cds{$g} transcript $g cds_length $cds_length transcript length ",length($transcript_seqs{$g})," tstart $tstart pstart $transcript_cds_start{$g} pend $transcript_cds_end{$g} tori $transcript_ori{$g}\n";
      print OUTFILE1 ">$g $transcript_cds_start{$g} $transcript_cds_end{$g}\n$transcript_seqs{$g}\n";
    }
  }
}

sub fix_start_stop_codon{
  my $cds_start_on_transcript=$_[0];
  my $cds_end_on_transcript=$_[1];
  my $transcript_seq=$_[2];
  my $first_codon=substr($transcript_seq,$cds_start_on_transcript,3);
  my $last_codon=substr($transcript_seq,$cds_end_on_transcript-3,3);
  print "DEBUG fixing start and stop starting at $cds_start_on_transcript $cds_end_on_transcript $first_codon $last_codon\n";
  if(valid_start($first_codon) && valid_stop($last_codon,$gcode)){
    print "DEBUG both codons are OK\n";
    return($cds_start_on_transcript,$cds_end_on_transcript);
  }elsif(valid_start($first_codon)){
    print "DEBUG stop broken\n";
    my $foundS=-1;
    for(my $i=$cds_start_on_transcript+3;$i<=length($transcript_seq)-3;$i+=3){
      if(valid_stop(substr($transcript_seq,$i,3),$gcode)){
        print "DEBUG found new stop codon downstream at $i, ",substr($transcript_seq,$i,3),"\n";
        $cds_end_on_transcript=$i+3;
        $foundS=$i;
        last;
      }
    }
    print "DEBUG failed to find new stop codon\n" if($foundS==-1);
  }elsif(valid_stop($last_codon,$gcode)){
    print "DEBUG start broken\n";
    my $foundU=-1;
    for(my $i=$cds_end_on_transcript-6;$i>=0;$i-=3){
      last if(valid_stop(substr($transcript_seq,$i,3),$gcode));
      $foundU=$i if(valid_start(substr($transcript_seq,$i,3)));
    }
    if($foundU>-1){
      print "DEBUG found new start codon upstream at $foundU\n";
      $cds_start_on_transcript=$foundU;
    }
  }else{#both broken, trust the frame look for a start
    print "DEBUG both start and stop broken\n";
    my $foundU=-1;
    for(my $i=$cds_start_on_transcript;$i>=0;$i-=3){
      last if(valid_stop(substr($transcript_seq,$i,3),$gcode));
      $foundU=$i if(valid_start(substr($transcript_seq,$i,3)));
    }
    if($foundU>-1){
      $cds_start_on_transcript=$foundU;
      print "DEBUG found new start codon upstream at $cds_start_on_transcript\n";
    }else{
      print "DEBUG failed to find new start codon, looking downstream\n";
      my $foundD=-1;
      for(my $i=$cds_start_on_transcript+3;$i<$cds_end_on_transcript-3;$i+=3){
        last if(valid_stop(substr($transcript_seq,$i,3),$gcode));
        if(valid_start(substr($transcript_seq,$i,3))){
          $foundD=$i;
          last;
        }
      }
      if($foundD>-1){
        print "DEBUG found new start codon upstream at $foundD\n";
        $cds_start_on_transcript=$foundD;
      }else{
        print "DEBUG failed to find new start codon\n";
      }
    }
    my $foundS=-1;
    for(my $i=$cds_start_on_transcript+3;$i<=length($transcript_seq)-3;$i+=3){
      if(valid_stop(substr($transcript_seq,$i,3),$gcode)){
        print "DEBUG found new stop codon downstream at $i, ",substr($transcript_seq,$i,3),"\n";
        $cds_end_on_transcript=$i+3;
        $foundS=$i;
        last;
      }
    }
    print "DEBUG failed to find new stop codon\n" if($foundS==-1);
  }
  return($cds_start_on_transcript,$cds_end_on_transcript);
}

sub valid_start{
  my $codon=uc($_[0]);
  if(length($codon)==3 && $codon eq "ATG"){
    return(1);
  }else{
    return(0);
  }
} 

sub valid_stop{
  my $codon=uc($_[0]);
  my $type=$_[1];
  if($type==0){
    if(length($codon)==3 && ($codon eq "TAG" || $codon eq "TAA" || $codon eq "TGA")){
      return(1);
    }else{
      return(0);
    }
  }elsif($type==1){
    if(length($codon)==3 && ($codon eq "AGA" || $codon eq "AGG" || $codon eq "TAA" || $codon eq "TAG" )){
      return(1);
    }else{
      return(0);
    }
  }else{
    return(0);
  }
}

sub check_in_frame_stops{
  my $in_frame_stop=-1;  
  my $cds_start_on_transcript=$_[0];
  my $cds_end_on_transcript=$_[1];
  my $transcript_seq=$_[2];
  print "DEBUG checking for in-frame stops $cds_start_on_transcript $cds_end_on_transcript\n";
  for($i=$cds_start_on_transcript;$i<$cds_end_on_transcript-3;$i+=3){
    if(valid_stop(substr($transcript_seq,$i,3),$gcode)){
      $in_frame_stop=$i;
      last;
    }
  }
  if($in_frame_stop>-1){
    print "DEBUG found in-frame stop at $in_frame_stop\n";
    $cds_end_on_transcript=$in_frame_stop+3;
  }
  return($cds_end_on_transcript);
}
