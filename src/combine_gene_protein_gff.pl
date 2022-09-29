#!/usr/bin/env perl
#This code adds gene and CDS record for every locus in GFF file
my $geneID="";
my @exons=();
my %loci;
my %protein_cds;
my %protein;
my %transcript_gff;
my %transcript_cds;
my @gene_records;
my %gene_record;
my $protID="";
my $dir="";
my %used_proteins;

#this is output of gffcompare -D -o protuniq ../GCF_000001735.4_TAIR10.1_genomic.fna.GCF_000001735.4_TAIR10.1_protein.faa.palign.gff
while(my $line=<STDIN>){#we just read in the whole file
  chomp($line);
  my @gff_fields=split(/\t/,$line);
  my @attributes=split(";",$gff_fields[8]);
  if($gff_fields[2] eq "gene"){
    if(not($protID eq "")){
      $protein_cds{$protID}=[@exons];
    }
    @exons=();
    $protID=substr($attributes[0],3);#this is protein name
    $protein{$protID}=$line;
  }elsif($gff_fields[2] eq "CDS"){
    push(@exons,$line);
  }
}
if(not($protID eq "")){
  $protein_cds{$protID}=[@exons];
}
@exons=();

#this is output of gffcompare -r protalign.gtf ../GCF_000001735.4_TAIR10.1_genomic.clean.fna.gtf -o protref
#here we load up all transcripts that matched proteins
open(FILE,$ARGV[0]);
while(my $line=<FILE>){
  chomp($line);
  my @gff_fields=split(/\t/,$line);
  my @attributes=split(";",$gff_fields[8]);
  if($gff_fields[2] eq "transcript"){
    $transcript_gff{$geneID}=[@exons] if(defined($transcript{$geneID}));
    if(defined($transcript_u{$geneID})){
      unless($#exons==-1){#not interested in single exon
        $transcript_gff_u{$geneID}=[@exons];
      }
    } 
    @exons=();
    $locID=substr($attributes[1],7);
    $geneID=substr($attributes[0],3);
    my $class_code="";
    my $protID="";
    foreach my $attr(@attributes){
      $class_code=substr($attr,11,1) if($attr =~ /^class_code=/);
      $protID=substr($attr,8) if($attr =~ /^cmp_ref=/);
    }
    if($class_code eq "k" || $class_code eq "=" || $class_code eq "j"){#equal intron chain or contains protein
      $transcript{$geneID}=$line;
      $transcript_cds{$geneID}=$protID;
      $transcripts_cds_loci{$locID}.="$geneID ";
    }elsif($class_code eq "u"){#no match to protein
      $transcript_u{$geneID}=$line;
      $transcripts_only_loci{$locID}.="$geneID ";
    }
  }elsif($gff_fields[2] eq "exon"){
    push(@exons,$line) if(defined($transcript{$geneID}) || defined($transcript_u{$geneID}));
  }
}
$transcript_gff{$geneID}=[@exons] if(not($geneID eq ""));
@exons=();

#process the loci
#print the header
print "##gff-version 3\n# EviAnn automated annotation\n";
for my $locus(keys %transcripts_cds_loci){
  my @output=();
  my $gene_feature="";
  my @transcripts_at_loci=split(/\s+/,$transcripts_cds_loci{$locus});
  my $locus_start=1000000000000;
  my $locus_end=0;
  my @gff_fields=split(/\t/,$transcript{$transcripts_at_loci[0]});
  my @attributes=split(";",$gff_fields[8]);
  $geneID=substr($attributes[1],7);#this is the XLOC
  my $parent=$geneID."-mRNA-";
  my $transcript_index=0;
  for my $t(@transcripts_at_loci){
    my @gff_fields_t=split(/\t/,$transcript{$t});
    my @attributes_t=split(";",$gff_fields_t[8]);
    my @cds_local=@{$protein_cds{$transcript_cds{$t}}};
    my $protID=$transcript_cds{$t};
    $used_proteins{$protID}=1;
    my @gff_fields_p=split(/\t/,$protein{$protID});
    my $start_cds=$gff_fields_p[3];
    my $end_cds=$gff_fields_p[4];
    my $transcript_start=$gff_fields_t[3];
    my $transcript_end=$gff_fields_t[4];
    my $transcript_cds_start_index=0;
    my $transcript_cds_end_index=$#{$transcript_gff{$t}};
    $locus_start=$gff_fields_t[3] if($gff_fields_t[3]<$locus_start);
    $locus_end=$gff_fields_t[4] if($gff_fields_t[4]>$locus_end);
    #here we fix CDS alignment problems
    #print "DEBUG $transcript{$t}\n";
    my $num_cds=-1;
    for (my $l=0;$l<=5 && not($num_cds==$#cds_local);$l++){
      $num_cds=$#cds_local;
      #print "$transcript_start $transcript_end $start_cds $end_cds\n";
      for(my $j=1;$j<=$#{$transcript_gff{$t}};$j++){
        my @gff_fields_curr=split(/\t/,${$transcript_gff{$t}}[$j]);
        my @gff_fields_prev=split(/\t/,${$transcript_gff{$t}}[$j-1]);
        my $num_cds=$#cds_local;
        if($start_cds>$gff_fields_prev[4] && $start_cds<$gff_fields_curr[3]){
        #we need to add another CDS entry to compensate for cds running into intron
          my $temp_cds=shift(@cds_local);#get the first cds to fix
          my @temp_cds_fields=split(/\t/,$temp_cds);
          $temp_cds_fields[3]=$gff_fields_curr[3];
          unshift(@cds_local,join("\t",@temp_cds_fields));
          $temp_cds_fields[4]=$gff_fields_prev[4];
          $temp_cds_fields[3]=$gff_fields_prev[4]-($gff_fields_curr[3]-$start_cds-1);
          unshift(@cds_local,join("\t",@temp_cds_fields));
        }
        if($end_cds>$gff_fields_prev[4] && $end_cds<$gff_fields_curr[3]){
#we need to add another CDS entry to compensate for cds running into intron
          my $temp_cds=pop(@cds_local);#get the first cds to fix
          my @temp_cds_fields=split(/\t/,$temp_cds);
          $temp_cds_fields[4]=$gff_fields_prev[4];
          push(@cds_local,join("\t",@temp_cds_fields));
          $temp_cds_fields[3]=$gff_fields_curr[3];
          $temp_cds_fields[4]=$gff_fields_curr[3]+($end_cds-$gff_fields_prev[4]-1);
          push(@cds_local,join("\t",@temp_cds_fields));
        }
      }
    $start_cds=(split(/\t/,@cds_local[0]))[3];
    $end_cds=(split(/\t/,@cds_local[-1]))[4];
    }
    my $start_cds_local=(split(/\t/,@cds_local[0]))[3];
    my $end_cds_local=(split(/\t/,@cds_local[-1]))[4];
#here we figure out which exons are UTR
    for(my $i=0;$i<=$#{$transcript_gff{$t}};$i++){
      my @gff_fields=split(/\t/,${$transcript_gff{$t}}[$i]);
      if($gff_fields[4]>=$start_cds_local && $gff_fields[3]<=$start_cds_local){
        $transcript_cds_start_index=$i;
        last;
      }
    }
    for(my $i=$#{$transcript_gff{$t}};$i>=0;$i--){
      my @gff_fields=split(/\t/,${$transcript_gff{$t}}[$i]);
      if($gff_fields[3]<=$end_cds_local && $gff_fields[4]>=$end_cds_local){
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
        $exon_start=$start_cds_local if($exon_start>$start_cds_local);
        $exon_end=$end_cds_local if($exon_end<$end_cds_local);
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
        push(@output,$gff_fields[0]."\tEviAnn\tfive_prime_UTR\t$gff_fields[3]\t".($start_cds_local-1)."\t".join("\t",@gff_fields[5..7])."\tID=$parent$transcript_index:5UTR:$i;Parent=$parent$transcript_index") if($start_cds_local>$gff_fields[3]);
      }else{
        push(@output,$gff_fields[0]."\tEviAnn\tfive_prime_UTR\t".join("\t",@gff_fields[3..7])."\tID=$parent$transcript_index:5UTR:$i;Parent=$parent$transcript_index");
      }
      $i++;
    }
#output cds from exons
    $i=1;
    if($transcript_cds_start_index<$transcript_cds_end_index){
      my @gff_fields=split(/\t/,${$transcript_gff{$t}}[$transcript_cds_start_index]);
      $i++;
      push(@output,$gff_fields[0]."\tEviAnn\tcds\t$start_cds_local\t".join("\t",@gff_fields[4..7])."\tID=$parent$transcript_index:cds:$i;Parent=$parent$transcript_index");
      for(my $j=$transcript_cds_start_index+1;$j<$transcript_cds_end_index;$j++){
        my @gff_fields=split(/\t/,${$transcript_gff{$t}}[$j]);
        push(@output,$gff_fields[0]."\tEviAnn\tcds\t".join("\t",@gff_fields[3..7])."\tID=$parent$transcript_index:cds:$i;Parent=$parent$transcript_index");
        $i++;
      }
      my @gff_fields=split(/\t/,${$transcript_gff{$t}}[$transcript_cds_end_index]);
      push(@output,$gff_fields[0]."\tEviAnn\tcds\t$gff_fields[3]\t$end_cds_local\t".join("\t",@gff_fields[5..7])."\tID=$parent$transcript_index:cds:$i;Parent=$parent$transcript_index");
    }else{#single exon
      my @gff_fields=split(/\t/,${$transcript_gff{$t}}[$transcript_cds_start_index]);
      push(@output,$gff_fields[0]."\tEviAnn\tcds\t$start_cds_local\t$end_cds_local\t".join("\t",@gff_fields[5..7])."\tID=$parent$transcript_index:cds:$i;Parent=$parent$transcript_index");
    }
    #output cds from protein alignment
    #for my $x(@cds_local){
    #  my @gff_fields=split(/\t/,$x);
    #  push(@output,$gff_fields[0]."\tEviAnn\tcds\t".join("\t",@gff_fields[3..7])."\tID=$parent$transcript_index:cds:$i;Parent=$parent$transcript_index");
    #  $i++;
    #}
#output 3'UTR
    $i=1;
    for(my $j=$transcript_cds_end_index;$j<=$last_j;$j++){
      my @gff_fields=split(/\t/,${$transcript_gff{$t}}[$j]);
      if($j==$transcript_cds_end_index){
        push(@output,$gff_fields[0]."\tEviAnn\tthree_prime_UTR\t".($end_cds_local+1)."\t".join("\t",@gff_fields[4..7])."\tID=$parent$transcript_index:3UTR:$i;Parent=$parent$transcript_index") if($end_cds_local<$gff_fields[4]);
      }else{
        push(@output,$gff_fields[0]."\tEviAnn\tthree_prime_UTR\t".join("\t",@gff_fields[3..7])."\tID=$parent$transcript_index:3UTR:$i;Parent=$parent$transcript_index");
      }
      $i++
    }
  }#end of transcripts loop
#now we know the locus start ans end, and we can output the gene record
  $gene_record{$gff_fields[0]." ".$locus_start}="$gff_fields[0]\tEviAnn\tgene\t$locus_start\t$locus_end\t".join("\t",@gff_fields[5..7])."\tID=$geneID;geneID=$geneID;type=protein_coding\n".join("\n",@output)."\n";
  push(@gene_records,$gff_fields[0]." ".$locus_start);
}

#finally output "intergenic" transcripts
#some of these may be completely messed up
for my $locus(keys %transcripts_only_loci){
  my @output=();
  my @transcripts_at_loci=split(/\s+/,$transcripts_only_loci{$locus});
  my @gff_fields=split(/\t/,$transcript_u{$transcripts_at_loci[0]});
  my $locus_start=1000000000000;
  my $locus_end=0;
  my $geneID=$locus."_lncRNA";
  my $parent=$geneID."-mRNA-";
  my $transcript_index=0;
  #here we first compute intron junction score as the number of distinct intron junctions over the number of total intron junctions
  my %distinct_intron_junctions=();
  my $total_intron_junctions=1;
  my $junction_score=0;
  if($#transcripts_at_loci>0){
    for my $t(@transcripts_at_loci){
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
  for my $t(@transcripts_at_loci){
    next if(not(defined($transcript_gff_u{$t})));
    my @gff_fields_t=split(/\t/,$transcript_u{$t});
    my @attributes_t=split(";",$gff_fields_t[8]);
    $locus_start=$gff_fields_t[3] if($gff_fields_t[3]<$locus_start);
    $locus_end=$gff_fields_t[4] if($gff_fields_t[4]>$locus_end);
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
    $gene_record{$gff_fields[0]." ".$locus_start}="$gff_fields[0]\tEviAnn\tgene\t$locus_start\t$locus_end\t".join("\t",@gff_fields[5..7])."\tID=$geneID;geneID=$geneID;type=lncRNA;junction_score=$junction_score;\n".join("\n",@output)."\n";
    push(@gene_records,$gff_fields[0]." ".$locus_start);
  }
}
#output unused proteins
foreach my $p(keys %protein){
  next if(defined($used_proteins{$p}));
  my @gff_fields_p=split(/\t/,$protein{$p});
  print STDERR "$gff_fields_p[0]\tEviAnn\t",join("\t",@gff_fields_p[2..$#gff_fields_p]),"\n";
  foreach my $cds(@{$protein_cds{$p}}){
    my @gff_fields_c=split(/\t/,$cds);
    print STDERR "$gff_fields_c[0]\tEviAnn\t",join("\t",@gff_fields_c[2..$#gff_fields_c]),"\n$gff_fields_c[0]\tEviAnn\texon\t",join("\t",@gff_fields_c[3..$#gff_fields_c]),"\n";
  }
}

#now we sort and output
my @gene_records_sorted=sort mysort @gene_records;
my %output;
foreach $g(@gene_records_sorted){
  if(not(defined($output{$g}))){
    print $gene_record{$g};
    $output{$g}=1;
  }
}

sub mysort{
my ($chroma,$coorda)=split(/\s+/,$a);
my ($chromb,$coordb)=split(/\s+/,$b);
return($chroma cmp $chromb || $coorda <=>$coordb);
}


  
