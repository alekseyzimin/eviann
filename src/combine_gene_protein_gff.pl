#!/usr/bin/env perl
#This code adds gene and CDS record for every locus in GFF file
my $geneID="";
my @exons=();
my %loci;
my %protein_cds;
my %transcript_gff;
my @gene_records;
my %gene_record;
my $protID="";
my $dir="";

#this is output of gffcompare -D -o protuniq ../GCF_000001735.4_TAIR10.1_genomic.fna.GCF_000001735.4_TAIR10.1_protein.faa.palign.gff
while(my $line=<STDIN>){#we just read in the whole file
  chomp($line);
  my @gff_fields=split(/\t/,$line);
  my @attributes=split(";",$gff_fields[8]);
  if($gff_fields[2] eq "transcript"){
    if(not($protID eq "")){
      #we need to fix the stop codon coord
      if($dir eq "+"){
        my $last_cds=pop(@exons);
        my @exon_fields=split(/\t/,$last_cds);
        push(@exons,join("\t",@exon_fields[0..3])."\t".($exon_fields[4]+3)."\t".join("\t",@exon_fields[5..$#exon_fields]));
      }else{
        my $first_cds=shift(@exons);
        my @exon_fields=split(/\t/,$first_cds);
        unshift(@exons,join("\t",@exon_fields[0..2])."\t".($exon_fields[3]-3)."\t".join("\t",@exon_fields[4..$#exon_fields]));
      }
      $protein_cds{$protID}=[@exons];
    }
    @exons=();
    $dir=$gff_fields[6];
    $geneID=substr($attributes[1],7);#this is the XLOC
    $protID=substr($attributes[0],3);#this is TCONS
    $protein{$protID}=$line;
    $loci{$geneID}.="$protID ";
  }else{
    push(@exons,$line);
  }
}
if(not($protID eq "")){
#we need to fix the stop codon coord
  if($dir eq "+"){
    my $last_cds=pop(@exons);
    my @exon_fields=split(/\t/,$last_cds);
    push(@exons,join("\t",@exon_fields[0..3])."\t".($exon_fields[4]+3)."\t".join("\t",@exon_fields[5..$#exon_fields]));
  }else{
    my $first_cds=shift(@exons);
    my @exon_fields=split(/\t/,$first_cds);
    unshift(@exons,join("\t",@exon_fields[0..2])."\t".($exon_fields[3]-3)."\t".join("\t",@exon_fields[4..$#exon_fields]));
  }
  $protein_cds{$protID}=[@exons];
}
@exons=();

#this is output of gffcompare -r protuniq.combined.gtf ../GCF_000001735.4_TAIR10.1_genomic.clean.fna.gtf -o protref
#here we load up all transcripts that matched proteins
open(FILE,$ARGV[0]);
while(my $line=<FILE>){
  chomp($line);
  my @gff_fields=split(/\t/,$line);
  my @attributes=split(";",$gff_fields[8]);
  if($gff_fields[2] eq "transcript"){
    $transcript_gff{$geneID}=[@exons] if(defined($transcript{$geneID}));
    if(defined($transcript_u{$geneID})){
      if($#exons==0){#not interested in single exon
        delete($transcript_gff_u{$geneID});
      }else{
        $transcript_gff_u{$geneID}=[@exons];
      }
    } 
    @exons=();
    $geneID=substr($attributes[0],3);
    if($attributes[6] eq "class_code=k" || $attributes[6] eq "class_code=="){#equal intron chain or contains protein
      $protID=substr($attributes[5],8);
      my @gff_fields_p=split(/\t/,$protein{$protID});
      if($gff_fields_p[4]-$gff_fields_p[3] > ($gff_fields[4]-$gff_fields[3])*0.25){#protein must match more than 25% of the transcript length
        $transcript{$geneID}=$line;
        $transcripts_for_prot{$protID}.="$geneID ";
      }
    }elsif($attributes[3] eq "class_code=u"){#no match to protein
      $lociID=substr($attributes[2],5);
      $transcript_u{$geneID}=$line;
      $transcripts_only_loci{$lociID}.="$geneID ";
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
for my $locus(keys %loci){
  my @output=();
  my $gene_feature="";
  my @proteins_at_loci=split(/\s+/,$loci{$locus});
  my @gff_fields=split(/\t/,$protein{$proteins_at_loci[0]});
  my @attributes=split(";",$gff_fields[8]);
  #figure out locus start and end, we have to go through all proteins and transcripts at this loci
  my $locus_start=1000000000000;
  my $locus_end=0;
  $geneID=substr($attributes[1],7);#this is the XLOC
  my $parent=$geneID."-mRNA-";
  my $transcript_index=0;
  for my $prot(@proteins_at_loci){
    my @gff_fields_p=split(/\t/,$protein{$prot});
    my @attributes_p=split(";",$gff_fields_p[8]);
    my $start_cds=(split(/\t/,${$protein_cds{$prot}}[0]))[3];
    my $end_cds=(split(/\t/,${$protein_cds{$prot}}[-1]))[4];
    $locus_start=$gff_fields_p[3] if($gff_fields_p[3]<$locus_start);
    $locus_end=$gff_fields_p[4] if($gff_fields_p[4]>$locus_end);
    if(defined($transcripts_for_prot{$prot})){
      my @transcripts=split(/\s+/,$transcripts_for_prot{$prot});
      #here we have to be careful.  I do not want to output exons that are all UTR. So first we read the CDSs for each transcript
      for my $t(@transcripts){
        my @gff_fields_t=split(/\t/,$transcript{$t});
        my @attributes_t=split(";",$gff_fields_t[8]);
        my $transcript_start=$gff_fields_t[3];
        my $transcript_end=$gff_fields_t[4];
        my $transcript_cds_start_index=0;
        my $transcript_cds_end_index=$#{$transcript_gff{$t}};
        #here we figure out which exons are UTR
        for(my $i=0;$i<=$#{$transcript_gff{$t}};$i++){
          my @gff_fields=split(/\t/,${$transcript_gff{$t}}[$i]);
          if($gff_fields[4]>$start_cds){
            $transcript_cds_start_index=$i;
            last;
          }
        }
        for(my $i=$#{$transcript_gff{$t}};$i>=0;$i--){
          my @gff_fields=split(/\t/,${$transcript_gff{$t}}[$i]);
          if($gff_fields[3]<$end_cds){
            $transcript_cds_end_index=$i;
            last;
          }
        }
        $locus_start=$transcript_start if($transcript_start<$locus_start);
        $locus_end=$transcript_end if($transcript_end>$locus_end);
        #output transcript
        $transcript_index++;
        push(@output,$gff_fields[0]."\tEviAnn\tmRNA\t$transcript_start\t$transcript_end\t".join("\t",@gff_fields_t[5..7])."\tID=$parent$transcript_index;Parent=$geneID;$attributes_p[2]");
        #output exons
        my $i=1;
        for my $x(@{$transcript_gff{$t}}){
          my @gff_fields=split(/\t/,$x);
          push(@output,$gff_fields[0]."\tEviAnn\texon\t".join("\t",@gff_fields[3..7])."\tID=$parent$transcript_index:exon:$i;Parent=$parent$transcript_index");
          $i++;
        }
        $i=1;
        #output 5'UTR
        for(my $j=0;$j<=$transcript_cds_start_index;$j++){
          my @gff_fields=split(/\t/,${$transcript_gff{$t}}[$j]);
          if($j==$transcript_cds_start_index){
            push(@output,$gff_fields[0]."\tEviAnn\tfive_prime_UTR\t$gff_fields[3]\t".($start_cds-1)."\t".join("\t",@gff_fields[5..7])."\tID=$parent$transcript_index:5UTR:$i;Parent=$parent$transcript_index") if($start_cds>$gff_fields[3]);
          }else{
            push(@output,$gff_fields[0]."\tEviAnn\tfive_prime_UTR\t".join("\t",@gff_fields[3..7])."\tID=$parent$transcript_index:5UTR:$i;Parent=$parent$transcript_index");
          }
          $i++;
        }
        #output cds
        $i=1;
        for my $x(@{$protein_cds{$prot}}){
          my @gff_fields=split(/\t/,$x);
          push(@output,$gff_fields[0]."\tEviAnn\tcds\t".join("\t",@gff_fields[3..7])."\tID=$parent$transcript_index:cds:$i;Parent=$parent$transcript_index");
          $i++;
        }
        #output 3'UTR
        $i=1;
        for(my $j=$transcript_cds_end_index;$j<=$#{$transcript_gff{$t}};$j++){
          my @gff_fields=split(/\t/,${$transcript_gff{$t}}[$j]);
          if($j==$transcript_cds_end_index){
            push(@output,$gff_fields[0]."\tEviAnn\tthree_prime_UTR\t".($end_cds+1)."\t".join("\t",@gff_fields[4..7])."\tID=$parent$transcript_index:3UTR:$i;Parent=$parent$transcript_index") if($end_cds<$gff_fields[4]);
          }else{
            push(@output,$gff_fields[0]."\tEviAnn\tthree_prime_UTR\t".join("\t",@gff_fields[3..7])."\tID=$parent$transcript_index:3UTR:$i;Parent=$parent$transcript_index");
          }
          $i++
        }
      }
    }else{
      $transcript_index++;
      push(@output,$gff_fields_p[0]."\tEviAnn\tmRNA\t".join("\t",@gff_fields_p[3..7])."\tID=$parent$transcript_index;Parent=$geneID;$attributes_p[2]");
      my $i=1;
      for my $x(@{$protein_cds{$prot}}){ 
        my @gff_fields=split(/\t/,$x);
        push(@output,$gff_fields[0]."\tEviAnn\texon\t".join("\t",@gff_fields[3..7])."\tID=$parent$transcript_index:exon:$i;Parent=$parent$transcript_index");
        $i++;
      }
      $i=1;
      for my $x(@{$protein_cds{$prot}}){
        my @gff_fields=split(/\t/,$x);
        push(@output,$gff_fields[0]."\tEviAnn\tcds\t".join("\t",@gff_fields[3..7])."\tID=$parent$transcript_index:cds:$i;Parent=$parent$transcript_index");
        $i++;
      }
    }
  }
#now we know the locus start ans end, and we can output the gene record
  $gene_record{$gff_fields[0]." ".$locus_start}="$gff_fields[0]\tEviAnn\tgene\t$locus_start\t$locus_end\t".join("\t",@gff_fields[5..7])."\tID=$geneID;geneID=$geneID;type=protein_coding\n".join("\n",@output)."\n";
  push(@gene_records,$gff_fields[0]." ".$locus_start);
}

#finally output "intergenic" transcripts
for my $locus(keys %transcripts_only_loci){
  my @output=();
  my @transcripts_at_loci=split(/\s+/,$transcripts_only_loci{$locus});
  my @gff_fields=split(/\t/,$transcript_u{$transcripts_at_loci[0]});
  my $locus_start=1000000000000;
  my $locus_end=0;
  my $geneID=$locus."_lncRNA";
  my $parent=$geneID."-mRNA-";
  my $transcript_index=0;
  for my $t(@transcripts_at_loci){
    my @gff_fields_t=split(/\t/,$transcript_u{$t});
    my @attributes_t=split(";",$gff_fields_t[8]);
    $locus_start=$gff_fields_t[3] if($gff_fields_t[3]<$locus_start);
    $locus_end=$gff_fields_t[4] if($gff_fields_t[4]>$locus_end);
    $transcript_index++;
    push(@output,join("\t",@gff_fields_t[0..1])."\tmRNA\t".join("\t",@gff_fields_t[3..7])."\tID=$parent$transcript_index;Parent=$geneID;$attributes_t[3]");
    my $i=1;
    for my $x(@{$transcript_gff_u{$t}}){
      my @gff_fields=split(/\t/,$x);
      push(@output,join("\t",@gff_fields[0..7])."\tID=$parent$transcript_index:exon:$i;Parent=$parent$transcript_index");
      $i++;
    }
  }
  $gene_record{$gff_fields[0]." ".$locus_start}="$gff_fields[0]\tEviAnn\tgene\t$locus_start\t$locus_end\t".join("\t",@gff_fields[5..7])."\tID=$geneID;geneID=$geneID;type=lncRNA\n".join("\n",@output)."\n";
  push(@gene_records,$gff_fields[0]." ".$locus_start);
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


  
