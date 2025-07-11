#!/bin/bash
#this code evaluates sensitivity and precision of genes and transcripts counting either a transcript or a CDS match as a match
function usage {
  echo "Usage: evaluate_gffcompare.sh [options]"
  echo "Options:"
  echo " -r FILE    reference annotation in GFF format"
  echo " -q FILE    query annotation in GFF format"
  echo " -g FILE    genome sequence file in fasta format"
}

while [[ $# > 0 ]]
do
    key="$1"

    case $key in
        -r|--reference)
            REF="$2"
            shift
            ;;
        -q|--query)
            QRY="$2"
            shift
            ;;
        -g|--genome)
            GENOME="$2"
            shift
            ;;
        -h|--help|-u|--usage)
            usage
            exit 255
            ;;
        *)
            echo "Unknown option $1"
            usage
            exit 1        # unknown option
            ;;
    esac
    shift
done

if [ ! -s $REF ] || [ ! -s $QRY ] || [ ! -s $GENOME ];then
  echo "ERROR: one of the input files is empty of does not exist"
  exit 1
fi

REF=$1
QRY=$2
GENOME=$3
REFL=`basename $REF`
QRYL=`basename $QRY`

#compute loci in the query
gffread -M --cluster-only -F --keep-genes $QRY > $QRYL.fix

#filter the refseq annotation keeping only protein coding and lnc_RNA genes
awk -F '\t' '{if($7 != "?") print }' $REF | \
  gffread -M -F -J -g $GENOME --ids <(perl -F'\t' -ane '{if(($F[2] eq "lnc_RNA" || $F[2] eq "mRNA" || $F[2] eq "transcript" || $F[2] eq "primary_transcript")){unless($F[8] =~/pseudo=true|exception=dicistronic gene|exception=trans-splicing|exception=RNA editing/){@f=split(";",$F[8]);print substr($f[0],3),"\n"}}}' $REF ) > $REFL.gene_transcript_lncRNA.gff

#produce table of locus gene transcript
perl -F'\t' -ane '{if($F[2] eq "lnc_RNA" || $F[2] eq "mRNA" || $F[2] eq "transcript" || $F[2] eq "primary_transcript"){@f=split(/;/,$F[8]);print "$locus ",substr($f[1],7)," ",substr($f[0],3),"\n";}elsif($F[2] eq "locus"){@f=split(/;/,$F[8]);$locus=substr($f[0],3);}}' $REFL.gene_transcript_lncRNA.gff > $REFL.locus_gene_transcript

perl -F'\t' -ane '{if($F[2] eq "lnc_RNA" || $F[2] eq "mRNA" || $F[2] eq "transcript" || $F[2] eq "primary_transcript"){@f=split(/;/,$F[8]);print "$gene ",substr($f[0],3),"\n";}elsif($F[2] eq "gene"){@f=split(/;/,$F[8]);$gene=substr($f[0],3);}}' $QRYL.fix > $QRYL.gene_transcript

#produce CDS-only reference
gffread -C $REFL.gene_transcript_lncRNA.gff | \
  perl -F'\t' -ane '{unless($F[2] eq "exon"){if($F[2] eq "CDS"){$F[2]="exon"}print join("\t",@F);}}' |\
  gffread -M > $REFL.gene_transcript_lncRNA.CDSasEXONS.gff

#produce CDS-only query
gffread -C $QRYL.fix|\
  perl -F'\t' -ane '{unless($F[2] eq "exon"){if($F[2] eq "CDS"){$F[2]="exon"}print join("\t",@F);}}' |\
  gffread -M > $QRYL.fix.CDSasEXONS.gff

#compute which transcript match
cat <(trmap -c '=' $REFL.gene_transcript_lncRNA.gff $QRYL.fix) |\
  perl -ane '{if($F[0]=~/^>/){$transc=substr($F[0],1)}else{$h{$F[5]}=$transc}}END{open(FILE,"'$REFL'.locus_gene_transcript");while($line=<FILE>){chomp($line); @f=split(/\s/,$line); print $line," $h{$f[2]}\n" if(defined $h{$f[2]});}}' > $REFL.matched_locus_gene_transcript

#compute which CDSs match
cat <(trmap -c '=' --strict-match $REFL.gene_transcript_lncRNA.CDSasEXONS.gff $QRYL.fix.CDSasEXONS.gff) |\
  perl -ane '{if($F[0]=~/^>/){$transc=substr($F[0],1)}else{$h{$F[5]}=$transc}}END{open(FILE,"'$REFL'.locus_gene_transcript");while($line=<FILE>){chomp($line); @f=split(/\s/,$line); print $line," $h{$f[2]}\n" if(defined $h{$f[2]});}}' > $REFL.matched_locus_gene_CDS

#produce matched qry loci
cat $REFL.matched_locus_gene_transcript $REFL.matched_locus_gene_CDS |\
  perl -ane '{$h{$F[3]}=1}END{open(FILE,"'$QRYL'.gene_transcript");while($line=<FILE>){chomp($line);@f=split(/\s/,$line);print $line,"\n" if(defined $h{$f[1]});}}' > $QRYL.matched_locus_transcriptORCDS

#compute statistcs
N_REF_TRANSCRIPTS=`awk '{print $3}' $REFL.locus_gene_transcript | wc -l`
N_REF_GENES=`awk '{print $2}' $REFL.locus_gene_transcript | sort|uniq |wc -l`
N_MATCHING_REF_TRANSCRIPTS=`awk '{print $3}' $REFL.matched_locus_gene_transcript |sort |uniq |wc -l`
N_MATCHING_REF_CDS=`awk '{print $3}' $REFL.matched_locus_gene_CDS |sort |uniq |wc -l`
N_MATCHING_REF_GENES=`cat $REFL.matched_locus_gene_transcript  $REFL.matched_locus_gene_CDS | awk '{print $2}' |sort |uniq |wc -l`
N_REF_CDS=`awk -F '\t' '{if($3=="transcript" ||$3=="mRNA"|| $3=="lnc_RNA" || $3=="primary_transcript") print}' $REFL.gene_transcript_lncRNA.CDSasEXONS.gff |wc -l`

N_QRY_TRANSCRIPTS=`awk -F '\t' '{if($3=="transcript" ||$3=="mRNA"|| $3=="lnc_RNA") print}' $QRYL.fix  |wc -l`
N_QRY_LOCI=`awk '{print $1}' $QRYL.gene_transcript |sort|uniq|wc -l`
N_MATCHING_QRY_LOCI=`awk '{print $1}' $QRYL.matched_locus_transcriptORCDS |sort| uniq | wc -l`
N_MATCHING_QRY_TRANSCRIPTS=`awk '{print $4}' $REFL.matched_locus_gene_transcript |sort |uniq |wc -l`
N_MATCHING_QRY_CDS=`awk '{print $4}' $REFL.matched_locus_gene_CDS |sort |uniq |wc -l`
N_QRY_CDS=`awk -F '\t' '{if($3=="transcript" ||$3=="mRNA"|| $3=="lnc_RNA") print}' $QRYL.fix.CDSasEXONS.gff |wc -l`

#compute Sn Pr and F1 and print
awk 'BEGIN{print "Transcript/CDS and gene evaluations:";
  print "Reference Transcripts: '$N_REF_TRANSCRIPTS' Reference genes: '$N_REF_GENES'";
  tsn=int('$N_MATCHING_REF_TRANSCRIPTS'/'$N_REF_TRANSCRIPTS'*1000+.5)/10;
  tpr=int('$N_MATCHING_QRY_TRANSCRIPTS'/'$N_QRY_TRANSCRIPTS'*1000+.5)/10;
  csn=int('$N_MATCHING_REF_CDS'/'$N_REF_CDS'*1000+.5)/10;
  cpr=int('$N_MATCHING_QRY_CDS'/'$N_QRY_CDS'*1000+.5)/10;
  gsn=int('$N_MATCHING_REF_GENES'/'$N_REF_GENES'*1000+.5)/10;
  gpr=int('$N_MATCHING_QRY_LOCI'/'$N_QRY_LOCI'*1000+.5)/10;

