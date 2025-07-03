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
gffread -M --cluster-only -F --keep-genes $QRY > $QRYL.fix
awk -F '\t' '{if($7 != "?") print }' $REF | \
gffread -M -F -J --keep-genes -g $GENOME --ids <(perl -F'\t' -ane '{if($F[2] eq "lnc_RNA" || $F[2] eq "mRNA" || $F[2] eq "transcript" || $F[2] eq "primary_transcript"){unless($F[8] =~/pseudo=true/){@f=split(";",$F[8]);print substr($f[0],3),"\n"}}}' $REF |grep -v 'rna-DA397_mgp01') > $REFL.gene_transcript_lncRNA
gffread -C $REFL.gene_transcript_lncRNA | perl -F'\t' -ane '{unless($F[2] eq "exon"){if($F[2] eq "CDS"){$F[2]="exon"}print join("\t",@F);}}' > $REFL.gene_transcript_lncRNA.CDSasEXONS.gff
gffread -C $QRYL.fix| perl -F'\t' -ane '{unless($F[2] eq "exon"){if($F[2] eq "CDS"){$F[2]="exon"}print join("\t",@F);}}' > $QRYL.fix.CDSasEXONS.gff
gffcompare -T -r $REFL.gene_transcript_lncRNA $QRYL.fix -o transc 2>&1
gffcompare -T --strict-match -e 0 -r $REFL.gene_transcript_lncRNA.CDSasEXONS.gff $QRYL.fix.CDSasEXONS.gff -o cds 1>/dev/null 2>&1
paste \
  <(cat <(perl -F'\t' -ane '{if($F[8]=~/class_code "="/){if($F[8]=~/transcript_id "(\S+)"/){print "$1\n";}}}' transc.annotated.gtf) <(perl -F'\t' -ane '{if($F[8]=~/class_code "="/){if($F[8]=~/transcript_id "(\S+)"/){print "$1\n";}}}' cds.annotated.gtf ) |sort|uniq|wc -l) \
  <(cat <(perl -F'\t' -ane '{if($F[8]=~/class_code "="/){if($F[8]=~/ref_gene_id "(\S+)"/){print "$1\n";}}}' transc.annotated.gtf) <(perl -F'\t' -ane '{if($F[8]=~/class_code "="/){if($F[8]=~/ref_gene_id "(\S+)"/){print "$1\n";}}}'  cds.annotated.gtf ) |sort|uniq|wc -l) \
  <(cat transc.annotated.gtf |perl -F'\t' -ane '{if($F[8]=~/class_code/){print "$1\n";}}'  |wc -l) \
  <(gffread -M $QRYL.fix| awk -F '\t' '{if($3=="locus") print}' |wc -l) \
  <(awk -F'\t' '{if($3=="transcript" || $3=="mRNA" || $3=="lnc_RNA" || $3=="primary_transcript") print;}' $REFL.gene_transcript_lncRNA | wc -l) \
  <(awk -F'\t' '{if($3=="gene") print;}' $REFL.gene_transcript_lncRNA | wc -l) | \
  awk 'BEGIN{print "Transcript/CDS and gene evaluations:"}{print "Reference Transcripts: "$5" Reference genes: "$6;tsn=int($1/$5*1000+.5)/10;tpr=int($1/$3*1000+.5)/10;gsn=int($2/$6*1000+.5)/10;gpr=int($2/$4*1000+.5)/10;print "Correct transcripts/CDS: "$1" Total transcripts: "$3"\nCorrect genes: "$2" Total genes: "$4"\nTranscript Sensitivity: "tsn" Precision: "tpr" F1: "int(2*tsn*tpr/(tsn+tpr)*10+0.5)/10"\nGene Sensitivity: "gsn" Precision: "gpr" F1: "int(2*gsn*gpr/(gsn+gpr)*10+.5)/10}'
PID=$$
(gffread -y prot.faa.$PID -g $GENOME $QRYL.fix && gffread -y ref.faa.$PID -g $GENOME $REFL.gene_transcript_lncRNA && cat <(ufasta one ref.faa.$PID |grep -v '^>' |sort |uniq ) <(cat prot.faa.$PID|ufasta one| grep -v '^>' |sort -S 10% |uniq) |sort |uniq -d |wc -l;ufasta one prot.faa.$PID |grep -v '^>' |sort|uniq|wc -l ;ufasta one ref.faa.$PID |grep -v '^>' |sort|uniq|wc -l) | \
  perl -ane '{print $F[0]," "}END{print "\n"}' |\
  awk '{print "Protein evaluations:\nCorrect proteins: "$1" Total unique proteins: "$2;sn=int($1/$3*1000+0.5)/10;pr=int($1/$2*1000+0.5)/10;print "Protein Sensitivity: "sn" Precision: "pr" F1: "int(2*sn*pr/(sn+pr)*10+0.5)/10}'
rm -f prot.faa.$PID ref.faa.$PID
rm -f transc.{loci,tracking}
rm -f cds.{loci,tracking}
#rm $REFL.gene_transcript_lncRNA $REFL.gene_transcript_lncRNA.CDSasEXONS.gff $QRYL.fix $QRYL.fix.CDSasEXONS.gff  
cat transc.stats

