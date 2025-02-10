#!/bin/bash
perl -F'\t' -ane '{if($F[2] eq "tRNA" && $F[5]>=57.3){chomp($F[8]);@f=split(/;/,$F[8]);$id=substr($f[0],3);$F[2]="gene";print join("\t",@F),";gene_biotype=tRNA\n";$F[2]="tRNA";print join("\t",@F[0..7]),"\tID=tRNA-$id;Parent=$id\n";$F[2]="exon";print join("\t",@F[0..7]),"\tID=tRNA-exon-$id;Parent=tRNA-$id\n";}}' 
