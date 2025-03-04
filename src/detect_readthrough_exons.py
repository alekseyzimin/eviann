#!/usr/bin/env python
#this code detects which exons in UTRs are possible readthroughs
import re
import sys

pbeg={}
pend={}
protein_gff=sys.argv[1]
with open(protein_gff,"r") as f:
    for line in f:
        if line[0] == "#": continue
        gff_fields=line.strip().split("\t")
        if len(gff_fields) !=9: continue 
        if gff_fields[2] == "gene":
            gff_attr=gff_fields[8].split(";")
            #print(gff_attr[0][3:])
            pbeg.update({gff_attr[0][3:]:int(gff_fields[3])})
            pend.update({gff_attr[0][3:]:int(gff_fields[4])})


for line in sys.stdin:
    if line[0] == "#": continue
    gff_fields=line.strip().split("\t")
    if len(gff_fields) !=9: continue
    gff_fields[3]=int(gff_fields[3])
    gff_fields[4]=int(gff_fields[4])
    if gff_fields[2] == "transcript":
        flag=False
        #print(gff_fields[8])
        x=re.search(r"class_code \"\S+\"",gff_fields[8])
        class_code=x.group()[12:-1]
        if re.search(class_code,"k=jcmnoiy") == None: continue
        x=re.search(r"transcript_id \"\S+\"",gff_fields[8])
        tid=x.group()[15:-1]
        x=re.search(r"cmp_ref \"\S+\"",gff_fields[8])
        protid=x.group()[9:-1]
        #print(tid+" "+protid+" "+class_code)
        tr=tid.split(".")
        transcript=".".join(tr[0:-1])+" "+tr[-1]
        #print(transcript)
        if ((tr[-1]=="5p" and pend[protid]<=gff_fields[4]) or (tr[-1]=="3p" and pbeg[protid]>=gff_fields[3])) and gff_fields[6] =="+": flag=True
        if ((tr[-1]=="5p" and pbeg[protid]>=gff_fields[3]) or (tr[-1]=="3p" and pend[protid]<=gff_fields[4])) and gff_fields[6] =="-": flag=True
    elif gff_fields[2] == "exon" and flag:
        if (((tr[-1]== "5p" and pend[protid]>gff_fields[3]) or (tr[-1] == "3p" and pbeg[protid]<gff_fields[4])) and gff_fields[6] == "+") or (((tr[-1]== "5p" and pbeg[protid]<gff_fields[4]) or (tr[-1] == "3p" and pend[protid]>gff_fields[3])) and gff_fields[6] == "-"):
            print(transcript+" "+str(gff_fields[3])+" "+str(gff_fields[4])+" "+gff_fields[6]+" "+str(pbeg[protid])+" "+str(pend[protid]))


