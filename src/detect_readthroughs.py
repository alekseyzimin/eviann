#!/usr/bin/env python

def check_readthroughs(tr):
    #for i in range(0,len(tr)): print("Check "+tr[i]+" exon "+str(startexon[tr[i]])+" "+str(endexon[tr[i]]))
    for skip in range(0,len(tr)):
        #print("Skipping "+str(skip))
        if skip==0: last_CDS_end=endCDS[tr[1]]
        else: last_CDS_end=endCDS[tr[0]]
        if skip==0: last_exon_end=endexon[tr[1]]
        else: last_exon_end=endexon[tr[0]]
        for i in range(0,len(tr)):
            if i==skip: continue
            #print(str(i)+" startexon "+str(startexon[tr[i]])+" last_exon_end "+str(last_exon_end)+" startCDS "+str(startCDS[tr[i]])+" last_CDS_end "+str(last_CDS_end))
            if startexon[tr[i]] > last_exon_end and startCDS[tr[i]] > last_CDS_end:
                #print("readthrough "+tr[skip])
                print(tr[skip])
                break
            if endCDS[tr[i]] > last_CDS_end: last_CDS_end=endCDS[tr[i]]
            if endexon[tr[i]] > last_exon_end: last_exon_end=endexon[tr[i]]

        
        

import sys
gene=""
endCDS={}
endexon={}
startCDS={}
startexon={}
transc=[]
for line in sys.stdin:
    if line[0] == "#": continue
    gff_fields=line.split("\t")
    gff_fields[3]=int(gff_fields[3])
    gff_fields[4]=int(gff_fields[4])
    if gff_fields[2] == "mRNA":
        attr=gff_fields[8].split(";");
        tid=attr[0][3:]
        gid=attr[1][6:-1]
        #print("Transcript id:"+tid+" Gene id:"+gid)
        endCDS.update({tid:0})
        endexon.update({tid:0})
        startCDS.update({tid:100000000000})
        startexon.update({tid:100000000000})
        if not gid==gene:
            if len(transc)>4: check_readthroughs(transc)
            transc.clear()
            gene=gid
        transc.append(tid)
    elif gff_fields[2] == "CDS":
        if startCDS[tid]>gff_fields[3]: startCDS[tid]=gff_fields[3]
        if endCDS[tid]<gff_fields[4]: endCDS[tid]=gff_fields[4]
    elif gff_fields[2] == "exon":
        if startexon[tid]>gff_fields[3]: startexon[tid]=gff_fields[3]
        if endexon[tid]<gff_fields[4]: endexon[tid]=gff_fields[4]
if len(transc)>4: check_readthroughs(transc)




