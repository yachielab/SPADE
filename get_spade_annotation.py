#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
import sys 
import codecs
from Bio import SeqIO

def get_spade_annotations(dir_name = ""):
    dirs = os.listdir(".")
    periodic_repeats = [] 
    spade_gbk_list = [f for f in dirs if "SPADE" in f]  
    for gbk in spade_gbk_list:
        for record in list(SeqIO.parse(codecs.open(gbk, 'r', 'utf-8', 'ignore'),"genbank")):
            gbk_id   = record.id
            CDS_list = []  
            for feature in record.features:
                if feature.type == "repeat_region":
                    if "note" in feature.qualifiers.keys():
                        if "SPADE" in feature.qualifiers["note"][0]:
                            start, end   = feature.location.start, feature.location.end
                            period       = feature.qualifiers["period"][0]
                            data_type    = feature.qualifiers["sequence_type"][0] 
                            repeat_motif = feature.qualifiers["unmasked_rpt_unit_seq"][0].replace(" ","_").upper()  
                            rptnum       = feature.qualifiers["rpt_num"][0]  
                            variable_repeat_motif = feature.qualifiers["rpt_unit_seq"][0].replace(" ","_").upper()  
                            periodicity = feature.qualifiers["periodicity_score"][0]  

                            periodic_repeats.append([gbk_id,str(start),str(end),period,periodicity,rptnum,data_type,repeat_motif,variable_repeat_motif,"N/A"])
                if feature.type == "CDS":
                    CDS_list.append(feature) 
           
            for CDS in CDS_list: 
                for repeat in periodic_repeats:
                    if int(repeat[1]) <= CDS.location.end and int(repeat[2]) >= CDS.location.start:
                        if repeat[-1] == "N/A":
                            repeat[-1] = CDS.qualifiers["product"][0]
                        else: 
                            repeat[-1] += ":" + CDS.qualifiers["product"][0]
    return periodic_repeats

def remove_overlap(periodic_repeats): 
    new_periodic_repeats  = []
    dna_s_e_list     = []
    protein_s_e_list = []
    for line in periodic_repeats: 
        start, end  = int(line[1]), int(line[2]) 
        seq_type    = line[-4] 
        if seq_type == "prot":
            protein_s_e_list.append([start, end]) 
            new_periodic_repeats.append(line) 
        else:
            dna_s_e_list.append([start,end,line]) 

    for d_set in dna_s_e_list: 
        flag = 0
        for p_set in protein_s_e_list:
            if d_set[0] < p_set[1] and d_set[1] > p_set[0]:
                flag = 1
            else: 
                pass 
        if flag == 0: 
           new_periodic_repeats.append(d_set[2]) 
        else: 
            pass 

    return new_periodic_repeats

def crispr_search(features, s, e):
    overlap = "None"
    for feature in features:
        if len(list(set(range(feature.location.start,feature.location.end)) & set(range(s,e)))) > 0:
            if feature.type == "repeat_region" and "rpt_family" in feature.qualifiers.keys():
                if feature.qualifiers["rpt_family"][0] == "CRISPR":
                    overlap = (1.0 * min([e,feature.location.end]) - max([s,feature.location.start])) / (max([e,feature.location.end]) - min([s,feature.location.start]))
                    break
    return overlap

if __name__ == "__main__":
    output = open("SPADE_annotations.tsv","w")
    output.write("\t".join(["NCBI RefSeqID","Start","End","Period","Periodicity","#Repeat unit","Sequence type","Repetitive motif","Repetitive motif masked variable positions","Overlapped CDS annotations"]) + "\n")
    for record in SeqIO.parse(sys.argv[1],"genbank"):
        periodic_repeats = get_spade_annotations()
        periodic_repeats = remove_overlap(periodic_repeats) 
        for row in periodic_repeats:
            output.write("\t".join(map(str,row)) + "\n")
    output.close()
