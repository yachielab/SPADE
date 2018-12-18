#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sys 
import numpy as np
import pycircos as pc
import matplotlib 
matplotlib.use("Agg") 
import matplotlib.pyplot as plt
from Bio import SeqIO

def remove_overlap(periodic_repeats): 
    new_periodic_repeats  = []
    dna_s_e_list          = []
    protein_s_e_list      = []
    for repeat in periodic_repeats: 
        start, end  = int(repeat[0]), int(repeat[1]) 
        seq_type    = repeat[2] 
        if seq_type == "prot":
            protein_s_e_list.append([start,end]) 
            new_periodic_repeats.append(repeat) 
        else:
            dna_s_e_list.append([start,end,repeat]) 

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

if __name__ == "__main__":
    record_parse = SeqIO.parse(sys.argv[1],"genbank")
    chuncos = pc.GENOME() #FFFECC
    chuncos.read_locus(record_parse,interspace=0, bottom=350, height=100, plot=True, requirement=lambda x: "" in x,circular=True, color_list=["#FFFECC"], features=True)
    kmer_count          = []
    dna_periodicity     = [] 
    protein_periodicity = [] 

    #read k-mer count
    for line in open(sys.argv[2]): 
        kmer_count.append(float(line.rstrip()))
    
    protein_periodicity = [0] * len(kmer_count) 
    dna_periodicity     = [0] * len(kmer_count) 
    
    periodic_repeats = [] 
    for key in chuncos.locus_dict.keys():
        for feat in chuncos.locus_dict[key]["features"]: 
            if feat.type == "repeat_region" and "note" in feat.qualifiers.keys() and "SPADE" in feat.qualifiers["note"][0]: 
                seq_type = feat.qualifiers["sequence_type"][0]
                periodic_repeats.append([feat.location.start, feat.location.end, seq_type]) 
    
    periodic_repeats = remove_overlap(periodic_repeats) 
    for repeat in periodic_repeats:
        start, end  = int(repeat[0]), int(repeat[1]) 
        seq_type    = repeat[2] 
        if seq_type == "prot":
            for i in range(start,end): 
                if kmer_count[i] > 0:
                    protein_periodicity[i] = kmer_count[i]  
                else:
                    protein_periodicity[i] = 1
        else:
            for i in range(start,end): 
                dna_periodicity[i] = kmer_count[i]

    for key in chuncos.locus_dict.keys():
        chuncos.plot_data(key, np.array(kmer_count), bottom=450, height=250, xaxes=False, yaxes=False, plot_style="fill", color1="#CCCCCC")
        chuncos.plot_data(key, np.array(dna_periodicity), bottom=450, height=250, xaxes=False, yaxes=False, plot_style="fill", color1="#531096", max_value=max(kmer_count))
        chuncos.plot_data(key, np.array(protein_periodicity), bottom=450, height=250, xaxes=False, yaxes=False, plot_style="fill", color1="#ff9300", max_value=max(kmer_count))

        #For CRISPR
        chuncos.plot_feature(feat_type="repeat_region", bottom=350, height=100, requirement=lambda x:("rpt_family" in x.qualifiers.keys() and "CRISPR" in x.qualifiers["rpt_family"][0]), color="#531096", expand=3, scatter=True)

        chuncos.plot_feature(feat_type="CDS", bottom=720, height=70, requirement=lambda x: x.strand==-1, color="#0494FF")
        chuncos.plot_feature(feat_type="CDS", bottom=800, height=70, requirement=lambda x: x.strand== 1, color="#FF7D76")
    
    chuncos.save(format_type="png",file_name=sys.argv[1].replace(".gb","_map")) 
