#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
import csv
import sys
import copy 
import time
import pickle
import argparse
import math
import shutil
import collections
import subprocess

__version__ = "1.0.0"

def savetxt(file_name,data,delimiter="\t",fmt=":.0f",header=""):
    if len(data.shape) == 1:
        with open(file_name,"w") as o:
            o.write(header) 
            for datum in data:
                line = "{" + fmt + "}"
                o.write(line.format(datum) + "\n") 
    else:
        with open(file_name,"w") as o:
            o.write(header) 
            for datum in data:
                line = delimiter.join(["{" + fmt + "}" for _ in range(len(datum))])
                o.write(line.format(*datum) + "\n") 

class SPADE(object):
    def __init__(self): 
        self.ps   = [] 
        self.num_threads = 1
        self.f_parse = None
        self.g_parse = None
        self.format_type = None

    def load(self, file_name):
        if file_name.split(".")[-1] == "fasta" or file_name.split(".")[-1] == "fna" or file_name.split(".")[-1] == "faa" or args.f=="fasta":
            self.f_parse   = SeqIO.parse(args.input,"fasta")
            self.format_type = "fasta"

        elif file_name.split(".")[-1] == "gbff" or file_name.split(".")[-1] == "gb" or args.f=="genbank":
            self.g_parse   = SeqIO.parse(args.input,"genbank")
            self.format_type = "genbank"
        
        else:
            try:
                self.g_parse = SeqIO.parse(args.input,"genbank")
            except: 
                try:
                    self.f_parse = SeqIO.parse(args.input,"fasta")
                except:
                    raise ValueError("Please input genbank or fasta format file")
    
    def run(self, num_threads, parameters):
        if self.format_type == "fasta":
            record_parse = self.f_parse
        else:
            record_parse = self.g_parse
        
        ps = [] 
        finished = 0 
        
        #Single process mode
        for record in record_parse:
            if self.g_parse != None and self.f_parse != None:
                record2 = self.f_parse.next()
                if len(set(str(record2.seq))) == 1:
                    seq = record2.seq
                    seq.Alphabet = Alphabet.DNAAlphabet()
                    record.seq = seq
        
            #Single process mode
            if num_threads < 2:
                if os.path.exists("./" + record.id) == False:
                    os.mkdir(record.id) 
                else:
                    pass 
                os.chdir(record.id)
                locus = LOCUS(str(record.seq),record,*(parameters+[self.format_type]))
                locus.all()
                os.chdir("../")

            #Multiprocess mode 
            else:
                ps.append(mp.Process(target=self.thread_run, args=(record,parameters)))     
            finished += 1
    
        if num_threads > 1:
            liveProcess  = ps[0:num_threads]
            finishedList = []
            presentindex = 0
            for p in liveProcess:
                p.start()
                presentindex += 1

            while 1:
                time.sleep(1)
                for i, p in enumerate(liveProcess):
                    if p.is_alive() or i in finishedList:
                        pass
                    else:
                        finishedList.append(i)
                        if presentindex < len(ps):
                            liveProcess.append(ps[presentindex])
                            liveProcess[-1].start()
                            presentindex += 1
                if len(finishedList) == len(liveProcess):
                    break	

    def thread_run(self,record,parameters):
        #Single process for multiprocess mode.
        if os.path.exists("./" + record.id) == False:
            os.mkdir(record.id) 
        else:
            pass 
        os.chdir(record.id)
        parameters.append(self.format_type)
        locus = LOCUS(str(record.seq),record,*parameters)
        locus.all()
        os.chdir("../")


class LOCUS(object):
    #Nm, Nr,  Nq
    def __init__(self, seq, record, k_size, w_size, g_size, tk, tp, tb, tm, tr, tq, pk_size, pw_size, pg_size, ptk, ptp, ptb, ptm, ptr, ptq, visualisation, delete, option_mafft, option_blastn, option_blastp, format_type):
        self.HRA_list      = []
        self.Id            = record.id 
        self.seq           = seq
        self.record        = record
        self.featnum       = 0
        self.k_size        = k_size
        self.w_size        = w_size 
        self.g_size        = g_size
        self.tm            = tm
        self.tr            = tr 
        self.tq            = tq
        self.tk            = tk
        self.tp            = tp
        self.tb            = tb
        self.pk_size       = pk_size 
        self.pw_size       = pw_size
        self.pg_size       = pg_size 
        self.ptk           = ptk
        self.ptp           = ptp
        self.ptb           = ptb
        self.ptm           = ptm
        self.ptq           = ptq
        self.ptr           = ptr
        self.option_mafft  = option_mafft 
        self.option_blastn = option_blastn
        self.option_blastp = option_blastp
        self.format_type   = format_type
        self.process       = "None"
        self.seqtype       = "nucl"
        self.visualisation = visualisation
        self.delete        = delete

        head_seq         = str(self.seq)[0:10000]
        ratio = (head_seq.count("A") + head_seq.count("T") + head_seq.count("G") + head_seq.count("C") + head_seq.count("N")) * 1.0 / len(head_seq)
        if ratio > 0.9 or self.seqtype == "nucl":
            self.seqtype = "nucl"
        else:
            self.seqtype = "prot"
            self.k_size  = self.pk_size
            self.w_size  = self.pw_size
            self.g_size  = self.pg_size
            self.tk      = self.ptk
            self.tp      = self.ptp
            self.tb      = self.ptb
            self.tm      = self.ptm
            self.tq      = self.ptq
            self.tr      = self.ptr
        #print self.seqtype
    
    def cumulative_kmer_count(self): 
        self.score_array, self.hra_range_list = kmer_count(self.seq, self.k_size, self.w_size, thresh=self.tk, gap=self.g_size, buf=self.tm, seqtype=self.seqtype)
        savetxt(self.Id + "_kmer_count.txt", self.score_array, delimiter="\t", fmt=":.0f") 
        #print("") 
    
    def find_protein_HRA(self): 
        i = 0
        output_dir_list = []
        for feat in self.record.features:
            if feat.type == "CDS" and "translation" in feat.qualifiers.keys():
                hra = HRA(i, "prot", self.pk_size, self.pw_size, self.pg_size, self.ptk, self.ptp, self.ptb, self.ptm, self.ptq, self.ptr, self.option_mafft, self.option_blastn, self.option_blastp)
                region_score_array, phra_range_list = kmer_count(feat.qualifiers["translation"][0], self.pk_size, self.pw_size, thresh=self.ptk, gap=self.pg_size, buf=self.pw_size,seqtype="prot")
                if len(phra_range_list) < 1:
                    pass 
                else:
                    for j in range(len(phra_range_list)):
                        hra = HRA(i+j, "prot", self.pk_size, self.pw_size, self.pg_size, self.ptk, self.ptp, self.ptb, self.ptm, self.ptq, self.ptr, self.option_mafft, self.option_blastn, self.option_blastp)
                        hra.hra_range          = phra_range_list[j]  
                        hra.region_score_array = region_score_array[hra.hra_range[2]:hra.hra_range[3]]
                        hra.feature            = feat
                        hra.strand             = feat.strand 
                        hra.region_seq         = str(feat.qualifiers["translation"][0][hra.hra_range[2]:hra.hra_range[3]]).upper()
                        if feat.location.start == 1 and feat.location.end == hra.seq_len:
                            hra.p_start = feat.location.parts[0].start.position + 3 * hra.hra_range[0]
                            hra.p_end   = feat.location.parts[0].start.position + 3 * hra.hra_range[1] 
                                
                        else:
                            hra.p_start     = feat.location.start.position + 3 * hra.hra_range[0]         
                            hra.p_end       = feat.location.start.position + 3 * hra.hra_range[1]
                        
                        hra.output_dir  = "_".join([hra.dtype,str(hra.p_start),str(hra.p_end)])
                        afeat           = copy.deepcopy(feat) 
                        afeat.location  = FeatureLocation(0,len(feat.qualifiers["translation"][0]))
                        hra.features    = [afeat]
                        hra.record      = SeqRecord(Seq(str(feat.qualifiers["translation"][0]), Alphabet.ProteinAlphabet())) 
                        hra.record.features.append(afeat) 
                        hra.w_size      = hra.hra_range[1] - hra.hra_range[0]  
                        if hra.output_dir not in output_dir_list:
                            output_dir_list.append(hra.output_dir) 
                            self.HRA_list.append(hra)
                        else:
                            pass
                i += 1

    def find_HRA(self):
        HRA.window_size = self.w_size
        HRA.seq_len = len(self.seq) 
        for i, hra_range in enumerate(self.hra_range_list):
            if self.seqtype == "nucl":
                hra = HRA(i, "nucl", self.k_size, self.w_size, self.g_size, self.tk, self.tp, self.tb, self.tm, self.tq, self.tr, self.option_mafft, self.option_blastn, self.option_blastp)
            else:
                hra = HRA(i, "prot", self.pk_size, self.pw_size, self.pg_size, self.ptk, self.ptp, self.ptb, self.ptm, self.ptq, self.ptr, self.option_mafft, self.option_blastn, self.option_blastp)
            hra.record, hra.features  = hra.extract(self.record,hra_range[2],hra_range[3])  
            hra.region_seq = str(hra.record.seq).upper() 
            hra.region_score_array = self.score_array[hra_range[2]:hra_range[3]]
            hra.output_dir = "_".join([hra.dtype,str(hra_range[0]),str(hra_range[1])])
            hra.hra_range  = hra_range
            hra.w_size = hra_range[1] - hra_range[0]
            if hra.w_size < self.w_size*0.1: 
                hra.w_size = self.w_size
            self.HRA_list.append(hra)

        for hra in self.HRA_list:
            path = os.getcwd()
            if self.process == "find_HRA":
                record_handle = open(hra.output_dir + "/" + hra.output_dir + ".gb","w")  
                SeqIO.write(hra.record, record_handle, "genbank") 
                record_handle.close() 
                subject_name = hra.output_dir + "/subject.fasta"
                subject = open(subject_name,"w")
                subject.write(">subject\n")
                subject.write(hra.region_seq + "\n")
                subject.close()                  
    
    def search_seed_all(self):
        new_hra_list = []
        rm_dir_list = []
        it = 0 
        for hra in self.HRA_list:
            path = os.getcwd()
            try:
                hra.kmer_period_matrix() 
                hra.search_peek()
                hra.search_seed()
                #hra.output_data()
                if len(hra.peak_period_set) > 0:
                    new_hra_list.append(hra)
                else:
                    del hra.period_matrix, hra.kmer_position_dict, hra.kmer_period_dict, hra.unit_seq_list, hra.sorted_kmer_list, hra.sumx_matrix
                    rm_dir_list.append(hra.output_dir)
            
            except Exception as e:
                print("Error in search seed all. dir", self.Id, hra.output_dir) 
                print(e) 
            it += 1
        
        del self.HRA_list
        self.HRA_list = new_hra_list
        new_hra_list = [] 
        for hra in self.HRA_list: 
            path = os.getcwd()
            if os.path.exists("./" + hra.output_dir) == False:
                os.mkdir(hra.output_dir) 
            else: 
                pass
            os.chdir(hra.output_dir)            
            hra.output_data()
            del hra.kmer_position_dict, hra.kmer_period_dict, hra.sorted_kmer_list
            os.chdir(path)
            new_hra_list.append(hra) 
        self.HRA_list = new_hra_list

    def mafft_all(self):
        """
        Run Mafft commands for all HRAs
        """
        mafft_coms = open("mafft_coms.sh","w") 
        for hra in self.HRA_list: 
            if 1:
                hra.mafft()
                mafft_coms.write(hra.mafft_com + "\n")
            else:
                print("Error in mafft_all. dir", self.Id, hra.output_dir)
                print(e)
        mafft_coms.close() 
        subprocess.call("bash mafft_coms.sh",shell=True)

    def decide_query_all(self):
        for hra in self.HRA_list:
            path = os.getcwd()
            try:
                os.chdir(hra.output_dir) 
                hra.decide_query() 
            except Exception as e:
                print("Error in decide_query_all. dir", self.Id, hra.output_dir)  
                print(e)
            os.chdir(path)
    
    def blast_all(self): 
        """
        Run Mafft commands for all HRAs
        """
        blast_coms = open("blast_coms.sh","w") 
        for hra in self.HRA_list: 
            try:
                hra.blast()
                blast_coms.write(hra.blast_com + "\n")
            except Exception as e:
                print("Error in blast_all. dir", self.Id, hra.output_dir) 
                print(e)
        blast_coms.close() 
        subprocess.call("bash blast_coms.sh",shell=True)

    def make_se_sets_all(self):
        new_hra_list = []
        rm_dir_list  = [] 
        for hra in self.HRA_list:
            path = os.getcwd()
            try:
                os.chdir(hra.output_dir) 
                hra.make_se_sets()
                hra.make_motif_array() 

                if len(hra.peak_period_set) >0:
                    new_hra_list.append(hra)
                else:
                    rm_dir_list.append(os.getcwd())
            except Exception as e:
                print("Error in make_se_sets_all. dir", self.Id, hra.output_dir)  
                print(e) 
            os.chdir(path) 
        
        self.HRA_list = new_hra_list

    def make_feature_all(self, make_file=True):
        spade_list    = [] 
        rm_dir_list   = [] 
        start_end_set = []
        for hra in self.HRA_list:
            path = os.getcwd()
            if 1:
            #try:
                os.chdir(hra.output_dir) 
                feats = hra.make_feature()
                if len(feats) > 0:
                    for feat in feats:
                        if [feat.location.start, feat.location.end, feat.qualifiers["periodicity_score"][0]] not in start_end_set:
                            spade_list.append(feat)
                            start_end_set.append([feat.location.start, feat.location.end, feat.qualifiers["periodicity_score"][0]])  
                        else:
                            pass 
                else:
                    rm_dir_list.append(os.getcwd())
                keys = list(hra.__dict__.keys())
                for key in keys:
                    if key != "thresh" and key != "peak_period_set" and key != "dtype" and key != "output_dir" and key != "k_size" and key != "strand":
                        del hra.__dict__[key]
            else:
            #except Exception as e:
                print("Error in make_feature_all. dir", self.Id, hra.output_dir) 
                print(e) 
            os.chdir(path)    
        
        self.record.features.extend(spade_list)
        self.record.features.sort(key=lambda x: x.location.start)
        if make_file == True:
            if self.format_type == "genbank":
                record_handle = open(self.record.id + "_SPADE.gb", "w") 
                SeqIO.write(self.record, record_handle, "genbank")
            
            if self.format_type == "fasta":
                record_handle = open(self.record.id + "_SPADE.gb", "w") 
                new_record    = SeqRecord(Seq(str(self.record.seq),Alphabet.DNAAlphabet()))
                new_record.id = self.record.id[0:16]
                self.record   = new_record
                self.record.features.extend(spade_list)
                self.record.features.sort(key=lambda x: x.location.start)
                SeqIO.write(self.record, record_handle, "genbank")
            record_handle.close()
        if self.delete == 1: 
            for apath in rm_dir_list:
                print(apath) 
                shutil.rmtree(apath)

        del self.record.features
        del self.record

    def visualisation_all(self):
        """
        make figure of each HRR and motif
        """
        for hra in self.HRA_list:
            path = os.getcwd()
            try:
                os.chdir(hra.output_dir) 
                hra.make_figure()
                hra.make_motif_logo()
            except Exception as e:
                print("Error in visualisation_all. dir", self.Id, hra.output_dir) 
                print(e)
            os.chdir(path)

    def all(self):
        self.cumulative_kmer_count()
        if self.process == "kmer_count":
            return 1
        else:
            self.find_HRA()
            self.find_protein_HRA() 
            self.search_seed_all()
            if self.process == "find_HRA":
                return 1     
            else:
                self.mafft_all()
                self.decide_query_all()
                self.blast_all()
                self.make_se_sets_all()
                self.make_feature_all()
                if self.visualisation == "Y":
                    self.visualisation_all() 

class HRA(object):
    record      = None
    seq_len     = None
    window_size = 1000
    def __init__(self, Id, dtype, ksize, wsize, gsize, tk, tp, tb, tm, tq, tr, optionm, optionbn, optionbp): 
        self.Id                 = Id
        self.output_dir         = ""
        self.hra_range          = []
        self.period_matrix      = []
        self.kmer_position_dict = {} 
        self.kmer_period_dict   = {} 
        self.region_score_array = [] 
        self.region_seq         = [] 
        self.peak_period_set    = [] 
        self.periodicity        = [] 
        self.dtype         = dtype 
        self.k_size        = ksize
        self.w_size        = wsize
        self.g_size        = gsize
        self.tk            = tk
        self.tp            = tp
        self.tb            = tb
        self.tm            = tm 
        self.tq            = tq 
        self.tr            = tr
        self.option_mafft  = optionm
        self.option_blastn = optionbn
        self.option_blastp = optionbp
        self.p_start       = 0
        self.p_end         = 0
        self.mattf_com     = ""
        self.blast_com     = ""
        self.feature       = None
        self.strand        = None

    def kmer_period_matrix(self):
        self.period_matrix, self.region_score_array, self.kmer_position_dict  = kmer_count_matrix(self.region_seq, self.k_size, self.w_size, seqtype=self.dtype)
 
    def search_peek(self):
        for kmer in self.kmer_position_dict:
            if len(self.kmer_position_dict[kmer][1]) > 1 and kmer != "":
                self.kmer_period_dict[kmer] = []
                self.kmer_period_dict[kmer].extend(list(map(lambda x,y: x - y, self.kmer_position_dict[kmer][1][1:], self.kmer_position_dict[kmer][1][0:-1])))
        
        #Position-Periodicity matrix integration with x axis  
        self.sumx_matrix      = np.sum(self.period_matrix[:,self.hra_range[0]-self.hra_range[2]:self.hra_range[1]-self.hra_range[2]], axis=1) / (self.hra_range[1]-self.hra_range[0]) 

        #Peak period detection
        maxIds = signal.argrelmax(self.sumx_matrix)[0].tolist()
        maxIds.append(np.argmax(self.sumx_matrix))
        maxIds = [Id for Id in maxIds if Id < HRA.window_size]
        maxIds = list(set(maxIds))
        maxIds.sort()
        maxIds.sort(key=lambda x: -1.0 * self.sumx_matrix[x])  
        if len(maxIds) > 0:
            maxId  = maxIds[0] 
            #Present version of SPADE care only the higherst peak. We will care and evaluate periodicity of oter peaks at Next verson's SPADE.
            #maxIds = [Id for Id in maxIds if self.sumx_matrix[Id] > 0.01 * self.sumx_matrix[maxId]]
            maxIds = [Id for Id in maxIds if Id > maxId and Id % maxId == 0] 
            maxIds = [maxId] + maxIds 
        self.peak_period_set = maxIds

        #sorting with first postion of each k-mer set
        self.sorted_kmer_list = list(self.kmer_period_dict.keys()) 
        self.sorted_kmer_list.sort(key=lambda x: self.kmer_position_dict[x][1][0]) 
         
    def check_periodicity(self,peak,limit,start,end):
        region_sum = np.sum(self.region_score_array[start:end])
        if region_sum == 0: 
            self.periodicity = 0 
        else:
            self.periodicity = np.sum(self.period_matrix[peak-int(round(0.2*peak)):peak+int(round(0.2*peak))+1,start:end]) * 1.0 /region_sum 
        if self.periodicity >= limit:
            return self.periodicity   
        else: 
            return False

    def search_seed(self):
        self.unit_poss_list   = []
        self.unit_seq_list    = [] 
        new_peak_period_set   = []
        self.seed_kmer_list   = []
        self.max_r_list       = [] 
        for period in self.peak_period_set[:1]:
            if self.check_periodicity(period,0.2,self.hra_range[0]-self.hra_range[2],self.hra_range[1]-self.hra_range[2]): 
                max_count  = 0
                seed_kmer  = ""
                for kmer in self.sorted_kmer_list:
                    count = 0
                    s_is = []
                    for j, p in enumerate(self.kmer_period_dict[kmer]): 
                        if abs(p-period) < 0.2 * period and period in self.kmer_period_dict[kmer]:
                            count += 1
                            s_is.append(j) 
                        else:
                            pass 

                    if count > max_count:
                        max_count       = count
                        seed_kmer       = kmer
                        seed_indexs     = s_is
                    else:
                        pass
                if max_count > 1 and seed_kmer not in self.seed_kmer_list:
                    seed_poss = [] 
                    for index in seed_indexs:
                        seed_poss.append(self.kmer_position_dict[seed_kmer][1][index])
                    seed_poss.append(self.kmer_position_dict[seed_kmer][1][index + 1])
                    self.seed_kmer_list.append(seed_kmer) 
                    unit_seqs, seed_poss = self.get_unit_seq_list(seed_poss, period)  
                    self.unit_seq_list.append(unit_seqs)
                    self.unit_poss_list.append(seed_poss)
                    new_peak_period_set.append(period)
                else:
                    self.periodicity = False 
        self.peak_period_set = new_peak_period_set
    
    def get_unit_seq_list(self, seed_poss, period):
        #Extraction of repeat unit sequence candidates based on Position-Periodicity Matrix
        r_list = []
        blank  =  0
        new_seed_poss = []
        seed_poss.sort()
        for seed_pos in seed_poss:
            r = 0
            blank = 0 
            while blank < 6  and  r < period  and seed_pos - r > 0:
                r += 1
                if self.region_score_array[seed_pos-r] < 2: 
                    blank += 1
                else:
                    blank = 0
            if r-blank == period:
                r_list.append(0)
            else:
                r_list.append(r-blank)

        max_r = max(r_list) 
        unit_seq_list = []
        for i,seed_pos in enumerate(seed_poss):
            r_num = seed_pos - max_r
            if r_num < 0:
                pass 
            #    r_num = seed_pos
            #    gaps = ""
            #    for j in range(abs(r_num)):
            #        gaps += "-"
            #    unit_seq = gaps + self.region_seq[0:period-abs(r_num)]
            else:
                if i < len(seed_poss) - 1 and (seed_poss[i+1] - seed_pos < 1.2 * period and seed_poss[i+1] - seed_pos > 0.8 * period):
                    unit_seq = self.region_seq[r_num:r_num + seed_poss[i+1] - seed_pos]
                else:
                    unit_seq = self.region_seq[r_num:r_num + period]
        
                new_seed_poss.append(r_num)
                unit_seq_list.append(unit_seq)
        
        self.max_r_list.append(max_r) 
        return unit_seq_list, new_seed_poss  

    def output_data(self, clear=1):
        #Data output for vidualization. 
        period = self.peak_period_set[0] 
        p_min = period - max(int(round(0.2*period)),5) if period - max(int(round(0.2*period)),5) > 0 else 0
        p_max = period + max(int(round(0.2*period)),5) if period + max(int(round(0.2*period)),5) < 1000 else 1000
        if p_max >= self.period_matrix.shape[0]: 
            p_max = self.period_matrix.shape[0]-1
        
        f = "ppm4vis.tsv"
        header = "\tPositoin\nPeriod\t"
        header += "\t".join(list(map(str,list(range(1,self.period_matrix.shape[1]+1))))) + "\n"
        savetxt(f, np.concatenate((np.array([list(range(p_min,p_max+1))]).T, self.period_matrix[p_min:p_max+1,:].astype(np.int)), axis=1), delimiter="\t", fmt=":.0f", header=header)
        #####
        #        Position
        #Pediod     1       2       3       4       5       6   
        #1          score   score   score   score   score   score
        #2          score   score   score   score   score   score
        #3          score   score   score   score   socre   score
        #####
        
        f = "pdist.tsv"
        header="Period\tPeriodicity\n"
        savetxt(f, np.concatenate((np.array([list(range(p_min,p_max+1))]).T,np.array([self.sumx_matrix[p_min:p_max+1]]).T),axis=1), delimiter="\t", fmt=":.2f", header=header)
        #####
        #Period Periodicity   
        #1        score  
        #2        score  
        #3        score  
        #####

        f = "ppm.tsv"
        header = "\tPositoin\nPeriod\t"
        header += "\t".join(list(map(str,list(range(1,self.period_matrix.shape[1]+1))))) + "\n"
        savetxt(f, np.concatenate((np.array([list(range(1,self.period_matrix.shape[0]+1))]).T, self.period_matrix.astype(np.int)), axis=1), delimiter="\t", fmt=":.0f", header=header)
        #####
        #        Position
        #Pediod     1       2       3       4       5       6   
        #1          score   score   score   score   score   score
        #2          score   score   score   score   score   score
        #3          score   score   score   score   socre   score
        #####
        
        f = "kmer.tsv"
        header = "Position\tScore\n"
        data   = [] 
        for n, score in enumerate(self.region_score_array):
            data.append([n+1,score]) 
        savetxt(f, np.array(data), delimiter="\t", fmt=":.0f", header=header) 
        #####
        #Position Score   
        #1        score  
        #2        score  
        #3        score  
        #####

        #Data output for motif detection
        subject_name = "subject.fasta"
        subject = open(subject_name,"w")
        subject.write(">subject\n")
        subject.write(self.region_seq + "\n")
        subject.close()  
        for i,period in enumerate(self.peak_period_set):
            fasta_name =  "unit_seq.fasta"
            fasta = open(fasta_name,"w") 
            for j, unit_seq in enumerate(self.unit_seq_list[i]):
                fasta.write(">motif_" + str(self.unit_poss_list[i][j]) + "\n")
                fasta.write(unit_seq + "\n")
            fasta.close() 
        else: 
            pass
 
    def mafft(self, Exec=0):
        mafft_coms = []
        for period in self.peak_period_set:
            if os.getcwd().split("/")[-1] == self.output_dir:   
                fasta_name =  "unit_seq.fasta"
            else: 
                fasta_name =  "./" + self.output_dir + "/unit_seq.fasta"
            
            if self.dtype == "nucl":
                mafft_coms.append("mafft {} --quiet --auto {} > {}".format(self.option_mafft, fasta_name, fasta_name.replace("unit_seq.fasta","align.unit_seq.fasta"))) 
            else: 
                mafft_coms.append("mafft {} --quiet --amino --auto {} > {}".format(self.option_mafft, fasta_name, fasta_name.replace("unit_seq.fasta","align.unit_seq.fasta")))
        self.mafft_com = "|".join(mafft_coms) 
        if Exec == 1:
            subprocess.call(self.mafft_com,shell=True)
        
    def decide_query(self):
        self.query_list             = []
        self.se_sets_list           = [] 
        self.aligned_positions_list = []
        for i,period in enumerate(self.peak_period_set[0:1]):
            fasta_name = "align.unit_seq.fasta" 
            fasta = open(fasta_name)
            seqs  = read_seq_data(fasta)
            fasta.close() 
            #Bit score calcuation using WebLogo package
            if self.dtype == "nucl":  	
                seqs.alphabet = std_alphabets["dna"] 
                data          = LogoData.from_seqs(seqs)	
                data.alphabet = std_alphabets["dna"]
                options       = LogoOptions()
                format        = LogoFormat(data, options)
                try:
                    fout = open("weblogo.txt","wb")
                    fout.write(txt_formatter(data, format).decode("utf-8"))    
                except:
                    fout = open("weblogo.txt","wb")
                    fout.write(txt_formatter(data, format))    
                fout.close()
            else:
                seqs.alphabet = std_alphabets["protein"] 
                data          = LogoData.from_seqs(seqs)	
                data.alphabet = std_alphabets["protein"] 
                options       = LogoOptions()
                format        = LogoFormat(data, options)
                try:
                    fout = open("weblogo.txt","wb")
                    fout.write(txt_formatter(data, format).decode("utf-8"))    
                except:
                    fout = open("weblogo.txt","wb")
                    fout.write(txt_formatter(data, format))    
                fout.close()
        
            fasta         = SeqIO.to_dict(SeqIO.parse(fasta_name,"fasta"))
            logo_result   = [line.split("\t") for line in open("weblogo.txt") if line[0] != "#"]
            bitscore_list = [float(elements[-4])/0.69 for elements in logo_result] 
            weight_list   = [float(elements[-1]) for elements in logo_result] 
            frequent_seq  = ""
            bit_weight_list = [] 
                    
            #High consensus letters extractoin from the bit score array.
            j = 0
            for char_list in zip(*(list(map(lambda x: str(x.seq).upper(), fasta.values())))):
                count_list = []
                char_set   = list(set(char_list))
                if "-" in char_set:
                    char_set.remove("-") 
                for char in char_set:
                    count_list.append(char_list.count(char))    
                
                if  weight_list[j] >= self.tq:
                    freq = max(count_list) * 1.0 / sum(count_list)
                    frequent_seq += char_set[count_list.index(max(count_list))]
                    bit_weight_list.append(bitscore_list[j] * weight_list[j])
                j += 1 
            query_candidates = find_candidates(bit_weight_list ,thresh=1.0, gap=self.tr, buf=0, min_score=self.tb)

            #If high score regions are separated by both sides. concat them.
            if query_candidates[0][2] == 0 and query_candidates[-1][3] == len(frequent_seq) and len(query_candidates) > 1: 
                query_seq = frequent_seq[query_candidates[-1][2]:query_candidates[-1][3]] + frequent_seq[query_candidates[0][2]:query_candidates[0][3]] 
            elif len(query_candidates) > 1:
                query_seq = ""
                for candidate in query_candidates:
                    query_seq += frequent_seq[candidate[2]:candidate[3]]
            else:    
                query_seq = frequent_seq[query_candidates[0][2]:query_candidates[0][3]]
            query = open("query.fasta","w")
            query.write(">query" + "\n") 
            query.write(query_seq + "\n") 
            query.close()
            se_sets = [[pos + query_candidates[0][2], pos + query_candidates[0][2] + len(query_seq)] for pos in self.unit_poss_list[i]]
            aligned_positions = [[0,len(query_seq)] for pos in self.unit_poss_list[i]]
            self.aligned_positions_list.append(aligned_positions) 
            self.se_sets_list.append(se_sets)
            self.query_list.append(query_seq) 
    
    def make_se_sets(self):
        self.variable_query_list = []
        self.repeat_num_list     = []
        for i,period in enumerate(self.peak_period_set[0:1]):
            self.variable_query_list.append(self.query_list[i])
            self.repeat_num_list.append(len(self.se_sets_list[i]))
            if len(self.query_list[i]) > self.k_size:
                fasta_name  = "align.unit_seq.fasta" 
                blast_file  = open("./blast.txt")
                subject_seq = open("./subject.fasta").readlines()[1].rstrip() 
                blast = [line.rstrip().split("\t") for line in blast_file]
                blast.sort(key=lambda x:(int(x[11]),float(x[-2])))
                
                if len(blast) > 2:
                    #If ovelap region between blast hits was detected, the hit represent less E-value compared to others were selected. 
                    fasta = open(fasta_name,"w")
                    blast = [elements for elements in blast if float(elements[-2]) <= 0.1 or float(elements[6])/float(elements[5]) >= 0.5]
                    consensus_motifs = []
                    true_motifs = [] 
                    pre_end = 0 
                    
                    new_blast = [] 
                    for elements in blast:
                        consensus_motif = ""
                        true_motif      = ""
     
                        if int(elements[9])-1 > 0:
                            for j in range(int(elements[9])-1): 
                                consensus_motif += "-"
                                if (int(elements[11])-1)-(int(elements[9])-1) + j < 0:
                                    true_motif += "-"
                                else:
                                    true_motif += subject_seq[(int(elements[11])-1)-(int(elements[9])-1) + j]
                        
                        for char_combi in zip(elements[1],elements[3]):
                            if char_combi[0] == "-":
                                pass 
                            else: 
                                consensus_motif += char_combi[1]
                                true_motif      += char_combi[1] 
                                
                        if len(self.query_list[i]) - int(elements[10]) > 0: 
                            for j in range(int(elements[10]), len(self.query_list[i])):
                                consensus_motif += "-"
                                if int(elements[12])-1+j >= len(subject_seq):
                                    true_motif += "-" 
                                else:
                                    true_motif      += subject_seq[int(elements[12])-1+j]
                                                       
                        if int(elements[11]) >= pre_end: 
                            consensus_motifs.append(consensus_motif)     
                            true_motifs.append(true_motif) 
                            pre_end    = int(elements[12])
                            pre_evalue = float(elements[-2]) 
                            new_blast.append(elements) 
                        
                    #The following code search the start position of repeat motif.
                    for s, char_list in enumerate(list(zip(*consensus_motifs))):
                        count_list = []  
                        char_set   = list(set(char_list))
                        for char in char_set:
                            count_list.append(char_list.count(char))
                        if "-" not in char_set or count_list[char_set.index("-")] * 1.0 / sum(count_list) < 0.5:
                            break
                    
                    #The following code search the end position of repeat motif.
                    for e, char_list in enumerate(reversed(list(zip(*consensus_motifs)))):
                        count_list = []  
                        char_set   = list(set(char_list))
                        for char in char_set:
                            count_list.append(char_list.count(char))
                        if "-" not in char_set or count_list[char_set.index("-")] * 1.0 / sum(count_list) < 0.5:
                            break
                    
                    #Threshold of Interspace or Periodic
                    if abs(period - (len(consensus_motif) - s - e)) <= 5:
                        s, e = 0, 0 
                    
                    #Cutting motif region from query sequence for blast
                    se_sets           = []
                    aligned_positions = [] 
                    new_motifs        = [] 
                    for n, consensus_motif in enumerate(consensus_motifs):
                        ms = int(blast[n][11])-1-(int(blast[n][9])-1) + s
                        me = int(blast[n][12])+(len(self.query_list[i]) - int(blast[n][10])) - e
                        se_sets.append([ms, me])  
                        aligned_positions.append([(int(blast[n][9])-1)-s,me - ms - (len(self.query_list[i]) - int(blast[n][10])) + e])
                        new_motif       = true_motifs[n][s:len(consensus_motif)-e]  
                        consensus_motif = consensus_motifs[n][s:len(consensus_motif)-e]  
                        if (int(blast[n][9])-1) - s <= 5:
                            pass     
                        else:
                            new_motif = consensus_motif[:(int(blast[n][9])-1)-s] + new_motif[(int(blast[n][9])-1)-s:]
                        
                        if (len(self.query_list[i]) - int(blast[n][10])) - e <= 5:
                            pass
                        else:
                            new_motif = new_motif[:e-(len(self.query_list[i]) - int(blast[n][10]) - 1)] + consensus_motif[e-(len(self.query_list[i]) - int(blast[n][10]) - 1):] 
                        new_motifs.append(new_motif)

                    flag = 0 
                    true_motifs = [] 
                    combi = list(zip(se_sets,new_motifs)) 
                    combi.sort(key=lambda x:x[0][0]) 
                    se_sets, new_motifs = list(zip(*combi))
                    se_sets = [list(se_set) for se_set in se_sets] 
                    for n, se_set in enumerate(se_sets):
                        if flag == 1:
                            pass 
                        else:
                            fasta.write(">motif_"+ str(se_set[0]) + "\n")
                            fasta.write(new_motifs[n] + "\n")
                            true_motifs.append(new_motifs[n])
                            start = se_set[0] 
                        if n < len(se_sets)-1 and start + len(new_motifs[n]) - new_motifs[n].count("-") > se_sets[n+1][0]:
                            flag = 1
                        else:
                            flag = 0 

                    #Save repeat motif and search the varible position in the repeat motif.
                    true_query = ""
                    variable_query = "" 
                    for char_list in zip(*true_motifs):
                        count_list = []
                        char_set   = list(set(char_list))
                        for char in char_set:
                            count_list.append(char_list.count(char))
                        #skip gap position
                        if char_set[count_list.index(max(count_list))] != "_" and char_set[count_list.index(max(count_list))] != "-":
                            char = char_set[count_list.index(max(count_list))] 
                            true_query += char 
                            if max(count_list) * 1.0 / (sum(count_list)-char_list.count("_")-char_list.count("-")) < 0.6:
                                variable_query += "*"
                            else: 
                                variable_query += char
                    
                    fasta.close()
                    self.query_list[i]             = true_query
                    self.variable_query_list[i]    = variable_query
                    self.se_sets_list[i]           = se_sets        
                    self.aligned_positions_list[i] = aligned_positions
                    self.repeat_num_list[i]        = len(true_motifs)
                blast_file.close()
    
    def blast(self, Exec=0):
        blast_coms = []
        for i, period in enumerate(self.peak_period_set):
            if os.getcwd().split("/")[-1] == self.output_dir:   
                query_name   = "./query.fasta"
                subject_name = "./subject.fasta" 
                blast_name   = "./blast.txt"
            else: 
                query_name   = "./" +  self.output_dir + "/query.fasta"
                subject_name = "./" +  self.output_dir + "/subject.fasta"
                blast_name   = "./" +  self.output_dir + "/blast.txt"
            
            if self.dtype == "nucl":
                #options   = '-strand plus -task blastn-short -penalty -2 -outfmt "6 qseqid qseq sseqid sseq pident qlen length mismatch gapopen qstart qend sstart send gaps evalue bitscore"'
                blast_com = "blastn -query {} -subject {} {} -out {}".format(query_name, subject_name, self.option_blastn, blast_name) 
            else:
                #options   = '-task blastp-short -outfmt "6 qseqid qseq sseqid sseq pident qlen length mismatch gapopen qstart qend sstart send gaps evalue bitscore"'
                blast_com = "blastp -query {} -subject {} {} -out {}".format(query_name, subject_name, self.option_blastp, blast_name) 
        
            blast_coms.append(blast_com)
        self.blast_com = "|".join(blast_coms)
        
        if Exec == 1:
            subprocess.call(blast_com,shell=True)
        
    def make_motif_array(self):
        self.range_list         = []
        self.repeat_type_list   = [] 
        self.periodicity_list   = [] 
        new_aligned_positions_list = [] 
        new_peak_period_set        = [] 
        new_query_list             = []
        new_variable_query_list    = [] 
        new_se_sets_list           = [] 
        new_repeat_num_list        = [] 
        for i,period in enumerate(self.peak_period_set[:1]):
            self.se_sets_list[i].sort()
            p_min = period - max(int(round(0.2*period)),5) if period - max(int(round(0.2*period)),5) > 0 else 0
            p_max = period + max(int(round(0.2*period)),5) if period + max(int(round(0.2*period)),5) < 1000 else 1000
            matrix = self.period_matrix[p_min:p_max+1,:].astype(np.int)
            score_array   = self.region_score_array
            se_sets       = self.se_sets_list[i]
            intensity     = 0
            if period < 6:
                se_sets = [[se_set[0], se_set[1]] for se_set in se_sets if np.sum(matrix[:,se_sets[0][0]:se_sets[-1][1]]) * 1.0 > 0]
            else:
                flag = 0
                new_se_sets = [] 
                for j, se_set in enumerate(se_sets[0:-1]):
                    if (period - int(round(0.2*period)) < se_sets[j+1][0] - se_set[0] and (period + int(round(0.2*period)) > se_sets[j+1][0] - se_set[0])) or (period - int(round(0.2*period)) < se_sets[j+1][1] - se_set[1] and (period + int(round(0.2*period)) > se_sets[j+1][1] - se_set[1])):
                        new_se_sets.append(se_set) 
                        flag = 1         
                    else: 
                        flag = 0 
                if flag == 1:
                    new_se_sets.append(se_sets[-1]) 
                else:
                    pass 
                se_sets = new_se_sets

            k = 0
            if len(se_sets) > 2:
                start = se_sets[0][0] if se_sets[0][0] >= self.hra_range[0]-self.hra_range[2] else self.hra_range[0]-self.hra_range[2]
                end   = se_sets[-1][1] if se_sets[-1][1] <= self.hra_range[1]-self.hra_range[2] else self.hra_range[1]-self.hra_range[2]
                if start < 0:
                    start = 0 

                if period < 6:
                    intensity = np.sum(matrix[:,start:end], axis=1) * 1.0
                else:
                    intensity = np.sum(matrix[matrix.shape[0]//2-int(round(0.2*period)):matrix.shape[0]//2+int(round(0.2*period))+1,start:end], axis=1) * 1.0 
                
                if np.sum(score_array[start:end]) == 0:
                    periodicity = 0
                else:
                    periodicity = np.sum(intensity) / np.sum(score_array[start:end])
                
                if periodicity < 0.5 and len(self.peak_period_set) > i + 1:
                    for add_period in self.peak_period_set[i+1:]:
                        if add_period % period == 0 or period % add_period == 0:
                            p_min = add_period - max(int(round(0.2*add_period)),5) if period - max(int(round(0.2*add_period)),5) > 0 else 0
                            p_max = add_period + max(int(round(0.2*add_period)),5) if period + max(int(round(0.2*add_period)),5) < 1000 else 1000
                            add_matrix = self.period_matrix[p_min:p_max+1,:].astype(np.int)
                            periodicity += np.sum(add_matrix[add_matrix.shape[0]/2-int(round(0.2*add_period)):add_matrix.shape[0]/2+int(round(0.2*add_period))+1,start:end]) / np.sum(score_array[start:end])
            else:
                periodicity = 0
            
            threash = self.tp
            if (periodicity >= threash) and len(se_sets) > 2 and (np.argmax(intensity) == len(intensity)//2 or period < 6):
                
                if period - len(self.query_list[i]) >  5:
                    repeat_type = "inter_space"
                else:
                    repeat_type = "tandem"
                
                new_peak_period_set.append(period)
                new_query_list.append(self.query_list[i]) 
                new_variable_query_list.append(self.variable_query_list[i])
                new_se_sets_list.append(se_sets)  
                new_aligned_positions_list.append(self.aligned_positions_list[i]) 
                new_repeat_num_list.append(self.repeat_num_list[i])   
                if self.dtype == "nucl":
                    self.range_list.append([self.hra_range[2] + se_sets[0][0], self.hra_range[2] + se_sets[-1][1]])
                else:
                    self.range_list.append([se_sets[0][0], se_sets[-1][1]])

                self.periodicity_list.append(periodicity) 
                self.repeat_type_list.append(repeat_type)
        
        self.peak_period_set        = new_peak_period_set
        self.query_list             = new_query_list
        self.variable_query_list    = new_variable_query_list
        self.se_sets_list           = new_se_sets_list
        self.aligned_positions_list = new_aligned_positions_list 
        self.repeat_num_list        = new_repeat_num_list
    
    def make_feature(self, make_gb=True, make_tsv=True):
        new_feat_list = []
        for i, period in zip([0],self.peak_period_set):
            #Present version of SPADE don't care abount the position of blast hits.
            #If repeating motifs are separated by space whose length is times longer than periodd, the repeat would be evaluated as 
            #two different repeats at next version of spade.
            se_sets = [] 
            for line in open("align.unit_seq.fasta"):
                if line[0] == ">":
                    s = int(line.rstrip().split("_")[-1]) 
                else:
                    e = s + len(line.rstrip()) 
                    e = s + period if e-s > period else e
                    se_sets.append([s,e]) 
            
            if se_sets[0][0] < 0: 
                se_sets[0][0] = 0 
            
            if self.dtype == "nucl" or self.feature == None:
                #For genome
                part_list_all = [] 
                for se_set in se_sets: 
                    part_list_all.append([self.hra_range[2] + se_set[0], self.hra_range[2] + se_set[1]]) 
                    
                #For partial genbank
                part_list = [] 
                for se_set in se_sets:
                    part_list.append([se_set[0], se_set[1]])
   
            else:
                #Present code dosen't adapt the protein periodic repeat in coding region which have splicing.
                #I will collect the following code as soon as possible.

                #For genome 
                part_list_all = [] 
                for se_set in se_sets: 
                    if self.feature.strand == 1: 
                        part_list_all.append([self.feature.location.start + 3  * (self.hra_range[2] + se_set[0]), self.feature.location.start + 3 * (self.hra_range[2] + se_set[1])])

                    else: 
                        part_list_all.append([self.feature.location.end - 3 * (self.hra_range[2] + se_set[1]), self.feature.location.end - 3 * (self.hra_range[2] + se_set[0])]) 

                if self.feature.strand == -1:
                    part_list_all.reverse()

                #For partial genbank 
                part_list = [] 
                for se_set in se_sets: 
                    part_list.append([se_set[0], se_set[1]])

            locations = []  
            for part in part_list_all:
                locations.append(FeatureLocation(part[0], part[-1]))  
            location = CompoundLocation(locations) 
            new_feat = SeqFeature(location , type="repeat_region", qualifiers={"dir_name":["",],"note":["",],"hrr_range":["",],"sequence_type":["",],"period":[], "rpt_type":[], "rpt_num":[], "periodicity_score":[], "rpt_unit_seq":[],"unmasked_rpt_unit_seq":[]})                
            new_feat.qualifiers["dir_name"][0]      = self.output_dir
            new_feat.qualifiers["note"][0]          = "Auto annotated by SPADE"
            new_feat.qualifiers["hrr_range"][0]     = ",".join(list(map(str,self.hra_range)))
            new_feat.qualifiers["sequence_type"][0] = self.dtype
            new_feat.qualifiers["period"].append(str(period)) 
            new_feat.qualifiers["rpt_type"].append(self.repeat_type_list[i])
            new_feat.qualifiers["rpt_num"].append(str(self.repeat_num_list[i])) 
            new_feat.qualifiers["periodicity_score"].append(str(self.periodicity_list[i]))
            new_feat.qualifiers["rpt_unit_seq"].append(str(self.variable_query_list[i]))
            new_feat.qualifiers["unmasked_rpt_unit_seq"].append(str(self.query_list[i]))

            locations = []  
            for part in part_list:
                locations.append(FeatureLocation(part[0], part[-1])) 
            location = CompoundLocation(locations) 
            new_feat_part            = SeqFeature(location , type="repeat_region")
            new_feat_part.qualifiers = new_feat.qualifiers 
            self.record.features.append(new_feat_part)
            new_feat_list.append(new_feat)
        
        if self.dtype == "nucl":
            for num in range(0,len(self.record.features[:-1*(i+1)])):
                s, e = self.record.features[num].location.start, self.record.features[num].location.end 
                strand = self.record.features[num].strand 
                start, end = s-self.hra_range[2], e-self.hra_range[2]
                if start < 0:
                    start = 0
                self.record.features[num].location = FeatureLocation(start,end,strand=strand) 
        
        self.record.features.sort(key=lambda x: x.location.start)
        if make_gb == True and len(new_feat_list) > 0: 
            record_handle  = open("repeat.gbk","w")
            self.record.id = self.record.id[0:16]
            SeqIO.write(self.record, record_handle, "genbank")  
        os.system("rm blast.txt subject.fasta query.fasta")
        return new_feat_list

    def extract(self, gb, start_pos, end_pos,  make_file=False):
        partial_gb = SeqRecord(Seq(str(gb.seq[start_pos:end_pos]),Alphabet.DNAAlphabet())) 
        partial_gb.id = gb.id
        partial_gb.description = gb.description
        feat_list = [] 
        for afeat in gb.features:
            flag  = 0
            for part in afeat.location.parts:
                start, end = part.start, part.end
                if (start <= end_pos and end >= start_pos and afeat.type != "source"):
                    flag = 1
            
            if flag == 1:
                partial_gb.features.append(copy.deepcopy(afeat))
                feat_list.append(copy.deepcopy(afeat))
            elif afeat.location.start > end_pos:
                break
        
        return partial_gb, feat_list
     
    def make_motif_logo(self):
        vs.motif_logo(self.dtype) 
        
    def make_figure(self):
        vs.load_data(self.dtype, 1, self.k_size, 10) 
    
    def save_region_record(self):
        pass 

if __name__ == "__main__":
    p = argparse.ArgumentParser(add_help=False)
    p.add_argument("-V", "--version", action="store_true", default=False)
    p.add_argument('-h', "--help", action='store_true', default=False) 

    p.add_argument("-in","--input", type=str, default="None", help="Input file name")
    p.add_argument("-f", type=str, default="auto", choices=("genbank","fasta","auto"))
    p.add_argument("-t", type=str, default="auto", choices=("auto","nucl","protein"))

    p.add_argument("-Nk", type=int, default=10)
    p.add_argument("-Nw", type=int, default=1000)
    p.add_argument("-Ns", type=int, default=20)  
    p.add_argument("-Nm", type=int, default=1000)
    p.add_argument("-Ng", type=int, default=200)
    p.add_argument("-Np", type=float, default=0.5)
    p.add_argument("-Nq", type=float, default=0.5)
    p.add_argument("-Nu", type=float, default=0.8)
    p.add_argument("-Nr", type=int, default=5)
    
    p.add_argument("-Pk", type=int, default=3)
    p.add_argument("-Pw", type=int, default=300)
    p.add_argument("-Ps", type=int, default=6)  
    p.add_argument("-Pm", type=int, default=300)
    p.add_argument("-Pg", type=int, default=50)
    p.add_argument("-Pp", type=float, default=0.3)
    p.add_argument("-Pq", type=float, default=0.5)
    p.add_argument("-Pu", type=float, default=0.8)
    p.add_argument("-Pr", type=int, default=5)
    
    p.add_argument("--mafft", type=str, default="--auto") 
    p.add_argument("--blastn", type=str, default='-strand plus -task blastn-short -penalty -2 -outfmt "6 qseqid qseq sseqid sseq pident qlen length mismatch gapopen qstart qend sstart send gaps evalue bitscore"') 
    p.add_argument("--blastp", type=str, default='-task blastp-short -outfmt "6 qseqid qseq sseqid sseq pident qlen length mismatch gapopen qstart qend sstart send gaps evalue bitscore"')

    p.add_argument("-n", "--num_threads", type=int, default=1) 
    p.add_argument("-v", type=str, default="Y", choices=("Y","N"))
    p.add_argument("-d", "--delete", action="store_true", default=False)
    args = p.parse_args()
    
    if args.help:
        import spade_help as sh
        print(sh.help_description)
        exit() 

    if args.version:
        print("SPADE: {}".format(__version__)) 
        exit()
    
    import numpy as np 
    import multiprocessing as mp
    import visualisation as vs
    from scipy import signal
    from Bio import Alphabet
    from Bio.Seq import Seq
    from Bio import SeqIO
    from Bio.SeqRecord import SeqRecord
    from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation
    from kmer_count import *
    from weblogolib import *
    spade = SPADE() 
    if args.input != "None":
        spade.load(args.input) 
    
    if args.delete:
        args.delete = 1 
 
    spade.run(args.num_threads,[args.Nk, args.Nw, args.Ng, args.Ns, args.Np, args.Nu, args.Nm, args.Nr, args.Nq, args.Pk, args.Pw, args.Pg, args.Ps, args.Pp, args.Pu, args.Pm, args.Pr, args.Pq, args.v, args.delete, args.mafft, args.blastn, args.blastp])
    finish = open("finish.txt","w")
    finish.write("fnish")
    finish.close() 
