#!/usr/bin/env python
# -*- coding: utf-8 -*-
import numpy as np 
import collections
from Bio import SeqIO 
import sys

def make_kmer_set(seq,window_size,k,pos):
    kmer_dict = collections.defaultdict(lambda: [0,[]])
    for i in range(window_size):
        target = seq[i:i+k]
        kmer_dict[target][0] += 1
        kmer_dict[target][1].append(pos + i) 
    return kmer_dict

def update(seq,kmer_dict,pre_kmer,new_kmer,k,window_size,pos):
    kmer_dict[pre_kmer][0] -= 1 
    if kmer_dict[pre_kmer][0] < 1:
        del kmer_dict[pre_kmer]
    kmer_dict[new_kmer][0] += 1
    kmer_dict[new_kmer][1].append(pos)


def update_pp(seq,kmer_dict,pre_kmer,new_kmer,k,window_size,pos):
    kmer_dict[pre_kmer][0] -= 1 
    kmer_dict[new_kmer][0] += 1
    kmer_dict[new_kmer][1].append(pos)

def kmer_count(seq, mer_size, window_size, thresh=5.0, gap=200, buf=1000, min_score=0, seqtype="nucl"):
    kmer_dict = make_kmer_set(seq[1:],window_size,mer_size,0) 
    score_array = np.zeros(len(seq)) 
    query = seq[0:mer_size]
    start = -1
    end   = 0
    space = 0 
    flag  = 0 
    candidate_list = []
    for i in range(len(seq)-mer_size): 
        if query not in kmer_dict or ((query.count("N") == len(query) or query.count("n") == len(query)) and seqtype == "nucl"):
            pass
        else: 
            n = 0
            for j in kmer_dict[query][1]:
                if j < i + 1:
                    n += 1
                else:
                    score_array[j:j+mer_size] += 1
            
            kmer_dict[query][0] = len(kmer_dict[query][1][n:])
            kmer_dict[query][1] = kmer_dict[query][1][n:]
            score_array[i:i+mer_size] += kmer_dict[query][0]

        score = score_array[i]
        if score > min_score: 
            if score >= thresh: 
                flag, space, end  = 1, 0, 0 
                if start == -1:
                    start = i
                    end = 0

            else:
                if start == -1:
                    start = i
                    end = 0
        else:
            space += 1 
            if start > 0 and end == 0: 
                if flag == 0:
                    start = -1
                else:
                    end   = i 
            

        if (space > gap or i == len(seq) - mer_size-1) and flag > 0:
            if end == 0:
                end = i + mer_size + 1 
            if start - buf <  0:
                s_buf = 0
                e_buf = end + buf
            else:
                s_buf = start - buf
                e_buf = end + buf 
            
            candidate_list.append([start,end,s_buf,min((e_buf,len(seq)))])
            flag, space = 0, 0   
            start, end  = -1, 0

        query = seq[i+1:i+1+mer_size]
        update(seq,kmer_dict,query,seq[i+1+window_size:i+1+mer_size+window_size],mer_size,window_size,i+1+window_size)
    return score_array, candidate_list


def kmer_count_matrix(seq, mer_size, window_size, max_window_size=5000, seqtype="nucl"):
    kmer_dict = make_kmer_set(seq,window_size,mer_size,0) 
    if window_size > max_window_size:
        limit_window_size = max_window_size
    else:
        limit_window_size = window_size

    try:
        score_matrix = np.zeros((limit_window_size,len(seq)))
    except:
        score_matrix = np.zeros((limit_window_size,len(seq)), dtype=np.int16)

    score_array  = np.zeros(len(seq)) 
    query = seq[0:mer_size]
    for i in range(len(seq)-mer_size): 
        if query not in kmer_dict or ((query.count("N") == len(query) or query.count("n") == len(query)) and seqtype == "nucl"):
            pass
        else: 
            n = 0
            m = 0 
            pre_j = i
            inter_space_sum = 0
            for j in kmer_dict[query][1]:
                if j < i + 1:
                    n += 1
                else:
                    if m == 0:
                        origin_space = j - pre_j
                    inter_space = j - pre_j
                    if inter_space < window_size:
                        score_array[i:i+mer_size] += 1
                        score_array[j:j+mer_size] += 1
                    if inter_space < limit_window_size and origin_space < limit_window_size:
                        score_matrix[inter_space,j:j+mer_size] += 1
                        score_matrix[origin_space,i:i+mer_size] += 1
                        #score_matrix[origin_space,j:j+mer_size] += 1
                        #score_matrix[inter_space,i:i+mer_size] += 1
                    m += 1
                    pre_j = j

        query = seq[i+1:i+1+mer_size]
        update_pp(seq,kmer_dict,query,seq[i+1+window_size:i+1+mer_size+window_size],mer_size,window_size,i+1+window_size)
    return score_matrix, score_array, kmer_dict

def find_candidates(score_list, thresh = 5.0, gap = 200, buf = 1000, min_score=0, end_cut=False):
    """
    The function return candidate repeats from kmer-counting result. 
    The "gap" means max range less than thresh hold value." The default size is 200.
    The "buf" means the side length of repeat region. The default size is 1000.
    """
    start = -1
    end = 0
    space = 0 
    flag = 0 
    candidate_list = []
    for i, score in enumerate(score_list):
        if score > min_score and start == -1:
            start = i
            end = 0
        
        if start > 0 and end == 0 and score <= min_score: 
            if flag == 0:
                start = -1
            else:
                end   = i 
        
        if score >= thresh: 
            flag, space, end  = 1, 0, 0 
        elif score < thresh and score > min_score:
            pass
        else:
            space += 1 

        if (space > gap or i == len(score_list) - 1) and flag > 0:
            if end == 0:
                end = i + 1 
            
            if end_cut == False and i == len(score_list) - 1: 
                end = i + 1 
                
            if start - buf <  0:
                s_buf = 0
                e_buf = end + buf
            else:
                s_buf = start - buf
                e_buf = end + buf 
            
            candidate_list.append([start,end,s_buf,min((e_buf,len(score_list)))])
            flag, space = 0, 0   
            start, end  = -1, 0
    
    return candidate_list

if __name__ == "__main__":
    record_parse = SeqIO.parse(sys.argv[1],sys.argv[2])
    for record in record_parse:
        score_array, HRRs = kmer_count(str(record.seq),10,1000,seqtype="nucl") 
        score_array = score_array.tolist()
        ouput_data  = zip(range(len(score_array)),score_array)
        np.savetxt("{}_scores".format(record.id) + ".txt", np.array(ouput_data), delimiter="\t", fmt="%.0f")
        with open("{}_HRRs".format(record.id) + ".txt","w") as o:
            for hrr in HRRs:
                o.write(",".join(list(map(str,hrr))) + "\n") 

        
