# -*- coding: utf-8 -*-

import numpy as np

def pvalue_split(res_folder, seeds, graph_nodes):
    '''
    Separate the pvalues for upregulated/downregulated nodes.
    '''
    pvalues_pos=[]
    pvalues_neg=[]
 
    f1=open(res_folder+"pvalues.txt")
    seq=f1.readline()
    #the pvalues need to be a parameter for the users  
    while(seq!=""):
        seq=seq.strip().split("\t")
        seq[1:]=np.array(seq[1:],float)
        if max(seq[1:4])>950:
            pvalues_pos.append(seq[0])
        if max(seq[4:])>950:
            pvalues_neg.append(seq[0])
        seq=f1.readline()

    pvalues_pos=pvalues_pos+seeds[0]+seeds[1]+seeds[2]
    pvalues_neg=pvalues_neg+seeds[3]+seeds[4]+seeds[5]
    #perfrom fishertest on the seeds and rwr nodes
    

    
    '''
    Separate the rwr_scores for upregulated/downregulated nodes.
    '''
    f1=open(res_folder+"rwr_scores.txt")
    rwr_pos={}
    rwr_neg={}
    seq=f1.readline()
    while(seq!=""):
        seq=seq.strip().split("\t")
        rwr_pos[seq[0]]=np.mean(np.array(seq[1:4],float))
        rwr_neg[seq[0]]=np.mean(np.array(seq[4:],float))
        seq=f1.readline()

    '''
    Force node assignment to be unique (either up- or down-regulated).
    '''
    common=list(set(pvalues_pos).intersection(pvalues_neg))
    remove_from_pos=[]
    remove_from_neg=[]
    for i in common:
        if rwr_pos[i]<rwr_neg[i]:
            remove_from_pos.append(i)
        if rwr_neg[i]<rwr_pos[i]:
            remove_from_neg.append(i)

    pvalues_pos=list(set(pvalues_pos).difference(remove_from_pos))
    pvalues_neg=list(set(pvalues_neg).difference(remove_from_neg))

    pvalues_pos=list(set(pvalues_pos).intersection(graph_nodes))
    pvalues_neg=list(set(pvalues_neg).intersection(graph_nodes))
    
    return(pvalues_pos, pvalues_neg)
    
 
 
    
