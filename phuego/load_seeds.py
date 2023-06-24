# -*- coding: utf-8 -*-

import numpy as np

def kinase_classification(pfam_domain_path):
    tyr_kinase=[]
    st_kinase=[]
    f1=open(pfam_domain_path)
    seq=f1.readline()
    seq=seq.strip().split("\t")
    tyr_kinase.extend(seq[1:])
    seq=f1.readline()
    seq=seq.strip().split("\t")
    st_kinase.extend(seq[1:])
    return tyr_kinase,st_kinase
    
def load_seeds(pfam_domain_path, sim_mean_std_path, sim_all_folder_path, 
               test_path, graph_nodes):
    '''
    Loading seed nodes.
    '''
    # zscores_global.
    f1 = open(sim_mean_std_path)
    seq=f1.readline()
    seq=f1.readline()
    zscores_global={}
    while(seq!=""):
        seq=seq.strip().split("\t")
        zscores_global[seq[0]]=np.array(seq[1].split("|"),dtype=float)
        seq=f1.readline()

    # seed nodes.
    f1 = open(test_path) 
    seq=f1.readline()
    ssim={}
    seeds_pos={}
    seeds_neg={}
    while(seq!=""):
        seq= seq.strip().split("\t")
        seq[0]=float(seq[0])
        if seq[0]>0.0:
            seeds_pos[seq[1]]=seq[0]
        if seq[0]<-0.0:
            seeds_neg[seq[1]]=-seq[0]
        seq=f1.readline()
    
    '''
    Classification and semantic similarity.
    '''
    # Seed node classification.
    tyr_kinase,st_kinase=kinase_classification(pfam_domain_path)
    tyr_pos=list((set(tyr_kinase).intersection(seeds_pos.keys())).intersection(set(graph_nodes)))
    tyr_neg=list((set(tyr_kinase).intersection(seeds_neg.keys())).intersection(set(graph_nodes)))
    st_pos= list((set(st_kinase).intersection(seeds_pos.keys())).intersection(set(graph_nodes)))
    st_neg=list((set(st_kinase).intersection(seeds_neg.keys())).intersection(set(graph_nodes)))
    all_the_rest_pos=list((set(seeds_pos.keys()).difference(set(tyr_pos+st_pos))).intersection(set(graph_nodes)))
    all_the_rest_neg=list((set(seeds_neg.keys()).difference(set(tyr_neg+st_neg))).intersection(set(graph_nodes)))
    seeds=[tyr_pos,st_pos,all_the_rest_pos,tyr_neg,st_neg,all_the_rest_neg]

    # Semantic similarity for each seed node.
    ssim={}
    for i in list(set(tyr_pos+tyr_neg+st_pos+st_neg+all_the_rest_pos+all_the_rest_neg)):
        ssim[i]={}
        fsim=open(sim_all_folder_path+'/'+i+"_all.txt")
        # The first readline trim away the header.
        seq_sim=fsim.readline()
        seq_sim=fsim.readline()
        while(seq_sim!=""):
            seq_sim=seq_sim.strip().split("\t")
            ssim[i][seq_sim[1]]=float(seq_sim[2])
            seq_sim=fsim.readline()
            
    '''
    Output.
    '''
    return seeds_pos,seeds_neg,seeds,zscores_global,ssim