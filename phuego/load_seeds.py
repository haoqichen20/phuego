# -*- coding: utf-8 -*-

import numpy as np

def layer_classification(layer_path):
    layer1=[]
    layer2=[]
    f1=open(layer_path)

    # Read first layer definition. 
    seq=f1.readline()
    # If not defined, return two empty layers.
    if not seq:
        return "",[],"",[]
    seq=seq.strip().split("\t")
    layer1_name = seq[0]
    layer1.extend(seq[1:])

    # Read second layer definition. 
    seq=f1.readline()
    # If not defined, return layer1 and empty layer2.
    if not seq:
        return layer1_name, layer1, "", []
    seq=seq.strip().split("\t")
    layer2_name = seq[0]
    layer2.extend(seq[1:])

    return layer1_name, layer1, layer2_name, layer2

def load_seeds(layer_path, sim_mean_std_path, sim_all_folder_path, 
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
        seq[1]=float(seq[1])
        if seq[1]>0.0:
            seeds_pos[seq[0]]=seq[1]
        if seq[1]<-0.0:
            seeds_neg[seq[0]]=-seq[1]
        seq=f1.readline()

    '''
    Layer classification and semantic similarity.
    '''
    # Classify seed nodes into layers.
    layer1_name,layer1,layer2_name,layer2=layer_classification(layer_path)
    if not layer1_name:
        all_the_rest_pos=list(set(seeds_pos.keys()).intersection(set(graph_nodes)))
        all_the_rest_neg=list(set(seeds_neg.keys()).intersection(set(graph_nodes)))
        seeds=[[],[],all_the_rest_pos,[],[],all_the_rest_neg]
        print("Input is retained as 1 layer.")
    elif not layer2_name:
        layer1_pos=list((set(layer1).intersection(seeds_pos.keys())).intersection(set(graph_nodes)))
        layer1_neg=list((set(layer1).intersection(seeds_neg.keys())).intersection(set(graph_nodes)))
        all_the_rest_pos=list((set(seeds_pos.keys()).difference(set(layer1_pos))).intersection(set(graph_nodes)))
        all_the_rest_neg=list((set(seeds_neg.keys()).difference(set(layer1_neg))).intersection(set(graph_nodes)))
        seeds=[layer1_pos,[],all_the_rest_pos,layer1_neg,[],all_the_rest_neg]
        print(f"Input is divided into 2 layers: {layer1_name} and all the rest.")
    else:
        layer1_pos=list((set(layer1).intersection(seeds_pos.keys())).intersection(set(graph_nodes)))
        layer1_neg=list((set(layer1).intersection(seeds_neg.keys())).intersection(set(graph_nodes)))
        layer2_pos=list((set(layer2).intersection(seeds_pos.keys())).intersection(set(graph_nodes)))
        layer2_neg=list((set(layer2).intersection(seeds_neg.keys())).intersection(set(graph_nodes)))
        all_the_rest_pos=list((set(seeds_pos.keys()).difference(set(layer1_pos+layer2_pos))).intersection(set(graph_nodes)))
        all_the_rest_neg=list((set(seeds_neg.keys()).difference(set(layer1_neg+layer2_neg))).intersection(set(graph_nodes)))
        seeds=[layer1_pos,layer2_pos,all_the_rest_pos,layer1_neg,layer2_neg,all_the_rest_neg]
        print(f"Input is divided into 3 layers: {layer1_name}, {layer2_name}, and all the rest.")

    # Semantic similarity for each seed node.
    ssim={}
    for i in list(set([node for _layer in seeds for node in _layer])):
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
    return seeds_pos,seeds_neg,seeds,layer1_name,layer2_name,zscores_global,ssim