# -*- coding: utf-8 -*-

import numpy as np
import igraph as ig
from .utils import fisher_test

def rwr_values(network, graph_nodes, ini_pos, ini_neg, seeds, seeds_pos, 
               seeds_neg, network_path, network_random_path, damping, res_folder):
    
    number_of_nodes=network.vcount()
    empirical_rwr=np.zeros((6,number_of_nodes),dtype=float)

    empirical_values={}
    pvalues={}
    for i in network.vs:
        pvalues[i["name"]]=np.zeros(6,dtype=int)
        empirical_values[i["name"]]=np.zeros(6,dtype=float)

    pos=False
    if ini_pos[0]!="False":
        pos=True
        to_delete=[]
        for j in network.vs.select(name_in=ini_pos):
            to_delete.append(j.index)

    neg=False
    if ini_neg[0]!="False":
        neg=True
        to_delete=[]
        for j in network.vs.select(name_in=ini_neg):
            to_delete.append(j.index)
    flag_pos=0
    flag_neg=0
    for i in enumerate(seeds):
        if len(i[1])>0:
            if pos==True and flag_pos==0 and i[0]<=2:
                network.delete_vertices(to_delete)
                flag_pos=1
            if neg==True and flag_neg==0 and i[0]>2:
                network.delete_vertices(to_delete)
                flag_neg=1

            number_of_nodes=network.vcount()
            reset_vertex=np.zeros(number_of_nodes)

            for j in network.vs.select(name_in=i[1]):
                if i[0]<=2:
                    reset_vertex[j.index]=seeds_pos[j["name"]]
                else:
                    reset_vertex[j.index]=seeds_neg[j["name"]]
            pagerank=np.array(network.personalized_pagerank(reset=reset_vertex,directed=False, damping=damping, weights='weight',implementation='prpack'))
            if flag_pos==1 and i[0]<=2:
                for j in to_delete:
                    pagerank=np.concatenate((pagerank[:j], [0], pagerank[j:]))

            if flag_neg==1 and i[0]>2:
                for j in to_delete:
                    pagerank=np.concatenate((pagerank[:j], [0], pagerank[j:]))

            empirical_rwr[i[0]]=pagerank

        if i[0]==2 and pos==True:
            network = ig.Graph.Read_Ncol(network_path, weights=True, directed=False)
            pos=False
            flag_pos=0
            to_delete=[]


    for i in enumerate(empirical_rwr.T):
        empirical_values[graph_nodes[i[0]]]=i[1]

    for ii in range(1000):
        network_random = ig.Graph.Read_Ncol(network_random_path+str(ii)+".txt", weights=True, directed=False)
        nodes=network_random.vs["name"]
        number_of_nodes=network_random.vcount()
        random_rwr=np.zeros((6,number_of_nodes),dtype=float)
        pos=False
        if ini_pos[0]!="False":
            to_delete=[]
            for j in network_random.vs.select(name_in=ini_pos):

                to_delete.append(j.index)
            pos=True
        neg=False
        if ini_neg[0]!="False":
            to_delete=[]
            for j in network_random.vs.select(name_in=ini_neg):
                to_delete.append(j.index)
            neg=True
        flag_pos=0
        flag_neg=0

        for jj in enumerate(seeds):
            if len(jj[1])>0:
                if pos==True and flag_pos==0 and jj[0]<=2:
                    network_random.delete_vertices(to_delete)
                    flag_pos=1
                if neg==True and flag_neg==0 and jj[0]>2:
                    network_random.delete_vertices(to_delete)
                    flag_neg=1

                number_of_nodes=network_random.vcount()
                reset_vertex=np.zeros(number_of_nodes)

                for j in network_random.vs.select(name_in=jj[1]):
                    if jj[0]<=2:
                        reset_vertex[j.index]=seeds_pos[j["name"]]
                    else:
                        reset_vertex[j.index]=seeds_neg[j["name"]]
                prandom=np.array(network_random.personalized_pagerank(reset=reset_vertex,directed=False, damping=damping, weights='weight',implementation='prpack'))

                if flag_pos==1 and jj[0]<=2:
                    for j in to_delete:
                        prandom=np.concatenate((prandom[:j], [0], prandom[j:]))

                if flag_neg==1 and jj[0]>2:
                    for j in to_delete:
                        prandom=np.concatenate((prandom[:j], [0], prandom[j:]))

                if jj[0]==2 and pos==True:
                    network_random = ig.Graph.Read_Ncol(network_random_path+str(ii)+".txt", weights=True, directed=False)
                    flag_pos=0
                    pos=False
                    to_delete=[]

                random_rwr[jj[0]]=prandom
            else:
                if jj[0]==2:
                    network_random = ig.Graph.Read_Ncol(network_random_path+str(ii)+".txt", weights=True, directed=False)
                    flag_pos=0
                    pos=False
                    to_delete=[]
        random_rwr=random_rwr.T
        for i in enumerate(nodes):
            pvalues[i[1]]=np.add(pvalues[i[1]],np.greater(empirical_values[i[1]],random_rwr[i[0]]))

    '''
    Output.
     '''
    #writing the results for global propagation
    f1=open(res_folder+"rwr_scores.txt","w")
    f1.write('uniprotid\tincreased_tyr\tincreased_kinases\tincreased_substrates\tdecreased_tyr\tdecreased_kinases\tdecreased_substrates\n')
    for i in empirical_values:
        f1.write(i+"\t"+"\t".join(map(str,empirical_values[i]))+"\n")
    f1.close()
    
    f1=open(res_folder+"pvalues.txt","w")
    f1.write('uniprotid\tincreased_tyr\tincreased_kinases\tincreased_substrates\tdecreased_tyr\tdecreased_kinases\tdecreased_substrates\n')
    for i in pvalues:
        f1.write(i+"\t"+"\t".join(map(str,pvalues[i]))+"\n")
    f1.close()
    
    f1=open(res_folder+"start_seeds.txt","w")
    f1.write('uniprotid\tLFC_value\n')
    for i in seeds[0]+seeds[1]+seeds[2]:
        f1.write(i+"\t"+str(seeds_pos[i])+"\n")

    for i in seeds[3]+seeds[4]+seeds[5]:
        f1.write(i+"\t"+str(-seeds_neg[i])+"\n")
    f1.close()
 
    # # Returned variables.
    # return empirical_values,pvalues


# Obtain the significant rwr nodes according to pvalues.txt, and split into increased/decreased nodes.
def pvalue_split(res_folder, seeds, graph_nodes, rwr_threshold,
                 fisher_threshold, fisher_geneset, uniprot_to_gene, geneset_path):
    '''
    Separate the pvalues for upregulated/downregulated nodes.
    '''
    pvalues_pos=[]
    pvalues_neg=[]
 
    f1=open(res_folder+"pvalues.txt")
    seq=f1.readline()
    seq=f1.readline()
    #the pvalues need to be a parameter for the users  
    while(seq!=""):
        seq=seq.strip().split("\t")
        seq[1:]=np.array(seq[1:],float)
        if max(seq[1:4]) > 1000 - rwr_threshold*1000:
            pvalues_pos.append(seq[0])
        if max(seq[4:]) > 1000 - rwr_threshold*1000:
            pvalues_neg.append(seq[0])
        seq=f1.readline()

    pvalues_pos=list(set(pvalues_pos+seeds[0]+seeds[1]+seeds[2]))
    pvalues_neg=list(set(pvalues_neg+seeds[3]+seeds[4]+seeds[5]))
    #perfrom fishertest on the seeds and rwr nodes
    
    fname = "increased_rwr_fisher"
    fisher_test(protein_list=pvalues_pos,
                starting_proteins=list(set(seeds[0]+seeds[1]+seeds[2])),
                fname=fname,
                threshold=fisher_threshold,
                component=fisher_geneset,
                path_def=res_folder,
                uniprot_to_gene=uniprot_to_gene,
                geneset_path=geneset_path,
                )
    fname = "decreased_rwr_fisher"
    fisher_test(protein_list=pvalues_neg,
            starting_proteins=list(set(seeds[3]+seeds[4]+seeds[5])),
            fname=fname,
            threshold=fisher_threshold,
            component=fisher_geneset,
            path_def=res_folder,
            uniprot_to_gene=uniprot_to_gene,
            geneset_path=geneset_path,
            )
    fname = "increased_seed_fisher"
    fisher_test(protein_list=list(set(seeds[0]+seeds[1]+seeds[2])),
            starting_proteins=list(set(seeds[0]+seeds[1]+seeds[2])),
            fname=fname,
            threshold=fisher_threshold,
            component=fisher_geneset,
            path_def=res_folder,
            uniprot_to_gene=uniprot_to_gene,
            geneset_path=geneset_path,
        )
    fname = "decreased_seed_fisher"
    fisher_test(protein_list=list(set(seeds[3]+seeds[4]+seeds[5])),
            starting_proteins=list(set(seeds[3]+seeds[4]+seeds[5])),
            fname=fname,
            threshold=fisher_threshold,
            component=fisher_geneset,
            path_def=res_folder,
            uniprot_to_gene=uniprot_to_gene,
            geneset_path=geneset_path,
        )

    
    '''
    Separate the rwr_scores for upregulated/downregulated nodes.
    '''
    f1=open(res_folder+"rwr_scores.txt")
    rwr_pos={}
    rwr_neg={}
    seq=f1.readline()
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
