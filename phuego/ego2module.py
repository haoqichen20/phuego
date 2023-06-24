# -*- coding: utf-8 -*-

from .utils import denoise_square
from .fishertest import fisher_test
import igraph as ig
import numpy as np
import os
import leidenalg as la
from itertools import combinations
from scipy.spatial.distance import jensenshannon


# Import network from phuego(), note that it's gic_raw, not gic
# Restructure output folders, removed the "test" variable.
# Delete folder_res, will use the same res_folder for write_modules()
# Damping in subnet.personalized_pagerank not touched.
# Import all necessary packages.
# Delete the commented code block.
# Moved supernodes_net.write() to last of function, 



def test_function(network, network_nodes, kde_cutoff, res_folder, uniprot_to_gene,
                  supernodes,all_nodes,direction,geneset_path,fisher_geneset, 
                  fisher_threshold,):
    
    for kde in kde_cutoff:        
        kde_net=network.induced_subgraph(network.vs.select(name_in=all_nodes[kde]),implementation='copy_and_delete')
        kde_nodes=kde_net.vs["name"]
        supernodes_net=ig.Graph()
        temp=[]
        for i in combinations(supernodes[kde].keys(),2):
            subnet=kde_net.induced_subgraph(kde_net.vs.select(name_in=supernodes[kde][i[0]]+supernodes[kde][i[1]]),implementation='copy_and_delete')
            subnet.vs.select(_degree=0).delete()
            subnet=denoise_square(subnet)
            nodes=subnet.vs["name"]
            n_nodes=len(nodes)
            index_net={k: v for v, k in enumerate(nodes)}
            components=subnet.connected_components()

            if i[0] in index_net and i[1] in index_net and components.membership[index_net[i[0]]]==components.membership[index_net[i[1]]]:
                reset_vertex=np.zeros(n_nodes)
                reset_vertex[index_net[i[0]]]=1.0
                pagerankA=np.array(subnet.personalized_pagerank(reset=reset_vertex,directed=False, damping=0.85, weights='weight',implementation='prpack'))

                reset_vertex=np.zeros(n_nodes)
                reset_vertex[index_net[i[1]]]=1.0
                pagerankB=np.array(subnet.personalized_pagerank(reset=reset_vertex,directed=False, damping=0.85, weights='weight',implementation='prpack'))
                js_dist=jensenshannon(pagerankA,pagerankB)

                if js_dist>0:
                    temp.append((i[0],i[1],js_dist))
        supernodes_net = ig.Graph.TupleList(temp,weights="weight",directed = False)
        #supernodes_net=denoise_square(supernodes_net)
        components = supernodes_net.connected_components(mode="strong")
        modules=[]
        mod={}
        for i in enumerate(components):
            if len(i[1])>=4:
                cc_net=supernodes_net.induced_subgraph(i[1],implementation='auto')
                nodes=cc_net.vs["name"]
                cluster = la.find_partition(cc_net, la.ModularityVertexPartition,weights='weight',n_iterations=-1,seed=42)
                for jj in cluster:
                    modules.append(cc_net.vs.select(jj)["name"])

            elif 2<=len(i[1])<4:
                sm_net=supernodes_net.induced_subgraph(i[1],implementation='auto')
                nodes=sm_net.vs["name"]
                modules.append(nodes)
        #print (modules)
        
        supernodes_net.write(f=res_folder+"/supernodes_net.txt",format="ncol")
        write_modules(clustering=modules,
                      nodes=supernodes[kde],
                      kde=kde,
                      direction=direction,
                      uniprot_to_gene=uniprot_to_gene,
                      res_folder=res_folder,
                      geneset_path=geneset_path,
                      fisher_geneset=fisher_geneset,
                      fisher_threshold=fisher_threshold,
                     )


# Remove the foldering structure here and replaced with names.
# Make fisher_test() input compatible with new fisher_test() function.
# Note that since this function is called within the for loop for kde_cutoff list, 
# the input is now a singular kde value, instead of kde_cutoff.

def write_modules(clustering,nodes, kde, direction, uniprot_to_gene, res_folder, 
                  geneset_path, fisher_geneset, fisher_threshold,):
    # base_folder=folder
    
    for j in enumerate(clustering):
        # if not os.path.exists(base_folder+"/module_"+str(j[0])):
        #     os.makedirs(base_folder+"/module_"+str(j[0]))
        # f1=open(base_folder+"/module_"+str(j[0])+"/cluster.txt","w")
        
        f1=open(res_folder+direction+"_module_"+str(j[0])+"_cluster_"+str(kde)+".txt", "w")
        fisher_proteins=set()
        for k in j[1]:
            f1.write(k+"\t"+"\t".join(nodes[k])+"\n")
            fisher_proteins.update(nodes[k])
        
        # fisher_folder=base_folder+"/module_"+str(j[0])+"/fisher/"
        # if not os.path.exists(fisher_folder):
        #     os.makedirs(fisher_folder)
        fname = direction+"_module_"+str(j[0])+"_fisher_"+str(kde)+"_"
        fisher_test(protein_list=fisher_proteins,
                    starting_proteins=j[1],
                    threshold=fisher_threshold,
                    component=fisher_geneset,
                    path_def=res_folder,
                    uniprot_to_gene=uniprot_to_gene,
                    geneset_path=geneset_path,
                    fname=fname,
                    )
            # fisher_proteins,0.05, ["P","C","F","R","K","RT","B","D"],fisher_folder,j[1],uniprot_to_gene,"intact")