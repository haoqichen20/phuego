# -*- coding: utf-8 -*-

from .utils import denoise_square
from .utils import fisher_test
import igraph as ig
import numpy as np
import leidenalg as la
from itertools import combinations
from scipy.spatial.distance import jensenshannon


def merge_egos(network, kde_cutoff, res_folder, uniprot_to_gene,
                  supernodes,all_nodes,direction,geneset_path,fisher_geneset, 
                  fisher_threshold,):
    for kde in kde_cutoff:        
        kde_net=network.induced_subgraph(network.vs.select(name_in=all_nodes[kde]),implementation='copy_and_delete')
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
        
        fname = direction+"_supernodes_net_"+str(kde)+".txt"
        supernodes_net.write(f=res_folder+fname,format="ncol")
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

def write_modules(clustering,nodes, kde, direction, uniprot_to_gene, res_folder, 
                  geneset_path, fisher_geneset, fisher_threshold, ):
    for j in enumerate(clustering):
        f1=open(res_folder+direction+"_module_"+str(j[0])+"_cluster_"+str(kde)+".txt", "w")
        fisher_proteins=set()
        for k in j[1]:
            f1.write(k+"\t"+"\t".join(nodes[k])+"\n")
            fisher_proteins.update(nodes[k])
        
        fname = direction+"_module_"+str(j[0])+"_fisher_"+str(kde)
        fisher_test(protein_list=fisher_proteins,
                    starting_proteins=j[1],
                    threshold=fisher_threshold,
                    component=fisher_geneset,
                    path_def=res_folder,
                    uniprot_to_gene=uniprot_to_gene,
                    geneset_path=geneset_path,
                    fname=fname,
                    )