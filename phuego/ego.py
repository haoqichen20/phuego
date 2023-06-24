# -*- coding: utf-8 -*-

from .utils import denoise_square
from .utils import calc_kde
from .fishertest import fisher_test
import numpy as np
import igraph
from scipy.spatial.distance import jensenshannon

def ego_filtering(network, pval, seeds, sim, zscores_global, kde_cutoff, 
                  direction, uniprot_to_gene, res_folder,
                  geneset_path, fisher_geneset, fisher_threshold):
    # Below: only add kde_cutoff to ego_friends input.
    subnet = network.induced_subgraph(list(network.vs.select(name_in=pval)),implementation='create_from_scratch')
    subnet.vs.select(_degree=0).delete()
    nodes=subnet.vs["name"]
    n_nodes=subnet.vcount()
    position_nodes={k: v for v, k in enumerate(nodes)}
    seed_nodes=list(set(seeds).intersection(set(nodes)))
    nodes_kde=ego_friends(subnet=subnet,
                          position_nodes=position_nodes,
                          seed_nodes=seed_nodes,
                          sim=sim,
                          zscores=zscores_global,
                          kde_cutoff=kde_cutoff,
                          )
    
    # Fisher test and write results.
    supernodes,all_nodes=write_results(nodes_kde=nodes_kde,
                                       seed_nodes=seed_nodes, 
                                       kde_cutoff=kde_cutoff, 
                                       direction=direction, 
                                       uniprot_to_gene=uniprot_to_gene, 
                                       res_folder=res_folder,
                                       geneset_path=geneset_path, 
                                       fisher_geneset=fisher_geneset, 
                                       fisher_threshold=fisher_threshold,
                                       )
    return supernodes,all_nodes

# Replaced [0.85, 0.9] with kde_cutoff, in two places.
# Did not touch damping.
def ego_friends(subnet,position_nodes,seed_nodes,sim,zscores,kde_cutoff):

    nodes_kde={}
    for i in kde_cutoff:
        nodes_kde[i]={}

    for i in enumerate(sorted(seed_nodes)):
        ego=subnet.induced_subgraph(subnet.neighborhood(vertices=position_nodes[i[1]], order=2, mode='all', mindist=0),implementation='create_from_scratch')
        nodes=ego.vs["name"]
        n_nodes=len(nodes)
        if n_nodes>5:
            position_ego={k: v for v, k in enumerate(nodes)}
            functional_nodes=[position_ego[i[1]]]

            for j in nodes:
                num=sim[i[1]][j]-float(zscores[j][0])
                den=float(zscores[j][1])
                if den>0 and (num/den)>1.645:
                    if i[1]!=j:
                        functional_nodes.append(position_ego[j])
            ego=ego.induced_subgraph(functional_nodes,implementation='copy_and_delete')
            ego=ego.induced_subgraph(ego.neighborhood(vertices=ego.vs.find(i[1]), order=2, mode='all', mindist=0),implementation='copy_and_delete')
            ego.vs.select(_degree=0).delete()
            nodes=ego.vs["name"]
            n_nodes=len(nodes)
            position_ego={k: v for v, k in enumerate(nodes)}
            if len(ego.connected_components(mode='strong'))>1 and i[1] in position_ego:
                ego=ego.induced_subgraph(ego.subcomponent(position_ego[i[1]], mode='all'),implementation='copy_and_delete')
                ego.vs.select(_degree=0).delete()
                nodes=ego.vs["name"]
                n_nodes=len(nodes)
                position_ego={k: v for v, k in enumerate(nodes)}

            if i[1] in position_ego and n_nodes>5:
                position_ego={k: v for v, k in enumerate(nodes)}
                first_shell=ego.neighbors(position_ego[i[1]], mode='all')

                distances=dict.fromkeys(first_shell,1)
                nodes_id=ego.vs.indices
                second_shell=list(set(nodes_id).difference(set(first_shell+[position_ego[i[1]]])))
                for j in second_shell:
                    distances[j]=2
                distances[position_ego[i[1]]]=0
                for j in ego.es:
                    node_A=j.tuple[0]
                    node_B=j.tuple[1]
                    dist=[distances[node_A],distances[node_B]]
                    if (dist[0]==1 and dist[1]==1) or (dist[0]==2 and dist[1]==2):
                        ego_sim= (sim[i[1]][nodes[node_A]]+sim[i[1]][nodes[node_B]])/2.0
                        ego.es[j.index]["weight"]=ego_sim
                    elif (dist[0]==1 and dist[1]==2):
                        ego_sim=sim[i[1]][nodes[node_B]]
                        ego.es[j.index]["weight"]=ego_sim
                    elif (dist[0]==2 and dist[1]==1):
                        ego_sim=sim[i[1]][nodes[node_A]]
                        ego.es[j.index]["weight"]=ego_sim
                ego=denoise_square(ego)

                reset_vertex=np.zeros(n_nodes)
                reset_vertex[position_ego[i[1]]]=1.0
                ego_rwr=np.array(ego.personalized_pagerank(reset=reset_vertex,directed=False, damping=0.85, weights='weight',implementation='prpack'))
                ssim=[]
                dist=[]

                for j in ego.vs:
                    reset_vertex=np.zeros(n_nodes)
                    reset_vertex[j.index]=1.0
                    ego_node=np.array(ego.personalized_pagerank(reset=reset_vertex,directed=False, damping=0.85,implementation='prpack',weights='weight'))
                    if i[1]==j["name"]:
                        dist.append(1.0)
                    else:
                        dist.append(1-jensenshannon(ego_node,ego_rwr))
                        #print (1000*np.log2(1+dist),1000*np.log2(1+1-jsd)
                    ssim.append(sim[i[1]][j["name"]])

                ssim=np.array(ssim)
                dist=np.array(dist)
                ssim=1000*np.log2(1+ssim)
                dist=1000*np.log2(1+dist)

                cdf_dist,grid_dist=calc_kde(dist)
                cdf_ssim,grid_ssim=calc_kde(ssim)

                for j in kde_cutoff:
                    position=np.searchsorted(cdf_dist,j)
                    interpolation=np.argmin([abs(j-cdf_dist[position]),abs(j-cdf_dist[position-1])])
                    grid_threshold_dist=grid_dist[position-interpolation]

                    position=np.searchsorted(cdf_ssim,j)
                    interpolation=np.argmin([abs(j-cdf_ssim[position]),abs(j-cdf_ssim[position-1])])
                    grid_threshold_ssim=grid_ssim[position-interpolation]

                    nodes_sim=[i[1]]
                    nodes_dist=[i[1]]
                    for k in enumerate(nodes):
                        if ssim[k[0]]>grid_threshold_ssim:
                            nodes_sim.append(k[1])
                        if dist[k[0]]>grid_threshold_dist:
                            nodes_dist.append(k[1])
                    nodes_kde[j][i[1]]=list(set(nodes_sim).intersection(nodes_dist))

    return nodes_kde

def write_results(nodes_kde, seed_nodes, kde_cutoff, direction, uniprot_to_gene, 
                  res_folder, geneset_path, fisher_geneset, fisher_threshold):
    all_nodes={}
    supernodes={}
    for j in kde_cutoff:
        all_nodes[j]=set()
        supernodes[j]={}

        # Variable temp is removed. It's just fisher_proteins as a set not a list.
        f2=open(res_folder+direction+"_sig"+"_cluster_"+str(j)+".txt", "w")
        fisher_proteins=set()
        for k in nodes_kde[j]:
            f2.write(k+"\t"+"\t".join(nodes_kde[j][k])+"\n")
            if len(nodes_kde[j][k])>=2:
                supernodes[j][k]=nodes_kde[j][k]
                fisher_proteins.update(nodes_kde[j][k])

        if fisher_proteins:
            # Create filename for fisher output.
            fname = direction+"_sig_fisher_"+str(j)+"_"
            fisher_test(protein_list=fisher_proteins,
                        starting_proteins=seed_nodes,
                        fname=fname,
                        threshold=fisher_threshold,
                        component=fisher_geneset,
                        path_def=res_folder,
                        uniprot_to_gene=uniprot_to_gene,
                        geneset_path=geneset_path,
                        )
            all_nodes[j]=list(fisher_proteins)
    return supernodes,all_nodes