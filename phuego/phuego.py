# -*- coding: utf-8 -*-

import igraph as ig
from .utils import load_gene_names
from .load_seeds import load_seeds
from .rwr_values import rwr_values
from .pvalue_splitting import pvalue_split
from .ego import ego_filtering
from .ego2module import test_function

# one liner function of phuego package.
def phuego(support_data_folder, res_folder, test_path, 
           fisher_geneset, fisher_threshold,
           ini_pos, ini_neg, damping, kde_cutoff,
           use_existing_rwr = False):
    """
    Create support data paths.
    """       
    # Networks
    network_path = support_data_folder + "networks/gic.txt"
    network_raw_path = support_data_folder + "networks/gic_raw.txt"
    network_random_path = support_data_folder + "networks/gic_random/"
    # Protein annotation, ID mapping
    pfam_domain_path = support_data_folder + "pfam_domains.txt"
    uniprot_to_gene_path = support_data_folder + "uniprot_to_gene.tab"
    # Semantic similarity
    sim_mean_std_path = support_data_folder + "gic_mean_std.txt"
    sim_all_folder_path = support_data_folder + "gic_sim/"
    # Genesets.
    geneset_path = support_data_folder + "proteome/"
    
    """
    Formatting user input.
    """
    # Force kde_cutoff input to be list.
    if(type(kde_cutoff) is not list):
           kde_cutoff = [ kde_cutoff ]
    # Force the fisher_geneset input to be list.
    if(type(fisher_geneset) is not list):
           fisher_geneset = [ fisher_geneset ]
    # ToDo: check if all folders end with /.
    # ToDo: check if folders/paths are compatible with windows.

    """
    load network.
    """
    network = ig.Graph.Read_Ncol(network_path, weights=True, directed=False)
    graph_nodes=network.vs["name"]
    
    network_raw = ig.Graph.Read_Ncol(network_raw_path, weights=True, directed=False)
    graph_nodes_raw = network_raw.vs["name"]
    
    """
    Load uniprot ID to gene name mapping.
    """
    uniprot_to_gene=load_gene_names(uniprot_to_gene_path)
    
    """
    load seeds.
    """
    seeds_pos,seeds_neg,seeds,zscores_global,ssim = load_seeds(
                     pfam_domain_path=pfam_domain_path,
                     sim_mean_std_path=sim_mean_std_path,
                     sim_all_folder_path=sim_all_folder_path,
                     test_path=test_path,
                     graph_nodes=graph_nodes,
                     )
    
    """
    Random walk with restart
    """
    # On genuine and randomized netowrk.
    # TIME CONSUMING: skip and use existing result by 'use_existing_rwr = True'
    if(use_existing_rwr):
           pass
    else:
           rwr_values(network=network, 
                      graph_nodes=graph_nodes, 
                      ini_pos=ini_pos,
                      ini_neg=ini_neg,
                      seeds=seeds,
                      seeds_pos=seeds_pos,
                      seeds_neg=seeds_neg,
                      network_path=network_path,
                      network_random_path=network_random_path,
                      damping=damping,
                      res_folder=res_folder,
                     )
    
    """
    Split the pvalues.
    """
    # pvalues is imported from file in res_folder.
    pvalues_pos, pvalues_neg = pvalue_split(res_folder, damping, seeds, 
              graph_nodes)
    
    """
    EGO.
    """
    if pvalues_pos:
           nodes_kde_pos, all_nodes_pos = ego_filtering(network=network,
                                   pval=pvalues_pos,
                                   seeds=seeds_pos,
                                   sim=ssim,
                                   zscores_global=zscores_global,
                                   kde_cutoff=kde_cutoff,
                                   direction="upregulated",
                                   uniprot_to_gene=uniprot_to_gene,
                                   res_folder=res_folder,
                                   geneset_path=geneset_path,
                                   fisher_geneset=fisher_geneset,
                                   fisher_threshold=fisher_threshold,
                                   )
    if pvalues_neg:
           nodes_kde_neg, all_nodes_neg = ego_filtering(network=network,
                                   pval=pvalues_neg,
                                   seeds=seeds_neg,
                                   sim=ssim,
                                   zscores_global=zscores_global,
                                   kde_cutoff=kde_cutoff,
                                   direction="downregulated",
                                   uniprot_to_gene=uniprot_to_gene,
                                   res_folder=res_folder,
                                   geneset_path=geneset_path,
                                   fisher_geneset=fisher_geneset,
                                   fisher_threshold=fisher_threshold,
                                   )
    
    """
    Module.
    """
    if all_nodes_pos:
           test_function(network=network_raw,
                         network_nodes=graph_nodes_raw,
                         kde_cutoff=kde_cutoff,
                         res_folder=res_folder,
                         uniprot_to_gene=uniprot_to_gene,
                         direction="upregulated",
                         supernodes=nodes_kde_pos,
                         all_nodes=all_nodes_pos,
                         geneset_path=geneset_path,
                         fisher_geneset=fisher_geneset,
                         fisher_threshold=fisher_threshold,
                         )
           
    if all_nodes_neg:
           test_function(network=network_raw,
                         network_nodes=graph_nodes_raw,
                         kde_cutoff=kde_cutoff,
                         res_folder=res_folder,
                         uniprot_to_gene=uniprot_to_gene,
                         direction="downregulated",
                         supernodes=nodes_kde_neg,
                         all_nodes=all_nodes_neg,
                         geneset_path=geneset_path,
                         fisher_geneset=fisher_geneset,
                         fisher_threshold=fisher_threshold,
                         )

    """
    Output CytoScape compatible network files.
    """
    number_of_nodes=network.vcount()
    number_of_genes = len(list(uniprot_to_gene.keys()))
    return(number_of_nodes, number_of_genes)