# -*- coding: utf-8 -*-

import igraph as ig
import sys
from .utils import load_gene_names
from .utils import add_trailing_slash
from .utils import convert_result
from .load_seeds import load_seeds
from .network_rwr import rwr_values
from .network_rwr import pvalue_split
from .ego import ego_filtering
from .ego2module import merge_egos
from .generate_net import generate_nets



# one liner function of phuego package.
def phuego(support_data_folder, res_folder, test_path, 
           fisher_geneset, fisher_threshold, fisher_background,
           ini_pos, ini_neg, damping, kde_cutoff,rwr_threshold,
           use_existing_rwr = False, convert2folder = False, 
           include_isolated_egos_in_KDE_net=False,
           net_format="ncol",):
       
    """
    Formatting user input.
    """
    # Force the folder paths to end with /.
    support_data_folder = add_trailing_slash(support_data_folder)
    res_folder = add_trailing_slash(res_folder)
        
    # Force damping to be [0.5, 0.95]
    if type(damping) is not float:
           sys.exit("damping should be a float number within range [0.5, 0.95]")
    if damping < 0.5 or damping > 0.95:
           sys.exit("damping should be a float number within range [0.5, 0.95]")
    # Force rwr_threshold to be [0.01, 0.1]
    if type(rwr_threshold) is not float:
           sys.exit("rwr_threshold should be a float number within range [0.01, 0.1]")
    if rwr_threshold > 0.1 or rwr_threshold < 0.01:
           sys.exit("rwr_threshold should be a float number within range [0.01, 0.1]")
    
    # Force kde_cutoff input to be list and within [0.5, 0.95]
    if type(kde_cutoff) is not list:
           kde_cutoff = [ kde_cutoff ]
    for kde in kde_cutoff:
           if kde < 0.5 or kde > 0.95:
                  sys.exit("kde_cutoff should be within range [0.5, 0.95]")

    # Force the fisher_geneset input to be list.
    if type(fisher_geneset) is not list:
           fisher_geneset = [ fisher_geneset ]
    for geneset in fisher_geneset:
           if geneset not in ["C","F","D","P","R","K","RT","B"]:
                  sys.exit("fisher genesets should be chosen from C,F,D,P,R,K,RT,B")

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
    if(fisher_background == "proteome"):
       geneset_path = support_data_folder + "proteome/"
    elif(fisher_background == "intact"):
       geneset_path = support_data_folder + "intact/"
    else:
       sys.exit("Please provide a valid fisher background, either proteome or intact.")

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
    pvalues_pos, pvalues_neg = pvalue_split(res_folder=res_folder,
                                          seeds=seeds,
                                          graph_nodes=graph_nodes,
                                          fisher_threshold=fisher_threshold,
                                          fisher_geneset=fisher_geneset,
                                          uniprot_to_gene=uniprot_to_gene,
                                          geneset_path=geneset_path,
                                          rwr_threshold=rwr_threshold,
                                          )
    
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
                                   direction="increased",
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
                                   direction="decreased",
                                   uniprot_to_gene=uniprot_to_gene,
                                   res_folder=res_folder,
                                   geneset_path=geneset_path,
                                   fisher_geneset=fisher_geneset,
                                   fisher_threshold=fisher_threshold,
                                   )
    
    """
    Module.
    """  
    try:
       if all_nodes_pos:
              merge_egos(network=network_raw,
                            kde_cutoff=kde_cutoff,
                            res_folder=res_folder,
                            uniprot_to_gene=uniprot_to_gene,
                            direction="increased",
                            supernodes=nodes_kde_pos,
                            all_nodes=all_nodes_pos,
                            geneset_path=geneset_path,
                            fisher_geneset=fisher_geneset,
                            fisher_threshold=fisher_threshold,
                            )
    except NameError:
       print(f"No node is significantly propagated during rwr of increased seed nodes. To confirm, length of pvalues_pos is {len(pvalues_pos)}")

    try:       
       if all_nodes_neg:
              merge_egos(network=network_raw,
                            kde_cutoff=kde_cutoff,
                            res_folder=res_folder,
                            uniprot_to_gene=uniprot_to_gene,
                            direction="decreased",
                            supernodes=nodes_kde_neg,
                            all_nodes=all_nodes_neg,
                            geneset_path=geneset_path,
                            fisher_geneset=fisher_geneset,
                            fisher_threshold=fisher_threshold,
                            )
    except NameError:
       print(f"No node is significantly propagated during rwr of decreased seed nodes. To confirm, length of pvalues_neg is {len(pvalues_neg)}")
           
    """
    Output CytoScape compatible network files.
    """
    generate_nets(res_folder=res_folder,
                  network=network_raw,
                  uniprot_to_gene=uniprot_to_gene,
                  kde_cutoff=kde_cutoff,
                  rwr_threshold=rwr_threshold,
                  include_isolated_egos_in_KDE_net=include_isolated_egos_in_KDE_net,
                  net_format=net_format,
                  )
    
    """
    Convert results if required.
    """
    if(convert2folder):
           convert_result(res_folder=res_folder,
                          kde_cutoff=kde_cutoff,
                          net_format=net_format,
                          )