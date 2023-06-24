# -*- coding: utf-8 -*-

import igraph as ig
import numpy as np
from .utils import load_gene_names
from .load_seeds import load_seeds
from .ego import ego_filtering
from .ego2module import test_function

# one liner function of phuego package.
def phuego_mock(support_data_folder, res_folder, test_path, 
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
           rwr_values_mock(network=network, 
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
    pvalues_pos, pvalues_neg = pvalue_split_mock(res_folder, damping, seeds, 
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

def rwr_values_mock(network, graph_nodes, ini_pos, ini_neg, seeds, seeds_pos, 
               seeds_neg, network_path, network_random_path, damping, res_folder):
    
	number_of_nodes=network.vcount()
	# rwr_random_values={}
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

	for ii in range(10):
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
		#st = time.time()

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
	f1=open(res_folder+str(damping)+"_rwr_scores.txt","w")
	for i in empirical_values:
		f1.write(i+"\t"+"\t".join(map(str,empirical_values[i]))+"\n")
	f1.close()

	f1=open(res_folder+str(damping)+"_pvalues.txt","w")
	for i in pvalues:
		f1.write(i+"\t"+"\t".join(map(str,pvalues[i]))+"\n")
	f1.close()

	f1=open(res_folder+str(damping)+"_start_seeds.txt","w")
	for i in seeds[0]+seeds[1]+seeds[2]:
		f1.write(i+"\t"+str(seeds_pos[i])+"\n")

	for i in seeds[3]+seeds[4]+seeds[5]:
		f1.write(i+"\t"+str(-seeds_neg[i])+"\n")
	f1.close()
 
	# # Returned variables.
	# return empirical_values,pvalues
 
 
 # -*- coding: utf-8 -*-

def pvalue_split_mock(res_folder, damping, seeds, graph_nodes):
    '''
    Separate the pvalues for upregulated/downregulated nodes.
    '''
    pvalues_pos=[]
    pvalues_neg=[]
 
    f1=open(res_folder+str(damping)+"_pvalues.txt")
    seq=f1.readline()
    while(seq!=""):
        seq=seq.strip().split("\t")
        seq[1:]=np.array(seq[1:],float)
        if max(seq[1:4])>9:
            pvalues_pos.append(seq[0])
        if max(seq[4:])>9:
            pvalues_neg.append(seq[0])
        seq=f1.readline()

    pvalues_pos=pvalues_pos+seeds[0]+seeds[1]+seeds[2]
    pvalues_neg=pvalues_neg+seeds[3]+seeds[4]+seeds[5]
    
    '''
    Separate the rwr_scores for upregulated/downregulated nodes.
    '''
    f1=open(res_folder+str(damping)+"_rwr_scores.txt")
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
    
