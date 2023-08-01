
import igraph as ig
import numpy as np
import pandas as pd
import os

def graph_to_df(G, seed, nodes_modules):
    edges = []
    # Network file as annotated edge list.
    for e in G.es:
        ProteinA = G.vs[e.source]["name"]
        ProteinA_GeneName = G.vs[e.source]["Gene_name"]
        ProteinB = G.vs[e.target]["name"]
        ProteinB_GeneName = G.vs[e.target]["Gene_name"]
        edge_dict = {
            "ProteinA": ProteinA,
            "ProteinB": ProteinB,
            "ProteinA_GeneName": ProteinA_GeneName,
            "ProteinB_GeneName": ProteinB_GeneName,
            "weight": e["weight"],
            "A_is_seed": ProteinA in seed,
            "B_is_seed": ProteinB in seed,
            "all_modules": []
        }
        in_any_module = False
        # Check if the edge is in the modules.
        for module_name, node_list in sorted(nodes_modules.items()):
            in_module = (ProteinA in node_list) and (ProteinB in node_list)
            edge_dict[f"is_{module_name}"] = in_module
            if in_module:
                in_any_module = True
                edge_dict["all_modules"].append(module_name.split("_")[1])
        edge_dict["inter_module"] = not in_any_module
        edges.append(edge_dict)
    df_edges = pd.DataFrame(edges)
    
    # Network node attribute file.
    nodes = []
    for v in G.vs:
        Protein = v["name"]
        Protein_GeneName = v["Gene_name"]
        node_dict = {
            "Protein": Protein,
            "Protein_GeneName": Protein_GeneName,
            "Protein_is_seed": Protein in seed,
        }
        for module_name, node_list in sorted(nodes_modules.items()):
            in_module = (Protein in node_list)
            node_dict[f"is_{module_name}"] = in_module
        nodes.append(node_dict)
    df_nodes = pd.DataFrame(nodes)
    return df_edges, df_nodes

def generate_nets(res_folder, network, uniprot_to_gene, kde_cutoff, rwr_threshold,
                  include_isolated_egos_in_KDE_net,net_format,):
    if os.path.isdir(res_folder):
        gene_label=[]
        for i in network.vs["name"]:
            gene_label.append(uniprot_to_gene.get(i,"Not_available"))
        network.vs["Gene_name"]=gene_label

        f1=open(res_folder+"start_seeds.txt")
        seq=f1.readline()
        seq=f1.readline()
        seeds_increase=[]
        seeds_decrease=[]
        while(seq!=""):
            seq= seq.strip().split("\t")
            seq[1]=float(seq[1])
            if seq[1]>0.0:
                seeds_increase.append(seq[0])
            if seq[1]<-0.0:
                seeds_decrease.append(seq[0])
            seq=f1.readline()
        f1.close()

        f1=open(res_folder+"pvalues.txt")
        pvalues_pos=[]
        pvalues_neg=[]
        seq=f1.readline()
        seq=f1.readline()
        while(seq!=""):
            seq=seq.strip().split("\t")
            seq[1:]=np.array(seq[1:],float)
            if max(seq[1:4]) > 1000 - rwr_threshold*1000:
                pvalues_pos.append(seq[0])
            if max(seq[4:]) > 1000 - rwr_threshold*1000:
                pvalues_neg.append(seq[0])
            seq=f1.readline()
        f1.close()
  
        """
        RWR NETWORKS -- down/upregulated.
        """
        pvalues_pos=pvalues_pos+seeds_increase
        pvalues_neg=pvalues_neg+seeds_decrease
        rwr_net_increased=network.induced_subgraph(pvalues_pos)
        rwr_net_increased["title"] ="RWR_increased_net"
        number_of_nodes=rwr_net_increased.vcount()
        rwr_net_increased.vs["Is_seed"]=list(np.full(number_of_nodes, False))
        for i in rwr_net_increased.vs:
            if i["name"] in seeds_increase:
                rwr_net_increased.vs[i.index]["Is_seed"]=True
        # Output the increased network.
        fname = res_folder+"rwr_increased."+net_format
        ig.write(rwr_net_increased,fname,format=net_format)


        rwr_net_decreased=network.induced_subgraph(pvalues_neg)
        rwr_net_decreased["title"] ="RWR_decreased_net"

        number_of_nodes=rwr_net_decreased.vcount()
        rwr_net_decreased.vs["Is_seed"]=list(np.full(number_of_nodes, False))
        for i in rwr_net_decreased.vs:
            if i["name"] in seeds_decrease:
                rwr_net_decreased.vs[i.index]["Is_seed"]=True
        # Output the decreased network.
        fname = res_folder+"rwr_decreased."+net_format
        ig.write(rwr_net_decreased,fname,format=net_format)


        for i in kde_cutoff:
            # The kde is a float. Convert to string for using in path.
            i = str(i)
            # Depending on the input, the sig_cluster file might not exist for increased seed nodes.
            if os.path.isfile(res_folder+"increased_sig_cluster_"+i+".txt"):
                """
                UPREGULATED -- SIGNATURE NETWORK.
                """
                f1=open(res_folder+"increased_sig_cluster_"+i+".txt")
                seq=f1.readline()
                nodes=set()
                while(seq!=""):
                    seq=seq.strip().split("\t")
                    if include_isolated_egos_in_KDE_net==False:
                        if len(seq)>2:
                            nodes.update(seq)
                    else:
                        nodes.update(seq)
                    seq=f1.readline()
                KDE_increased=rwr_net_increased.induced_subgraph(nodes)
                KDE_increased["title"] ="KDE_increased_net"
                number_of_nodes=KDE_increased.vcount()
                ##write the net
                fname = res_folder+"KDE_increased_"+i+"."+net_format
                ig.write(KDE_increased,fname,format=net_format)
   
                """
                UPREGULATED -- MODULE NETWORKS
                """ 
                files = os.listdir(res_folder)
                module_files = [file_name for file_name in files 
                                if ("increased_module_" in file_name) and 
                                    (("_cluster_"+i+".txt") in file_name)]
                # Get the module name, such as "module_0".
                modules = []
                for file_name in module_files:
                    module = file_name.split("_")[2]
                    module = "module_" + module
                    modules.append(module)

                all_nodes=set()
                nodes_module=dict()
                for module, module_file in zip(modules, module_files): 
                    KDE_increased.vs[module]=list(np.full(number_of_nodes, False))
                    f1 = open(res_folder+module_file)
                    seq=f1.readline()
                    nodes_module[module]=set()
                    while(seq!=""):
                        seq=seq.strip().split("\t")
                        nodes_module[module].update(seq)
                        seq=f1.readline()
                    # Label the KDE_increased vertices attribute with module name.
                    # This vertice attribute will be inherited to module_net when subgraph is induced.
                    KDE_increased.vs[module]=np.in1d(KDE_increased.vs["name"], list(nodes_module[module]))
                    all_nodes.update(nodes_module[module])
                #write the net
                module_net=KDE_increased.induced_subgraph(all_nodes)
                module_net["title"] ="Module_increased_net"
                fname = res_folder+"module_net_increased_"+i+"."+net_format
                ig.write(module_net,fname,format=net_format)
                
                # Create the dataframe for annotated csv output of module network.
                seed = seeds_increase + seeds_decrease
                df_module_net_edges, df_module_net_nodes = graph_to_df(module_net, seed, nodes_module)
                df_module_net_edges.to_csv(res_folder+"module_net_increased_edgelist_"+i+".csv")
                df_module_net_nodes.to_csv(res_folder+"module_net_increased_nodes_attribute_"+i+".csv")
            else:
                pass

            # Depending on the input, the sig_cluster file might not exist for decreased seed nodes.
            if os.path.isfile(res_folder+"decreased_sig_cluster_"+i+".txt"):
                """
                DOWNREGULATED -- SIGNATURE NETWORK
                """
                f1=open(res_folder+"decreased_sig_cluster_"+i+".txt")
                seq=f1.readline()
                nodes=set()
                while(seq!=""):
                    seq=seq.strip().split("\t")
                    if include_isolated_egos_in_KDE_net==False:
                        if len(seq)>2:
                            nodes.update(seq)
                    else:
                        nodes.update(seq)
                    seq=f1.readline()
                KDE_decreased=rwr_net_decreased.induced_subgraph(nodes)
                KDE_decreased["title"] ="KDE_decreased_net"
                number_of_nodes=KDE_decreased.vcount()
                ##write the net
                fname = res_folder+"KDE_decreased_"+i+"."+net_format
                ig.write(KDE_decreased,fname,format=net_format)

                """
                DOWNREGULATED -- MODULE NETWORKS
                """ 
                # Get the module file names (for downregulated and a specific kde_cutoff)
                files = os.listdir(res_folder)
                module_files = [file_name for file_name in files 
                                if ("decreased_module_" in file_name) and 
                                    (("_cluster_"+i+".txt") in file_name)]
                # Get the module name, such as "module_0".
                modules = []
                for file_name in module_files:
                    module = file_name.split("_")[2]
                    module = "module_" + module
                    modules.append(module)
                
                all_nodes=set()
                nodes_module=dict()
                for module, module_file in zip(modules, module_files):
                    KDE_decreased.vs[module]=list(np.full(number_of_nodes, False))
                    f1=open(res_folder+module_file)
                    seq=f1.readline()
                    nodes_module[module]=set()
                    while(seq!=""):
                        seq=seq.strip().split("\t")
                        nodes_module[module].update(seq)
                        seq=f1.readline()
                    all_nodes.update(nodes_module[module])
                    KDE_decreased.vs[module]=np.in1d(KDE_decreased.vs["name"], list(nodes_module[module]))

                #write the net
                module_net=KDE_decreased.induced_subgraph(all_nodes)
                module_net["title"] ="Module_decreased_net"
                fname = res_folder+"module_net_decreased_"+i+"."+net_format
                ig.write(module_net,fname,format=net_format)
                
                # Create the dataframe for annotated csv output of module network.
                seed = seeds_increase + seeds_decrease
                df_module_net_edges, df_module_net_nodes = graph_to_df(module_net, seed, nodes_module)
                df_module_net_edges.to_csv(res_folder+"module_net_decreased_edgelist_"+i+".csv")
                df_module_net_nodes.to_csv(res_folder+"module_net_decreased_nodes_attribute_"+i+".csv")
            else:
                pass
