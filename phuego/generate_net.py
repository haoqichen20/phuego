
import igraph as ig
from os import path
import numpy as np
import sys,os
import os.path

def generate_nets(res_folder, network, uniprot_to_gene, kde_cutoff,
                  include_isolated_egos_in_KDE_net,net_format,):
    if os.path.isdir(res_folder):
        gene_label=[]
        for i in network.vs["name"]:
            gene_label.append(uniprot_to_gene.get(i,"Not_available"))
        network.vs["Gene_name"]=gene_label

        f1=open(res_folder+"start_seeds.txt")
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
        while(seq!=""):
            seq=seq.strip().split("\t")
            seq[1:]=np.array(seq[1:],float)
            if max(seq[1:4])>950:
                pvalues_pos.append(seq[0])
            if max(seq[4:])>950:
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

        """
        UPREGULATED
        """
        for i in kde_cutoff:
            """
            SIGNATURE NETWORK.
            """
            # The kde is a float. Convert to string for using in path.
            i = str(i)
            f1=open(res_folder+"upregulated_sig_cluster_"+i+".txt")
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
            MODULE NETWORKS
            """
            files = os.listdir(res_folder)
            module_files = [file_name for file_name in files 
                             if ("upregulated_module_" in file_name) and 
                                (("_cluster_"+i+".txt") in file_name)]
            # Get the module name, such as "module_0".
            modules = []
            for file_name in module_files:
                module = file_name.split("_")[2]
                module = "module_" + module
                modules.append(module)

            all_nodes=set()
            for module, module_file in zip(modules, module_files): 
                KDE_increased.vs[module]=list(np.full(number_of_nodes, False))
                f1 = open(res_folder+module_file)
                seq=f1.readline()
                nodes_module=set()
                while(seq!=""):
                    seq=seq.strip().split("\t")
                    nodes_module.update(seq)
                    seq=f1.readline()
                KDE_increased.vs[module]=np.in1d(KDE_increased.vs["name"], list(nodes_module))
                all_nodes.update(nodes_module)
            #write the net
            module_net=KDE_increased.induced_subgraph(all_nodes)
            module_net["title"] ="Module_increased_net"
            fname = res_folder+"module_net_increased_"+i+"."+net_format
            ig.write(module_net,fname,format=net_format)

        """
        DOWNREGULATED
        """
        for i in kde_cutoff:
            """
            SIGNATURE NETWORK
            """
            # The kde is a float. Convert to string for using in path.
            i = str(i)
            f1=open(res_folder+"downregulated_sig_cluster_"+i+".txt")
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
            MODULE NETWORKS
            """
            # Get the module file names (for downregulated and a specific kde_cutoff)
            files = os.listdir(res_folder)
            module_files = [file_name for file_name in files 
                             if ("downregulated_module_" in file_name) and 
                                (("_cluster_"+i+".txt") in file_name)]
            # Get the module name, such as "module_0".
            modules = []
            for file_name in module_files:
                module = file_name.split("_")[2]
                module = "module_" + module
                modules.append(module)
            
            all_nodes=set()
            for module, module_file in zip(modules, module_files):
                KDE_decreased.vs[module]=list(np.full(number_of_nodes, False))
                f1=open(res_folder+module_file)
                seq=f1.readline()
                nodes_module=set()
                while(seq!=""):
                    seq=seq.strip().split("\t")
                    nodes_module.update(seq)
                    seq=f1.readline()
                all_nodes.update(nodes_module)
                KDE_decreased.vs[module]=np.in1d(KDE_decreased.vs["name"], list(nodes_module))

            #write the net
            module_net=KDE_decreased.induced_subgraph(all_nodes)
            module_net["title"] ="Module_decreased_net"
            fname = res_folder+"module_net_decreased_"+i+"."+net_format
            ig.write(module_net,fname,format=net_format)
