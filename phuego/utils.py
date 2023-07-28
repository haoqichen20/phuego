#!/usr/bin/env python

import os
import sys
import shutil
import numpy as np
import math
from sklearn.neighbors import KernelDensity


def load_gene_names(uniprot_to_gene_path):
	uniprot_to_gene={}
	gene_to_uniprot={}
	localization={}
	f1=open(uniprot_to_gene_path,"r")
	seq=f1.readline()
	while(seq!=""):
		seq=seq.strip().split("\t")
		uniprot_to_gene[seq[0]]=seq[1].split(";")[0].strip()
		seq=f1.readline()
	return uniprot_to_gene


def denoise_square(G):
	weight=G.strength(G.vs, mode='all', loops=False, weights='weight')
	for i in G.es():
		node_A=i.tuple[0]
		node_B=i.tuple[1]
		den=math.sqrt(weight[node_A]*weight[node_B])
		num=i["weight"]
		G.es[i.index]["weight"]=num/den
	return (G)


def calc_kde(vector):
	obs = len(vector)
	sigma = np.std(vector, ddof=1)
	IQR = (np.percentile(vector, q=75) - np.percentile(vector, q=25)) / 1.3489795003921634
	sigma = min(sigma, IQR)
	if sigma > 0:
		bw= sigma * (obs * 3 / 4.0) ** (-1 / 5)
	else:
		IQR = (np.percentile(vector, q=99) - np.percentile(vector, q=1)) / 4.6526957480816815
		if IQR > 0:
			bw = IQR * (obs * 3 / 4.0) ** (-1 / 5)
	if bw:
		kde = KernelDensity(kernel = 'gaussian', bandwidth=bw).fit(vector.reshape(-1,1))
		grid = np.linspace(min(vector)-1,max(vector)+1,len(vector)*100).reshape(-1,1)
		log_dens = kde.score_samples(grid)
		pdf=np.exp(log_dens)
		grid=grid.ravel()
		normalization=sum(pdf.ravel())
		cdf=np.cumsum(pdf)/normalization
	else:
		onevec=np.ones(len(vector))
		return onevec,onevec

	return cdf,grid


# Check if folder path end with forward slash, if not, add it.
def add_trailing_slash(folder_path):
    if not folder_path.endswith('/'):
        folder_path = folder_path + '/'
    return folder_path


def convert_result(res_folder, kde_cutoff, net_format):
    # Create direction folder and move files.
    files = os.listdir(res_folder)
    for direction in ["decreased", "increased"]:
        # Get file names with direction flag.
        # rwr result files without direction flag will remain in res_folder.
        direction_files = [file_name for file_name in files if direction in file_name]
        
        # Create direction master folder.
        direction_rf = res_folder+direction+"/"
        if os.path.exists(direction_rf):
            print(f"Result sub-folder {direction} already exists from previous run, will be deleted now (along with files)")
            shutil.rmtree(direction_rf)
        os.mkdir(direction_rf)
        # Move rwr network file.
        shutil.move(src=res_folder+"rwr_"+direction+"."+net_format, 
                    dst=direction_rf+"rwr_net."+net_format)
        
        # Create seed fisher folder.
        direction_seed_fisher_rf = direction_rf+"seed_fisher/"
        os.mkdir(direction_seed_fisher_rf)
        # Move all fisher files for seeds.
        srcf = [file_name for file_name in direction_files if "_seed_fisher_" in file_name]
        for file in srcf:
            new_file = file.split("_")[-1]
            shutil.move(src=res_folder+file, dst=direction_seed_fisher_rf+new_file)
            
        # Create rwr fisher folder.
        direction_rwr_fisher_rf = direction_rf+"rwr_fisher/"
        os.mkdir(direction_rwr_fisher_rf)
        # Move all fisher files for rwr.
        srcf = [file_name for file_name in direction_files if "_rwr_fisher_" in file_name]
        for file in srcf:
            new_file = file.split("_")[-1]
            shutil.move(src=res_folder+file, dst=direction_rwr_fisher_rf+new_file)            
        
        
        # Create subfolders and move files.
        for kde in kde_cutoff:
            ######## Create the kde-level folder ########
            direction_KDE_rf = direction_rf+"KDE_"+str(kde)+"/"
            os.mkdir(direction_KDE_rf)
            # Move KDE_egos.txt.
            srcf = [file_name for file_name in direction_files if "_sig_cluster_" in file_name and str(kde) in file_name]
            # there should only be one _sig_cluster_ file per kde x direction.
            if(len(srcf)>1):
                sys.exit(f"More than one _sig_cluster_ file is present for '{direction}'. KDE is {kde}")
            for file in srcf:
                new_file = "KDE_egos.txt"
                shutil.move(src=res_folder+file, dst=direction_KDE_rf+new_file)
            
            ######## Create fisher folder ########
            direction_KDE_fisher_rf = direction_KDE_rf + "fisher/"
            os.mkdir(direction_KDE_fisher_rf)
            # Move all fisher files for signature.
            srcf = [file_name for file_name in direction_files if "_sig_fisher_" in file_name and str(kde) in file_name]
            for file in srcf:
                new_file = file.split("_")[-1]
                shutil.move(src=res_folder+file, dst=direction_KDE_fisher_rf+new_file)
                
            ######## Create network folder ########
            direction_KDE_networks_rf = direction_KDE_rf + "networks/"
            os.mkdir(direction_KDE_networks_rf)
            # Move supernodes_net.
            srcf = [file_name for file_name in direction_files if "_supernodes_net_" in file_name and str(kde) in file_name]
            # There should only be one supernodes_net per kde x direction.
            if(len(srcf)>1):
                sys.exit(f"More than one supernodes_net is present for '{direction}'. KDE is {kde}")
            for file in srcf:
                new_file = "supernodes_net.txt"
                shutil.move(src=res_folder+file, dst=direction_KDE_networks_rf+new_file)
            # Move sig_net and module_net
            srcf = [file_name for file_name in direction_files if net_format in file_name and str(kde) in file_name]
            for file in srcf:
                new_file = file.replace("_"+direction+"_"+str(kde), "")
                shutil.move(src=res_folder+file, dst=direction_KDE_networks_rf+new_file)
            # Move module_net csv file.
            srcf = [file_name for file_name in direction_files if "module_net_" in file_name and ".csv" in file_name and str(kde) in file_name]
            for file in srcf:
                # new_file = "module_net_annotated_edgelist.csv"
                new_file = file.replace("_"+direction, "")
                new_file = new_file.replace("_"+str(kde), "")
                shutil.move(src=res_folder+file, dst=direction_KDE_networks_rf+new_file)
            
            ######## Create module folders ########
            direction_KDE_modules_rf = direction_KDE_rf + "modules/"
            os.mkdir(direction_KDE_modules_rf)
            # Collect the module names for the regulation direction.
            direction_module_names = set()
            for file_name in direction_files:
                if("module_" in file_name and "cluster_"+str(kde) in file_name):
                    module = file_name.split("_")[2]
                    module = "module_" + module
                    direction_module_names.add(module)
            
            # Module folder has extra levels for each module. 
            for module in direction_module_names:
                # Module folder.
                direction_KDE_modules_module_rf = direction_KDE_modules_rf+module+"/"
                os.mkdir(direction_KDE_modules_module_rf)
                # Move module_egos file.
                srcf = [file_name for file_name in direction_files if module in file_name and "_cluster_" in file_name and str(kde) in file_name]
                for file in srcf:
                    new_file = "module_egos.txt"
                    shutil.move(res_folder+file, direction_KDE_modules_module_rf+new_file)
                
                # Module fisher folder.
                direction_KDE_modules_module_fisher_rf = direction_KDE_modules_rf+module+"/fisher/"
                os.mkdir(direction_KDE_modules_module_fisher_rf)
                # Move module fisher files.
                srcf = [file_name for file_name in direction_files if module in file_name and "_fisher_" in file_name and str(kde) in file_name]
                for file in srcf:
                    new_file = file.split("_")[-1]
                    shutil.move(res_folder+file, direction_KDE_modules_module_fisher_rf+new_file)