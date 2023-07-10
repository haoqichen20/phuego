# -*- coding: utf-8 -*-

import os
import sys
import shutil
import glob

# Reorganizing output files into the folder structure of the original version.
def convert_result(res_folder, kde_cutoffs, fisher_backgrounds, test_name):
    """
    DOWN - making folders.
    """
    down_rf = res_folder+"downregulated/" 
    if os.path.exists(down_rf):
        print("Folder already exists from previous run, will be deleted now (along with files)")
        shutil.rmtree(down_rf)
    os.mkdir(down_rf)
        
    # Creating next level folder.
    down_cluster = down_rf+"cluster/"
    down_fisher = down_rf+"fisher/"
    down_networks = down_rf+"networks/"
    down_modules = down_rf+"results_modules/"
    os.mkdir(down_cluster)
    os.mkdir(down_fisher)
    os.mkdir(down_networks)
    os.mkdir(down_modules)

    # Creating the kde-level folder.
    for kde_cutoff in kde_cutoffs:
        os.mkdir(down_cluster+str(kde_cutoff))
        os.mkdir(down_networks+str(kde_cutoff))

    # Fisher folder has extra levels.
    for fisher_bg in fisher_backgrounds:
        os.mkdir(down_fisher+fisher_bg+"/")
        for kde_cutoff in kde_cutoffs:
            os.mkdir(down_fisher+fisher_bg+"/"+str(kde_cutoff))


    # Module folder has extra levels for each module. First get number of modules.
    # Note that number of modules is dependent on the kde cutoff.
    files = os.listdir(res_folder)
    down_files = [file_name for file_name in files if "downregulated_" in file_name]

    down_module_names = dict()
    for kde in kde_cutoffs:
        down_module_names[kde] = set()
        for file_name in down_files:
            if("module_" in file_name and str(kde) in file_name):
                module = file_name.split("_")[2]
                module = "module_" + module
                down_module_names[kde].add(module)
        # Make first layer folder.
        os.mkdir(down_modules+str(kde))
        for module in down_module_names[kde]:
            # Make second/third layer folder.
            os.mkdir(down_modules+str(kde)+"/"+module)
            os.mkdir(down_modules+str(kde)+"/"+module+"/fisher/")
            
    """
    DOWN - moving files.
    """
    # Move into first level.
    files = os.listdir(res_folder)
    srcf = [file_name for file_name in files if "downregulated_" in file_name]
    for file in srcf:
        newfile = file.replace("downregulated_", "")
        shutil.move(src=res_folder+file, 
                    dst=down_rf+newfile)

    # Move into second level.
    files = os.listdir(down_rf)
    labels = ["sig_cluster_", "sig_fisher_", "module_"]
    dst_folders = [down_cluster, down_fisher, down_modules]
    rep_labels = ["","","module_"]
    for label, dst_folder, rep_label in zip(labels, dst_folders, rep_labels):
        srcf = [file_name for file_name in files if label in file_name]
        for file in srcf:
            newfile = file.replace(label, rep_label)
            shutil.move(src=down_rf+file, dst=dst_folder+newfile)
            
    # Move into third level.
    # Signatures.
    for kde in kde_cutoffs:
        shutil.move(src=down_cluster+str(kde)+".txt", 
                    dst=down_cluster+str(kde)+"/"+test_name)
        
    # Networks.

    # Fishers.
    # Need to provide fisher_background to filename and use the commented line instead.
    files = os.listdir(down_fisher)
    for fisher_bg in fisher_backgrounds:
        bg_folder = down_fisher+fisher_bg+"/"
        for kde in kde_cutoffs:
            kde_folder = bg_folder+str(kde)+"/"
            srcf = [file_name for file_name in files if str(kde) in file_name]
            # srcf = [file_name for file_name in files if fisher_bg in file_name and str(kde) in file_name]
            for file in srcf:
                # newfile = file.replace(fisher_bg+"_"+str(kde)+"_", "")
                newfile = file.replace(str(kde)+"__", "")
                shutil.move(src=down_fisher+file, dst=kde_folder+newfile)

    # Modules.
    files = os.listdir(down_modules)
    for kde in kde_cutoffs:
        for module in down_module_names[kde]:
            shutil.move(src=down_modules+module+"_cluster_"+str(kde)+".txt",
                        dst=down_modules+str(kde)+"/"+module+"/"+"cluster.txt")
            # The cluster files are already moved.
            # Fisher_background need to be added here when package is updated.
            srcf = [file_name for file_name in files if (module in file_name) and (str(kde) in file_name) and ("fisher" in file_name)]
            for file in srcf:
                newfile = file.split("_")[-1]
                shutil.move(src=down_modules+file,
                            dst=down_modules+str(kde)+"/"+module+"/fisher/"+newfile)
    
    """
    UP - making folders.
    """
    up_rf = res_folder+"upregulated/" 
    if os.path.exists(up_rf):
        print("Folder already exists from previous run, will be deleted now (along with files)")
        shutil.rmtree(up_rf)
    os.mkdir(up_rf)
        
    # Creating next level folder.
    up_cluster = up_rf+"cluster/"
    up_fisher = up_rf+"fisher/"
    up_networks = up_rf+"networks/"
    up_modules = up_rf+"results_modules/"
    os.mkdir(up_cluster)
    os.mkdir(up_fisher)
    os.mkdir(up_networks)
    os.mkdir(up_modules)

    # Creating the kde-level folder.
    for kde_cutoff in kde_cutoffs:
        os.mkdir(up_cluster+str(kde_cutoff))
        os.mkdir(up_networks+str(kde_cutoff))

    # Fisher folder has extra levels.
    for fisher_bg in fisher_backgrounds:
        os.mkdir(up_fisher+fisher_bg+"/")
        for kde_cutoff in kde_cutoffs:
            os.mkdir(up_fisher+fisher_bg+"/"+str(kde_cutoff))


    # Module folder has extra levels for each module. First get number of modules.
    # Note that number of modules is dependent on the kde cutoff.
    files = os.listdir(res_folder)
    up_files = [file_name for file_name in files if "upregulated_" in file_name]

    up_module_names = dict()
    for kde in kde_cutoffs:
        up_module_names[kde] = set()
        for file_name in up_files:
            if("module_" in file_name and str(kde) in file_name):
                module = file_name.split("_")[2]
                module = "module_" + module
                up_module_names[kde].add(module)
        # Make first layer folder.
        os.mkdir(up_modules+str(kde))
        for module in up_module_names[kde]:
            # Make second/third layer folder.
            os.mkdir(up_modules+str(kde)+"/"+module)
            os.mkdir(up_modules+str(kde)+"/"+module+"/fisher/")
            
    """
    UP - moving files.
    """
    # Move into first level.
    files = os.listdir(res_folder)
    srcf = [file_name for file_name in files if "upregulated_" in file_name]
    for file in srcf:
        newfile = file.replace("upregulated_", "")
        shutil.move(src=res_folder+file, 
                    dst=up_rf+newfile)

    # Move into second level.
    files = os.listdir(up_rf)
    labels = ["sig_cluster_", "sig_fisher_", "module_"]
    dst_folders = [up_cluster, up_fisher, up_modules]
    rep_labels = ["","","module_"]
    for label, dst_folder, rep_label in zip(labels, dst_folders, rep_labels):
        srcf = [file_name for file_name in files if label in file_name]
        for file in srcf:
            newfile = file.replace(label, rep_label)
            shutil.move(src=up_rf+file, dst=dst_folder+newfile)
            
    # Move into third level.
    # Signatures.
    for kde in kde_cutoffs:
        shutil.move(src=up_cluster+str(kde)+".txt", 
                    dst=up_cluster+str(kde)+"/"+test_name)
        
    # Networks.

    # Fishers.
    # Need to provide fisher_background to filename and use the commented line instead.
    files = os.listdir(up_fisher)
    for fisher_bg in fisher_backgrounds:
        bg_folder = up_fisher+fisher_bg+"/"
        for kde in kde_cutoffs:
            kde_folder = bg_folder+str(kde)+"/"
            srcf = [file_name for file_name in files if str(kde) in file_name]
            # srcf = [file_name for file_name in files if fisher_bg in file_name and str(kde) in file_name]
            for file in srcf:
                # newfile = file.replace(fisher_bg+"_"+str(kde)+"_", "")
                newfile = file.replace(str(kde)+"__", "")
                shutil.move(src=up_fisher+file, dst=kde_folder+newfile)

    # Modules.
    files = os.listdir(up_modules)
    for kde in kde_cutoffs:
        for module in up_module_names[kde]:
            shutil.move(src=up_modules+module+"_cluster_"+str(kde)+".txt",
                        dst=up_modules+str(kde)+"/"+module+"/"+"cluster.txt")
            # The cluster files are already moved.
            # Fisher_background need to be added here when package is updated.
            srcf = [file_name for file_name in files if (module in file_name) and (str(kde) in file_name) and ("fisher" in file_name)]
            for file in srcf:
                newfile = file.split("_")[-1]
                shutil.move(src=up_modules+file,
                            dst=up_modules+str(kde)+"/"+module+"/fisher/"+newfile)