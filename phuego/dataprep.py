# -*- coding: utf-8 -*-

import requests
import hashlib
import tarfile
import os
import sys

def dataprep(support_data_folder, need_fisher=True, need_gic_sim=True, need_networks=True, 
             remove_zip_file = False):
    
    # Input retrieved from Zenodo. Should be static unless files are changed by developer.
    project_url = "https://zenodo.org/api/records/8040791"
    fisher_md5 = "d9f871680c2cc07bbaf35a2b42083522"
    gic_sim_md5 = "fa608665e4e6ce2241f5b5dcf105c552"
    networks_md5 = "6c249ac5c049815eb7a3fadcbb5c149c"
    
    # Send a GET request to the download URL
    response = requests.get(project_url)
    if response.status_code == 200:
        data = response.json()
        fisher_url = data["files"][0]["links"]["self"]
        gic_sim_url = data["files"][1]["links"]["self"]
        networks_url = data["files"][2]["links"]["self"]
    else:
        sys.exit("Cannot access Zenodo. Please troubleshoot.")
        
    if(need_fisher):
        download_dataset(fisher_url, support_data_folder, fisher_md5, remove_zip_file)
    if(need_gic_sim):
        download_dataset(gic_sim_url, support_data_folder, gic_sim_md5, remove_zip_file)
    if(need_networks):
        download_dataset(networks_url, support_data_folder, networks_md5, remove_zip_file)

def download_dataset(download_url, support_data_folder, file_md5, remove_zip_file):
    # Send a GET request to the download URL
    response = requests.get(download_url)

    if response.status_code == 200:
        # Extract the filename from the response headers
        filename = response.headers.get("Content-Disposition").split("=")[1]
        
        # Save the dataset to a file
        with open(support_data_folder+filename, "wb") as file:
            file.write(response.content)
        
        # Check if md5_checksum matches. If so decompress file.
        md5_checksum = calculate_md5(support_data_folder+filename)
        if(file_md5 == md5_checksum):
            decompress_tar_gz(support_data_folder+filename, support_data_folder)
            print(f"MD5 checksum match, dataset '{filename}' downloaded successfully!")
        else:
            sys.exit("MD5 checksum doesn't match, download incomplete.")

        # Remove the tar.gz file after decompression, if required by the user.    
        if(remove_zip_file):
            os.remove(support_data_folder+filename)
            print(f"{filename} has been removed after decompression.")
    else:
        sys.exit("Cannot access Zenodo. Please troubleshoot.")
        

def calculate_md5(file_path):
    md5_hash = hashlib.md5()
    with open(file_path, "rb") as file:
        # Read the file in chunks to conserve memory
        for chunk in iter(lambda: file.read(4096), b""):
            md5_hash.update(chunk)
    return md5_hash.hexdigest()


def decompress_tar_gz(file_path, output_dir):
    with tarfile.open(file_path, "r:gz") as tar:
        tar.extractall(path=output_dir)