# -*- coding: utf-8 -*-

from .utils_CLI import load_test_example
from .utils_CLI import dataprep
from .phuego import phuego
from .phuego_mock import phuego_mock
import click

__version__ = "1.0.0"

'''
This is the CLI tool for phuego.
'''
@click.command()
# Folders.
@click.option("--support_data_folder", "-sf", type=str, required=False, 
              help="Folder that stored support data. Directories separated by /")
@click.option("--res_folder", "-rf", type=str, required=False, 
              help="Folder that stored output result. Directories separated by /")
@click.option("--test_path", "-tpath", default="", type=str, required=False, 
              help="The path of user input test file. If defined, -mock and -test should not be used.")
@click.option("--ini_path", "-ipath", default="", type=str, required=False, 
              help="The path of a csv file, that specify two sets of nodes to be removed from reference network during propagation.")
@click.option("--convert2folder/--dont_convert2folder", default=True, required=False, 
              help="Should phuego organize result files into a folder structure?")


# Support dataset download.
@click.option("--need_fisher", is_flag=True, required=False,
              help="Should phuego download the fisher geneset from zenodo?")
@click.option("--need_gic_sim", is_flag=True, required=False, 
              help="Should phuego download the semantic similarity from zenodo?")
@click.option("--need_networks", is_flag=True, required=False, 
              help="Should phuego download the networks from zenodo?")
@click.option("--remove_zip_file", is_flag=True, required=False, 
              help="Should phuego remove the zip files after downloading?")

# Mock/test.
@click.option("--run_mock", is_flag=True, required=False, 
              help="Should phuego perform a mock test run using test data? If defined, -test and -tpath would be ignored.")
@click.option("--run_test", is_flag=True, required=False, 
              help="Should phuego perform a full test run using test data? If defined, -tpath would be ignored.")


# network propagation parameters.
@click.option("--use_existing_rwr", "-ru", is_flag=True, required=False, 
              help="Should phuego reuse existing network propagation results?")
@click.option("--damping", "-d", default=0.85, type=float, required=False, 
              help="Damping factor for random walk with restart algorithm. Float number within range [0.5, 0.95]")
@click.option("--rwr_threshold", "-rt", default=0.05, type=float, required=False, 
              help="Threshold of significance for propagated nodes. Float number within range [0.01, 0.1]")

# ego propagation and fisher test.
@click.option("--kde_cutoff", "-k", default=[0.85], type=float, multiple=True,
              required=False, 
              help="KDE cutoff value for removing less similar nodes in ego network. A float number within range [0.5, 0.95]. Multiple numbers can be provided at the same time by reusing the argument flag.")
@click.option("--fisher_geneset", "-fg", default=["K"], type=str, multiple=True,
              required=False, 
              help="Abbreviation of genesets to be tested when annotating modules of propagation, can be one of 'C,F,D,P,R,K,RT,B'. Multiple genesets can be provided at the same time by reusing the argument flag. Refer to documents to learn what each abbreviation stands for.")
@click.option("--fisher_threshold", "-ft", default=0.05, type=float, 
              required=False, 
              help="Threshold of significance of geneset enrichment analysis") 
@click.option("--fisher_background", "-fb", default="intact", type=str, 
              required=False, 
              help="Geneset background of geneset enrichment analysis") 

# network output
@click.option("--include_isolated_egos_in_kde_net", "-ie", is_flag=True, required=False, 
              help="Should we include isolated nodes in the output network?") 
@click.option("--net_format", "-nf", default="graphml", type=str, required=False, 
              help="file format of output network. Can be one of 'edgelist', 'pajek', 'ncol', 'lgl', 'graphml', 'dimacs', 'gml', 'dot', 'leda'")

# Versioning.
@click.option('--version', '-v', is_flag=True, 
              help='Print version to stdout')

def main(support_data_folder, res_folder, test_path, convert2folder, use_existing_rwr, run_test, run_mock, 
         need_fisher, need_gic_sim, need_networks, remove_zip_file,
         damping, fisher_geneset, fisher_threshold, fisher_background, kde_cutoff, ini_path, rwr_threshold, 
         include_isolated_egos_in_kde_net, net_format, version) -> None:
       
    # Print the version number.
    if version:
        bold_package_name = click.style("phuego", bold=True)
        click.echo(f"{bold_package_name} ({__version__})")
        return

    # Download support dataset.
    if(need_fisher | need_gic_sim | need_networks):
       dataprep(support_data_folder=support_data_folder, 
              need_fisher=need_fisher, 
              need_gic_sim=need_gic_sim, 
              need_networks=need_networks,
              remove_zip_file=remove_zip_file)
       

    # Click turns the multiple input into tuple. Make it as list now.
    fisher_geneset = list(fisher_geneset)
    kde_cutoff = list(kde_cutoff)
    # ini_pos = list(ini_pos)
    # ini_neg = list(ini_neg)
        
    # Create ini_pos and ini_neg, by default or by reading the ini_path file. 
    if(ini_path == ""):
       ini_pos = ["False"]
       ini_neg = ["False"]
    else:
       f1=open(ini_path)
       ini_pos = f1.readline().strip()
       ini_pos = ini_pos.split(",")
       ini_neg = f1.readline().strip()
       ini_neg = ini_neg.split(",")

    # Run phuego.
    if(run_mock):
           test_path, test_df = load_test_example()
           print("Run phuego_mock with test dataset, whose first few lines are: \n",test_df.head())
           phuego_mock(
              support_data_folder=support_data_folder,
              res_folder=res_folder,
              test_path=test_path,
              fisher_geneset=fisher_geneset,
              fisher_threshold=fisher_threshold,
              fisher_background=fisher_background,
              ini_pos=ini_pos,
              ini_neg=ini_neg,
              damping=damping,
              kde_cutoff=kde_cutoff,
              use_existing_rwr=use_existing_rwr,
              convert2folder=convert2folder,
              include_isolated_egos_in_KDE_net=include_isolated_egos_in_kde_net,
              net_format=net_format,
              )
    elif(run_test):
           test_path, test_df = load_test_example()
           print("Run phuego with test dataset, whose first few lines are: \n",test_df.head())
           phuego(
              support_data_folder=support_data_folder,
              res_folder=res_folder,
              test_path=test_path,
              fisher_geneset=fisher_geneset,
              fisher_threshold=fisher_threshold,
              fisher_background=fisher_background,
              ini_pos=ini_pos,
              ini_neg=ini_neg,
              damping=damping,
              rwr_threshold=rwr_threshold,
              kde_cutoff=kde_cutoff,
              use_existing_rwr=use_existing_rwr,
              convert2folder=convert2folder,
              include_isolated_egos_in_KDE_net=include_isolated_egos_in_kde_net,
              net_format=net_format,
              )
    elif(test_path):
           print("Run phuego with user input dataset")
           phuego(
              support_data_folder=support_data_folder,
              res_folder=res_folder,
              test_path=test_path,
              fisher_geneset=fisher_geneset,
              fisher_threshold=fisher_threshold,
              fisher_background=fisher_background,
              ini_pos=ini_pos,
              ini_neg=ini_neg,
              damping=damping,
              rwr_threshold=rwr_threshold,
              kde_cutoff=kde_cutoff,
              use_existing_rwr=use_existing_rwr,
              convert2folder=convert2folder,
              include_isolated_egos_in_KDE_net=include_isolated_egos_in_kde_net,
              net_format=net_format,
              )
    else:
           print("Nothing is provided for running phuego.")