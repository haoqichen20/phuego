# -*- coding: utf-8 -*-

from .utils_CLI import load_test_example
from .utils_CLI import PythonLiteralOption
from .phuego import phuego
from .phuego_mock import phuego_mock
import click


__version__ = "1.1.0"


'''
This is the CLI tool for phuego.
'''

# changes the default parameters to -h and --help instead of just --help
CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])

@click.command(context_settings=CONTEXT_SETTINGS)
# Folders.
@click.option("--support_data_folder", "-sf", type=str, 
              help="Folder that stored support data. Directories separated by /")
@click.option("--res_folder", "-rf", type=str, 
              help="Folder that stored output result. Directories separated by /")
@click.option("--test_path", "-tpath", default="", type=str, 
              help="The path of user input test file. If defined, -mock and -test should not be used.")
@click.option("--ini_path", "-ipath", default="", type=str, 
              help="The path of a csv file, that specify two sets of nodes to be"
                   " removed from reference network during propagation.")
@click.option("--convert2folder/--dont_convert2folder", default=True, 
              help="Should phuego organize result files into a folder structure?")
# Mock/test.
@click.option("--run_mock", is_flag=True, 
              help=("Should phuego perform a mock test run using test data?"
                    " If defined, -test and -tpath would be ignored."))
@click.option("--run_test", is_flag=True, 
              help=("Should phuego perform a full test run using test data?"
                    " If defined, -tpath would be ignored."))
# network propagation parameters.
@click.option("--use_existing_rwr", "-ru", is_flag=True, 
              help="Should phuego reuse existing network propagation results?")
@click.option("--damping_seed_propagation", "-ds", default=0.85, type=float, 
              help=("Damping factor of pagerank algorithm for seeds propagation,"
                    " Float number within range [0.5, 0.95]"))
@click.option("--rwr_threshold", "-rt", default=0.05, type=float, 
              help=("Threshold of significance for propagated nodes."
                    " Float number within range [0.01, 0.1]"))
# ego propagation and fisher test.
@click.option("--damping_ego_decomposition", "-de", default=0.85, type=float,
              help=("Damping factor of pagerank algorithm for ego decomposition,"
                    " Float number within range [0.5, 0.95]"))
@click.option("--damping_module_detection", "-dm", default=0.85, type=float,
              help=("Damping factor of pagerank algorithm for module detection,"
                    " Float number within page [0.5, 0.95]"))
@click.option("--kde_cutoff", "-k", cls=PythonLiteralOption, default="[0.85]", 
               help=("KDE cutoff value for removing nodes in ego network that are less similar to seed."
                     " Provided as a double quoted list of float number(s). Example: -k \"[0.5, 0.75]\"."
                     " Value should be within range [0.5,0.95]"))
@click.option("--fisher_geneset", "-fg", cls=PythonLiteralOption, default="['K']", 
               help=("Genesets to be tested when annotating modules of propagation."
                     " Provided as a quoted list of string(s). Example: -fg \"['C','F']\"."
                     " Value should be a subset of 'C,F,D,P,R,K,RT,B'."
                     " Refer to documentation to learn about the geneset abbreviation."))
@click.option("--fisher_threshold", "-ft", default=0.05, type=float,
              help="Threshold of significance of geneset enrichment analysis") 
@click.option("--fisher_background", "-fb", default="intact", type=str,
              help="Geneset background of geneset enrichment analysis") 
# network output
@click.option("--include_isolated_egos_in_kde_net", "-ie", is_flag=True, 
              help="Should we include isolated nodes in the output network?") 
@click.option("--net_format", "-nf", default="graphml", type=str, 
              help=("file format of output network. Can be one of"
                    " 'edgelist', 'pajek', 'ncol', 'lgl', 'graphml', 'dimacs', 'gml', 'dot', 'leda'"))
# Versioning.
@click.option('--version', '-v', is_flag=True, 
              help='Print version to stdout')
def main(support_data_folder, res_folder, test_path, convert2folder, use_existing_rwr, run_test, run_mock, 
         damping_seed_propagation, damping_ego_decomposition, damping_module_detection, 
         fisher_geneset, fisher_threshold, fisher_background, kde_cutoff, ini_path, rwr_threshold, 
         include_isolated_egos_in_kde_net, net_format, version) -> None:
    """
    phuEGO documentation: https://phuego.readthedocs.io/en/latest/index.html
    """

    """
    Print the version number.
    """  
    if version:
        bold_package_name = click.style("phuego", bold=True)
        click.echo(f"{bold_package_name} ({__version__})")
        return


    """
    User input formatting and assertion.
    """
    # Assert that kde_cutoff is a list and contains only floats.
    assert (
       isinstance(kde_cutoff, list) and all(isinstance(item, float) for item in kde_cutoff)
    ), (
       "kde_cutoff must be a list of floats.")
 
    # Assert that fisher_geneset is a list and contains only the allowed geneset label.
    allowed_labels = {"C", "F", "D", "P", "R", "K", "RT", "B"}
    assert (
       isinstance(fisher_geneset, list) and all(item in allowed_labels for item in fisher_geneset)
    ), (
       "Variable must be a list containing only specific strings: C, F, D, P, R, K, RT, B"
       )

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


    """
    Run phuego.
    """
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
              damping_seed_propagation=damping_seed_propagation,
              damping_ego_decomposition=damping_ego_decomposition,
              damping_module_detection=damping_module_detection,
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
              damping_seed_propagation=damping_seed_propagation,
              damping_ego_decomposition=damping_ego_decomposition,
              damping_module_detection=damping_module_detection,
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
              damping_seed_propagation=damping_seed_propagation,
              damping_ego_decomposition=damping_ego_decomposition,
              damping_module_detection=damping_module_detection,
              rwr_threshold=rwr_threshold,
              kde_cutoff=kde_cutoff,
              use_existing_rwr=use_existing_rwr,
              convert2folder=convert2folder,
              include_isolated_egos_in_KDE_net=include_isolated_egos_in_kde_net,
              net_format=net_format,
              )
    else:
           print("Nothing is provided for running phuego.")