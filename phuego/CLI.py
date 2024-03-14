# -*- coding: utf-8 -*-

from .utils_CLI import load_test_example
from .utils_CLI import PythonLiteralOption
from .utils_CLI import GroupedOptions
from .phuego import phuego
from .phuego_mock import phuego_mock
import click

VERSION = '1.1.0'

# changes the default parameters to -h and --help instead of just --help
CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])

@click.group(context_settings=CONTEXT_SETTINGS)
@click.version_option(version=VERSION, prog_name='phuEGO', message='%(prog)s version %(version)s')
def cli():
    """phuEGO documentation: https://phuego.readthedocs.io/en/latest/index.html"""
    pass

def register_command(command_func):
    """Decorator to register commands with the CLI."""
    cli.add_command(command_func)
    return command_func

def grouped_options(**groups):
    """Decorator to attach grouped options metadata to a command function."""
    def wrapper(func):
        func.grouped_options = groups
        return func
    return wrapper



@register_command
@click.command(cls=GroupedOptions)
@grouped_options(
    PATHS=['support_data_folder','res_folder'],
    SEED_PROPAGATION=['damping_seed_propagation','rwr_threshold'],
    EGO_DECOMPOSITION_CLUSTERING=['damping_ego_decomposition','kde_cutoff', 'damping_module_detection'],
    FISHER_TEST=['fisher_geneset', 'fisher_threshold', 'fisher_background'],
    OUTPUT=['net_format', 'convert2folder']
)
# Paths.
@click.option("--support_data_folder", "-sf", type=str, 
              help="Folder that stored support data. Directories separated by /")
@click.option("--res_folder", "-rf", type=str, 
              help="Folder that stored output result. Directories separated by /")
# Seed propagation.
@click.option("--damping_seed_propagation", "-ds", default=0.85, type=float, 
              help=("Damping factor of pagerank algorithm for seeds propagation,"
                    " Float number within range [0.5, 0.95]"))
@click.option("--rwr_threshold", "-rt", default=0.05, type=float, 
              help=("Threshold of significance for propagated nodes."
                    " Float number within range [0.01, 0.1]"))
# Ego decomposition and clustering.
@click.option("--damping_ego_decomposition", "-de", default=0.85, type=float,
              help=("Damping factor of pagerank algorithm for ego decomposition,"
                    " Float number within range [0.5, 0.95]"))
@click.option("--kde_cutoff", "-k", cls=PythonLiteralOption, default="[0.85]", 
               help=("KDE cutoff value for removing nodes in ego network that are less similar to seed."
                     " Provided as a double quoted list of float number(s). Example: -k \"[0.5, 0.75]\"."
                     " Value should be within range [0.5,0.95]"))
@click.option("--damping_module_detection", "-dm", default=0.85, type=float,
              help=("Damping factor of pagerank algorithm for module detection,"
                    " Float number within page [0.5, 0.95]"))
# Fisher test.
@click.option("--fisher_geneset", "-fg", cls=PythonLiteralOption, default="['K']", 
               help=("Genesets to be tested when annotating modules of propagation."
                     " Provided as a quoted list of string(s). Example: -fg \"['C','F']\"."
                     " Value should be a subset of 'C,F,D,P,R,K,RT,B'."
                     " Refer to documentation to learn about the geneset abbreviation."))
@click.option("--fisher_threshold", "-ft", default=0.05, type=float,
              help="Threshold of significance of geneset enrichment analysis") 
@click.option("--fisher_background", "-fb", default="intact", type=str,
              help="Geneset background of geneset enrichment analysis") 
# Output configuration.
@click.option("--net_format", "-nf", default="graphml", type=str, 
              help=("file format of output network. Can be one of"
                    " 'edgelist', 'pajek', 'ncol', 'lgl', 'graphml', 'dimacs', 'gml', 'dot', 'leda'"))
@click.option("--convert2folder/--dont_convert2folder", default=True, 
              help="Should phuego organize result files into a folder structure?")
def run_mock(support_data_folder, res_folder, 
    damping_seed_propagation, rwr_threshold, 
    damping_ego_decomposition, kde_cutoff, damping_module_detection, 
    fisher_geneset, fisher_threshold, fisher_background, 
    net_format, convert2folder):
    """
    Perform a mock run with test dataset.
    """
    test_path, test_df = load_test_example()
    print("Run phuego_mock with test dataset, whose first few lines are: \n",test_df.head())

    # Mock run does not allow removing nodes from reference network.
    ini_pos = ["False"]
    ini_neg = ["False"]
    # Mock run does not allow reusing existing options.
    use_existing_rwr = False
    # Mock run does not allow excluding isolated egos.
    include_isolated_egos_in_kde_net = True

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


# @register_command
# @click.command(cls=GroupedOptions)
# @grouped_options(
#     Paths=[],
#     Seed_propagation=[],
#     Ego_decomposition=[],
#     Module_detection=[],
#     Output=[]
# )
# def run_test():
#     """
#     Perform a test run with test dataset.
#     """
#     test_path, test_df = load_test_example()
#     print("Run phuego with test dataset, whose first few lines are: \n",test_df.head())
#     phuego(
#         support_data_folder=support_data_folder,
#         res_folder=res_folder,
#         test_path=test_path,
#         fisher_geneset=fisher_geneset,
#         fisher_threshold=fisher_threshold,
#         fisher_background=fisher_background,
#         ini_pos=ini_pos,
#         ini_neg=ini_neg,
#         damping_seed_propagation=damping_seed_propagation,
#         damping_ego_decomposition=damping_ego_decomposition,
#         damping_module_detection=damping_module_detection,
#         rwr_threshold=rwr_threshold,
#         kde_cutoff=kde_cutoff,
#         use_existing_rwr=use_existing_rwr,
#         convert2folder=convert2folder,
#         include_isolated_egos_in_KDE_net=include_isolated_egos_in_kde_net,
#         net_format=net_format,
#         )

# @register_command
# @click.command(cls=GroupedOptions)
# @grouped_options(
#     Paths=[],
#     Seed_propagation=[],
#     Ego_decomposition=[],
#     Module_detection=[],
#     Output=[]
# )
# def run_phos():
#     """
#     Run phuEGO by dividing input into three layers: tyrosine kinase, ser/thr kinase, other proteins. 
#     Recommended for phosphoproteomics dataset.
#     """
#     print("Run phuego with user input dataset")
#     phuego(
#         support_data_folder=support_data_folder,
#         res_folder=res_folder,
#         test_path=test_path,
#         fisher_geneset=fisher_geneset,
#         fisher_threshold=fisher_threshold,
#         fisher_background=fisher_background,
#         ini_pos=ini_pos,
#         ini_neg=ini_neg,
#         damping_seed_propagation=damping_seed_propagation,
#         damping_ego_decomposition=damping_ego_decomposition,
#         damping_module_detection=damping_module_detection,
#         rwr_threshold=rwr_threshold,
#         kde_cutoff=kde_cutoff,
#         use_existing_rwr=use_existing_rwr,
#         convert2folder=convert2folder,
#         include_isolated_egos_in_KDE_net=include_isolated_egos_in_kde_net,
#         net_format=net_format,
#         )


# @register_command
# @click.command(cls=GroupedOptions)
# @grouped_options(
#     Paths=[],
#     Seed_propagation=[],
#     Ego_decomposition=[],
#     Module_detection=[],
#     Output=[]
# )
# def run_custom():
#     """
#     Run phuEGO by dividing input into user-defined layers (no more than three).
#     """


# @register_command
# @click.command(cls=GroupedOptions)
# @grouped_options(
#     Paths=[],
#     Seed_propagation=[],
#     Ego_decomposition=[],
#     Module_detection=[],
#     Output=[]
# )
# def run_sc():
#     """
#     Run phuEGO by dividing input into three layers: receptors, transcription factors, and others.
#     Recommended for single cell RNA-seq dataset.
#     """




if __name__ == '__main__':
    cli()