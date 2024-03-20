# -*- coding: utf-8 -*-

from .utils_CLI import load_test_example
from .utils_CLI import PythonLiteralOption
from .utils_CLI import GroupedOptions
from .phuego import phuego
from .phuego_mock import phuego_mock
import click

VERSION = '1.2.0'

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
    PATHS=['support_folder','result_folder'],
    SEED_PROPAGATION=['damping_seed','rwr_threshold'],
    EGO_DECOMPOSITION_CLUSTERING=['damping_ego','kde_cutoff', 'damping_module'],
    FISHER_TEST=['fisher_geneset', 'fisher_threshold', 'fisher_background'],
    OUTPUT=['net_format', 'convert2folder']
)
# Input/Output.
@click.option("--support_folder", "-sf", type=str, 
              help="Folder that stores support data.")
@click.option("--result_folder", "-rf", type=str, 
              help="Folder that stores phuEGO result.")
# Seed propagation.
@click.option("--damping_seed", "-ds", default=0.85, type=float, 
              help=("Damping factor of pagerank algorithm for seeds propagation,"
                    " Float number within range [0.5, 0.95]"
                    " [Default: 0.85]"))
@click.option("--rwr_threshold", "-rt", default=0.05, type=float, 
              help=("Threshold of significance for propagated nodes."
                    " Float number within range [0.01, 0.1]"
                    " [Default: 0.05]"))
# Ego decomposition and clustering.
@click.option("--damping_ego", "-de", default=0.85, type=float,
              help=("Damping factor of pagerank algorithm for ego decomposition,"
                    " Float number within range [0.5, 0.95]"
                    " [Default: 0.85]"))
@click.option("--kde_cutoff", "-k", cls=PythonLiteralOption, default="[0.85]", 
               help=("KDE cutoff value for removing nodes in ego network that are less similar to seed."
                     " Provided as a double quoted list of float number(s). Example: -k \"[0.5, 0.75]\"."
                     " Value should be within range [0.5,0.95]"
                     " [Default: [0.85]]"))
@click.option("--damping_module", "-dm", default=0.85, type=float,
              help=("Damping factor of pagerank algorithm for module detection,"
                    " Float number within page [0.5, 0.95]"
                    " [Default: 0.85]"))
# Fisher test.
@click.option("--fisher_geneset", "-fg", cls=PythonLiteralOption, default="['K']", 
               help=("Genesets to be tested when annotating modules of propagation."
                     " Provided as a quoted list of string(s). Example: -fg \"['C','F']\"."
                     " Value should be a subset of 'C,F,D,P,R,K,RT,B'."
                     " Refer to documentation to learn about the geneset abbreviation."
                     " [Default: \"['K']\"]"))
@click.option("--fisher_threshold", "-ft", default=0.05, type=float,
              help=("Threshold of significance of geneset enrichment analysis"
                    " [Default: 0.05]")) 
@click.option("--fisher_background", "-fb", default="intact", type=str,
              help=("Geneset background of geneset enrichment analysis"
                    " [Default: 'intact']"))
# Output configuration.
@click.option("--net_format", "-nf", default="graphml", type=str, 
              help=("file format of output network. Can be one of"
                    " 'edgelist', 'pajek', 'ncol', 'lgl', 'graphml', 'dimacs', 'gml', 'dot', 'leda'"
                    " [Default: 'graphml']"))
@click.option("--convert2folder/--dont_convert2folder", default=True, 
              help=("Should phuego organize result files into a folder structure?"
                    " [Default: --convert2folder]"))
def mock(support_folder, result_folder, 
    damping_seed, rwr_threshold, 
    damping_ego, kde_cutoff, damping_module, 
    fisher_geneset, fisher_threshold, fisher_background, 
    net_format, convert2folder):
    """
    Perform a mock run with test dataset. p value of seed propagation is estimated against 10 randomisation thus unreliable.
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
        support_data_folder=support_folder,
        res_folder=result_folder,
        test_path=test_path,
        fisher_geneset=fisher_geneset,
        fisher_threshold=fisher_threshold,
        fisher_background=fisher_background,
        ini_pos=ini_pos,
        ini_neg=ini_neg,
        damping_seed_propagation=damping_seed,
        damping_ego_decomposition=damping_ego,
        damping_module_detection=damping_module,
        kde_cutoff=kde_cutoff,
        use_existing_rwr=use_existing_rwr,
        convert2folder=convert2folder,
        include_isolated_egos_in_KDE_net=include_isolated_egos_in_kde_net,
        net_format=net_format,
        )



@register_command
@click.command(cls=GroupedOptions)
@grouped_options(
    PATHS=['support_folder','result_folder'],
    SEED_PROPAGATION=['use_existing_rwr','damping_seed','rwr_threshold'],
    EGO_DECOMPOSITION_CLUSTERING=['damping_ego','kde_cutoff', 'damping_module'],
    FISHER_TEST=['fisher_geneset', 'fisher_threshold', 'fisher_background'],
    OUTPUT=['net_format', 'convert2folder']
)
# Input/Output.
@click.option("--support_folder", "-sf", type=str, 
              help="Folder that stores support data.")
@click.option("--result_folder", "-rf", type=str, 
              help="Folder that stores phuEGO result.")
# Seed propagation.
@click.option("--use_existing_rwr", "-ru", is_flag=True, 
              help=("Should phuego reuse existing network propagation results?"
                    " [Flag option. If unused, default is False]"))
@click.option("--damping_seed", "-ds", default=0.85, type=float, 
              help=("Damping factor of pagerank algorithm for seeds propagation,"
                    " Float number within range [0.5, 0.95]"
                    " [Default: 0.85]"))
@click.option("--rwr_threshold", "-rt", default=0.05, type=float, 
              help=("Threshold of significance for propagated nodes."
                    " Float number within range [0.01, 0.1]"
                    " [Default: 0.05]"))
# Ego decomposition and clustering.
@click.option("--damping_ego", "-de", default=0.85, type=float,
              help=("Damping factor of pagerank algorithm for ego decomposition,"
                    " Float number within range [0.5, 0.95]"
                    " [Default: 0.85]"))
@click.option("--kde_cutoff", "-k", cls=PythonLiteralOption, default="[0.85]", 
               help=("KDE cutoff value for removing nodes in ego network that are less similar to seed."
                     " Provided as a double quoted list of float number(s). Example: -k \"[0.5, 0.75]\"."
                     " Value should be within range [0.5,0.95]"
                     " [Default: [0.85]]"))
@click.option("--damping_module", "-dm", default=0.85, type=float,
              help=("Damping factor of pagerank algorithm for module detection,"
                    " Float number within page [0.5, 0.95]"
                    " [Default: 0.85]"))
# Fisher test.
@click.option("--fisher_geneset", "-fg", cls=PythonLiteralOption, default="['K']", 
               help=("Genesets to be tested when annotating modules of propagation."
                     " Provided as a quoted list of string(s). Example: -fg \"['C','F']\"."
                     " Value should be a subset of 'C,F,D,P,R,K,RT,B'."
                     " Refer to documentation to learn about the geneset abbreviation."
                     " [Default: \"['K']\"]"))
@click.option("--fisher_threshold", "-ft", default=0.05, type=float,
              help=("Threshold of significance of geneset enrichment analysis"
                    " [Default: 0.05]")) 
@click.option("--fisher_background", "-fb", default="intact", type=str,
              help=("Geneset background of geneset enrichment analysis"
                    " [Default: 'intact']"))
# Output configuration.
@click.option("--net_format", "-nf", default="graphml", type=str, 
              help=("file format of output network. Can be one of"
                    " 'edgelist', 'pajek', 'ncol', 'lgl', 'graphml', 'dimacs', 'gml', 'dot', 'leda'"
                    " [Default: 'graphml']"))
@click.option("--convert2folder/--dont_convert2folder", default=True, 
              help=("Should phuego organize result files into a folder structure?"
                    " [Default: --convert2folder]"))
def test(support_folder, result_folder, 
    damping_seed, rwr_threshold, use_existing_rwr,
    damping_ego, kde_cutoff, damping_module, 
    fisher_geneset, fisher_threshold, fisher_background, 
    net_format, convert2folder):
    """
    Perform a full run with test dataset.
    """
    test_path, test_df = load_test_example()
    print("Run phuego with test dataset, whose first few lines are: \n",test_df.head())

    # Test run does not allow removing nodes from reference network.
    ini_pos = ["False"]
    ini_neg = ["False"]
    # Test run does not allow excluding isolated egos.
    include_isolated_egos_in_kde_net = True
    # Test run does not allow defining layers.
    layer_division = "phos"
    layer_def_path = ""
    # Test run does not allow defining other semantic similarity.
    semsim = "gic"
    
    phuego(
        support_data_folder=support_folder,
        res_folder=result_folder,
        test_path=test_path,
        fisher_geneset=fisher_geneset,
        fisher_threshold=fisher_threshold,
        fisher_background=fisher_background,
        layer_division=layer_division,
        layer_definition_path=layer_def_path,
        ini_pos=ini_pos,
        ini_neg=ini_neg,
        damping_seed_propagation=damping_seed,
        damping_ego_decomposition=damping_ego,
        damping_module_detection=damping_module,
        rwr_threshold=rwr_threshold,
        semsim=semsim,
        kde_cutoff=kde_cutoff,
        use_existing_rwr=use_existing_rwr,
        convert2folder=convert2folder,
        include_isolated_egos_in_KDE_net=include_isolated_egos_in_kde_net,
        net_format=net_format,
        )



@register_command
@click.command(cls=GroupedOptions)
@grouped_options(
    PATHS=['support_folder','result_folder','test_path'],
    SEED_PROPAGATION=['layer_division','layer_def_path','ini_path','use_existing_rwr','damping_seed','rwr_threshold','semsim'],
    EGO_DECOMPOSITION_CLUSTERING=['damping_ego','kde_cutoff', 'damping_module'],
    FISHER_TEST=['fisher_geneset', 'fisher_threshold', 'fisher_background'],
    OUTPUT=['net_format', 'convert2folder', 'include_isolated_egos_in_kde_net']
)
# Input/Output.
@click.option("--support_folder", "-sf", type=str, 
              help="Folder that stores support data.")
@click.option("--result_folder", "-rf", type=str, 
              help="Folder that stores phuEGO result.")
@click.option("--test_path", "-tpath", type=str, 
              help="Path to user test file")
# Seed propagation.
@click.option("--layer_division","-ld", default="phos", type=str,
              help=("""How to divide the seed nodes into layers.
                     \b\n'phos': divide seed nodes into tyrosine kinases, ser/thr kinases and others;
                     \b\n'custom': use -ldpath to provide a text file with user-defined layer division;
                     \b\n'one': analyze the input data without layer division. [Default: \"phos\"]"""))
@click.option("--layer_def_path", "-ldpath", default="", type=str,
              help=("A text file with user-defined layer division."
                    " Required if -ld == \"custom\"."
                    " Refer to documentation to learn about formatting."
                    " [Default: \"\"]"))
@click.option("--ini_path", "-ipath", default="", type=str, 
              help=("The path of a csv file, that specify two sets of nodes to be"
                    " removed from reference network during propagation."
                    " [Default: \"\"]"))
@click.option("--use_existing_rwr", "-ru", is_flag=True, 
              help=("Should phuego reuse existing seeds network propagation results?"
                    " [Flag option. If unused, default is False]"))
@click.option("--damping_seed", "-ds", default=0.85, type=float, 
              help=("Damping factor of pagerank algorithm for seeds propagation,"
                    " Float number within range [0.5, 0.95]"
                    " [Default: 0.85]"))
@click.option("--rwr_threshold", "-rt", default=0.05, type=float, 
              help=("Threshold of significance for propagated nodes."
                    " Float number within range [0.01, 0.1]"
                    " [Default: 0.05]"))
@click.option("--semsim", "-ss", default="gic", type=str, 
              help=("The type of semantic similarity to use. If using slim version"
                    " of premade supporting dataset, only \"gic\" is available."
                    " [Default: \"gic\"]"))
# Ego decomposition and clustering.
@click.option("--damping_ego", "-de", default=0.85, type=float,
              help=("Damping factor of pagerank algorithm for ego decomposition,"
                    " Float number within range [0.5, 0.95]"
                    " [Default: 0.85]"))
@click.option("--kde_cutoff", "-k", cls=PythonLiteralOption, default="[0.85]", 
               help=("KDE cutoff value for removing nodes in ego network that are less similar to seed."
                     " Provided as a double quoted list of float number(s). Example: -k \"[0.5, 0.75]\"."
                     " Value should be within range [0.5,0.95]"
                     " [Default: [0.85]]"))
@click.option("--damping_module", "-dm", default=0.85, type=float,
              help=("Damping factor of pagerank algorithm for module detection,"
                    " Float number within page [0.5, 0.95]"
                    " [Default: 0.85]"))
# Fisher test.
@click.option("--fisher_geneset", "-fg", cls=PythonLiteralOption, default="['K']", 
               help=("Genesets to be tested when annotating modules of propagation."
                     " Provided as a quoted list of string(s). Example: -fg \"['C','F']\"."
                     " Value should be a subset of 'C,F,D,P,R,K,RT,B'."
                     " Refer to documentation to learn about the geneset abbreviation."
                     " [Default: \"['K']\"]"))
@click.option("--fisher_threshold", "-ft", default=0.05, type=float,
              help=("Threshold of significance of geneset enrichment analysis"
                    " [Default: 0.05]")) 
@click.option("--fisher_background", "-fb", default="intact", type=str,
              help=("Geneset background of geneset enrichment analysis"
                    " [Default: 'intact']"))
# Output configuration.
@click.option("--net_format", "-nf", default="graphml", type=str, 
              help=("file format of output network. Can be one of"
                    " 'edgelist', 'pajek', 'ncol', 'lgl', 'graphml', 'dimacs', 'gml', 'dot', 'leda'"
                    " [Default: 'graphml']"))
@click.option("--convert2folder/--dont_convert2folder", default=True, 
              help=("Should phuego organize result files into a folder structure?"
                    " [Default: --convert2folder]"))
@click.option("--include_isolated_egos_in_kde_net", "-ie", is_flag=True, 
              help=("Should we include isolated nodes in the output network?"
                    " [Flag option. If unused, default is False]")) 
def main(support_folder, result_folder, test_path, 
    layer_division, layer_def_path, ini_path, damping_seed, 
    rwr_threshold, use_existing_rwr, semsim,
    damping_ego, kde_cutoff, damping_module, 
    fisher_geneset, fisher_threshold, fisher_background, 
    net_format, convert2folder, include_isolated_egos_in_kde_net):
    """
    Run phuEGO with user input.
    """
    print("Run phuEGO with user input dataset")

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

    # Assert if the layer_division and the layer_def_path is provided properly. 
    # Including checkin g if the file exist (can be accessed).


    phuego(
        support_data_folder=support_folder,
        res_folder=result_folder,
        test_path=test_path,
        fisher_geneset=fisher_geneset,
        fisher_threshold=fisher_threshold,
        fisher_background=fisher_background,
        layer_division=layer_division,
        layer_definition_path=layer_def_path,
        ini_pos=ini_pos,
        ini_neg=ini_neg,
        damping_seed_propagation=damping_seed,
        damping_ego_decomposition=damping_ego,
        damping_module_detection=damping_module,
        rwr_threshold=rwr_threshold,
        semsim=semsim,
        kde_cutoff=kde_cutoff,
        use_existing_rwr=use_existing_rwr,
        convert2folder=convert2folder,
        include_isolated_egos_in_KDE_net=include_isolated_egos_in_kde_net,
        net_format=net_format,
        )



if __name__ == '__main__':
    cli()