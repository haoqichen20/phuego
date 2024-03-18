[![Read the Docs](https://img.shields.io/readthedocs/phuego)](https://phuego.readthedocs.io/en/latest/index.html)
[![PyPI - Version](https://img.shields.io/pypi/v/phuego)](https://pypi.org/project/phuego/)
[![PyPI - Python Version](https://img.shields.io/pypi/pyversions/phuego)](https://pypi.org/project/phuego/)
[![PyPI - Wheel](https://img.shields.io/pypi/wheel/phuego)](https://pypi.org/project/phuego/)

# phuEGO: a network-based method to reconstruct active signalling pathways

phuEGO is a network-based method to reconstruct active signalling pathways from phosphoproteomics datasets. It combines three-layer network propagation with ego network decomposition to provide small networks comprising active functional signalling modules. PhuEGO boosts the signal-to-noise ratio from global phosphoproteomics datasets, enriches the resulting networks for functional phosphosites and allows the improved comparison and integration across datasets.

Refer to [documentation](https://phuego.readthedocs.io/en/latest/) for details on installation and executation.
The package is currently distributed with [PyPI](https://pypi.org/project/phuego/). 
Check the [Zenodo](https://zenodo.org/records/8094690) page for the premade support dataset, which contains processed reference network with semantic similarity,
or the [Nextflow pipeline](https://github.com/haoqichen20/phuego_support) to generate your customized support dataset.


## Citation

Please cite phuEGO if you use it in your analysis.

> **phuEGO: A network-based method to reconstruct active signalling pathways from phosphoproteomics datasets** <br> _Girolamo Giudice, Haoqi Chen, Evangelia Petsalaki_ <br>
> bioRxiv, 2023 <br>
> doi: [10.1101/2023.08.07.552249](https://doi.org/10.1101/2023.08.07.552249) <br>


## Contributors

The algorithm and initial scripts of phuEGO are developed by Girolamo Giudice ([@girolamogiudice](https://github.com/girolamogiudice)) and Evangelia Petsalaki at [EMBL-EBI](https://www.ebi.ac.uk/).

The Python package and command line interface are developed by Haoqi Chen ([@haoqichen20](https://github.com/haoqichen20)) at [EMBL-EBI](https://www.ebi.ac.uk/).
