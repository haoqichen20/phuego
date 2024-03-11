Welcome to phuego's documentation!
==================================

**phuEGO** is a network-based method to reconstruct active signalling
pathways from phosphoproteomics datasets. It combines three-layer
network propagation with ego network decomposition to provide small
networks comprising active functional signalling modules. phuEGO boosts
the signal-to-noise ratio from global phosphoproteomics datasets,
enriches the resulting networks for functional phosphosites and allows
the improved comparison and integration across datasets.

The Python package offers a command line interface (CLI) that implement
the method. Please follow the :ref:`documentation` below **step-by-step** for
setting up and running your analysis.


.. _documentation:

Documentation
-------------

.. toctree::
   :maxdepth: 2

   installation
   support_data
   mock_test_run
   input_preparation
   usage_cli
   output
   usage_api


Changelog
---------

For a detailed history of changes, see the
`Changelog <https://github.com/haoqichen20/phuego/blob/master/CHANGELOG.md>`__.


.. _citation:

Citation
--------

Please cite phuEGO if you use it in your analysis.

   **phuEGO: A network-based method to reconstruct active signalling
   pathways from phosphoproteomics datasets** *Girolamo Giudice, Haoqi
   Chen, Evangelia Petsalaki* bioRxiv, 2023 doi:
   `10.1101/2023.08.07.552249 <https://doi.org/10.1101/2023.08.07.552249>`__

.. code:: bibtex

   @article{giudice_phuego_2023,
       title = {{phuEGO}: {A} network-based method to reconstruct active signalling pathways from phosphoproteomics datasets},
       url = {https://www.biorxiv.org/content/early/2023/08/07/2023.08.07.552249},
       doi = {10.1101/2023.08.07.552249},
       abstract = {Signalling networks are critical for virtually all cell functions. Our current knowledge of cell signalling has been summarised in signalling pathway databases, which, while useful, are highly biassed towards well-studied processes, and don’t capture context specific network wiring or pathway cross-talk. Mass spectrometry-based phosphoproteomics data can provide a more unbiased view of active cell signalling processes in a given context, however, it suffers from low signal-to-noise ratio and poor reproducibility across experiments. Methods to extract active signalling signatures from such data struggle to produce unbiased and interpretable networks that can be used for hypothesis generation and designing downstream experiments. Here we present phuEGO, which combines three-layer network propagation with ego network decomposition to provide small networks comprising active functional signalling modules. PhuEGO boosts the signal-to-noise ratio from global phosphoproteomics datasets, enriches the resulting networks for functional phosphosites and allows the improved comparison and integration across datasets. We applied phuEGO to five phosphoproteomics data sets from cell lines collected upon infection with SARS CoV2. PhuEGO was better able to identify common active functions across datasets and to point to a subnetwork enriched for known COVID-19 targets. Overall, phuEGO provides a tool to the community for the improved functional interpretation of global phosphoproteomics datasets.Competing Interest StatementThe authors have declared no competing interest.},
       journal = {bioRxiv},
       author = {Giudice, Girolamo and Chen, Haoqi and Petsalaki, Evangelia},
       year = {2023},
       note = {Publisher: Cold Spring Harbor Laboratory
   \_eprint: https://www.biorxiv.org/content/early/2023/08/07/2023.08.07.552249.full.pdf},
   }


Contributors
------------

The algorithm and initial scripts of phuEGO are developed by Girolamo
Giudice (`@girolamogiudice <https://github.com/girolamogiudice>`__) and
Evangelia Petsalaki at `EMBL-EBI <https://www.ebi.ac.uk/>`__.

The Python package and CLI application are developed by Haoqi Chen
(`@haoqichen20 <https://github.com/haoqichen20>`__) at
`EMBL-EBI <https://www.ebi.ac.uk/>`__.