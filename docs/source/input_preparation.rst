Input preparation
=================

The main input of phuEGO is a tab-delimited text file with two columns of no headers.
The first column contains protein uniprot IDs and the second column contains
:ref:`aggregated <aggregation>` log2-fold change (LFC) or other scores.
The score dictates the probability by which the random walk process will restart
from the associated node, subject to influence by the scores of other nodes
within the same layer. User can :ref:`view the example input data <example_input>`.

.. _aggregation:

Aggregating phosphosites (**Active development**)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A phosphoproteomics experiment typically produced a list of phosphosites with 
LFC between test and control samples. To perform network propagation
using protein-protein interaction network, the phosphosites need to be aggregated
into a protein list. In the :ref:`publication <citation>`, this is performed by
filtering phosphosites on `functional scores <https://www.nature.com/articles/s41587-019-0344-3>`__, 
and aggregating using maximal LFC. 

The authors of phuEGO are currently developing additional functions for this step, which will be released
in near-future update of the package. Besides functional score, 
predicted signs of regulation of the phosphosites on protein functions will also be used for aggregation.

The users can also develop their own method of phosphosite aggregation.

.. container::

   The package provides the following Python function to perform this step:

      .. code-block:: python

         pass