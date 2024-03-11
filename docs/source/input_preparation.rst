Input preparation
=================

The main input of phuEGO is a text file where first column are protein uniprot
IDs and second columns are :ref:`aggregated <aggregation>` LFC or other scores
The two columns should be separated by tab and with no header.
The score dictate the probability by which the random walk process will restart
from the associated node, subject to influence by the scores of other nodes
within the same layer. User can :ref:`view the example input data <example_input>`.

.. _aggregation:

Aggregating phosphosites
~~~~~~~~~~~~~~~~~~~~~~~~

A phosphoproteomics experiment typically produced a list of phosphosites with 
log2-fold change between test and control samples. To perform network propagation
using protein-protein interaction network, the phosphosites need to be aggregated
into a protein list. In the :ref:`publication <citation>`, this is done with the 
following function:


.. container::

   The package provides the following Python function to perform this step:

      .. code-block:: python

         pass


Phosphosite functional score and SELPHI2.0 sign prediction (**Beta**)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

How phosphosrylation regulate protein function is a topic being actively studied.
Signs of log-fold change of phosphosites does not always correponds to signs of
regulation. On the other hand, some phosphorylation sites play a more important
role than others. To take these into consideration, the users can make use of 
`functional score predictions <url_to_Ochoa_paper>`__ and 
`signs predictions <url_to_SELPHI2.0_paper>`__ for phosphosite integration.

.. container::

   The package provides the following Python function to perform this step:

      .. code-block:: python

         pass