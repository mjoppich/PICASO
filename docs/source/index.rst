**********************************
Welcome to PICASO's documentation!
**********************************

.. role:: raw-html(raw)
   :format: html



**PICASO** is a Python library for :raw-html:`<strong>P</strong>rofiling <strong>I</strong>ntegrative <strong>C</strong>ommunities of <strong>A</strong>ggregated <strong>S</strong>ingle-cell <strong>O</strong>mics data (<strong>PICASO</strong>)`.

Check out the :doc:`usage` section for further information, including
how to :ref:`installation` the project.

.. note::

   This project is under active development.


Abstract
========

Various single-cell modalities covering transcriptomics, epigenetic and spatio-temporal changes in health and disease phenotypes are used in an exploratory way to understand biological systems at single-cell resolution. However, the vast amount of such single-cell data is not systematically linked to existing biomedical data. Networks have previously been used to represent harmonized biomedical data. Integrating various resources of biomedical data in networks  has recently received increasing attention. These aggregated networks can provide additional insight into the biology of complex human diseases at cell-type level, however, lack inclusion of single cell expression data. Here, we present the PICASO framework, which incorporates single-cell gene expression data as an additional layer to represent associations between cell types, disease phenotypes, drugs and genes. The PICASO network includes several standardized biomedical databases such as STRING, Uniprot, GeneOntology, Reactome, OmniPath and OpenTargets. Using multiple cell type-specific instances of the framework, each annotated and scored with their respective expression data, comparisons between disease states can be made by computing respective sub-networks and comparing the expression scores between conditions. Ultimately, these group-specific networks will allow the identification of relevant genes, processes and potentially druggable targets, as well as the comparison of different measured groups and thus the identification of group-specific communities and interactions.

If you find this work helpful for your analysis, `please cite <https://www.biorxiv.org/content/10.1101/2024.08.28.610120v1>`_ :

`PICASO: Profiling Integrative Communities of Aggregated Single-cell Omics data<https://www.biorxiv.org/content/10.1101/2024.08.28.610120v1>`_
`Markus Joppich, Rafael Kramann, Sikander Hayat<https://www.biorxiv.org/content/10.1101/2024.08.28.610120v1>`_
`bioRxiv 2024.08.28.610120; doi: https://doi.org/10.1101/2024.08.28.610120 <https://www.biorxiv.org/content/10.1101/2024.08.28.610120v1>`_


.. image:: ../picaso_framework.png
    :align: center
    :alt: PICASO graphical abstract

Contents
--------

.. toctree::
    :maxdepth: 1
    :titlesonly:

    usage
    examples
    api
