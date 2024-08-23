# PICASO

## Abstract

Various single-cell modalities covering transcriptomics, epigenetic and spatio-temporal changes in health and disease phenotypes are used in an exploratory way to understand biological systems at single-cell resolution. However, the vast amount of such single-cell data is not systematically linked to existing biomedical data. Networks have previously been used to represent harmonized biomedical data. Integrating various resources of biomedical data in networks  has recently received increasing attention. These aggregated networks can provide additional insight into the biology of complex human diseases at cell-type level, however, lack inclusion of single cell expression data. Here, we present the PICASO framework, which incorporates single-cell gene expression data as an additional layer to represent associations between cell types, disease phenotypes, drugs and genes. The PICASO network includes several standardized biomedical databases such as STRING, Uniprot, GeneOntology, Reactome, OmniPath and OpenTargets. Using multiple cell type-specific instances of the framework, each annotated and scored with their respective expression data, comparisons between disease states can be made by computing respective sub-networks and comparing the expression scores between conditions. Ultimately, these group-specific networks will allow the identification of relevant genes, processes and potentially druggable targets, as well as the comparison of different measured groups and thus the identification of group-specific communities and interactions.

## Collaborators

[Collaboration with @HayatLab](https://github.com/hayatlab)

# Setup

## Creating environment

    conda create --name regnetworks -c conda-forge python=3.11
    conda activate regnetworks

    # To run in jupyter
    conda install -c conda-forge python-igraph ipykernel
    python -m ipykernel install --user --name=regnetworks

    pip install scanpy networkx scikit-learn natsort seaborn matplotlib leidenalg markov_clustering

    # For the AIDescriptor class
    pip install huggingface_hub llama-cpp-python 'transformers[tf-cpu]'

## Using PICASO

At the moment users have to add the path to the PICASO-repository to their python environment. This can be achieved by extending the python path:

    import sys
    sys.path.insert(0, <path-to-PICASO-folder>)

After doing so you can use the PICASO framework by importing all functions from the kgraph-module:

    from PICASO.kgraph import *

Please also consider our notebooks for usage examples of the PICASO framework:

- Creating the base PICASO network: [notebook](./blob/main/scripts/create_basic_knowledgegraph.ipynb)

- KPMP snRNA-seq data preparation: [notebook](./blob/main/kpmp/process_snrna.ipynb)

- KPMP PICASO data preparation: [notebook](./blob/main/kpmp/kpmp_celltype_zone_prepare.ipynb)

- KPMP PICASO analysis: [notebook](./blob/main/scripts/kpmp_celltype_zone_diff_analysis_objectified.ipynb)


- Myocardial Infarction PICASO analysis: [notebook](./blob/main/scripts/mi_celltype_zone_diff_analysis.ipynb)



