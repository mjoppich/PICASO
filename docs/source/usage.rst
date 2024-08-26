Usaging PICASO
==============

.. _installation:

Installation
------------

At the moment users have to add the path to the PICASO-repository to their python environment. This can be achieved by extending the python path:

.. code-block:: python
    import sys
    sys.path.insert(0, <path-to-PICASO-folder>)

After doing so you can use the PICASO framework by importing all functions from the kgraph-module:

.. code-block:: python
    from PICASO.kgraph import *

Examples
--------

Please also consider our notebooks for usage examples of the PICASO framework:

- Creating the base PICASO network: [notebook](./blob/main/scripts/create_basic_knowledgegraph.ipynb)

- KPMP snRNA-seq data preparation: [notebook](./blob/main/kpmp/process_snrna.ipynb)

- KPMP PICASO data preparation: [notebook](./blob/main/kpmp/kpmp_celltype_zone_prepare.ipynb)

- KPMP PICASO analysis: [notebook](./blob/main/scripts/kpmp_celltype_zone_diff_analysis_objectified.ipynb)

- Myocardial Infarction PICASO analysis: [notebook](./blob/main/scripts/mi_celltype_zone_diff_analysis.ipynb)
