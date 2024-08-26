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

- Creating the base :ref:`PICASO network </PICASO/scripts/create_basic_knowledgegraph.ipynb>`

More examples to come.
