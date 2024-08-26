************
Using PICASO
************

.. _installation:

Installation
============


At the moment users have to add the path to the PICASO-repository to their python environment. This can be achieved by extending the python path:

.. code-block:: python

    import sys
    sys.path.insert(0, <path-to-PICASO-folder>)

After doing so you can use the PICASO framework by importing all functions from the kgraph-module:

.. code-block:: python

    from PICASO.kgraph import *

Please check out our :doc:`examples` on how to use PICASO.