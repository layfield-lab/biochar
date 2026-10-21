Installation
============

Requirements
------------

- Python 3.9 or later
- `RDKit <https://www.rdkit.org>`_ ≥ 2023.9
- NumPy ≥ 1.24
- SciPy ≥ 1.10
- NetworkX ≥ 3.1

conda (recommended)
-------------------

``biochar`` is available on **conda-forge** and installs all dependencies,
including RDKit, in a single command:

.. code-block:: bash

   conda install -c conda-forge biochar

To install into a fresh environment:

.. code-block:: bash

   conda create -n biochar -c conda-forge biochar
   conda activate biochar

From PyPI
---------

.. code-block:: bash

   pip install biochar

.. note::
   On some platforms the ``rdkit`` wheel from PyPI may not be available.
   If the install fails, use the conda route above.

Development install
-------------------

Clone the repository and install in editable mode with all dev dependencies:

.. code-block:: bash

   git clone https://github.com/layfield-lab/biochar.git
   cd biochar
   pip install -e ".[dev]"

Verify the installation by running the test suite:

.. code-block:: bash

   pytest tests/ -v
