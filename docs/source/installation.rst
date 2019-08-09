.. _installation:

Installation
============

There are several ways to obtain SplAdder. You can directly install it from the Python package index
using `pip` or you can clone the git repository and set it up yourself.

Install from PyPi
^^^^^^^^^^^^^^^^^

The installation from Pypi is very straightforward. You can install the latest version of SplAdder
using `pip`::
    
    pip install spladder

If you would like to get a specific version, you can do::

    pip install spladder=2.2.0

Install from source
^^^^^^^^^^^^^^^^^^^

If you would like to get the latest changes available on GitHub, you can clone the repository::

    git clone https://github.com/ratschlab/spladder.git

This will create a directory named `spladder` inside your current directory. You can then install
SplAdder locally by first changing into the `spladder` directory and then typing::

    python setup.py install

If you are interested in using a specific branch, e.g., `development`, you would need to change into
that branch first, before you install. You can achieve this with::
    
    git checkout development
    python setup.py install
