Get started
===========

Installing Python
-----------------

The script has been built using Python 3.7.7. Since it doesn't rely on very extensive Python features I would expect it to run just fine on any version of Python3.
If you do not already have Python installed I recommend you `download Anaconda <https://www.anaconda.com>`_.


Installing dependencies
-----------------------

The script uses several dependencies, once Python is installed you can install these by running

.. code-block:: bash

    pip install --user --upgrade numpy pandas matplotlib statsmodels scipy


Downloading script
------------------

Download the latest version of the script from `GitHub <https://github.com/matthijsz/KinshipCorrelationGenerator>`_.

Reformatting pedigree
---------------------

Start by generating the reformatted pedigree file, if this has not already been done for your sample. Note you will only have to do this once!
See :ref:`Generating a reformatted pedigree file`.

Generating correlations
-----------------------

Finally, make sure the reformatted pedigree, script, and your data file are in the same directory (not all of which is strictly necessary, but makes it easier). Then generate correlations, for example:

.. code-block:: bash

    python KinshipCorrelationGenerator.py --data my_input.csv --outprefix my_output

See :ref:`Computing correlations`.


