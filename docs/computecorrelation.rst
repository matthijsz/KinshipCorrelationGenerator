Computing correlations
======================

Creating the input data file
----------------------------

First you will need to create an input data file. This file should be a comma-seperated file with a header and the following columns:

* * :code:FISNumber`: Personal identifier
* :code:`age`: Age

As well as any additional phenotype columns. Any column other than FISNumber and age will be treated as a phenotype, unless otherwise specified using :code:`--exclude`.

The following names for phenotype columns are not allowed:
:code:`FISNumber`, :code:`sex`, :code:`age`, :code:`Source`, :code:`index`. As these names are used internally.

Generating weighted correlations
--------------------------------

By default this script will assign weights to each pair of observations (within kinship within phenotype) depending on the number of occurrences of the personal identifier of both members of that pair to prevent bias from larger families.
Pairs where both phenotypic values are missing are excluded from weight calculations.
Other options for dealing with larger families are available, see :ref:`Script arguments`.

To generate the kinsip correlations run the following:

.. code-block:: bash

    python KinshipCorrelationGenerator.py --data mydata.csv --outprefix my_output

This will generate 2 files: :code:`my_output_Fam_correlations.csv`, which contains 1 column per phenotype containing all the phenotypic correlations, and :code:`my_output_Fam_N.csv`
which contains the sum of weights by default.

To generate the extended kinsip correlations run the following:

.. code-block:: bash

    python KinshipCorrelationGenerator.py --data mydata.csv --outprefix my_output --extended


Common options
--------------

By default the script calculates weighted Pearson correlation, you can change this to a weighted Spearman correlation by specifying method.

.. code-block:: bash

    python KinshipCorrelationGenerator.py --data mydata.csv --outprefix my_output --extended --method spearman

To linearly correct your phenotypes for age before calculating the correlations, you can specify correct.

.. code-block:: bash

    python KinshipCorrelationGenerator.py --data mydata.csv --outprefix my_output --extended --correct age

Correct accepts many combinations of input, for example if you want to correct for age and sex and their interaction use:

.. code-block:: bash

    python KinshipCorrelationGenerator.py --data mydata.csv --outprefix my_output --extended --correct age+sex+age*sex

You can also use custom covariates as long as they are present with the same name in your input file. If you do so make sure to also add custom covariates to :code:`--exclude`!

.. code-block:: bash

    python KinshipCorrelationGenerator.py --data mydata.csv --outprefix my_output --extended --correct age+bmi --exclude bmi

By default, no correlation is computed (and NA returned) when there are less than 30 complete pairs. You can change this to 15 (for example) as follows:

.. code-block:: bash

    python KinshipCorrelationGenerator.py --data mydata.csv --outprefix my_output --extended --min_n 15



