Computing correlations
======================

Creating the input data file
----------------------------

First you will need to create an input data file. This file should be a comma-seperated file with a header and the following columns:

* :code:`FISNumber`: Personal identifier
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

Generating cross-trait correlations
-----------------------------------

To generate cross-trait (ie. bivariate) correlations of the combinations of your phenotypes (i.e. twin1 trait1 - twin1 trait2, mother trait1, son trait2, etc.) you can use the bivar option:

.. code-block:: bash

    python KinshipCorrelationGenerator.py --data mydata.csv --outprefix my_output --bivar

By adding this argument the script will now save two excel (.xlsx) files, one for the correlations, one for the sum of weights (or N depeding on your other :ref:`Script arguments`) named :code:`my_output_bivar_Fam_correlations.xlsx` and :code:`my_output_bivar_Fam_N.xlsx`.
These tables will have 1 sheet per kinship (1 for MZM, 1 for MZF, etc), and these sheets contain the full correlation matrix of all phenotypes included in your input file. The diagonal of these matrices are the standard within-phenotype kinship correlations (same as the standard output), the off-diagonal are the cross-trait cross-kinship correlations.
Note these tables are not symmetrical, as they represent different things. Columns are phenotypes in ID_0, and rows are phenotypes in ID_1. As a specific example: in the MotherSon correlation table, column x, row y represents mother trait x and son trait y. Conversiley, column y row x represents the correlation of son's trait x with mother's trait y.

.. note::
    The calculation for values on the diagonal in each sheet of the bivariate table is identical to that of the output without the :code:`--bivar` option. Therefore if this option is enabled no csv files are stored, only the excel files.


Other common options
--------------------

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



