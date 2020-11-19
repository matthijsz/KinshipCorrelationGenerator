Generating a reformatted pedigree file
======================================

Creating the input pedigree file
--------------------------------

Generating a reformatted pedigree file requires your existing pedigree file to have a specific format.
First, the file must be uncompressed, and comma-separated. The header-line must be as follows:

.. code-block:: bash

    !FamID,ID,Father,Mother,Gender,Twincode,DZtwincode,TwinHousehold3,SibHousehold2,SpouseHousehold3

Each of these columns explained in more detail below:

* :code:`!FamID`: Family identifier, unique for each family (definition of `family` is up to the user)
* :code:`ID`: Personal identifier, unique for each individual
* :code:`Father`: Personal identifier (:code:`ID`) of the father.
* :code:`Mother`: Personal identifier (:code:`ID`) of the mother.
* :code:`Gender`: Gender (M=Male, F=Female)
* :code:`TwinCode`: MZ Twin code, unique identifier of each MZ twin pair. (First MZ twin pair in the data is TwinCode 1, second pair is TwinCode 2, etc)
* :code:`DZtwincode`: DZ Twin code, unique identifier of each DZ twin pair.
* :code:`TwinHousehold3`: Twin pair code, unique identifier of each MZ DZ pair within a family.
* :code:`SibHousehold2`: Sibling household pair, unique identifier for all sets of siblings.
* :code:`SpouseHousehold3`: Spouse household pair, unique identifier for each spousal pair.


Running the code
----------------

Running the following code:

.. code-block:: bash

    python KinshipCorrelationGenerator.py --pedigree mypedigree.csv

Will generate a reformatted pedigree file :code:`reformatted_pedigree.pickle`. This reformatted pedigree file will allow you to generate twin- sibling- spouse- and parent-offspring correlations.

If you would like to generate an extended pedigree, use the following:

.. code-block:: bash

    python KinshipCorrelationGenerator.py --pedigree mypedigree.csv --extended

Will generate a reformatted extended pedigree file :code:`reformatted_extended_pedigree.pickle`. This file also contains twin- sib- spouse- and parent-offspring pairs. Therefore if you have a reformatted extended pedigree, you do not need a regular reformatted pedigree.

.. note::

    It is recommended to generate a reformatted pedigree file for the entire pedigree. The script can handle missing data, and reformatting the pedigree takes considerably longer than generating correlations.
    Therefore this reformatting only needs to be done once for an entire pedigree, and can subsequently be used for any number of correlations.


.. warning::
    Generating a pedigree file will always overwrite any existing file named :code:`reformatted_pedigree.pickle` or :code:`reformatted_extended_pedigree.pickle`

    If you have existing reformatted pedigrees in your directory that you would like to continue using, please rename this before generating a new one!


