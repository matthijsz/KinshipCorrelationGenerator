Script arguments
================

-h, --help
^^^^^^^^^^

Prints a short help summary with overview of arguments

--morehelp
^^^^^^^^^^

Print more help on a specific function, or specific functions. Pass the argument name, or argument names (comma seperated), or :code:`all` for a more detailed description of the whole script and its arguments.

--data
^^^^^^

Path to the datafile. This datafile should at least contain the columns :code:`FISNumber` and :code:`age`. See :ref:`Creating the input data file`.

--outprefix
^^^^^^^^^^^

Prefix to use for output files.

--extended
^^^^^^^^^^

Add calculation of extended-pedigree correlations (UncleAunt, AuntNephew, NephewNiece, etc) when used to compute correlations. Generates an extended reformatted pedigree when together with :ref:`--pedigree`

--pedigree
^^^^^^^^^^

Path to the raw pedigree file, only used when generating a new reformatted pedigree file. See :ref:`Generating a reformatted pedigree file`.

--method
^^^^^^^^

Method for computing correlation, should be Pearson, or Spearman, default is Pearson.

--correct
^^^^^^^^^

Formula to correct phenotypes before computing correlations. Defaults to no correction

--exclude
^^^^^^^^^

Phenotypes for which no correlations should be calculated. For example, when custom covariates are used. Identifier column, age and sex columns are always excluded.

--raw_n
^^^^^^^

Return an additional .csv file containing the raw number of samples in each kinship in addition to the weighted N file.

--randomsample
^^^^^^^^^^^^^^

Use only 1 pair per family, instead of calculating weighted correlations.

--use_repeated_families
^^^^^^^^^^^^^^^^^^^^^^^

Include all participants from larger families, i.e. don't drop or weigh for duplicate samples or larger families.


--longitudinal
^^^^^^^^^^^^^^

When the input data is of longitudinal nature, this script can perform family-based selection, such that the difference in the year of survey completion among family-members is minimal.
If you use this option the input data should contain the columns :code:`index` (within-subject number of the survey, starting at 1 for fist survey), and :code:`Source` (string label of the survey). Additionally, you should also add the argument :ref:`--surveycompletion`.

--surveycompletion
^^^^^^^^^^^^^^^^^^

Dat file with survey completion years. Should contain the columns :code:`FISNumber` (identical to datafile), :code:`Source` (identical to datafile), and :code:`invjr` (year of survey completion).

Additonal non-argument settings
^^^^^^^^^^^^^^^^^^
The top of the Python file contains some additional settings you can tweak.
* :code:`upper_boundary`: Upper age boundary for inclusion in :code:`--longitudinal`. default = 110
* :code:`lower_boundary`: Lower age boundary for inclusion in :code:`--longitudinal`. default = 0
* :code:`check_cutoff_drops`: save an Age-cutoff_drops.txt file detailing subjects dropped by the cutoffs. default = False
* :code:`seed`: seed for random selection of subjects with :code:`--randomsample`. default = 1415926536
* :code:`explore_plot`: Generate scatterplots of each correlation. default = False:
* :code:`save_separate_data`: Store the raw data used for each correlation in separate files. default = False
* :code:`parallel`: Generate reformatted pedigree using multiple processing threads (Currently not working). default = False


