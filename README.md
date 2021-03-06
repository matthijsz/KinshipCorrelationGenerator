[![Documentation Status](https://readthedocs.org/projects/kinshipcorrelationgenerator/badge/?version=latest)](https://kinshipcorrelationgenerator.readthedocs.io/en/latest/?badge=latest)

# Patch notes

##### 25-03-2021 (v1.2.1)
 - Added argument `min_n` which was in documentation but not implemented in code.
 
##### 20-03-2021 (v1.2.0)
 - Fixed [Issue#1](https://github.com/matthijsz/weightedcorr/issues/1) in weighted correlation
 
##### 17-03-2021 (v1.1.3)
 - Fixed minor issue that would cause an error after succesful operation.

##### 17-03-2021 (v1.1.2)
 - Fixed issue that caused `randomsample` and `raw_n` to not work properly.
 - Added some lines to make it easier to work with straight from Python.
 - Added `openpyxl` to the list of dependencies to install in documentation.
 
##### 24-11-2020 (v1.1.1)
 - Fixed an issue that would cause a new reformatted extended pedigree to be saved to the wrong directory.
 - Fixed an issue that would cause columns in `--exlucde` to still be included.

##### 23-11-2020 (v1.1.0)
 - Added `--bivar` option. See [Generating cross-trait-correlations](https://kinshipcorrelationgenerator.readthedocs.io/en/latest/computecorrelation.html#generating-cross-trait-correlations) in the documentation.
 - Rearranged arguments so the order they are printed in with `-h` and `--morehelp` are in line with documentation.

##### 19-11-2020 (v1.0.0)
 - Initial release

# Description

This is the main page for the KinShipCorrelationGenerator. It can be used to quickly, and efficiently calculate phenotypic correlations in different familial pairs, (e.g. twin-correlations, parent-offspring correlations, etc) in large complex pedigrees.
By default correlations are weighted for multiple occurrences of individuals to correct for larger families. 

See the [documentation](https://kinshipcorrelationgenerator.readthedocs.io/) for more information.

# Usage

See [documentation](https://kinshipcorrelationgenerator.readthedocs.io/).

