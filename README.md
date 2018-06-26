# testing_selac
Code for profiling selac R package

## Sourcing the R_profile_test_likelihood.R file in a recent version of Rstudio should: 
- install selac 1.6.1 and its dependencies in a separate library subdirectory
- install profvis in your home library if you haven't already
- generate Rprof files for each of the four selac tests

## Caveates:
- this is a draft script that has only been tested on Ubuntu bionic
- the install script assumes you have the selac source package at "~/git/selac_1.6.1.tar.gz"
