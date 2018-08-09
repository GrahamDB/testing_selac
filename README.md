# testing_selac
Code for profiling selac R package
## R_profile_helper.R
### To use: 
* source the script, which will install profvis in your home library if you haven't already
* call the function 'run_ecoli_profile_mode', specifying:
- the profile test
- the random number seed
- the codon model (in general, just use "selac")
- the nucleotide model (one of "GTR", "UNREST", or "JC")
- the commit short hash for the version of selac to benchmark (ie dd94866)
- the gamma model (one of "none", "quadrature", "lognormal", or "median")
- the number of cores (for line profiling use 1)
### run_ecoli_profile_mode profile tests
#### SHORTTEST
#### TEST
#### SHORT
#### SHORTTESTHMM
#### SHORTHMM
#### FASTHMMTEST
Single shot HMM likelihood evaluation over 10 codon triplets drawn from the E coli KOSI07 gene.
#### HMMEVAL50
Single shot HMM likelihood evaluation over 50 codon triplets drawn from the E coli KOSI07 gene.
#### HMMEVALFULL
Single shot HMM likelihood evaluation for the entire E coli KOSI07 gene.
#### FASTHMMDEBUG
In this mode profiling and the outer tryCatch is turned off, but is otherwise identical to FASTHMMTEST.
#### FASTHMMSINGLEDEBUG
In this mode lapply is used instead of mclapply to allow errors to be examined in context. This mode is
otherwise identical to FASTHMMDEBUG

## R_profile_test_likelihood.R
### Sourcing the R_profile_test_likelihood.R file in a recent version of Rstudio should: 
- install selac 1.6.1 and its dependencies in a separate library subdirectory
- install profvis in your home library if you haven't already
- generate Rprof files for each of the four selac tests

### Caveates:
- this is a draft script that has only been tested on Ubuntu bionic
- the install script assumes you have the selac source package at "~/git/selac_1.6.1.tar.gz"
