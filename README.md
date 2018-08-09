# testing_selac
Code for profiling selac R package
## shell functions
Source useful.txt to load functions.  
- The wait_and_run polls every hour to see if a benchmark was running 
at some point in the last half hour.  Continues to loop until no such benchmarks are found.
Useful to run batches in separate terminals while avoiding overlaps.
### Script to install specific versions of selac if missing
./check_ecoli_profile_short.sh 2053 abcdef0 b012345 [...] 
- This takes 3 minutes per missing version, and 5 seconds per previously checked version.
- This short script installs and test runs all of the listed selac revisions, skipping 
any version that has both been previously installed and returned a finite likelihood 
value for the test for the given seed (2053 is chosen to be consistant).  
- hmm_benchmark_battery and hmm_longbenchmark_battery automatically run this for you.
### HMM quick batch benchmark for 5c98a1f revision, seeds 3010 through 3015
I expect this to take roughly an hour for that revision; preliminary estimates are four hours for ab3e84e
source useful.txt; wait_and_run hmm_benchmark_battery 5c98a1f
### HMM long batch benchmark for 5c98a1f revision, seeds 3010 through 3015
I expect this to take roughly 18 hours for that revision; preliminary estimates are 17 hours for ab3e84e
source useful.txt; wait_and_run hmm_longbenchmark_battery 5c98a1f
### HMM quick batch benchmark for 5c98a1f revision, seeds 3040 through 3045
source useful.txt; wait_and_run hmm_benchmark_battery 5c98a1f 3040

## shell scripts

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
