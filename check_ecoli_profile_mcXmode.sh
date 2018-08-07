#!/bin/bash
echo "$0" "$@"
prof_mode=$1
shift
echo args "$@"
sole_SEED=$1
shift
echo args "$@"
num_cores=$1
shift
echo args "$@"
gamma_mode="none"
nuc_mode="GTR"
codon_mode="selac"
sys_args="mode=\"${prof_mode}\", nuc.model=\"${nuc_mode}\", gamma.type=\"${gamma_mode}\", nCores=${num_cores}"
out_file="log_profile_ecoli${prof_mode}.txt"
profile_run () {
  ( date;
time R -q --no-save -e 'source("R_profile_helper.R", keep.source = T); run_ecoli_profile_mode(seed='$1', '"$sys_args"', ref="'"$2"'")';
date ;) |& tee -a "${out_file}"
sleep 5;
}
while [[ "$#" -gt 0 ]]; do
  profile_run "${sole_SEED}" "$1"
  shift 1
done
