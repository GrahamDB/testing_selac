#!/bin/bash
echo "$0" "$@"
prof_mode=$1
shift 1
echo args "$@"
start=$1
stop=${2:-$start}
shift 
shift
echo args "$@"
gamma_mode="${1:-quadrature}"
[[ -z ${gamma_mode} ]] && gamma_mode="quadrature"
shift 1
echo args "$@"
nuc_mode="${1:-UNREST}"
[[ -z ${nuc_mode} ]] && nuc_mode="UNREST"
shift 1
echo args "$@"
codon_mode="selac"
sys_args="mode=\"${prof_mode}\", nuc.model=\"${nuc_mode}\", gamma.type=\"${gamma_mode}\", nCores=1"
out_file="log_profile_ecoli${prof_mode}.txt"
seq $start $stop | while read iter_seed; do
( date;
time R -q --no-save -e 'source("R_profile_helper.R", keep.source = T); run_ecoli_profile_mode(seed='$iter_seed', '"$sys_args"', ref="dd94866")';
date ;) |& tee -a "${out_file}";
sleep 30;
( date;
time R -q --no-save -e 'source("R_profile_helper.R", keep.source = T); run_ecoli_profile_mode(seed='$iter_seed', '"$sys_args"', ref="744c6c8")';
date ;) |& tee -a "${out_file}";
sleep 30;
( date;
time R -q --no-save -e 'source("R_profile_helper.R", keep.source = T); run_ecoli_profile_mode(seed='$iter_seed', '"$sys_args"')';
date ;) |& tee -a "${out_file}"
sleep 30;
done
