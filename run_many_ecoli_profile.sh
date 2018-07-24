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
profile_run () {
  ( date;
time R -q --no-save -e 'source("R_profile_helper.R", keep.source = T); run_ecoli_profile_mode(seed='$1', '"$sys_args"', ref="'"$2"'")';
date ;) |& tee -a "${out_file}"
sleep 30;
}
profile_run_seq () {
  profile_run "$1" "$2"
  profile_run "$1" "$3"
  profile_run "$1" "$4"
}
profile_run_combo () {
  seed_offset="$1"
  rA="dd94866"
  rB="d60d4c6"
  rC="38a4c36"
  profile_run_seq "$((seed_offset+0))" "$rA" "$rB" "$rC"
  profile_run_seq "$((seed_offset+3))" "$rB" "$rC" "$rA"
  profile_run_seq "$((seed_offset+4))" "$rC" "$rA" "$rB"
  profile_run_seq "$((seed_offset+1))" "$rA" "$rC" "$rB"
  profile_run_seq "$((seed_offset+5))" "$rC" "$rB" "$rA"
  profile_run_seq "$((seed_offset+2))" "$rB" "$rA" "$rC"
}
seq $start 6 $stop | while read iter_seed; do profile_run_combo "$iter_seed"; done

#( date;
#time R -q --no-save -e 'source("R_profile_helper.R", keep.source = T); run_ecoli_profile_mode(seed='$iter_seed', '"$sys_args"', ref="dd94866")';
#date ;) |& tee -a "${out_file}";
#sleep 30;
#( date;
#time R -q --no-save -e 'source("R_profile_helper.R", keep.source = T); run_ecoli_profile_mode(seed='$iter_seed', '"$sys_args"', ref="744c6c8")';
#date ;) |& tee -a "${out_file}";
#sleep 30;
#( date;
#time R -q --no-save -e 'source("R_profile_helper.R", keep.source = T); run_ecoli_profile_mode(seed='$iter_seed', '"$sys_args"')';
#date ;) |& tee -a "${out_file}"
#sleep 30;
#done
