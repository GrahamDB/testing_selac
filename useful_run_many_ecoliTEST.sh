#!/bin/bash
start=$1
stop=$2

shift 2
seq $start $stop | while read iter_seed; do
( date;
time R -q --no-save -e 'source("R_profile_helper.R", keep.source = T); run_full_selac_optimize(seed='$iter_seed', nCores =1, ref="dd94866")';
date ;) |& tee -a log_profile_ecoliTEST.txt;
sleep 5m;
( date;
time R -q --no-save -e 'source("R_profile_helper.R", keep.source = T); run_test_ecoli_optimize(seed='$iter_seed', nCores =1, ref="744c6c8")';
date ;) |& tee -a log_profile_ecoliTEST.txt;
sleep 5m;
( date;
time R -q --no-save -e 'source("R_profile_helper.R", keep.source = T); run_test_ecoli_optimize(seed='$iter_seed', nCores =1)';
date ;) |& tee -a log_profile_ecoliTEST.txt
sleep 5m;
done
