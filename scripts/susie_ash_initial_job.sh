#!/bin/bash
username=apm2217
./src/mm_jobman.sh \
 ./commands_to_submit.txt \
 -c 4 -m 8 \
 --job-size 72 \
 --parallel-commands 9 \
 --mount gao-851725442056/ftp_fgc_xqtl/interactive_analysis/$username/susie-ash-data:/home/$username/data \
 --mount gao-851725442056/ftp_fgc_xqtl/analysis_result/interactive_analysis/$username/simulation_output:/home/$username/output \
 --mountOpt "mode=r" "mode=rw" \
 --cwd "/home/$username/data" \
 --image ghcr.io/cumc/pecotmr_docker:latest \
 --entrypoint "source /usr/local/bin/_activate_current_env.sh" \
 --env ENV_NAME=pecotmr \
 --imageVolSize 10 \
 --opcenter 44.222.241.133 \
 --no-fail-fast
