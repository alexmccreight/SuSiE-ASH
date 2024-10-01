#!/bin/bash
username=apm2217
./src/mm_jobman.sh \
 ./commands_to_submit.txt \
 -c 16 -m 32 \
 --job-size 45 \
 --parallel-commands 9 \
 --mount statfungen/ftp_fgc_xqtl/interactive_analysis/$username/susie-ash-data:/home/$username/data \
 --mount statfungen/ftp_fgc_xqtl/analysis_result/interactive_analysis/$username/eqtl_output:/home/$username/output \
 --mountOpt "mode=r" "mode=rw" \
 --cwd "/home/$username/data" \
 --image ghcr.io/cumc/pecotmr_docker:latest \
 --entrypoint "source /usr/local/bin/_activate_current_env.sh" \
 --env ENV_NAME=pecotmr \
 --imageVolSize 10 \
 --opcenter 44.222.241.133 \
 --no-fail-fast
