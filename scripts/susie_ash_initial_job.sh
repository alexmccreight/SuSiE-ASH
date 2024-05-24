#!/bin/bash
username=apm2217
./src/mm_jobman.sh \
 ./commands_to_submit.txt \
 -c 2 -m 16 \
 --job-size 2 \
 --parallel-commands 2 \
 --mount statfungen/ftp_fgc_xqtl/interactive_analysis/$username/susie-ash-data:/home/$username/data \
 --mount statfungen/ftp_fgc_xqtl/analysis_result/interactive_analysis/$username/susie-ash-output:/home/$username/output \
 --mountOpt "mode=r" "mode=rw" \
 --cwd "/home/$username/data" \
 --image ghcr.io/cumc/pecotmr_docker:latest \
 --entrypoint "source /usr/local/bin/_activate_current_env.sh" \
 --env ENV_NAME=pecotmr \
 --imageVolSize 10 \
 --opcenter 23.22.157.8 \
 --no-fail-fast
