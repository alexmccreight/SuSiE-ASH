#!/bin/bash
username=apm2217
./src/mm_jobman.sh \
 ./commands_to_submit.txt \
 -c 4 -m 128 \
 --job-size 22 \
 --parallel-commands 1 \
 --mount statfungen/ftp_fgc_xqtl/interactive_sessions/$username/susie-ash-data/UKBB_data:/home/$username/data \
 --mount statfungen/ftp_fgc_xqtl/interactive_sessions/$username/LD_blocks:/home/$username/output \
 --mountOpt "mode=rw" "mode=rw" \
 --cwd "/home/$username/data" \
 --image ghcr.io/cumc/pecotmr_docker:latest \
 --entrypoint "source /usr/local/bin/_activate_current_env.sh" \
 --env ENV_NAME=pecotmr \
 --imageVolSize 10 \
 --opcenter 44.222.241.133 \
 --no-fail-fast
