#!/bin/bash
username=apm2217
./src/mm_jobman.sh \
 ./commands_to_submit.txt \
 -c 1 -m 16 \
 --job-size 1 \
 --parallel-commands 1 \
 --mount statfungen/ftp_fgc_xqtl/interactive_analysis/$username/susie-ash-data:/home/$username/data \
 --mount statfungen/ftp_fgc_xqtl/analysis_result/interactive_analysis/$username/eqtl_output2:/home/$username/output \
 --mountOpt "mode=r" "mode=rw" \
 --cwd "/home/$username/data" \
 --imageVolSize 10 \
 --opcenter 44.222.241.133 \
 --no-fail-fast
