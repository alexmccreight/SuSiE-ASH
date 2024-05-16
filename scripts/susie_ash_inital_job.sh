#!/bin/bash

username=apm2217
./src/susie_ash_initial_job.sh \
    .
-c 2 -m 16 \
--job-size 100 \
--mount your-bucket-name:/home/$username/data \
--mount your-bucket-name:/home/$username/output \
--mountOpt "mode=r" "mode=rw" \
--cwd "/home/$username/data" \
--image your-docker-image \
--entrypoint "source /usr/local/bin/_activate_current_env.sh" \
--env ENV_NAME=pecotmr \
--imageVolSize 10 \
--opcenter 23.22.157.8 \
--no-fail-fast


username=apm2217
./src/susie_ash_initial_job.sh \
 ./PATHWAY/commands_to_submit.txt \
 -c 2 -m 16 \
 --job-size 100 \
 --mount statfungen/ftp_fgc_xqtl:/home/$username/data \
-  statfungen/ftp_fgc_xqtl/sos_cache/$username:/home/$username/.sos \
-  statfungen/ftp_fgc_xqtl/analysis_result/finemapping_twas:/home/$username/output \
 --mountOpt "mode=r" "mode=rw" "mode=rw" \
 --cwd "/home/$username/data" \
 --image ghcr.io/cumc/pecotmr_docker:latest \
 --entrypoint "source /usr/local/bin/_activate_current_env.sh" \
 --env ENV_NAME=pecotmr \
 --imageVolSize 10 \
 --opcenter 23.22.157.8 \
 --download "statfungen/ftp_fgc_xqtl/ROSMAP/genotype/analysis_ready/geno_by_chrom/:/home/$username/input/" \
 --download-include "ROSMAP_NIA_WGS.leftnorm.bcftools_qc.plink_qc.1.*" \
 --ebs-mount "/home/$username/input=60" \
 --no-fail-fast  