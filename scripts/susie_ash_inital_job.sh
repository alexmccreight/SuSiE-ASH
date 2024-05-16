#!/bin/bash

username=apm2217

-c 2 -m 16 \
--job-size 100 \
--mount your-bucket-name:/home/$username/data \
--mount your-bucket-name:/home/$username/output \
--mountOpt "mode=r" "mode=rw" \
--cwd "/home/$username/data" \
--image your-docker-image \
--entrypoint "source /usr/local/bin/_activate_current_env.sh" \
--env ENV_NAME=your-env-name \
--imageVolSize 10 \
--opcenter 23.22.157.8 \
--no-fail-fast