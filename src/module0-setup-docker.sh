#!/usr/bin/env bash
# run docker at working dir
docker run --name imputation \
    -v ${PWD}:/mnt/ \
    -it thehung92phuyen/biotools:v4.0 bash