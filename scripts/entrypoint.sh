#!/bin/bash

. /kb/deployment/user-env.sh

python3 ./scripts/prepare_deploy_cfg.py ./deploy.cfg ./work/config.properties

if [ -f ./work/token ] ; then
  export KB_AUTH_TOKEN=$(<./work/token)
fi

if [ $# -eq 0 ] ; then
  sh ./scripts/start_server.sh
elif [ "${1}" = "test" ] ; then
  echo "Run Tests"
  make test
elif [ "${1}" = "async" ] ; then
  sh ./scripts/run_async.sh
elif [ "${1}" = "init" ] ; then
  echo "Initialize module"
  # Copy reference data to /data
  cp -r ./data/* /data/
  # Sanity check: ensure a key file exists (e.g., family_centroids.npz)
  if [ -f /data/family_centroids.npz ]; then
    touch /data/__READY__
    echo "Reference data initialized and __READY__ file created."
  else
    echo "Sanity check failed: family_centroids.npz not found in /data."
    exit 1
  fi
elif [ "${1}" = "bash" ] ; then
  bash
elif [ "${1}" = "report" ] ; then
  export KB_SDK_COMPILE_REPORT_FILE=./work/compile_report.json
  make compile
else
  echo Unknown
fi
