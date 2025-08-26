#!/bin/bash
script_dir=$(dirname "$(readlink -f "$0")")
export KB_DEPLOYMENT_CONFIG=$script_dir/../deploy.cfg
export KB_AUTH_TOKEN=`cat /kb/module/work/token`
export PYTHONPATH=$script_dir/../lib:$PATH:$PYTHONPATH
cd $script_dir/../test

# Use pytest instead of nose for better compatibility
python -m pytest --tb=short --cov=kbase_protein_query_module --cov-report=html:/kb/module/work/test_coverage --cov-report=term-missing -v .
