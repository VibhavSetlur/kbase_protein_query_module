#!/bin/bash
script_dir=$(dirname "$(readlink -f "$0")")
export KB_DEPLOYMENT_CONFIG=$script_dir/../deploy.cfg
export KB_AUTH_TOKEN=`cat /kb/module/work/token`
echo "Removing temp files..."
rm -rf /kb/module/work/tmp/*
echo "...done removing temp files."
export PYTHONPATH=$script_dir/../lib:$PATH:$PYTHONPATH
cd $script_dir/../test
python -m pytest --tb=short --cov=kbase_protein_query_module --cov-report=html:/kb/module/work/test_coverage --cov-report=term-missing -v .
