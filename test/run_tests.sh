#!/bin/bash
script_dir=$(dirname "$(readlink -f "$0")")
export KB_DEPLOYMENT_CONFIG=$script_dir/../deploy.cfg
export KB_AUTH_TOKEN=`cat /kb/module/work/token`
echo "Removing temp files..."
rm -rf /kb/module/work/tmp/*
echo "...done removing temp files."
export PYTHONPATH=$script_dir/../lib:$PATH:$PYTHONPATH
export KPQM_TEST_FAST=1
# Ensure callback URL is set so code paths depending on it don't attempt network
export SDK_CALLBACK_URL=${SDK_CALLBACK_URL:-}
cd $script_dir/../test
echo "Current directory: $(pwd)"
echo "Python path: $PYTHONPATH"
echo "Test files found:"
find . -name "test_*.py" | head -10
echo "Running pytest..."
python -m pytest --tb=short --cov=kbase_protein_query_module --cov-report=html:/kb/module/work/test_coverage --cov-report=term-missing -v .
