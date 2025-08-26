#!/bin/bash
script_dir=$(dirname "$(readlink -f "$0")")
export KB_DEPLOYMENT_CONFIG=$script_dir/../deploy.cfg
export PYTHONPATH=$script_dir/../lib:$PATH:$PYTHONPATH
uwsgi --master --processes 5 --threads 5 --http :5000 --wsgi-file $script_dir/../lib/kbase_protein_query_module/kbase_protein_query_moduleServer.py
