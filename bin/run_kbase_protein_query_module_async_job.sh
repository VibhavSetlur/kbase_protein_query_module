#!/bin/bash
script_dir=$(dirname "$(readlink -f "$0")")
export PYTHONPATH=$script_dir/../lib:$PATH:$PYTHONPATH
python -u $script_dir/../lib/kbase_protein_query_module/kbase_protein_query_moduleServer.py $1 $2 $3
