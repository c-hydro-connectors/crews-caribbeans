#!/bin/bash

#Export conda path (on idroclima the standard python is 2.7)
export PATH=/home/silvestro/fp_virtualenv_python3_hyde/bin:$PATH
source activate wradlib

python3 crews_preprocess_mcm_data.py -settings_file crews_preprocess_mcm_data.json -time "2020-10-30 12:00"
