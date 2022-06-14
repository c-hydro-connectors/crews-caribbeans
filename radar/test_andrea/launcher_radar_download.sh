#!/bin/bash

#Export conda path (on idroclima the standard python is 2.7)
export PATH=/home/silvestro/fp_virtualenv_python3_hyde/bin:$PATH
source activate wradlib

python3 crews_import_radar_data_FLA.py -settings_file crews_import_radar_data.json -time "2021-08-13 00:00"

