#!/bin/bash

#Export conda path (on idroclima the standard python is 2.7)
export PATH=/home/silvestro/fp_virtualenv_python3_hyde/bin:$PATH
source activate wradlib

#python3 crews_merge_to_remote.py -settings_file crews_merge_to_remote_ghe.json -time "2021-08-21 12:00"

python3 crews_merge_to_remote.py -settings_file crews_merge_to_remote_ghe_historical.json -time "2021-08-10 00:00"


python3 crews_merge_to_remote.py -settings_file crews_merge_to_remote_gsmap.json -time "2021-08-10 00:00"

python3 crews_merge_to_remote.py -settings_file crews_merge_to_remote_gsmap_historical.json -time "2021-08-10 00:00"

