#!/bin/bash

#Export conda path (on idroclima the standard python is 2.7)
export PATH=/home/silvestro/fp_virtualenv_python3_hyde/bin:$PATH
source activate wradlib

python3 crews_import_satellite_data.py -settings_file crews_import_satellite_data_gsmap_historical.json -time "2021-08-21 12:00"

python3 crews_import_satellite_data.py -settings_file crews_import_satellite_data_gsmap.json -time "2022-03-17 10:00"



python3 crews_import_satellite_data.py -settings_file crews_import_satellite_data_ghe_historical.json -time "2021-08-21 12:00"

python3 crews_import_satellite_data.py -settings_file crews_import_satellite_data_ghe.json -time "2021-08-21 12:00"


