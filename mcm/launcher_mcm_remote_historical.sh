#!/bin/bash -e

#-----------------------------------------------------------------------------------------
# Script information
script_name='Modified Conditional Merging - OPERATIONAL'
script_version="1.0.0"
script_date='2021/11/03'

#-----------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------
# Script settings
library_path=/home/silvestro/fp_libs/fp-hyde/
virtualenv_path=/home/silvestro/fp_virtualenv_python3_hyde/bin
virtualenv_name=fp_virtualenv_python3_hyde_libraries

settings_file=mcm_operational_config_remote.json
settings_file_griso=griso_historical_config.json

# Get time information (-u to get gmt time)
time_start="2013-12-23"
time_end="2013-12-26"

#-----------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------
# Activate virtualenv
export PATH=$PATH:$virtualenv_path
source activate $virtualenv_name 

# Add path to pythonpath
export PYTHONPATH=$PYTHONPATH:$library_path

#-----------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------
# Script start
echo " ==================================================================================="
echo " ==> "$script_name" (Version: "$script_version" Release_Date: "$script_date")"

date_now=$time_start #"$(date --date="$time_start 1 hour ago" '+%Y-%m-%d %H:%M')"

while [ "$date_now" != "$time_end" ]; do

date_now=$(date -I -d "$date_now + 1 day")

for hour in {0..23..1}; do
 
echo " "
echo " ==> Computing time : $date_now $hour:00"
echo " "

python3 $library_path/apps/model/mcm/hyde_data_dynamic_modified_conditional_merging.py -settings_file $settings_file -time "$date_now $hour:00" || echo " ==> Failed"

echo " "
echo " ==> Computing time : $date_now $hour:30"
echo " "

python3 $library_path/apps/model/mcm/hyde_data_dynamic_modified_conditional_merging.py -settings_file $settings_file -time "$date_now $hour:30" || echo " ==> Failed"

done
done

#-----------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------
# Info script end
echo " ==> "$script_name" (Version: "$script_version" Release_Date: "$script_date")"
echo " ==> ... END"
echo " ==> Bye, Bye"
echo " ==================================================================================="
# ----------------------------------------------------------------------------------------
