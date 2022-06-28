#!/bin/bash -e

#-----------------------------------------------------------------------------------------
# Script information
script_name='LAUNCHER - RADAR POLAR VOLUME DOWNLOADER'
script_version="1.0.0"
script_date='2022/06/14'

#-----------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------
# Script settings
system_library_folder='/home/idrologia/libraries/'
lock_folder='/home/idrologia/op_chain/lock/'
file_lock_init=true
#-----------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------
# Get file information
virtualenv_folder=$system_library_folder'python/'
virtualenv_name='wradlib'

script_folder=$system_library_folder'crews-caribbeans'
script_file=$script_folder'/radar/crews_import_radar_data_polar_volume.py'
settings_file=/home/idrologia/op_chain/preprocessing/radar/crews_import_radar_polar_volume.json

# Get information (-u to get gmt time)
time_now="2021-09-30 23:00"

year=${time_now:0:4}
month=${time_now:5:2}
day=${time_now:8:2}
hour=${time_now:11:2}
#-----------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------
# Activate virtualenv
export PATH=$virtualenv_folder/bin:$PATH
source activate $virtualenv_name

# Add path to pythonpath
export PYTHONPATH="${PYTHONPATH}:$script_folder"

# Add additional bins to path
#export PATH=$cdo_folder:$PATH
#-----------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------
# Preliminary check for forecast file existence
if [ $type == "nwp_gfs-det" ];
    then 
    for f in $nwp_gfs_files; do 
        [ -e "$f" ] && echo " --> CHECK FORECAST FILE... OK!" || { echo "CHECK FORECAST FILE... ERROR! Forecast file does not exist!"; exit 1; }
        break
        done
fi
#-----------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------
    # Get lock information
    file_lock_start_raw="downloader_${type}_lock_%YYYY%MM%DD_%HH_START.txt"
    file_lock_end_raw="downloader_${type}_lock_%YYYY%MM%DD_%HH_END.txt"

    folder_lock_raw=$lock_folder'preprocessing/'
    
    # ----------------------------------------------------------------------------------------
    # Define path data
    folder_lock_def=${folder_lock_raw/"%YYYY"/$year}
    folder_lock_def=${folder_lock_def/"%MM"/$month}
    folder_lock_def=${folder_lock_def/"%DD"/$day}

    file_lock_start_def=${file_lock_start_raw/"%YYYY"/$year}
    file_lock_start_def=${file_lock_start_def/"%MM"/$month}
    file_lock_start_def=${file_lock_start_def/"%DD"/$day}
    file_lock_start_def=${file_lock_start_def/"%HH"/$hour}

    file_lock_end_def=${file_lock_end_raw/"%YYYY"/$year}
    file_lock_end_def=${file_lock_end_def/"%MM"/$month}
    file_lock_end_def=${file_lock_end_def/"%DD"/$day}
    file_lock_end_def=${file_lock_end_def/"%HH"/$hour}

    # Create folder(s)
    if [ ! -d "$folder_lock_def" ]; then
    	mkdir -p $folder_lock_def
    fi
    # ----------------------------------------------------------------------------------------	

#-----------------------------------------------------------------------------------------
# Info script start
echo " ==================================================================================="
echo " ==> "$script_name" (Version: "$script_version" Release_Date: "$script_date")"
echo " ==> START ${dom} ..."

# Execution pid
execution_pid=$$

#-----------------------------------------------------------------------------------------
# File lock definition
path_file_lock_def_start=$folder_lock_def/$file_lock_start_def 
path_file_lock_def_end=$folder_lock_def/$file_lock_end_def

# Init lock conditions
echo " ====> INITILIZE LOCK FILES ... "
if $file_lock_init; then
# Delete lock files
if [ -f "$path_file_lock_def_start" ]; then
   rm "$path_file_lock_def_start"
fi
if [ -f "$path_file_lock_def_end" ]; then
   rm "$path_file_lock_def_end"
fi
echo " ====> INITILIZE LOCK FILES ... DONE!"
else
echo " ====> INITILIZE LOCK FILES ... SKIPPED!"
fi
#-----------------------------------------------------------------------------------------  

#-----------------------------------------------------------------------------------------
# Run check
if [ -f $path_file_lock_def_start ] && [ -f $path_file_lock_def_end ]; then   
#-----------------------------------------------------------------------------------------
# Process completed
echo " ===> EXECUTION SKIPPED! ALL DATA WERE PROCESSED DURING A PREVIOUSLY RUN"
#-----------------------------------------------------------------------------------------

elif [ -f $path_file_lock_def_start ] && [ ! -f $path_file_lock_def_end ]; then
#-----------------------------------------------------------------------------------------
# Process running condition
echo " ===> EXECUTION SKIPPED! SCRIPT IS STILL RUNNING ... WAIT FOR PROCESS END"
#-----------------------------------------------------------------------------------------

elif [ ! -f $path_file_lock_def_start ] && [ ! -f $path_file_lock_def_end ]; then
     
 #-----------------------------------------------------------------------------------------
 # Lock File START
 time_step=$(date +"%Y-%m-%d %H:%S")
 echo " ================================ " >> $path_file_lock_def_start
 echo " ==== EXECUTION START REPORT ==== " >> $path_file_lock_def_start
 echo " "
 echo " ==== PID:" $execution_pid >> $path_file_lock_def_start
 echo " ==== Algorithm: $script_name" >> $path_file_lock_def_start
 echo " ==== RunTime: $time_step" >> $path_file_lock_def_start
 echo " ==== ExecutionTime: $time_now" >> $path_file_lock_def_start
 echo " ==== Status: RUNNING" >> $path_file_lock_def_start
 echo " "
 echo " ================================ " >> $path_file_lock_def_start

# Run python script (using setting and time)
if python3 $script_file -settings_file $settings_file -time "$time_now" 
then
    time_step=$(date +"%Y-%m-%d %H:%S")
    echo " ============================== " >> $path_file_lock_def_end
    echo " ==== EXECUTION END REPORT ==== " >> $path_file_lock_def_end
    echo " "
    echo " ==== PID:" $execution_pid >> $path_file_lock_def_end
    echo " ==== Algorithm: $script_name" >> $path_file_lock_def_end
    echo " ==== RunTime: $time_step" >> $path_file_lock_def_end
    echo " ==== ExecutionTime: $time_now" >> $path_file_lock_def_end
    echo " ==== Status: COMPLETED" >>  $path_file_lock_def_end
    echo " "
    echo " ============================== " >> $path_file_lock_def_end
else
    rm $path_file_lock_def_start
    echo " ===> EXECUTION FAILED! SCRIPT HAS CRASHED!"
    exit 1
fi

else 
#-----------------------------------------------------------------------------------------
# Exit unexpected mode
echo " ===> EXECUTION FAILED! SCRIPT ENDED FOR UNKNOWN LOCK FILES CONDITION!"
#-----------------------------------------------------------------------------------------

fi

# Info script end
echo " ==> "$script_name" (Version: "$script_version" Release_Date: "$script_date")"
echo " ==> ... END"
echo " ==> Bye, Bye"
echo " ==================================================================================="
# ----------------------------------------------------------------------------------------
    

