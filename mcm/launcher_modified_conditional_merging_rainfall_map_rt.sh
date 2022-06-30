#!/bin/bash -e

#-----------------------------------------------------------------------------------------
# Script information
script_name='LAUNCHER - MERGE REMOTE PRODUCTS'
script_version="1.0.0"
script_date='2022/06/14'

Help()
{
   # Display Help
   echo "---------------------------------------------------"
   echo "Download multiple satellite product for historicla analysis"
   echo
   echo "Syntax: ./launcher_merge_remote_rt.sh --satellite {sat_name} --radar {radar_type}"
   echo "{sat_name} available: ghe, gsmap"
   echo "{radar_type} available: polar_volume, rainfall_map"
   echo "---------------------------------------------------"
}
#-----------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------
# Script settings
system_library_folder='/home/silvestro/MSPG/libraries/'
lock_folder='/home/silvestro/MSPG/op_chain/lock/'
file_lock_init=true
#-----------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------
# Managing inputs
if [ $1 == "--satellite" ];
    then type=$2; 
    echo " --> Use satellite product: "$type
elif  [ $1 == "-h" ]; 
    then Help
    exit
else
    echo " --> ERROR! Incorrect parameter/s provided. Use flag -h for help!"
    exit
fi
if [ $3 == "--radar" ];
    then radar_type=$4;
    echo " --> Use radar type: "$radar_type
else
    radar_type="rainfall_map";
    echo " --> --radar setting not correctly specified. rainfall_map product will be considered!"
fi
#-----------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------
# Get file information
virtualenv_folder='/home/silvestro/fp_virtualenv_python3_hyde/'
virtualenv_name='fp_virtualenv_python3_hyde_libraries'

script_folder=$system_library_folder'prism'
script_file=$script_folder'/mcm/modified_conditional_merging.py'
settings_file=/home/silvestro/MSPG/op_chain/processing/modified_conditional_merging_${type}_${radar_type}_rt.json

# Get information (-u to get gmt time)
time_now=$(date -u +"%Y-%m-%d %H:%M")

year=${time_now:0:4}
month=${time_now:5:2}
day=${time_now:8:2}
hour=${time_now:11:2}
mins=${time_now:14:2}

if [ $mins -ge 30 ]
then
time_now=$(date -u +"%Y-%m-%d %H:30")
else
time_now=$(date -u +"%Y-%m-%d %H:00")
fi
#-----------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------
# Activate virtualenv
export PATH=$virtualenv_folder/bin:$PATH
source activate $virtualenv_name

# Add path to pythonpath
export PYTHONPATH="${PYTHONPATH}:$script_folder":$system_library_folder

# Add additional bins to path
#export PATH=$cdo_folder:$PATH
#-----------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------
    # Get lock information
    file_lock_start_raw="download_${type}_lock_%YYYY%MM%DD_%HH_START.txt"
    file_lock_end_raw="download_${type}_lock_%YYYY%MM%DD_%HH_END.txt"

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
    

