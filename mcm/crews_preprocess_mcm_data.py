"""
CREWS Caribbean - Preprocess MCM data
__date__ = '20220228'
__version__ = '1.0.0'
__author__ =
        'Flavio Pignone (flavio.pignone@cimafoundation.org',
        'Andrea Libertino (andrea.libertino@cimafoundation.org',
__library__ = 'crews'
General command line:
### python crews_download_gauge_data.py -time "YYYY-MM-DD HH:MM"
Version(s):
20220120 (1.0.0) --> Beta release for CREWS Caribbean
"""
# -------------------------------------------------------------------------------------

# -------------------------------------------------------------------------------------
from drops2 import sensors
from drops2.utils import DropsCredentials
import numpy as np
import pandas as pd
import datetime as dt
import xarray as xr
import logging, json, warnings, os, time
import netrc
from argparse import ArgumentParser
import rioxarray
import rasterio as rio

# -------------------------------------------------------------------------------------

# -------------------------------------------------------------------------------------
def main():

    # -------------------------------------------------------------------------------------
    # Version and algorithm information
    alg_name = 'CREWS Caribbean - Import radar data '
    alg_version = '1.0.0'
    alg_release = '2022-01-20'
    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # Get algorithm settings
    alg_settings, alg_time = get_args()

    # Set algorithm settings
    data_settings = read_file_json(alg_settings)

    # Set system settings
    warnings.filterwarnings('ignore')
    drops_settings = data_settings["data"]["dynamic"]["gauge"]["drops2"]
    if not all([drops_settings['DropsUser'], drops_settings['DropsPwd']]):
        netrc_handle = netrc.netrc()
        try:
            drops_settings['DropsUser'], _, drops_settings['DropsPwd'] = netrc_handle.authenticators(drops_settings['DropsAddress'])
        except:
            logging.error('--> Valid netrc authentication file not found in home directory! Generate it or provide user and password in the settings!')
            raise FileNotFoundError('Verify that your .netrc file exists in the home directory and that it includes proper credentials!')

    # Set algorithm logging
    os.makedirs(data_settings['data']['log']['folder'], exist_ok=True)
    set_logging(logger_file=os.path.join(data_settings['data']['log']['folder'], data_settings['data']['log']['filename']))
    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # Info algorithm
    logging.info(' ============================================================================ ')
    logging.info(' ==> START ... ')
    logging.info(' ')

    # Time algorithm information
    start_time = time.time()
    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # Get algorithm settings
    date_run = dt.datetime.strptime(alg_time, "%Y-%m-%d %H:%M")
    start_run = date_run - data_settings['data']['dynamic']['time']['steps_observed'] * pd.Timedelta(data_settings['data']['dynamic']['time']['time_frequency'])
    steps_run = pd.date_range(start_run,date_run,closed='right',freq=data_settings['data']['dynamic']['time']['time_frequency'])
    drops_settings['frequency'] = data_settings['data']['dynamic']['time']['time_frequency']

    for time_now in steps_run:
        missing_map = False
        logging.info("-> Computing step: " + time_now.strftime("%Y-%m-%d %H:%M"))
        template_filled = fill_template_time(data_settings["template"], time_now)
        try: 
            time_now_start = time_now - pd.Timedelta(data_settings['data']['dynamic']['time']['time_frequency'])
            logging.info("--> Preparing rain gauge data...")
            df_pluvio = dload_drops_gauges(time_now_start, time_now, drops_settings)
            df_pluvio.loc["value"] = df_pluvio.loc["value"] * data_settings['data']['dynamic']['gauge']['scale_factor']
            logging.info("---> Writing rain gauge output...")
            folder_out_gauge = data_settings['data']['dynamic']['gauge']['output']['folder'].format(**template_filled)
            filename_out_gauge = data_settings['data']['dynamic']['gauge']['output']['filename'].format(**template_filled)
            os.makedirs(folder_out_gauge, exist_ok=True)
            df_pluvio.T.to_csv(os.path.join(folder_out_gauge, filename_out_gauge), index=False, header=True)
            logging.info("--> Preparing rain gauge data...DONE")
        except:
            logging.warning("---> WARNING! No gauge at " + time_now_start.strftime("%Y-%m-%d %H:%M"))
            folder_out_gauge = data_settings['data']['dynamic']['gauge']['output']['folder'].format(**template_filled)
            filename_out_gauge = data_settings['data']['dynamic']['gauge']['output']['filename'].format(**template_filled)
            
        logging.info("--> Preparing radar data...")
        logging.info("---> Reading radar maps...")
        radar_steps_in = pd.date_range(time_now_start,time_now,closed='right',freq=data_settings['data']['dynamic']['radar']["local"]['time_frequency'])
        first_step = True
        counter=0
        for time_now_rad in radar_steps_in:
            template_filled = fill_template_time(data_settings["template"], time_now_rad)
            file_now_radar = os.path.join(data_settings['data']['dynamic']['radar']['local']["folder"],data_settings['data']['dynamic']['radar']['local']["filename"]).format(**template_filled)

            try:
                da_radar = rioxarray.open_rasterio(file_now_radar)
                counter = counter + 1
                da_radar_values = rioxarray.open_rasterio(file_now_radar).values
                da_radar_values[da_radar_values > data_settings['data']['dynamic']['radar']['local']["max_valid"]] = np.nan

                if first_step is True:
                    rad_rain_value = da_radar_values * data_settings['data']['dynamic']['radar']['scale_factor']
                    first_step = False
                else:
                    rad_rain_value += da_radar_values * data_settings['data']['dynamic']['radar']['scale_factor']
            
            except rio.errors.RasterioIOError:
                # If a radar map is missing it skips the time step
                logging.warning("---> WARNING! Radar map " + time_now_rad.strftime("%Y-%m-%d %H:%M") + " not found!")
                pass
                #missing_map = True
                #logging.warning("---> WARNING! Radar map " + time_now.strftime("%Y-%m-%d %H:%M") + " is missing! Skip to next time step")
                #break
             

 #       if missing_map is True:
        if counter<=data_settings['data']['dynamic']['radar']['min_step_acceptable']:
            continue
        else:
            da_radar.values = np.nan_to_num(rad_rain_value, copy=True, nan=-9999.0)
            da_radar.rio.set_nodata(-9999.0, inplace=True)
            logging.info("---> Writing radar output...")
            template_filled = fill_template_time(data_settings["template"], time_now)
            folder_out_radar = data_settings['data']['dynamic']['radar']['output']['folder'].format(**template_filled)
            filename_out_radar = data_settings['data']['dynamic']['radar']['output']['filename'].format(**template_filled)
            os.makedirs(folder_out_gauge, exist_ok=True)
            da_radar.rio.to_raster(os.path.join(folder_out_radar, filename_out_radar), compress="DEFLATE", dtype="float32")
            logging.info("--> Preparing radar data...DONE")

    # -------------------------------------------------------------------------------------
    # Info algorithm
    time_elapsed = round(time.time() - start_time, 1)

    logging.info(' ')
    logging.info(' ==> ' + alg_name + ' (Version: ' + alg_version + ' Release_Date: ' + alg_release + ')')
    logging.info(' ==> TIME ELAPSED: ' + str(time_elapsed) + ' seconds')
    logging.info(' ==> ... END')
    logging.info(' ==> Bye, Bye')
    logging.info(' ============================================================================ ')
    # -------------------------------------------------------------------------------------

# -------------------------------------------------------------------------------------
# Function for downloading rain fauge data with drops2
def dload_drops_gauges(time_start, time_end,drops_settings):
    logging.info("---> Querying drops dataset...")
    DropsCredentials.set(drops_settings["DropsAddress"], drops_settings["DropsUser"], drops_settings["DropsPwd"])
    sensors_list_P = sensors.get_sensor_list(drops_settings["DropsSensor"],
                                             geo_win=(drops_settings['lon_left'], drops_settings['lat_bottom'], drops_settings['lon_right'], drops_settings['lat_top']),
                                             group=drops_settings["DropsGroup"])
    date_from = time_start.strftime("%Y%m%d%H%M")
    date_to = time_end.strftime("%Y%m%d%H%M")

    df_pluvio = sensors.get_sensor_data(drops_settings["DropsSensor"], sensors_list_P, date_from, date_to, as_pandas=True)
    df_pluvio = df_pluvio.dropna('columns', how='all')

    df_pluvio_now = pd.DataFrame(index=['name','lon','lat','value'], columns = df_pluvio.columns)

    for col in df_pluvio.columns:
        series = df_pluvio[col]
        series[series<0] = np.nan
        series = series.dropna('rows', how='all')

        df_pluvio_now.loc['lat', col] = [i for i in sensors_list_P.list if i.id == col][0].lat
        df_pluvio_now.loc['lon', col] = [i for i in sensors_list_P.list if i.id == col][0].lng
        df_pluvio_now.loc['name', col] = [i for i in sensors_list_P.list if i.id == col][0].name

        if (len(series>0)) and (date_to in series.keys()):
            df_pluvio_now.loc['value',col] = series.resample(drops_settings['frequency'], closed='right', label='right').agg(pd.DataFrame.sum)[date_to]
            # df_pluvio_now.loc['value',col] = series.resample(drops_settings['frequency'], closed='right', label='right')
        else:
            df_pluvio_now.loc['value', col] = np.nan

    return df_pluvio_now.dropna('columns', how='any')
# -------------------------------------------------------------------------------------

# -------------------------------------------------------------------------------------
def fill_template_time(template_empty, time_now):
    template_filled = {}
    for key in template_empty.keys():
        template_filled[key] = time_now.strftime(template_empty[key])
    return template_filled


# -------------------------------------------------------------------------------------
# Method to read file json
def read_file_json(file_name):

    env_ws = {}
    for env_item, env_value in os.environ.items():
        env_ws[env_item] = env_value

    with open(file_name, "r") as file_handle:
        json_block = []
        for file_row in file_handle:

            for env_key, env_value in env_ws.items():
                env_tag = '$' + env_key
                if env_tag in file_row:
                    env_value = env_value.strip("'\\'")
                    file_row = file_row.replace(env_tag, env_value)
                    file_row = file_row.replace('//', '/')

            # Add the line to our JSON block
            json_block.append(file_row)

            # Check whether we closed our JSON block
            if file_row.startswith('}'):
                # Do something with the JSON dictionary
                json_dict = json.loads(''.join(json_block))
                # Start a new block
                json_block = []

    return json_dict
# -------------------------------------------------------------------------------------

# -------------------------------------------------------------------------------------
# Method to get script argument(s)
def get_args():
    parser_handle = ArgumentParser()
    parser_handle.add_argument('-settings_file', action="store", dest="alg_settings")
    parser_handle.add_argument('-time', action="store", dest="alg_time")
    parser_values = parser_handle.parse_args()

    if parser_values.alg_settings:
        alg_settings = parser_values.alg_settings
    else:
        alg_settings = 'configuration.json'

    if parser_values.alg_time:
        alg_time = parser_values.alg_time
    else:
        alg_time = None

    return alg_settings, alg_time
# -------------------------------------------------------------------------------------

# -------------------------------------------------------------------------------------
# Method to set logging information
def set_logging(logger_file='log.txt', logger_format=None):
    if logger_format is None:
        logger_format = '%(asctime)s %(name)-12s %(levelname)-8s ' \
                        '%(filename)s:[%(lineno)-6s - %(funcName)20s()] %(message)s'

    # Remove old logging file
    if os.path.exists(logger_file):
        os.remove(logger_file)

    # Set level of root debugger
    logging.root.setLevel(logging.INFO)

    # Open logging basic configuration
    logging.basicConfig(level=logging.INFO, format=logger_format, filename=logger_file, filemode='w')

    # Set logger handle
    logger_handle_1 = logging.FileHandler(logger_file, 'w')
    logger_handle_2 = logging.StreamHandler()
    # Set logger level
    logger_handle_1.setLevel(logging.INFO)
    logger_handle_2.setLevel(logging.INFO)
    # Set logger formatter
    logger_formatter = logging.Formatter(logger_format)
    logger_handle_1.setFormatter(logger_formatter)
    logger_handle_2.setFormatter(logger_formatter)
    # Add handle to logging
    logging.getLogger('').addHandler(logger_handle_1)
    logging.getLogger('').addHandler(logger_handle_2)
# -------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------
def kilometers2degrees(kilometer, radius=6371):
    import math
    return kilometer / (2.0 * radius * math.pi / 360.0)

# ----------------------------------------------------------------------------

# ----------------------------------------------------------------------------
if __name__ == '__main__':
    main()

# ----------------------------------------------------------------------------

# ----------------------------------------------------------------------------
