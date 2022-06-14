"""
CREWS Caribbean - Import satellite data
__date__ = '20220304'
__version__ = '1.0.0'
__author__ =
        'Flavio Pignone (flavio.pignone@cimafoundation.org',
        'Andrea Libertino (andrea.libertino@cimafoundation.org',
__library__ = 'crews'
General command line:
### python crews_import_satellite_data.py -time "YYYY-MM-DD HH:MM"
Version(s):
20220304 (1.0.0) --> Beta release for CREWS Caribbean
"""
# -------------------------------------------------------------------------------------

# -------------------------------------------------------------------------------------
import os, logging, netrc, json, time
import pytz
import pandas as pd
import numpy as np
import xarray as xr
import datetime as dt
from argparse import ArgumentParser
from drops2 import coverages
from drops2.utils import DropsCredentials
import rioxarray


def main():
    # -------------------------------------------------------------------------------------
    # Version and algorithm information
    alg_name = 'CREWS Caribbean - Import satellite data '
    alg_version = '1.0.0'
    alg_release = '2022-03-04'
    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # Get algorithm settings
    alg_settings, alg_time = get_args()

    # Set algorithm settings
    data_settings = read_file_json(alg_settings)
    drops_settings = data_settings["data"]["dynamic"]["drops2"]

    # Set algorithm logging
    os.makedirs(data_settings['data']['log']['folder'], exist_ok=True)
    set_logging(logger_file=os.path.join(data_settings['data']['log']['folder'], data_settings['data']['log']['filename']))
    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # Get time settings
    date_run = dt.datetime.strptime(alg_time, "%Y-%m-%d %H:%M")
    start_run = date_run - data_settings['data']['dynamic']['time']['steps_observed'] * pd.Timedelta(
        data_settings['data']['dynamic']['time']['time_frequency'])
    steps_run = pd.date_range(start_run, date_run, closed='right',
                              freq=data_settings['data']['dynamic']['time']['time_frequency'])
    drops_settings['frequency'] = data_settings['data']['dynamic']['time']['time_frequency']
    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # Create directories
    #template_time_step = fill_template_time(data_settings['template'], date_run)
    #output_folder = data_settings['data']['dynamic']['outcome']['folder'].format(**template_time_step)
    #output_file = os.path.join(output_folder, data_settings['data']['dynamic']['outcome']['filename']).format(
    #    **template_time_step)
    #os.makedirs(output_folder, exist_ok=True)
    # -------------------------------------------------------------------------------------
    # Info algorithm
    logging.info(' ============================================================================ ')
    logging.info(' ==> START ... ')
    logging.info(' ')

    logging.info(" --> Time now : " + alg_time)

    # Time algorithm information
    start_time = time.time()
    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    logging.info(" --> Set up drops2 connection...")
    drops_settings = data_settings['data']['dynamic']['drops2']
    if not all([drops_settings['DropsUser'], drops_settings['DropsPwd']]):
        netrc_handle = netrc.netrc()
        try:
            drops_settings['DropsUser'], _, drops_settings['DropsPwd'] = netrc_handle.authenticators(
                drops_settings['DropsAddress'])
        except:
            logging.error(
                ' --> Netrc authentication file not found in home directory! Generate it or provide user and password in the settings!')
            raise FileNotFoundError(
                'Verify that your .netrc file exists in the home directory and that it includes proper credentials!')


    DropsCredentials.set(drops_settings['DropsAddress'], drops_settings['DropsUser'], drops_settings['DropsPwd'])
    logging.info(" --> Set up drops2 connection...DONE")

    for date_to in steps_run:
        logging.info("-> Computing step: " + date_to.strftime("%Y-%m-%d %H:%M"))
        date_from = date_to - pd.Timedelta(data_settings['data']['dynamic']['time']['time_frequency'])
        try:
            # Download data
            data_id = drops_settings['DropsDataId']

            logging.info(" ---> Time window for file search : " + date_from.strftime("%Y%m%d%H%M") + " - " + date_to.strftime("%Y%m%d%H%M"))
            model_dates_raw = coverages.get_dates(data_id, date_from.strftime("%Y%m%d%H%M"), date_to.strftime("%Y%m%d%H%M"))
            model_dates = [i for i in model_dates_raw if i <= pytz.utc.localize(date_to) and i >= pytz.utc.localize(date_from)]
            logging.info(" ---> Found " + str(len(model_dates)) + " files...")

            for date_now in model_dates:
                template_filled = fill_template_time(data_settings["template"], date_now)
                out_folder_time = data_settings["data"]["dynamic"]["outcome"]["folder"].format(**template_filled)
                os.makedirs(out_folder_time, exist_ok=True)
                file_out = os.path.join(out_folder_time,data_settings["data"]["dynamic"]["outcome"]["filename"]).format(**template_filled)
                if os.path.isfile(file_out) is False or data_settings["flags"]["overwrite_existing"] is True:
                    logging.info(" ---> Download step : " + date_now.strftime("%Y-%m-%d %H:%M") + "...")
                    data_drops = coverages.get_data(data_id, date_now, drops_settings['DropsVarId'], drops_settings['DropsLevel'])
                    if not data_settings["data"]["dynamic"]["satellite"]["flip_lat"]:
                        yaxis = data_drops.latitude.values[::-1,1]
                    else:
                        yaxis = data_drops.latitude.values[:, 1]
                    xaxis = data_drops.longitude.values[1, :]
                    rain_values = np.squeeze(np.flipud(data_drops[drops_settings['DropsVarId']].values))
                    if not data_settings["data"]["dynamic"]["satellite"]["additional_nulls"] is None:
                        for null_val in data_settings["data"]["dynamic"]["satellite"]["additional_nulls"]:
                            rain_values[rain_values==null_val] = np.nan
                    if data_settings["data"]["dynamic"]["satellite"]["xaxis_0-360"]:
                        rain_values = np.concatenate((rain_values[:,xaxis>180],rain_values[:,xaxis<180]), axis=1)
                        xaxis = np.concatenate((xaxis[xaxis > 180], xaxis[xaxis < 180]), axis=0)
                        xaxis[xaxis > 180] = -360 + xaxis[xaxis > 180]
                    da_rain_tif = xr.DataArray(rain_values, dims=["y", "x"], coords={"x": xaxis,"y": yaxis}).rio.write_crs("epsg:4326",inplace=True)
                    if data_settings["flags"]["crop_bbox"]:
                        mask_lon = (da_rain_tif.x >= data_settings["data"]["dynamic"]["satellite"]["bbox"]["lon_left"]) & (da_rain_tif.x <= data_settings["data"]["dynamic"]["satellite"]["bbox"]["lon_right"])
                        mask_lat = (da_rain_tif.y >= data_settings["data"]["dynamic"]["satellite"]["bbox"]["lat_bottom"]) & (da_rain_tif.y <= data_settings["data"]["dynamic"]["satellite"]["bbox"]["lat_top"])
                        da_rain_tif = da_rain_tif.where(mask_lon & mask_lat, drop=True)
                    da_rain_tif.rio.to_raster(file_out, compress="DEFLATE", dtype="float32")
                else:
                    logging.info(" ---> Step : " + date_now.strftime("%Y-%m-%d %H:%M") + "... Already downloaded!")
                    continue
        except:
            logging.info(" ---> Step : " + date_now.strftime("%Y-%m-%d %H:%M") + "... PROBLEM!")
            continue
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
if __name__ == '__main__':
    main()