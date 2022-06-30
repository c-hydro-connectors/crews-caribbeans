"""
CREWS Caribbean - Import satellite data
__date__ = '20220309'
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
    alg_name = 'CREWS Caribbean - Merge radar and satellite GHE '
    alg_version = '1.0.0'
    alg_release = '2022-03-09'
    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # Get algorithm settings
    alg_settings, alg_time = get_args()

    # Set algorithm settings
    data_settings = read_file_json(alg_settings)
    
    # Set algorithm logging
    os.makedirs(data_settings['data']['log']['folder'], exist_ok=True)
    set_logging(logger_file=os.path.join(data_settings['data']['log']['folder'], data_settings['data']['log']['filename']))
    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # Get time settings
    date_run = dt.datetime.strptime(alg_time, "%Y-%m-%d %H:%M")
    start_run = date_run - data_settings['data']['dynamic']['time']['steps_observed'] * pd.Timedelta(data_settings['data']['dynamic']['time']['time_frequency'])
    steps_run = pd.date_range(start_run, date_run, closed='right',freq=data_settings['data']['dynamic']['time']['time_frequency'])
    
    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # Create directories
    template_time_step = fill_template_time(data_settings['template'], date_run)

    output_folder = data_settings['data']['dynamic']['outcome']['folder'].format(**template_time_step)
    output_file = os.path.join(output_folder, data_settings['data']['dynamic']['outcome']['filename']).format(**template_time_step)
    os.makedirs(output_folder, exist_ok=True)
    
    # -------------------------------------------------------------------------------------
    # Info algorithm
    logging.info(' ============================================================================ ')
    logging.info(' ==> START ... ')
    logging.info(' ')
    logging.info(" --> Time now : " + alg_time)

    # Time algorithm information
    start_time = time.time()
    # -------------------------------------------------------------------------------------

    griglia =  rioxarray.open_rasterio(data_settings['data']['static']['grid'])     
    
    # pesi per unione
    griglia_lon,griglia_lat = np.meshgrid(griglia.x.values,griglia.y.values)
    center_x = np.min(griglia.x.values) + ((np.max(griglia.x.values)-np.min(griglia.x.values))/2)
    center_y = np.min(griglia.y.values) + ((np.max(griglia.y.values) - np.min(griglia.y.values)) / 2)
    distance = ((griglia_lon-center_x)**2+(griglia_lat-center_y)**2)**.5
    distance_km= distance/0.009                                  
    pesi_rad=1-(distance_km*0.004)
    pesi_rad[pesi_rad<0]=0   
    temp_rad=np.zeros(griglia.shape)
    temp_rad[0]=pesi_rad
    pesi_radar = griglia.copy()
    pesi_radar.values[0].astype(np.float32)
    pesi_radar.values=temp_rad
    pesi_sat=distance_km*0.004
    pesi_sat[pesi_sat>1]=1
    temp_sat=np.zeros(griglia.shape)
    temp_sat[0]=pesi_sat
    pesi_satellite = griglia.copy()
    pesi_satellite.values[0].astype(np.float32)
    pesi_satellite.values=temp_sat   
        
    for time_now in steps_run:
        logging.info("-> Computing step: " + time_now.strftime("%Y-%m-%d %H:%M"))
        template_filled = fill_template_time(data_settings["template"], time_now)
        output_folder = data_settings['data']['dynamic']['outcome']['folder'].format(**template_filled)
        os.makedirs(output_folder, exist_ok=True)
        
        ###################################################################
        #read radar data
        logging.info("--> Preparing radar data...")
        logging.info("---> Reading radar maps...")
        file_now_radar = os.path.join(data_settings['data']['dynamic']['input_rad']["folder"],data_settings['data']['dynamic']['input_rad']["filename"]).format(**template_filled)
        try:
            da_radar = rioxarray.open_rasterio(file_now_radar)
            pesi_op_radar = pesi_radar.copy()
        except:
            logging.warning("---> WARNING! Radar map not exist" + time_now.strftime("%Y-%m-%d %H:%M"))      
            da_radar = griglia.copy()
            da_radar.values=np.zeros(griglia.shape)
            pesi_op_radar = 0

        # rigriglio
        da_radar_grid = da_radar.reindex({'band': np.ones(1), 'x': griglia.x.values, 'y': griglia.y.values},
                                                     method='nearest')

        ###################################################################
        #read saellite data
        logging.info("--> Preparing satellite data...")
        logging.info("---> Reading satellite maps...")
        
        # satelliti risoluzione 1h e prodotto finale risoluzione 30min
        if (data_settings['data']['dynamic']['input_sat']["resolution"]=="60min" and data_settings['data']['dynamic']['time']['time_frequency']=="30min"):
            # istatne XX:30                
            if time_now.minute==30: 
                time_now_new = time_now + pd.Timedelta(data_settings['data']['dynamic']['time']['time_frequency'])
                template_filled_new = fill_template_time(data_settings["template"], time_now_new)        
                file_now_sat = os.path.join(data_settings['data']['dynamic']['input_sat']["folder"],data_settings['data']['dynamic']['input_sat']["filename"]).format(**template_filled_new)     
                
            # istatne XX:00
            else:
                file_now_sat = os.path.join(data_settings['data']['dynamic']['input_sat']["folder"],data_settings['data']['dynamic']['input_sat']["filename"]).format(**template_filled)
                  
            try:
                da_satellite = rioxarray.open_rasterio(file_now_sat)
                da_satellite = da_satellite/2
                pesi_op_satellite = pesi_satellite.copy()
            except:
                logging.warning("---> WARNING! Satellite map not exist" + time_now.strftime("%Y-%m-%d %H:%M") + " : " + file_now_sat)
                logging.info("FILE : " + file_now_sat)
                da_satellite = griglia.copy()
                da_satellite.values=np.zeros(griglia.shape)
                pesi_op_satellite = 0
        
        # satelliti risoluzione 30min e prodotto finale risoluzione 30min                                        
        elif data_settings['data']['dynamic']['input_sat']["resolution"]=="30min" and data_settings['data']['dynamic']['time']['time_frequency']=="30min":
            file_now_sat = os.path.join(data_settings['data']['dynamic']['input_sat']["folder"],data_settings['data']['dynamic']['input_sat']["filename"]).format(**template_filled)         
            try:
                da_satellite = rioxarray.open_rasterio(file_now_sat)
                da_satellite = da_satellite/2
                pesi_op_satellite = pesi_satellite.copy()
            except:
                logging.warning("---> WARNING! Satellite map not exist" + time_now.strftime("%Y-%m-%d %H:%M") + " : " + file_now_sat)
                da_satellite = griglia.copy()
                da_satellite.values=np.zeros(griglia.shape)
                pesi_op_satellite = 0

        else:
            logging.error("--> ERROR! Only satellite products with 30min or 60min temporal resolution are supported!")
            raise NotImplementedError
        
        '''                                
        # satelliti risoluzione ??min e prodotto finale risoluzione ??min
        else:
                
                file_now_sat = os.path.join(data_settings['data']['dynamic']['input_sat']["folder"],data_settings['data']['dynamic']['input_sat']["filename"]).format(**template_filled)  
                try:
                        da_satellite = rioxarray.open_rasterio(file_now_sat)
                except rio.errors.RasterioIOError:
                        logging.warning("---> WARNING! Satellite map not exist" + time_now.strftime("%Y-%m-%d %H:%M") 
                        missing_map = True
                        da_satellite = np.zeros(1470*1440)#################################
        '''
                                        
        # rigriglio
        da_satellite_grid = da_satellite.reindex({'band':np.ones(1),'x':griglia.x.values,'y':griglia.y.values},method='nearest')


        ###################################################################
        # DATAFUSION                                       
                                           
                                        
        remote = griglia.copy()
        #remote.values[0]=da_satellite_grid[0]*pesi_sat+da_radar[0]*pesi_rad                           
        remote.values[0].astype(np.float32)
        remote.values=np.float32(da_satellite_grid*pesi_op_satellite+da_radar_grid*pesi_op_radar)
                                     
        remote.rio.to_raster(os.path.join(data_settings['data']['dynamic']['outcome']["folder"],data_settings['data']['dynamic']['outcome']["filename"]).format(**template_filled), compress="DEFLATE", dtype="float32")
        logging.info("--> Preparing remote data...DONE")
                                        
       
 

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
