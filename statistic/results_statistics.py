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
import os, logging, netrc, json, time, warnings

import matplotlib.pyplot as plt
import pytz
import pandas as pd
import numpy as np
import xarray as xr
import datetime as dt
from argparse import ArgumentParser
from drops2 import coverages
from drops2.utils import DropsCredentials
from drops2 import sensors
import rioxarray
import rasterio as rio
import matplotlib.pyplot as pl


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
    drops_settings = data_settings["data"]["dynamic"]["gauge"]["drops2"]
    if not all([drops_settings['DropsUser'], drops_settings['DropsPwd']]):
        netrc_handle = netrc.netrc()
        try:
            drops_settings['DropsUser'], _, drops_settings['DropsPwd'] = netrc_handle.authenticators(drops_settings['DropsAddress'])
        except:
            logging.error('--> Valid netrc authentication file not found in home directory! Generate it or provide user and password in the settings!')
            raise FileNotFoundError('Verify that your .netrc file exists in the home directory and that it includes proper credentials!')
    drops_settings['frequency'] = data_settings['data']['dynamic']['time']['time_frequency']

    n_nearest = data_settings["data"]["static"]["number_nearest_cells_sampling"]
    
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

    
    
    logging.info("--> Preparing data...")
    griglia =  rioxarray.open_rasterio(data_settings['data']['static']['grid']) 
    lon_grid,lat_grid = np.meshgrid(griglia.x.values,griglia.y.values)
     
    logging.info("--> Preparing rain gauge data...")
    df_pluvio_anagrafica = dload_drops_gauges(steps_run[0], date_run, drops_settings,'false')
    '''
    df_pluvio_data = pd.DataFrame(index=['time'],columns=['name'])   
    df_sensor = pd.DataFrame(index=['name','lon','lat','value'], columns = df_pluvio.columns)   
    '''
    
    
    #inizializzo df 
    df_gauge_finale = pd.DataFrame(index=steps_run,columns=df_pluvio_anagrafica.columns.values)
    df_radar_finale = pd.DataFrame(index=steps_run,columns=df_pluvio_anagrafica.columns.values)
    df_satellite_finale = pd.DataFrame(index=steps_run,columns=df_pluvio_anagrafica.columns.values)
    
    pluvio_validi=0
    for time_now in steps_run:   
        # inizializzo vettori per le mezzore

        template_filled = fill_template_time(data_settings["template"], time_now)
        
        ###################################################################
        #read gauge data
        
        time_now_start = time_now - pd.Timedelta(data_settings['data']['dynamic']['time']['time_frequency'])
        logging.info("--> Preparing rain gauge data...")
        df_pluvio = dload_drops_gauges(time_now_start, time_now, drops_settings,'true')
        pluvio_validi=np.max([pluvio_validi,np.shape(df_pluvio)[1]]) 
        
        PointValueGauge=np.zeros(np.shape(df_pluvio)[1])
        PointValueRadar=np.zeros(np.shape(df_pluvio)[1])
        PointValueSatellite=np.zeros(np.shape(df_pluvio)[1])
        
        for index,col in enumerate(df_pluvio.columns):
            PointValueGauge[index]=df_pluvio[col].value
        
        df_gauge_finale.loc[time_now,df_pluvio.columns.values]=PointValueGauge
        
    
        
        
        ###################################################################
        #read radar data
        logging.info("--> Preparing radar data...")
        logging.info("---> Reading radar maps...")
        file_now_radar = os.path.join(data_settings['data']['dynamic']['input_rad']["folder"],data_settings['data']['dynamic']['input_rad']["filename"]).format(**template_filled)
        try:
            da_radar = rioxarray.open_rasterio(file_now_radar)
        except:
            logging.warning("---> WARNING! Radar map not exist" + time_now.strftime("%Y-%m-%d %H:%M"))      
            da_radar = griglia.copy()
            da_radar.values=np.zeros(griglia.shape)

        # valori su pluvio
        for index,col in enumerate(df_pluvio.columns):
            info = df_pluvio[col]
            '''          
            PointValueRadar[index] = da_radar.reindex({'band': np.ones(1), 'x': info.lon, 'y': info.lat},method='nearest')          
            '''
            abslat = np.abs(lat_grid-info.lat)
            abslon= np.abs(lon_grid-info.lon)
            #c = np.maximum(abslon,abslat)
            dist = np.sqrt(abslat**2+abslon**2)
            latlon_idx = dist.flatten().argsort()[:n_nearest]
            vec_neighbour = da_radar.values.flat[latlon_idx]
            PointValueRadar[index] = vec_neighbour[np.argmin(np.abs(vec_neighbour-info["value"]))]
            # latlon_idx = np.argmin(c)
            #PointValueRadar[index] = da_radar.values.flat[latlon_idx]
        df_radar_finale.loc[time_now,df_pluvio.columns.values]=PointValueRadar  
            
            
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
            except:
                logging.warning("---> WARNING! Satellite map not exist" + time_now.strftime("%Y-%m-%d %H:%M"))
                da_satellite = griglia.copy()
                da_satellite.values=np.zeros(griglia.shape)
        
        # satelliti risoluzione 30min e prodotto finale risoluzione 30min                                        
        elif data_settings['data']['dynamic']['input_sat']["resolution"]=="30min" and data_settings['data']['dynamic']['time']['time_frequency']=="30min":
            file_now_sat = os.path.join(data_settings['data']['dynamic']['input_sat']["folder"],data_settings['data']['dynamic']['input_sat']["filename"]).format(**template_filled)         
            try:
                da_satellite = rioxarray.open_rasterio(file_now_sat)
                da_satellite = da_satellite/2
            except:
                logging.warning("---> WARNING! Satellite map not exist" + time_now.strftime("%Y-%m-%d %H:%M"))
                da_satellite = griglia.copy()
                da_satellite.values=np.zeros(griglia.shape)

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
                                        
        #rigriglio
        da_satellite_grid = da_satellite.reindex({'band':np.ones(1),'x':griglia.x.values,'y':griglia.y.values},method='nearest')
        
        # valori su pluvio
        for index,col in enumerate(df_pluvio.columns):
            info = df_pluvio[col]
            '''          
            PointValueRadar[index] = da_radar.reindex({'band': np.ones(1), 'x': info.lon, 'y': info.lat},method='nearest')          
            '''
            abslat = np.abs(lat_grid-info.lat)
            abslon= np.abs(lon_grid-info.lon)
            dist = np.sqrt(abslat ** 2 + abslon ** 2)
            latlon_idx = dist.flatten().argsort()[:n_nearest]
            vec_neighbour = da_satellite_grid.values.flat[latlon_idx]
            PointValueSatellite[index] = vec_neighbour[np.argmin(np.abs(vec_neighbour - info["value"]))]
            #c = np.maximum(abslon,abslat)
            #latlon_idx = np.argmin(c)
            #PointValueSatellite[index] = da_satellite_grid.values.flat[latlon_idx]
       
        df_satellite_finale.loc[time_now,df_pluvio.columns.values]=PointValueSatellite

    df_satellite_finale = df_satellite_finale.fillna(0)
    df_radar_finale = df_radar_finale.fillna(0)
    per_asse=np.ceil((np.max([df_gauge_finale.max().max(),df_radar_finale.max().max(),df_satellite_finale.max().max()]))*1.3)
    
    fig,ax = pl.subplots(1,1)
    ax.set_aspect('equal')
    for index,step in enumerate(df_gauge_finale.index.values):
        ax.scatter(df_gauge_finale.iloc[index],df_radar_finale.iloc[index],color='blue')  
        ax.scatter(df_gauge_finale.iloc[index],df_satellite_finale.iloc[index],color='red')
        ax.scatter(df_radar_finale.iloc[index],df_satellite_finale.iloc[index],color='black')
    ax.set_aspect('equal')
    #ax.set_xlim(0,per_asse)
    #ax.set_ylim(0,per_asse)
    pl.savefig(os.path.join(output_folder, steps_run[-1].strftime("%Y%m%d%H%M") + "_cumulative_single.png"))
    pl.close()
    
    fig = pl.figure(figsize=(20, 10))
    ax = fig.add_subplot(221)
    for index,step in enumerate(df_gauge_finale.index.values):
        ax.scatter(df_gauge_finale.iloc[index],df_radar_finale.iloc[index],color='blue')
    ax.set_aspect('equal')
    ax.set_xlabel(['Gauge (' + str(pluvio_validi) + ') - 30mins accumulation' ])
    ax.set_ylabel('Radar - 30mins accumulation')
    ax.set_xlim(0,per_asse)
    ax.set_ylim(0,per_asse)
    ax = fig.add_subplot(222)
    for index,step in enumerate(df_gauge_finale.index.values):
        ax.scatter(df_gauge_finale.iloc[index],df_satellite_finale.iloc[index],color='blue')
    ax.set_aspect('equal')
    ax.set_xlabel('Gauge - 30mins accumulation')
    ax.set_ylabel('Satellite - 30mins accumulation')
    ax.set_xlim(0,per_asse)
    ax.set_ylim(0,per_asse)  
    ax = fig.add_subplot(223)
    for index,step in enumerate(df_gauge_finale.index.values):
        ax.scatter(df_radar_finale.iloc[index],df_satellite_finale.iloc[index],color='blue')
    ax.set_aspect('equal')
    ax.set_xlabel('Radar - 30mins accumulation')
    ax.set_ylabel('Satellite - 30mins accumulation')
    ax.set_xlim(0,per_asse)
    ax.set_ylim(0,per_asse)
    pl.savefig(os.path.join(output_folder, steps_run[-1].strftime("%Y%m%d%H%M") + "_cumulative_panels.png"))
    pl.close()
    
    
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

    
    
    # scatterplot
    
    #data
    
    
# -------------------------------------------------------------------------------------
# Function for downloading rain fauge data with drops2
def dload_drops_gauges(time_start, time_end,drops_settings,removeNAN):
    logging.info("---> Querying drops dataset...")
    DropsCredentials.set(drops_settings["DropsAddress"], drops_settings["DropsUser"], drops_settings["DropsPwd"])
    sensors_list_P = sensors.get_sensor_list(drops_settings["DropsSensor"],
                                             geo_win=(drops_settings['lon_left'], drops_settings['lat_bottom'], drops_settings['lon_right'], drops_settings['lat_top']),
                                             group=drops_settings["DropsGroup"])
    date_from = time_start.strftime("%Y%m%d%H%M")
    date_to = time_end.strftime("%Y%m%d%H%M")

    df_pluvio = sensors.get_sensor_data(drops_settings["DropsSensor"], sensors_list_P, date_from, date_to, as_pandas=True)
    if removeNAN=='true':
        df_pluvio = df_pluvio.dropna('columns', how='all')

    df_pluvio_now = pd.DataFrame(index=['name','lon','lat','value'], columns = df_pluvio.columns)

    for col in df_pluvio.columns:
        series = df_pluvio[col]
        series[series<0] = np.nan
        if removeNAN=='true':
            series = series.dropna('rows', how='all')

        df_pluvio_now.loc['lat', col] = [i for i in sensors_list_P.list if i.id == col][0].lat
        df_pluvio_now.loc['lon', col] = [i for i in sensors_list_P.list if i.id == col][0].lng
        df_pluvio_now.loc['name', col] = [i for i in sensors_list_P.list if i.id == col][0].name

        if (len(series>0)) and (date_to in series.keys()):
            df_pluvio_now.loc['value',col] = series.resample(drops_settings['frequency'], closed='right', label='right').agg(pd.DataFrame.sum)[date_to]
        else:
            df_pluvio_now.loc['value', col] = np.nan
    return df_pluvio_now
    '''if removeNAN=='true':
        return df_pluvio_now.dropna('columns', how='any')
    else:
        return df_pluvio_now'''
    
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