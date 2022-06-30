"""
CREWS Caribbean - Import radar data
__date__ = '20220120'
__version__ = '1.0.0'
__author__ =
        'Lauro Rossi (lauro.rossi@cimagoundation.org,
        'Flavio Pignone (flavio.pignone@cimafoundation.org',
        'Andrea Libertino (andrea.libertino@cimafoundation.org',
__library__ = 'crews'
General command line:
### python crews_import_radar_data.py -time "YYYY-MM-DD HH:MM"
Version(s):
20220120 (1.0.0) --> Beta release for CREWS Caribbean
"""
# -------------------------------------------------------------------------------------

# -------------------------------------------------------------------------------------
# Complete library
from pykdtree.kdtree import KDTree
import pandas as pd
import wradlib as wrl
import matplotlib.pyplot as pl
import numpy as np
import warnings
import os, glob
import datetime
import json, logging
import datetime as dt
from argparse import ArgumentParser
import wradlib.georef as georef
import xarray as xr
import pickle
import rioxarray
import time as Time

from pyresample.geometry import SwathDefinition
from pyresample.kd_tree import resample_nearest


    #try:
    #    get_ipython().magic("matplotlib inline")
    #except:
    #    pl.ion()
    #'''

def main(): #main(path, in_dir, out_dir, plot, out_flag):

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
    proj_wgs84 = georef.epsg_to_osr(4326)

    # Set algorithm logging
    os.makedirs(data_settings['data']['log']['folder'], exist_ok=True)
    set_logging(logger_file=os.path.join(data_settings['data']['log']['folder'], data_settings['data']['log']['filename']))
    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # Get algorithm settings
    settings_file, alg_time = get_args()
    date_run = dt.datetime.strptime(alg_time, "%Y-%m-%d %H:%M")
    start_run = date_run - dt.timedelta(hours=data_settings['data']['dynamic']['time']['time_observed_period_h'])

    #date_range = pd.date_range(start_run, date_run, freq='1H')
    date_range = pd.date_range(start_run, date_run, freq=data_settings['data']['dynamic']['time']['radar_frequency'])
    # -------------------------------------------------------------------------------------
    # Loop through dates
    start_time = Time.time()
    for time_now in date_range:
        logging.info(" --> ###################################################")
        logging.info(" --> Computing time : " + time_now.strftime("%Y-%m-%d %H:%M"))
        template_filled = fill_template_time(data_settings["template"], time_now)

        folder_in_now = data_settings["data"]["dynamic"]["input"]["folder"].format(**template_filled)
        file_in_now_raw = data_settings["data"]["dynamic"]["input"]["filename"].format(**template_filled)

        folder_out_now = data_settings["data"]["dynamic"]["output"]["folder"].format(**template_filled)
        os.makedirs(folder_out_now, exist_ok=True)
        folder_ancillary_now = data_settings["data"]["dynamic"]["ancillary"]["folder"].format(**template_filled)
        os.makedirs(folder_ancillary_now, exist_ok=True)

        
        #Set variable "WRADLIB_DATA"
        os.environ["WRADLIB_DATA"] = folder_in_now
        
        try:
            file_in_now = glob.glob(os.path.join(folder_in_now,file_in_now_raw))[0]
            logging.info(" ---> File " + os.path.basename(file_in_now) + " found!")
            filename = wrl.util.get_wradlib_data_file(file_in_now)
            rbdict = wrl.io.read_rainbow(filename)
            presenceZ=1
     
        except:
            logging.error("ERROR! File " + file_in_now_raw + " not found")
            #raise FileNotFoundError
            presenceZ=0

        if presenceZ>0:
            
            #Conversion of the maps (to be revised)
            data_rain_bin = rbdict['product']['data']['radarpicture']['datamap']['data']
            step = (float(rbdict['product']['data']['radarpicture']['@max']) - float(rbdict['product']['data']['radarpicture']['@min']))/254
            data_rain = (data_rain_bin * step) + float(rbdict['product']['data']['radarpicture']['@min'])
            data_rain[data_rain<0]=0

            lon = rbdict['product']['data']['sensorinfo']['lon']
            lat = rbdict['product']['data']['sensorinfo']['lat']
            alt = rbdict['product']['data']['sensorinfo']['alt']

            ######### griglia proiezione e plot #########
            logging.info(" ---> Gridding and reprojecting")
            radar_location = (float(lon), float(lat), float(alt)) # (lon, lat, alt) in decimal degree and meters
            lat_max = float(rbdict['product']['data']['radarpicture']['projection']['@lat_ul'])
            lat_min = float(rbdict['product']['data']['radarpicture']['projection']['@lat_lr'])
            lon_max = float(rbdict['product']['data']['radarpicture']['projection']['@lon_lr'])
            lon_min = float(rbdict['product']['data']['radarpicture']['projection']['@lon_ul'])

            x = np.linspace(lon_min, lon_max, num=int(rbdict['product']['data']['radarpicture']['projection']['@size_x']))
            y = np.linspace(lat_min, lat_max, num=int(rbdict['product']['data']['radarpicture']['projection']['@size_y']))

            logging.info(" --> Save output and closing...")
            if data_settings["data"]["dynamic"]["output"]["format"] == "tif":
                file_out_now = folder_out_now + '/' + data_settings["data"]["dynamic"]["output"]["filename_tif"].format(**template_filled)
                #da_rain_bin_tif = xr.DataArray(np.flipud(data_rain_bin), dims=["y","x"], coords={"x":x,"y":y}).rio.write_crs("epsg:4326", inplace=True)
                da_rain_tif = xr.DataArray(np.flipud(data_rain), dims=["y","x"], coords={"x":x,"y":y}).rio.write_crs("epsg:4326", inplace=True)
                da_rain_tif.rio.to_raster(file_out_now, compress="DEFLATE", dtype="float32")
                #da_rain_bin_tif.rio.to_raster(file_out_now + '_bin', compress="DEFLATE", driver='GTiff', dtype="float32")
            elif data_settings["data"]["dynamic"]["output"]["format"] == "netcdf":
                file_out_now = folder_out_now + '/' + data_settings["data"]["dynamic"]["output"]["filename_nc"].format(**template_filled)
                xr.DataArray(np.flipud(data_rain), dims=["lat","lon"], coords={"lon":x,"lat":y}).to_netcdf(os.path.join(folder_out_now, file_out_now))
            else:
                logging.error("ERROR! Only tif and netcdf output types supported!")
                raise NotImplementedError

    # -------------------------------------------------------------------------------------
    # Info algorithm
    time_elapsed = round(Time.time() - start_time, 1)

    logging.info(' ')
    logging.info(' ==> ' + alg_name + ' (Version: ' + alg_version + ' Release_Date: ' + alg_release + ')')
    logging.info(' ==> TIME ELAPSED: ' + str(time_elapsed) + ' seconds')
    logging.info(' ==> ... END')
    logging.info(' ==> Bye, Bye')
    logging.info(' ============================================================================ ')
    # -------------------------------------------------------------------------------------

# -------------------------------------------------------------------------------------



'''controlli ulteriori per togliere rumore sono:
1) fare verifica sulla verticale direzione (z) per dBZ. Se supera valore soglia di XXX, allora si rimuove il valore definendolo clutter.
Il problema è che fare i controlli sulla verticale con le elevazioni che partono ad angoli diversi, ma soprattutto hanno aperture diverse. 
Forse la soluzione è grigliare tutto e farlo a questo punto.
           
2) verifica su variabile V (velocità radiale)
Se la velocità è compresa tre -0.2 e 0.2, abbianta al controllo precedente, si considerano i pixel come ground clutter.
Anche qui si può passare prima dalla grigliatura, anche se perde un po' di senso essendo velocità radiale, ma vabbè
           
Quindi ci sarebbe da aggiungere la lettura del *V.vol che è identica a *dBZ.vol stesse keys, senza fare la parte di process data.
Rigrigliare tutto per applicare 1 e 2
Fare solo ora la ZR trasnform
salvare prodotto finale (che potrà essere la seconda elevazione o combinazione delle varie)
'''
# -------------------------------------------------------------------------------------

# -------------------------------------------------------------------------------------
def plot_polar (azi,r,data,title):
    azi = np.radians (azi)
    #azi = np.radians (np.arange (0,484,1)/484*360)
    rmesh, thetamesh = np.meshgrid(r[0:20],azi)
    fig, ax = pl.subplots(subplot_kw={'projection': 'polar'})
    ax.set_title(title)
    ax.set_theta_zero_location("N")
    ax.set_theta_direction(-1)
    pm = ax.contourf(thetamesh,rmesh,data[:,0:20],10)
    pl.colorbar(pm)
    #ax.contourf(thetamesh,rmesh,data)
    # display the Polar plot
    pl.show()
# -------------------------------------------------------------------------------------

# -------------------------------------------------------------------------------------
def simple_plot_1map(data1,r,azi,title1):
    fig = pl.figure(figsize=(8, 8))
    ax, pm = wrl.vis.plot_ppi(data1, r=r, az=azi, fig=fig, proj='cg')
    ax.set_title(title1)
# -------------------------------------------------------------------------------------

# -------------------------------------------------------------------------------------
def simple_plot_2maps(data1, title1, data2, title2):
    fig = pl.figure(figsize=(12, 8))
    ax = fig.add_subplot(121)
    ax, pm = wrl.vis.plot_ppi(data1, ax=ax)
    ax.set_title(title1)
    ax = fig.add_subplot(122)
    ax, pm = wrl.vis.plot_ppi(data2, ax=ax)
    ax.set_title(title2)
# -------------------------------------------------------------------------------------

# -------------------------------------------------------------------------------------
def simple_plot_3maps(data1, title1, data2, title2, data3, title3):
    fig = pl.figure(figsize=(12, 4))
    ax = fig.add_subplot(131)
    ax, pm = wrl.vis.plot_ppi(data1, ax=ax)
    pm.set_clim(-20, 60)
    cb = pl.colorbar(pm, shrink=0.8)
    ax.set_title(title1)
    ax = fig.add_subplot(132)
    ax, pm = wrl.vis.plot_ppi(data2, ax=ax)
    pm.set_clim(-20, 60)
    cb = pl.colorbar(pm, shrink=0.8)
    ax.set_title(title2)
    ax = fig.add_subplot(133)
    ax, pm = wrl.vis.plot_ppi(data3, ax=ax)
    pm.set_clim(-20, 60)
    cb = pl.colorbar(pm, shrink=0.8)
    ax.set_title(title3)
# -------------------------------------------------------------------------------------

# -------------------------------------------------------------------------------------
def plot_beams(fig, data, mybeams, sub=111):
    ax = fig.add_subplot(sub)
    labelsize=13
    for beam in range(mybeams.start, mybeams.stop):
        pl.plot(data[beam], label="{0} deg".format(beam))
    pl.grid()
    pl.text(0.99, 0.88, "Reflectivity along beams",
            horizontalalignment='right',
            transform = ax.transAxes, fontsize="large")
    pl.xlabel("range (km)", fontsize="large")
    pl.ylabel("Reflectivity (dBZ)", fontsize="large")
    pl.legend(loc="upper left")
    ax.tick_params(axis='x', labelsize=labelsize)
    ax.tick_params(axis='y', labelsize=labelsize)
    pl.ylim(-20,80)
    pl.xlim(0,250)
# -------------------------------------------------------------------------------------

# -------------------------------------------------------------------------------------
def plot_pia(fig, pia, sub=111, title=None):
    ax = fig.add_subplot(sub)
    labelsize=13
    pl.plot(pia.T)
    pl.grid()
    pl.ylim(0,30)
    pl.ylabel("PIA (dB)", fontsize="large")
    pl.xlabel("range (km)", fontsize="large")
    pl.text(0.01, 0.88, title,
            transform = ax.transAxes, fontsize="large")
    ax.tick_params(axis='x', labelsize=labelsize)
    ax.tick_params(axis='y', labelsize=labelsize)
    pl.xlim(0,250)
# -------------------------------------------------------------------------------------

# -------------------------------------------------------------------------------------
def mod_kramer_attenuation(data):
    # Generalised Kraemer procedure: adding additional constraints
    # The function wradlib.atten.correct_attenuation_constrained() allows us to pass any kind of constraint function or lists of constraint functions via the argument constraints. The arguments of these functions are passed via a nested list as argument constraint_args. For example, Jacobi et al., 2016 suggested to constrain both the corrected reflectivity (by a maximum of 59 dBZ) and the resulting path-intgrated attenuation PIA (by a maximum of 20 dB):
    # data= data*2.5

    pia_mkraemer = wrl.atten.correct_attenuation_constrained(
        data,
        a_max=1.67e-4,
        a_min=2.33e-5,
        n_a=100,
        b_max=0.7,
        b_min=0.65,
        n_b=6,
        gate_length=1.,
        constraints=[wrl.atten.constraint_dbz,
                     wrl.atten.constraint_pia],
        constraint_args=[[59.0], [20.0]])

    return (pia_mkraemer)
# -------------------------------------------------------------------------------------

# -------------------------------------------------------------------------------------
def clmap_gabella(data):
# Clutter detection using the Gabella approach

    clmap = wrl.clutter.filter_gabella(data,
                                            wsize=5,
                                            thrsnorain=0.,
                                            tr1=6.,
                                            n_p=8,
                                            tr2=1.3)

    return clmap
        #Plot data with annotation
# -------------------------------------------------------------------------------------

# -------------------------------------------------------------------------------------
def plot_ppi (data, r, azi, sensorname, lon, lat, elev_angle, date, time, unit, path, out_dir, out_flag):
        fig = pl.figure(figsize=(8,8))
        cgax, pm = wrl.vis.plot_ppi(data, r=r, az=azi, fig=fig, proj='cg')
        pm.set_clim(-20, 60)
        title = '{0} ({1}E {2}N)\nPPI at {3}° {4} {5}\n'.format(sensorname, lon, lat,
                                                              elev_angle, date, time)
        caax = cgax.parasites[0]
        t = pl.title(title, fontsize=12)
        t.set_y(1.1)
        cbar = pl.gcf().colorbar(pm, pad=0.075)
        caax.set_xlabel('x_range [km]')
        caax.set_ylabel('y_range [km]')
        pl.text(1.0, 1.05, 'azimuth', transform=caax.transAxes, va='bottom',
                ha='right')
        cbar.set_label('reflectivity [' + unit + ']')
        pl.tight_layout()
        if out_flag== "show":
            pl.show()
        elif out_flag== "save":
            os.system("mkdir -p "+os.path.join(path,out_dir))
            pl.savefig(os.path.join(path,out_dir,"PPI"+elev_angle+date+time+".png"))
#        pl.tight_layout()
# -------------------------------------------------------------------------------------

# -------------------------------------------------------------------------------------
def fill_template_time(template_empty, time_now):
    template_filled = {}
    for key in template_empty.keys():
        template_filled[key] = time_now.strftime(template_empty[key])
    return template_filled
# -------------------------------------------------------------------------------------

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
def kilometers2degrees(kilometer, radius=6371):
    import math
    return kilometer / (2.0 * radius * math.pi / 360.0)

# ----------------------------------------------------------------------------
# =======================================================

if __name__ == '__main__':
    main()
