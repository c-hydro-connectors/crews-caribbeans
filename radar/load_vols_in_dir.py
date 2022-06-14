import wradlib as wrl
import matplotlib.pyplot as pl
import numpy as np
import warnings
import os
import datetime

warnings.filterwarnings('ignore')
    #try:
    #    get_ipython().magic("matplotlib inline")
    #except:
    #    pl.ion()
    #'''

def plot_polar (azi,r,data,title):
    azi = np.radians (azi)
    #azi = np.radians (np.arange (0,484,1)/484*360)
    rmesh, thetamesh = np.meshgrid(r[0:20],azi)
    fig, ax = pl.subplots(subplot_kw={'projection': 'polar'})
    ax.set_title(title)
    ax.set_theta_zero_location("N")
    ax.set_theta_direction(-1)
    ax.contourf(thetamesh,rmesh,data[:,0:20],10)
    #ax.contourf(thetamesh,rmesh,data)
    # display the Polar plot
    pl.show()

def simple_plot_1map(data1,r,azi,title1):
    fig = pl.figure(figsize=(8, 8))
    ax, pm = wrl.vis.plot_ppi(data1, r=r, az=azi, fig=fig, proj='cg')
    ax.set_title(title1)

def simple_plot_2maps(data1, title1, data2, title2):
    fig = pl.figure(figsize=(12, 8))
    ax = fig.add_subplot(121)
    ax, pm = wrl.vis.plot_ppi(data1, ax=ax)
    ax.set_title(title1)
    ax = fig.add_subplot(122)
    ax, pm = wrl.vis.plot_ppi(data2, ax=ax)
    ax.set_title(title2)

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

def load_vols_in_dir(path, in_dir, out_dir, plot, out_flag):

    #Set variable "WRADLIB_DATA"
    os.environ["WRADLIB_DATA"] = path
    #Load Rainbow file

    # get list of all files
    for root, dirs, files in os.walk(os.path.join(path,in_dir)):
         files.sort()
         #for each file in the folder in the in_dir folder
         for file in files:
            #consider only files that contains "dBZ"
            #if ("dBZ.vol" in file):
            if ("dBZ.vol" in file) and (file[10]=='0') and (file[11]=='0'): # only every hour
                filename = wrl.util.get_wradlib_data_file(os.path.join(path,in_dir,file))
                rbdict = wrl.io.read_rainbow(filename)

                #Get total scans
                tot_scans = len(rbdict['volume']['scan']['slice'])

                #for i in range (tot_scans):
                print (file )
                for i in range(1,2):

                    #Get azimuthal data
                    azi = rbdict['volume']['scan']['slice'][i]['slicedata']['rayinfo'][0]['data']
                    azidepth =  float(rbdict['volume']['scan']['slice'][i]['slicedata']['rayinfo'][0]['@depth'])
                    azirange = float(rbdict['volume']['scan']['slice'][i]['slicedata']['rayinfo'][0]['@rays'])
                    azirange = 360.0 #wrong parameter for 2hn and 3rd override
                    azi =  azi * azirange / 2 ** azidepth
                    #Get reflectivity data
                    elev_angle = str(rbdict['volume']['scan']['slice'][i]['posangle'])
                    data = rbdict['volume']['scan']['slice'][i]['slicedata']['rawdata']['data']
                    datadepth = float(rbdict['volume']['scan']['slice'][i]['slicedata']['rawdata']['@depth'])
                    datamin = float(rbdict['volume']['scan']['slice'][i]['slicedata']['rawdata']['@min'])
                    datamax = float(rbdict['volume']['scan']['slice'][i]['slicedata']['rawdata']['@max'])
                    data = datamin + data * (datamax - datamin) / 2 ** datadepth

                    #Get annotation data
                    unit = rbdict['volume']['scan']['slice'][i]['slicedata']['rawdata']['@type']
                    time = rbdict['volume']['scan']['slice'][i]['slicedata']['@time']
                    date = rbdict['volume']['scan']['slice'][i]['slicedata']['@date']
                    lon = rbdict['volume']['radarinfo']['@lon']
                    lat = rbdict['volume']['radarinfo']['@lat']
                    sensorname = rbdict['volume']['radarinfo']['name']

                    #simple_plot_1map(data, "Reflectivity PPI "+elev_angle+"\n"+date+time)

                    # Create range array
                    stoprange = float(rbdict['volume']['scan']['slice'][i]['stoprange'])

                    # Attention! stoprange for the 1st elevation is incorrect: need to be overridden
                    oDate = datetime.datetime.strptime(date, '%Y-%m-%d')
                    rangestep = float(rbdict['volume']['scan']['slice'][i]['rangestep'])

                    if oDate.strftime('%Y') == '2016':
                        if i == 0: stoprange = 399
                        if i == 1: stoprange = 249
                        if i == 2: stoprange = 249
                        if i == 3: stoprange = 149
                        if i == 4: stoprange = 50
                    elif oDate.strftime('%Y') == '2021':
                        if i == 0: stoprange = 398.5

                    r = np.arange(0, stoprange, rangestep)

                    print (elev_angle, date)

                    if plot:
                        data[data<10]=-20
                        plot_polar(azi,r,data,"PPI"+elev_angle+"_"+date+time)
                        #plot_ppi(data, r, azi, sensorname, lon, lat, elev_angle, date, time, unit, path, out_dir,out_flag) # out_flag can be "save" or "plot"
                        #simple_plot_1map(data,r,azi,"PPI"+elev_angle+"_"+date+time+"\n")
                        #Clutter detection using the Gabella approach
                        #clmap = clmap_gabella(data)
                        #simple_plot_2maps(data, "Reflectivity", clmap, "Clutter map")
                        #mybeams = slice(55,59)

                        #solo per testare l'attenuazione
                        #data = data * 2.5

                        pia_mkraemer = mod_kramer_attenuation(data)
                        data_corrected = data+pia_mkraemer

                        #simple_plot_3maps(data, "Reflect. orig.", pia_mkraemer, "Attenuation", data_corrected, "Reflect. corr. (orig+att)")

                        ## this is to plot 1-D attenuation graph
                        #fig = pl.figure(figsize=(10, 6))
                        #plot_beams(fig, data, mybeams, 211)
                        #plot_pia(fig, pia_mkraemer[mybeams], 212,
                        #     "PIA according to modified Kraemer")

                    print()
# =======================================================

if __name__ == '__main__':

    # these parameters will be put in a .json file

    #root path
    path = "/home/silvestro/crews/data"
    #input folder
    in_dir = 'source/radar_vol_pol'
    #output folder
    out_dir= 'outcome/radar'

    load_vols_in_dir(path, in_dir, out_dir, plot=True, out_flag="save")
