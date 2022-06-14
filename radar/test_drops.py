import drops2
from drops2 import sensors, coverages
from drops2.utils import DropsCredentials
import numpy as np
import pandas as pd
import datetime as dt

DropsCredentials.set("http://caribbean.mydewetra.cimafoundation.org/dds/rest", "admin", "geoserver")
sensor_classes = sensors.get_sensor_classes()
print(sensor_classes)
sensor_class='PLUVIOMETRO'
sensors_list_P = sensors.get_sensor_list(sensor_class, group="Dewetra%Caraibi")
date_from="202102260000"
date_to="202102270000"
df_pluvio = sensors.get_sensor_data(sensor_class, sensors_list_P, date_from, date_to, aggr_time=3600, as_pandas=True)
