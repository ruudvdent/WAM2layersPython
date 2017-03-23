# -*- coding: utf-8 -*-
"""
Created on Thu Dec 22 08:07:59 2016

@author: vijverbe
"""

from ecmwfapi import ECMWFDataServer
     
server = ECMWFDataServer()


grid = "1.125/1.125"

#server.retrieve({
#    'dataset'   : "era5_test",
#    'stream'    : "oper/enda",   # 'oper' specifies the high resolution daily data, as opposed to monthly means, wave, eda edmm, etc.
#    'type'      : "an",    # We want instantaneous parameters, which are archived as type Analysis ('an') as opposed to forecast (fc)
#    'levtype'   : "sfc",   # Surface level, as opposed to pressure level (pl) or model level (ml)
#    'param'     : "26",   # For parameter codes see the ECMWF parameter database at http://apps.ecmwf.int/codes/grib/param-db
#    'grid'      : grid,     # The spatial resolution in ERA5 is 31 km globally on a Gaussian grid. Here we us lat/long with 0.25 degrees, which is approximately the equivalent of 31km.
#    'time'      : "12:00:00",  # ERA5 provides hourly analysis
#    'date'      : "2016-01-01/to/2016-01-05",
#    'target'    : "era5_test_2016-01-01to2016-01-05_hourly.grib" # Default output format is GRIB
# })
# 
server.retrieve({
    "class": "ei",
    "dataset": "era5_test",
    "date": "2016-01-01",
    "expver": "1",
    "grid": grid,
    "levtype": "sfc",
    "param": "26",
    "step": "0",
    "stream": "oper",
    "time": "12:00:00",
    "type": "an",
    "target": "lakemask.nc",
    'format' : "netcdf"
})
#%%
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import matplotlib as cm
import numpy as np
path = '/nobackup/users/vijverbe/Download_scripts/lakemask.nc'

lakemask = np.squeeze(Dataset(path, mode='r')['cl'][:,:,:])
latitude = Dataset(path, mode='r')['latitude'][:]
longitude = Dataset(path, mode='r')['longitude'][:]



#%% ################### GRIDCELL PLOT
from mpl_toolkits.basemap import Basemap



plt.figure() # nice projections are 'cea', 'cyl', 'eck4' and 'robin', however only 'cea' and 'cyl' allows for not showing the entire world
map = Basemap(projection='cyl',lat_0=0,lon_0=0,
              llcrnrlon=-180, llcrnrlat=-78,
              urcrnrlon=180, urcrnrlat=78)
map.drawcoastlines(linewidth=0.25)
map.drawcountries(linewidth=0.25)

x_shift,y_shift = np.meshgrid(longitude,latitude)
x = x_shift - 0.75
y = y_shift + 0.75

lol = map.pcolormesh(x,y,lakemask, cmap=cm.cm.coolwarm,     latlon=True, vmin=0,vmax=1)
map.colorbar(lol,location='right',pad="5%", label = '(-)')

plt.savefig('lakemask_ERA5.PNG', format='PNG', dpi=200)
