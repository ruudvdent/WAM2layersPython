#!/usr/bin/python
from ecmwfapi import ECMWFDataServer
import numpy as np

server = ECMWFDataServer()

#{
#    "url"   : "https://api.ecmwf.int/v1",
#    "key"   : "776c3a66a6aac2386985b1b1479e6f13",
#    "email" : "r.j.vanderent@tudelft.nl"
#}

for y in range(2010,2011):
    year = np.str(y)
    server.retrieve({
        'dataset' : "interim",
        'date'    : year+"0101/to/"+year+"1231",
        'time'    : "00:00:00/06:00:00/12:00:00/18:00:00",
        'stream'  : "oper",
        'grid'    : "1.5/1.5",
        'step'    : "0",
        'levtype' : "sfc",
        'type'    : "an",
        'class'   : "ei",
        'param'   : "136.128",
        #'area'    : "79.5/-180/-57/178.5",
        'format'  : "netcdf",
        'target'  : year+"-tcw.nc"
        })
        
for y in range(2010,2011):
    year = np.str(y)
    server.retrieve({
        'dataset' : "interim",
        'date'    : year+"0101/to/"+year+"1231",
        'time'    : "00:00:00/06:00:00/12:00:00/18:00:00",
        'stream'  : "oper",
        'grid'    : "1.5/1.5",
        'step'    : "0",
        'levtype' : "sfc",
        'type'    : "an",
        'class'   : "ei",
        'param'   : "137.128",
        #'area'    : "79.5/-180/-57/178.5",
        'format'  : "netcdf",
        'target'  : year+"-tcwv.nc"
        })

for y in range(2010,2011):
    year = np.str(y)
    server.retrieve({
        'dataset' : "interim",
        'date'    : year+"0101/to/"+year+"1231",
        'time'    : "00:00:00/06:00:00/12:00:00/18:00:00",
        'stream'  : "oper",
        'grid'    : "1.5/1.5",
        'step'    : "0",
        'levtype' : "sfc",
        'type'    : "an",
        'class'   : "ei",
        'param'   : "71.162",
        #'area'    : "79.5/-180/-57/178.5",
        'format'  : "netcdf",
        'target'  : year+"-ewvf.nc"
        })
        
for y in range(2010,2011):
    year = np.str(y)
    server.retrieve({
        'dataset' : "interim",
        'date'    : year+"0101/to/"+year+"1231",
        'time'    : "00:00:00/06:00:00/12:00:00/18:00:00",
        'stream'  : "oper",
        'grid'    : "1.5/1.5",
        'step'    : "0",
        'levtype' : "sfc",
        'type'    : "an",
        'class'   : "ei",
        'param'   : "88.162",
        #'area'    : "79.5/-180/-57/178.5",
        'format'  : "netcdf",
        'target'  : year+"-eclwf.nc"
        })
        
for y in range(2010,2011):
    year = np.str(y)
    server.retrieve({
        'dataset' : "interim",
        'date'    : year+"0101/to/"+year+"1231",
        'time'    : "00:00:00/06:00:00/12:00:00/18:00:00",
        'stream'  : "oper",
        'grid'    : "1.5/1.5",
        'step'    : "0",
        'levtype' : "sfc",
        'type'    : "an",
        'class'   : "ei",
        'param'   : "90.162",
        #'area'    : "79.5/-180/-57/178.5",
        'format'  : "netcdf",
        'target'  : year+"-ecfwf.nc"
        })
        
for y in range(2010,2011):
    year = np.str(y)
    server.retrieve({
        'dataset' : "interim",
        'date'    : year+"0101/to/"+year+"1231",
        'time'    : "00:00:00/06:00:00/12:00:00/18:00:00",
        'stream'  : "oper",
        'grid'    : "1.5/1.5",
        'step'    : "0",
        'levtype' : "sfc",
        'type'    : "an",
        'class'   : "ei",
        'param'   : "72.162",
        #'area'    : "79.5/-180/-57/178.5",
        'format'  : "netcdf",
        'target'  : year+"-nwvf.nc"
        })
        
for y in range(2010,2011):
    year = np.str(y)
    server.retrieve({
        'dataset' : "interim",
        'date'    : year+"0101/to/"+year+"1231",
        'time'    : "00:00:00/06:00:00/12:00:00/18:00:00",
        'stream'  : "oper",
        'grid'    : "1.5/1.5",
        'step'    : "0",
        'levtype' : "sfc",
        'type'    : "an",
        'class'   : "ei",
        'param'   : "89.162",
        #'area'    : "79.5/-180/-57/178.5",
        'format'  : "netcdf",
        'target'  : year+"-nclwf.nc"
        })
        
for y in range(2010,2011):
    year = np.str(y)
    server.retrieve({
        'dataset' : "interim",
        'date'    : year+"0101/to/"+year+"1231",
        'time'    : "00:00:00/06:00:00/12:00:00/18:00:00",
        'stream'  : "oper",
        'grid'    : "1.5/1.5",
        'step'    : "0",
        'levtype' : "sfc",
        'type'    : "an",
        'class'   : "ei",
        'param'   : "91.162",
        #'area'    : "79.5/-180/-57/178.5",
        'format'  : "netcdf",
        'target'  : year+"-ncfwf.nc"
        })

        
