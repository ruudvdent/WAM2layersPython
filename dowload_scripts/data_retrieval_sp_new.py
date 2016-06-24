#!/usr/bin/python
from ecmwfapi import ECMWFDataServer
import numpy as np

server = ECMWFDataServer()

# Lines below are hardcopied into api.py
#{
#    URL   : "https://api.ecmwf.int/v1",
#    KEY   : "776c3a66a6aac2386985b1b1479e6f13",
#    EMAIL : "r.j.vanderent@tudelft.nl"
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
        'param'   : "134.128",
        # 'area'    : "79.5/-180/-57/178.5",
        'format'  : "netcdf",
        'target'  : year+"-sp.nc"
        })