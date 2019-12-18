#!/usr/bin/env python
from ecmwfapi import ECMWFDataServer
from ecmwfapi import *

# server = ECMWFDataServer()

# server = ECMWFDataServer(url="https://api.ecmwf.int/v1",key="d72fc7d368a78b613af37d3c471cc7ab",email="milicak@itu.edu.tr")
server = ECMWFService("mars", url="https://api.ecmwf.int/v1",key="d72fc7d368a78b613af37d3c471cc7ab",email="milicak@itu.edu.tr")

param_r = "228.128"

date = '2019-12-02'

for param_string, param in zip(["_r"],
                               [param_r]):

    server.execute({
        "class": "tr",
        "date": date,
        "expver": "9001",
        "levtype": "sfc",
        "param": param,
        "step": "0/3/6/9/12",
        "stream": "enfo",
        "format": "netcdf",
        "time": "00:00:00",
        "type": "cf",
    },
        date+param_string + "IFS_Turkey.nc",
    )


