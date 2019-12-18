#!/usr/bin/env python
from ecmwfapi import ECMWFDataServer
from ecmwfapi import *

# server = ECMWFDataServer()

# server = ECMWFDataServer(url="https://api.ecmwf.int/v1",key="d72fc7d368a78b613af37d3c471cc7ab",email="milicak@itu.edu.tr")
server = ECMWFService("mars", url="https://api.ecmwf.int/v1",key="d72fc7d368a78b613af37d3c471cc7ab",email="milicak@itu.edu.tr")

param_lw = "175.128"

for param_string, param in zip(["_lw"],
                               [param_lw]):

    server.execute({
        "class": "mc",
        "date": "2019-11-29",
        "expver": "0001",
        "levtype": "sfc",
        "param": param,
        "step": "0/3/6/9/12",
        "stream": "oper",
        # "format": "netcdf",
        "time": "00:00:00",
        "type": "fc",
    },
        "2019-11-29"+param_string + "IFS_Turkey.grib",
    )


