#!/usr/bin/env python
from ecmwfapi import ECMWFDataServer
from ecmwfapi import *

# server = ECMWFDataServer()

# server = ECMWFDataServer(url="https://api.ecmwf.int/v1",key="d72fc7d368a78b613af37d3c471cc7ab",email="milicak@itu.edu.tr")
server = ECMWFService("mars", url="https://api.ecmwf.int/v1",key="d72fc7d368a78b613af37d3c471cc7ab",email="milicak@itu.edu.tr")

param_u, param_v, param_t = "165.128", "166.128", "167.128"

for param_string, param in zip(["_u", "_v", "_t"],
                               [param_u, param_v, param_t]):

    server.execute({
        "class": "tr",
        "date": "2019-11-29",
        "expver": "9001",
        "levtype": "sfc",
        "param": param,
        "step": "0/3/6/9/12",
        "stream": "enfo",
        "format": "netcdf",
        "time": "00:00:00",
        "type": "cf",
    },
        "2019-11-29"+param_string + "IFS_Turkey.nc",
    )


