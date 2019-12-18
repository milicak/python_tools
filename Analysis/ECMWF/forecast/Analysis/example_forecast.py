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
        "class": "od",
        "date": "2019-11-29",
        "expver": "1",
        "grid": "0.125/0.125",
        "levtype": "sfc",
        "param": param,
        "step": "0/1/2/3/4/5/6/7/8/9/10/11/12/13/14/15/16/17/18/19/20/21/22/23/24/25/26/27/28/29/30/31/32/33/34/35/36/37/38/39/40/41/42/43/44/45/46/47/48/49/50/51/52/53/54/55/56/57/58/59/60/61/62/63/64/65/66/67/68/69/70/71/72/73/74/75/76/77/78/79/80/81/82/83/84/85/86/87/88/89/90/93/96/99/102/105/108/111/114/117/120/123/126/129/132/135/138/141/144/150/156/162/168/174/180/186/192/198/204/210/216/222/228/234/240",
        "stream": "oper",
        "format": "netcdf",
        "time": "00:00:00",
        "type": "fc",
    },
        "target"+"2019-11-29"+param_string + ".nc",
    )


