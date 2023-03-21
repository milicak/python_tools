import numpy as np
import os
import xarray as xr
import subprocess
from natsort import natsorted

class Regridding:
    def __init__(self,inpath,outpath,maskpath,variables,resolution):
        self.inpath = inpath
        self.outpath=outpath
        self.maskpath=maskpath
        self.variables=variables
        self.res=resolution

    def regriddingMain(self,area):
        os.chdir(self.inpath)
        self.regridOutput(area)

    def regridOutput(self, area):
        print('regridding output')
        vars = self.variables
        res = self.res
        box = (',').join([str(f) for f in area]).replace('[', '').replace(']','')
        mask = os.path.join(self.maskpath, 'mask%s_%s.nc' % (int(res * 100000), box))
        ncs = natsorted([os.path.join(self.inpath, f) for f in os.listdir(self.inpath) if f.endswith('ous.nc')])
        ncs = ncs[791:]

        if not os.path.exists(mask):
            print ("getting mask")
            print ('python {exe} rast mask {nc} {msk} --dx {res} --dy {res} --topology {topo} --box {bbox}'.format(
                    bbox=box, exe=regridderPath,
                    nc=ncs[0], msk=mask, res=res,
                    topo=topology))
            subprocess.call(
                'python {exe} rast mask {nc} {msk} --dx {res} --dy {res} --topology {topo} --box {bbox}'.format(
                    bbox=box, exe=regridderPath,
                    nc=ncs[0], msk=mask, res=res,
                    topo=topology), shell=True)
        # n = 0
        # startTs = 0
        n = 791
        startTs = 791
        for nc in ncs:
            n += len(xr.open_dataset(nc).time)

            for var in vars:
                outfile = '{outpath}/REG_{startTs}_{n}_0_{var}_{res}.nc'.format(outpath=self.outpath,startTs=startTs, n=n, var=var, res=res)

                cmd = 'python {exe} rast var {ncs} {msk}  {out} --var {var}'.format(
                    exe=regridderPath,
                    ncs=nc, msk=mask, out=outfile, var=var)

                subprocess.call(cmd, shell=True)
            startTs = n

# [minLon,minLat,maxLon,maxLat]
box=[28.25, 40.575, 29.35, 41.226]# bounding box
# box=[26.0, 40.0, 30.0, 41.5]# bounding box

inpath='/work/opa/mi19918/Projects/uTSS_SHYFEM/work/OUT'
# outpath='/work/opa/mi19918/Projects/uTSS_SHYFEM/work/daily'
outpath='/work/opa/mi19918/Projects/uTSS_SHYFEM/work/tmp3'
maskpath='/users_home/opa/mi19918//tools/regridUnst'
regridderPath='/users_home/opa/mi19918/python_models/regularize_unstructured/toolz/unstruct'
variables=['u_velocity','v_velocity','water_level']
variables=['water_level']
# resolution=0.01
resolution=0.001
topology='element_index' #  element_index for shyfem, tri for ww3
os.makedirs(maskpath,exist_ok=True)
os.makedirs(outpath,exist_ok=True)

Regridding(inpath,outpath,maskpath,variables,resolution).regriddingMain(box)
