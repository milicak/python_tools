import sys
sys.path.append('/users_home/opa/mi19918/python_models/opendrift/')
from datetime import timedelta
import numpy as np
from opendrift.readers.unstructured import shyfem
from opendrift.models.oceandrift import OceanDrift

o = OceanDrift(loglevel=20)  # Set loglevel to 0 for debug information
o.set_config('general:coastline_action', 'previous')
o.set_config('drift:advection_scheme', 'runge-kutta4')

shyfem = shyfem.Reader('/work/opa/mi19918/Projects/uTSS_SHYFEM/work/out_0729/dnm2.nc')
o.add_reader(shyfem)
print(shyfem)


# Seed elements at defined positions, depth and time
# number of floats
N = 1000
z = -15*np.random.uniform(0, 1, N)
# radius is in meter
# o.seed_elements(lon=29.0, lat=40.97, radius=1000, number=N,
                # z=z, time=shyfem.start_time)

o.seed_elements(lon=29.0, lat=40.97, radius=1000, number=N,
# o.seed_elements(lon=29.0, lat=40.8, radius=1000, number=N,
                time=shyfem.start_time)

# o.run(time_step=1800, duration=timedelta(hours=24))
# o.run(time_step=1800, duration=timedelta(hours=24),outfile='/work/opa/mi19918/Projects/uTSS_SHYFEM/work/out_0729/opendrift.nc')
o.run(time_step=1800, duration=timedelta(hours=72),outfile='/work/opa/mi19918/Projects/uTSS_SHYFEM/work/out_0729/opendrift.nc')
# o.animation(color='z', markersize=5, corners=[26.0, 30.0, 40.0, 41.5], filename='/users_home/opa/mi19918/shyfem_floats.gif')



o2 = OceanDrift(loglevel=20)  # Set loglevel to 0 for debug information
o2.set_config('general:coastline_action', 'previous')
o2.set_config('drift:advection_scheme', 'runge-kutta4')
o2.add_reader(shyfem)
o2.set_config('drift:horizontal_diffusivity', 10)
o2.seed_elements(lon=29.0, lat=40.97, radius=1000, number=N,
# o2.seed_elements(lon=29.0, lat=40.8, radius=1000, number=N,
                time=shyfem.start_time)
o2.run(time_step=1800,duration=timedelta(hours=72),outfile='/work/opa/mi19918/Projects/uTSS_SHYFEM/work/out_0729/opendrift_hormix.nc')

o2.animation(compare=[o], legend=['Diffusion' , 'No diffusion'],
             legend_loc='upper center', 
             corners=[26.0, 30.0, 40.0, 41.5],fast=True,
            filename='/users_home/opa/mi19918/shyfem_floats2.gif')
