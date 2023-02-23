import os
from os import path
import pandas as pd
import time
start = time.time()

# Input CMEMS User and Password (case sensitive)
# and desired directory to write data to.
# NOTE THAT YOUR PASSWORD CANNOT END WITH AN AMPERSAND - THIS WILL CAUSE THIS TO FAIL
USER = 'milicak1'
PASSWORD = 'mhmt@MI1067'
out_dir = '/Volumes/A1/workdir/milicak/datasets/glorys/'
new_name = "test"
# Set Lat/lon bounds
min_lon = str(-104)
max_lon = str(57)
min_lat = str(-30)
max_lat = str(81)


def name_from_day(day_string):
    return f'GLORYS_REANALYSIS_{day_string}.nc'


def get_days_in_year(year):
    all_days = pd.date_range(f'{year}-01-01', f'{year}-12-31')
    day_strings = [d.strftime('%Y-%m-%d') for d in all_days]
    return day_strings


def download_day(day_string):
    t1 = f'{day_string} 00:00:00'
    t2 = f'{day_string} 11:59:59'
    new_name = name_from_day(day_string)
    command = f'python -m motuclient --motu http://my.cmems-du.eu/motu-web/Motu --service-id GLOBAL_MULTIYEAR_PHY_001_030-TDS --product-id cmems_mod_glo_phy_my_0.083_P1D-m --longitude-min {min_lon} --longitude-max {max_lon} --latitude-min {min_lat} --latitude-max {max_lat} --date-min {t1} --date-max {t2} --depth-min 0.493 --depth-max 5727.918 --variable so --variable thetao --variable uo --variable siconc --variable vo --variable zos --out-dir {out_dir} --out-name {new_name} --user {USER} --pwd {PASSWORD}'
    os.system(command)


def download_year(year):
    for day in get_days_in_year(year):
        ntries = 0
        while not path.isfile(path.join(out_dir, name_from_day(day))):
            if ntries > 3:
                print(f'Failed to download {day}')
                break
            download_day(day)
            ntries += 1


if __name__ == '__main__':
    download_year(2017)
