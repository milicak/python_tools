import os
import datetime

OUTPUT_DIRECTORY = '/okyanus/users/milicak/python_tools/Analysis/opendrift/Analysis'
OUTPUT_FILENAME = 'data_MS_forecast.nc'

service_id = 'MEDSEA_ANALYSISFORECAST_PHY_006_013-TDS'
product_id = 'med-cmcc-cur-an-fc-h'
lon_min = '-17.2917'
lon_max = '36.2917'
lat_min = '30.1875'
lat_max = '45.9792'

# now = datetime.datetime.now()
# dd/mm/YY H:M:S
# dt_string = now.strftime("%d/%m/%Y %H:%M:%S")

today = datetime.date.today()
# 5 day into future
endday = today+datetime.timedelta(days=5)
# YY/mm/dd
d1 = today.strftime("%Y-%m-%d")
d2 = endday.strftime("%Y-%m-%d")
quotes = '"'
date_str = quotes + d1 + ' 00:00:00' + quotes
date_end = quotes + d2 + ' 00:00:00' + quotes

# date_str = '"2022-06-24 00:00:00"'
# date_end = '"2022-06-29 00:00:00"'
depth_min = '1.0182'
depth_max = '1.0183'

vars_all = ' --variable uo --variable vo '

cmnd = 'python -m motuclient --motu http://nrt.cmems-du.eu/motu-web/Motu --service-id ' + \
        service_id + ' --product-id ' + product_id + ' --longitude-min ' + lon_min + \
        ' --longitude-max ' + lon_max + ' --latitude-min ' + lat_min + \
        ' --latitude-max ' + lat_max + ' --date-min ' + date_str + \
        ' --date-max ' + date_end + ' --depth-min ' + depth_min + \
        ' --depth-max ' + depth_max + vars_all + \
        ' --out-dir ' + OUTPUT_DIRECTORY + ' --out-name ' + OUTPUT_FILENAME +\
        ' --user "milicak1" --pwd "mhmt@MI1067"'

os.system(cmnd)

