import os
import datetime

OUTPUT_DIRECTORY = '/okyanus/users/milicak/python_tools/Analysis/opendrift/Analysis'
OUTPUT_FILENAME = 'data_BS_forecast.nc'

service_id = 'BLKSEA_ANALYSISFORECAST_PHY_007_001-TDS'
product_id = 'bs-cmcc-cur-an-fc-h'
lon_min = '27.370073318481445'
lon_max = '41.96255111694336'
lat_min = '40.86015319824219'
lat_max = '46.80463409423828'

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
depth_min = '2.5'
depth_max = '2.5012000000000003'

cmnd = 'python -m motuclient --motu http://nrt.cmems-du.eu/motu-web/Motu --service-id ' + \
        service_id + ' --product-id ' + product_id + ' --longitude-min ' + lon_min + \
        ' --longitude-max ' + lon_max + ' --latitude-min ' + lat_min + \
        ' --latitude-max ' + lat_max + ' --date-min ' + date_str + \
        ' --date-max ' + date_end + ' --depth-min ' + depth_min + \
        ' --depth-max ' + depth_max + ' --variable uo --variable vo ' + \
        ' --out-dir ' + OUTPUT_DIRECTORY + ' --out-name ' + OUTPUT_FILENAME +\
        ' --user "milicak1" --pwd "mhmt@MI1067"'

os.system(cmnd)

