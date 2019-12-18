python -m motuclient --motu http://nrt.cmems-du.eu/motu-web/Motu --service-id
BLKSEA_ANALYSIS_FORECAST_PHYS_007_001-DGF --product-id sv03-bs-cmcc-tem-an-fc-d
--date-min "2019-12-10 00:00:00" --date-max "2019-12-10 00:00:00" --out-dir
/okyanus/users/milicak/python_tools/Analysis/ECMWF/forecast/Analysis --out-name
dnm.nc --user milicak1 --pwd 'mhmt@MI1067'


python -m motuclient --motu http://nrt.cmems-du.eu/motu-web/Motu --service-id BLKSEA_ANALYSIS_FORECAST_PHYS_007_001-TDS --product-id sv03-bs-cmcc-sal-an-fc-h --longitude-min 27.370073318481445 --longitude-max 31.96255111694336 --latitude-min 40.86015319824219 --latitude-max 42.80463409423828 --date-min "2019-12-01 00:30:00" --date-max "2019-12-05 11:30:00" --depth-min 2.5 --depth-max 2.5012000000000003 --variable vosaline --out-dir /okyanus/users/milicak/python_tools/Analysis/ECMWF/forecast/Analysis --out-name dnm2.nc --user milicak1 --pwd 'mhmt@MI1067'
