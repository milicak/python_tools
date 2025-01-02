import cdsapi
import urllib

c = cdsapi.Client()

root_folder = '/okyanus/users/milicak/dataset/ERA5/global_daily/'

variables = [
             '2m_temperature',
             'total_precipitation'
         ]


# choose years to request (this should overlap with the `period_hist` chosen for the cm data)
# this is chosen shorter for demonstration purposes
years =  list(map(str, range(1992, 2011)))

# choose months to request
months = list(map(str, range(1, 13)))

# choose a variable (must be a valid ERA5 CDS API name)
variable = "total_precipitation"
# variable = "2m_temperature"

# choose a required statistic (valid names given in the application description above)
statistic = "daily_mean"

# choose an area (should be the same as above)
area = {"lat": [35, 45], "lon": [25, 40]}

# Loop over years and months
filenames_for_cleanup= []
for yr in years:
    print(f"----- Requesting year: {yr} -----")
    for mn in months:
        result = c.service(
            "tool.toolbox.orchestrator.workflow",
            params={
                 "realm": "user-apps",
                 "project": "app-c3s-daily-era5-statistics",
                 "version": "master",
                 "kwargs": {
                     "dataset": "reanalysis-era5-single-levels",
                     "product_type": "reanalysis",
                     "variable": variable,
                     "statistic": statistic,
                     "year": yr,
                     "month": mn,
                     "time_zone": "UTC+00:0",
                     "frequency": "1-hourly",
                     "grid": "0.1/0.1",
                     "area": area,
                     },
            "workflow_name": "application"
        })


        # filename = f"{root_folder}/era5_{variable}_{statistic}_{yr}_{mn}.nc"
        filename = f"{root_folder}/era5_{variable}_{statistic}_{yr}_{str(mn).zfill(2)}.nc"
        url = result[0]['location']

        # Download nc file
        urllib.request.urlretrieve(url, filename)
        # Append filename to list of filenames to cleanup
        filenames_for_cleanup.append(filename)

