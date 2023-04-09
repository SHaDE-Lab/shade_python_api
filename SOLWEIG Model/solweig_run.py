import fnmatch
import json
import os

import fiona
import pandas as pd
import rasterio
import rasterio as rio
import requests
from datetime import date
import solweig_mrt as sol
import building_mask
import rasterio.mask

def convert_datetime(datetime):
    # ex '2023-01-29 18:00:00' in UTC
    timetuple = datetime.timetuple()
    year = timetuple.tm_year
    month = timetuple.tm_mon
    day = timetuple.tm_mday
    doy = timetuple.tm_yday
    hour = timetuple.tm_hour
    minu = timetuple.tm_min
    return year, month, day, doy, hour, minu

def open_files():
    # dsm,vegdsm,dem,res,trans,svf, svfN, svfW, svfE, svfS, svfveg, svfNveg,svfEveg, svfSveg,svfWveg, svfaveg, svfEaveg, svfSaveg, svfWaveg,svfNaveg,walls,dirwalls
    svf_walls_paths = 'tempe_camp_svfs/'
    files_list = [svf_walls_paths + files for files in os.listdir(svf_walls_paths)]

    DSM = rio.open(fnmatch.filter(files_list, '*_dsm_al.tif')[0])
    Vegdsm = rio.open(fnmatch.filter(files_list, '*_cdsm_al.tif')[0])
    Dem = rio.open(fnmatch.filter(files_list, '*_dem_al.tif')[0])
    Svf = rio.open(fnmatch.filter(files_list, '*svfbu_al.tif')[0])
    SvfN = rio.open(fnmatch.filter(files_list, '*svfN_al.tif')[0])
    SvfW = rio.open(fnmatch.filter(files_list, '*svfW_al.tif')[0])
    SvfE = rio.open(fnmatch.filter(files_list, '*svfE_al.tif')[0])
    SvfS = rio.open(fnmatch.filter(files_list, '*svfS_al.tif')[0])
    Svfveg = rio.open(fnmatch.filter(files_list, '*svfveg_al.tif')[0])
    SvfNveg = rio.open(fnmatch.filter(files_list, '*svfNveg_al.tif')[0])
    SvfEveg = rio.open(fnmatch.filter(files_list, '*svfEveg_al.tif')[0])
    SvfSveg = rio.open(fnmatch.filter(files_list, '*svfSveg_al.tif')[0])
    SvfWveg = rio.open(fnmatch.filter(files_list, '*svfWveg_al.tif')[0])
    Svfaveg = rio.open(fnmatch.filter(files_list, '*svfaveg_al.tif')[0])
    SvfEaveg = rio.open(fnmatch.filter(files_list, '*svfEaveg_al.tif')[0])
    SvfSaveg = rio.open(fnmatch.filter(files_list, '*svfSaveg_al.tif')[0])
    SvfWaveg = rio.open(fnmatch.filter(files_list, '*svfWaveg_al.tif')[0])
    SvfNaveg = rio.open(fnmatch.filter(files_list, '*svfNaveg_al.tif')[0])
    Walls = rio.open(fnmatch.filter(files_list, '*_walls_al.tif')[0])
    Dirwalls = rio.open(fnmatch.filter(files_list, '*_dirwalls_al.tif')[0])
    return DSM, Vegdsm, Dem, Svf, SvfN, SvfW, SvfE, SvfS,Svfveg,SvfNveg, SvfEveg, SvfSveg, SvfWveg, Svfaveg, SvfEaveg, SvfSaveg, SvfWaveg, SvfNaveg, Walls, Dirwalls

#TODO - Modify this function when scheduling aws jobs to accept date and time
def run_solweig():
    # Get Tomorrow's Dates
    presentday = date.today()

    tomorrow = presentday
    # gets time at today at 21 UTC  2pm AZ time
    tomorrow_ts = pd.Timestamp(tomorrow).replace(hour=21, minute=0, second=0, microsecond=0)
    print(tomorrow_ts)
    # Geting data from API
    r = requests.get('https://api.oikolab.com/weather',
                     params={'param': ['wind_speed', 'temperature', 'relative_humidity', 'surface_solar_radiation'],
                             # Only Include Tomorrow's Data
                             'start': tomorrow_ts.date(),
                             'end': tomorrow_ts.date(),
                             # Tempe, AZ Location
                             'lat': 33.29,
                             'lon': -112.42,
                             # Plug in API Key From Discord (will put in .env later)
                             'api-key': 'f56043c50c174dea9628ec8f280e8a11'}
                     )
    print("Response Status:", r.status_code)
    print("Response Reason:", r.reason)
    print("Response Reason:", r.content)


    # Processing Response into Panda DataFrame
    weather_data = json.loads(r.json()['data'])
    df = pd.DataFrame(index=pd.to_datetime(weather_data['index'],
                                           unit='s'),
                      data=weather_data['data'],
                      columns=weather_data['columns'])

    #open files
    DSM, Vegdsm, Dem, Svf, SvfN, SvfW, SvfE, SvfS,Svfveg,SvfNveg, SvfEveg, SvfSveg, SvfWveg, Svfaveg, SvfEaveg, SvfSaveg, SvfWaveg, SvfNaveg, Walls, Dirwalls = open_files()

    # read the rasters
    dsm = DSM.read(1)
    vegdsm = Vegdsm.read(1)
    dem = Dem.read(1)
    svf = Svf.read(1)
    svfN = SvfN.read(1)
    svfW = SvfW.read(1)
    svfE = SvfE.read(1)
    svfS = SvfS.read(1)
    svfveg = Svfveg.read(1)
    svfNveg = SvfNveg.read(1)
    svfEveg = SvfEveg.read(1)
    svfSveg = SvfSveg.read(1)
    svfWveg = SvfWveg.read(1)
    svfaveg = Svfaveg.read(1)
    svfEaveg = SvfEaveg.read(1)
    svfSaveg = SvfSaveg.read(1)
    svfWaveg = SvfWaveg.read(1)
    svfNaveg = SvfNaveg.read(1)
    walls = Walls.read(1)
    dirwalls = Dirwalls.read(1)

    # API and hardcoded data
    # time in UTC (MST: UTC - 7:00)
    datetime = tomorrow_ts
    Ws = df.loc[datetime]['wind_speed (m/s)']
    Ta = df.loc[datetime]['temperature (degC)']
    RH = df.loc[datetime]['relative_humidity (0-1)'] * 100
    radG = df.loc[datetime]['surface_solar_radiation (W/m^2)']
    location = {'latitude': 33.29, 'longitude': -112.42}

    year, month, day, doy, hour, minu = convert_datetime(datetime)
    tzone = 'UTC'

    res = 1
    trans = 3

    if DSM.bounds == Dem.bounds == Vegdsm.bounds:  # simple sanity check here to make sure rasters are aligned
        rez = sol.Solweig_2021a_calc(dsm, vegdsm, dem, res, trans, svf, svfN, svfW, svfE, svfS, svfveg, svfNveg, svfEveg,
                                     svfSveg, svfWveg, svfaveg, svfEaveg, svfSaveg, svfWaveg, svfNaveg, walls, dirwalls,
                                     location, tzone, year, month, day, doy, hour, minu, Ws, Ta, RH,
                                     radG, Twater=15.0, ani=0, cyl=1, usevegdem=1, onlyglobal=1, elvis=0, landcover=0,
                                     lc_grid=None)
    mrt = rez['Tmrt']

    # if hour < 10:
    #     root = '0' + str(hour) + '00' + str(year) + '_' + str(month) + '_' + str(day) + '_0' + str(hour) + '00'
    # else:
    #     root = str(hour) + '00' +  str(year) + str(month) + str(day) + '_' + str(hour) + '00'

    root = 'output/' + str(year) + '-' + str(month) + '-' + str(day) + '-' + str(hour) + '00'

    rt = rio.open(root+'_mrt.tif', 'w', driver='GTiff', height=mrt.shape[0], width=mrt.shape[1], count=1,
                  crs=DSM.crs, transform=DSM.transform, dtype=mrt.dtype)
    rt.write(mrt, 1)
    rt.close

    # maskshpfn = 'Maps/Tempe_MaskedBuildingsRoads.shp'
    # rasterfn = root + '_mrt.tif'
    # newrasterfn = 'output/masked_' + str(year) + '-' + str(month) + '-' + str(day) + '-' + str(hour) + '00'
    #
    # building_mask.adjustWithMask(maskshpfn, rasterfn, newrasterfn)

run_solweig()
