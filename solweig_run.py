import fnmatch
import json
import os
from datetime import datetime, timezone, timedelta
from os import environ, path

import pandas as pd
import rasterio as rio
import requests
from dotenv import load_dotenv

import solweig_mrt as sol

# Use the .env file to get the API key for the OIKOLAB API
basedir = path.abspath(path.dirname(__file__))
load_dotenv(path.join(basedir, ".env"))

# General Config
OIKOLAB_KEY = environ.get("OIKOLAB_KEY")

DEFAULT_WIND_SPEED = 3.0
DEFAULT_TEMPERATURE = 30.0
DEFAULT_RELATIVE_HUMIDITY = 37.0
DEFAULT_SOLAR_RADIATION = 1360.0

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
    # dsm,vegdsm,dem,res,trans,svf, svfN, svfW, svfE, svfS, svfveg, svfNveg,svfEveg, svfSveg,svfWveg, svfaveg,
    # svfEaveg, svfSaveg, svfWaveg,svfNaveg,walls,dirwalls
    svf_walls_paths = os.path.join(os.getcwd(), '2020_svf_s/')
    files_list = [svf_walls_paths + files for files in os.listdir(svf_walls_paths)]

    DSM = rio.open(fnmatch.filter(files_list, '*updated_dsm.tif')[0])
    Vegdsm = rio.open(fnmatch.filter(files_list, '*updated_cdsm.tif')[0])
    Dem = rio.open(fnmatch.filter(files_list, '*updated_dem.tif')[0])
    Svf = rio.open(fnmatch.filter(files_list, '*updated_bu.tif')[0])
    SvfN = rio.open(fnmatch.filter(files_list, '*updated_N.tif')[0])
    SvfW = rio.open(fnmatch.filter(files_list, '*updated_W.tif')[0])
    SvfE = rio.open(fnmatch.filter(files_list, '*updated_E.tif')[0])
    SvfS = rio.open(fnmatch.filter(files_list, '*updated_S.tif')[0])
    Svfveg = rio.open(fnmatch.filter(files_list, '*updated_veg.tif')[0])
    SvfNveg = rio.open(fnmatch.filter(files_list, '*updated_Nveg.tif')[0])
    SvfEveg = rio.open(fnmatch.filter(files_list, '*updated_Eveg.tif')[0])
    SvfSveg = rio.open(fnmatch.filter(files_list, '*updated_Sveg.tif')[0])
    SvfWveg = rio.open(fnmatch.filter(files_list, '*updated_Wveg.tif')[0])
    Svfaveg = rio.open(fnmatch.filter(files_list, '*updated_aveg.tif')[0])
    SvfEaveg = rio.open(fnmatch.filter(files_list, '*updated_Eaveg.tif')[0])
    SvfSaveg = rio.open(fnmatch.filter(files_list, '*updated_Saveg.tif')[0])
    SvfWaveg = rio.open(fnmatch.filter(files_list, '*updated_Waveg.tif')[0])
    SvfNaveg = rio.open(fnmatch.filter(files_list, '*updated_Naveg.tif')[0])
    Walls = rio.open(fnmatch.filter(files_list, '*updated_walls.tif')[0])
    Dirwalls = rio.open(fnmatch.filter(files_list, '*updated_dirwalls.tif')[0])
    return DSM, Vegdsm, Dem, Svf, SvfN, SvfW, SvfE, SvfS, Svfveg, SvfNveg, SvfEveg, SvfSveg, SvfWveg, Svfaveg, SvfEaveg, SvfSaveg, SvfWaveg, SvfNaveg, Walls, Dirwalls

def close_files( DSM, Vegdsm, Dem, Svf, SvfN, SvfW, SvfE, SvfS, Svfveg, SvfNveg, SvfEveg, SvfSveg, SvfWveg, Svfaveg, SvfEaveg, SvfSaveg, SvfWaveg, SvfNaveg, Walls, Dirwalls):
    DSM.close()
    Vegdsm.close()
    Dem.close()
    Svf.close()
    SvfN.close()
    SvfW.close()
    SvfE.close()
    SvfS.close()
    Svfveg.close()
    SvfNveg.close()
    SvfEveg.close()
    SvfSveg.close()
    SvfWveg.close()
    Svfaveg.close()
    SvfEaveg.close()
    SvfSaveg.close()
    SvfWaveg.close()
    SvfNaveg.close()
    Walls.close()
    Dirwalls.close()

def run_solweig(time_to_run):
    # Geting data from API
    r = requests.get('https://api.oikolab.com/weather',
                     params={'param': ['wind_speed', 'temperature', 'relative_humidity', 'surface_solar_radiation'],
                             # The day plus or minus one day for buffer 
                             'start': time_to_run.date() - timedelta(days=1),
                             'end': time_to_run.date() + timedelta(days=1),
                             # Tempe, AZ Location
                             'lat': 33.29,
                             'lon': -112.42,
                             'api-key': OIKOLAB_KEY}
                     )
    datetime = time_to_run.replace(tzinfo=None)

    if r.status_code != 200:
        print(f"Error: {r.status_code} {r.json()}", flush=True) 
        Ws = DEFAULT_WIND_SPEED
        Ta = DEFAULT_TEMPERATURE
        RH = DEFAULT_RELATIVE_HUMIDITY
        radG = DEFAULT_SOLAR_RADIATION
    else: 
        weather_data = json.loads(r.json()['data'])
        df = pd.DataFrame(index=pd.to_datetime(weather_data['index'],
                                            unit='s'),
                        data=weather_data['data'],
                        columns=weather_data['columns'])
         # Convert to no timezone to use as key
        # Convert to no timezone
        df.tz_localize(None)
        Ws = df.loc[datetime]['wind_speed (m/s)']
        Ta = df.loc[datetime]['temperature (degC)']
        # expected to be a percentage between 0% and 100%
        RH = df.loc[datetime]['relative_humidity (0-1)'] * 100
        radG = df.loc[datetime]['surface_solar_radiation (W/m^2)']
   
    # open files
    DSM, Vegdsm, Dem, Svf, SvfN, SvfW, SvfE, SvfS, Svfveg, SvfNveg, SvfEveg, SvfSveg, SvfWveg, Svfaveg, SvfEaveg, SvfSaveg, SvfWaveg, SvfNaveg, Walls, Dirwalls = open_files()

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

    # Close files
    close_files(DSM, Vegdsm, Dem, Svf, SvfN, SvfW, SvfE, SvfS, Svfveg, SvfNveg, SvfEveg, SvfSveg, SvfWveg, Svfaveg, SvfEaveg, SvfSaveg, SvfWaveg, SvfNaveg, Walls, Dirwalls)
   
    location = {'latitude': 33.29, 'longitude': -112.42}

    year, month, day, doy, hour, minu = convert_datetime(time_to_run)
    tzone = 'UTC'

    res = 1
    trans = 3

    if DSM.bounds == Dem.bounds == Vegdsm.bounds:  # simple sanity check here to make sure rasters are aligned
        mrt = sol.Solweig_2021a_calc(dsm, vegdsm, dem, res, trans, svf, svfN, svfW, svfE, svfS, svfveg, svfNveg,
                                     svfEveg,
                                     svfSveg, svfWveg, svfaveg, svfEaveg, svfSaveg, svfWaveg, svfNaveg, walls, dirwalls,
                                     location, tzone, year, month, day, doy, hour, minu, Ws, Ta, RH,
                                     radG, Twater=15.0, ani=0, cyl=1, usevegdem=1, onlyglobal=1, elvis=0, landcover=0,
                                     lc_grid=None)

    root = 'output/' + time_to_run.tz_convert(None).strftime('%Y-%m-%d-%H00')

    rt = rio.open(root + '_mrt.tif', 'w', driver='GTiff', height=mrt.shape[0], width=mrt.shape[1], count=1,
                  crs=DSM.crs, transform=DSM.transform, dtype=mrt.dtype, compress='lzw')
    rt.write(mrt, 1)
    rt.close

if __name__ == '__main__':
    today = datetime.now(timezone.utc)
    print(f"Running SOLWEIG for ${today}")
    # gets the timestamp but zeros out the minutes, seconds, and microseconds
    today_ts = pd.Timestamp(today).replace(minute=0, second=0, microsecond=0)
    run_solweig(today_ts)
