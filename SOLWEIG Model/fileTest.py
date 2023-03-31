from solweig_mrt import *
import requests
import json
import pandas as pd
from datetime import datetime, timedelta
import rasterio as rio

import building_mask

# Get Tomorrow's Dates 
presentday = datetime.now()

tomorrow = presentday

tomorrow_ts = pd.Timestamp(tomorrow).replace(hour=0,minute=0,second=0,microsecond=0)

#Geting data from API
r = requests.get('https://api.oikolab.com/weather',
                 params={'param': ['wind_speed','temperature','relative_humidity','surface_solar_radiation'],
                            #Only Include Tomorrow's Data
                        'start': tomorrow_ts.date(),
                         'end': tomorrow_ts.date(),
                            #Tempe, AZ Location
                         'lat': 33.427204,
                         'lon': -111.939896,
                            #Plug in API Key From Discord (will put in .env later)
                         'api-key': 'f56043c50c174dea9628ec8f280e8a11'}
                 )
print("Response Status:", r)

#Processing Response into Panda DataFrame
weather_data = json.loads(r.json()['data'])
df = pd.DataFrame(index=pd.to_datetime(weather_data['index'],
                                       unit='s'),
                  data=weather_data['data'],
                  columns=weather_data['columns'])





# API and hardcoded data
# time in UTC (MST: UTC - 7:00)
datetime = '2023-01-29 18:00:00'
ws = df.loc[datetime]['wind_speed (m/s)']
tmp = df.loc[datetime]['temperature (degC)']
rH = df.loc[datetime]['relative_humidity (0-1)']*100
rad = df.loc[datetime]['surface_solar_radiation (W/m^2)']
loc = {'latitude': 33.427204,
        'longitude': -111.939896}

print(loc)

tzone = 'US/Mountain'
year = 2023
month = 1
day = 30
hour = 11
min = 0

res = 1
trans = 3
twater = 15.0
ani = 0
cyl = 1
usevegdem = 1
onlyglobal = 1
elvis = 0
landcover = 0
lc_grid = None



# dsm,vegdsm,dem,res,trans,svf, svfN, svfW, svfE, svfS, svfveg, svfNveg,svfEveg, svfSveg,svfWveg, svfaveg, svfEaveg, svfSaveg, svfWaveg,svfNaveg,walls,dirwalls
from PIL import Image
import numpy

Image.MAX_IMAGE_PIXELS = 4811370788


tif_files = ['tempe_campus_dsm_al.tif', 'tempe_campus_cdsm_al.tif','tempe_campus_dem_al.tif','svf_al.tif','svfN_al.tif','svfW_al.tif', 'svfE_al.tif', 'svfS_al.tif','svfveg_al.tif','svfNveg_al.tif','svfEveg_al.tif','svfSveg_al.tif','svfWveg_al.tif', 'svfaveg_al.tif', 'svfEveg_al.tif', 'svfSveg_al.tif', 'svfWveg_al.tif', 'svfNveg_al.tif','tempe_campus_walls_al.tif', 'tempe_campus_dirwalls_al.tif']
svfs = []
raw = []

for svf in tif_files:
    # im = Image.open(svf)
    im = rio.open(svf)
    raw.append(im)
    arr = im.read(1)
    svfs.append(arr)



result = Solweig_2021a_calc(dsm=svfs[0], vegdsm=svfs[1], dem=svfs[2], res=res, trans=trans, svf=svfs[3], svfN=svfs[4], svfW=svfs[5], svfE=svfs[6], svfS=svfs[7], svfveg=svfs[8], svfNveg=svfs[9], svfEveg=svfs[10], svfSveg=svfs[11], svfWveg=svfs[12], svfaveg=svfs[13], svfEaveg=svfs[14], svfSaveg=svfs[15], svfWaveg=svfs[16], svfNaveg=svfs[17], walls=svfs[18], dirwalls=svfs[19], location=loc, tzone=tzone, year=year, month=month, day=day, doy=day, hour=hour, minu=min, Ws=ws, Ta=tmp, RH=rH, radG=rad ,Twater=twater, ani=ani, cyl=cyl, usevegdem=usevegdem, onlyglobal=onlyglobal, elvis=elvis, landcover=landcover, lc_grid=lc_grid)
mrt = result['Tmrt']
out_path = 'mrt_output'
if hour < 10:
    root = '0'+str(hour)+'00'+'/'''+str(year)+str(month)+str(day)+'_'+str(hour)+'00'

else:
    root = str(hour)+'00'+'/'''+str(year)+str(month)+str(day)+'_'+str(hour)+'00'

rt = rio.open('mrt2202.tif', 'w', driver = 'GTiff', height = mrt.shape[0], width = mrt.shape[1], count = 1,
                      crs= raw[0].crs,transform=raw[0].transform,dtype=mrt.dtype)
rt.write(mrt,1)

rt.close()

#Not Exactly sure how root works but I think this is clsoe
maskshpfn = 'Maps/Tempe_buildings_WGS84.shp'
rasterfn = root + '.tif'
newrasterfn = 'masked_' + root + '.tif'

building_mask(maskshpfn, rasterfn, newrasterfn)