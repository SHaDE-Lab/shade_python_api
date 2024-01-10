####################################################################################################
# This scripts calculates the mean radiant temperature for a set of digital surface model          #
# The main function was adopted from Fredrik Lindberg, fredrikl@gvc.gu.se                          #
# Goteborg Urban Climate Group                                                                     #
# Gothenburg University                                                                            #
# No landcover input is required and isotropy is assumed                                           #
#                                                                                                  #
#                                                                                                  #
####################################################################################################


# import numpy as np
# from solweig_cylindric_wedge_svfalfa import*
# from solweig_daylen import*
# from solweig_gvf_sunonsurface import *
# from solweig_k_rads import *
# from solweig_l_rads import * 
# from solweig_tswavedelay import *
# from solweig_wallshadow import *
# from clearness_index import *
# from diffusefraction import *
# from sunpos import altitude
# from sunpos import altitude
# from sunpos import azimuth
# from sunpos import zenith
# from sunpos import *
# import rasterio as rio


def Solweig_2021a_calc(dsm,vegdsm,dem,res,trans,svf, svfN, svfW, svfE, svfS, svfveg, svfNveg,svfEveg, svfSveg,svfWveg, svfaveg, svfEaveg, svfSaveg, svfWaveg,svfNaveg,walls,dirwalls,location,tzone,year,month,day,doy,hour,minu,Ws,Ta, RH, radG,Twater,ani,cyl,usevegdem,onlyglobal,elvis,landcover,lc_grid):
    
    from sunpos import altitude
    from sunpos import azimuth
    from sunpos import zenith
    from sunpos import sun_distance
    from sunpos import maxalt
    from sunpos import timestamp
    from clearness_index import clearnessindex
    from diffusefraction import diffusefraction
    from solweig_cylindric_wedge_svfalfa import cylindric_wedge
    from solweig_daylen import daylen
    from solweig_gvf_sunonsurface import sunonsurface_2018a,gvf_2018a
    from solweig_k_rads import Kvikt_veg,Kup_veg_2015a,Kside_veg_v2019a
    from solweig_l_rads import Lvikt_veg,Lside_veg_v2015a
    from solweig_tswavedelay import TsWaveDelay_2015a
    from solweig_wallshadow import shadowingfunction_wallheight_23
    #from umep_shadow_working import shadow_withveg_wall
    #from umep_shadow import shadow_withveg_only
    #from Solweig_2015a_metdata_noload import Solweig_2015a_metdata_noload
    import numpy as np
    import math
    import datetime
    import calendar
    import pytz
    import pysolar.solar as ps
    #from __future__ import division

    # albedo_b =None, absK=None, absL=None,ewall=None, Fside=None, Fup=None, Fcyl=None,, 
    #, dectime, altmax=None, dirwalls=None, walls=None, cyl=None
    # TgK=None, Tstart=None, alb_grid=None, emis_grid=None, TgK_wall=None, Tstart_wall=None,TmaxLST=None
    #TmaxLST_wall=None, first=None, second=None, svfbuveg=None,firstdaytime, timeadd, timestepdec=None, Tgmap1=None, 
    #Tgmap1E=None, Tgmap1S=None, Tgmap1W=None, Tgmap1N=None, CI=None, TgOut1, diffsh=None, ani=None
    #radD=None, radI=None altitude, azimuth, zen, jday, 
    
    

    # Variable definitions:
    # dsm = digital surface model
    # scale = height to pixel size (2m pixel gives scale = 0.5)
    # header = ESRI Ascii Grid header
    # sizey,sizex = no. of pixels in x and y
    # svf,svfN,svfW,svfE,svfS = SVFs for building and ground
    # svfveg,svfNveg,svfEveg,svfSveg,svfWveg = Veg SVFs blocking sky
    # svfaveg,svfEaveg,svfSaveg,svfWaveg,svfNaveg = Veg SVFs blocking buildings
    # vegdem = Vegetation canopy DSM
    # vegdem2 = Vegetation trunk zone DSM
    # albedo_b = buildings
    # absK = human absorption coefficient for shortwave radiation
    # absL = human absorption coefficient for longwave radiation
    # ewall = Emissivity of building walls
    # Fside = The angular factors between a person and the surrounding surfaces
    # Fup = The angular factors between a person and the surrounding surfaces
    # altitude = Sun altitude (degree)
    # azimuth = Sun azimuth (degree)
    # zen = Sun zenith angle (radians)
    # jday = day of year
    # usevegdem = use vegetation scheme
    # onlyglobal = calculate dir and diff from global
    # buildings = Boolena grid to identify building pixels
    # location = geographic location
    # height = height of measurements point
    # psi = 1 - Transmissivity of shortwave through vegetation
    # landcover = use landcover scheme !!!NEW IN 2015a!!!
    # sensorheight = Sensorheight of wind sensor
    # lc_grid = grid with landcoverclasses
    # lc_class = table with landcover properties
    # dectime = decimal time
    # altmax = maximum sun altitude
    # dirwalls = aspect of walls
    # walls = one pixel row outside building footprints
    # cyl = consider man as cylinder instead of cude
	
	
    tmp = svf + svfveg - 1.
    tmp[tmp < 0.] = 0.
    # %matlab crazyness around 0
    svfalfa = np.arcsin(np.exp((np.log((1. - tmp)) / 2.)))
	
    lat = location['latitude']
    lon = location['longitude']
	
	
	
    year = year
    month = month
    day = day
  
    hour = hour
    minu = minu
    time = timestamp(year,month,day,hour,minu,tzone)
    Ta = Ta
    RH = RH
    radG = radG
	
    alt= altitude(lat,lon,time)  #90. - sun['zenith'] 
    azmt= azimuth (lat,lon,time)#sun['azimuth']
    zen=  zenith (lat,lon,time)* (np.pi/180.)
	
	

    Twater = Twater
    Ws = Ws
    
    # # metdata = np.zeros((1, 24)) - 999

    # # metdata[0, 0] = year
    # # metdata[0, 1] = doy
    # # metdata[0, 2] = hour
    # # metdata[0, 3] = minu
    # # metdata[0, 11] = Ta
    # # metdata[0, 10] = RH
    # # metdata[0, 14] = radG
    # # metdata[0, 21] = radD
    # # metdata[0, 22] = radI
    # # metdata[0, 9] = Ws
	
    # # YYYY, altitude, azimuth, zen, jday, leafon, dectime, altmax = Solweig_2015a_metdata_noload(metadata,lat,lon,tzone)
    
    # # #psi = leafon * self.trans
    # # #psi[leafon == 0] = 0.5
	
    # # P = -999.0
	# # # %Creating vectors from meteorological input
    # # DOY = metdata[:, 1]
    # # hours = metdata[:, 2]
    # # minu = metdata[:, 3]
    # # Ta = metdata[:, 11]
    # # RH = metdata[:, 10]
    # # radG = metdata[:, 14]
    # # radD = metdata[:, 21]
    # # radI = metdata[:, 22]
    # # P = metdata[:, 12]
    # # Ws = metdata[:, 9]
	
    #some constants here :: change before running if need be 
    albedo_b = 0.20
    albedo_g = 0.15
    ewall = 0.90
    eground = 0.95
    absK = 0.70
    absL = 0.95
	
    P = -999.0
    
	
    #timeconstants
    timestepdec = 0
    timeadd = 0.
    timeaddE = 0.
    timeaddS = 0.
    timeaddW = 0.
    timeaddN = 0.
    firstdaytime = 1.
	
    #constants for a standing person 
    Fside = 0.22 #change 
    Fup = 0.06
    Fcyl = 0.28
    height = 1.1
	
	#building footprints from DSM and DEM
	
    buildings = dsm - dem
    buildings[buildings < 2.] = 1.
    buildings[buildings >= 2.] = 0.
	
	
	###########
    sizex = dsm.shape[0]  # rows
    sizey = dsm.shape[1]  # cols
    rows = dsm.shape[0]
    cols = dsm.shape[1]

    trunkratio = 25/100.
    psi = trans/100.0
	
    vegdsm2= vegdsm * trunkratio
    vegmax = vegdsm.max()
    amaxvalue = dsm.max() - dsm.min()
    amaxvalue = np.maximum(amaxvalue, vegmax)

    # Elevation vegdsms if buildingDSM includes ground heights
    vegdem = vegdsm + dsm
    vegdem[vegdem == dsm] = 0
    vegdem2 = vegdsm2 + dsm
    vegdem2[vegdem2 == dsm] = 0

    # Bush separation
    bush = np.logical_not((vegdem2*vegdem))*vegdem

	
	# ###############
    # vegdsm2= vegdsm * trunkratio
	
	# # amaxvalue
    # vegmax = vegdsm.max()
    # amaxvalue = dsm.max() - dsm.min()
    # amaxvalue = np.maximum(amaxvalue, vegmax)

    # # Elevation vegdsms if buildingDSM includes ground heights
    # vegdem = vegdsm + dsm
    # vegdem[vegdem == dsm] = 0
    # vegdem2 = vegdsm2 + dsm       #modified to deal with the absence of trunk information
    # vegdem2[vegdem2 == dsm] = 0  #modified to deal with the absence of trunk information
	
	# # Bush separation
    # bush = np.logical_not((vegdem2*vegdem))*vegdem
	
    svfbuveg =  svf - (1. - svfveg) * (1. - psi) 
	
	# %Initialization of maps
    Knight = np.zeros((rows, cols))
    Tgmap1 = np.zeros((rows, cols))
    Tgmap1E = np.zeros((rows, cols))
    Tgmap1S = np.zeros((rows, cols))
    Tgmap1W = np.zeros((rows, cols))
    Tgmap1N = np.zeros((rows, cols))
    TgOut1 = np.zeros((rows, cols))
	
    TgK = Knight + 0.37
    Tstart = Knight - 3.41
    alb_grid = Knight + albedo_g
    emis_grid = Knight + eground
    TgK_wall = 0.37
    Tstart_wall = -3.41
    TmaxLST = 15.
    TmaxLST_wall = 15.
	
	# %Parameterisarion for Lup
    #if not height:
    height = 1.1
	
    diffsh = None #condiering isotropic sky

    # %Radiative surface influence, Rule of thumb by Schmid et al. (1990).
    first = np.round(height)
    if first == 0.:
      first = 1.
    second = np.round((height * 20.))
    
    #############wall calculations################
    #walls = findwalls (dsm,3.0)
    
    #dirwalls = filter1Goodwin_as_aspect_v3(walls,res,dsm)
    

    # # # Core program start # # #
    # Instrument offset in degrees
    t = 0.

    # Stefan Bolzmans Constant
    SBC = 5.67051e-8

    # Find sunrise decimal hour - new from 2014a
    _, _, _, SNUP = daylen(doy, location['latitude'])

    # Vapor pressure
    ea = 6.107 * 10 ** ((7.5 * Ta) / (237.3 + Ta)) * (RH / 100.)

    # Determination of clear - sky emissivity from Prata (1996)
    msteg = 46.5 * (ea / (Ta + 273.15))
    esky = (1 - (1 + msteg) * np.exp(-((1.2 + 3.0 * msteg) ** 0.5))) + elvis  # -0.04 old error from Jonsson et al.2006

    if alt > 0: # # # # # # DAYTIME # # # # # #
        # Clearness Index on Earth's surface after Crawford and Dunchon (1999) with a correction
        #  factor for low sun elevations after Lindberg et al.(2008)
        I0, CI, Kt, I0et, CIuncorr = clearnessindex(zen, doy, Ta, RH , radG, location, P)
        if (CI > 1) or (CI == np.inf):
            CI = 1

        # Estimation of radD and radI if not measured after Reindl et al.(1990)
        if onlyglobal == 1:
            I0, CI, Kt, I0et, CIuncorr = clearnessindex(zen, doy, Ta, RH , radG, location, P)
            if (CI > 1) or (CI == np.inf):
             CI = 1

        radI, radD = diffusefraction(radG, alt, Kt, Ta, RH)
        print('rad diffuse done')
        print(radI, radD)
        print(datetime.datetime.now())
		
		#metdata = np.zeros((1, 24)) - 999
        metdata =np.zeros((1,24)) -999
        metdata[0, 0] = year
        metdata[0, 1] = doy
        metdata[0, 2] = hour
        metdata[0, 3] = minu
        metdata[0, 11] = Ta
        metdata[0, 10] = RH
        metdata[0, 14] = radG
        metdata[0, 21] = radD
        metdata[0, 22] = radI
        metdata[0, 9] = Ws
	
        ##YYYY, altitude, azimuth, zen, jday, leafon, dectime, altmax = Solweig_2015a_metdata_noload(metdata,lat,lon,time) use when doy =221
        dectime= doy+hour / 24 + minu / (60*24.)
        altmax= maxalt(lat,doy,year) #80.15 # to be calculated 
        print('rad dectime done')
        print(datetime.datetime.now())
    
    #psi = leafon * self.trans
    #psi[leafon == 0] = 0.5
	
	# %Creating vectors from meteorological input
        DOY = metdata[:, 1]
        hours = metdata[:, 2]
        minu = metdata[:, 3]
        Ta = metdata[:, 11]
        RH = metdata[:, 10]
        radG = metdata[:, 14]
        radD = metdata[:, 21]
        radI = metdata[:, 22]
        P = metdata[:, 12]
        Ws = metdata[:, 9]

        # Diffuse Radiation
        # Anisotropic Diffuse Radiation after Perez et al. 1993
        if ani == 1:
            patchchoice = 1
            zenDeg = zen*(180/np.pi)
            lv = Perez_v3(zenDeg, azmt, radD, radI, jday, patchchoice)   # Relative luminance
            aniLum = np.zeros((rows, cols))
            for idx in range(0,145):
                aniLum = aniLum + diffsh[:,:,idx] * lv[0][idx][2]     # Total relative luminance from sky into each cell

            dRad = aniLum * radD   # Total diffuse radiation from sky into each cell
        else:
            dRad = radD * svfbuveg
            lv = 0
        #print(dRad)

        # Shadow  images
        if usevegdem == 1:
            vegsh, sh, vbshvegsh, wallsh, wallsun, wallshve, facesh, facesun = shadowingfunction_wallheight_23(dsm,vegdem, vegdem2, azmt, alt, res, amaxvalue, bush, walls, (dirwalls * np.pi / 180.)) 
            #shade = shadow_withveg_only(dsm,vegdsm,azmt,alt,res)            
            shadow = sh - (1 - vegsh) * (1 - psi)
             #shade = sh_tot
        #sh-(1-vegsh)*(1-psi)
        else:
            sh, wallsh, wallsun, facesh, facesun = shadowingfunction_wallheight_13(dsm, azmt, alt, res, walls, dirwalls * np.pi / 180.)
            shadow = sh
        print('on the fly shadow done')
        print(datetime.datetime.now())
        # # # Surface temperature parameterisation during daytime # # # #
        # new using max sun alt.instead of  dfm
        Tgamp = (TgK * altmax - Tstart) + Tstart
        Tgampwall = (TgK_wall * altmax - (Tstart_wall)) + (Tstart_wall)
        Tg = Tgamp * np.sin((((dectime - np.floor(dectime)) - SNUP / 24) / (TmaxLST / 24 - SNUP / 24)) * np.pi / 2) + Tstart # 2015 a, based on max sun altitude
        Tgwall = Tgampwall * np.sin((((dectime - np.floor(dectime)) - SNUP / 24) / (TmaxLST_wall / 24 - SNUP / 24)) * np.pi / 2) + (Tstart_wall) # 2015a, based on max sun altitude

        if Tgwall < 0:  # temporary for removing low Tg during morning 20130205
            # Tg = 0
            Tgwall = 0

        # New estimation of Tg reduction for non - clear situation based on Reindl et al.1990
        radI0, _ = diffusefraction(I0, alt, 1., Ta, RH)
        corr = 0.1473 * np.log(90 - (zen / np.pi * 180)) + 0.3454  # 20070329 correction of lat, Lindberg et al. 2008
        CI_Tg = (radI / radI0) + (1 - corr)
        if (CI_Tg > 1) or (CI_Tg == np.inf):
            CI_Tg = 1
        Tg = Tg * CI_Tg  # new estimation
        Tgwall = Tgwall * CI_Tg
        if landcover == 1:
            Tg[Tg < 0] = 0  # temporary for removing low Tg during morning 20130205

        # # # # Ground View Factors # # # #
        gvfLup, gvfalb, gvfalbnosh, gvfLupE, gvfalbE, gvfalbnoshE, gvfLupS, gvfalbS, gvfalbnoshS, gvfLupW, gvfalbW,\
        gvfalbnoshW, gvfLupN, gvfalbN, gvfalbnoshN, gvfSum, gvfNorm = gvf_2018a(wallsun, walls, buildings, res, shadow, first,
                second, dirwalls, Tg, Tgwall, Ta, emis_grid, ewall, alb_grid, SBC, albedo_b, rows, cols,
                                                                 Twater, lc_grid, landcover)
        print ('gvf done')
        print(datetime.datetime.now())
        # # # # Lup, daytime # # # #
        # Surface temperature wave delay - new as from 2014a
        Lup, timeaddnotused, Tgmap1 = TsWaveDelay_2015a(gvfLup, firstdaytime, timeadd, timestepdec, Tgmap1)
        LupE, timeaddnotused, Tgmap1E = TsWaveDelay_2015a(gvfLupE, firstdaytime, timeadd, timestepdec, Tgmap1E)
        LupS, timeaddnotused, Tgmap1S = TsWaveDelay_2015a(gvfLupS, firstdaytime, timeadd, timestepdec, Tgmap1S)
        LupW, timeaddnotused, Tgmap1W = TsWaveDelay_2015a(gvfLupW, firstdaytime, timeadd, timestepdec, Tgmap1W)
        LupN, timeaddnotused, Tgmap1N = TsWaveDelay_2015a(gvfLupN, firstdaytime, timeadd, timestepdec, Tgmap1N)
        print('l rad done')
        print(datetime.datetime.now())
        # # For Tg output in POIs
        TgTemp = Tg * shadow + Ta
        TgOut, timeadd, TgOut1 = TsWaveDelay_2015a(TgTemp, firstdaytime, timeadd, timestepdec, TgOut1) #timeadd only here v2021a

        # Building height angle from svf
        F_sh = cylindric_wedge(zen, svfalfa, rows, cols)  # Fraction shadow on building walls based on sun alt and svf
        F_sh[np.isnan(F_sh)] = 0.5

        # # # # # # # Calculation of shortwave daytime radiative fluxes # # # # # # #
        Kdown = radI * shadow * np.sin(alt * (np.pi / 180)) + dRad + albedo_b * (1 - svfbuveg) * (radG * (1 - F_sh) + radD * F_sh) # *sin(altitude(i) * (pi / 180))

        #Kdown = radI * shadow * np.sin(altitude * (np.pi / 180)) + radD * svfbuveg + albedo_b * (1 - svfbuveg) * \
                            #(radG * (1 - F_sh) + radD * F_sh) # *sin(altitude(i) * (pi / 180))

        Kup, KupE, KupS, KupW, KupN = Kup_veg_2015a(radI, radD, radG, alt, svfbuveg, albedo_b, F_sh, gvfalb,
                    gvfalbE, gvfalbS, gvfalbW, gvfalbN, gvfalbnosh, gvfalbnoshE, gvfalbnoshS, gvfalbnoshW, gvfalbnoshN)

        Keast, Ksouth, Kwest, Knorth, KsideI, KsideD = Kside_veg_v2019a(radI, radD, radG, shadow, svfS, svfW, svfN, svfE,
                    svfEveg, svfSveg, svfWveg, svfNveg, azmt, alt, psi, t, albedo_b, F_sh, KupE, KupS, KupW,
                    KupN, cyl, lv, ani, diffsh, rows, cols)
        print('k rad done')
        print(datetime.datetime.now())

        firstdaytime = 0

    else:  # # # # # # # NIGHTTIME # # # # # # # #
        CI=1 ##may change here

        Tgwall = 0
        # CI_Tg = -999  # F_sh = []

        # Nocturnal K fluxes set to 0
        Knight = np.zeros((rows, cols))
        Kdown = np.zeros((rows, cols))
        Kwest = np.zeros((rows, cols))
        Kup = np.zeros((rows, cols))
        Keast = np.zeros((rows, cols))
        Ksouth = np.zeros((rows, cols))
        Knorth = np.zeros((rows, cols))
        KsideI = np.zeros((rows, cols))
        KsideD = np.zeros((rows, cols))
        F_sh = np.zeros((rows, cols))
        Tg = np.zeros((rows, cols))
        shadow = np.zeros((rows, cols))

        # # # # Lup # # # #
        Lup = SBC * emis_grid * ((Knight + Ta + Tg + 273.15) ** 4)
        if landcover == 1:
            Lup[lc_grid == 3] = SBC * 0.98 * (Twater + 273.15) ** 4  # nocturnal Water temp

        LupE = Lup
        LupS = Lup
        LupW = Lup
        LupN = Lup

        # # For Tg output in POIs
        TgOut = Ta + Tg

        I0 = 0
        timeadd = 0
        firstdaytime = 1

    # # # # Ldown # # # #
    Ldown = (svf + svfveg - 1) * esky * SBC * ((Ta + 273.15) ** 4) + (2 - svfveg - svfaveg) * ewall * SBC *  ((Ta + 273.15) ** 4) + (svfaveg - svf) * ewall * SBC * ((Ta + 273.15 + Tgwall) ** 4) + (2 - svf - svfveg) * (1 - ewall) * esky * SBC * ((Ta + 273.15) ** 4)  # Jonsson et al.(2006)
    # Ldown = Ldown - 25 # Shown by Jonsson et al.(2006) and Duarte et al.(2006)

    #if CI < 0.95:  # non - clear conditions
    #    c = 1 - CI
    #    Ldown = Ldown * (1 - c) + c * ((svf + svfveg - 1) * SBC * ((Ta + 273.15) ** 4) + (2 - svfveg - svfaveg) *ewall * SBC * ((Ta + 273.15) ** 4) + (svfaveg - svf) * ewall * SBC * ((Ta + 273.15 + Tgwall) ** 4) +(2 - svf - svfveg) * (1 - ewall) * SBC * ((Ta + 273.15) ** 4))  # NOT REALLY TESTED!!! BUT MORE CORRECT?

    # # # # Lside # # # #
    Least, Lsouth, Lwest, Lnorth = Lside_veg_v2015a(svfS, svfW, svfN, svfE, svfEveg, svfSveg, svfWveg, svfNveg,
                    svfEaveg, svfSaveg, svfWaveg, svfNaveg, azmt, alt, Ta, Tgwall, SBC, ewall, Ldown,
                                                      esky, t, F_sh, CI, LupE, LupS, LupW, LupN)

    # # # # Calculation of radiant flux density and Tmrt # # # #
    if cyl == 1 and ani == 1:  # Human body considered as a cylinder with Perez et al. (1993)
        Sstr = absK * ((KsideI + KsideD) * Fcyl + (Kdown + Kup) * Fup + (Knorth + Keast + Ksouth + Kwest) * Fside) + absL * (Ldown * Fup + Lup * Fup + Lnorth * Fside + Least * Fside + Lsouth * Fside + Lwest * Fside)
    elif cyl == 1 and ani == 0: # Human body considered as a cylinder with isotropic all-sky diffuse
	    Sstr = absK * (KsideI * Fcyl + (Kdown + Kup) * Fup + (Knorth + Keast + Ksouth + Kwest) * Fside) + absL * (Ldown * Fup + Lup * Fup + Lnorth * Fside + Least * Fside + Lsouth * Fside + Lwest * Fside)
    else: # Human body considered as a standing cube
	    Sstr = absK * ((Kdown + Kup) * Fup + (Knorth + Keast + Ksouth + Kwest) * Fside) +absL * (Ldown * Fup + Lup * Fup + Lnorth * Fside + Least * Fside + Lsouth * Fside + Lwest * Fside)
    print('calculating mrt')
    print(datetime.datetime.now())

    Tmrt = np.sqrt(np.sqrt((Sstr / (absL * SBC)))) - 273.2
    
    #if write==True:
    # mrt = rio.open(out, 'w', driver = 'GTiff', height = a.shape[0], width = a.shape[1], count = 1,
	#					  crs= DSM.crs,transform=DSM.transform,dtype=dirwalls.dtype)
	#					  
	# vf.write(dirwalls,1)

	# vf.close
    print('mrt done')
    print(datetime.datetime.now())
    

    return {'Tmrt': Tmrt,'Kdown': Kdown, 'Kup':Kup,'Ldown': Ldown,'Lup': Lup,'shadow':shadow,
	'Keast':Keast, 'Ksouth':Ksouth, 'Kwest':Kwest, 'Knorth':Knorth, 
	'Least':Least,'Lsouth':Lsouth, 'Lwest':Lwest, 'Lnorth':Lnorth}
	
	#'shade':shade'Tg': Tg, ea, esky, I0, CI, shadow, firstdaytime, timestepdec, \
    #timeadd, Tgmap1, Tgmap1E, Tgmap1S, Tgmap1W, Tgmap1N, Keast, Ksouth, Kwest, Knorth, Least, \
   #Lsouth, Lwest, Lnorth, KsideI, TgOut1, TgOut, radI, radD}