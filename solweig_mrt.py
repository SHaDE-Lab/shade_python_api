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

def Solweig_2021a_calc(dsm, vegdsm, dem, res, trans, svf, svfN, svfW, svfE, svfS, svfveg, svfNveg, svfEveg, svfSveg,
                       svfWveg, svfaveg, svfEaveg, svfSaveg, svfWaveg, svfNaveg, walls, dirwalls, location, tzone, year,
                       month, day, doy, hour, minu, Ws, Ta, RH, radG, Twater, ani, cyl, usevegdem, onlyglobal, elvis,
                       landcover, lc_grid):
    from sunpos import altitude
    from sunpos import azimuth
    from sunpos import zenith
    from sunpos import maxalt
    from sunpos import timestamp
    from clearness_index import clearnessindex
    from diffusefraction import diffusefraction
    from solweig_cylindric_wedge_svfalfa import cylindric_wedge
    from solweig_daylen import daylen
    from solweig_gvf_sunonsurface import gvf_2018a
    from solweig_k_rads import Kup_veg_2015a, Kside_veg_v2019a
    from solweig_l_rads import Lside_veg_v2015a
    from solweig_tswavedelay import TsWaveDelay_2015a
    from solweig_wallshadow import shadowingfunction_wallheight_23
    import numpy as np
    import datetime
    import gc

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

    time = timestamp(year, month, day, hour, minu, tzone)

    alt = altitude(lat, lon, time)  # 90. - sun['zenith']
    azmt = azimuth(lat, lon, time)  # sun['azimuth']
    zen = zenith(lat, lon, time) * (np.pi / 180.)

    albedo_b = 0.20
    albedo_g = 0.15
    ewall = 0.90
    eground = 0.95
    absK = 0.70
    absL = 0.95

    P = -999.0

    # timeconstants
    timestepdec = 0
    timeadd = 0.
    firstdaytime = 1.

    # constants for a standing person
    Fside = 0.22  # change
    Fup = 0.06
    Fcyl = 0.28

    # building footprints from DSM and DEM

    buildings = dsm - dem
    buildings[buildings < 2.] = 1.
    buildings[buildings >= 2.] = 0.

    ###########
    rows = dsm.shape[0]
    cols = dsm.shape[1]

    trunkratio = 25 / 100.
    psi = trans / 100.0

    vegdsm2 = vegdsm * trunkratio
    vegmax = vegdsm.max()
    amaxvalue = dsm.max() - dsm.min()
    amaxvalue = np.maximum(amaxvalue, vegmax)

    # Elevation vegdsms if buildingDSM includes ground heights
    vegdem = vegdsm + dsm
    vegdem[vegdem == dsm] = 0
    vegdem2 = vegdsm2 + dsm
    vegdem2[vegdem2 == dsm] = 0

    # Bush separation
    bush = np.logical_not((vegdem2 * vegdem)) * vegdem

    svfbuveg = svf - (1. - svfveg) * (1. - psi)

    # %Initialization of maps
    print('initializing maps', flush=True)
    Knight = np.zeros((rows, cols))

    TgK = Knight + 0.37
    Tstart = Knight - 3.41
    alb_grid = Knight + albedo_g
    emis_grid = Knight + eground
    TgK_wall = 0.37
    Tstart_wall = -3.41
    TmaxLST = 15.
    TmaxLST_wall = 15.

    # %Parameterisarion for Lup
    height = 1.1

    diffsh = None  # condiering isotropic sky

    # %Radiative surface influence, Rule of thumb by Schmid et al. (1990).
    first = np.round(height)
    if first == 0.:
        first = 1.
    second = np.round((height * 20.))

    #############wall calculations################
    # walls = findwalls (dsm,3.0)

    # dirwalls = filter1Goodwin_as_aspect_v3(walls,res,dsm)

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

    if alt > 0:  # # # # # # DAYTIME # # # # # #
        print('daytime', flush=True)
        # Clearness Index on Earth's surface after Crawford and Dunchon (1999) with a correction
        #  factor for low sun elevations after Lindberg et al.(2008)
        I0, CI, Kt, I0et, CIuncorr = clearnessindex(zen, doy, Ta, RH, radG, location, P)
        if (CI > 1) or (CI == np.inf):
            CI = 1

        # Estimation of radD and radI if not measured after Reindl et al.(1990)
        if onlyglobal == 1:
            I0, CI, Kt, I0et, CIuncorr = clearnessindex(zen, doy, Ta, RH, radG, location, P)
            if (CI > 1) or (CI == np.inf):
                CI = 1

        radI, radD = diffusefraction(radG, alt, Kt, Ta, RH)
        print('rad diffuse done', flush=True)
        print(radI, radD)
        print(datetime.datetime.now())

        dectime = doy + hour / 24 + minu / (60 * 24.)
        altmax = maxalt(lat, doy, year)  # 80.15 # to be calculated
        print('rad dectime done', flush=True)
        print(datetime.datetime.now())

        # Diffuse Radiation
        # Anisotropic Diffuse Radiation after Perez et al. 1993
        dRad = radD * svfbuveg
        lv = 0
        print('on the fly shadow starting', flush=True)
        # Shadow  images
        vegsh, sh, wallsun = shadowingfunction_wallheight_23(dsm, vegdem, vegdem2, azmt, alt, res, amaxvalue, bush,
                                                             walls, (dirwalls * np.pi / 180.))
        del dsm, vegdem, vegdem2, bush
        shadow = sh - (1 - vegsh) * (1 - psi)

        print('on the fly shadow done', flush=True)
        print(datetime.datetime.now())

        # # # Surface temperature parameterisation during daytime # # # #
        # new using max sun alt.instead of
        Tgamp = (TgK * altmax - Tstart) + Tstart
        del TgK
        Tgampwall = (TgK_wall * altmax - Tstart_wall) + Tstart_wall
        del TgK_wall
        Tg = Tgamp * np.sin((((dectime - np.floor(dectime)) - SNUP / 24) / (
                TmaxLST / 24 - SNUP / 24)) * np.pi / 2) + Tstart  # 2015 a, based on max sun altitude
        del Tgamp, Tstart
        Tgwall = Tgampwall * np.sin((((dectime - np.floor(dectime)) - SNUP / 24) / (
                TmaxLST_wall / 24 - SNUP / 24)) * np.pi / 2) + Tstart_wall  # 2015a, based on max sun altitude
        del Tgampwall, Tstart_wall, TmaxLST_wall, SNUP
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
        del radI0
        gc.collect()
        # # # # Ground View Factors # # # #
        print('gvf starting', flush=True)
        gvfLup, gvfalb, gvfalbnosh, gvfLupE, gvfalbE, gvfalbnoshE, gvfLupS, gvfalbS, gvfalbnoshS, gvfLupW, gvfalbW, \
            gvfalbnoshW, gvfLupN, gvfalbN, gvfalbnoshN = gvf_2018a(wallsun, walls, buildings, res,
                                                                                    shadow, first,
                                                                                    second, dirwalls, Tg, Tgwall, Ta,
                                                                                    emis_grid, ewall, alb_grid, SBC,
                                                                                    albedo_b, rows, cols,
                                                                                    Twater, lc_grid, landcover)
        del emis_grid, alb_grid
        print('gvf done', flush=True)
        print(datetime.datetime.now())
        # # # # Lup, daytime # # # #
        # Surface temperature wave delay - new as from 2014a
        print('l rad starting', flush=True)
        Lup = TsWaveDelay_2015a(gvfLup, firstdaytime, timeadd, timestepdec)
        del gvfLup
        LupE = TsWaveDelay_2015a(gvfLupE, firstdaytime, timeadd, timestepdec)
        del gvfLupE
        LupS = TsWaveDelay_2015a(gvfLupS, firstdaytime, timeadd, timestepdec)
        del gvfLupS
        LupW = TsWaveDelay_2015a(gvfLupW, firstdaytime, timeadd, timestepdec)
        del gvfLupW
        LupN = TsWaveDelay_2015a(gvfLupN, firstdaytime, timeadd, timestepdec)
        del gvfLupN
        print('l rad done', flush=True)
        print(datetime.datetime.now())
        # # For Tg output in POIs
        # TgOut, timeadd, TgOut1 = TsWaveDelay_2015a(TgTemp, firstdaytime, timeadd, timestepdec, TgOut1) #timeadd only here v2021a

        # Building height angle from svf
        F_sh = cylindric_wedge(zen, svfalfa, rows, cols)  # Fraction shadow on building walls based on sun alt and svf
        del svfalfa, zen
        F_sh[np.isnan(F_sh)] = 0.5

        print('k rad starting', flush=True)
        # # # # # # # Calculation of shortwave daytime radiative fluxes # # # # # # #
        Kdown = radI * shadow * np.sin(alt * (np.pi / 180)) + dRad + albedo_b * (1 - svfbuveg) * (
                radG * (1 - F_sh) + radD * F_sh)  # *sin(altitude(i) * (pi / 180))
        del dRad
        Kup, KupE, KupS, KupW, KupN = Kup_veg_2015a(radI, radD, radG, alt, svfbuveg, albedo_b, F_sh, gvfalb,
                                                    gvfalbE, gvfalbS, gvfalbW, gvfalbN, gvfalbnosh, gvfalbnoshE,
                                                    gvfalbnoshS, gvfalbnoshW, gvfalbnoshN)
        del svfbuveg, gvfalbE, gvfalbS, gvfalbW, gvfalbN, gvfalbnoshE, gvfalbnoshS, gvfalbnoshW, gvfalbnoshN
        Keast, Ksouth, Kwest, Knorth, KsideI, KsideD = Kside_veg_v2019a(radI, radD, radG, shadow, svfS, svfW, svfN,
                                                                        svfE,
                                                                        svfEveg, svfSveg, svfWveg, svfNveg, azmt, alt,
                                                                        psi, t, albedo_b, F_sh, KupE, KupS, KupW,
                                                                        KupN, cyl, lv, ani, diffsh, rows, cols)
        del radD, KupE, KupS, KupW, KupN,shadow, albedo_b, diffsh, rows, cols, lv
        print('k rad done', flush=True)
        print(datetime.datetime.now())

    else:  # # # # # # # NIGHTTIME # # # # # # # #
        print('nighttime', flush=True)

        CI = 1  # may change here

        Tgwall = 0
        # CI_Tg = -999  # F_sh = []

        # Nocturnal K fluxes set to 0
        # omitted to save memory

        # # # # Lup # # # #
        Lup = SBC * emis_grid * ((Ta + 273.15) ** 4)
        del emis_grid
        if landcover == 1:
            Lup[lc_grid == 3] = SBC * 0.98 * (Twater + 273.15) ** 4  # nocturnal Water temp
        del lc_grid

    # # # # Ldown # # # #
    # Ldown = (svf + svfveg - 1) * esky * SBC * ((Ta + 273.15) ** 4) + (2 - svfveg - svfaveg) * ewall * SBC *  ((Ta + 273.15) ** 4) + (svfaveg - svf) * ewall * SBC * ((Ta + 273.15 + Tgwall) ** 4) + (2 - svf - svfveg) * (1 - ewall) * esky * SBC * ((Ta + 273.15) ** 4)  # Jonsson et al.(2006)
    # Ldown = Ldown - 25 # Shown by Jonsson et al.(2006) and Duarte et al.(2006)
    # In-place operations and avoiding unnecessary intermediates
    Ldown = svf + svfveg
    Ldown *= esky * SBC * ((Ta + 273.15) ** 4)

    Ldown += (2 - svfveg - svfaveg) * ewall * SBC * ((Ta + 273.15) ** 4)
    Ldown += (svfaveg - svf) * ewall * SBC * ((Ta + 273.15 + Tgwall) ** 4)
    Ldown += (2 - svf - svfveg) * (1 - ewall) * esky * SBC * ((Ta + 273.15) ** 4)

    # # # # Lside # # # #
    if alt > 0:
        Least, Lsouth, Lwest, Lnorth = Lside_veg_v2015a(svfS, svfW, svfN, svfE, svfEveg, svfSveg, svfWveg, svfNveg,
                                                        svfEaveg, svfSaveg, svfWaveg, svfNaveg, azmt, alt, Ta, Tgwall,
                                                        SBC,
                                                        ewall, Ldown,
                                                        esky, t, F_sh, CI, LupE, LupS, LupW, LupN)
    else:
        Least, Lsouth, Lwest, Lnorth = Lside_veg_v2015a(svfS, svfW, svfN, svfE, svfEveg, svfSveg, svfWveg, svfNveg,
                                                        svfEaveg, svfSaveg, svfWaveg, svfNaveg, azmt, alt, Ta, Tgwall,
                                                        SBC,
                                                        ewall, Ldown,
                                                        esky, t, Lup, CI, Lup, Lup, Lup, Lup)

    # # # # Calculation of radiant flux density and Tmrt # # # #
    # if cyl == 1 and ani == 1:  # Human body considered as a cylinder with Perez et al. (1993)
    #     Sstr = absK * ((KsideI + KsideD) * Fcyl + (Kdown + Kup) * Fup + (Knorth + Keast + Ksouth + Kwest) * Fside) + absL * (Ldown * Fup + Lup * Fup + Lnorth * Fside + Least * Fside + Lsouth * Fside + Lwest * Fside)
    # elif cyl == 1 and ani == 0: # Human body considered as a cylinder with isotropic all-sky diffuse
    #     Sstr = absK * (KsideI * Fcyl + (Kdown + Kup) * Fup + (Knorth + Keast + Ksouth + Kwest) * Fside) + absL * (Ldown * Fup + Lup * Fup + Lnorth * Fside + Least * Fside + Lsouth * Fside + Lwest * Fside)
    # else: # Human body considered as a standing cube
    #     Sstr = absK * ((Kdown + Kup) * Fup + (Knorth + Keast + Ksouth + Kwest) * Fside) +absL * (Ldown * Fup + Lup * Fup + Lnorth * Fside + Least * Fside + Lsouth * Fside + Lwest * Fside)
    print('calculating mrt', flush=True)

    print(datetime.datetime.now())
    if alt > 0:
        Sstr = absK * (KsideI * Fcyl + (Kdown + Kup) * Fup + (Knorth + Keast + Ksouth + Kwest) * Fside) + absL * (
                Ldown * Fup + Lup * Fup + Lnorth * Fside + Least * Fside + Lsouth * Fside + Lwest * Fside)
    else:  # night time all the k is 0 so we can skip it
        Sstr = absK * (absL * (
                Ldown * Fup + Lup * Fup + Lnorth * Fside + Least * Fside + Lsouth * Fside + Lwest * Fside))
    Tmrt = np.sqrt(np.sqrt((Sstr / (absL * SBC)))) - 273.2


    print('mrt done', flush=True)
    print(datetime.datetime.now())

    return Tmrt
