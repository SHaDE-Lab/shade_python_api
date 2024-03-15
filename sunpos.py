############################################################################################################################################
# This scripts calculates sun position and time parameters for the mean radiant temperature and shade estimations using SOLWEIG            #
# maxalt and sun distance source :                                                                                                         #
# https://github.com/UMEP-dev/UMEP/blob/6fc4a1adaf5b66096b4d3dd558b867d44a04496c/Utilities/SEBESOLWEIGCommonFiles/sun_distance.py          #
# Goteborg Urban Climate Group                                                                                                             #
# Gothenburg University                                                                                                                    #
#                                                                                                                                          #
#                                                                                                                                          #
#                                                                                                                                          #
############################################################################################################################################


import numpy as np
import datetime
import pytz
import pysolar.solar as ps


# import pvlib as pv


def timestamp(year, month, day, hour, minute, timezone):
    # get appropiate timezone here : https://gist.github.com/heyalexej/8bf688fd67d7199be4a1682b3eec7568

    timezone = pytz.timezone(timezone)

    return datetime.datetime(year, month, day, hour, minute, tzinfo=timezone)


def azimuth(lat, lon, date):
    az = ps.get_azimuth(lat, lon, date)

    return round(az, 4)


def altitude(lat, lon, date):
    alt = ps.get_altitude(lat, lon, date)
    return round(alt, 4)


def zenith(lat, lon, date):
    alt = ps.get_altitude(lat, lon, date)
    zen = 90 - alt
    return round(zen, 4)


def sun_distance(jday):
    """
    #% Calculatesrelative earth sun distance
    #% with day of year as input.
    #% Partridge and Platt, 1975

    """
    b = 2. * np.pi * jday / 365.
    D = np.sqrt((1.00011 + np.dot(0.034221, np.cos(b)) + np.dot(0.001280, np.sin(b)) + np.dot(0.000719, np.cos(
        (2. * b))) + np.dot(0.000077, np.sin((2. * b)))))
    return D


# def  altitude(lat,lon,date):
#	alt = ps.get_altitude(lat,lon,date)
#	return round(alt,4)


def maxalt(lat, doy, year):
    # leap year check and declination cal

    if ((year % 400 == 0) or ((year % 4 == 0) and (year % 100 != 0))):
        dec = 23.45 * np.sin(((2 * np.pi) / 366) * (doy + 284))
    else:
        dec = 23.45 * np.sin(((2 * np.pi) / 365) * (doy + 284))

    # if the declination and the latitude are on the same side of the hemisphere
    if (lat > 0 and dec > 0) or (lat < 0 and dec < 0):
        if lat > dec:
            maxalt = 90 - (lat - dec)
        else:
            maxalt = 90 - (dec - lat)

            # if the latitude and the declination are no on the same side of the hemisphere
    else:
        maxalt = 90 - (lat + dec)

    # correct alt if greater than 90.
    if maxalt > 90:
        maxalt = 180 - maxalt

    return maxalt