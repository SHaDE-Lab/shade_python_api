##############################################################################################################################################
# surface temperature delay functions for Solweig                                                                                            #
# source :https://github.com/UMEP-dev/UMEP/tree/SuPy-QGIS3/SOLWEIG/SOLWEIGpython                                                             #  
# Goteborg Urban Climate Group                                                                                                               #
# Gothenburg University                                                                                                                      #
#                                                                                                                                            #
#                                                                                                                                            #
#                                                                                                                                            #
##############################################################################################################################################

import numpy as np


def TsWaveDelay_2015a(gvfLup, firstdaytime, timeadd, timestepdec):
    Tgmap1 = np.zeros_like(gvfLup)
    if firstdaytime == 1:  # "first in morning"
        Tgmap1 = gvfLup

    if timeadd >= (59 / 1440):  # more or equal to 59 min
        weight1 = np.exp(-33.27 * timeadd)  # surface temperature delay function - 1 step
        Tgmap1 = gvfLup * (1 - weight1) + Tgmap1 * weight1
        Lup = Tgmap1
        if timestepdec > (59 / 1440):
            timeadd = timestepdec
        else:
            timeadd = 0
    else:
        timeadd = timeadd + timestepdec
        weight1 = np.exp(-33.27 * timeadd)  # surface temperature delay function - 1 step
        Lup = (gvfLup * (1 - weight1) + Tgmap1 * weight1)

    return Lup
