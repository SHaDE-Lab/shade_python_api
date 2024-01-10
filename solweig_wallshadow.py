##############################################################################################################################################
# Wall and shadowing functions for Solweig                                                                                                     #
# source: https://github.com/UMEP-dev/UMEP/tree/SuPy-QGIS3/Utilities/SEBESOLWEIGCommonFiles                                                  #  
# source: https://github.com/UMEP-dev/UMEP/tree/SuPy-QGIS3/WallHeight                                                                        # 
# Goteborg Urban Climate Group                                                                                                               #
# Gothenburg University                                                                                                                      #
#                                                                                                                                            #
#                                                                                                                                            #
#                                                                                                                                            #
##############################################################################################################################################

from __future__ import division
import numpy as np
import rasterio as rio
import scipy.ndimage.interpolation as sc
import math



# import matplotlib.pylab as plt

def findwalls(dsm, walllimit,outfolder):
    
    # This function identifies walls based on a DSM and a wall-height limit
    # Walls are represented by outer pixels within building footprints
    #
    # Fredrik Lindberg, Goteborg Urban Climate Group
    # fredrikl@gvc.gu.se
    # 20150625
    
    DSM = rio.open(dsm) 
    
    a = DSM.read(1)
    
    root = dsm[dsm.find('box'):dsm.find('_dsm.tif')]
    
    out = outfolder + root+ '_walls' + '.tif'

    col = a.shape[0]
    row = a.shape[1]
    walls = np.zeros((col, row))
    domain = np.array([[0, 1, 0], [1, 0, 1], [0, 1, 0]])
    for i in np.arange(1, row-1):
        for j in np.arange(1, col-1):
            dom = a[j-1:j+2, i-1:i+2]
            # walls[j, i] = np.min(dom[np.where(domain == 1)])  # new 20171006
            walls[j, i] = np.max(dom[np.where(domain == 1)])  # new 20171006

    # walls = a-walls  # new 20171006
    walls = np.copy(walls - a)  # new 20171006
    walls[(walls < walllimit)] = 0

    walls[0:walls .shape[0], 0] = 0
    walls[0:walls .shape[0], walls .shape[1] - 1] = 0
    walls[0, 0:walls .shape[0]] = 0
    walls[walls .shape[0] - 1, 0:walls .shape[1]] = 0
    
    vf = rio.open(out, 'w', driver = 'GTiff', height = a.shape[0], width = a.shape[1], count = 1,crs= DSM.crs,transform=DSM.transform,dtype=walls.dtype)
    vf.write(walls,1)
    vf.close

    return walls


def filter1Goodwin_as_aspect_v3(res,dsm,outfolder):
    """
    tThis function applies the filter processing presented in Goodwin et al (2010) but instead for removing
    linear fetures it calculates wall aspect based on a wall pixels grid, a dsm (a) and a res factor
    Fredrik Lindberg, 2012-02-14
    fredrikl@gvc.gu.se
    Translated: 2015-09-15
    :param walls:
    :param scale:
    :param a:
    :return: dirwalls
    """
    DSM = rio.open(dsm) 
    
    a = DSM.read(1)
    
    root = dsm[dsm.find('box'):dsm.find('_dsm.tif')]
	
	
    
    out = outfolder + root+ '_dirwalls' + '.tif'

    row = a.shape[0]
    col = a.shape[1]
    walls=findwalls(dsm, 3,outfolder)
    filtersize = np.floor((res + 0.0000000001) * 9)
    if filtersize <= 2:
        filtersize = 3
    else:
        if filtersize != 9:
            if filtersize % 2 == 0:
                filtersize = filtersize + 1

    filthalveceil = int(np.ceil(filtersize / 2.))
    filthalvefloor = int(np.floor(filtersize / 2.))

    filtmatrix = np.zeros((int(filtersize), int(filtersize)))
    buildfilt = np.zeros((int(filtersize), int(filtersize)))

    filtmatrix[:, filthalveceil - 1] = 1
    n = filtmatrix.shape[0] - 1
    buildfilt[filthalveceil - 1, 0:filthalvefloor] = 1
    buildfilt[filthalveceil - 1, filthalveceil: int(filtersize)] = 2

    y = np.zeros((row, col))  # final direction
    z = np.zeros((row, col))  # temporary direction
    x = np.zeros((row, col))  # building side
    walls[walls > 0] = 1

    for h in range(0, 180):  # =0:1:180 #%increased resolution to 1 deg 20140911
        filtmatrix1temp = sc.rotate(filtmatrix, h, order=1, reshape=False, mode='nearest')  # bilinear
        filtmatrix1 = np.round(filtmatrix1temp)
        # filtmatrix1temp = sc.imrotate(filtmatrix, h, 'bilinear')
        # filtmatrix1 = np.round(filtmatrix1temp / 255.)
        # filtmatrixbuildtemp = sc.imrotate(buildfilt, h, 'nearest')
        filtmatrixbuildtemp = sc.rotate(buildfilt, h, order=0, reshape=False, mode='nearest')  # Nearest neighbor
        # filtmatrixbuild = np.round(filtmatrixbuildtemp / 127.)
        filtmatrixbuild = np.round(filtmatrixbuildtemp)
        index = 270 - h
        if h == 150:
            filtmatrixbuild[:, n] = 0
        if h == 30:
            filtmatrixbuild[:, n] = 0
        if index == 225:
            # n = filtmatrix.shape[0] - 1  # length(filtmatrix);
            filtmatrix1[0, 0] = 1
            filtmatrix1[n, n] = 1
        if index == 135:
            # n = filtmatrix.shape[0] - 1  # length(filtmatrix);
            filtmatrix1[0, n] = 1
            filtmatrix1[n, 0] = 1

        for i in range(int(filthalveceil) - 1, row - int(filthalveceil) - 1):  # i=filthalveceil:sizey-filthalveceil
            for j in range(int(filthalveceil) - 1,
                           col - int(filthalveceil) - 1):  # (j=filthalveceil:sizex-filthalveceil
                if walls[i, j] == 1:
                    wallscut = walls[i - filthalvefloor:i + filthalvefloor + 1,
                               j - filthalvefloor:j + filthalvefloor + 1] * filtmatrix1
                    dsmcut = a[i - filthalvefloor:i + filthalvefloor + 1, j - filthalvefloor:j + filthalvefloor + 1]
                    if z[i, j] < wallscut.sum():  # sum(sum(wallscut))
                        z[i, j] = wallscut.sum()  # sum(sum(wallscut));
                        if np.sum(dsmcut[filtmatrixbuild == 1]) > np.sum(dsmcut[filtmatrixbuild == 2]):
                            x[i, j] = 1
                        else:
                            x[i, j] = 2

                        y[i, j] = index

    y[(x == 1)] = y[(x == 1)] - 180
    y[(y < 0)] = y[(y < 0)] + 360

    grad, asp = get_ders(a, res)

    y = y + ((walls == 1) * 1) * ((y == 0) * 1) * (asp / (math.pi / 180.))

    dirwalls = y
    
    vf = rio.open(out, 'w', driver = 'GTiff', height = a.shape[0], width = a.shape[1], count = 1,crs= DSM.crs,transform=DSM.transform,dtype=dirwalls.dtype)

    vf.write(dirwalls,1)

    vf.close

    return dirwalls


def cart2pol(x, y, units='deg'):
    radius = np.sqrt(x**2 + y**2)
    theta = np.arctan2(y, x)
    if units in ['deg', 'degs']:
        theta = theta * 180 / np.pi
    return theta, radius


def get_ders(dsm, res):
    # dem,_,_=read_dem_grid(dem_file)
    dx = 1/res
    # dx=0.5
    fy, fx = np.gradient(dsm, dx, dx)
    asp, grad = cart2pol(fy, fx, 'rad')
    grad = np.arctan(grad)
    asp = asp * -1
    asp = asp + (asp < 0) * (np.pi * 2)
    return grad, asp


def shadowingfunction_wallheight_23(a, vegdem, vegdem2, azimuth, altitude, scale, amaxvalue, bush, walls, aspect):
    """
    This function calculates shadows on a DSM and shadow height on building
    walls including both buildings and vegetion units.
    New functionallity to deal with pergolas, August 2021
    
    INPUTS:
    a = DSM
    vegdem = Vegetation canopy DSM (magl)
    vegdem2 = Trunkzone DSM (magl)
    azimuth and altitude = sun position
    scale= scale of DSM (1 meter pixels=1, 2 meter pixels=0.5)
    walls= pixel row 'outside' buildings. will be calculated if empty
    aspect=normal aspect of walls
    
    OUTPUT:
    sh=ground and roof shadow
    wallsh = height of wall that is in shadow
    wallsun = hieght of wall that is in sun
    original Matlab code:
    Fredrik Lindberg 2013-08-14
    fredrikl@gvc.gu.se
    :param a:
    :param vegdem:
    :param vegdem2:
    :param azimuth:
    :param altitude:
    :param scale:
    :param amaxvalue:
    :param bush:
    :param walls:
    :param aspect:
    :return:
    """

    # conversion
    degrees = np.pi/180.
    azimuth *= degrees
    altitude *= degrees
    
    # measure the size of the image
    sizex = np.shape(a)[0]
    sizey = np.shape(a)[1]
    
    # initialise parameters
    dx = 0
    dy = 0
    dz = 0
    temp = np.zeros((sizex, sizey))
    tempvegdem = np.zeros((sizex, sizey))
    tempvegdem2 = np.zeros((sizex, sizey))
    templastfabovea = np.zeros((sizex, sizey))
    templastgabovea = np.zeros((sizex, sizey))
    bushplant = bush > 1
    sh = np.zeros((sizex, sizey)) #shadows from buildings
    vbshvegsh = np.copy(sh) #vegetation blocking buildings
    vegsh = np.add(np.zeros((sizex, sizey)), bushplant, dtype=float) #vegetation shadow
    f = np.copy(a)
    shvoveg = np.copy(vegdem) # for vegetation shadowvolume
    # g = np.copy(sh)
    wallbol = (walls > 0).astype(float)

    # other loop parameters
    pibyfour = np.pi/4
    threetimespibyfour = 3*pibyfour
    fivetimespibyfour = 5*pibyfour
    seventimespibyfour = 7*pibyfour
    sinazimuth = np.sin(azimuth)
    cosazimuth = np.cos(azimuth)
    tanazimuth = np.tan(azimuth)
    signsinazimuth = np.sign(sinazimuth)
    signcosazimuth = np.sign(cosazimuth)
    dssin = np.abs(1/sinazimuth)
    dscos = np.abs(1/cosazimuth)
    tanaltitudebyscale = np.tan(altitude)/scale

    index = 0

    # new case with pergola (thin vertical layer of vegetation), August 2021
    dzprev = 0

    # main loop
    while (amaxvalue >= dz) and (np.abs(dx) < sizex) and (np.abs(dy) < sizey):
        if ((pibyfour <= azimuth) and (azimuth < threetimespibyfour)) or ((fivetimespibyfour <= azimuth) and (azimuth < seventimespibyfour)):
            dy = signsinazimuth * index
            dx = -1 * signcosazimuth * np.abs(np.round(index / tanazimuth))
            ds = dssin
        else:
            dy = signsinazimuth * np.abs(np.round(index * tanazimuth))
            dx = -1 * signcosazimuth * index
            ds = dscos

        # note: dx and dy represent absolute values while ds is an incremental value
        dz = (ds * index) * tanaltitudebyscale
        tempvegdem[0:sizex, 0:sizey] = 0
        tempvegdem2[0:sizex, 0:sizey] = 0
        temp[0:sizex, 0:sizey] = 0
        templastfabovea[0:sizex, 0:sizey] = 0.
        templastgabovea[0:sizex, 0:sizey] = 0.
        absdx = np.abs(dx)
        absdy = np.abs(dy)
        xc1 = int((dx+absdx)/2)
        xc2 = int(sizex+(dx-absdx)/2)
        yc1 = int((dy+absdy)/2)
        yc2 = int(sizey+(dy-absdy)/2)
        xp1 = -int((dx-absdx)/2)
        xp2 = int(sizex-(dx+absdx)/2)
        yp1 = -int((dy-absdy)/2)
        yp2 = int(sizey-(dy+absdy)/2)

        tempvegdem[xp1:xp2, yp1:yp2] = vegdem[xc1:xc2, yc1:yc2] - dz
        tempvegdem2[xp1:xp2, yp1:yp2] = vegdem2[xc1:xc2, yc1:yc2] - dz
        temp[xp1:xp2, yp1:yp2] = a[xc1:xc2, yc1:yc2]-dz

        f = np.fmax(f, temp) #Moving building shadow
        shvoveg = np.fmax(shvoveg, tempvegdem) # moving vegetation shadow volume
        sh[f > a] = 1
        sh[f <= a] = 0   
        fabovea = (tempvegdem > a).astype(int)   #vegdem above DEM
        gabovea = (tempvegdem2 > a).astype(int)   #vegdem2 above DEM
        
        #new pergola condition
        templastfabovea[xp1:xp2, yp1:yp2] = vegdem[xc1:xc2, yc1:yc2]-dzprev
        templastgabovea[xp1:xp2, yp1:yp2] = vegdem2[xc1:xc2, yc1:yc2]-dzprev
        lastfabovea = templastfabovea > a
        lastgabovea = templastgabovea > a
        dzprev = dz
        vegsh2 = np.add(np.add(np.add(fabovea, gabovea, dtype=float),lastfabovea, dtype=float),lastgabovea, dtype=float)
        vegsh2[vegsh2 == 4] = 0.
        # vegsh2[vegsh2 == 1] = 0. # This one is the ultimate question...
        vegsh2[vegsh2 > 0] = 1.

        # vegsh2 = fabovea - gabovea #old without pergolas
        # vegsh = np.max([vegsh, vegsh2], axis=0) #old without pergolas

        vegsh = np.fmax(vegsh, vegsh2)
        vegsh[vegsh*sh > 0] = 0    
        vbshvegsh = np.copy(vegsh) + vbshvegsh # removing shadows 'behind' buildings

        # # vegsh at high sun altitudes # Not needed when pergolas are included
        # if index == 0:
        #     firstvegdem = np.copy(tempvegdem) - np.copy(temp)
        #     firstvegdem[firstvegdem <= 0] = 1000
        #     vegsh[firstvegdem < dz] = 1
        #     vegsh *= (vegdem2 > a)
        #     vbshvegsh = np.zeros((sizex, sizey))

        # # Bush shadow on bush plant # Not needed when pergolas are included
        # if np.max(bush) > 0 and np.max(fabovea*bush) > 0:
        #     tempbush = np.zeros((sizex, sizey))
        #     tempbush[int(xp1):int(xp2), int(yp1):int(yp2)] = bush[int(xc1):int(xc2), int(yc1):int(yc2)] - dz
        #     g = np.max([g, tempbush], axis=0)
        #     g = bushplant * g
    
        index += 1

    # Removing walls in shadow due to selfshadowing
    azilow = azimuth - np.pi/2
    azihigh = azimuth + np.pi/2
    if azilow >= 0 and azihigh < 2*np.pi:    # 90 to 270  (SHADOW)
        facesh = np.logical_or(aspect < azilow, aspect >= azihigh).astype(float) - wallbol + 1    # TODO check
    elif azilow < 0 and azihigh <= 2*np.pi:    # 0 to 90
        azilow = azilow + 2*np.pi
        facesh = np.logical_or(aspect > azilow, aspect <= azihigh) * -1 + 1    # (SHADOW)
    elif azilow > 0 and azihigh >= 2*np.pi:    # 270 to 360
        azihigh -= 2 * np.pi
        facesh = np.logical_or(aspect > azilow, aspect <= azihigh)*-1 + 1    # (SHADOW)

    sh = 1-sh
    vbshvegsh[vbshvegsh > 0] = 1
    vbshvegsh = vbshvegsh-vegsh
    
    # if np.max(bush) > 0: # Not needed when pergolas are included
    #     g = g-bush
    #     g[g > 0] = 1
    #     g[g < 0] = 0
    #     vegsh = vegsh-bushplant+g
    #     vegsh[vegsh < 0] = 0

    vegsh[vegsh > 0] = 1
    shvoveg = (shvoveg-a) * vegsh    #Vegetation shadow volume
    vegsh = 1-vegsh
    vbshvegsh = 1-vbshvegsh
    
    # wall shadows
    shvo = f - a   # building shadow volume
    facesun = np.logical_and(facesh + (walls > 0).astype(float) == 1, walls > 0).astype(float)
    wallsun = np.copy(walls-shvo)
    wallsun[wallsun < 0] = 0
    wallsun[facesh == 1] = 0    # Removing walls in "self"-shadow
    wallsh = np.copy(walls-wallsun)

    wallshve = shvoveg * wallbol
    wallshve = wallshve - wallsh
    wallshve[wallshve < 0] = 0
    id = np.where(wallshve > walls)
    wallshve[id] = walls[id]
    wallsun = wallsun-wallshve    # problem with wallshve only
    id = np.where(wallsun < 0)
    wallshve[id] = 0
    wallsun[id] = 0
    
    return vegsh, sh, vbshvegsh, wallsh, wallsun, wallshve, facesh, facesun
	












