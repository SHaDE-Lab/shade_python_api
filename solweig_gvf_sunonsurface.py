##############################################################################################################################################
# Sunonsurface and ground view factor functions to estimate upwelling longwave rad and  for Solweig                                          #
# source :https://github.com/UMEP-dev/UMEP/tree/SuPy-QGIS3/SOLWEIG/SOLWEIGpython                                                             #  
# Goteborg Urban Climate Group                                                                                                               #
# Gothenburg University                                                                                                                      #
#                                                                                                                                            #
#                                                                                                                                            #
#                                                                                                                                            #
##############################################################################################################################################


import numpy as np


def sunonsurface_2018a(azimuthA, res, buildings, shadow, sunwall, first, second, aspect, walls, Tg, Tgwall, Ta,
                       emis_grid, ewall, alb_grid, SBC, albedo_b, Twater, lc_grid, landcover):
    # This version of sunonsurfce goes with SOLWEIG 2015a. It also simulates
    # Lup and albedo based on landcover information and shadow patterns.
    # Fredrik Lindberg, fredrikl@gvc.gu.se

    sizex = np.shape(walls)[0]
    sizey = np.shape(walls)[1]

    # sizex=size(buildings,1);sizey=size(buildings,2);
    wallbol = (walls > 0) * 1
    sunwall[sunwall > 0] = 1  # test 20160910

    # conversion into radians
    azimuth = azimuthA * (np.pi / 180)

    # loop parameters
    index = 0
    Lup = SBC * emis_grid * (Tg * shadow + Ta + 273.15) ** 4 - SBC * emis_grid * (Ta + 273.15) ** 4  # +Ta
    if landcover == 1:
        Tg[lc_grid == 3] = Twater - Ta  # Setting water temperature

    Lwall = SBC * ewall * (Tgwall + Ta + 273.15) ** 4 - SBC * ewall * (Ta + 273.15) ** 4  # +Ta
    albshadow = alb_grid * shadow

    weightsumsh = np.zeros((sizex, sizey))
    weightsumwall = np.zeros((sizex, sizey))
    first = np.round(first * res)
    if first < 1:
        first = 1
    second = np.round(second * res)
    weightsumLupsh = np.zeros((sizex, sizey))
    weightsumalbsh = np.zeros((sizex, sizey))
    weightsumalbnosh = np.zeros((sizex, sizey))
    weightsumalbwallnosh = np.zeros((sizex, sizey))

    # other loop parameters
    pibyfour = np.pi / 4
    threetimespibyfour = 3 * pibyfour
    fivetimespibyfour = 5 * pibyfour
    seventimespibyfour = 7 * pibyfour
    sinazimuth = np.sin(azimuth)
    cosazimuth = np.cos(azimuth)
    tanazimuth = np.tan(azimuth)
    signsinazimuth = np.sign(sinazimuth)
    signcosazimuth = np.sign(cosazimuth)

    # The Shadow casting algorithm
    for n in np.arange(0, second):
        if (pibyfour <= azimuth and azimuth < threetimespibyfour) or (
                fivetimespibyfour <= azimuth and azimuth < seventimespibyfour):
            dy = signsinazimuth * index
            dx = -1 * signcosazimuth * np.abs(np.round(index / tanazimuth))
        else:
            dy = signsinazimuth * abs(round(index * tanazimuth))
            dx = -1 * signcosazimuth * index

        absdx = np.abs(dx)
        absdy = np.abs(dy)

        xc1 = ((dx + absdx) / 2)
        xc2 = (sizex + (dx - absdx) / 2)
        yc1 = ((dy + absdy) / 2)
        yc2 = (sizey + (dy - absdy) / 2)

        xp1 = -((dx - absdx) / 2)
        xp2 = (sizex - (dx + absdx) / 2)
        yp1 = -((dy - absdy) / 2)
        yp2 = (sizey - (dy + absdy) / 2)

        tempbu = np.zeros_like(buildings)
        tempbu[int(xp1):int(xp2), int(yp1):int(yp2)] = buildings[int(xc1):int(xc2),
                                                       int(yc1):int(yc2)]  # moving building
        tempsh = np.zeros_like(shadow)
        tempsh[int(xp1):int(xp2), int(yp1):int(yp2)] = shadow[int(xc1):int(xc2), int(yc1):int(yc2)]  # moving shadow

        tempLupsh = np.zeros_like(Lup)
        tempLupsh[int(xp1):int(xp2), int(yp1):int(yp2)] = Lup[int(xc1):int(xc2), int(yc1):int(yc2)]  # moving Lup/shadow

        tempalbsh = np.zeros_like(albshadow)
        tempalbsh[int(xp1):int(xp2), int(yp1):int(yp2)] = albshadow[int(xc1):int(xc2),
                                                          int(yc1):int(yc2)]  # moving Albedo/shadow
        tempalbnosh = np.zeros_like(alb_grid)
        tempalbnosh[int(xp1):int(xp2), int(yp1):int(yp2)] = alb_grid[int(xc1):int(xc2), int(yc1):int(yc2)]  # moving Albedo

        tempbu = np.min([buildings, tempbu], axis=0)  # utsmetning av buildings

        weightsumsh += tempsh * tempbu

        weightsumLupsh += tempLupsh * tempbu

        weightsumalbsh += tempalbsh * tempbu

        weightsumalbnosh += tempalbnosh * tempbu

        tempwallsun = np.zeros_like(sunwall)
        tempwallsun[int(xp1):int(xp2), int(yp1):int(yp2)] = sunwall[int(xc1):int(xc2),
                                                            int(yc1):int(yc2)]  # moving buildingwall insun image
        tempb = tempwallsun * tempbu
        tempbwall = tempbu * -1 + 1
        tempbub = (tempb > 0) * 1
        tempbubwall = (tempbwall > 0) * 1
        weightsumwall += tempbub
        weightsumalbwallnosh += tempbubwall * albedo_b

        ind = 1
        if (n + 1) <= first:
            weightsumwall_first = weightsumwall / ind
            weightsumsh_first = weightsumsh / ind
            wallsuninfluence_first = weightsumwall_first > 0
            weightsumLwall_first = weightsumwall / ind * Lwall
            weightsumLupsh_first = weightsumLupsh / ind

            weightsumalbwall_first = weightsumwall / ind  * albedo_b
            weightsumalbsh_first = weightsumalbsh / ind
            weightsumalbwallnosh_first = weightsumalbwallnosh / ind  # *albedo_b
            weightsumalbnosh_first = weightsumalbnosh / ind
            wallinfluence_first = weightsumalbwallnosh_first > 0
            #         gvf1=(weightsumwall+weightsumsh)/first;
            #         gvf1(gvf1>1)=1;
            ind += 1
        index += 1

    wallsuninfluence_second = weightsumwall > 0
    wallinfluence_second = weightsumalbwallnosh > 0
    # gvf2(gvf2>1)=1;

    # Removing walls in shadow due to selfshadowing
    azilow = azimuth - np.pi / 2
    azihigh = azimuth + np.pi / 2
    if azilow >= 0 and azihigh < 2 * np.pi:  # 90 to 270  (SHADOW)
        facesh = (np.logical_or(aspect < azilow, aspect >= azihigh).astype(float) - wallbol + 1)
    elif azilow < 0 and azihigh <= 2 * np.pi:  # 0 to 90
        azilow = azilow + 2 * np.pi
        facesh = np.logical_or(aspect > azilow, aspect <= azihigh) * -1 + 1  # (SHADOW)    # check for the -1
    elif azilow > 0 and azihigh >= 2 * np.pi:  # 270 to 360
        azihigh = azihigh - 2 * np.pi
        facesh = np.logical_or(aspect > azilow, aspect <= azihigh) * -1 + 1  # (SHADOW)

    # removing walls in self shadowing
    keep = (weightsumwall == second) - facesh
    keep[keep == -1] = 0
    weightsumwall[keep == 1] = 0

    # gvf from shadow and Lup
    gvfLup1 = ((weightsumLwall_first + weightsumLupsh_first) / (first + 1)) * wallsuninfluence_first + \
              (weightsumLupsh_first) / (first) * (wallsuninfluence_first * -1 + 1)
    gvfLup2 = ((weightsumwall * Lwall + weightsumLupsh) / (second + 1)) * wallsuninfluence_second + \
              (weightsumLupsh) / (second) * (wallsuninfluence_second * -1 + 1)

    # gvf from shadow and albedo
    gvfalb1 = ((weightsumalbwall_first + weightsumalbsh_first) / (first + 1)) * wallsuninfluence_first + \
              (weightsumalbsh_first) / (first) * (wallsuninfluence_first * -1 + 1)
    gvfalb2 = ((weightsumwall * albedo_b + weightsumalbsh) / (second + 1)) * wallsuninfluence_second + \
              (weightsumalbsh) / (second) * (wallsuninfluence_second * -1 + 1)

    # gvf from albedo only
    gvfalbnosh1 = ((weightsumalbwallnosh_first + weightsumalbnosh_first) / (first + 1)) * wallinfluence_first + \
                  (weightsumalbnosh_first) / (first) * (wallinfluence_first * -1 + 1)  #
    gvfalbnosh2 = ((weightsumalbwallnosh + weightsumalbnosh) / (second)) * wallinfluence_second + \
                  (weightsumalbnosh) / (second) * (wallinfluence_second * -1 + 1)

    # Weighting
    gvfLup = (gvfLup1 * 0.5 + gvfLup2 * 0.4) / 0.9
    gvfLup += ((SBC * emis_grid * (Tg * shadow + Ta + 273.15) ** 4) - SBC * emis_grid * (Ta + 273.15) ** 4) * (
            buildings * -1 + 1)
    gvfalb = (gvfalb1 * 0.5 + gvfalb2 * 0.4) / 0.9
    gvfalb += alb_grid * (buildings * -1 + 1) * shadow
    gvfalbnosh = (gvfalbnosh1 * 0.5 + gvfalbnosh2 * 0.4) / 0.9
    gvfalbnosh *= buildings + alb_grid * (buildings * -1 + 1)

    return gvfLup, gvfalb, gvfalbnosh


def gvf_2018a(wallsun, walls, buildings, res, shadow, first, second, dirwalls, Tg, Tgwall, Ta, emis_grid, ewall,
              alb_grid, SBC, albedo_b, rows, cols, Twater, lc_grid, landcover):
    azimuthA = np.arange(5, 359, 20)  # Search directions for Ground View Factors (GVF)

    #### Ground View Factors ####
    gvfLupE = 0
    gvfLupS = 0
    gvfLupW = 0
    gvfLupN = 0
    gvfalbE = 0
    gvfalbS = 0
    gvfalbW = 0
    gvfalbN = 0
    gvfalbnoshE = 0
    gvfalbnoshS = 0
    gvfalbnoshW = 0
    gvfalbnoshN = 0

    #  sunwall=wallinsun_2015a(buildings,azimuth(i),shadow,psi(i),dirwalls,walls);
    sunwall = (wallsun / walls * buildings) == 1  # new as from 2015a
    # 5 25 45 65 85 105 125 145 165 185 205 225 245 265 285 305 325 345
    for j in np.arange(0, azimuthA.__len__()):
        gvfLupi, gvfalbi, gvfalbnoshi = sunonsurface_2018a(azimuthA[j], res, buildings, shadow, sunwall,
                                                                    first,
                                                                    second, dirwalls * np.pi / 180, walls, Tg, Tgwall,
                                                                    Ta,
                                                                    emis_grid, ewall, alb_grid, SBC, albedo_b, Twater,
                                                                    lc_grid, landcover)
        if (azimuthA[j] >= 0) and (azimuthA[j] < 180):
            gvfLupE += gvfLupi
            gvfalbE += gvfalbi
            gvfalbnoshE += gvfalbnoshi

        if (azimuthA[j] >= 90) and (azimuthA[j] < 270):
            gvfLupS += gvfLupi
            gvfalbS += gvfalbi
            gvfalbnoshS += gvfalbnoshi

        if (azimuthA[j] >= 180) and (azimuthA[j] < 360):
            gvfLupW += gvfLupi
            gvfalbW += gvfalbi
            gvfalbnoshW += gvfalbnoshi

        if (azimuthA[j] >= 270) or (azimuthA[j] < 90):
            gvfLupN += gvfLupi
            gvfalbN += gvfalbi
            gvfalbnoshN += gvfalbnoshi

    bias = SBC * emis_grid * (Ta + 273.15) ** 4


    gvfLup = gvfLupN + gvfLupS + gvfLupE + gvfLupW
    gvfLup /= azimuthA.__len__()
    gvfLup += bias

    gvfalb = gvfalbN + gvfalbS + gvfalbE + gvfalbW
    gvfalb /= azimuthA.__len__()

    gvfalbnosh = gvfalbnoshN + gvfalbnoshS + gvfalbnoshE + gvfalbnoshW
    gvfalbnosh /= azimuthA.__len__()

    gvfLupE /= azimuthA.__len__() / 2
    gvfLupE += bias
    gvfLupS /= azimuthA.__len__() / 2
    gvfLupS += bias
    gvfLupW /= azimuthA.__len__() / 2
    gvfLupW += bias
    gvfLupN /= azimuthA.__len__() / 2
    gvfLupN += bias

    gvfalbE /= (azimuthA.__len__() / 2)
    gvfalbS /= (azimuthA.__len__() / 2)
    gvfalbW /= (azimuthA.__len__() / 2)
    gvfalbN /= (azimuthA.__len__() / 2)

    gvfalbnoshE /= (azimuthA.__len__() / 2)
    gvfalbnoshS /= (azimuthA.__len__() / 2)
    gvfalbnoshW /= (azimuthA.__len__() / 2)
    gvfalbnoshN /= (azimuthA.__len__() / 2)

    return gvfLup, gvfalb, gvfalbnosh, gvfLupE, gvfalbE, gvfalbnoshE, gvfLupS, gvfalbS, gvfalbnoshS, gvfLupW, gvfalbW, gvfalbnoshW, gvfLupN, gvfalbN, gvfalbnoshN
