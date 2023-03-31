##############################################################################################################################################
# Fucntions to estimate all longwave radiations for Solweig                                                                                 #
# source :https://github.com/UMEP-dev/UMEP/tree/SuPy-QGIS3/SOLWEIG/SOLWEIGpython                                                             #  
# Goteborg Urban Climate Group                                                                                                               #
# Gothenburg University                                                                                                                      #
#                                                                                                                                            #
#                                                                                                                                            #
#                                                                                                                                            #
##############################################################################################################################################


import numpy as np

def Lvikt_veg(svf,svfveg,svfaveg,vikttot):

    # Least
    viktonlywall=(vikttot-(63.227*svf**6-161.51*svf**5+156.91*svf**4-70.424*svf**3+16.773*svf**2-0.4863*svf))/vikttot

    viktaveg=(vikttot-(63.227*svfaveg**6-161.51*svfaveg**5+156.91*svfaveg**4-70.424*svfaveg**3+16.773*svfaveg**2-0.4863*svfaveg))/vikttot

    viktwall=viktonlywall-viktaveg

    svfvegbu=(svfveg+svf-1)  # Vegetation plus buildings
    viktsky=(63.227*svfvegbu**6-161.51*svfvegbu**5+156.91*svfvegbu**4-70.424*svfvegbu**3+16.773*svfvegbu**2-0.4863*svfvegbu)/vikttot
    viktrefl=(vikttot-(63.227*svfvegbu**6-161.51*svfvegbu**5+156.91*svfvegbu**4-70.424*svfvegbu**3+16.773*svfvegbu**2-0.4863*svfvegbu))/vikttot
    viktveg=(vikttot-(63.227*svfvegbu**6-161.51*svfvegbu**5+156.91*svfvegbu**4-70.424*svfvegbu**3+16.773*svfvegbu**2-0.4863*svfvegbu))/vikttot
    viktveg=viktveg-viktwall

    return viktveg,viktwall,viktsky,viktrefl



def Lside_veg_v2015a(svfS,svfW,svfN,svfE,svfEveg,svfSveg,svfWveg,svfNveg,svfEaveg,svfSaveg,svfWaveg,svfNaveg,azimuth,altitude,Ta,Tw,SBC,ewall,Ldown,esky,t,F_sh,CI,LupE,LupS,LupW,LupN):

    # This m-file is the current one that estimates L from the four cardinal points 20100414
    
    #Building height angle from svf
    svfalfaE=np.arcsin(np.exp((np.log(1-svfE))/2))
    svfalfaS=np.arcsin(np.exp((np.log(1-svfS))/2))
    svfalfaW=np.arcsin(np.exp((np.log(1-svfW))/2))
    svfalfaN=np.arcsin(np.exp((np.log(1-svfN))/2))
    
    vikttot=4.4897
    aziW=azimuth+t
    aziN=azimuth-90+t
    aziE=azimuth-180+t
    aziS=azimuth-270+t
    
    F_sh = 2*F_sh-1  #(cylindric_wedge scaled 0-1)
    
    c=1-CI
    Lsky_allsky = esky*SBC*((Ta+273.15)**4)*(1-c)+c*SBC*((Ta+273.15)**4)
    
    ## Least
    [viktveg, viktwall, viktsky, viktrefl] = Lvikt_veg(svfE, svfEveg, svfEaveg, vikttot)
    
    if altitude > 0:  # daytime
        alfaB=np.arctan(svfalfaE)
        betaB=np.arctan(np.tan((svfalfaE)*F_sh))
        betasun=((alfaB-betaB)/2)+betaB
        # betasun = np.arctan(0.5*np.tan(svfalfaE)*(1+F_sh)) #TODO This should be considered in future versions
        if (azimuth > (180-t))  and  (azimuth <= (360-t)):
            Lwallsun=SBC*ewall*((Ta+273.15+Tw*np.sin(aziE*(np.pi/180)))**4)*\
                viktwall*(1-F_sh)*np.cos(betasun)*0.5
            Lwallsh=SBC*ewall*((Ta+273.15)**4)*viktwall*F_sh*0.5
        else:
            Lwallsun=0
            Lwallsh=SBC*ewall*((Ta+273.15)**4)*viktwall*0.5
    else: #nighttime
        Lwallsun=0
        Lwallsh=SBC*ewall*((Ta+273.15)**4)*viktwall*0.5
    
    Lsky=((svfE+svfEveg-1)*Lsky_allsky)*viktsky*0.5
    Lveg=SBC*ewall*((Ta+273.15)**4)*viktveg*0.5
    Lground=LupE*0.5
    Lrefl=(Ldown+LupE)*(viktrefl)*(1-ewall)*0.5
    Least=Lsky+Lwallsun+Lwallsh+Lveg+Lground+Lrefl
    
    # clear alfaB betaB betasun Lsky Lwallsh Lwallsun Lveg Lground Lrefl viktveg viktwall viktsky
    
    ## Lsouth
    [viktveg,viktwall,viktsky,viktrefl]=Lvikt_veg(svfS,svfSveg,svfSaveg,vikttot)
    
    if altitude>0: # daytime
        alfaB=np.arctan(svfalfaS)
        betaB=np.arctan(np.tan((svfalfaS)*F_sh))
        betasun=((alfaB-betaB)/2)+betaB
        # betasun = np.arctan(0.5*np.tan(svfalfaS)*(1+F_sh))
        if (azimuth <= (90-t))  or  (azimuth > (270-t)):
            Lwallsun=SBC*ewall*((Ta+273.15+Tw*np.sin(aziS*(np.pi/180)))**4)*\
                viktwall*(1-F_sh)*np.cos(betasun)*0.5
            Lwallsh=SBC*ewall*((Ta+273.15)**4)*viktwall*F_sh*0.5
        else:
            Lwallsun=0
            Lwallsh=SBC*ewall*((Ta+273.15)**4)*viktwall*0.5
    else: #nighttime
        Lwallsun=0
        Lwallsh=SBC*ewall*((Ta+273.15)**4)*viktwall*0.5

    Lsky=((svfS+svfSveg-1)*Lsky_allsky)*viktsky*0.5
    Lveg=SBC*ewall*((Ta+273.15)**4)*viktveg*0.5
    Lground=LupS*0.5
    Lrefl=(Ldown+LupS)*(viktrefl)*(1-ewall)*0.5
    Lsouth=Lsky+Lwallsun+Lwallsh+Lveg+Lground+Lrefl
    
    # clear alfaB betaB betasun Lsky Lwallsh Lwallsun Lveg Lground Lrefl viktveg viktwall viktsky
    
    ## Lwest
    [viktveg,viktwall,viktsky,viktrefl]=Lvikt_veg(svfW,svfWveg,svfWaveg,vikttot)
    
    if altitude>0: # daytime
        alfaB=np.arctan(svfalfaW)
        betaB=np.arctan(np.tan((svfalfaW)*F_sh))
        betasun=((alfaB-betaB)/2)+betaB
        # betasun = np.arctan(0.5*np.tan(svfalfaW)*(1+F_sh))
        if (azimuth > (360-t))  or  (azimuth <= (180-t)):
            Lwallsun=SBC*ewall*((Ta+273.15+Tw*np.sin(aziW*(np.pi/180)))**4)*\
                viktwall*(1-F_sh)*np.cos(betasun)*0.5
            Lwallsh=SBC*ewall*((Ta+273.15)**4)*viktwall*F_sh*0.5
        else:
            Lwallsun=0
            Lwallsh=SBC*ewall*((Ta+273.15)**4)*viktwall*0.5
    else: #nighttime
        Lwallsun=0
        Lwallsh=SBC*ewall*((Ta+273.15)**4)*viktwall*0.5

    Lsky=((svfW+svfWveg-1)*Lsky_allsky)*viktsky*0.5
    Lveg=SBC*ewall*((Ta+273.15)**4)*viktveg*0.5
    Lground=LupW*0.5
    Lrefl=(Ldown+LupW)*(viktrefl)*(1-ewall)*0.5
    Lwest=Lsky+Lwallsun+Lwallsh+Lveg+Lground+Lrefl
    
    # clear alfaB betaB betasun Lsky Lwallsh Lwallsun Lveg Lground Lrefl viktveg viktwall viktsky
    
    ## Lnorth
    [viktveg,viktwall,viktsky,viktrefl]=Lvikt_veg(svfN,svfNveg,svfNaveg,vikttot)
    
    if altitude>0: # daytime
        alfaB=np.arctan(svfalfaN)
        betaB=np.arctan(np.tan((svfalfaN)*F_sh))
        betasun=((alfaB-betaB)/2)+betaB
        # betasun = np.arctan(0.5*np.tan(svfalfaN)*(1+F_sh))
        if (azimuth > (90-t))  and  (azimuth <= (270-t)):
            Lwallsun=SBC*ewall*((Ta+273.15+Tw*np.sin(aziN*(np.pi/180)))**4)*\
                viktwall*(1-F_sh)*np.cos(betasun)*0.5
            Lwallsh=SBC*ewall*((Ta+273.15)**4)*viktwall*F_sh*0.5
        else:
            Lwallsun=0
            Lwallsh=SBC*ewall*((Ta+273.15)**4)*viktwall*0.5
    else: #nighttime
        Lwallsun=0
        Lwallsh=SBC*ewall*((Ta+273.15)**4)*viktwall*0.5

    Lsky=((svfN+svfNveg-1)*Lsky_allsky)*viktsky*0.5
    Lveg=SBC*ewall*((Ta+273.15)**4)*viktveg*0.5
    Lground=LupN*0.5
    Lrefl=(Ldown+LupN)*(viktrefl)*(1-ewall)*0.5
    Lnorth=Lsky+Lwallsun+Lwallsh+Lveg+Lground+Lrefl

    # clear alfaB betaB betasun Lsky Lwallsh Lwallsun Lveg Lground Lrefl viktveg viktwall viktsky
    
    return Least,Lsouth,Lwest,Lnorth