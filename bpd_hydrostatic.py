# -*- coding: utf-8 -*-
"""
Created on Thu May 13 10:27:29 2021

@author: t.hoerbrand
"""


import numpy as np
import pandas as pd
import iapws
from plotnine import *



# ----------------------------------------------------  
# functions
# ----------------------------------------------------

def hydrostatic_bpdc(dmin, dmax, steps, p_bar = 1.01325, method = "iapws"):
    '''
    Function to calculate the hydrostatic boiling point depth curve, similar to:
    Haas Jr., J.L., 1971. The effect of salinity on the maximum thermal gradient 
    of a hydrothermal system at hydrostatic pressure. Econ. Geol. 66, 940–946.

    Parameters
    ----------
    dmin : float
        mimimum depth for bpd [m]
    dmax : float
        maximum depth for bpd [m]
    steps : integer
        number of steps for integral
    p_bar : float, optional
        surface pressure of well [bara]. The default is 1.01325.
    method : string
        Method to calculate brine density. The default is "iapws".

    Returns
    -------
    df : pd.DataFrame
        

    '''
    g = 9.81
    depth = np.linspace(dmin, dmax, steps)
    depth_diff = np.diff(depth)

    if method == "iapws":
        pressure = np.array([p_bar*0.1]) # conversion to MPa for iaps
        tsat = iapws.iapws97._TSat_P(pressure[0]) #conversion to K for iapws  
        rho = iapws.iapws97.IAPWS97_PT(p_bar, tsat).rho # calculate density of first step
        density = np.array([rho])
        
        for i in depth_diff:
            # calculate new pressure for this step
            dp = rho*g*i*1e-5*0.1 # pressure difference in MPa
            p = pressure[-1] + dp
            # calculate new temperature for this step
            saturation_temp = iapws.iapws97._TSat_P(p) 
            # calculate new density for next step
            rho = iapws.iapws97.IAPWS97_PT(p, saturation_temp).rho
            # save pressure, temperature and density
            pressure = np.append(pressure, p)
            tsat = np.append(tsat, saturation_temp)
            density = np.append(density, rho)
        
        pressure = pressure*10
        tsat = tsat-273.15
        
        df = pd.DataFrame(list(zip(depth, pressure, tsat, density)), columns = ["depth_m", "p_bar", "t_c", "rho"])
        
    return df



# ----------------------------------------------------  
# data and bpd from ICWallis
# ---------------------------------------------------- 


# load real data
heating_37days = pd.read_csv(r'C:/Users/t.hoerbrand/Documents/github/bpd_curve/data/Data-Temp-Heating37days.csv')
heating_37days['pressure_bara'] = heating_37days.pres_barg - 1.01325

# calculate bpd curve for pure water (iapws)
heating_37days['pressure_mpa'] = heating_37days.pressure_bara * 0.1  # convert pressure to MPa for ipaws
pressure = heating_37days['pressure_mpa'].tolist()
tsat = []
for p in pressure:
    saturation_temp = iapws.iapws97._TSat_P(p) - 273.15  # calculate saturation temp in Kelvin & convert to degC
    tsat.append(saturation_temp)
heating_37days['tsat_degC'] = tsat


    
# ----------------------------------------------------  
# calculate hydrostatic pressure
# ---------------------------------------------------- 
    

hyd = hydrostatic_bpdc(min(heating_37days["depth_m"]), 
                           max(heating_37days["depth_m"]), 
                           1000, 
                           p_bar = min(heating_37days["pressure_bara"]))


# ----------------------------------------------------  
# compare results
# ---------------------------------------------------- 
ggplot(heating_37days, aes(y = "depth_m")) +\
    geom_line(aes(x="tsat_degC", colour = "'well pressure + iapws'")) + \
    geom_line(aes(x = "temp_degC", colour = "'well data only'")) + \
    geom_line(aes(x = "t_c", y = "depth_m", colour = "'hydrostatic + iapws'"), data = hyd) + \
    scale_y_reverse() + \
    labs(x = "T [°C]", y = "depth [m]", colour = "legend")
    
ggplot(heating_37days, aes(y = "pressure_bara")) +\
    geom_line(aes(x="tsat_degC", colour = "'well pressure + iapws'")) + \
    geom_line(aes(x = "temp_degC", colour = "'well data only'")) + \
    geom_line(aes(x = "t_c", y = "p_bar", colour = "'hydrostatic + iapws'"), data = hyd) + \
    scale_y_reverse() + \
    labs(x = "T [°C]", y = "depth [m]", colour = "legend")
    
ggplot(heating_37days, aes(y = "depth_m")) + \
    geom_line(aes(x="pressure_bara", colour = "'well pressure + iapws'")) + \
    geom_line(aes(x = "pressure_bara", colour = "'well data only'")) + \
    geom_line(aes(x = "p_bar", y = "depth_m", colour = "'hydrostatic + iapws'"), data = hyd) + \
    scale_y_reverse() + \
    labs(x = "T [°C]", y = "depth [m]", colour = "legend")

        
        
        