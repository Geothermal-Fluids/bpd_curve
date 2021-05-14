# -*- coding: utf-8 -*-
"""
Created on Thu May 13 10:27:29 2021

@author: t.hoerbrand
"""


import numpy as np
import pandas as pd
import iapws
import pint as pint
from plotnine import *
from collections import namedtuple
from warnings import warn

# setup pint
u = pint.UnitRegistry()


# ----------------------------------------------------  
# functions
# ----------------------------------------------------

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

def hydrostatic_bpdc(depth, p0=None, method="iapws", units='common'):
    '''
    Function to calculate the hydrostatic boiling point depth curve, similar to:
    Haas Jr., J.L., 1971. The effect of salinity on the maximum thermal gradient 
    of a hydrothermal system at hydrostatic pressure. Econ. Geol. 66, 940–946.

    Parameters
    ----------
    depth : array-like
        An array of depths at which to evaluate the function. 
        Depth intervals should not be bigger than 1 m!
    p : float, optional
        surface pressure of well. The default is 1.01325 bar.
    method : string
        Method to calculate brine density. The default is "iapws".
    units : string or tuple
        String units for depth, pressure, temperature, and density
        respectively. The following are allowed:
        - 'common' to use metres, bar, deg Celsius, and kg/m**3.
        - 'SI' to use metres, Pascal, Kelvin, and kg/m**3.
        - 'Imperial' to use feet, psia, deg Fahrenheit, and g/cm**3.
        - tuple like ('m', 'atm', 'degC', 'g/cm**3') for custom units.

    Returns
    -------
    namedtuple
    '''
    u = UnitRegistry()

    g = 9.81 * u.m / u.s**2

    # Assign units to everything.
    if units == 'SI':
        u_d, u_p, u_t, u_rho = u.m, u.pa, u.K, u.kg/u.m**3
    if units == 'common':
        u_d, u_p, u_t, u_rho = u.m, u.bar, u('degC'), u.kg/u.m**3
    elif units == 'Imperial':
        u_d, u_p, u_t, u_rho = u.ft, u.psi, u('degF'), u.g/u.cm**3
    else:
        u_d, u_p, u_t, u_rho = list(map(u, units))
    
    # Override units with pint Quantity, if possible.
    # And compute the diff array.
    if isinstance(depth, pint.Quantity):
        u_d = depth.units
    else:
        depth = np.asanyarray(depth) * u_d
    depth_diff = np.diff(depth)
    if max(depth_diff.magnitude) > 1:
        warn("Step size of depth interval > 1. This leads to an error propagation in calculated temperatures. Consider reducing the interval size.")
                
    # Override units with pint Quantity, if possible.
    # And deal with not getting a p0 pressure.
    if isinstance(p0, pint.Quantity):
        u_p = p0.units
    elif p0 is None:
        p0 = 101325 * u_p
    else:
        p0 = p0 * u_p  
    p0 = p0.to('MPa')  # Use MPa for calculations (required by IAPWS).

    # Compute using IAPWS option.
    if method == "iapws":
        pressure = np.atleast_1d(p0)
        tsat = iapws.iapws97._TSat_P(p0.m) * u.K
        rho = iapws.iapws97.IAPWS97_PT(p0.m, tsat.m).rho * u.kg / u.m**3
        density = np.atleast_1d(rho)
        
        for i in depth_diff:
            # Calculate new pressure for this step.
            new_p = pressure[-1] + rho * g * i
            pressure = np.append(pressure, new_p)  # Has units.
            
            # Calculate new temperature for this step.
            t = iapws.iapws97._TSat_P(new_p.m) * u.K
            tsat = np.append(tsat, t)
            
            # Calculate new density for next step.
            rho = iapws.iapws97.IAPWS97_PT(new_p.m, t.m).rho * u.kg / u.m**3
            density = np.append(density, rho)
        
        # Finalize units.
        pressure = pressure.to(u_p)
        tsat = tsat.to(u_t)
        density = density.to(u_rho)
        
    # Return a namedtuple to the user to retain units.
    BPD = namedtuple('BPD', ['depth', 'pressure', 'tsat', 'rho'])
    return BPD(depth, pressure, tsat, density)


# ----------------------------------------------------  
# data and bpd from ICWallis
# ---------------------------------------------------- 


# load real data
df = pd.read_csv(r'C:/Users/t.hoerbrand/Documents/github/bpd_curve/data/Data-Temp-Heating37days.csv')
df['pressure_bara'] = df.pres_barg - 1.01325

# calculate bpd curve for pure water (iapws)
df['pressure_mpa'] = df.pressure_bara * 0.1  # convert pressure to MPa for ipaws
pressure = df['pressure_mpa'].tolist()
tsat = []
for p in pressure:
    saturation_temp = iapws.iapws97._TSat_P(p) - 273.15  # calculate saturation temp in Kelvin & convert to degC
    tsat.append(saturation_temp)
df['tsat_degC'] = tsat


    
# ----------------------------------------------------  
# calculate hydrostatic pressure
# ---------------------------------------------------- 
    

hyd = hydrostatic_bpdc(depth=df['depth_m'], 
                       p0=df["pressure_bara"].min(),
                       method='iapws',
                       units=['m', 'bar', 'degC', 'kg/m**3'],
                      )

# retain units and convert to DataFrame
hyd_df = pd.DataFrame(hyd).T  # I don't really know why we have to transpose here.
hyd_df.columns = columns=hyd._fields

# remove units and convert to DataFrame
hyd_df = pd.DataFrame(hyd._asdict())


# ----------------------------------------------------  
# compare results
# ---------------------------------------------------- 
gd = ggplot(df, aes(y = "depth_m")) +\
    geom_line(aes(x="tsat_degC", colour = "'well pressure + iapws'")) + \
    geom_line(aes(x = "temp_degC", colour = "'well data only'")) + \
    geom_line(aes(x = "tsat", y = "depth", colour = "'hydrostatic + iapws'"), data = hyd_df) + \
    scale_y_reverse() + \
    labs(x = "T [°C]", y = "depth [m]", colour = "legend")
    
gp = ggplot(df, aes(y = "pressure_bara")) +\
    geom_line(aes(x="tsat_degC", colour = "'well pressure + iapws'")) + \
    geom_line(aes(x = "temp_degC", colour = "'well data only'")) + \
    geom_line(aes(x = "tsat", y = "pressure", colour = "'hydrostatic + iapws'"), data = hyd_df) + \
    scale_y_reverse() + \
    labs(x = "T [°C]", y = "pressure [bara]", colour = "legend")
    
ggplot(df, aes(y = "depth_m")) + \
    geom_line(aes(x="pressure_bara", colour = "'well pressure + iapws'")) + \
    geom_line(aes(x = "pressure_bara", colour = "'well data only'")) + \
    geom_line(aes(x = "pressure", y = "depth", colour = "'hydrostatic + iapws'"), data = hyd_df) + \
    scale_y_reverse() + \
    labs(x = "T [°C]", y = "depth [m]", colour = "legend")

        
ggsave(plot = gd, filename = "C:/Users/t.hoerbrand/Documents/temp/gd.png", width = 6, height = 9)
ggsave(plot = gp, filename = "C:/Users/t.hoerbrand/Documents/temp/gp.png", width = 6, height = 9)



        