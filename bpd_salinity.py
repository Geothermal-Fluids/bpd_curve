# -*- coding: utf-8 -*-
"""
Created on Thu May 13 10:27:29 2021

@author: t.hoerbrand
"""


import numpy as np
import pandas as pd
import iapws
import pint as pint
from pint import UnitRegistry
from plotnine import *
import phreeqpy.iphreeqc.phreeqc_com as phreeqc_mod # install .com module for 64 bit from IPHREEQC
from collections import namedtuple
from warnings import warn

# setup pint
u = pint.UnitRegistry()


# ----------------------------------------------------  
# functions
# ----------------------------------------------------

def make_chemical_model(chem_ions, chem_gas_mol, t, p, chem_units = "mg/L", pH = 7):
    '''  
    Function to create a PHREEQC model for a brine to calculate the boiling point and density.

    Parameters
    ----------
    chem_ions_mol : pd.DataFrame
        Contains the ionic composition of the fluid with the chemical units defined in chem_units.
    chem_gas_mol : pd.DataFrame
        Contains the gas amount of the individual gases in mol/kgw.
    t : float
        Temperature [°C]
    p : float
        Pressure [bara]
    chem_units : string, optional
        Units in which the chemical composition is supplied. 
        Options are "mg/L" (equal to ppm), "mol/L" (molarity) and "mol/kgw" (molality). 
        The default is "mg/L" which is equal to ppm.
    pH : float, optional
        DESCRIPTION. The default is 7.
        
    Returns
    -------
    chemical model string
    
    '''
    p = p/1.01325 # convert to atm for PHREEQC
    
    chemical_model = """
        TITLE 
        SOLUTION 1
        """ + \
        "\n units " + chem_units + \
        "\n temp " + str(t) + \
        "\n pH " + str(pH) + \
        "\n pressure " + str(p) + \
        "\n Na " + str(chem_ions["Na"][0]) + \
        "\n K " + str(chem_ions["Na"][0]) + \
        "\n Ca " + str(chem_ions["Na"][0]) + \
        "\n S(6) " + str(chem_ions["Na"][0]) + \
        "\n HCO3 " + str(chem_ions["Na"][0]) + \
        "\n Cl " + str(chem_ions["Cl"][0]) + \
        "\n EQUILIBRIUM_PHASES"   + \
        "\n Mtg(g) 10 " + str(chem_gas_mol["CH4"][0]) + \
        "\n Ntg(g) 10 " + str(chem_gas_mol["N2"][0]) + \
        "\n Hdg(g) 10 " + str(chem_gas_mol["H2"][0]) + \
        "\n CO2(g) 10 " + str(chem_gas_mol["CO2"][0]) + """
        SELECTED_OUTPUT
        -reset false
        USER_PUNCH
        -headings CO2 N2 CH4 H2O H2 rho t
        10 PUNCH 10^SI("CO2(g)") * 1.01325 /PR_PHI("CO2(g)")
        20 PUNCH 10^SI("Ntg(g)") * 1.01325 /PR_PHI("Ntg(g)")
        30 PUNCH 10^SI("Mtg(g)") * 1.01325 /PR_PHI("Mtg(g)")
        40 PUNCH 10^SI("H2O(g)") * 1.01325 /PR_PHI("H2O(g)")
        50 PUNCH 10^SI("Hdg(g)") * 1.01325 /PR_PHI("Hdg(g)")
        60 PUNCH RHO*1000
        70 PUNCH TC
        END
        """
        
    return chemical_model
      
  

def run_chemcial_model(chemical_model, database_path):  
    """
    

    Parameters
    ----------
    chemical_model : string
        PHREEQC input string, generated with "make_chemical_model"
    database_path : string
        Path to thermodynamic database used with phreeqc, e.g. r"C:/phreeqc/database/phreeqc.dat"

    Returns
    -------
    output : tuple
        results of PHREEQC simulation (selected output)
        Contains header in first row [0], intialization in second row [1] and results in third row [2]

    """
    
    # Initialize IPhreeqc
    phreeqc = phreeqc_mod.IPhreeqc()
    
    # Load thermodynamic database provided by phreeqc
    phreeqc.load_database(database_path)
    
    # Run chemical model
    phreeqc.run_string(chemical_model)  

    # Get selected output    
    #components = phreeqc.get_component_list()
    output = phreeqc.get_selected_output_array()
    
    return output


def get_chemical_properties(output):
    """
    Format simulation results from "run_chemical_model".

    Parameters
    ----------
    output : tuple
        Results of function "run_chemical_model".

    Returns
    -------
    partial_pressure : float
        Sum of partial pressures [atm].
    rho_phreeqc : float
        Denisty in [kg/m**3].
    t_phreeqc : float
        Temperature in [°C]

    """
    
    # Indices / output format provided with the USER_PUNCH block in function "make_chemical_model"
    partial_pressure = sum(output[2][0:5]) # why the hell is it not 0:4 here? Whatever it works...
    partial_pressure = partial_pressure
    rho_phreeqc = output[2][5]
    t_phreeqc = output[2][6]
    
    return partial_pressure, rho_phreeqc, t_phreeqc
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


# All units "pinted" except for Temperature. 
def bpdc_salinity(depth, p0, chem_ions, chem_gas_mol, database_path, tolerance = 0.01, chem_units = "mol/kgw", units='common'):
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
    p0 = p0.to('MPa')  

    # Compute using PHREEQC option.
    pressure = np.atleast_1d(p0)
    # select starting temperature for iterate_chemical_models (reduces number of iterations)
    t_model = 100
    #print(chem_ions, chem_gas_mol, t_model, p0.to('atm').m, database_path, tolerance, chem_units)

    tsat = iapws.iapws97._TSat_P(p0.m) - 273.15
    
    
    # calculate density wwith chemical model
    chemical_model = make_chemical_model(chem_ions, chem_gas_mol, tsat, 
                                             pressure[-1].m/1.01325, chem_units = chem_units, pH = 7) 
    output = run_chemcial_model(chemical_model, database_path)
    pp, rho, t_dont_save = get_chemical_properties(output)
    
    #tsat = tsat * u.K
    #print(tsat)
    tsat = np.atleast_1d(tsat)
    rho = rho * u.kg / u.m**3
    density = np.atleast_1d(rho)
    
    for i in depth_diff:
        # Calculate new pressure for this step.
        new_p = pressure[-1] + rho * g * i
        pressure = np.append(pressure, new_p)  # Has units.
        
        # select starting temperature for iterate_chemical_models (reduces number of iterations)
        if 'tsat' in locals():
            t_model = tsat[-1]
        else:
            t_model = 100
            
        
        
        # Calculate new temperature for this step.
        #t = t * u.K
        t = iapws.iapws97._TSat_P(new_p.m) -273.15
        tsat = np.append(tsat, t)
        print(tsat)
        
        # iterate chemical model until partial pressure of model is equal to pressure
        chemical_model = make_chemical_model(chem_ions, chem_gas_mol, t, 
                                             pressure[-1].m/1.01325, chem_units = chem_units, pH = 7) 
        output = run_chemcial_model(chemical_model, database_path)
        pp, rho, t_dont_save = get_chemical_properties(output)
        rho = rho * u.kg / u.m**3
 
        
        # Calculate new density for next step.
        density = np.append(density, rho)
    
    # Finalize units.
    pressure = pressure.to(u_p)
    #tsat = tsat.to(u_t)
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
# calculate hydrostatic bpd
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
# calculate hydrostatic bpd with salinity correction
# ---------------------------------------------------- 

chem_ions_mol = pd.DataFrame([[0.5, 0.5 , 0.5, 0.5, 0.5, 0.5]], columns = ['Na', 'K', 'Ca', 'Cl', 'SO4', 'HCO3'])
chem_gas_mol = pd.DataFrame([[0.01, 0.01, 0.001, 0.001]], columns = ['CH4', 'CO2', 'N2', 'H2'])




sal_curve = bpdc_salinity(depth = df["depth_m"],
                p0=df["pressure_bara"].min(),
                chem_ions = chem_ions_mol,
                chem_gas_mol = chem_gas_mol,
                database_path = r"C:/phreeqc/database/pitzer.dat",
                tolerance = 0.1,
                chem_units = "mol/kgw")

sal_df = pd.DataFrame(sal_curve._asdict())



# ----------------------------------------------------  
# compare results
# ---------------------------------------------------- 

ggplot(df, aes(y = "depth_m")) +\
    geom_line(aes(x = "tsat", y = "depth", colour = "'hydrostatic + iapws'"), data = hyd_df) + \
    geom_line(aes(x = "tsat", y = "depth", colour = "'salinity corrected'"), data = sal_df) + \
    scale_y_reverse() + \
    labs(x = "T [°C]", y = "depth [m]", colour = "legend")



        