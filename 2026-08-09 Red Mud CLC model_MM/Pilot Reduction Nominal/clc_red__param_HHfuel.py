# -*- coding: utf-8 -*-
"""
Created on Mon Mar  17 09:40:48 2025

@author: d4mcbrid
"""
import clc_des__param_L as dp
# Reactor Conditions
m_in = 6.94e-7*2*2*1.5*10000*1.565*1.75*1*1                                        # inlet mass flowrate, kg/s, calculated at reactor conditions
G = m_in/dp.Ac                    # inlet mass flux, kg/m²/s
P = 1.01325                                            # inlet reactor pressure, bar
T_in = 850+273.                                     # inlet temperature, K
trace = 1e-4                                        # trace component definition to avoid numerical errors
yis = [0.11, 0.18, 0.44, 0.14, trace,0.1299]#      # inlet mole fractions
#[CH4, H2, CO, CO2, H2O, inert]

# Mass transfer parameters
rho = 0.262                                         # density of mixture, kg/m³
mu = 4.02e-5                                      # dynamic viscosity of mixture, kg/m/s, McCabe, Smith, and Harriott
dij2 = 3.61e-5                                    # diffusivity of CH4-H2O, m²/s
tau = 3.0                                           # particle tortuosity
rho_s = 2800                                        # particle density (pure solid, not bed), red mud, kg/m³

# Heat transfer parameters
cp_f = 40.93#39.91                                        # heat capacity of mixture, J/mol/K
cpc = cp_f                                          # heat capacity of mixture in the particle, J/mol/K
l_e0 = 0.01                                         # static contribution effect of thermal conductivity, J/m/s/K
cp_s = 1230                                         # oxygen carrier specific heat capacity, J/kg/K
l_i1 = 5                                            # thermal conductivity of Fe2O3, W/m/K (Garcia Labiano, 2005)
l_i2 = 2                                            # thermal conductivity of SiO2, W/m/K
l_i3 = 1                                            # thermal conductivity of Na2O, W/m/K (highly uncertain)
l_i4 = 4                                            # thermal conductivity of TiO2, W/m/K
l_i5 = 7                                            # thermal conductivity of Al2O3, W/m/K
l_i = 3.30                                          # thermal conductivity of the OC (combined)


Dei = [7.97e-07, 3.02e-06, 8.25e-07, 6.4e-07,1.02e-06,8.21e-07]# effective diffusivity coefficients for each species, m²/s


# Reaction rate parameters


a0 = 4620                                            # initial specific surface area of OC, m²/kg
n_AE = 2                                         # reaction rate constant

dH1 = 120530                                        # heat of reaction, J/mol
dH2 = -47697                                        # heat of reaction, J/mol
dH3 = -14495                                        # heat of reaction, J/mol