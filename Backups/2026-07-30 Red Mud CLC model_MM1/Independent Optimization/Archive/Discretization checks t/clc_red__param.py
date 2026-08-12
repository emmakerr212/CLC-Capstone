# -*- coding: utf-8 -*-
"""
Created on Mon Mar  17 09:40:48 2025

@author: d4mcbrid
"""
import clc_des__param_L as dp
# Reactor Conditions
m_in = 6.94e-7*2*2*1.5*10000*1.565*1.75*1*1#2*2*2*2*2*2#2.6975e-7*1#1.35e-7#6.95e-7#2.6975e-7                                          # inlet mass flowrate, kg/s, calculated at reactor conditions
G = m_in/dp.Ac#0.0001767#m_in/dp.Ac#0.00017671458676442585                      # inlet mass flux, kg/m²/s
P = 1.01325                                            # inlet reactor pressure, bar
T_in = 850+273.                                     # inlet temperature, K
trace = 1e-4                                        # trace component definition to avoid numerical errors
yis = [0.085, 0.14, 0.335, 0.14, trace, 0.3-trace]#[0.07,0.12,0.28,0.14,trace,0.65-trace]#[0.085, 0.14, 0.335, 0.14, trace, 0.3-trace]#[0.07,0.12,0.28,0.14,trace,0.65-trace]#[0.1, 0.16, 0.385, 0.14, trace, 0.215-trace]#[0.07,0.12,0.28,0.14,trace,0.65-trace]#[0.1, 0.16, 0.385, 0.14, trace, 0.215-trace]#[0.085, 0.14, 0.335, 0.14, trace, 0.3-trace]#[0.085, 0.14, 0.335, 0.14, trace, 0.3-trace]#[0.05,0.1,0.25,0.2,0.1,0.3]#[0.2, 0.25, 0.335, 0.14, trace, 0.075-trace]#[0.085, 0.14, 0.335, 0.14, trace, 0.3-trace]#[0.085, 0.14, 0.335, 0.14, trace, 0.3-trace]#[0.2, 0.25, 0.335, 0.14, trace, 0.075-trace]  #[0.085, 0.14, 0.335, 0.14, trace, 0.3-trace]            # inlet mole fractions
#[CH4, H2, CO, CO2, H2O, inert]

# Mass transfer parameters
rho = 0.278#0.1079                                          # density of mixture, kg/m³
mu = 4.194e-5                                      # dynamic viscosity of mixture, kg/m/s, McCabe, Smith, and Harriott
dij2 = 3.61e-5#3.1278e-5                                    # diffusivity of CH4-H2O, m²/s
tau = 3.0                                           # particle tortuosity
rho_s = 2800                                        # particle density (pure solid, not bed), red mud, kg/m³

# Heat transfer parameters
cp_f = 39.91                                        # heat capacity of mixture, J/mol/K
cpc = cp_f                                          # heat capacity of mixture in the particle, J/mol/K
l_e0 = 0.01                                         # static contribution effect of thermal conductivity, J/m/s/K
cp_s = 1230                                         # oxygen carrier specific heat capacity, J/kg/K
l_i1 = 5                                            # thermal conductivity of Fe2O3, W/m/K (Garcia Labiano, 2005)
l_i2 = 2                                            # thermal conductivity of SiO2, W/m/K
l_i3 = 1                                            # thermal conductivity of Na2O, W/m/K (highly uncertain)
l_i4 = 4                                            # thermal conductivity of TiO2, W/m/K
l_i5 = 7                                            # thermal conductivity of Al2O3, W/m/K
l_i = 3.30                                          # thermal conductivity of the OC (combined)


Dei = [8.57e-07, 3.01e-06, 8.24e-07, 6.45e-07,1.02e-06,8.23e-07]#[9.28e-7,3.2e-6,8.92e-7,6.98e-7,1.1e-6,8.91e-7]#[8.45e-7,3.09e-6,8.67e-7,6.98e-7,1.1e-6,8.87e-7]# [1.16e-7,8.27e-7,3.06e-6,6.89e-7,1.09e-6,8.85e-7]#[1.57e-7, 8.27e-7, 3.06e-6, 5.03e-7, 6.42e-7, 8.36e-07]   # effective diffusivity coefficients for each species, m²/s


# Reaction rate parameters
# R1 (m1): CH4 + 12Fe2O3 -> CO2 + 2H2O + 8Fe3O4
# R4 (smr): CH4 + H2O <-> 3H2 + CO
# R5 (wgs): CO + H2O <-> H2 + CO2

a0 = 4620                                            # initial specific surface area of OC, m²/kg
n_AE = 2                                         # reaction rate constant

dH1 = 120530                                        # heat of reaction, J/mol
dH2 = -47697                                        # heat of reaction, J/mol
dH3 = -14495                                        # heat of reaction, J/mol