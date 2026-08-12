# -*- coding: utf-8 -*-
"""
Created on Fri Jun  3 14:08:48 2022

@author: kmtoffol
"""
import clc_des__param_L as dp

# Reactor Conditions
m_in = 1*0.9*1.574e-6*2*2*2*4*2*10*2*2*2*4*2*2.05*2*4*2*1*0.5*0.75*0.75*0.5*4*1.5*1*0.5*0.5                # inlet mass flowrate of air, kg/s
G = m_in/dp.Ac                    # inlet mass flux, kg/m²/s
P = 1.01325                   # inlet reactor pressure, bar
T_air = 1000           # inlet temperature, K
yi = 0.18        # inlet fraction of oxygen in the air stream

# Mass transfer parameters
rho = 0.3529                 # density of air at 1000 K and 1.01325 bar, kg/m³ -> CHANGE if pressure is different
mu = 4.153e-5               # dynamic viscosity of  air at 1000 K (not much impact from pressure), kg/m/s
dij = 1.61e-4              # binary diffusivity coefficient of oxygen in nitrogen at 1000 K and 1 atm, m²/s, Hirschfelder eqn
rho_s = 2800                # particle density (pure solid, not bed), red mud, kg/m³

# Heat transfer parameters
cpf = 1144.805            # heat capacity of air, J/kg/K
cp_f = 33.15              # heat capacity of air, J/mol/K
l_e0 = 0.01                 # static contribution effect of thermal conductivity, J/m/s/K
cp_s = 1230                 # oxygen carrier specific heat capacity, J/kg/K
dH = -467658.2          # heat of reaction for oxidation at 1000 K, J/mol of O2
l_i = 3.30                   # thermal conductivity of solid, W/m/K (weighted based on info from Garcia Labiano, 2005)




# Reaction Parameters
a0 = 4620                    # initial specific surface area of OC, m²/kg
