# -*- coding: utf-8 -*-
"""
Created on Thu May  15 14:46:45 2025

@author: d4mcbrid
"""

import numpy as np
import pyomo.environ as pyo
import pyomo.dae as pdae
import clc_des__param_L as dp
import clc_gen__param as gp
import clc_ox__param as op
import time
import matplotlib.pyplot as plt
import os
import csv
from pyomo.util.infeasible import log_infeasible_constraints
import pandas as pd
import scipy.integrate
import logging
from X_avg_calc import _X_avg_

# CLC Oxidation Pyomo Model (Multiscale)

# Substitution for scaling: y=z/L (so y=[0,1], z=[0,1]L); Also, t=[0,1], tf is total time (similarly multiplied)

def oxidation_pyo(init, disc, nd, t_solve, cFe2O3):

    [nt, nint, mcp, nr] = disc
    ny = nint*mcp

    # Define initial guesses
    def _C_ic_(m,i,k):
        return init['C'][int(round(i*ny,0)),int(round(k*nt,0))]

    def _T_ic_(m,i,k):
        return init['T'][int(round(i*ny,0)),int(round(k*nt,0))]

    def _Cc_ic_(m,i,j,k):  
        return init['Cc'][int(round(i*ny,0)),int(round(j/m.r.last()*nr,0)),int(round(k*nt,0))]

    def _Tc_ic_(m,i,j,k):
        return init['Tc'][int(round(i*ny,0)),int(round(j/m.r.last()*nr,0)),int(round(k*nt,0))]

    def _X_ic_(m,i,j,k):
        return init['X'][int(round(i*ny,0)),int(round(j/m.r.last()*nr,0)),int(round(k*nt,0))]
    
    def _dCcdr2_init(m,i,j,k):
        return -50


    # Define model (for optimization)
    m = pyo.ConcreteModel()

    ###############################################################################
    # Define parameters
    ###############################################################################

    # Note: reorganizing such that parameters are retrieved from an external source, to assist with consistency if I need to make changes
    # Reactor Design
    m.tf = pyo.Param(initialize=t_solve)     #             # interval duration, s
    m.L = pyo.Param(initialize=dp.L)                           # reactor length, m
    m.LRM = pyo.Param(initialize=3.15)                      # length of red mud section of the bed, m (note the actual RM section is dictated by the discretization elements, having inert section at beginning/end helped with model convergence)
    m.Lq = pyo.Param(initialize=(m.L-m.LRM)/2)              # length of quartz sand section of the bed, m
    m.Di = pyo.Param(initialize=dp.Di)                      # reactor diameter, m
    m.Dp = pyo.Param(initialize=dp.Dp)                      # particle diameter, m
    m.Ac = pyo.Param(initialize=np.pi*m.Di**2/4.0)          # reactor cross-sectional area, m²

    # Reactor Conditions
    m.P = pyo.Param(initialize=op.P)                        # inlet reactor pressure, bar
    m.T_air = pyo.Param(initialize=op.T_air)                # inlet temperature, K
    m.Tref = pyo.Param(initialize=op.T_air)                     # reference temperature for certain calculations
    m.T_in = pyo.Param(initialize=op.T_air)

    # Flexible Parameters
    m.rho = pyo.Param(initialize=op.rho)                    # density of air at 1000 K and 1 bar, kg/m³ -> CHANGE if pressure is different
    m.eb = pyo.Param(initialize=dp.eb)                      # reactor bed porosity
    m.ec = pyo.Param(initialize=dp.ec)                      # particle porosity
    m.dpore = pyo.Param(initialize=dp.dpore)                # pore diameter, m
    m.cFe2O30_1 = pyo.Param(initialize=cFe2O3)                  # Fe2O3 weight fraction in OC (MeO vs support material), kg Fe2O3/kg OC
    m.cFe2O30_2 = pyo.Param(initialize=cFe2O3/gp.mw_Fe2O3*1000)   # Fe2O3 concentration, mol Fe2O3/kg OC

    # General parameters
    m.mw_air = pyo.Param(initialize=gp.mw_air)              # molar mass of air, g/mol
    m.mw_Fe3O4 = pyo.Param(initialize=gp.mw_Fe3O4)                # molar mass of Fe3O4, g/mol
    m.mw_Fe2O3 = pyo.Param(initialize=gp.mw_Fe2O3)              # molar mass of Fe2O3, g/mol
    m.Rg = pyo.Param(initialize=gp.Rg)                      # universal gas constant, J/mol/K
    m.Rg2 = pyo.Param(initialize=gp.Rg2)                    # universal gas constant, g*m²/s²/K/mol
    m.Rgb = pyo.Param(initialize=gp.Rgb)                    # universal gas constant, m³*bar/mol/K

    # Mass transfer parameters
    m.mu = pyo.Param(initialize=op.mu)                      # dynamic viscosity of air at 1000 K, kg/m/s
    m.dij = pyo.Param(initialize=op.dij)                    # binary diffusivity coefficient of oxygen in nitrogen at 1000 K, m²/s
    m.tau = pyo.Param(initialize=dp.tau)                    # particle tortuosity
    m.rho_s = pyo.Param(initialize=op.rho_s)                # particle density (pure solid, not bed), Red mud, kg/m³

    # Heat transfer parameters
    m.cpf = pyo.Param(initialize=op.cpf)                    # heat capacity of air, J/kg/K
    m.cp_f = pyo.Param(initialize=m.cpf*m.mw_air/1000.0)    # heat capacity of air, J/mol/K
    m.l_e0 = pyo.Param(initialize=op.l_e0)                  # static contribution effect of thermal conductivity, J/m/s/K
    m.cp_s = pyo.Param(initialize=op.cp_s)                  # oxygen carrier specific heat capacity, J/kg/K
    m.dH = pyo.Param(initialize=op.dH)                      # heat of reaction for oxidation, J/mol of O2
    m.l_i = pyo.Param(initialize=op.l_i)                    # thermal conductivity of Red mud, W/m/K

    # Reaction Parameters
    m.a0 = pyo.Param(initialize=op.a0)                      # initial specific surface area of OC, m²/kg

    # Define radial, temporal, and axial domains
    m.r = pdae.ContinuousSet(bounds=(0.0,m.Dp/2))
  
    m.y = pdae.ContinuousSet(initialize=[0,0.015,0.03,0.05,0.075,0.1,0.125,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1])
    m.t = pdae.ContinuousSet(initialize=[0,0.005,0.01,0.02,0.03,0.04,0.05,0.08,0.1,0.12,0.15,0.18,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.7,1])
    m.hr = pyo.Param(initialize=m.r.last()/nr)
    m.hy = pyo.Param(initialize=m.y.last()/ny)

    # Define state variables
    m.C = pyo.Var(m.y, m.t, initialize=_C_ic_,domain=pyo.NonNegativeReals)              # Concentration in fluid phase, mol/m³
    m.T = pyo.Var(m.y, m.t, initialize=_T_ic_,domain=pyo.NonNegativeReals)              # Temperature in fluid phase, K
    m.Cc = pyo.Var(m.y, m.r, m.t, initialize=_Cc_ic_,domain=pyo.NonNegativeReals)       # Concentration in the particle, mol/m³
    m.Tc = pyo.Var(m.y, m.r, m.t, initialize=_Tc_ic_,domain=pyo.NonNegativeReals)       # Temperature in the particle, K
    m.X = pyo.Var(m.y, m.r, m.t, initialize=_X_ic_,bounds=(0.01,0.99))       # Solid conversion in the particle

    # Define derivatives
    m.dCdt = pdae.DerivativeVar(m.C, wrt=m.t)
    m.dTdt = pdae.DerivativeVar(m.T, wrt=m.t)
    m.dCcdt = pdae.DerivativeVar(m.Cc, wrt=m.t)
    m.dTcdt = pdae.DerivativeVar(m.Tc, wrt=m.t)
    m.dXdt = pdae.DerivativeVar(m.X, wrt=m.t)
    m.dCdy = pdae.DerivativeVar(m.C, wrt=m.y)
    m.dCdy2 = pdae.DerivativeVar(m.C, wrt=(m.y,m.y))
    m.dTdy = pdae.DerivativeVar(m.T, wrt=m.y)
    m.dTdy2 = pdae.DerivativeVar(m.T, wrt=(m.y,m.y))
    m.dCcdr = pdae.DerivativeVar(m.Cc, wrt=m.r)
    m.dCcdr2 = pdae.DerivativeVar(m.Cc, wrt=(m.r,m.r),initialize=_dCcdr2_init) 
    m.dTcdr = pdae.DerivativeVar(m.Tc, wrt=m.r)
    m.dTcdr2 = pdae.DerivativeVar(m.Tc, wrt=(m.r,m.r))

    # Differentiate system
    disc_fd = pyo.TransformationFactory('dae.finite_difference')
    disc_oc = pyo.TransformationFactory('dae.collocation')

    disc_fd.apply_to(m, wrt=m.t, nfe=nt, scheme='BACKWARD')
    disc_fd.apply_to(m, wrt=m.r, nfe=nr, scheme='CENTRAL')
    disc_oc.apply_to(m, wrt=m.y, nfe=nint, ncp=mcp, scheme='LAGRANGE-RADAU')



    def _dCcdr2_cent(m, i, j, k):
        if j == m.r.first() or j == m.r.last():
            if j == m.r.first():
                return m.dCcdr2[i,j,k] == (m.Cc[i,m.r.next(m.r.next(j)),k]-2*m.Cc[i,m.r.next(j),k]+m.Cc[i,m.r.first(),k])/m.hr**2
            if j == m.r.last():
                return m.dCcdr2[i,j,k] == (m.Cc[i,m.r.last(),k]-2*m.Cc[i,m.r.prev(j),k]+m.Cc[i,m.r.prev(m.r.prev(j)),k])/m.hr**2
        else:
            return pyo.Constraint.Skip
    m.dCcdr2_centDiff = pyo.Constraint(m.y, m.r, m.t, rule=_dCcdr2_cent)

    def _dTcdr2_cent(m, i, j, k):
        if j == m.r.first() or j == m.r.last():
            if j == m.r.first():
                return m.dTcdr2[i,j,k] == (m.Tc[i,m.r.at(3),k]-2*m.Tc[i,m.r.at(2),k]+m.Tc[i,m.r.at(1),k])/m.hr**2
            if j == m.r.last():
                return m.dTcdr2[i,j,k] == (m.Tc[i,m.r.at(-1),k]-2*m.Tc[i,m.r.at(-2),k]+m.Tc[i,m.r.at(-3),k])/m.hr**2
        else:
            return pyo.Constraint.Skip
    m.dTcdr2_centDiff = pyo.Constraint(m.y, m.r, m.t, rule=_dTcdr2_cent)


    ###############################################################################
    # Bulk mass transfer variables
    ###############################################################################

    # Inlet mass flux; kg/m2/s

    m.Gair = pyo.Param(initialize=op.G)
    m.Ginert = pyo.Param(initialize=0)
  
    m.Dmi = pyo.Param(initialize=m.dij)                          # diffusion coefficient of oxygen in the mixture, m²/s
    m.Sc = pyo.Param(initialize=m.mu/(m.rho*m.Dmi))              # Schmidt number of oxygen
    m.av = pyo.Param(initialize=6*(1-m.eb)/m.Dp)                 # external particle surface area per unit volume, 1/m

    # Total inlet mass flux; kg/m2/s
    def _G_(m):
        return m.Gair + m.Ginert
    m.G = pyo.Expression(rule=_G_)

    # Gas velocity, m/s
    def _yi_(m):
        return op.yi*m.Gair/m.G
    m.yi = pyo.Expression(rule=_yi_)

    # Gas velocity, m/s
    def _vel_(m, k):
        return m.G/m.rho
    m.vel = pyo.Expression(m.t, rule=_vel_)

    # Volumetric flowrate (total), m³/s
    def _Q_(m, k):
        return m.G*m.Ac/m.rho
    m.Q = pyo.Expression(m.t, rule=_Q_)

    # Particle Reynolds number based on initial superficial velocity
    def _ReP_(m, k):
        return m.vel[k]*m.Dp*m.rho/((1-m.eb)*m.mu)
    m.ReP = pyo.Expression(m.t, rule=_ReP_)

    # Reynolds number based on initial superficial velocity
    def _Re_(m, k):
        return m.vel[k]*m.Dp*m.rho/m.mu
    m.Re = pyo.Expression(m.t, rule=_Re_)

    # Mass transfer coefficient of oxygen, m/s
    def _kc_(m, k):
        return 0.357*m.Sc**(-2/3)*m.Re[k]**(-0.359)*m.G/(m.rho*m.eb)
    m.kc = pyo.Expression(m.t, rule=_kc_)

    # Effective axial dispersion coefficient of oxygen [combined with Peclet number], m²/s
    def _Dax_(m, k):
        return m.vel[k]*m.Dp*(m.eb/(m.tau*m.ReP[k]*m.Sc)+0.45/(1+7.3/(m.ReP[k]*m.Sc)))
    m.Dax = pyo.Expression(m.t, rule=_Dax_)

    def _yout_(m, k):
        return m.C[m.y.last(),k]/(m.C[m.y.last(),k]+(1-op.yi)*op.P/gp.Rgb/op.T_air/nd[0])
    m.yout = pyo.Expression(m.t, rule=_yout_)    

    ###############################################################################
    # Bulk heat transfer variables
    ###############################################################################

    # Effective axial thermal conductivity, W/m/K
    def _l_ax_(m, k):
        return (m.l_e0+0.7*m.cpf*(m.vel[k]/(1-m.eb)*m.Dp*m.rho))
    m.l_ax = pyo.Expression(m.t, rule=_l_ax_)


    # Heat transfer coefficient between bulk and OC, W/m²/K
    def _hf_(m, k):
        return 1.37*0.357*m.Re[k]**(-0.359)*m.Sc**(-2/3)*m.cpf*m.G/m.eb
    m.hf = pyo.Expression(m.t, rule=_hf_)


    ###############################################################################
    # OC mass transfer variables
    ###############################################################################

    m.Dki = pyo.Expression(expr=m.dpore/3*(8*m.Rg2*m.Tref/3.146/m.mw_air)**0.5)     # Knudson diffusivity coefficient, m²/s
    m.Dei = pyo.Expression(expr=m.Dmi*m.Dki*m.eb/m.tau/(m.Dmi+m.Dki))               # effective diffusivity coefficient, m²/s
    m.k_ox = pyo.Param(initialize=2.282e-7)
    m.order = pyo.Param(initialize=1)

    # Reaction rate
    

    def _rxn_(m, i, j, k):
        if i < m.Lq/m.L or i > (m.Lq+m.LRM)/m.L: # no reaction occurs
            return 0
        else:
            return m.a0*m.k_ox*m.Cc[i,j,k]**m.order*m.cFe2O30_1*2*(0.99-m.X[i,j,k])*(-pyo.log(1.000000000000000001-m.X[i,j,k]))**0.5
    m.r_ox = pyo.Expression(m.y, m.r, m.t, rule=_rxn_)

    # Previous reduction stage conversion calculation
    zs_int = np.array([ 0.02    , 0.03    , 0.036667,
           0.05    , 0.058333, 0.075   , 0.083333, 0.1     , 0.108333,
           0.125   , 0.133333, 0.15    , 0.166667, 0.2     , 0.216667,
           0.25    , 0.266667, 0.3     , 0.316667, 0.35    , 0.358333,
           0.375   , 0.383333, 0.4     , 0.416667, 0.45    , 0.466667,
           0.5     , 0.516667, 0.55    , 0.566667, 0.6     , 0.616667,
           0.65    , 0.666667, 0.7     , 0.716667, 0.75    , 0.766667,
           0.8     , 0.816667, 0.85    , 0.866667, 0.9     , 0.916667,
           0.95    , 0.966667    ]) # axial discrete elements

    r_int = np.array([0,0.000833,0.001667,0.0025]) # radial discrete elements

    X_red = _X_avg_(filename_X,zs_int,r_int,dp.Dp/2,dp.L) # average conversion from reduction stage (used to estimate cFe3O4  = cFe2O3 (mol/kg)*X_red)

    # Conversion rate for OC particle
    def _conv_(m, i, j, k):
        return m.dXdt[i,j,k]*m.cFe2O30_2*X_red == m.tf*4*m.r_ox[i,j,k]*nd[0]
    m.x_conv = pyo.Constraint(m.y, m.r, m.t, rule=_conv_)


    ###############################################################################
    # OC heat transfer variables
    ###############################################################################

    m.l_s = pyo.Expression(expr=m.l_i)     # weighted thermal conductivity of OC, W/m/K


    ###############################################################################
    # PBR mass balance
    ###############################################################################

    def _PBR_mass(m, i, k):
        return 1/m.tf*m.eb*m.dCdt[i,k] + m.Q[k]/m.Ac*(1/m.L)*m.dCdy[i,k] == m.Dax[k]*(1/(m.L**2))*m.dCdy2[i,k] + m.kc[k]*m.av*(m.Cc[i,m.r.last(),k]-m.C[i,k])
    m.PBR_mass = pyo.Constraint(m.y, m.t, rule=_PBR_mass)

    # PBR mass balance BC; feed stream condition at inlet
    def _BC1_(m, k):
        if k == m.t.first():
            return (m.C[0,m.t.first()] - trace*op.P/gp.Rgb/op.T_air/nd[0])**2 <= m.C[0,m.t.first()]*1e-4
        else:
            return m.C[0,k] == m.yi*m.P/(m.Rgb*m.T_in)/nd[0]
    m.BC1 = pyo.Constraint(m.t, rule=_BC1_)
    # PBR mass balance BC; mass insulation at reactor outlet
    def _BC2_(m, k):
        return (m.C[m.y.last(),k] - m.C[m.y.at(-2),k])**2 <= 1e-2*m.C[m.y.at(-2),k]
    m.BC2 = pyo.Constraint(m.t, rule=_BC2_)


    ###############################################################################
    # PBR energy balance [possibly move up apparently]
    ###############################################################################

    def _PBR_energy(m, i, k):
        return m.eb*m.cp_f*m.P/(m.Rgb*m.T_air)*(1/m.tf)*m.dTdt[i,k] + m.cp_f*m.Q[k]*m.P/(m.Rgb*m.T_air)/m.Ac*(1/m.L)*m.dTdy[i,k] == m.l_ax[k]*((1/m.L)**2)*m.dTdy2[i,k] + m.hf[k]*m.av*(m.Tc[i,m.r.last(),k]-m.T[i,k])
    m.PBR_energy = pyo.Constraint(m.y, m.t, rule=_PBR_energy)

    # PBR energy balance BC; heat feed stream conditions at inlet
    def _BC3_(m, k):
      return (m.T_air/nd[1] - m.T[m.y.at(1),k])**2 <= 1e-6
    m.BC3 = pyo.Constraint(m.t, rule=_BC3_)


    # PBR energy balance BC; heat insulation condition at reactor's outlet stream
    def _BC4_(m, k):
        return (m.T[m.y.last(),k] - m.T[m.y.at(-2),k])**2 <= 1e-2*m.T[m.y.at(-2),k]
    m.BC4 = pyo.Constraint(m.t, rule=_BC4_)


    ###############################################################################
    # OC Particle mass balance
    ###############################################################################

    def _OC_mass(m, i, j, k):
        if j == m.r.first():
            return m.ec*m.dCcdt[i,j,k] == m.tf * (3*m.Dei*(m.dCcdr2[i,j,k]) - m.ec*m.rho_s*m.r_ox[i,j,k])
        return m.ec*j**2*m.dCcdt[i,j,k] == m.tf * (m.Dei*(2*j*m.dCcdr[i,j,k]+m.dCcdr2[i,j,k]*j**2) - j**2*m.ec*m.rho_s*m.r_ox[i,j,k])
    m.OC_mass = pyo.Constraint(m.y, m.r, m.t, rule=_OC_mass)

    # OC mass balance BC; mass insulation
    def _BC5_(m, i, k):
        return m.dCcdr[i,m.r.first(),k] == 0
    m.BC5 = pyo.Constraint(m.y, m.t, rule=_BC5_)

    # OC mass balance BC; mass transfer between reactor bulk phase and OC particle
    def _BC6_(m, i, k):
        return -m.Dei*m.dCcdr[i,m.r.last(),k] == m.kc[k]*(m.Cc[i,m.r.last(),k]-m.C[i,k])
    m.BC6 = pyo.Constraint(m.y, m.t, rule=_BC6_)


    ###############################################################################
    # OC Particle energy balance
    ###############################################################################

    def _OC_energy(m, i, j, k):
        if j == m.r.first():
            return ((1-m.ec)*m.rho_s*m.cp_s + m.ec*m.P/(m.Rgb*m.T_air)*m.cp_f)*m.dTcdt[i,j,k] == m.tf * (3*m.l_s*m.dTcdr2[i,j,k] - m.ec*m.rho_s*m.r_ox[i,j,k]*m.dH*nd[0]/nd[1])
        return ((1-m.ec)*m.rho_s*m.cp_s + m.ec*m.P/(m.Rgb*m.T_air)*m.cp_f)*j**2*m.dTcdt[i,j,k] == m.tf * (m.l_s*(2*j*m.dTcdr[i,j,k]+m.dTcdr2[i,j,k]*j**2) - j**2*m.ec*m.rho_s*m.r_ox[i,j,k]*m.dH*nd[0]/nd[1])
    m.OC_energy = pyo.Constraint(m.y, m.r, m.t, rule=_OC_energy)

    # OC energy balance BC; heat insulation
    def _BC7_(m, i, k):
        return m.dTcdr[i,m.r.first(),k] == 0
    m.BC7 = pyo.Constraint(m.y, m.t, rule=_BC7_)

    # OC energy balance BC; heat transfer between reactor bulk phase and OC particle
    def _BC8_(m, i, k):
        return -m.l_s*m.dTcdr[i,m.r.last(),k] == m.hf[k]*(m.Tc[i,m.r.last(),k]-m.T[i,k])
    m.BC8 = pyo.Constraint(m.y, m.t, rule=_BC8_)
    ###########################################################################
    # Pilot scale up characteristics [for now left as just calculating the values - convert to constraints for optimization]
    ###########################################################################
    # pressure drop
    def _Ergun_(m,k):
        return 3.1*(150*(1-m.eb)**2/(m.eb**3)*m.mu*m.vel[k]/(m.Dp**2)+1.75*(1-m.eb)/(m.eb**3)*m.rho*m.vel[k]**2/(m.Dp))
     
    m.dP = pyo.Expression(m.t,rule=_Ergun_) # pressure drop, Pa
    
    m.umf = pyo.Var(initialize=0.5,domain=pyo.NonNegativeReals) # minimum fluidization velocity
    def _umf_(m):
        return 1305.556*m.rho*m.umf*m.Dp/m.mu + 30.789*(m.rho*m.umf*m.Dp/m.mu)**1.657 + 19.474*(m.rho*m.umf*m.Dp/m.mu)**2-m.rho*m.Dp**3*(m.rho_s*(1-m.ec)-m.rho)*9.81/(m.mu**2) == 0
    m.u_mf = pyo.Constraint(rule=_umf_)
    
    m.LD = pyo.Expression(expr=m.L/m.Di)
   
    def Q_out(m):# Average heat extracted via outlet gas (rough approximation using left reimann sums)
        value = 0
        for t in range(len(m.t)):
            if t == 0:
                value += 0
            else:
                value += (m.t.at(t+1)-m.t.at(t))*m.Q[m.t.at(t)]*(m.C[m.y.last(),m.t.at(t)]*nd[0]+(1-op.yi)*op.P/gp.Rgb/op.T_air)*m.cp_f*(m.T[m.y.last(),m.t.at(t)]*nd[1]-op.T_air)
        return value
    def Q_in(m): # Average heat remaining in the packed bed at the end of oxidation stage (rough approximation using left reimann sums)
        value = 0
        for z in range(len(m.y)):
            if z == 0:
                value += 0
            else:
                value += m.L*(m.y.at(z+1)-m.y.at(z))*m.rho_s*(1-m.eb)*(1-m.ec)*m.Ac*m.cp_s*(m.Tc[m.y.at(z),m.r.last(),m.t.last()]*nd[1]-op.T_air)
        return value/m.tf
    m.Q_in = pyo.Expression(rule=Q_in)
    m.Q_out = pyo.Expression(rule=Q_out)
    
    # integrated heats (more accurate)
    
    def Q_out_int(m,k):  # Average heat extracted via outlet gas
       return  m.tf*m.Q[k]*(m.C[m.y.last(),k]*nd[0]+(1-m.yi)*op.P/gp.Rgb/op.T_air)*m.cp_f*(m.T[m.y.last(),k]*nd[1]-op.T_air)
   
    def Q_in_int(m,i):# Average heat remaining in the packed bed at the end of oxidation stage
        return m.L*m.rho_s*(1-m.eb)*(1-m.ec)*m.Ac*m.cp_s*(m.Tc[i,m.r.last(),m.t.last()]*nd[1]-op.T_air)
    
    m.Q_out_int = pdae.integral.Integral(m.t,wrt=m.t,rule=Q_out_int)
    m.Q_in_int = pdae.integral.Integral(m.y,wrt=m.y,rule=Q_in_int)

    return m

if __name__ == "__main__":
     # Import Yin's data for model fitting/validation
    wdir = os.path.dirname(os.path.realpath(__file__))


    tfreq = 20                  # Frequency of time points in discretization, s
    
    
    tf = 1400                   # Reduction stage simulation time, s
    ncon = 1      # Number of flux steps (control actions) per interval
    
    
    cFe2O3 = 0.4928#                                                     # Fe2O3 concentration (mass fraction)
    

    
    # Define nondimensionalizing coefficients for concentration and temperature
    Cnd = op.yi*op.P/(gp.Rgb*op.T_air)
    Tnd = op.T_air
    nd = [Cnd,Tnd]
     
    # Discretization
 
    nr = 3                          # Number of radial nodes (particle)
    nint = 25                     # Number of axial elements (bulk reactor)
    mcp = 2                      # Total number of orthogonal collocation points sent to pyomo
    nt = 21        # Number of temporal nodes
    ntcon = int(nt/ncon)            # Number of nodes in each step of the control profile
    ny = nint*mcp                   # Total number of axial points
    l_t = tf/nt                 # Length of each time node which must be solved by solve_ivp
    disc = [nt, nint, mcp, nr]

    trace = 1e-4
    yi = op.yi #inlet mole fractions
    
    # Initial Conditions
    C0 = trace*op.P/gp.Rgb/op.T_air    # initial fluid-phase concentrations of compounds, mol/m³
    T0 = op.T_air                                          # initial temperature at reactor inlet, K
    X0 = 0.5                                          # initialization of conversion of Fe3O4 to Fe2O3
    
    # Define initial guesses (and initial conditions for t=0)
    init = {'C': np.full([ny+1,nt+1],0,dtype=float),
            'T': np.full([ny+1,nt+1],0,dtype=float),
            'Xp': np.full([ny+1,nt+1],0,dtype=float),
            'Cc': np.full([ny+1,nr+1,nt+1],0,dtype=float),
            'Tc': np.full([ny+1,nr+1,nt+1],0,dtype=float),
            'X': np.full([ny+1,nr+1,nt+1],0,dtype=float)}
    
    # Define matrices to collect data
    plantsolves = []
    ts = np.full([nt+1],0,dtype=float)
    tTs = ts[nt:]
    Cs = np.full([ny+1,len(ts)],0,dtype=float)
    youts = np.full([len(ts)],0,dtype=float)
    youtsdry = np.full([len(ts)],0,dtype=float)
    Ts = np.full([ny+1,len(ts)],0,dtype=float)
    Ccs = np.full([ny+1,nr+1,len(ts)],0,dtype=float)
    Tcs = np.full([ny+1,nr+1,len(ts)],0,dtype=float)
    Xs = np.full([ny+1,nr+1,len(ts)],0,dtype=float)
    Gs = np.full([len(ts)],0,dtype=float)
    Ssets = np.full([len(tTs)],0,dtype=float)
    SCO2s = np.full([len(ts)],0,dtype=float)
    Gsets = np.full([len(tTs)],0,dtype=float)
    GCO2s = np.full([len(ts)],0,dtype=float)
    GCH4s = np.full([len(ts)],0,dtype=float)
    r_O2 = np.full([ny+1,len(ts)],0,dtype=float)
    zs = np.full([ny+1],0,dtype=float)
    tps = [k*l_t for k in range(nt,int((tf)/l_t)+1)]
    

    plantopt = pyo.SolverFactory('ipopt')#, executable='C:\Jaime\MSYS2\home\ja2mccor\ipopt')#pyo.SolverFactory('ipopt', executable='C:\software\cygwin\home\ipopt')#'C:\Dana\ipopt')
    plantopt.options['halt_on_ampl_error'] = 'yes'
    plantopt.options['linear_solver'] = 'ma97'

    plantopt.options['ma97_u'] = 0.1

    
    # Define initial guesses
    for k in range(nt+1):
        for i in range(ny+1):
            init['C'][i,k] = C0/nd[0]
            init['T'][i,k] = T0/nd[1]
            init['Xp'][i,k] = X0
            for j in range(nr+1):
                init['Cc'][i,j,k] = C0/nd[0]
                init['Tc'][i,j,k] = T0/nd[1]
                init['X'][i,j,k] = X0
    
    tt = 0
    filename_X = 'XPilot_red_base.xlsx'
    Xinit  = pd.read_excel(filename_X, index_col=[0,1,2]) # 0: y, 1: r, 2: t
    # calculate the average conversion of the OC at the end of reduction
    
   
    X_red_old = 0.5 #placeholder

    plant = oxidation_pyo(init, disc, nd, tf, cFe2O3)

    X_red_avg = np.full([ny+1,nr+1],0,dtype=float)
    # Fix initial conditions for plant
    for i in range(ny+1):
        plant.C[plant.y.at(i+1),plant.t.first()].fix(init['C'][i,0])
        plant.T[plant.y.at(i+1),plant.t.first()].fix(init['T'][i,0])
        for j in range(nr+1):
            plant.Cc[plant.y.at(i+1),plant.r.at(j+1),plant.t.first()].fix(init['Cc'][i,j,0])
            plant.Tc[plant.y.at(i+1),plant.r.at(j+1),plant.t.first()].fix(init['Tc'][i,j,0])
            if Xinit.loc[plant.y.at(i+1)].loc[plant.r.at(j+1)].loc[1].val > 0.99:
                plant.X[plant.y.at(i+1),plant.r.at(j+1),plant.t.first()].fix(0.01)
            else:
                plant.X[plant.y.at(i+1),plant.r.at(j+1),plant.t.first()].fix(1-Xinit.loc[plant.y.at(i+1)].loc[plant.r.at(j+1)].loc[1].val)
            X_red_avg[i,j] = Xinit.loc[plant.y.at(i+1)].loc[plant.r.at(j+1)].loc[1].val


    zs_int = np.array([0.02    , 0.03    , 0.036667,
           0.05    , 0.058333, 0.075   , 0.083333, 0.1     , 0.108333,
           0.125   , 0.133333, 0.15    , 0.166667, 0.2     , 0.216667,
           0.25    , 0.266667, 0.3     , 0.316667, 0.35    , 0.358333,
           0.375   , 0.383333, 0.4     , 0.416667, 0.45    , 0.466667,
           0.5     , 0.516667, 0.55    , 0.566667, 0.6     , 0.616667,
           0.65    , 0.666667, 0.7     , 0.716667, 0.75    , 0.766667,
           0.8     , 0.816667, 0.85    , 0.866667, 0.9     , 0.916667,
           0.95    , 0.966667     ])
    r_int = np.array([0,0.000833,0.001667,0.0025])

    X_red = _X_avg_(filename_X,zs_int,r_int,dp.Dp/2,3.1)

    print(X_red)
    start = time.time()
    # Solve plant model



#    
    plantresults = plantopt.solve(plant,tee=True)
    plantsolves += [plantresults['Solver']]

    # collect results
    for k in range(nt+1):

       ts[k] = plant.t.at(k+1)*pyo.value(plant.tf)
       youts[k] = plant.yout[plant.t.at(k+1)]()

       for i in range(ny+1):
           zs[i] = plant.y.at(i+1)*dp.L
           Cs[i,k] = plant.C[plant.y.at(i+1),plant.t.at(k+1)].value*nd[0]
           Ts[i,k] = plant.T[plant.y.at(i+1),plant.t.at(k+1)].value*nd[1]
           r_O2[i,k] = pyo.value(plant.r_ox[plant.y.at(i+1),plant.r.last(),plant.t.at(k+1)])
           for j in range(nr+1):
               Ccs[i,j,k] = plant.Cc[plant.y.at(i+1),plant.r.at(j+1),plant.t.at(k+1)].value*nd[0]
               Tcs[i,j,k] = plant.Tc[plant.y.at(i+1),plant.r.at(j+1),plant.t.at(k+1)].value*nd[1]
               Xs[i,j,k] = plant.X[plant.y.at(i+1),plant.r.at(j+1),plant.t.at(k+1)].value
               
    l1_CO2 = youts[:]
    l2 = ts
    l1_X = Xs[:,-1,-1]
    l3 = zs
    df_l = pd.DataFrame({'length':l3,'final conversion':l1_X})
    df_CO2 = pd.DataFrame({' fraction':l1_CO2,'time':l2})
    df_CO2.to_excel('results_time_ox.xlsx')
    df_l.to_excel('results_length_ox.xlsx')
    foldername = "Pilot_ox_base"
    
    zs_int = np.array([0.02    , 0.03    , 0.036667,
           0.05    , 0.058333, 0.075   , 0.083333, 0.1     , 0.108333,
           0.125   , 0.133333, 0.15    , 0.166667, 0.2     , 0.216667,
           0.25    , 0.266667, 0.3     , 0.316667, 0.35    , 0.358333,
           0.375   , 0.383333, 0.4     , 0.416667, 0.45    , 0.466667,
           0.5     , 0.516667, 0.55    , 0.566667, 0.6     , 0.616667,
           0.65    , 0.666667, 0.7     , 0.716667, 0.75    , 0.766667,
           0.8     , 0.816667, 0.85    , 0.866667, 0.9     , 0.916667,
           0.95    , 0.966667     ])

    r_int = np.array([0,0.000833,0.001667,0.0025])
    # Calculating average converison after oxidation
    filename_X_ox = 'X_ox.xlsx'
    ii = 0
    writer = pd.ExcelWriter(filename_X_ox)
    df = pd.DataFrame(columns=['i','j','k','val'])
    for i in plant.X.keys():
        val = pyo.value(plant.X[i])
        #print("X[{}] = {}".format(i,val))
        df.loc[ii] = [i[0],i[1],i[2],val]
        ii += 1
    df.to_excel(writer,sheet_name='X',index=None)
    writer.close()
    X_ox = _X_avg_(filename_X_ox,zs_int,r_int,dp.Dp/2,3.1)
    
    data = {'value':[pyo.value(plant.vel[plant.t.last()]),pyo.value(plant.umf*plant.Dp*plant.rho/plant.mu),round(pyo.value(9.81*plant.Dp**3*plant.rho*(plant.rho_s-plant.rho)/(plant.mu**2))),pyo.value(plant.LD),tf/60,plant.Dp/(plant.Di/2),pyo.value(plant.dP[plant.t.last()])/101325*100,pyo.value(plant.Di/plant.Dp)],'restriction':[pyo.value(plant.umf),'<350','<5000000','2-6','15 min','0.04-0.25','<8%','>2']}
    print(pd.DataFrame(data,index=[ 'Gas velocity','Remf','Ar','L/D','Cycle time','dp/R0','Pressure drop %','D/dp']))
    print('min temp: ',np.min(Ts))  
    print('max temp: ',np.max(Ts)) 
    print('min conc: ',np.min(Cs))  
    print('max conc: ',np.max(Cs)) 
    print('average final conversion: ',X_ox)#np.average(Xs[5:-1,-1,-1]))
    
    save = False# True
    x = zs[3:-1]-zs[3]      # RM section (change indices if discretization or LRM changes)
    y = Xs[3:-1,-1,-1]       # conversion in the RM section
    I2 = scipy.integrate.simpson(y=y,x=x) 
    if not os.path.exists(foldername):
        os.makedirs(foldername)
    if save:
         ts_save = np.asarray(ts)
         zs_save = np.array(zs)
         np.save(foldername+"\\times.npy",ts_save)
         np.save(foldername+"\\length.npy",zs_save)
         np.save(foldername+"\\Cs.npy",Cs)
         np.save(foldername+"\\Ts.npy",Ts)
         np.save(foldername+"\\Ccs.npy",Ccs)
         np.save(foldername+"\\Tcs.npy",Tcs)
         np.save(foldername+"\\Xs.npy",Xs)
         np.save(foldername+'\\youts',youts)
         np.save(foldername+'\\youtsdry',youtsdry)
         with open(foldername+'\\output_file.txt', 'w') as f:
             f.write(str(pd.DataFrame(data,index=['Gas velocity','Remf','Ar','L/D','Cycle time','dp/R0','Pressure drop %','D/dp'])))
             f.write("\nAverage OC conversion: "+str(X_ox))
             f.write('\nMin/Max conc: '+str(np.min(Cs))+'/'+str(np.max(Cs)))
             f.write('\nMin/Max temp: '+str(np.min(Ts))+'/'+str(np.max(Ts)))
             f.write('\nExtracted heat (Reimann/trapezoid), J: '+str(pyo.value(plant.Q_out*plant.tf))+'/'+str(pyo.value(plant.Q_out_int)))
             f.write('\nAccumulated heat (Reimann/trapezoid), J: '+str(pyo.value(plant.Q_in*plant.tf))+'/'+str(pyo.value(plant.Q_in_int)))
             f.write('\ntotal heat (Reimann/trapezoid): '+str((pyo.value(plant.Q_out*plant.tf) + pyo.value(plant.Q_in*plant.tf))/1000) +'/'+str((pyo.value(plant.Q_out_int)+pyo.value(plant.Q_in_int))/1000)+' kJ')
             f.write('\n cFe2O3: '+str(cFe2O3))
             f.write('\nDesign Variables:')
             f.write('\n  tf: '+str(pyo.value(plant.tf)))
             f.write('\n  yO2: '+str(pyo.value(plant.yi)))
             f.write('\n  G (kg/m2/s): '+str(pyo.value(plant.G)))
             f.write('\n  Di: '+str(pyo.value(plant.Di)))
             f.write('\n  L: 3.1 (rounded)')
             

    
             
    plt.figure()
    plt.plot([0.,3.25],[Ts[0,1],Ts[0,1]],'--r')
    plt.plot(zs,Ts[:,-1])
    plt.xlabel('Reactor length (m)')
    plt.ylabel('Final Temperature')
    
    plt.figure()
    plt.plot(ts,Cs[-1,:])                                                         
    plt.xlabel('Time (s)')
    plt.ylabel('outlet O2 concentration (mol/m3)')
    plt.show()
    
    plt.figure()
    plt.plot(zs,Cs[:,-1])
    plt.xlabel('Reactor length (m)')
    plt.ylabel('O2 concentration (mol/m3)')
    plt.show()
    
    plt.figure()
    plt.plot(zs,Xs[:,-1,-1],label='Final')
    plt.plot(zs,Xs[:,-1,0],label='Initial')
    plt.xlabel('Reactor length (m)')
    plt.ylabel('Surface conversion')
    plt.legend()
    
    plt.figure()
    plt.plot(ts,Xs[4,-1,:],label='entrance')
    plt.plot(ts,Xs[22,-1,:],label='middle')
    plt.plot(ts,Xs[-4,-1,:],label='exit')
    plt.xlabel('Time (s)')
    plt.ylabel('Surface conversion')
    plt.legend()
    plt.show()
    
    plt.figure()
    plt.plot([0,3.25],[Cs[0,1],Cs[0,1]],'--r',label='inlet O2')
    plt.plot(zs,Cs[:,-1],label='-1')
    plt.plot(zs,Cs[:,-2],label='-2')
    plt.plot(zs,Cs[:,-3],label='-3')
    plt.plot(zs,Cs[:,-5],label='-5')
    plt.plot(zs,Cs[:,-6],label='-6')
    plt.plot(zs,Cs[:,-7],label='-7')
    plt.plot(zs,Cs[:,-8],label='-8')
    plt.plot(zs,Cs[:,-9],label='-9')
    plt.plot(zs,Cs[:,-10],label='-10')
    plt.plot(zs,Cs[:,-11],label='-11')
    plt.plot(zs,Cs[:,-12],label='-12')
    plt.plot(zs,Cs[:,-13],label='-13')
    plt.plot(zs,Cs[:,-14],label='-14')
    plt.xlabel('Reactor length (m)')
    plt.ylabel('O2 concentration (mol/m3)')
    plt.legend()
    plt.show()
    logging.basicConfig(level=logging.INFO)
    

