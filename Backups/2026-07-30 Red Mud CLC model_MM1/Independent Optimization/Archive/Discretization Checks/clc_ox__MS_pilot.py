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

def oxidation_pyo(init, disc, nd, cFe2O3, filename_X, t_red):

    [nt, nint, mcp, nr] = disc
    ny = nint*mcp
    
    # filename = 'C_ox_init.xlsx'
    # C = pd.read_excel(filename, index_col=[0,1]) # 0: y, 1: t
    # filename = 'Cc_ox_init.xlsx'
    # Cc = pd.read_excel(filename, index_col=[0,1,2]) # 0:y, 1: r, 2: t
    # filename = 'T_ox_init.xlsx'
    # T = pd.read_excel(filename, index_col=[0,1]) # 0: y, 1: t
    # filename = 'Tc_ox_init.xlsx'
    # Tc  = pd.read_excel(filename, index_col=[0,1,2]) # 0: y, 1: r, 2: t
    # filename = 'X_ox_init.xlsx'
    # X  = pd.read_excel(filename, index_col=[0,1,2]) # 0: y, 1: r, 2: t

    # Define initial guesses
    def _C_ic_(m,i,k):
        return init['C'][int(round(i*ny,0)),int(round(k*nt,0))]
       # return C.loc[i].loc[k].val

    def _T_ic_(m,i,k):
       return init['T'][int(round(i*ny,0)),int(round(k*nt,0))]
       #return T.loc[i].loc[k].val

    def _Cc_ic_(m,i,j,k):  
       return init['Cc'][int(round(i*ny,0)),int(round(j/m.r.last()*nr,0)),int(round(k*nt,0))]
      # return Cc.loc[i].loc[j].loc[k].val

    def _Tc_ic_(m,i,j,k):
        return init['Tc'][int(round(i*ny,0)),int(round(j/m.r.last()*nr,0)),int(round(k*nt,0))]
        #return Tc.loc[i].loc[j].loc[k].val

    def _X_ic_(m,i,j,k):
       return init['X'][int(round(i*ny,0)),int(round(j/m.r.last()*nr,0)),int(round(k*nt,0))]
      # return X.loc[i].loc[j].loc[k].val
    
    def _dCcdr2_init(m,i,j,k):
        return -50


    # Define model (for optimization)
    m = pyo.ConcreteModel()

    ###############################################################################
    # Define parameters
    ###############################################################################

    # Note: reorganizing such that parameters are retrieved from an external source, to assist with consistency if I need to make changes
    # Reactor Design
    m.tf = pyo.Var(initialize=1400,bounds=(900,3600))#pyo.Param(initialize=t_solve)     #     pyo.Var(initialize=t_solve,bounds=(500,3000))#          # interval duration, s
    m.L = pyo.Var(initialize=3.1, bounds=(2.38,7.38))#pyo.Param(initialize=L)#1.75)#dp.L)            # reactor length, m
    m.Di = pyo.Var(initialize=1.5,bounds=(0.99,2.47))#pyo.Param(initialize=D)                       # reactor diameter, m
    m.Dp = pyo.Param(initialize=dp.Dp)                      # particle diameter, m
    m.Ac = pyo.Expression(expr=np.pi*m.Di**2/4.0)          # reactor cross-sectional area, m²

    # Reactor Conditions
    m.P = pyo.Param(initialize=op.P)                        # inlet reactor pressure, bar
    m.T_air = pyo.Param(initialize=op.T_air)                # inlet temperature, K
    m.Tref = pyo.Param(initialize=op.T_air)                     # reference temperature for certain calculations
    m.T_in = pyo.Param(initialize=op.T_air)

    # Flexible Parameters
    m.rho = pyo.Param(initialize=op.rho)                    # density of air at 700 K and 1 bar, kg/m³ -> CHANGE if pressure is different
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
    m.mu = pyo.Param(initialize=op.mu)                      # dynamic viscosity of air at 700 K, kg/m/s
    m.dij = pyo.Param(initialize=op.dij)                    # binary diffusivity coefficient of oxygen in nitrogen at 700 K, m²/s
    m.tau = pyo.Param(initialize=dp.tau)                    # particle tortuosity
    m.rho_s = pyo.Param(initialize=op.rho_s)                # particle density (pure solid, not bed), Ni on Al2O3, kg/m³

    # Heat transfer parameters
    m.cpf = pyo.Param(initialize=op.cpf)                    # heat capacity of air, J/kg/K
    m.cp_f = pyo.Param(initialize=m.cpf*m.mw_air/1000.0)    # heat capacity of air, J/mol/K
    m.l_e0 = pyo.Param(initialize=op.l_e0)                  # static contribution effect of thermal conductivity, J/m/s/K
    m.cp_s = pyo.Param(initialize=op.cp_s)                  # oxygen carrier specific heat capacity, J/kg/K
    m.dH = pyo.Param(initialize=op.dH)                      # heat of reaction for oxidation, J/mol of O2
    m.l_i = pyo.Param(initialize=op.l_i)                    # thermal conductivity of Ni, W/m/K
    trace = 1e-4
    # Reaction Parameters
    m.a0 = pyo.Param(initialize=op.a0)                      # initial specific surface area of OC, m²/kg

    # Define radial, temporal, and axial domains
    m.r = pdae.ContinuousSet(bounds=(0.0,m.Dp/2))
   # m.r = pdae.ContinuousSet(initialize=[0,0.00125,0.0017,0.002,0.00225,0.0025])
    #m.t = pdae.ContinuousSet(bounds=(0.0,1.0))
    #Note I increased no. of discrete points for oxidation local SA to improve resolution of graphs
  #  m.t = pdae.ContinuousSet(initialize=[0,0.002,0.005,0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.125,0.15,0.175,0.2,0.225,0.25,0.275,0.3,0.325,0.35,0.375,0.4,0.425,0.45,0.475,0.5,0.525,0.55,0.6,0.625,0.65,0.675,0.7,0.725,0.75,0.775,0.8,0.825,0.85,0.875,0.9,0.95,1])
   # m.t = pdae.ContinuousSet(initialize=[0,0.002,0.005,0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.125,0.15,0.175,0.2,0.225,0.25,0.275,0.3,0.325,0.35,0.375,0.4,0.425,0.45,0.475,0.5,0.6,0.65,0.7,0.75,0.8,0.9,1])
   # m.t = pdae.ContinuousSet(bounds=(0,1),initialize=[0,0.002,0.005,0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.125,0.15,0.175,0.2,0.225,0.25,0.275,0.3,0.325,0.35,0.375,0.4,0.425,0.45,0.475,0.5,0.6,0.7,0.8,0.9,1])
    #               base    m.t = pdae.ContinuousSet(bounds=(0,1),initialize=[0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.9,1])
   #m.y = pdae.ContinuousSet(bounds=(0.0,1.0))
  #            base       m.y = pdae.ContinuousSet(bounds=(0,1),initialize=[0,0.01,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1])
   # m.y = pdae.ContinuousSet(initialize=[0,0.01,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1])
    m.y = pdae.ContinuousSet(initialize=[0,0.02,0.04,0.06,0.08,0.1,0.12,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1])
    m.t = pdae.ContinuousSet(initialize=[0,0.002,0.005,0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.125,0.15,0.175,0.2,0.225,0.25,0.275,0.3,0.325,0.35,0.375,0.4,0.425,0.45,0.475,0.5,0.525,0.55,0.6,0.65,0.7,0.75,0.8,0.9,1])
    m.hr = pyo.Param(initialize=m.r.last()/nr)
    m.hy = pyo.Param(initialize=m.y.last()/ny)

    # Define state variables
    m.C = pyo.Var(m.y, m.t, initialize=_C_ic_,domain=pyo.NonNegativeReals)              # Concentration in fluid phase, mol/m³
    m.T = pyo.Var(m.y, m.t, initialize=_T_ic_,domain=pyo.NonNegativeReals)              # Temperature in fluid phase, K
    m.Cc = pyo.Var(m.y, m.r, m.t, initialize=_Cc_ic_,domain=pyo.NonNegativeReals)       # Concentration in the particle, mol/m³
    m.Tc = pyo.Var(m.y, m.r, m.t, initialize=_Tc_ic_,domain=pyo.NonNegativeReals)       # Temperature in the particle, K
    m.X = pyo.Var(m.y, m.r, m.t, initialize=_X_ic_,bounds=(0.01,0.99))#(1e-8,0.999))         # Solid conversion in the particle

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
    m.dCcdr2 = pdae.DerivativeVar(m.Cc, wrt=(m.r,m.r))#,initialize=_dCcdr2_init) #pyo.Expression(m.y,m.r,m.t,expr=0)# pdae.DerivativeVar(m.Cc, wrt=(m.r,m.r),initialize=_dCcdr2_init)
    m.dTcdr = pdae.DerivativeVar(m.Tc, wrt=m.r)
    m.dTcdr2 = pdae.DerivativeVar(m.Tc, wrt=(m.r,m.r))

    # Differentiate system
    disc_fd = pyo.TransformationFactory('dae.finite_difference')
    disc_oc = pyo.TransformationFactory('dae.collocation')
    #disc_fd.apply_to(m, wrt=m.t, nfe=nt, scheme='BACKWARD')
    disc_fd.apply_to(m, wrt=m.t, nfe=nt, scheme='BACKWARD')
    disc_fd.apply_to(m, wrt=m.r, nfe=nr, scheme='CENTRAL')
    disc_oc.apply_to(m, wrt=m.y, nfe=nint, ncp=mcp, scheme='LAGRANGE-RADAU')
 #   disc_oc.apply_to(m, wrt=m.r, nfe=5, ncp=1, scheme='LAGRANGE-RADAU')


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
    m.m_air = pyo.Var(initialize=0.401,bounds=(0.2*0.401,10*0.401)) # kg/s
    m.Gair = pyo.Expression(expr=m.m_air/m.Ac)#(initialize=op.G,bounds=(0.5*op.G,4*op.G))#pyo.Var(initialize=op.G,bounds=(0.06,10))
    m.Ginert = pyo.Param(initialize=0)
  
    m.Dmi = pyo.Param(initialize=m.dij)                          # diffusion coefficient of oxygen in the mixture, m²/s
    m.Sc = pyo.Param(initialize=m.mu/(m.rho*m.Dmi))              # Schmidt number of oxygen
    m.av = pyo.Param(initialize=6*(1-m.eb)/m.Dp)                 # external particle surface area per unit volume, 1/m

    # Total inlet mass flux; kg/m2/s
    def _G_(m):
        return m.Gair + m.Ginert
    m.G = pyo.Expression(rule=_G_)

    # # Gas velocity, m/s
    # def _yi_(m):
    #     return op.yi*m.Gair/m.G
    m.yi = pyo.Var(initialize=op.yi,bounds=(0.08,0.21))
   # m.yi = pyo.Var(initialize=op.yi,bounds=(0.01,0.21))

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
    m.k_ox = pyo.Param(initialize=2.706e-7)
    m.order = pyo.Param(initialize=1)#pyo.Var(initialize=1,bounds=(0.1,4))

    # Reaction rate
    

    def _rxn_(m, i, j, k):
        # if i < m.Lq/m.L or i > (m.Lq+m.LRM)/m.L: # no reaction occurs
        #     return 0
        # #return m.a0*m.k_ox*m.Cc[i,j,k]**m.order*m.cFe2O30_1*2*(1-m.X[i,j,k])*(m.log_term[i,j,k])**0.5#(1-m.X[i,j,k])
        # else:
        return m.a0*m.k_ox*m.Cc[i,j,k]**m.order*m.cFe2O30_1*2*(0.99-m.X[i,j,k])*(-pyo.log(1.000000000000000001-m.X[i,j,k]))**0.5#(1-m.X[i,j,k])
           # return m.a0*m.k_ox*m.Cc[i,j,k]**m.order*m.cFe2O30_1*2*m.poly_log[i,j,k]#(1-m.X[i,j,k])
    m.r_ox = pyo.Expression(m.y, m.r, m.t, rule=_rxn_)

    # Previous reduction stage conversion calculation

    # zs_int = np.array([ 0, 0.003333, 0.01, 0.023333, 0.05, 0.066667, 0.1, 0.116667, 0.15,
    #                    0.166667, 0.2, 0.216667, 0.25, 0.266667, 0.3, 0.316667, 0.35,
    #                    0.358333, 0.375, 0.383333, 0.4, 0.416667, 0.45, 0.466667, 0.5, 
    #                    0.516667, 0.55, 0.566667, 0.6, 0.616667, 0.65, 0.666667, 0.7, 0.716667,
    #                    0.75, 0.766667, 0.8, 0.816667, 0.85, 0.866667, 0.9, 0.916667, 0.95,
    #                    0.966667, 1])## np.full(len(m.y),0)
    zs_int = np.array([0, 0.006667, 0.02, 0.026667, 0.04, 0.046667, 0.06, 0.066667, 0.08, 0.086667, 0.1, 0.106667, 0.12, 0.13, 0.15, 0.166667, 0.2, 0.216667, 0.25, 0.266667, 0.3, 0.316667, 0.35, 0.358333, 0.375, 0.383333, 0.4, 0.416667, 0.45, 0.466667, 0.5, 0.516667, 0.55, 0.566667, 0.6, 0.616667, 0.65, 0.666667, 0.7, 0.716667, 0.75, 0.766667, 0.8, 0.816667, 0.85, 0.866667, 0.9, 0.916667, 0.95, 0.966667, 1])

    r_int = np.array([0,0.000833,0.001667,0.0025])#np.full(len(m.r),0)
    # for i in range(len(zs_int)):
    #     zs_int[i] = m.y.at(i+1)
    # for j in range(len(r_int)):
    #     r_int[j] = m.r.at(j+1)
    X_red = _X_avg_(filename_X,zs_int,r_int,dp.Dp/2,dp.L)

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
    #return m.eb*m.Dax[k]*(1/m.L)*m.dCdy[n,0,k] == (m.Q[k]*m.C[n,0,k]-m.yis[n]*m.Q[k]*m.P/(m.Rgb*m.T_in)/nd[0])/m.Ac
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
      return (m.T_air/nd[1] - m.T[m.y.at(1),k])**2 <= 1e-6# m.eb*m.l_ax[k]*(1/m.L)*m.dTdy[0,k] == (m.T[0,k]-m.T_air/nd[1])*m.cp_f*m.Q[k]*m.P/(m.Rgb*m.T_air*m.Ac)
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
    
    def _Ergun_(m,k):
        return m.L*(150*(1-m.eb)**2/(m.eb**3)*m.mu*m.vel[k]/(m.Dp**2)+1.75*(1-m.eb)/(m.eb**3)*m.rho*m.vel[k]**2/(m.Dp))
     
    m.dP = pyo.Expression(m.t,rule=_Ergun_)
    m.umf = pyo.Var(initialize=0.5,domain=pyo.NonNegativeReals)
    def _umf_(m):
        return 1305.556*m.rho*m.umf*m.Dp/m.mu + 30.789*(m.rho*m.umf*m.Dp/m.mu)**1.657 + 19.474*(m.rho*m.umf*m.Dp/m.mu)**2-m.rho*m.Dp**3*(m.rho_s-m.rho)*9.81/(m.mu**2) == 0
    m.u_mf = pyo.Constraint(rule=_umf_)
    m.LD = pyo.Expression(expr=m.L/m.Di)
    def OC_conversion(m,i):
        return m.X[i,m.r.last(),m.t.last()]
    m.final_OC_conversion = pdae.integral.Integral(m.y,wrt=m.y,rule=OC_conversion)
    def Q_out(m): # Average heat extracted via outlet gas
        value = 0
        for t in range(len(m.t)):
            if t == 0:
                value += 0
            else:
                value += (m.t.at(t+1)-m.t.at(t))*m.Q[m.t.at(t)]*(m.C[m.y.last(),m.t.at(t)]*nd[0]+(1-m.yi)*op.P/gp.Rgb/op.T_air)*m.cp_f*(m.T[m.y.last(),m.t.at(t)]*nd[1]-op.T_air)
        return value
    def Q_in(m): # Average heat remaining in the packed bed at the end of oxidation stage
        value = 0
        for z in range(len(m.y)):
            if z == 0:
                value += 0
            else:
                value += m.L*(m.y.at(z+1)-m.y.at(z))*m.rho_s*(1-m.eb)*(1-m.ec)*m.Ac*m.cp_s*(m.Tc[m.y.at(z),m.r.last(),m.t.last()]*nd[1]-op.T_air)
        return value/m.tf
    m.Q_in = pyo.Expression(rule=Q_in)
    m.Q_out = pyo.Expression(rule=Q_out)
    
    # integrated heats
    def Q_out_int(m,k):
       return  m.Q[k]*(m.C[m.y.last(),k]*nd[0]+(1-m.yi)*op.P/gp.Rgb/op.T_air)*m.cp_f*(m.T[m.y.last(),k]*nd[1]-op.T_air)
   
    def Q_in_int(m,i):
        return m.L*m.rho_s*(1-m.eb)*(1-m.ec)*m.Ac*m.cp_s*(m.Tc[i,m.r.last(),m.t.last()]*nd[1]-op.T_air)
    
    m.Q_out_int = pdae.integral.Integral(m.t,wrt=m.t,rule=Q_out_int)
    m.Q_in_int = pdae.integral.Integral(m.y,wrt=m.y,rule=Q_in_int)
  
    ###########################################################################
    # Design Constraints
    ###########################################################################
    
    # Velocity less than fluidization
    def _velocity_constraint_(m,k):
        return m.vel[k] <= m.umf
    m.velocity_constraint = pyo.Constraint(m.t,rule=_velocity_constraint_)
    
    # Reactor length-to-diameter ratio
    def _LD_H_(m):
        return m.L/m.Di <= 6
    m.LDH = pyo.Constraint(rule=_LD_H_)
    
    def _LD_L_(m):
        return m.L/m.Di >= 2
    m.LDL = pyo.Constraint(rule=_LD_L_)
    
    def _max_dP_(m,k):
        return m.dP[k] <= 100000 # Pa
    m.maxdP = pyo.Constraint(m.t,rule=_max_dP_)
    
    def _max_T(m,i,k): # maximum temperture constraint
        return m.T[i,k] <= (1100+273)/nd[1]
    m.max_T = pyo.Constraint(m.y,m.t,rule=_max_T)
    
    # Ensure near full OC conversion
    
    def X_final_last(m,i):
        return m.X[i,m.r.last(),m.t.last()] >= 0.95
    m.OC_conversion = pyo.Constraint(m.y,rule=X_final_last)
   
    ###########################################################################
    # Economics
    ###########################################################################
    # Number of cycles in one year
    m.t_op = pyo.Param(initialize=8000) # number of operating hours in one year
    m.t_purge = pyo.Param(initialize=0.25) # total length of purge stages in one cycle, h
    m.t_ox = pyo.Expression(expr = m.tf/3600) # length of the oxidation stage, h
    m.t_red = pyo.Param(initialize=t_red/3600) # length of the reduction stage, h
    m.n_cycles = pyo.Expression(expr = m.t_op/(m.t_purge+m.t_ox+m.t_red)) # number of cycles in one year
    
    # Cost Parameters
    m.Pe = pyo.Param(initialize=0.157/3.6e6)                    # cost of electricity, $/J
    m.eta_fan = pyo.Param(initialize=0.75)                      # fan efficiency
    m.eta_elec = pyo.Param(initialize=0.4)                      # electrical efficiency factor
    m.eta_heat = pyo.Param(initialize=0.7)                      # heat efficiency factor (assumed)
    m.i = pyo.Param(initialize=0.0795)                          # interest rate for annualization factor
    m.n = pyo.Param(initialize=20)                              # reactor lifetime, years
    m.CEPCI25 = pyo.Param(initialize=796.5)                     # CEPCI 2025
    m.CEPCI94 = pyo.Param(initialize=368.1)                     # CEPCI 1994
    m.MS = pyo.Param(initialize=1050)                           # M&S index 1994
    m.CAD_conversion = pyo.Param(initialize=1.37)               # conversion factor for $USD to $CAD
    m.Fp = pyo.Param(initialize=1)                              # pressure factor (1 is for up to 50 psi)
    m.Fm = pyo.Param(initialize=3.67)                           # material factor (1 is for carbon steel, 3.67 is for stainless steel)
    m.Fc = pyo.Param(initialize=m.Fp*m.Fm)                      # combined material/pressure factor    
    m.cRM = pyo.Param(initialize=520)                           # cost of red mud OC, $USD(2014)/tn
    m.USD_yr_conv = pyo.Param(initialize=1.39)                  # conversion factor for $USD(2014) to $USD(2025)
    m.tn_kg_conv = pyo.Param(initialize=907.185)                # conversion factor, kg/US tn
    m.OC_replacement = pyo.Param(initialize=75)                 # number of cycles for stable OC operation
   # m.n_OC_changes = pyo.Param(initialize=39)                   # number of OC changes in 1 year, based on Hong et al. attrition rate
    m.n_OC_changes = pyo.Param(initialize=1)                   # number of OC changes in 1 year, assumed

    # Flow rate cost
    m.Fan_cost = pyo.Expression(expr=m.Pe*m.Q[m.t.last()]*m.dP[m.t.last()]/m.eta_fan*m.tf*m.n_cycles) # cost to overcome pressure drop, $/year
    
    # Heat profit
    m.Q_total = pyo.Expression(expr=(m.Q_in+m.Q_out)*m.tf) # total energy produced, J
    m.heat_profit = pyo.Expression(expr=m.Q_total*m.Pe*m.eta_elec*m.eta_heat*m.n_cycles) # profit from heat, $/year
    
    # Reactor capital cost
    m.AF = pyo.Expression(expr=(m.i*(1+m.i)**m.n)/((1+m.i)**m.n-1)) # annualization factor
    m.D_ft = pyo.Expression(expr=m.Di*3.28) # reactor diameter, ft
    m.L_ft = pyo.Expression(expr=m.L*3.28) # reactor length, ft
    m.capital = pyo.Expression(expr=(m.CEPCI25/m.CEPCI94)*(m.MS/280)*101.9*m.D_ft**1.066*m.L_ft**0.802*(2.18+m.Fc)*m.CAD_conversion) # total capital cost, $
    m.capital_annualized = pyo.Expression(expr=m.capital*m.AF)  # annualized capital cost, $/year
    
    # OC replacement cost
    m.cRM_OC = pyo.Expression(expr=m.cRM*m.USD_yr_conv/m.tn_kg_conv) # cost of OC, $CAD/kg
    # m.OC_cost = pyo.Expression(expr=m.Ac*m.L*(1-m.ec)*(1-m.eb)*m.rho_s*m.cRM_OC*m.n_cycles/m.OC_replacement) # cost for OC replacement, $/year based on assumed number of stable cycles
    m.OC_cost = pyo.Expression(expr=m.Ac*m.L*(1-m.ec)*(1-m.eb)*m.rho_s*m.cRM_OC*m.n_OC_changes*m.CAD_conversion) # cost for OC replacement, $/year based on attrition rate

    # total annualized cost
    m.obj = pyo.Objective(expr=-m.Fan_cost+m.heat_profit-m.capital_annualized-m.OC_cost, sense=pyo.maximize)
    
    return m

if __name__ == "__main__":
     # Import Yin's data for model fitting/validation
    wdir = os.path.dirname(os.path.realpath(__file__))


    tfreq = 20#50#200#50#75#50#22.5                    # Frequency of time points in discretization, s
    
    
    tf = 1400#500#1000#2000#3600#7000#8200                      # Reduction stage simulation time, s
    ncon = 1#int(t_pred/t_samp)       # Number of flux steps (control actions) per interval
    
    
    cFe2O3 = 0.4928                                                      # Fe2O3 concentration (mass fraction)
    

    
    # Define nondimensionalizing coefficients for concentration and temperature
    Cnd = op.yi*op.P/(gp.Rgb*op.T_air)
    Tnd = op.T_air
    nd = [Cnd,Tnd]
     
    # Discretization
 
    nr = 3#6                          # Number of radial nodes (particle)
    nint = 25#22#22#14#20#10                        # Number of axial elements (bulk reactor)
    mcp = 2#3#5#4                         # Total number of orthogonal collocation points sent to pyomo
    nt = 38#21# 34# 38#46#38
    # nt = 46#38# int(tf/tfreq)          # Number of temporal nodes
    ntcon = int(nt/ncon)            # Number of nodes in each step of the control profile
    ny = nint*mcp                   # Total number of axial points
    l_t = tf/nt                 # Length of each time node which must be solved by solve_ivp
    disc = [nt, nint, mcp, nr]

    trace = 1e-4
    yi = op.yi# [1-5*trace, trace, trace, trace, trace, trace] #inlet mole fractions
    
    # Initial Conditions
    C0 = trace*op.P/gp.Rgb/op.T_air    # initial fluid-phase concentrations of compounds, mol/m³
    T0 = op.T_air                                          # initial temperature at reactor inlet, K
    X0 = 0.5#0.01                                           # initial conversion of NiO to Ni
    
    # Define initial guesses (and initial conditions for t=0)
    init = {'C': np.full([ny+1,nt+1],0,dtype=float),
            'T': np.full([ny+1,nt+1],0,dtype=float),
            'Xp': np.full([ny+1,nt+1],0,dtype=float),
            'Cc': np.full([ny+1,nr+1,nt+1],0,dtype=float),
            'Tc': np.full([ny+1,nr+1,nt+1],0,dtype=float),
            'X': np.full([ny+1,nr+1,nt+1],0,dtype=float)}
    
    # Define matrices to collect data
    plantsolves = []
    ts = np.full([nt+1],0,dtype=float)#[k*l_t for k in range(int(tf/l_t)+1)]
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
    
    optfile='ipopt1.opt'
    #plantopt = pyo.SolverFactory('ipopt', executable='C:\Dana\ipopt')'
    plantopt = pyo.SolverFactory('ipopt')#, executable='C:\Jaime\MSYS2\home\ja2mccor\ipopt')#pyo.SolverFactory('ipopt', executable='C:\software\cygwin\home\ipopt')#'C:\Dana\ipopt')
    plantopt.options['halt_on_ampl_error'] = 'yes'
    #plantopt.options["option_file_name"] = optfile
    plantopt.options['linear_solver'] = 'ma97'
    # plantopt.options['acceptable_tol'] = 1e-5
    # plantopt.options['tol'] = 1e-5
    #plantopt.options['ma86_u'] = 0.5
   # plantopt.options['ma57_pivtol'] = 0.0001    
    plantopt.options['ma97_u'] = 0.5
    #plantopt.options['ma77_u'] = 0.00001
  #  plantopt.options['nlp_scaling_method'] = 'gradient-based'
    # plantopt.options['bound_push'] = 1e-6
    # plantopt.options['bound_frac'] = 1e-6    
    # plantopt.options['mu_strategy'] = 'adaptive'
    # plantopt.options['mu_oracle'] = 'quality-function'
    # plantopt.options['expect_infeasible_problem'] = 'yes'
    # plantopt.options['print_level'] = 5
    # plantopt.options['hessian_approximation'] = 'limited-memory'
    # plantopt.options['derivative_test'] = 'second-order'
    # plantopt.options['derivative_test_print_all'] = 'yes'
    
    
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
    filename_X = 'X_red_baseline_feas_newdisc.xlsx'
    Xinit  = pd.read_excel(filename_X, index_col=[0,1,2]) # 0: y, 1: r, 2: t
    # calculate the average conversion of the OC at the end of reduction
    
   
    X_red_old = 0.5 #placeholder
  #  X_red = 0.5 #placeholder
    # Define plant model
    plant = oxidation_pyo(init, disc, nd, cFe2O3,'X_red_baseline_feas_newdisc.xlsx',1300)

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
  
  
    # zs_int = np.array([ 0, 0.003333, 0.01, 0.023333, 0.05, 0.066667, 0.1, 0.116667, 0.15,
    #                    0.166667, 0.2, 0.216667, 0.25, 0.266667, 0.3, 0.316667, 0.35,
    #                    0.358333, 0.375, 0.383333, 0.4, 0.416667, 0.45, 0.466667, 0.5, 
    #                    0.516667, 0.55, 0.566667, 0.6, 0.616667, 0.65, 0.666667, 0.7, 0.716667,
    #                    0.75, 0.766667, 0.8, 0.816667, 0.85, 0.866667, 0.9, 0.916667, 0.95,
    #                    0.966667, 1])## np.full(len(m.y),0)
    zs_int = np.array([0, 0.006667, 0.02, 0.026667, 0.04, 0.046667, 0.06, 0.066667, 0.08, 0.086667, 0.1, 0.106667, 0.12, 0.13, 0.15, 0.166667, 0.2, 0.216667, 0.25, 0.266667, 0.3, 0.316667, 0.35, 0.358333, 0.375, 0.383333, 0.4, 0.416667, 0.45, 0.466667, 0.5, 0.516667, 0.55, 0.566667, 0.6, 0.616667, 0.65, 0.666667, 0.7, 0.716667, 0.75, 0.766667, 0.8, 0.816667, 0.85, 0.866667, 0.9, 0.916667, 0.95, 0.966667, 1])

    r_int = np.array([0,0.000833,0.001667,0.0025])#np.full(len(m.r),0)

    X_red = _X_avg_(filename_X,zs_int,r_int,dp.Dp/2,dp.L)
    # print(X_red_old)
    print(X_red)
    start = time.time()
    # Solve plant model


    
    # plant.tf.fix(1575)
    # plant.Gair.fix(0.06)
    # plant.yi.fix(0.01)
#    
   # plant.k_ox = 2.706e-7#6e-7#6.7e-7#1.07e-7#1.842e-6#6.7e-7#1.07e-7#6e-7#1.842e-06#1.07e-7
    plant.OC_conversion.deactivate()
    plant.m_air.fix(0.401)
    plant.yi.fix(0.18)
    plant.tf.fix(1400)
    plant.L.fix(3.1)
    plant.Di.fix(1.5)
    plant.obj.deactivate()
    plantresults = plantopt.solve(plant,tee=True)
    plantsolves += [plantresults['Solver']]
    print('End of run 1')
    
    plant.m_air.unfix()
    plant.yi.unfix()
    plant.tf.unfix()
    #plant.L.unfix()
    #plant.Di.unfix()
    plant.obj.activate()
    # #plant.OC_conversion.activate()
    plantresults = plantopt.solve(plant,tee=True)
    plantsolves += [plantresults['Solver']]
    plant.OC_conversion.activate()
    plant.L.unfix()
    plant.Di.unfix()
    plantresults = plantopt.solve(plant,tee=True)
    plantsolves += [plantresults['Solver']]
##     #Compile data

    for k in range(nt+1):
      # Gs = plant.G.value
    #   SCO2s[k] =pyo.value(plant.SCO2[plant.t.at(k+1)])
       ts[k] = plant.t.at(k+1)*pyo.value(plant.tf)
       youts[k] = plant.yout[plant.t.at(k+1)]()

       #youtsdry[k] = plant.youtdry[plant.t.at(k+1)]()
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
    foldername = "Oxidation1"#'O2 Concentration''-'+str(op.yi)#str(round(pyo.value(plant.Gair),3))#+str(tf)+"-"+str(round(pyo.value(plant.Gair),3))+'final_design'
    
    # Calculating average converison after oxidation
    filename_X_ox = 'X_oxidation2.xlsx'
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
    X_ox = _X_avg_(filename_X_ox,zs_int,r_int,dp.Dp/2,pyo.value(plant.L))
    
    data = {'value':[pyo.value(plant.vel[plant.t.last()]),pyo.value(plant.umf*plant.Dp*plant.rho/plant.mu),round(pyo.value(9.81*plant.Dp**3*plant.rho*(plant.rho_s-plant.rho)/(plant.mu**2))),pyo.value(plant.LD),tf/60,plant.Dp/(plant.Di/2),pyo.value(plant.dP[plant.t.last()])/101325*100,pyo.value(plant.Di/plant.Dp)],'restriction':[pyo.value(plant.umf),'<350','<5000000','2-6','15 min','0.04-0.25','<8%','>2']}
    print(pd.DataFrame(data,index=[ 'Gas velocity','Remf','Ar','L/D','Cycle time','dp/R0','Pressure drop %','D/dp']))
    print('min temp: ',np.min(Ts))  
    print('max temp: ',np.max(Ts)) 
    print('min conc: ',np.min(Cs))  
    print('max conc: ',np.max(Cs)) 
    print('average conversion increase: ',np.average(Xs[:,-1,-1]-Xs[:,-1,0]))
    print('average final conversion: ',X_ox)#np.average(Xs[5:-1,-1,-1]))
    print('total heat: ',pyo.value(plant.Q_out)/1000+pyo.value(plant.Q_in)/1000,' kW')
    save = True

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
             f.write('\nExtracted heat: '+str(pyo.value(plant.Q_out)/1000))
             f.write('\nAccumulated heat: '+str(pyo.value(plant.Q_in)/1000))
             f.write('\ntotal heat: '+str(pyo.value(plant.Q_out)/1000+pyo.value(plant.Q_in)/1000)+' kW')
             f.write('\n cFe2O3: '+str(cFe2O3))
             f.write('\nDesign Variables:')
             f.write('\n  tf: '+str(pyo.value(plant.tf)))
             f.write('\n  yO2: '+str(pyo.value(plant.yi)))
             f.write('\n  G: '+str(pyo.value(plant.G)))
             f.write('\n  Di: '+str(pyo.value(plant.Di)))
             f.write('\n  L: '+str(pyo.value(plant.L)))
             f.write('\nObjective Values:')
             f.write('\n  Flow rate: '+str(pyo.value(plant.Fan_cost)))
             f.write('\n  Heat: '+str(pyo.value(plant.heat_profit)))
             f.write('\n  Reactor CAPEX: '+str(pyo.value(plant.capital_annualized)))
             f.write('\n  OC Cost: '+str(pyo.value(plant.OC_cost)))
             f.write('\n  TAC: '+str(pyo.value(plant.obj)))
             
    # filename_X_ox = 'X_ox.xlsx'
    # ii = 0
    # writer = pd.ExcelWriter(filename_X_ox)
    # df = pd.DataFrame(columns=['i','j','k','val'])
    # for i in plant.X.keys():
    #     val = pyo.value(plant.X[i])
    #     #print("X[{}] = {}".format(i,val))
    #     df.loc[ii] = [i[0],i[1],i[2],val]
    #     ii += 1
    # df.to_excel(writer,sheet_name='X',index=None)
    # writer.close()
    # X_ox = _X_avg_(filename_X_ox,zs_int,r_int,dp.Dp/2,dp.L)
    
    
    
    
             
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
    
    # Variables
    # total_vars = len(list(plant.component_data_objects(pyo.Var)))
    # active_vars = len([v for v in plant.component_data_objects(pyo.Var) if not v.fixed])
    # fixed_vars = len([v for v in plant.component_data_objects(pyo.Var) if v.fixed])
    
    # # Constraints
    # total_cons = len(list(plant.component_data_objects(pyo.Constraint)))
    # active_cons = len(list(plant.component_data_objects(pyo.Constraint, active=True)))
    # eq_cons = len([c for c in plant.component_data_objects(pyo.Constraint, active=True)
    #                if c.equality])
    
    # print("Total variables:", total_vars)
    # print("Free variables:", active_vars)
    # print("Fixed variables:", fixed_vars)
    # print("Total constraints:", total_cons)
    # print("Active constraints:", active_cons)
    # print("Equality constraints:", eq_cons)
    
    result_vals = {'value':[pyo.value(plant.C[plant.y.at(4),plant.t.at(21)]),pyo.value(plant.T[plant.y.at(4),plant.t.at(21)]),pyo.value(plant.Cc[plant.y.at(4),plant.r.last(),plant.t.at(21)]),pyo.value(plant.Tc[plant.y.at(4),plant.r.last(),plant.t.at(21)]),pyo.value(plant.X[plant.y.at(4),plant.r.last(),plant.t.at(21)]),\
                            pyo.value(plant.dCdt[plant.y.at(4),plant.t.at(21)]),pyo.value(plant.dTdt[plant.y.at(4),plant.t.at(21)]),pyo.value(plant.dCcdt[plant.y.at(4),plant.r.last(),plant.t.at(21)]),pyo.value(plant.dTcdt[plant.y.at(4),plant.r.last(),plant.t.at(21)]),pyo.value(plant.dXdt[plant.y.at(4),plant.r.last(),plant.t.at(21)]),\
                            pyo.value(plant.dCdy[plant.y.at(4),plant.t.at(21)]),pyo.value(plant.dTdy[plant.y.at(4),plant.t.at(21)]),pyo.value(plant.dCcdr[plant.y.at(4),plant.r.last(),plant.t.at(21)]),pyo.value(plant.dTcdr[plant.y.at(4),plant.r.last(),plant.t.at(21)]),\
                            pyo.value(plant.dCdy2[plant.y.at(4),plant.t.at(21)]),pyo.value(plant.dTdy2[plant.y.at(4),plant.t.at(21)]),pyo.value(plant.dCcdr2[plant.y.at(4),plant.r.last(),plant.t.at(21)]),pyo.value(plant.dTcdr2[plant.y.at(4),plant.r.last(),plant.t.at(21)]),pyo.value(plant.r_ox[plant.y.at(4),plant.r.last(),plant.t.at(21)])]\
                   }
    data_frame = pd.DataFrame(result_vals,index=['C','T','Cc','Tc','X','dCdt','dTdt','dCcdt','dTcdt','dXdt','dCdy','dTdy','dCcdr','dTcdr','dCdy2','dTdy2','dCcdr2','dTcdr2','r_ox'])
    data_frame.to_csv('out_nr3.csv',index=False)
'''
    # Plot results
    c = 0
    ts = ts
    legz0 = " z=0"
    legz1 = " z=1"
    legz2 = " z=L/2"
    legzL = " z=L"
    
    print("Analytical ln")
    c += 1
    fig = plt.figure(c)
    plt.plot(ts,np.transpose(Xs[11,-1,:]),'b-', ts,np.transpose(Xs[int(ny//2),-1,:]),'g-', ts,np.transpose(Xs[29,-1,:]),'r-')
    plt.xlabel("Time (s)")
    plt.title("analytical kinetics")
    plt.ylabel("Conversion")
    plt.legend(["ms"+legz0, "ms"+legz2, "ms"+legzL])
    
    c += 1
    fig = plt.figure(c)
    plt.plot(ts,youts[:],label ="Multiscale simulation")
    plt.plot(tO2_Liu, C_O2_Liu, 'o', label="Liu et al")
    plt.title("analytical kinetics")
    plt.xlabel("Time (s)")
    plt.ylabel("Outlet fraction ")
    plt.legend()
    
 
    
    c += 1
    fig = plt.figure(c)
    for i in range(11,30):
        plt.plot(ts,Xs[i,-1,:],label=plant.y.at(i))
    plt.legend(loc='upper center', bbox_to_anchor=(0.5, 1.5),fancybox=True, shadow=True, ncol=4)
    plt.title('conversion at different points in RM')
    
    c += 1
    fig = plt.figure(c)
    for k in range(0,25):
        plt.plot(zs,Cs[:,k],label= plant.t.at(1+k)*tf)
    plt.legend(loc='upper center', bbox_to_anchor=(0.5, 1.5),fancybox=True, shadow=True, ncol=4)
    plt.show()
    #print('optimized k times ', factor)
    print(time.time()-start)
    c += 1
    fig = plt.figure(c)
    plt.plot(zs,r_CO[:,6])
    plt.xlabel('packed bed length (m)')
    plt.ylabel('reaction rate')
    plt.title('A2 200 s')
    plt.show()
#    #print(plant.obj.value())
#    #log_infeasible_constraints(plant)

'''
'''
filename = 'C.xlsx'
ii = 0
writer = pd.ExcelWriter(filename,engine='xlsxwriter')
df = pd.DataFrame(columns=['n','i','k','val'])
for i in plant.C.keys():
    val = pyo.value(plant.C[i])
   # print("C[{}] = {}".format(i,val))
    df.loc[ii] = [i[0],i[1],i[2],val]
    ii += 1
df.to_excel(writer,sheet_name='C',index=None)
writer.close()

filename = 'Cc.xlsx'
ii = 0
writer = pd.ExcelWriter(filename,engine='xlsxwriter')
df = pd.DataFrame(columns=['n','i','j','k','val'])
for i in plant.Cc.keys():
    val = pyo.value(plant.Cc[i])
   # print("Cc[{}] = {}".format(i,val))
    df.loc[ii] = [i[0],i[1],i[2],i[3],val]
    ii += 1
df.to_excel(writer,sheet_name='Cc',index=None)
writer.close()

filename = 'T.xlsx'
ii = 0
writer = pd.ExcelWriter(filename,engine='xlsxwriter')
df = pd.DataFrame(columns=['i','k','val'])
for i in plant.T.keys():
    val = pyo.value(plant.T[i])
    #print("T[{}] = {}".format(i,val))
    df.loc[ii] = [i[0],i[1],val]
    ii += 1
df.to_excel(writer,sheet_name='T',index=None)
writer.close()

filename = 'Tc.xlsx'
ii = 0
writer = pd.ExcelWriter(filename,engine='xlsxwriter')
df = pd.DataFrame(columns=['i','j','k','val'])
for i in plant.Tc.keys():
    val = pyo.value(plant.Tc[i])
   # print("Tc[{}] = {}".format(i,val))
    df.loc[ii] = [i[0],i[1],i[2],val]
    ii += 1
df.to_excel(writer,sheet_name='Tc',index=None)
writer.close()

filename = 'X.xlsx'
ii = 0
writer = pd.ExcelWriter(filename,engine='xlsxwriter')
df = pd.DataFrame(columns=['i','j','k','val'])
for i in plant.X.keys():
    val = pyo.value(plant.X[i])
    #print("X[{}] = {}".format(i,val))
    df.loc[ii] = [i[0],i[1],i[2],val]
    ii += 1
df.to_excel(writer,sheet_name='X',index=None)
writer.close()

filename = 'dCdt.xlsx'
ii = 0
writer = pd.ExcelWriter(filename,engine='xlsxwriter')
df = pd.DataFrame(columns=['n','i','k','val'])
for i in plant.dCdt.keys():
    val = pyo.value(plant.dCdt[i])
    #print("dCdt[{}] = {}".format(i,val))
    df.loc[ii] = [i[0],i[1],i[2],val]
    ii += 1
df.to_excel(writer,sheet_name='dCdt',index=None)
writer.close()

filename = 'dTdt.xlsx'
ii = 0
writer = pd.ExcelWriter(filename,engine='xlsxwriter')
df = pd.DataFrame(columns=['i','k','val'])
for i in plant.dTdt.keys():
    val = pyo.value(plant.dTdt[i])
    #print("T[{}] = {}".format(i,val))
    df.loc[ii] = [i[0],i[1],val]
    ii += 1
df.to_excel(writer,sheet_name='dTdt',index=None)
writer.close()

filename = 'dCcdt.xlsx'
ii = 0
writer = pd.ExcelWriter(filename,engine='xlsxwriter')
df = pd.DataFrame(columns=['n','i','j','k','val'])
for i in plant.dCcdt.keys():
    val = pyo.value(plant.dCcdt[i])
   # print("Cc[{}] = {}".format(i,val))
    df.loc[ii] = [i[0],i[1],i[2],i[3],val]
    ii += 1
df.to_excel(writer,sheet_name='dCcdt',index=None)
writer.close()

filename = 'dTcdt.xlsx'
ii = 0
writer = pd.ExcelWriter(filename,engine='xlsxwriter')
df = pd.DataFrame(columns=['i','j','k','val'])
for i in plant.dTcdt.keys():
    val = pyo.value(plant.dTcdt[i])
   # print("Tc[{}] = {}".format(i,val))
    df.loc[ii] = [i[0],i[1],i[2],val]
    ii += 1
df.to_excel(writer,sheet_name='dTcdt',index=None)
writer.close()

filename = 'dXdt.xlsx'
ii = 0
writer = pd.ExcelWriter(filename,engine='xlsxwriter')
df = pd.DataFrame(columns=['i','j','k','val'])
for i in plant.dXdt.keys():
    val = pyo.value(plant.dXdt[i])
    #print("X[{}] = {}".format(i,val))
    df.loc[ii] = [i[0],i[1],i[2],val]
    ii += 1
df.to_excel(writer,sheet_name='dXdt',index=None)
writer.close()

filename = 'dCdy.xlsx'
ii = 0
writer = pd.ExcelWriter(filename,engine='xlsxwriter')
df = pd.DataFrame(columns=['n','i','k','val'])
for i in plant.dCdy.keys():
    val = pyo.value(plant.dCdy[i])
    #print("dCdt[{}] = {}".format(i,val))
    df.loc[ii] = [i[0],i[1],i[2],val]
    ii += 1
df.to_excel(writer,sheet_name='dCdy',index=None)
writer.close()

filename = 'dCdy2.xlsx'
ii = 0
writer = pd.ExcelWriter(filename,engine='xlsxwriter')
df = pd.DataFrame(columns=['n','i','k','val'])
for i in plant.dCdy2.keys():
    val = pyo.value(plant.dCdy2[i])
    #print("dCdt[{}] = {}".format(i,val))
    df.loc[ii] = [i[0],i[1],i[2],val]
    ii += 1
df.to_excel(writer,sheet_name='dCdy2',index=None)
writer.close()

filename = 'dTdy.xlsx'
ii = 0
writer = pd.ExcelWriter(filename,engine='xlsxwriter')
df = pd.DataFrame(columns=['i','k','val'])
for i in plant.dTdy.keys():
    val = pyo.value(plant.dTdy[i])
    #print("T[{}] = {}".format(i,val))
    df.loc[ii] = [i[0],i[1],val]
    ii += 1
df.to_excel(writer,sheet_name='dTdy',index=None)
writer.close()

filename = 'dTdy2.xlsx'
ii = 0
writer = pd.ExcelWriter(filename,engine='xlsxwriter')
df = pd.DataFrame(columns=['i','k','val'])
for i in plant.dTdy2.keys():
    val = pyo.value(plant.dTdy2[i])
    #print("T[{}] = {}".format(i,val))
    df.loc[ii] = [i[0],i[1],val]
    ii += 1
df.to_excel(writer,sheet_name='dTdy2',index=None)
writer.close()

filename = 'dCcdr.xlsx'
ii = 0
writer = pd.ExcelWriter(filename,engine='xlsxwriter')
df = pd.DataFrame(columns=['n','i','j','k','val'])
for i in plant.dCcdr.keys():
    val = pyo.value(plant.dCcdr[i])
   # print("Cc[{}] = {}".format(i,val))
    df.loc[ii] = [i[0],i[1],i[2],i[3],val]
    ii += 1
df.to_excel(writer,sheet_name='dCcdr',index=None)
writer.close()

filename = 'dCcdr2.xlsx'
ii = 0
writer = pd.ExcelWriter(filename,engine='xlsxwriter')
df = pd.DataFrame(columns=['n','i','j','k','val'])
for i in plant.dCcdr2.keys():
    val = pyo.value(plant.dCcdr2[i])
   # print("Cc[{}] = {}".format(i,val))
    df.loc[ii] = [i[0],i[1],i[2],i[3],val]
    ii += 1
df.to_excel(writer,sheet_name='dCcdr2',index=None)
writer.close()

filename = 'dTcdr.xlsx'
ii = 0
writer = pd.ExcelWriter(filename,engine='xlsxwriter')
df = pd.DataFrame(columns=['i','j','k','val'])
for i in plant.dTcdr.keys():
    val = pyo.value(plant.dTcdr[i])
   # print("Tc[{}] = {}".format(i,val))
    df.loc[ii] = [i[0],i[1],i[2],val]
    ii += 1
df.to_excel(writer,sheet_name='dTcdr',index=None)
writer.close()

filename = 'dTcdr2.xlsx'
ii = 0
writer = pd.ExcelWriter(filename,engine='xlsxwriter')
df = pd.DataFrame(columns=['i','j','k','val'])
for i in plant.dTcdr2.keys():
    val = pyo.value(plant.dTcdr2[i])
   # print("Tc[{}] = {}".format(i,val))
    df.loc[ii] = [i[0],i[1],i[2],val]
    ii += 1
df.to_excel(writer,sheet_name='dTcdr2',index=None)
writer.close()

'''
#X_ = np.zeros(len(Xs[int(ny//2),-1,:]))
#X_TS = np.zeros(len(Xs[int(ny//2),-1,:]))
##    for i in range(1,len(X_)+1):
##        X_[i-1] = (-pyo.log(1-Xs[int(ny//2),-1,i-1]))**0.5
##        X_TS[i-1] = pyo.value(plant.A2[0.478766,0.00075,plant.t[i]])
##        
##    c += 1
##    plt.figure(c)
##    plt.plot(ts,X_,label="analytical")
##    plt.plot(ts,X_TS,label="Taylor Series approx.")
##    plt.legend()
##    plt.show()