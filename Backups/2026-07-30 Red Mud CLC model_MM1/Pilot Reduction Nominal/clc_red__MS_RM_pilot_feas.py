# -*- coding: utf-8 -*-
"""
Created on Mon Mar  17 09:36:45 2025

@author: d4mcbrid
"""
# This model version considers three separate sections of the reactor
# Quartz sand -> red mud -> quartz sand 

import numpy as np
import pyomo.environ as pyo
import pyomo.dae as pdae
import clc_des__param_L as dp
import clc_gen__param as gp
import clc_red__param as rp
import time
import matplotlib.pyplot as plt
import os
import csv
from pyomo.util.infeasible import log_infeasible_constraints
import pandas as pd
import scipy.integrate
from X_avg_calc import _X_avg_

        
# 0: CH4, 1: H2, 2: CO, 3: CO2, 4: H2O, 5: inert
# CLC Reduction Pyomo Model (Multiscale)

# Substitution for scaling: y=z/L (so y=[0,1], z=[0,1]L); Also, t=[0,1], tf is total time (similarly multiplied)

def reduction_pyo(init, disc, nd, t_solve, yis, cFe2O3):
  
    [nt, nint, mcp, nr] = disc
    ny = nint*mcp

    # If the model is struggling, it can help to initialize the states with a previous solution
    # I used pandas/excel to save some solutions at the end of this code file
    
    # filename = 'C.xlsx'
    # C = pd.read_excel(filename, index_col=[0,1,2]) # 0: el, 1: y, 2: t
    # filename = 'Cc.xlsx'
    # Cc = pd.read_excel(filename, index_col=[0,1,2,3]) # 0: el, 1: y, 2: r, 3: t
    # filename = 'T.xlsx'
    # T = pd.read_excel(filename, index_col=[0,1]) # 0: y, 1: t
    # filename = 'Tc.xlsx'
    # Tc  = pd.read_excel(filename, index_col=[0,1,2]) # 0: y, 1: r, 2: t
    # filename = 'X.xlsx'
    # X  = pd.read_excel(filename, index_col=[0,1,2]) # 0: y, 1: r, 2: t
    # filename = 'dCdt.xlsx'
    # dCdt = pd.read_excel(filename, index_col=[0,1,2]) # 0: el, 1: y, 2: t
    # filename = 'dTdt.xlsx'
    # dTdt = pd.read_excel(filename, index_col=[0,1]) # 0: y, 1: t
    # filename = 'dCcdt.xlsx'
    # dCcdt = pd.read_excel(filename, index_col=[0,1,2,3]) # 0: el, 1: y, 2: r, 3: t
    # filename = 'dTcdt.xlsx'
    # dTcdt  = pd.read_excel(filename, index_col=[0,1,2]) # 0: y, 1: r, 2: t
    # filename = 'dXdt.xlsx'
    # dXdt  = pd.read_excel(filename, index_col=[0,1,2]) # 0: y, 1: r, 2: t
    # filename = 'dCdy.xlsx'
    # dCdy = pd.read_excel(filename, index_col=[0,1,2]) # 0: el, 1: y, 2: t
    # filename = 'dCdy2.xlsx'
    # dCdy2 = pd.read_excel(filename, index_col=[0,1,2]) # 0: el, 1: y, 2: t
    # filename = 'dTdy.xlsx'
    # dTdy = pd.read_excel(filename, index_col=[0,1]) # 0: y, 1: t
    # filename = 'dTdy2.xlsx'
    # dTdy2 = pd.read_excel(filename, index_col=[0,1]) # 0: y, 1: t
    # filename = 'dCcdr.xlsx'
    # dCcdr = pd.read_excel(filename, index_col=[0,1,2,3]) # 0: el, 1: y, 2: r, 3: t
    # filename = 'dCcdr2.xlsx'
    # dCcdr2 = pd.read_excel(filename, index_col=[0,1,2,3]) # 0: el, 1: y, 2: r, 3: t
    # filename = 'dTcdr.xlsx'
    # dTcdr  = pd.read_excel(filename, index_col=[0,1,2]) # 0: y, 1: r, 2: t
    # filename = 'dTcdr2.xlsx'
    # dTcdr2  = pd.read_excel(filename, index_col=[0,1,2]) # 0: y, 1: r, 2: t
    
    

    # Define initial guesses
    def _C_ic_(m,n,i,k):
        return init['C'][n,int(round(i*ny,0)),int(round(k*nt,0))]
        #return C.loc[n].loc[i].loc[k].val

    def _T_ic_(m,i,k):
        return init['T'][int(round(i*ny,0)),int(round(k*nt,0))]
        #return T.loc[i].loc[k].val

    def _Cc_ic_(m,n,i,j,k):
        return init['Cc'][n,int(round(i*ny,0)),int(round(j/m.r.last()*nr,0)),int(round(k*nt,0))]
        #return Cc.loc[n].loc[i].loc[j].loc[k].val

    def _Tc_ic_(m,i,j,k):
        return init['Tc'][int(round(i*ny,0)),int(round(j/m.r.last()*nr,0)),int(round(k*nt,0))]
        #return Tc.loc[i].loc[j].loc[k].val

    def _X_ic_(m,i,j,k):
        return init['X'][int(round(i*ny,0)),int(round(j/m.r.last()*nr,0)),int(round(k*nt,0))]
        #return X.loc[i].loc[j].loc[k].val
    
    def _dCdt_ic_(m,n,i,k):
        return init['C'][n,int(round(i*ny,0)),int(round(k*nt,0))]
        #return dCdt.loc[n].loc[i].loc[k].val
    
    def _dTdt_ic_(m,i,k):
        return init['T'][int(round(i*ny,0)),int(round(k*nt,0))]
        #return dTdt.loc[i].loc[k].val
    
    def _dCcdt_ic_(m,n,i,j,k):
        return init['Cc'][n,int(round(i*ny,0)),int(round(j/m.r.last()*nr,0)),int(round(k*nt,0))]
        #return dCcdt.loc[n].loc[i].loc[j].loc[k].val

    def _dTcdt_ic_(m,i,j,k):
        return init['Tc'][int(round(i*ny,0)),int(round(j/m.r.last()*nr,0)),int(round(k*nt,0))]
        #return dTcdt.loc[i].loc[j].loc[k].val
    
    def _dXdt_ic_(m,i,j,k):
        return init['X'][int(round(i*ny,0)),int(round(j/m.r.last()*nr,0)),int(round(k*nt,0))]
        #return dXdt.loc[i].loc[j].loc[k].val
    
    def _dCdy_ic_(m,n,i,k):
        return init['C'][n,int(round(i*ny,0)),int(round(k*nt,0))]
        #return dCdy.loc[n].loc[i].loc[k].val
    
    def _dCdy2_ic_(m,n,i,k):
        return init['C'][n,int(round(i*ny,0)),int(round(k*nt,0))]
        #return dCdy2.loc[n].loc[i].loc[k].val
    
    def _dTdy_ic_(m,i,k):
        return init['T'][int(round(i*ny,0)),int(round(k*nt,0))]
        #return dTdy.loc[i].loc[k].val
    
    def _dTdy2_ic_(m,i,k):
        return init['T'][int(round(i*ny,0)),int(round(k*nt,0))]
        #return dTdy2.loc[i].loc[k].val
    
    def _dCcdr_ic_(m,n,i,j,k):
        return init['Cc'][n,int(round(i*ny,0)),int(round(j/m.r.last()*nr,0)),int(round(k*nt,0))]
        #return dCcdr.loc[n].loc[i].loc[j].loc[k].val
    
    def _dCcdr2_ic_(m,n,i,j,k):
        return init['Cc'][n,int(round(i*ny,0)),int(round(j/m.r.last()*nr,0)),int(round(k*nt,0))]
        #return dCcdr2.loc[n].loc[i].loc[j].loc[k].val
    
    def _dTcdr_ic_(m,i,j,k):
        return init['Tc'][int(round(i*ny,0)),int(round(j/m.r.last()*nr,0)),int(round(k*nt,0))]
        #return dTcdr.loc[i].loc[j].loc[k].val
    
    def _dTcdr2_ic_(m,i,j,k):
        return init['Tc'][int(round(i*ny,0)),int(round(j/m.r.last()*nr,0)),int(round(k*nt,0))]
        #return dTcdr2.loc[i].loc[j].loc[k].val
    
    # Define model (for optimization)
    m = pyo.ConcreteModel()

    ###############################################################################
    # Define parameters
    ###############################################################################

    # Note: reorganizing such that parameters are retrieved from an external source, to assist with consistency if I need to make changes
    # Reactor Design
    m.tf = pyo.Param(initialize=t_solve)                    # interval duration, s
    m.L = pyo.Param(initialize=dp.L)                        # reactor length, m
    m.LRM = pyo.Param(initialize=3.15)                      # length of red mud section of the bed, m (note the actual RM section is dictated by the discretization elements)
    m.Lq = pyo.Param(initialize=(m.L-m.LRM)/2)              # length of quartz sand section of the bed, m
    m.Di = pyo.Param(initialize=dp.Di)                      # reactor diameter, m
    m.Dp = pyo.Param(initialize=dp.Dp)                      # particle diameter, m
    m.Ac = pyo.Param(initialize=np.pi*m.Di**2/4.0)          # reactor cross-sectional area, m²

    # Reactor Conditions
    m.P = pyo.Param(initialize=rp.P)                        # inlet reactor pressure, bar
    m.T_in = pyo.Param(initialize=rp.T_in)                  # inlet temperature, K
    m.Tref = pyo.Param(initialize=rp.T_in)                  # reference temperature for certain calculations

    # Flexible Parameters
    m.rho = pyo.Param(initialize=rp.rho)                    # density of methane, kg/m³
    m.eb = pyo.Param(initialize=dp.eb)                      # reactor bed porosity
    m.ec = pyo.Param(initialize=dp.ec)                      # particle porosity
    m.dpore = pyo.Param(initialize=dp.dpore)                # pore diameter, m
    m.cFe2O30_1 = pyo.Param(initialize=cFe2O3)                  # Fe2O3 weight fraction in OC (MeO vs support material), kg Fe2O3/kg OC
    m.cFe2O30_2 = pyo.Param(initialize=cFe2O3/gp.mw_Fe2O3*1000)   # Fe2O3 concentration, mol Fe2O3/kg OC

    # General parameters
    m.mw_Fe3O4 = pyo.Param(initialize=gp.mw_Fe3O4)                # molar mass of Fe3O4, g/mol
    m.mw_Fe2O3 = pyo.Param(initialize=gp.mw_Fe2O3)              # molar mass of Fe2O3, g/mol
    m.Rg = pyo.Param(initialize=gp.Rg)                      # universal gas constant, J/mol/K
    m.Rgb = pyo.Param(initialize=gp.Rgb)                    # universal gas constant, m³*bar/mol/K

    # Mass transfer parameters
    m.mu = pyo.Param(initialize=rp.mu)                      # dynamic viscosity of mixture, kg/m/s
    m.dij = pyo.Param(initialize=rp.dij2)                   # binary diffusivity coefficient of CH4-H2O, m²/s
    m.tau = pyo.Param(initialize=dp.tau)                    # particle tortuosity
    m.rho_s = pyo.Param(initialize=rp.rho_s)                # particle density (pure solid, not bed), Red Mud, kg/m³

    # Heat transfer parameters
    m.cp_f = pyo.Param(initialize=rp.cp_f)                  # heat capacity of fuel mixture, J/mol/K
    m.cpc = pyo.Param(initialize=rp.cpc)                    # heat capacity of fuel in the particle, J/mol/K
    m.l_e0 = pyo.Param(initialize=rp.l_e0)                  # static contribution effect of thermal conductivity, J/m/s/K
    m.cp_s = pyo.Param(initialize=rp.cp_s)                  # oxygen carrier specific heat capacity, J/kg/K
    m.l_i = pyo.Param(initialize=rp.l_i)                    # thermal conductivity of red mud, W/m/K (weighted)
  

    # Reaction Parameters
    m.a0 = pyo.Param(initialize=rp.a0)                      # initial specific surface area of OC, m²/kg
    m.n_AE = pyo.Param(initialize=rp.n_AE)                  # reaction rate constant, Deng


    # Define components set; and radial, temporal, and axial domains
    m.el = pyo.Set(initialize=[0,1,2,3,4,5])#couldbeentirelywrong; 0: CH4, 1: H2, 2: CO, 3: CO2, 4: H2O, 5: inert
    m.r = pdae.ContinuousSet(bounds=(0.0,m.Dp/2))
    m.t = pdae.ContinuousSet(initialize=[0,0.02,0.04,0.06,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.6,0.7,0.8,0.9,1])
    m.y = pdae.ContinuousSet(initialize=[0,0.015,0.03,0.05,0.075,0.1,0.125,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1])

    m.hr = pyo.Param(initialize=m.r.last()/nr)
    m.hy = pyo.Param(initialize=m.y.last()/ny)

    # Define state variables
    m.C = pyo.Var(m.el, m.y, m.t, initialize=_C_ic_)  # Concentration in fluid phase, mol/m³
    m.T = pyo.Var(m.y, m.t, initialize=_T_ic_)             # Temperature in fluid phase, K
    m.Cc = pyo.Var(m.el, m.y, m.r, m.t, initialize=_Cc_ic_)#,domain=pyo.NonNegativeReals) # Concentration in the particle, mol/m³
    m.Tc = pyo.Var(m.y, m.r, m.t, initialize=_Tc_ic_)#,bounds=(0.9,1.1))       # Temperature in the particle, K
    m.X = pyo.Var(m.y, m.r, m.t, initialize=_X_ic_,bounds=(0.01,1))         # Solid conversion in the particle

    # Define derivatives
    m.dCdt = pdae.DerivativeVar(m.C, wrt=m.t)#,initialize=_dCdt_ic_)
    m.dTdt = pdae.DerivativeVar(m.T, wrt=m.t)#,initialize=_dTdt_ic_)
    m.dCcdt = pdae.DerivativeVar(m.Cc, wrt=m.t)#,initialize=_dCcdt_ic_)
    m.dTcdt = pdae.DerivativeVar(m.Tc, wrt=m.t)#,initialize=_dTcdt_ic_)
    m.dXdt = pdae.DerivativeVar(m.X, wrt=m.t)#,initialize=_dXdt_ic_)
    m.dCdy = pdae.DerivativeVar(m.C, wrt=m.y)#,initialize=_dCdy_ic_)
    m.dCdy2 = pdae.DerivativeVar(m.C, wrt=(m.y,m.y))#,initialize=_dCdy2_ic_)
    m.dTdy = pdae.DerivativeVar(m.T, wrt=m.y)#,initialize=_dTdy_ic_)
    m.dTdy2 = pdae.DerivativeVar(m.T, wrt=(m.y,m.y))#,initialize=_dTdy2_ic_)
    m.dCcdr = pdae.DerivativeVar(m.Cc, wrt=m.r)#,initialize=_dCcdr_ic_)
    m.dCcdr2 = pdae.DerivativeVar(m.Cc, wrt=(m.r,m.r))#,initialize=_dCcdr2_ic_)
    m.dTcdr = pdae.DerivativeVar(m.Tc, wrt=m.r)#,initialize=_dTcdr_ic_)
    m.dTcdr2 = pdae.DerivativeVar(m.Tc, wrt=(m.r,m.r))#,initialize=_dTcdr2_ic_)

    # Differentiate system
    disc_fd = pyo.TransformationFactory('dae.finite_difference')
    disc_oc = pyo.TransformationFactory('dae.collocation')
    disc_fd.apply_to(m, wrt=m.t, nfe=nt, scheme='BACKWARD')
    disc_fd.apply_to(m, wrt=m.r, nfe=nr, scheme='CENTRAL')
    disc_oc.apply_to(m, wrt=m.y, nfe=nint, ncp=mcp, scheme='LAGRANGE-RADAU')



    def _dCcdr2_cent(m, n, i, j, k):
        if j == m.r.first() or j == m.r.last():
            if j == m.r.first():
                return m.dCcdr2[n,i,j,k] == (m.Cc[n,i,m.r.at(int(round(j/(m.Dp/2)*nr,0))+2+1),k]-2*m.Cc[n,i,m.r.at(int(round(j/(m.Dp/2)*nr,0))+1+1),k]+m.Cc[n,i,j,k])/m.hr**2
            if j == m.r.last():
                return m.dCcdr2[n,i,j,k] == (m.Cc[n,i,j,k]-2*m.Cc[n,i,m.r.at(int(round(j/(m.Dp/2)*nr,0))-1+1),k]+m.Cc[n,i,m.r.at(int(round(j/(m.Dp/2)*nr,0))-2+1),k])/m.hr**2
        else:
            return pyo.Constraint.Skip
    m.dCcdr2_centDiff = pyo.Constraint(m.el, m.y, m.r, m.t, rule=_dCcdr2_cent)

    def _dTcdr2_cent(m, i, j, k):
        if j == m.r.first() or j == m.r.last():
            if j == m.r.first():
                return m.dTcdr2[i,j,k] == (m.Tc[i,m.r.at(int(round(j/m.r.last()*nr,0))+2+1),k]-2*m.Tc[i,m.r.at(int(round(j/m.r.last()*nr,0))+1+1),k]+m.Tc[i,j,k])/m.hr**2
            if j == m.r.last():
                return m.dTcdr2[i,j,k] == (m.Tc[i,j,k]-2*m.Tc[i,m.r.at(int(round(j/m.r.last()*nr,0))-1+1),k]+m.Tc[i,m.r.at(int(round(j/m.r.last()*nr,0))-2+1),k])/m.hr**2
        else:
            return pyo.Constraint.Skip
    m.dTcdr2_centDiff = pyo.Constraint(m.y, m.r, m.t, rule=_dTcdr2_cent)


    ###############################################################################
    # Bulk mass transfer variables
    ###############################################################################

    # Inlet compositions
    def _yis_(m, n):
        return yis[n]
    m.yis = pyo.Param(m.el, initialize=_yis_)

    # Outlet compositions
    def _yout_(m, n, k):
        return m.C[n,m.y.last(),k]/sum(m.C[nn,m.y.last(),k] for nn in m.el)
    m.yout = pyo.Expression(m.el, m.t, rule=_yout_)
    
    def _youtdry_(m,n,k): # dry outlet fraction (no water)
        return m.C[n,m.y.last(),k]/(m.C[0,m.y.last(),k]+m.C[1,m.y.last(),k]+m.C[2,m.y.last(),k]+m.C[3,m.y.last(),k]+m.C[5,m.y.last(),k])
    m.youtdry = pyo.Expression(m.el,m.t,rule=_youtdry_)
    # Inlet methane mass flux; kg/m2/s
    m.G = pyo.Param(initialize=rp.G*0.9) #note rp.G*0.9 is the baseline flow rate

    # 0: CH4, 1: H2, 2: CO, 3: CO2, 4: H2O, 5: inert
    def _MWs_(m, n):
        if n == 0:
            return gp.mw_CH4
        elif n == 1:
            return gp.mw_H2
        elif n == 2:
            return gp.mw_CO
        elif n == 3:
            return gp.mw_CO2
        elif n == 4:
            return gp.mw_H2O
        elif n == 5:
            return gp.mw_N2
    m.MWs = pyo.Param(m.el, initialize=_MWs_)
   
  
    def _y_profiles_(m, n, i, k):
        return m.C[n,i,k]/sum(m.C[nn,i,k] for nn in m.el)
    m.y_profiles = pyo.Expression(m.el, m.y, m.t, rule=_y_profiles_)


    m.Dmi = pyo.Expression(expr=m.dij)                          # diffusion coefficient of the mixture, m²/s
    m.Sc = pyo.Expression(expr=m.mu/(m.rho*m.Dmi))              # Schmidt number 
    m.av = pyo.Expression(expr=6*(1-m.eb)/m.Dp)                 # external particle surface area per unit volume, 1/m
    m.trace = pyo.Param(initialize=1e-4)
    def _y0_(m,n):
        if n == 0 :
            return m.trace
        elif n == 1:
            return m.trace
        elif n == 2:
            return m.trace
        elif n == 3:
            return m.trace
        elif n == 4:
            return m.trace
        elif n == 5:
            return 1-5*m.trace
    m.y0 = pyo.Param(m.el, initialize=_y0_)
    m.CT = pyo.Param(initialize=rp.P/(gp.Rgb*rp.T_in))
    # Gas velocity, m/s
    def _vel_(m, k):
        return m.G/m.rho
    m.vel = pyo.Expression(m.t, rule=_vel_)

    # Volumetric flowrate (total), m³/s
    def _Q_(m, k):
        return m.G*m.Ac/m.rho
    m.Q = pyo.Expression(m.t, rule=_Q_)

    # CO2 Selectivity
#    def _SCO2_(m, k): # mol/m3 * m3/s ; kg/m2/s -> kg/s -> g/s -> mol/s
#        return (m.C[3,m.y.last(),k]*m.Q[k]) / (m.yis[0]*m.G[k]*m.Ac*1000/gp.mw_CH4 )
#    m.SCO2 = pyo.Expression(m.t, rule=_SCO2_)
    def _SCO2_(m, k): # should work out to the same as CO2 outlet concentration / CH4 inlet concentration
        return m.C[3,m.y.last(),k] / m.C[0,m.y.first(),k]
    m.SCO2 = pyo.Expression(m.t, rule=_SCO2_)

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



    ###############################################################################
    # Bulk heat transfer variables
    ###############################################################################

    # Effective axial thermal conductivity, W/m/K
    def _l_ax_(m, k):
        return m.l_e0+0.7*(m.cp_f*1000/25.6)*(m.vel[k]/(1-m.eb)*m.Dp*m.rho) #25.6
    m.l_ax = pyo.Expression(m.t, rule=_l_ax_)

    def _hf_(m, k):
        return 1.37*0.357*m.Re[k]**(-0.359)*m.Sc**(-2/3)*(m.cp_f*1000/25.6)*m.G/m.eb # note m.cp_f*1000/25.6 converts heat capacity to be per kg instead of per mol, update if fuel composition changes!
    m.hf = pyo.Expression(m.t, rule=_hf_)


    ###############################################################################
    # OC mass transfer variables
    ###############################################################################

    def _Dei_(m, n):
        return rp.Dei[n]
    m.Dei = pyo.Param(m.el, initialize=_Dei_)


    ###############################################################################
    # Reaction rates
    ###############################################################################
    # 0: CH4, 1: H2, 2: CO, 3: CO2, 4: H2O, 5: inert
    # R1 (m1): CH4 + 12Fe2O3 -> CO2 + 2H2O + 8Fe3O4
    # R2: CO + 3Fe2O3 -> CO2 + 2Fe3O4
    # R3: H2 + 3Fe2O3 -> H2O + 2Fe3O4

    
    m.km1 = pyo.Param(initialize=6.676e-8)

    m.km2 = pyo.Param(initialize=2.45e-7)
    m.order2 = pyo.Param(initialize =  1)
    m.order1 = pyo.Param(initialize=1)
    m.km3 = pyo.Param(initialize=3.97e-7)
    m.order3 = pyo.Param(initialize=1)

    m.dH1 = pyo.Param(initialize=rp.dH1)
    m.dH2 = pyo.Param(initialize=rp.dH2)
    m.dH3 = pyo.Param(initialize=rp.dH3)

   
    # Here, not non-dimensionalizing Ps, but non-dimensionalizing reaction rates
    def _Ps_(m, n, i, j, k): #(bar)
        return m.Cc[n,i,j,k]*m.Rgb*m.Tc[i,j,k]*nd[0]
    m.Ps = pyo.Expression(m.el, m.y, m.r, m.t, rule=_Ps_)
    
    def _rm1_(m, i, j, k):
       return m.a0*m.km1*m.n_AE*(1-m.X[i,j,k])*(-pyo.log(1.000001-m.X[i,j,k]))**0.5*m.Cc[0,i,j,k]**m.order1*m.cFe2O30_1
    m.rm1 = pyo.Expression(m.y, m.r, m.t, rule=_rm1_)
    
    def _rm2_(m, i, j, k):
        return m.a0*m.km2*m.n_AE*(1-m.X[i,j,k])*(-pyo.log(1.000001-m.X[i,j,k]))**0.5*m.Cc[2,i,j,k]**m.order2*m.cFe2O30_1
    m.rm2 = pyo.Expression(m.y, m.r, m.t, rule=_rm2_)
    
    def _rm3_(m, i, j, k):
        return m.a0*m.km3*m.n_AE*(1-m.X[i,j,k])*(-pyo.log(1.000001-m.X[i,j,k]))**0.5*m.Cc[1,i,j,k]**m.order3*m.cFe2O30_1
    m.rm3 = pyo.Expression(m.y, m.r, m.t, rule=_rm3_)

    # # Overall reaction rate
    def _rtot_(m, n, i, j, k):
        if i < m.Lq/m.L or i > (m.Lq+m.LRM)/m.L: # no reaction occurs
            return 0
        elif n == 0: # CH4
            return -m.rm1[i,j,k] 
        elif n == 1: # H2
            return -m.rm3[i,j,k]
        elif n == 2: # CO
            return -m.rm2[i,j,k]
        elif n == 3: # CO2
            return m.rm1[i,j,k] + m.rm2[i,j,k]
        elif n == 4: # H2O
            return 2*m.rm1[i,j,k]+m.rm3[i,j,k]
        elif n == 5: # Inert
            return 0
    m.rtot = pyo.Expression(m.el, m.y, m.r, m.t, rule=_rtot_)

    # Conversion rate for OC particle
    def _conv_(m, i, j, k):
        if i < m.Lq/m.L or i > (m.Lq+m.LRM)/m.L:
            return m.dXdt[i,j,k] == 0
        else:
            return m.dXdt[i,j,k]*m.cFe2O30_2 == m.tf* (12*m.rm1[i,j,k]+3*m.rm2[i,j,k]+3*m.rm3[i,j,k])*nd[0]
    m.x_conv = pyo.Constraint(m.y, m.r, m.t, rule=_conv_)
    
    # Overall heat of reaction (J/mol)
    def _dH_(m, i, j, k):
        if i < m.Lq/m.L or i > (m.Lq+m.LRM)/m.L:
            return 0
        else:
            return (m.rm1[i,j,k]*m.dH1 + m.rm2[i,j,k]*m.dH2 + m.rm3[i,j,k]*m.dH3) * nd[0] / nd[1]
    m.dH = pyo.Expression(m.y, m.r, m.t, rule=_dH_)


    ###############################################################################
    # OC heat transfer variables
    ###############################################################################

    m.l_s = pyo.Expression(expr = m.l_i)     # weighted thermal conductivity of OC, W/m/K


    ###############################################################################
    # PBR mass balance
    ###############################################################################

    def _PBR_mass(m, n, i, k):
        return 1/m.tf*m.eb*m.dCdt[n,i,k] + m.Q[k]/m.Ac*(1/m.L)*m.dCdy[n,i,k] == m.Dax[k]*(1/(m.L**2))*m.dCdy2[n,i,k] + m.kc[k]*m.av*(m.Cc[n,i,m.r.last(),k]-m.C[n,i,k])

    m.PBR_mass = pyo.Constraint(m.el, m.y, m.t, rule=_PBR_mass)

    # PBR mass balance BC; feed stream condition at inlet
    def _BC1_(m, n, k):
        if k == m.t.first():
            return (m.C[n,0,m.t.first()] - m.y0[n]*rp.P/gp.Rgb/rp.T_in/nd[0])**2 <= m.C[n,0,m.t.first()]*1e-4
        else:
            return m.C[n,0,k] == m.yis[n]*m.P/(m.Rgb*m.T_in)/nd[0]

    m.BC1 = pyo.Constraint(m.el, m.t, rule=_BC1_)

    # PBR mass balance BC; mass insulation at reactor outlet
    def _BC2_(m, n, k):
        return (m.C[n,m.y.last(),k] - m.C[n,m.y.at(-2),k])**2 <= 1e-2*m.C[n,m.y.at(-2),k]
    m.BC2 = pyo.Constraint(m.el, m.t, rule=_BC2_)


    ###############################################################################
    # PBR energy balance [possibly move up apparently]
    ###############################################################################

    def _PBR_energy(m, i, k):
        return m.eb*m.cp_f*m.P/(m.Rgb*m.T_in)*(1/m.tf)*m.dTdt[i,k] + m.cp_f*m.Q[k]*m.P/(m.Rgb*m.T_in)/m.Ac*(1/m.L)*m.dTdy[i,k] == m.l_ax[k]*((1/m.L)**2)*m.dTdy2[i,k] + m.hf[k]*m.av*(m.Tc[i,m.r.last(),k]-m.T[i,k])
    m.PBR_energy = pyo.Constraint(m.y, m.t, rule=_PBR_energy)

    # PBR energy balance BC; heat feed stream conditions at inlet
    def _BC3_(m, k):
       return (m.T_in/nd[1] - m.T[m.y.at(1),k])**2 <= 1e-6
    m.BC3 = pyo.Constraint(m.t, rule=_BC3_)
    
        
    # PBR energy balance BC; heat insulation condition at reactor's outlet stream
    def _BC4_(m, k):
        return (m.T[m.y.last(),k] - m.T[m.y.at(-2),k])**2 <= 1e-2*m.T[m.y.at(-2),k]
    m.BC4 = pyo.Constraint(m.t, rule=_BC4_)


    ###############################################################################
    # OC Particle mass balance
    ###############################################################################

    def _OC_mass(m, n, i, j, k):
        if j == m.r.first():
            return m.ec*m.dCcdt[n,i,j,k] == m.tf * (3*m.Dei[n]*(m.dCcdr2[n,i,j,k]) + m.ec*m.rho_s*m.rtot[n,i,j,k])
        return m.ec*j**2*m.dCcdt[n,i,j,k] == m.tf * (m.Dei[n]*(2*j*m.dCcdr[n,i,j,k]+m.dCcdr2[n,i,j,k]*j**2) + j**2*m.ec*m.rho_s*m.rtot[n,i,j,k])
    m.OC_mass = pyo.Constraint(m.el, m.y, m.r, m.t, rule=_OC_mass)
    
    # OC mass balance BC; mass insulation
    def _BC5_(m, n, i, k):
        return m.dCcdr[n,i,m.r.first(),k] == 0
    m.BC5 = pyo.Constraint(m.el, m.y, m.t, rule=_BC5_)

    # OC mass balance BC; mass transfer between reactor bulk phase and OC particle
    def _BC6_(m, n, i, k):
        return -m.Dei[n]*m.dCcdr[n,i,m.r.last(),k] == m.kc[k]*(m.Cc[n,i,m.r.last(),k]-m.C[n,i,k])
    m.BC6 = pyo.Constraint(m.el, m.y, m.t, rule=_BC6_)


    ###############################################################################
    # OC Particle energy balance
    ###############################################################################

    def _OC_energy(m, i, j, k):
        if j == m.r.first():
            return ((1-m.ec)*m.rho_s*m.cp_s + m.ec*m.P/(m.Rgb*m.T_in)*m.cp_f)*m.dTcdt[i,j,k] == m.tf * (3*m.l_s*m.dTcdr2[i,j,k] - m.ec*m.rho_s*m.dH[i,j,k])
        return ((1-m.ec)*m.rho_s*m.cp_s + m.ec*m.P/(m.Rgb*m.T_in)*m.cp_f)*j**2*m.dTcdt[i,j,k] == m.tf * (m.l_s*(2*j*m.dTcdr[i,j,k]+m.dTcdr2[i,j,k]*j**2) - j**2*m.ec*m.rho_s*m.dH[i,j,k])
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
        return(3.1*(150*(1-m.eb)**2/(m.eb**3)*m.mu*m.vel[k]/(m.Dp**2)+1.75*(1-m.eb)/(m.eb**3)*m.rho*m.vel[k]**2/(m.Dp))) #3.1 is LRM
    m.dP = pyo.Expression(m.t,rule=_Ergun_) #pressure drop
    m.umf = pyo.Var(initialize=0.5) # solves for minimum fluidization velocity
    def _umf_(m):
        return 1305.556*m.rho*m.umf*m.Dp/m.mu + 30.789*(m.rho*m.umf*m.Dp/m.mu)**1.657 + 19.474*(m.rho*m.umf*m.Dp/m.mu)**2-m.rho*m.Dp**3*(m.rho_s*(1-m.ec)-m.rho)*9.81/(m.mu**2) == 0
    m.u_mf = pyo.Constraint(rule=_umf_)
    m.LD = pyo.Expression(expr=m.L/m.Di)
    def _fuel_conv_(m,n,k):
        return 1-m.C[n,m.y.last(),k]/(m.yis[n]*101325/(8.314*m.T_in)/nd[0])
    m.fuel_conv = pyo.Expression(m.el,m.t,rule=_fuel_conv_)
   
    
    def Q_out(m): # Average heat extracted via outlet gas (left reimann sum)
        value = 0
        for t in range(len(m.t)):
            if t == 0:
                value += 0
            else:
                value += (m.t.at(t+1)-m.t.at(t))*m.Q[m.t.at(t)]*((m.C[0,m.y.last(),m.t.at(t)]+m.C[1,m.y.last(),m.t.at(t)]+m.C[2,m.y.last(),m.t.at(t)]+m.C[3,m.y.last(),m.t.at(t)]+m.C[4,m.y.last(),m.t.at(t)]+m.C[5,m.y.last(),m.t.at(t)])*nd[0])*m.cp_f*(m.T[m.y.last(),m.t.at(t)]*nd[1]-rp.T_in)
        return value
    def Q_in(m): # Average heat remaining in the packed bed at the end of oxidation stage (left reimann sum)
        value = 0
        for z in range(len(m.y)):
            if z == 0:
                value += 0
            else:
                value += m.L*(m.y.at(z+1)-m.y.at(z))*m.rho_s*(1-m.eb)*(1-m.ec)*m.Ac*m.cp_s*(m.Tc[m.y.at(z),m.r.last(),m.t.last()]*nd[1]-rp.T_in)
        return value/m.tf
    m.Q_in = pyo.Expression(rule=Q_in)
    m.Q_out = pyo.Expression(rule=Q_out)
    
    def Q_out_int(m,k):  # Average heat extracted via outlet gas (trapezoidal integration in pyomo for better accuracy)
        return m.tf*m.Q[k]*((m.C[0,m.y.last(),k]+m.C[1,m.y.last(),k]+m.C[2,m.y.last(),k]+m.C[3,m.y.last(),k]+m.C[4,m.y.last(),k]+m.C[5,m.y.last(),k])*nd[0])*m.cp_f*(m.T[m.y.last(),k]*nd[1]-rp.T_in)
    m.Q_out_int = pdae.integral.Integral(m.t,wrt=m.t,rule=Q_out_int)
    
    def Q_in_int(m,i):# Average heat remaining in the packed bed at the end of oxidation stage(trapezoidal integration in pyomo for better accuracy)
        return m.L*m.rho_s*(1-m.eb)*(1-m.ec)*m.Ac*m.cp_s*(m.Tc[i,m.r.last(),m.t.last()]*nd[1]-rp.T_in)
    m.Q_in_int = pdae.integral.Integral(m.y,wrt=m.y,rule=Q_in_int)
    
    m.fuel_conversion_int = pdae.integral.Integral(m.el,m.t,wrt=m.t,rule=_fuel_conv_) #average fuel conversion for each component integrated over time
    return m
if __name__ == "__main__":

    tfreq = 50                   # Frequency of time points in discretization, s
    tf = 1300                   # Reduction stage simulation time, s

        
    cFe2O3 = 0.4928                                                 # Fe2O3 concentration (mass fraction)
    

    
 
    # Define nondimensionalizing coefficients for concentration and temperature
    Cnd = rp.yis[0]*rp.P/(gp.Rgb*rp.T_in)
    Tnd = rp.T_in
    nd = [Cnd,Tnd]
     
    # Discretization
    nel = 6
    nr = 3                          # Number of radial nodes (particle)
    nint = 25                    # Number of axial elements (bulk reactor)
    mcp = 2                         # Total number of orthogonal collocation points sent to pyomo
    nt = 18         # Number of temporal nodes
    ny = nint*mcp                   # Total number of axial points
    l_t = tf/nt                 # Average length of each time node
    disc = [nt, nint, mcp, nr]

    trace = 1e-4
    yis = rp.yis                     #inlet mole fractions
    
    # Initial Conditions
    y0 = [trace, trace, trace, trace, trace, 1.-5*trace]    # initial mole fractions of different compounds
    C0 = [y0[n]*rp.P/gp.Rgb/rp.T_in for n in range(nel)]    # initial fluid-phase concentrations of compounds, mol/m³
    T0 = rp.T_in                                          # initial temperature at reactor inlet, K
    X0 = 0.01                                           # initial conversion of Fe2O3 to Fe3O4
   

    
    # Define initial guesses (and initial conditions for t=0)
    init = {'C': np.full([nel,ny+1,nt+1],0,dtype=float),
            'T': np.full([ny+1,nt+1],0,dtype=float),
            'Xp': np.full([ny+1,nt+1],0,dtype=float),
            'Cc': np.full([nel,ny+1,nr+1,nt+1],0,dtype=float),
            'Tc': np.full([ny+1,nr+1,nt+1],0,dtype=float),
            'X': np.full([ny+1,nr+1,nt+1],0,dtype=float)}
    
    # Define matrices to collect data
    plantsolves = []
    ts = [k*l_t for k in range(int(tf/l_t)+1)]
    tTs = ts[nt:]
    Cs = np.full([nel,ny+1,len(ts)],0,dtype=float)
    youts = np.full([nel,len(ts)],0,dtype=float)
    youtsdry = np.full([nel,len(ts)],0,dtype=float)
    Ts = np.full([ny+1,len(ts)],0,dtype=float)
    Ccs = np.full([nel,ny+1,nr+1,len(ts)],0,dtype=float)
    Tcs = np.full([ny+1,nr+1,len(ts)],0,dtype=float)
    Xs = np.full([ny+1,nr+1,len(ts)],0,dtype=float)
    Gs = np.full([len(ts)],0,dtype=float)
    Ssets = np.full([len(tTs)],0,dtype=float)
    SCO2s = np.full([len(ts)],0,dtype=float)
    Gsets = np.full([len(tTs)],0,dtype=float)
    GCO2s = np.full([len(ts)],0,dtype=float)
    GCH4s = np.full([len(ts)],0,dtype=float)
    zs = np.full([ny+1],0,dtype=float)
    tps = [k*l_t for k in range(nt,int((tf)/l_t)+1)]

    rtot_CO2 = np.full([ny+1,len(ts)],0,dtype=float)
    rtot_CH4 = np.full([ny+1,len(ts)],0,dtype=float)
    rtot_H2O = np.full([ny+1,len(ts)],0,dtype=float)
    rtot_1 = np.full([ny+1,len(ts)],0,dtype=float)
    rtot_2 = np.full([ny+1,len(ts)],0,dtype=float)
    rtot_3 = np.full([ny+1,len(ts)],0,dtype=float)
    

    plantopt = pyo.SolverFactory('ipopt')#, executable='C:\Jaime\MSYS2\home\ja2mccor\ipopt')#pyo.SolverFactory('ipopt')#, executable='C:\Dana\ipopt')
    plantopt.options['linear_solver'] = 'ma97'
    plantopt.options['ma97_u'] = 0.0001#0.0001##0.0001
    
    
    # Define initial guesses
    for k in range(nt+1):
        for i in range(ny+1):
            for n in range(nel):
                init['C'][n,i,k] = C0[n]/nd[0]
            init['T'][i,k] = T0/nd[1]
            init['Xp'][i,k] = X0
            for j in range(nr+1):
                for n in range(nel):
                    init['Cc'][n,i,j,k] = C0[n]/nd[0]
                init['Tc'][i,j,k] = T0/nd[1]
                init['X'][i,j,k] = X0
    
    tt = 0
    
    
    
    # Define plant model
    plant = reduction_pyo(init, disc, nd, tf, yis, cFe2O3)

    
    # Fix initial conditions for plant
    for i in range(ny+1):
        for n in range(nel):
            plant.C[n,plant.y.at(i+1),plant.t.first()].fix(init['C'][n,i,0])
        plant.T[plant.y.at(i+1),plant.t.first()].fix(init['T'][i,0])
        for j in range(nr+1):
            for n in range(nel):
                plant.Cc[n,plant.y.at(i+1),plant.r.at(j+1),plant.t.first()].fix(init['Cc'][n,i,j,0])
            plant.Tc[plant.y.at(i+1),plant.r.at(j+1),plant.t.first()].fix(init['Tc'][i,j,0])
            plant.X[plant.y.at(i+1),plant.r.at(j+1),plant.t.first()].fix(init['X'][i,j,0])

    start = time.time()
    # Solve model
  
    plantresults = plantopt.solve(plant,tee=True)
    plantsolves += [plantresults['Solver']]

##     #Compile data
    print(time.time()-start)
    
    foldername = "Pilot_red_base"
    save = False#True
    if not os.path.exists(foldername):
        os.makedirs(foldername)
    
    
    for k in range(nt+1):
       Gs = plant.G.value
       SCO2s[k] =pyo.value(plant.SCO2[plant.t.at(k+1)])
       ts[k] = tf*plant.t.at(k+1)
       for n in range(nel):
           youts[n,k] = plant.yout[n,plant.t.at(k+1)]()

           youtsdry[n,k] = plant.youtdry[n,plant.t.at(k+1)]()
       for i in range(ny+1):
           zs[i] = plant.y.at(i+1)*dp.L
           rtot_CO2[i,k] = pyo.value(plant.rtot[3,plant.y.at(i+1),plant.r.last(),plant.t.at(k+1)])
           rtot_H2O[i,k] = pyo.value(plant.rtot[4,plant.y.at(i+1),plant.r.last(),plant.t.at(k+1)])
           rtot_1[i,k] = pyo.value(plant.rm1[plant.y.at(i+1),plant.r.last(),plant.t.at(k+1)])
           rtot_2[i,k] = pyo.value(plant.rm2[plant.y.at(i+1),plant.r.last(),plant.t.at(k+1)])
           rtot_3[i,k] = pyo.value(plant.rm3[plant.y.at(i+1),plant.r.last(),plant.t.at(k+1)])
           for n in range(nel):
               Cs[n,i,k] = plant.C[n,plant.y.at(i+1),plant.t.at(k+1)].value*nd[0]
           Ts[i,k] = plant.T[plant.y.at(i+1),plant.t.at(k+1)].value*nd[1]
           for j in range(nr+1):
               for n in range(nel):
                   Ccs[n,i,j,k] = plant.Cc[n,plant.y.at(i+1),plant.r.at(j+1),plant.t.at(k+1)].value*nd[0]
               Tcs[i,j,k] = plant.Tc[plant.y.at(i+1),plant.r.at(j+1),plant.t.at(k+1)].value*nd[1]
               Xs[i,j,k] = plant.X[plant.y.at(i+1),plant.r.at(j+1),plant.t.at(k+1)].value
    l1_CO2 = youtsdry[3,:]
    l1_CH4 = youtsdry[0,:]
    l1_CO = youtsdry[2,:]
    l1_H2 = youtsdry[1,:]
    l1_X = Xs[:,-1,-1]
    l1_T = Ts[:,-1]
    l2 = ts
    l3 = zs
    l3_CO2 = Cs[3,:,-1]
    l3_N2 = Cs[5,:,-1]
    l3_CO = Cs[2,:,-1]
    l3_H2 = Cs[1,:,-1]
    l3_H2O = Cs[4,:,-1]
    l3_CH4 = Cs[0,:,-1]
    df_t = pd.DataFrame({'time':l2,'CO2 fraction':l1_CO2, 'CH4 fraction':l1_CH4,'CO fraction':l1_CO,'H2 fraction':l1_H2})
    df_l = pd.DataFrame({'length':l3,'final conversion':l1_X,'Midtime temp':l1_T,'final CO2 concentration':l3_CO2,'final N2 concentration':l3_N2,'final CO concentration':l3_CO,'final H2 concentration':l3_H2,'final CH4 concentration':l3_CH4,'final H2O concentration':l3_H2O})
    df_t.to_excel('results_time.xlsx')
    df_l.to_excel('results_length.xlsx')


    data = {'value':[[round(pyo.value(plant.fuel_conversion_int[0]),3),round(pyo.value(plant.fuel_conversion_int[1]),3),round(pyo.value(plant.fuel_conversion_int[2]),3)],pyo.value(plant.vel[plant.t.last()]),pyo.value(plant.umf*plant.Dp*plant.rho/plant.mu),round(pyo.value(9.81*plant.Dp**3*plant.rho*(plant.rho_s-plant.rho)/(plant.mu**2))),pyo.value(plant.LD),tf/60,plant.Dp/(plant.Di/2),pyo.value(plant.dP[plant.t.last()]),pyo.value(plant.Di/plant.Dp)],'restriction':['>0.85',pyo.value(plant.umf),'<350','<5000000','2-6','15 min','0.04-0.25','<100000','>2']}
    print(pd.DataFrame(data,index=['Average conversion (CH4, H2, CO)', 'Gas velocity','Remf','Ar','L/D','Cycle time','dp/R0','Pressure drop %','D/dp']))
         
    filename = 'X'+foldername+'.xlsx'
    ii = 0
    writer = pd.ExcelWriter(filename)
    df = pd.DataFrame(columns=['i','j','k','val'])
    for i in plant.X.keys():
        val = pyo.value(plant.X[i])
        #print("X[{}] = {}".format(i,val))
        df.loc[ii] = [i[0],i[1],i[2],val]
        ii += 1
    df.to_excel(writer,sheet_name='X',index=None)
    writer.close()
    
    # calculate average conversion
    zs_int = np.array([ 0.02    , 0.03    , 0.036667,
           0.05    , 0.058333, 0.075   , 0.083333, 0.1     , 0.108333,
           0.125   , 0.133333, 0.15    , 0.166667, 0.2     , 0.216667,
           0.25    , 0.266667, 0.3     , 0.316667, 0.35    , 0.358333,
           0.375   , 0.383333, 0.4     , 0.416667, 0.45    , 0.466667,
           0.5     , 0.516667, 0.55    , 0.566667, 0.6     , 0.616667,
           0.65    , 0.666667, 0.7     , 0.716667, 0.75    , 0.766667,
           0.8     , 0.816667, 0.85    , 0.866667, 0.9     , 0.916667,
           0.95    , 0.966667])
    r_int = np.array([0,0.000833,0.001667,0.0025])
    X_average = _X_avg_(filename,zs_int,r_int,dp.Dp/2,3.1)
   
   # Model checks to make sure there's no unreasonable results
    print('min/max conc',np.min(Cs),np.max(Cs))
    print('min/max T',np.min(Ts),np.max(Ts))
    print(Xs[:,-1,-1])
    
    print('CH4 max (theoretical): 0.9214','model:',np.max(Cs[0]))
    print('H2 max (theoretical): 1.51','model:',np.max(Cs[1]))
    print('CO max (theoretical): 3.63','model:',np.max(Cs[2]))
    print('CO2 max (theoretical): 6.07','model:',np.max(Cs[3]))
    print('H2O max (theoretical): 3.36','model:',np.max(Cs[4]))
    x = zs[3:-1]-zs[3]      # RM section (change indices if discretization or LRM changes)
    y = Xs[3:-1,-1,-1]       # conversion in the RM section
    I2 = scipy.integrate.simpson(y,x) 
    print('average conversion of the OC: '+ str(X_average))
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
              f.write(str(pd.DataFrame(data,index=['Average conversion (CH4, H2, CO)', 'Gas velocity','Remf','Ar','L/D','Cycle time','dp/R0','Pressure drop %','D/dp'])))
              f.write("\nAverage OC conversion: "+str(X_average))
              f.write('\nMin/Max conc: '+str(np.min(Cs))+'/'+str(np.max(Cs)))
              f.write('\nMin/Max temp: '+str(np.min(Ts))+'/'+str(np.max(Ts)))
              f.write('\nExtracted heat (Reimann/trapezoid), J: '+str(pyo.value(plant.Q_out*plant.tf))+'/'+str(pyo.value(plant.Q_out_int)))
              f.write('\nAccumulated heat (Reimann/trapezoid), J: '+str(pyo.value(plant.Q_in*plant.tf))+'/'+str(pyo.value(plant.Q_in_int)))
              f.write('\ntotal heat (Reimann/trapezoid): '+str((pyo.value(plant.Q_out*plant.tf) + pyo.value(plant.Q_in*plant.tf))/1000) +'/'+str((pyo.value(plant.Q_out_int)+pyo.value(plant.Q_in_int))/1000)+' kJ')
              f.write('\n fuel comp: '+str(rp.yis))
              f.write('\n cFe2O3: '+str(cFe2O3))
              f.write('\nDesign Variables:')
              f.write('\n  tf: '+str(pyo.value(plant.tf)))
              f.write('\n  G (kg/m2/s): '+str(pyo.value(plant.G)))
              f.write('\n  L: 3.1 (rounded)')
              f.write('\n  Di: '+str(pyo.value(plant.Di)))
              f.write('\n Solve details:')
              
    plt.figure()
    plt.plot(zs,Ts[:,-1],label='end')
    plt.plot(zs,Ts[:,10],label='middle')
    plt.plot(zs,Ts[:,1],label='start')
    plt.xlabel('Length (m)')
    plt.ylabel('Temperature (K)')
    plt.legend()
    
    plt.figure()
    plt.plot(ts,Ts[-1,:])
    plt.xlabel('Time (s)')
    plt.ylabel('Outlet temperature (K)')
    
    plt.figure()
    plt.plot(ts,youtsdry[0],label='CH4')
    plt.plot(ts,youtsdry[1],label='H2')
    plt.plot(ts,youtsdry[2],label='CO')
    plt.plot(ts,youtsdry[3],label='CO2')
    plt.xlabel('Time (s)')
    plt.ylabel('Outlet fraction')
    plt.legend()
    plt.show()
    
         
