# -*- coding: utf-8 -*-
"""
Created on Mon Mar  17 09:36:07 2022

@author: d4mcbrid
"""

import numpy as np

# # Reactor Design
L = 3.25              # reactor length, m
Di = 1.5#75e-3#15e-3                # reactor diameter, m
Dp = 5e-3             # particle diameter, m

# # Flexible Parameters
eb = 0.376#0.41#0.37               # reactor bed porosity, Jeschar
ec = 0.55#0.7                # particle porosity
dpore = 21.74e-9          # pore diameter, m
tau = 3.0               # particle tortuosity

#L = 0.4*2*2*1.5*1.2              # reactor length, m
#Di = 12e-3*2*2*2*2*2*2                # reactor diameter, m

Ac = np.pi*Di**2/4.0    # reactor cross-sectional area, m²


LRM = L-0.1