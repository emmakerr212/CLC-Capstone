# -*- coding: utf-8 -*-
"""
Created on Mon Mar  17 09:36:07 2022

@author: d4mcbrid
"""

import numpy as np

# # Reactor Design
L = 3.25              # reactor length, m  (note: this is longer than the actual length. I've used a small section of inert at inlet/outlet as this helped with model convergence. Actual length of reactor is determined by LRM and discrete elements)
Di = 1.5              # reactor diameter, m
Dp = 5e-3             # particle diameter, m

# # Flexible Parameters
eb = 0.376              # reactor bed porosity, Jeschar
ec = 0.55               # particle porosity
dpore = 21.74e-9          # pore diameter, m
tau = 3.0               # particle tortuosity



Ac = np.pi*Di**2/4.0    # reactor cross-sectional area, m²


LRM = L-0.1 # LRM is actually ~3.1 m based on the discretization