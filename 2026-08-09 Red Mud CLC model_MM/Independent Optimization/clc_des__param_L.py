# -*- coding: utf-8 -*-
"""
Created on Mon Mar  17 09:36:07 2022

@author: d4mcbrid
"""

import numpy as np

# # Reactor Design
L = 3.25              # reactor length, m (note: this is a placeholder since L/D are decision variables in the optimizaion)
Di = 1.5             # reactor diameter, m
Dp = 5e-3             # particle diameter, m

# # Flexible Parameters
eb = 0.376              # reactor bed porosity, Mueller
ec = 0.55            # particle porosity
dpore = 21.74e-9          # pore diameter, m
tau = 3.0               # particle tortuosity


Ac = np.pi*Di**2/4.0    # reactor cross-sectional area, m²


LRM = L-0.1