# -*- coding: utf-8 -*-
"""
Created on Wed Nov  1 01:33:31 2023

@author: Paul Villavicencio Fernandez
"""

# %reset -f


##############################################################################
# IMPORT LIBRARIES
##############################################################################
import numpy as np
import openseespy.opensees as ops
import matplotlib.pyplot as plt
import os
from units import SI_m_units as U

##############################################################################
# Create model
##############################################################################
ops.wipe()

ops.model('basic', '-ndm', 2 ,'-ndf', 3)

# Define geometry for model
g=9.807*U.m/U.sec**2
P=80*U.kN
m=P/g*U.kg

# Geometry
ops.node(1,0.0,0.0)
ops.node(2,0.0,0.0)

# Applying kinematic boundary condition
ops.fix(1,1,1,1);
ops.fix(2,0,0,1);

# Define material models



