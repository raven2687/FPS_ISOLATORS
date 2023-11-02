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
from scipy.io import loadmat

##############################################################################
# Load Register
##############################################################################
Reg_Concepcion = loadmat('CONCEPCION_MAULE_2010.mat')

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

#---------------------------------------------------
# Define material models
kvc=2.1e6*U.kgf/U.cm**2
xi=0.02
cv=2*xi*np.sqrt(kvc/m)
ops.uniaxialMaterial('Elastic',1,kvc,cv);
ops.uniaxialMaterial('Elastic',2,0.0);

#---------------------------------------------------
# Define frictional model
mu1=0.012;
Fr_Model=['Coulomb','VelDependent','VelPressureDep','VelDepMultiLinear']
frnTag=[1,2,3,4];

#---------------------------------------------------
# To do: begin loop
#---------------------------------------------------
ops.frictionModel(Fr_Model[0],frnTag[0], mu1)

#---------------------------------------------------
# Applying static vertical load
tsTag=1
ops.timeSeries('Linear',tsTag)
ptTag=1
ops.pattern('Plain',ptTag,tsTag)
nodeTag=2
ops.load(nodeTag,0.,-P,0,)

#----------------------------------------------
# Start of analysis generation
#----------------------------------------------

# create the DOF numberer
ops.numberer('Plain')
# create the system of equation
ops.system('BandGeneral')
# create the constraint handler
ops.constraints('Plain')
# create the convergence test
ops.test('NormDispIncr',1e-12,100)
# create the solution algorithm
ops.algorithm('Newton');
# create the integration scheme
ops.integrator('LoadControl',0.1)
# create the analysis object
ops.analysis('Static')
#----------------------------------------------
# End of analysis generation
#----------------------------------------------

# perform the gravity load analysis, requires 10 steps to reach the load level
ops.analyze(10)
# set the gravity loads to be constant & reset the time in the domain
ops.loadConst('-time',0.0);





