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
from OpenSees_Solvers import Solve_System

##############################################################################
# Load Register
##############################################################################
Reg_Concepcion = loadmat('CONCEPCION_MAULE_2010.mat')
dt=Reg_Concepcion['Dt']['Ch1'][0][0][0][0]
Acc1=Reg_Concepcion['Acc']['Ch1'][0][0][0]*U.cm/U.sec**2
Acc2=Reg_Concepcion['Acc']['Ch2'][0][0][0]*U.cm/U.sec**2

npts=Acc1.size
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

# --------------------------------
# Perform an eigenvalue analysis
# --------------------------------
lambda_i = ops.eigen('-fullGenLapack', 2)
omega=np.sqrt(lambda_i)
Ti=2*np.pi/omega
fi=1/Ti

# ------------------------------
# Start of model generation
# ------------------------------

# Define dynamic loads
# --------------------
# set time series to be passed to uniform excitation
Factor = 1
tsTag=2
ops.timeSeries('Path', tsTag,'-dt',dt, '-values',*Acc1, '-factor', Factor) 
tsTag=3
ops.timeSeries('Path', tsTag,'-dt',dt, '-values',*Acc2, '-factor', Factor) 
# Define where and how (pattern tag, dof) acceleration is applied
ops.pattern("UniformExcitation", 2, 1, "-accel", 2)
ops.pattern("UniformExcitation", 3, 2, "-accel", 3)

# calculate the Rayleigh damping factors for nodes & elements
alphaM=0.05;    # mass proportional damping;       D = alphaM*M
betaK=0.0;      # stiffness proportional damping;  D = betaK*Kcurrent
betaKinit=0.0;  # stiffness proportional damping;  D = beatKinit*Kinit
betaKcomm=0.0;  # stiffness proportional damping;  D = betaKcomm*KlastCommit
ops.rayleigh(alphaM, betaK, betaKinit, betaKcomm)

# ------------------------------
# Start of analysis generation
# ------------------------------
# create the DOF numberer
ops.numberer('Plain')
# create the system of equation
ops.system('BandGeneral')
# create the constraint handler
ops.constraints('Plain')
# create the convergence test
Tol=1.0e-12
ops.test('NormDispIncr',Tol,100)

# create the solution algorithm
ops.algorithm('Newton');
# create the integration scheme
ops.integrator('Newmark',0.5,0.25)
# create the analysis object
ops.analysis('Transient')


ok = ops.analyze(npts, dt)

if (ok != 0):
    print("analysis FAILED")
else:
    print("analysis SUCCESSFUL")
    
    
u2=ops.nodeDisp(2,1)

# ------------------------------
# End of analysis generation
# ------------------------------
dtAna=dt; dtMin=1.0e-8; dtMax=dtAna; tFinal=npts*dt
tCurrent = ops.getTime(); ok=0; timeu2=[0.0]; u2=[0.0]
Solve_System(dt,dtAna,dtMin,dtMax,ok,tCurrent,tFinal,Tol,timeu2,u2)












