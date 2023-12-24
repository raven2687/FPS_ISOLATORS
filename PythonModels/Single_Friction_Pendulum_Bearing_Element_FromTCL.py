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
# Create folder to save results
##############################################################################
Carpeta='Frictional_Models_1GDL';		# Output folder
if not os.path.exists(Carpeta):
    os.makedirs(Carpeta)

##############################################################################
# Load Register
##############################################################################
Reg_Concepcion = loadmat('CONCEPCION_MAULE_2010.mat')
dt=Reg_Concepcion['Dt']['Ch1'][0][0][0][0]
Acc1=Reg_Concepcion['Acc']['Ch1'][0][0][0]*U.cm/U.sec**2
Acc2=Reg_Concepcion['Acc']['Ch2'][0][0][0]*U.cm/U.sec**2
npts=Acc1.size
Time=np.arange(0,npts*dt,dt)

##############################################################################
# Create model
##############################################################################
ops.wipe()
ops.model('basic', '-ndm', 2 ,'-ndf', 3)

# Define geometry for model
g=9.807*U.m/U.sec**2
P=160*U.Tonf
m=P/g*U.kg
masa=[m, m, 0]
# Geometry
ops.node(1,0.0,0.0)
ops.node(2,0.0,0.0,'-mass',*masa)

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
mu1=0.10;
Fr_Model=['Coulomb','VelDependent','VelPressureDep','VelDepMultiLinear']
frnTag=[1,2,3,4];

#---------------------------------------------------
# To do: begin loop
#---------------------------------------------------
ops.frictionModel(Fr_Model[0],frnTag[0], mu1)
delta_1=1*U.mm
kInit=P/delta_1
Reff=160*U.cm
ops.element('singleFPBearing',1,1,2,frnTag[0],Reff,kInit,'-P',1,'-Mz',2,'-orient',0,1,0,-1,0,0)
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


Case=Fr_Model[0]
ops.recorder('Node','-file', Carpeta+'/Node_Dsp_'+str(Case)+'.out', '-time', '-node', 2, '-dof', 1, 2, 3, 'disp')
ops.recorder('Node','-file', Carpeta+'/Node_Vel_'+str(Case)+'.out', '-time', '-node', 2, '-dof', 1, 2, 3, 'vel')
ops.recorder('Node','-file', Carpeta+'/Node_Acc_'+str(Case)+'.out', '-time', '-node', 2, '-dof', 1, 2, 3, 'accel')
ops.recorder('Node','-file', Carpeta+'/Node_AbsAcc_'+str(Case)+'.out', '-timeSeries', 2, 3, '-time', '-node', 1, 2, '-dof', 1, 2, 'accel')
ops.recorder('Node','-file', Carpeta+'/RBase_'+str(Case)+'.out', '-time', '-node', 1, '-dof', 1, 2, 3, 'reaction')

ops.recorder('Element','-file', Carpeta+'/Elmt_Frc_'+str(Case)+'.out', '-time', '-ele', 1, 'force')
ops.recorder('Element','-file', Carpeta+'/Elmt_Def_'+str(Case)+'.out', '-time', '-ele', 1, 'basicDeformation')
ops.recorder('Element','-file', Carpeta+'/Elmt_N_'+str(Case)+'.out', '-time', '-ele', 1, 'frictionModel','normalForce')
ops.recorder('Element','-file', Carpeta+'/Elmt_Vel_'+str(Case)+'.out', '-time', '-ele', 1, 'frictionModel', 'vel')
ops.recorder('Element','-file', Carpeta+'/Elmt_Ff_'+str(Case)+'.out', '-time', '-ele', 1, 'frictionModel', 'frictionForce')
ops.recorder('Element','-file', Carpeta+'/Elmt_COF_'+str(Case)+'.out', '-time', '-ele', 1, 'frictionModel', 'COF')

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


# ------------------------------
# End of analysis generation
# ------------------------------
dtAna=dt; dtMin=1.0e-8; dtMax=dtAna; tFinal=npts*dt
tCurrent = ops.getTime(); ok=0; timeu2=[0.0]; u2=[0.0]
Solve_System(dt,dtAna,dtMin,dtMax,ok,tCurrent,tFinal,Tol)

# ------------------------------
# End of analysis generation
# ------------------------------
Node_Dsp=np.loadtxt(Carpeta+'/Node_Dsp_'+str(Case)+'.out')
Node_Vel=np.loadtxt(Carpeta+'/Node_Vel_'+str(Case)+'.out')
Node_Acc=np.loadtxt(Carpeta+'/Node_Acc_'+str(Case)+'.out')
Node_AbsAcc=np.loadtxt(Carpeta+'/Node_AbsAcc_'+str(Case)+'.out')
RBase=np.loadtxt(Carpeta+'/RBase_'+str(Case)+'.out')
    
Elmt_Frc=np.loadtxt(Carpeta+'/Elmt_Frc_'+str(Case)+'.out')
Elmt_Def=np.loadtxt(Carpeta+'/Elmt_Def_'+str(Case)+'.out')
Elmt_N=np.loadtxt(Carpeta+'/Elmt_N_'+str(Case)+'.out')
Elmt_Vel=np.loadtxt(Carpeta+'/Elmt_Vel_'+str(Case)+'.out')
Elmt_Ff=np.loadtxt(Carpeta+'/Elmt_Ff_'+str(Case)+'.out')
Elmt_COF=np.loadtxt(Carpeta+'/Elmt_COF_'+str(Case)+'.out')
    
fig = plt.figure()
sub1=plt.subplot(3,1,1)
sub1.set_title('Input: Concepción Maule 2010',fontsize='small')
sub1.set_ylabel('Acc \n'+r'$[cm/seg^2]$',fontsize='small')
sub1.set_xlabel('Time [seg]',fontsize='small')
sub1.tick_params(labelsize='small')
sub1.grid(True)
sub1.plot(Time,Acc1)

sub2=plt.subplot(3,1,2)
sub2.set_ylabel('Axial \n Force \n'+r'$[kN]$',fontsize='small')
sub2.set_xlabel('Time [seg]',fontsize='small')
sub2.tick_params(labelsize='small')
sub2.grid(True)
sub2.plot(Elmt_N[:,0],Elmt_N[:,1]/1000)

sub3=plt.subplot(3,1,3)
sub3.set_ylabel('Displ. \n'+r'$[mm]$',fontsize='small')
sub3.set_xlabel('Time [seg]',fontsize='small')
sub1.tick_params(labelsize='small')
sub3.grid(True)
sub3.plot(Node_Dsp[:,0],Node_Dsp[:,2])



fig2 = plt.figure()
sub4=plt.subplot(5, 4,(13,15))
sub4.set_xlabel('Longitudinal displacement'+r'$[mm]$',fontsize='small')
sub4.set_ylabel('Long. Force \n'+r'$[kN]$',fontsize='small')
sub4.tick_params(labelsize='small')
sub4.grid(True)


sub5=plt.subplot(5, 4,16)
sub5.set_ylabel('Long. displacement'+r'$[mm]$',fontsize='small')
sub5.set_xlabel('Lateral displacement'+r'$[mm]$',fontsize='small')
sub5.tick_params(labelsize='small')
sub5.grid(True)


sub6=plt.subplot(5, 4,(17,19))
sub6.set_xlabel('Lateral displacement'+r'$[mm]$',fontsize='small')
sub6.set_ylabel('Lateral Force \n'+r'$[kN]$',fontsize='small')
sub6.tick_params(labelsize='small')
sub6.grid(True)


sub7=plt.subplot(5, 4,20)
sub7.set_ylabel('Long. force'+r'$[kN]$',fontsize='small')
sub7.set_xlabel('Lateral force'+r'$[kN]$',fontsize='small')
sub7.tick_params(labelsize='small')
sub7.grid(True)










