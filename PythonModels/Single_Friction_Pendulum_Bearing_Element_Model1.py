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

plt.close('all')

##############################################################################
# Create folder to save results
##############################################################################
def Plot_max(Xvalues,Yvalues,Label,Elem,size):
    Ymax=np.max(abs(Yvalues))
    Ypos=np.argmax(abs(Yvalues))
    Xmax=Xvalues[Ypos]
    Ymax=Yvalues[Ypos]
    Elem.text(Xmax,Ymax,Label+str(np.char.mod('%3.3f',Ymax)),verticalalignment='center')
    Elem.scatter(Xmax,Ymax,size)
    return Ymax,Xmax
##############################################################################
# Create folder to save results
##############################################################################
Carpeta='Frictional_Models_1GDL_Model1';		# Output folder
if not os.path.exists(Carpeta):
    os.makedirs(Carpeta)

##############################################################################
# Sinusoidal Input
##############################################################################
g=9.807*U.m/U.sec**2
GMdirection=1;			    # ground-motion direction
GMSineAccAmpl=0.4*g;	    # sine ground-motion acceleration amplitude
TPeriodSine=1;	            # period of input sine wave
DurationSine=3;	            # duration of input sine wave
DtAnalysis=0.01;	        # time-step Dt for lateral analysis
TmaxAnalysis=10;	        # maximum duration of ground-motion analysis
IDloadTag=400;	            # for uniformSupport excitation
DtGround=0.005;	            # time-step Dt for input grond motion
omegaSine=2*np.pi/TPeriodSine;
vel0=GMSineAccAmpl*(-1)/omegaSine;
    
##############################################################################
# Create model
##############################################################################
ops.wipe()
ops.model('basic', '-ndm', 2 ,'-ndf', 3)

# Define model loads
m=100*U.Ton;
P=m*g;
# Geometry
ops.node(1,0.0,0.0)
ops.node(2,0.0,0.0)
ops.mass(2,*[m, m,0.0])

# Applying kinematic boundary condition
ops.fix(1,1,1,1);
ops.fix(2,0,0,1);

#---------------------------------------------------
# Define material models

Ev=2.1e6*U.kgf/U.cm**2;
Av=(250*U.mm)**2*np.pi/4;
Lv=50*U.mm

kvc=Ev*Av/Lv;

xi=0.00
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
Reff=155*U.cm
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
Tol=1.0e-0
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
    

# ------------------------------
# Plot 1
# ------------------------------
fig = plt.figure()
sub1=plt.subplot(3,2,1)
sub1.set_title('Input: Concepción Maule 2010 N-S',fontsize='small')
sub1.set_ylabel('Acc \n'+r'$[g]$',fontsize='small')
sub1.set_xlabel('Time [seg]',fontsize='small')
sub1.tick_params(labelsize='small')
sub1.grid(True)
sub1.plot(Time,Acc1/g,label=r'$\ddot{u}_{g1}$',color='k', linewidth=0.5)
sub1.legend()
Ymax,Xmax=Plot_max(Time,Acc1/g,r'$a_{max}$=',sub1,12)
sub1.set_xlim((0, tFinal))
sub1.set_ylim((-1.2*Ymax, 1.2*Ymax))



sub2=plt.subplot(3,2,2)
sub2.set_title('Input: Concepción Maule 2010 E-O',fontsize='small')
sub2.set_ylabel('Acc \n'+r'$[g]$',fontsize='small')
sub2.set_xlabel('Time [seg]',fontsize='small')
sub2.tick_params(labelsize='small')
sub2.grid(True)
sub2.plot(Time,Acc2/g,label=r'$\ddot{u}_{g2}$',color='k', linewidth=0.5)
sub2.set_xlim((0, tFinal))
sub2.legend()
Ymax,Xmax=Plot_max(Time,Acc2/g,r'$a_{max}$=',sub2,12)
sub2.set_xlim((0, tFinal))
sub2.set_ylim((-1.2*Ymax, 1.2*Ymax))



sub3=plt.subplot(3,2,(3,4))
sub3.set_ylabel('Axial \n Force \n'+r'$[kN]$',fontsize='small')
sub3.set_xlabel('Time [seg]',fontsize='small')
sub3.tick_params(labelsize='small')
sub3.grid(True)
sub3.plot(Elmt_N[:,0],Elmt_N[:,1]/1000)
Ymax,Xmax=Plot_max(Elmt_N[:,0],Elmt_N[:,1]/1000,r'$P_{max}$=',sub3,12)
sub3.set_xlim((0, tFinal))
sub3.set_ylim((-1.2*Ymax, 1.2*Ymax))

sub4=plt.subplot(3,2,(5,6))
sub4.set_ylabel('Displ. \n'+r'$[mm]$',fontsize='small')
sub4.set_xlabel('Time [seg]',fontsize='small')
sub4.tick_params(labelsize='small')
sub4.grid(True)
sub4.plot(Node_Dsp[:,0],Node_Dsp[:,1]*1000, label=r'$u_1$')
sub4.plot(Node_Dsp[:,0],Node_Dsp[:,2]*1000, label=r'$u_2$')
sub4.legend()
Ymax,Xmax=Plot_max(Node_Dsp[:,0],Node_Dsp[:,1]*1000,r'$\delta_{max}$=',sub4,12)
sub4.set_xlim((0, tFinal))
sub4.set_ylim((-1.2*Ymax, 1.2*Ymax))

# ------------------------------
# GRAFICO 2
# ------------------------------

an = np.linspace(0, 2 * np.pi, 100)
rd=1.2*max([abs(Node_Dsp[:,1]).max(),abs(Node_Dsp[:,2]).max()])*1000
rf=1.2*max([abs(Elmt_Frc[:,1]).max(),abs(Elmt_Frc[:,2]).max()])/1000

fig2 = plt.figure()
sub4=plt.subplot(2, 4,(1,3))
sub4.set_xlabel('Longitudinal displacement'+r'$[mm]$',fontsize='small')
sub4.set_ylabel('Long. Force \n'+r'$[kN]$',fontsize='small')
sub4.tick_params(labelsize='small')
sub4.grid(True)
sub4.plot(Node_Dsp[:,1]*1000,Elmt_Frc[:,4]/1000)


sub5=plt.subplot(2, 4,4)
sub5.set_ylabel('Long. displacement'+r'$[mm]$',fontsize='small')
sub5.set_xlabel('Lateral displacement'+r'$[mm]$',fontsize='small')
sub5.tick_params(labelsize='small')
sub5.grid(True)
sub5.axis('equal')
sub5.plot(rd * np.cos(an), rd * np.sin(an),c='gray')
sub5.plot(Node_Dsp[:,1]*1000,Node_Dsp[:,2]*1000)


sub6=plt.subplot(2, 4,(5,7))
sub6.set_xlabel('Lateral displacement'+r'$[mm]$',fontsize='small')
sub6.set_ylabel('Lateral Force \n'+r'$[kN]$',fontsize='small')
sub6.tick_params(labelsize='small')
sub6.grid(True)
sub6.plot(Node_Dsp[:,2]*1000,Elmt_Frc[:,5]/1000)

sub7=plt.subplot(2, 4,8)
sub7.set_ylabel('Long. force'+r'$[kN]$',fontsize='small')
sub7.set_xlabel('Lateral force'+r'$[kN]$',fontsize='small')
sub7.tick_params(labelsize='small')
sub7.grid(True)
sub7.axis('equal')
sub7.plot(rf * np.cos(an), rf * np.sin(an))
sub7.plot(Elmt_Frc[:,1]/1000,Elmt_Frc[:,2]/1000)







