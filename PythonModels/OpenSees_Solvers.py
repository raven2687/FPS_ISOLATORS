# -*- coding: utf-8 -*-
"""
Created on Thu Nov  2 03:16:07 2023

@author: paulv
"""

##############################################################################
# IMPORT LIBRARIES
##############################################################################
import numpy as np
import openseespy.opensees as ops
import matplotlib.pyplot as plt
import os
from units import SI_m_units as U
from scipy.io import loadmat


def Solve_System(dt,dtAna,dtMin,dtMax,ok,tCurrent,tFinal,Tol,timeu2,u2):
    test = {1:'NormDispIncr', 2: 'RelativeEnergyIncr', 4: 'RelativeNormUnbalance',5: 'RelativeNormDispIncr', 6: 'NormUnbalance'}
    algorithm = {1:'KrylovNewton', 2: 'SecantNewton' , 4: 'RaphsonNewton',5: 'PeriodicNewton', 6: 'BFGS', 7: 'Broyden', 8: 'NewtonLineSearch'}
    while ok==0 and tCurrent<tFinal:    
        ok=ops.analyze(1,dtAna)
        if ok!=0:
            for i in test:    
                for j in algorithm:
                    if j<4:
                        ops.algorithm(algorithm[j], '-initial')
                    else:
                        ops.algorithm(algorithm[j])
                    print('test: ',test[i],' Algorithm: ',algorithm[j],' ok=',ok)
                    ops.test(test[i], Tol, 1000)
                    while ok!=0 and dtAna/2>=dtMin:
                        ok = ops.analyze(1,dtAna)
                        if ok!=0:
                            dtAna=dtAna/2
                            print('REDUCING time step size (dtNew=',dtAna,')')
                        else:
                            tCurrent = ops.getTime()
                            ops.algorithm('Newton')
                            print('t=',tCurrent,'sec')
                            timeu2.append(tCurrent)
                            u2.append(ops.nodeDisp(2,1))
                            if dtAna*2<=dtMax:
                                dtAna=dtAna*2
                                print('INCREASING time step size (dtNew=',dtAna,')')
                    if ok!=0:
                        print('Change Algorithm')
                        dtAna=dt
                        
        else:
            tCurrent = ops.getTime()
            print('t=',tCurrent,'sec')
            timeu2.append(tCurrent)
            u2.append(ops.nodeDisp(2,1))
            if dtAna*2<=dtMax:
                dtAna=dtAna*2
                print('INCREASING time step size (dtNew=',dtAna,')')
       
    if ok!=0:
        print(f'Model Failed (time={tCurrent})')
    else:
        print('Response-history analysis completed')
    return u2,timeu2