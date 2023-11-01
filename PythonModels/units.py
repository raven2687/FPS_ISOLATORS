# -*- coding: utf-8 -*-
"""
Created on Sun Jan  8 17:13:08 2023
UNITS LIBRARY
@author: paulv
"""
class SI_m_units:
    
    # length
    m=1.0
    
    # force 
    N=1.00
        
    # mass
    kg=1.0
    
    # time
    sec=1.00
    
    # LENGTH UNITS TRANSFORMATIONS
    cm=m/100
    mm=m/1000
    Km=m*1000
    
    # gravity
    g = 9.81*m/sec**2
    
    # English Units
    inch=2.54*cm
    ft=12*inch
    
    # FORCE UNITS TRANSFORMATIONS
    kgf=g*N
    kN=1000*N
    Tonf=1000*kgf
    lbf=kgf/2.2
    kipf=1000*lbf
    
    # MASS TRANSFORMATIONS
    lb=kg/2.2
    kip=1000*lb
    Ton=1000*kg
    
    #Pressure
    psi=1*lbf/inch**2
    MPa=N/mm**2
    
class SI_mm_units:
    
    # length
    mm=1.0
    
    # force 
    N=1.00
        
    # mass
    Ton=1.0
    
    # time
    sec=1.00

    # LENGTH UNITS TRANSFORMATIONS
    m=mm*1000
    cm=m/100
    Km=m*1000
    
    # gravity
    g = 9.81*m/sec**2
    
    # English Units
    inch=2.54*cm
    ft=12*inch
    
    # FORCE UNITS TRANSFORMATIONS
    kgf=g*N
    
    Tonf=1000*kgf
    lbf=kgf/2.2
    kipf=1000*lbf
    
    # MASS TRANSFORMATIONS
    kg=Ton/1000
    lb=kg/2.2
    kip=1000*lb
    
    
    #Pressure
    psi=1*lbf/inch**2
    MPa=N/mm**2    

