'''
Dark Photon in Medium Effects
(Calculates absorption of dark photons in
 various materials)
'''
import matplotlib
from matplotlib import pyplot as plt
import numpy as np
from numpy import pi as pi
from numpy import sqrt, cos, sin, exp
from numpy import random as rand
from matplotlib import ticker, cm
import math
import scipy as sp
from scipy import special as spc

element_list = ["Xenon","Argon","Helium","Sulfur Hexafluoride"]

#element = "Xenon" #Must be an element of element list
#density = 0.005 #g/cm^3

def f1f2Value(element,energy):
    '''
    Determine the atomic scattering factor
    given an element and an energy 

    Parameters
    ----------
    element (string) : Must be part of element list
    energy (float): Energy of dark photon in MeV 

    Returns
    -------
    f1: Real part of atomic scattering factor
    f2: Imaginary part of atomis scattering factor
    '''
    
    
    if element != "Sulfur Hexafluoride":
        f1Array = np.array([])
        f2Array = np.array([])
        EnergyArray = np.array([]) #MeV
        
        filename = element +"OpticalProperties.csv"
        file = open(filename)
        for line in file:
            line = line.split(',')
            EnergyArray = np.append(EnergyArray,float(line[0])/1e6)
            f1Array = np.append(f1Array,float(line[1]))
            f2Array = np.append(f2Array,float(line[2]))
        file.close()
        
        f1val = np.interp(energy, EnergyArray, f1Array)
        f2val = np.interp(energy, EnergyArray, f2Array)
            
    else:
        f1Array1 = np.array([])
        f2Array1 = np.array([])
        EnergyArray1 = np.array([]) #MeV
    
        f1Array2 = np.array([])
        f2Array2 = np.array([])
        EnergyArray2 = np.array([]) #MeV
        
        filename1 = "SulfurOpticalProperties.csv"
        filename2 = "FluorineOpticalProperties.csv"
        file1 = open(filename1)
        file2 = open(filename2)
        
        for line in file1:
            line = line.split(',')
            EnergyArray1 = np.append(EnergyArray1,float(line[0])/1e6)
            f1Array1 = np.append(f1Array1,float(line[1]))
            f2Array1 = np.append(f2Array1,float(line[2]))
        for line in file2:
            line = line.split(',')
            EnergyArray2 = np.append(EnergyArray2,float(line[0])/1e6)
            f1Array2 = np.append(f1Array2,float(line[1]))
            f2Array2 = np.append(f2Array2,float(line[2]))
        
        file1.close()
        file2.close()
        f1val1 = np.interp(energy,EnergyArray1,f1Array1)
        f1val2 = np.interp(energy,EnergyArray2,f1Array2)
        f2val1 = np.interp(energy,EnergyArray1,f2Array1)
        f2val2 = np.interp(energy,EnergyArray2,f2Array2)
        
        f1val = f1val1 + 6*f1val2
        f2val = f2val1 + 6*f2val2
        
    
    
    return(f1val,f2val)

def RefIndex(element, density, energy):
    '''
    Determine the refractive index of
    a given element at a specific energy
    
    Parameters
    element (string) : element of the detector
    density (float): density of target mass g/cm^3
    energy (float): energy of dark photon in MeV

    Returns
    nrefre = real part of refractive index
    nrefim = imaginary part of refractive index
    '''
    r0 = 2.82*1e-15 #m
    hc = 4.14e-21 * 3e8 #MeV m
    
    atomic_masses = {"Xenon": 131.2, "Argon":39.94,
                         "Helium":4, "Sulfur Hexafluoride": 146.06}#g/mol
    
    element_mass = atomic_masses[element]
    num_dens = (density/element_mass) * 6.02e23 * (100)**3 #m^{-3}
    
    f1val,f2val = f1f2Value(element,energy)
    
    nrefre = 1 - r0/(2 *pi) * (hc/energy)**2 * num_dens * f1val
    nrefim = -r0/(2 *pi) * (hc/energy)**2 * num_dens * f2val
    
    if element == "Helium" and energy > 0.004:
        filename = "HeliumIndexRefraction1kgm3.csv"
        EnergyArray = np.array([])
        deltaArray = np.array([])
        betaArray = np.array([])
        
        file = open(filename)
        for line in file:
            line = line.split(',')
            EnergyArray = np.append(EnergyArray,float(line[0])*1e-6)
            deltaArray = np.append(deltaArray, float(line[1]))
            betaArray = np.append(betaArray, float(line[2]))
            
        file.close()
        
        nrefre = 1 - (density*1e3)*np.interp(energy, EnergyArray, deltaArray)
        nrefim = -1* (density*1e3)* np.interp(energy, EnergyArray, betaArray)
    
    return(nrefre, nrefim)

def PiT(element,density,energy):
    '''
    Determine the polarization tensor of
    a given element at a specific energy
    
    Parameters
    element (string) : element of the detector
    density (float): density of target mass g/cm^3
    energy (float): energy of dark photon in MeV

    Returns
    PiTre = real part of polarization tensor (MeV)
    PiTim = imaginary part of polarization tensor (MeV)
    '''
    nrefre,nrefim = RefIndex(element,density,energy)
    nref2re = nrefre**2 - nrefim**2
    nref2im = 2 * nrefre * nrefim
    
    PiTre = energy**2 * (1 - nref2re)
    PiTim = energy**2 * -nref2im
    
    return(PiTre, PiTim)

def Eps2Scaling(element,density,energy,mA):
    '''
    Determine the polarization tensor of
    a given element at a specific energy
    
    Parameters
    element (string) : element of the detector
    density (float): density of target mass g/cm^3
    energy (float): energy of dark photon in MeV
    mA (float): mass of dark photon in MeV
    
    Returns
    Scale: Scaling of the effective mixing squared
    '''
    PiTre,PiTim = PiT(element, density, energy)
    
    Scale = mA**4 / ((mA**2 - PiTre)**2 + PiTim**2)
    return(Scale)
