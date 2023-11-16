'''
Rate of dark photon absorption in a specific experiment
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

import Dark_Matter_Distribution_and_J_Factors as DMMod
import Dark_Photon_in_Medium_Effects as DPMod

f0 = 0.1
mchi = 1 #keV
mA = 0.29 #eV
sigmav = DMMod.SigmavOpt(f0,mchi*1e-6)
epsilon = 2e-13

Volume = 10 #m^3
element = "Helium"
density = 0.0002 #g cm^{-3}

mass = Volume * density * (100)**3 * (1e-6) #ton

MeV_to_inv_cm = 5.06e10
mchiGeV = mchi * 1e-6
mchiMeV = mchi * 1e-3
mAMeV = mA * 1e-6


#Determine Rate for Single Set of parameters
'''
Phi = DMMod.FluxDiff(f0, sigmav, mchiGeV) #cm^{-2} s^{-1}

eps2eff = epsilon**2 * DPMod.Eps2Scaling(element,density,mchiMeV,mAMeV)

PiTre,PiTim = DPMod.PiT(element, density, mchiMeV) # MeV^2

Rate = eps2eff * Volume * Phi * (PiTim/mchiMeV) * (100)**3 * MeV_to_inv_cm #s^{-1}

print(Rate/Volume, "s^{-1} m^{-3}")
print(Rate*3.1e7/mass, "yr^{-1} ton^{-1}")
'''

#Determine Exclusions Given a Desired Rate
'''
Rate_Des = 2 *mass /(3.1e7) #2 events per ton per year

mA_vals = np.logspace(-2,2,300) #eV
eps_true_vals = np.array([])

for mA in mA_vals:
    print('mA', mA)
    mAMeV = mA * 1e-6
    
    Phi = DMMod.FluxDiff(f0, sigmav, mchiGeV) #cm^{-2} s^{-1}

    eps2eff = epsilon**2 * DPMod.Eps2Scaling(element,density,mchiMeV,mAMeV)

    PiTre,PiTim = DPMod.PiT(element, density, mchiMeV) # MeV^2

    Rate = eps2eff * Volume * Phi * (PiTim/mchiMeV) * (100)**3 * MeV_to_inv_cm #s^{-1}
    
    eps_des = np.sqrt(Rate_Des/Rate) * epsilon
    
    eps_true_vals = np.append(eps_true_vals,eps_des)
    
fig = plt.figure()
plt.plot(mA_vals, eps_true_vals)
plt.xscale('log')
plt.yscale('log')
'''
