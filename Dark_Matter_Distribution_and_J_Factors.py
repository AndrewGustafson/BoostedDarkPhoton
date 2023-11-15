'''
Dark Matter Distributions and J Factors
(Calculates Distribution of annihilating
 dark matter in galaxy and the corresponding
 modern day J-Factor)
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

kpc_to_cm = 3.084e21
cm_to_kpc = 1/kpc_to_cm

def rho_NFW(r):
    '''
    

    Parameters
    ----------
    r : Distance from Galactic Center (kpc)

    returns
    dens : density at r (GeV/cm^3)

    '''
    rho_s = 0.184 #GeV cm^-3
    r_s = 24.42 #kpc
    dens = rho_s *(r_s/r)* (1 + r/r_s)**(-2)
    return(dens)

r_vals = np.linspace(1e-5,60,3000) #kpc
dr = r_vals[1] - r_vals[0]

massMW = 4*pi*sum(rho_NFW(r_vals)*r_vals**2 * dr)*kpc_to_cm**3 #GeV
#print(mass)
rho2Integral = 4*pi*sum(rho_NFW(r_vals)**2 * r_vals**2 * dr)*kpc_to_cm**3 #GeV cm^-3
#print(rho2Integral)

def f_diff(f0,sigmav, mchi):
    '''
    Calculate the fraction of chi compared to the NFW density
        Used in diffusion limit
    Parameters
    ----------
    f0 : Initial fraction of NFW made by \chi
    sigmav : annihilation cross section (cm^3 s^-1)
    mchi: mass of \chi (GeV)

    Returns
    -------
    f: current fraction of \chi compared to NFW DM
    '''
    t0 = 13.6e9 * 3.1e7 #age of universe in seconds
    
    f = (1/f0 + (sigmav * t0 *rho2Integral)/(mchi*massMW))**(-1)
    return(f)

def rho_diff(r, f0,sigmav,mchi):
    '''
    Calculate the density of chi in the diffusion limit
    Parameters
    ----------
    r: distance from galactic center (kpc)
    f0 : Initial fraction of NFW made by \chi
    sigmav : annihilation cross section (cm^3 s^-1)
    mchi: mass of \chi (GeV)

    Returns
    -------
    dens: density of \chi at position r (GeV cm^-3)
    '''
    f = f_diff(f0,sigmav,mchi)
    rhoTot = rho_NFW(r) #GeV cm^-3
    dens = f*rhoTot #GeV cm^-3
    return(dens)

def rho_SIE(r,f0,sigmav,mchi):
    '''
    Calculate the density of chi in the spatially
        independent evolution
    Parameters
    ----------
    r: distance from galactic center (kpc)
    f0 : Initial fraction of NFW made by \chi
    sigmav : annihilation cross section (cm^3 s^-1)
    mchi: mass of \chi (GeV)

    Returns
    -------
    dens: density of \chi at position r (GeV cm^-3)
    '''
    t0 = 13.6e9 * 3.1e7 #age of universe in seconds
    numerator = f0 * rho_NFW(r) * mchi
    denominator = mchi + f0 * rho_NFW(r) * sigmav * t0
    dens = numerator/denominator
    return(dens)

def JFactorDiff(f0,sigmav,mchi, num = int(4e6)):
    '''
    Calculate the J Factor for chi in diffusion limit
    Parameters
    ----------
    f0 : Initial fraction of NFW made by \chi
    sigmav : annihilation cross section (cm^3 s^-1)
    mchi: mass of \chi (GeV)

    Returns
    -------
    J: J-Factor in GeV^2 cm^{-5}
    '''
    rsol = 8.3 #kpc
    xvals = 60*rand.random(num)#np.linspace(0,60,xnum) #kpc
    costhetavals = 1-2*rand.random(num)#np.linspace(-1,1,cosnum)
    xrange = 60
    costhetarange = 2
    diff_element = xrange * costhetarange / num
    
    
    rvals = np.sqrt(rsol**2 + xvals**2 - 2*rsol* xvals*costhetavals)
    J = 2*pi*diff_element*np.sum((rho_diff(rvals,f0,sigmav,mchi)**2))\
        *kpc_to_cm #GeV^2 cm^-5
    return(J)

def JFactorSIE(f0,sigmav,mchi,num = int(4e6)):
    '''
    Calculate the J Factor for chi with
        Spatially Independent Evolution
    Parameters
    ----------
    f0 : Initial fraction of NFW made by \chi
    sigmav : annihilation cross section (cm^3 s^-1)
    mchi: mass of \chi (GeV)

    Returns
    -------
    J: J-Factor in GeV^2 cm^{-5}
    '''
    rsol = 8.3 #kpc
    xvals = 60*rand.random(num)#np.linspace(0,60,xnum) #kpc
    costhetavals = 1-2*rand.random(num)#np.linspace(-1,1,cosnum)
    xrange = 60
    costhetarange = 2
    diff_element = xrange * costhetarange / num
    
    
    rvals = np.sqrt(rsol**2 + xvals**2 - 2*rsol* xvals*costhetavals)
    J = 2*pi*diff_element*np.sum((rho_SIE(rvals,f0,sigmav,mchi)**2))\
        *kpc_to_cm #GeV^2 cm^-5
    return(J)

def FluxDiff(f0,sigmav,mchi):
    '''
    Calculate the flux of dark photons in the diffusion limit

    f0 : Initial fraction of NFW made by \chi
    sigmav : annihilation cross section (cm^3 s^-1)
    mchi: mass of \chi (GeV)

    Returns
    -------
    Phi: Flux in cm^{-2} s^{-1}
    '''
    
    J = JFactorDiff(f0,sigmav,mchi) #GeV^{2} cm^{-5}
    Phi = (1/(4*pi)) * sigmav/(4*mchi**2) * (2/3) * J #cm^{-2} s^{-1}
    
    return(Phi)
    
def FluxSIE(f0,sigmav,mchi):
    '''
    Calculate the flux of dark photons in
        spatially independent evolution

    f0 : Initial fraction of NFW made by \chi
    sigmav : annihilation cross section (cm^3 s^-1)
    mchi: mass of \chi (GeV)

    Returns
    -------
    Phi: Flux in cm^{-2} s^{-1}
    '''
    
    J = JFactorSIE(f0,sigmav,mchi) #GeV^{2} cm^{-5}
    Phi = (1/(4*pi)) * sigmav/(4*mchi**2) * (2/3) * J #cm^{-2} s^{-1}
    
    return(Phi)

def SigmavOpt(f0,mchi):
    '''
    Calculate the annihilation cross section which produces
        the largest flux of dark photons

    f0 : Initial fraction of NFW made by \chi
    mchi: mass of \chi (GeV)

    Returns
    -------
    sigmav: cross section cm^{-3} s^{-1}
    '''
    t0 = 13.6e9 * 3.1e7 #age of universe in seconds
    sigmav = mchi/(t0 * f0) * (massMW/rho2Integral)
    return(sigmav)

frac = 0.1
mchi = 1e-6 #GeV
sigmav_vals = np.logspace(-30,-20,30)

'''
conv = False
while conv == False:
    print('num',num)
    newJ = JFactorSIE(1,0,1e-5,num)
    print('New J',newJ*cm_to_kpc)
    
    if abs(newJ-prevJ)/newJ < 0.005:
        conv = True
    else:
        num = num*2
        prevJ = newJ
'''    


'''
SIE_Flux_Vals = np.array([])
Diff_Flux_Vals = np.array([])
for sv in sigmav_vals:
    print('sigma v',sv)
    SIE_Flux_Vals = np.append(SIE_Flux_Vals,FluxSIE(frac,sv,mchi))
    Diff_Flux_Vals = np.append(Diff_Flux_Vals,FluxDiff(frac,sv,mchi))

fig = plt.figure()

plt.plot(sigmav_vals,SIE_Flux_Vals,label = "SIE")
plt.plot(sigmav_vals, Diff_Flux_Vals,label = "Diff limit")

plt.xscale('log')
plt.yscale('log')
plt.legend()
plt.xlabel('$\sigma$v [$cm^3 s^{-1}$]')
plt.ylabel('$|Phi$ [$cm^{-2} s^{-1}$]')
'''
