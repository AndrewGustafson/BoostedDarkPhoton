# BoostedDarkPhoton
Code to determine the rate of dark photon absorptions in different detector materials.

## Dark_Matter_Distribution_and_J_Factors.py
Code that determines the present-day distribution of dark matter given the initial fraction, mass, and annihilation cross section. It also takes these same values and computes the flux of dark photons at Earth.

## Dark_Photon_in_Medium_Effects.py
Code that determines the effective mixing parameter of dark photons in different media. It reads in optical properties of different materials via included data files.

## Dark_Photon_Rate.py
Determines the rate of dark photon absorption events in a detector. Inputs are the material and size of the detector, along with the parameters of the theory. Included (but commented out) is lines of code that produces a curve in kinetic mixing vs dark photon mass parameter space that corresponds to a desired rate of events

## Data files
The following files contain the atomic scattering factors as a function of energy.

ArgonOpticalProperties.csv

XenonOpticalProperties.csv

HeliumOpticalProperties.csv

SulfurOpticalProperties.csv

FluorineOpticalProperties.csv

The following file has the values of 1 - refractive index for helium with a density of 1 kilogram per cubic meter

HeliumIndexRefraction1kgm3.csv
