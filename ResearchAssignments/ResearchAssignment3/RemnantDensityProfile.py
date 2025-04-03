# %%
# Import Modules
import numpy as np
from astropy import units as u
from astropy import constants as const

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

from ReadFile import Read
from CenterOfMass2 import CenterOfMass

# %% [markdown]
# # Plan for Research Project:
# 
# # Part 1: Find the halo denisty profile of the merger remnant.
# Note: Using snap 801 for MW and M31
# 
# Steps:
# 1) Extract halo information of MW and M31 and add them together
# 2) Sum mass enclosed in each shell
# 3) Divide mass by spherical volume of the shell
# 4) Store the density value in an array
# 5) Plot that density as a function of radius of the shell

# %%
# This is no longer needed
"""
def CombineFilesDM(file1, file2):
    '''
    Combine DM data from two files together and make an array out of them
    Inputs:
        file1: string
            Name of file 1
        file2: string
            Name of file 2
    Outputs:
        DMMerger: numpy array
            Array of DM data of merged galaxies/files
    '''
    # Extract info for MW and M31
    MWt, MWp, MWd = Read("MW_801.txt")
    M31t, M31p, M31d = Read("M31_801.txt")

    # Index to halo particles only
    MWindex = np.where(MWd['type'] == 1)
    M31index = np.where(M31d['type'] == 1)

    # Get halo particle info
    MW_DM = MWd[MWindex]
    M31_DM = M31d[M31index]

    # Add halo particles together from galaxies
    '''
    Indexes:
        0: type
        1: mass
        2-4: x,y,z
        5-7: vx, vy, vz
    '''    
    Merger = np.append(MW_DM, M31_DM)

    return Merger

# Get details about Merger remnant
Merger = CombineFilesDM("MW_801.txt", "M31_801.txt")
"""

# %%
# Get file info
filename = "Merger_801.txt"
Remnant_t, Remnant_p, Remnant_d = Read(filename)

# %%
# get COM of remnant
Remnant = CenterOfMass(filename, 1)

# Error tolerance
delta = 0.01
# Decrease volume factor
voldec = 2.0

# This is the center of mass of the remnant
Remnant_COMp = Remnant.COM_P(delta, voldec)

# %%
def Integrator(rmax, dr):
    '''
    Make spheres and sum masses inside
    Inputs:
        rmax: float
            Maximum radius to be considered in kpc
        dr: float
            Added radius for each shell
    Outputs:
        radii: numpy array
            Radius of sphere
        masses: numpy array
            Total mass in each sphere
    '''
    # Determine particle distance from COM
    xnew = abs(Remnant_COMp[0].value - Remnant_d['x'])
    ynew = abs(Remnant_COMp[1].value - Remnant_d['y'])
    znew = abs(Remnant_COMp[2].value - Remnant_d['z'])

    # Find magnitude of particle distance
    rnew = np.sqrt(xnew**2 + ynew**2 + znew**2)
    
    # Initial radius (if zero return nan for density)
    r0 = 0.1
    # Initialize radius and mass array
    radii = np.array([])
    masses = np.array([])

    # Loop over shells 
    while r0 < rmax:
        # Add radius to radii array
        radii = np.append(radii, r0)

        # Create index where particles are within sphere
        index = np.where(rnew < r0)

        # Take masses of the particles in sphere
        tot_mass = np.sum(Remnant_d['m'][index])
        if tot_mass == 'nan':
            tot_mass == 0
        # Add new mass value to array
        masses = np.append(masses, tot_mass)

        # Create new shell
        r0 += dr

    return radii, masses

# %%
def density(rmax, dr):
    ''' 
    Convert mass values into spherical density values 
    Inputs:
        rmax: float
            Maximum radius to be considered
        dr: float
            Radius step size
    Outputs:
        radii: numpy array
            Array of radii
        densities: numpy array
            Array of densities at radius values
    '''
    # Integrate masses with respect to radius
    radii, masses = Integrator(rmax, dr)

    # Convert radii to sphere volume
    sphere_vol = (4/3) * np.pi * radii**3

    # Convert masses and volumes to densities
    densities = (masses)/sphere_vol

    return radii, densities

# %%
Halo_density = density(3000, 5)

# %%
# Initialize Plots
fig, ax = plt.subplots()
# Plotting particle masses at given radius
ax.plot(Halo_density[0], Halo_density[1], 'b', linewidth = 2, label = 'Dark Matter')

# Formatting
ax.set_yscale('log')
ax.grid()
ax.set(title='Merger Remnant Density Profile', xlabel='Radius (kpc)', ylabel='Density (Msun*10e10/kpc^3)')
ax.legend(loc='best')


