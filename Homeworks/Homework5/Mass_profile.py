#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 18 18:54:07 2025

@author: alessandro
"""

import numpy as np
import matplotlib.pyplot as plt
import astropy.units as u
from astropy.constants import G
from ReadFile import Read
from CenterOfMass import CenterOfMass

class MassProfile:
    
    def __init__(self, galaxy, snap):
        '''
        Reconstruct the file to be read
        
        PARAMETERS
        ----------
        galaxy: string with galaxy name
        snap: snapshot number (e.g. 0, 1, etc.)
        '''
        # Add a string of the filenumber to the value “000”
        ilbl = '000' + str(snap)
        # Remove all but the last 3 digits
        ilbl = ilbl[-3:]
        # Name the file
        self.filename = "%s_"%(galaxy) + ilbl + '.txt'
        
        # Get file info
        self.time, self.total, self.data = Read(self.filename)
        # Get x,y,z positions of the particles
        self.x = self.data['x'] * u.kpc
        self.y = self.data['y'] * u.kpc
        self.z = self.data['z'] * u.kpc
        # Get mass of the particles
        self.m = self.data['m']
        
        self.gname = galaxy
        
    def MassEnclosed(self, ptype, radii):
        '''
        Calculate the mass enclosed for each particle type within a galaxy at each radius
        Inputs:
            ptype: integer
                Type of particle user is interested in
                    Type 1 is dark matter
                    Type 2 is disk stars
                    Type 3 is bulge stars
            radii: array
                Array of radii for which the mass enclosed will be calculated at
        Outputs:
            masses: array
                Mass enclosed for each radius in inputed array
        '''
        # Finding center of mass of galaxy using disk stars
        COM = CenterOfMass(self.filename, 2)
        COM_p = COM.COM_P(.1)
        
        # Indexing particles of interest
        index = np.where(self.data['type'] == ptype)
        # Obtain new x, y, z
        x_part = self.data['x'][index] * u.kpc
        y_part = self.data['y'][index] * u.kpc
        z_part = self.data['z'][index] * u.kpc
        mass_part = self.data['m'][index]
        
        # Find relative distance from COM
        part_dist = np.sqrt((x_part - COM_p[0])**2 + (y_part - COM_p[1])**2 + (z_part - COM_p[2])**2)
        
        # Initializing mass array
        masses = np.array([])
        radii *= u.kpc
        
        # Loop over larger and larger radii
        for radius in radii:
            # Get particle info enclosed in radius
            part_enclosed = np.where(part_dist <= radius)
            mass = mass_part[part_enclosed]
            # Add to masses array
            masses = np.append(masses, np.sum(mass))
        # Units, units, units
        masses = masses * 1e10 * u.Msun
        return masses
    
    
    def MassEnclosedTotal(self, radii):
        '''
        Get total mass of all particles within each radius from center of mass
        Inputs:
            radii: array
                user-inputed array of radii where total mass will be calculated within
        Outputs:
            masses: array
                Total enclosed mass of all particles within each radius
        '''
        # Find masses of particles for each type for radii array
        dm_mass = self.MassEnclosed(1, radii)
        disk_mass = self.MassEnclosed(2, radii)
        # Ignore bulge for M33
        if self.gname != "M33":
            bulge_mass = self.MassEnclosed(3, radii)
            total_masses = dm_mass + disk_mass + bulge_mass
        else: 
            total_masses = dm_mass + disk_mass
        
        return total_masses
    
    
    def HernquistMass(self, r, m_halo, h_a):
        '''
        Function that defines the Hernquist 1990 mass profile 
        Inputs:
            r: array
                distance from COM in kpc
            a: astropy quantity
                scale radius of the Hernquist profile in kpc
            m: float
                total mass in units of 1e12 Msun 
            
        Ouputs:
            mass:  astropy quantity
                total mass within the input radius r in Msun
        '''
        # Constructing Equation
        a = m_halo * u.Msun
        b = r**2/(h_a + r)**2
        # Hernquist Profile and units
        masses = a*b
        
        return masses
    
    def CircularVelocity(self, ptype, radii):
        '''
        Calculate circular speed of particles in km/s
        Inputs:
            ptype: integer
                Type of particle user is interested in
                    Type 1 is dark matter
                    Type 2 is disk stars
                    Type 3 is bulge stars
            radii: array
                Array of radii for which the velocity will be calculated at
        Outputs:
            circular_speeds: array
                circular speed of particles at each radius
        '''
        # Convert G units
        G.to(u.kpc*u.km**2/u.s**2/u.Msun)
        # Find masses of particles
        M = self.MassEnclosed(ptype, radii)
        # Calculating circular velocity 
        circular_speeds = np.sqrt(G*M/radii)
        
        return np.round(circular_speeds, 2)
    
    def CircularVelocityTotal(self, radii):
        '''
        ...
        Input: 
            radii: array
                Array of radii for which the velocity will be calculated at
        Output:
            v_circ_total: array
                Circular speed of all particles at radius
        '''
        # Convert G units
        G.to(u.kpc*u.km**2/u.s**2/u.Msun)
        # Find masses of particles
        M = self.MassEnclosedTotal(radii)
        # Calculating circular velocity 
        circular_speeds = np.sqrt(G*M/radii)
        
        return np.round(circular_speeds, 2)
    
    def HernquistVCirc(self, r, m_halo, h_a):
        '''
        Function that defines the Hernquist 1990 mass profile 
        Inputs:
            r: array
                distance from COM in kpc
            a: astropy quantity
                scale radius of the Hernquist profile in kpc
            m_halo: float
                total halo mass in units of 1e12 Msun 
            
        Ouputs:
            hern_speed: astropy quantity
                Speed of particles according to Hernquist mass profile
        '''
        # Convert G units
        G.to(u.kpc*u.km**2/u.s**2/u.Msun)
        # Get Hernquist Mass
        M = self.HernquistMass(r, m_halo, h_a)
        # Apply mass to velocity equation
        hern_speed = np.sqrt(G*M/r)
        return np.round(hern_speed, 2)
        
# Array of radii
radii = np.linspace(0.1, 30, 100)
# Get galaxy info
MW = MassProfile('MW', 000)
M31 = MassProfile('M31', 000)
M33 = MassProfile('M33', 000)
# Get halo masses
# MW halo
MW_halo_index = np.where(MW.data['type'] == 1)
MW_halo_tot_mass = np.sum(MW.data['m'][MW_halo_index]) * 1e10
# M31 halo
M31_halo_index = np.where(M31.data['type'] == 1)
M31_halo_tot_mass = np.sum(M31.data['m'][M31_halo_index]) * 1e10
# M33 halo
M33_halo_index = np.where(M33.data['type'] == 1)
M33_halo_tot_mass = np.sum(M33.data['m'][M33_halo_index]) * 1e10

# Masses
'''
Milky Way
'''
# Get particle masses for MW
MW_halo = MW.MassEnclosed(1, radii)
MW_disk = MW.MassEnclosed(2, radii)
MW_bulge = MW.MassEnclosed(3, radii)
# Get total mass info
MW_total = MW.MassEnclosedTotal(radii)
# Get Hernquist info
MW_Hern = MW.HernquistMass(radii, MW_halo_tot_mass, 60)


# Initialize Plots
fig1, ax1 = plt.subplots()
# Plotting particle masses at given radius
ax1.plot(radii, MW_halo, 'b', linewidth = 2, label = 'Dark Matter')
ax1.plot(radii, MW_disk, color='red', linewidth = 2, label = 'Disk')
ax1.plot(radii, MW_bulge, 'g', linewidth = 2, label = 'Bulge')
# Plotting total mass at given radius
ax1.plot(radii, MW_total, color='black', linewidth = 2, label = 'Total Mass')
# Hernquist fit at given radius
ax1.plot(radii, MW_Hern, color='magenta', linewidth = 2, linestyle = 'dotted', label = 'Hernquist Profile (a=60 kpc)')

# Formatting
ax1.set_yscale('log')
ax1.grid()
ax1.set(title='Milky Way Mass Profiles', xlabel='Radius (kpc)', ylabel='Mass (Msun)')
ax1.legend(loc='best')


'''
M31
'''
# Get particle masses for M31
M31_halo = M31.MassEnclosed(1, radii)
M31_disk = M31.MassEnclosed(2, radii)
M31_bulge = M31.MassEnclosed(3, radii)
# Get total mass info
M31_total = M31.MassEnclosedTotal(radii)
# Get Hernquist info
M31_Hern = M31.HernquistMass(radii, M31_halo_tot_mass, 60)


# Initialize Plots
fig2, ax2 = plt.subplots()
# Plotting particle masses at given radius
ax2.plot(radii, M31_halo, 'b', linewidth = 2, label = 'Dark Matter')
ax2.plot(radii, M31_disk, color='red', linewidth = 2, label = 'Disk')
ax2.plot(radii, M31_bulge, 'g', linewidth = 2, label = 'Bulge')
# Plotting total mass at given radius
ax2.plot(radii, M31_total, color='black', linewidth = 2, label = 'Total Mass')
# Hernquist fit at given radius
ax2.plot(radii, M31_Hern, color='magenta', linewidth = 2, linestyle = 'dotted', label = 'Hernquist Profile (a=60 kpc)')

# Formatting
ax2.set_yscale('log')
ax2.grid()
ax2.set(title='M31 Mass Profiles', xlabel='Radius (kpc)', ylabel='Mass (Msun)')
ax2.legend(loc='best')

'''
M33
'''
# Get particle masses for M33
M33_halo = M33.MassEnclosed(1, radii)
M33_disk = M33.MassEnclosed(2, radii)
M33_bulge = M33.MassEnclosed(3, radii)
# Get total mass info
M33_total = M33.MassEnclosedTotal(radii)
# Get Hernquist info
M33_Hern = M33.HernquistMass(radii, M33_halo_tot_mass, 25)


# Initialize Plots
fig3, ax3 = plt.subplots()
# Plotting particle masses at given radius
ax3.plot(radii, M33_halo, 'b', linewidth = 2, label = 'Dark Matter')
ax3.plot(radii, M33_disk, color='red', linewidth = 2, label = 'Disk')
# Plotting total mass at given radius
ax3.plot(radii, M33_total, color='black', linewidth = 2, label = 'Total Mass')
# Hernquist fit at given radius
ax3.plot(radii, M33_Hern, color='magenta', linewidth = 2, linestyle = 'dotted', label = 'Hernquist Profile (a=25 kpc)')

# Formatting
ax3.set_yscale('log')
ax3.grid()
ax3.set(title='M33 Mass Profiles', xlabel='Radius (kpc)', ylabel='Mass (Msun)')
ax3.legend(loc='best')


# Velocities 
'''
Milky Way
'''
# Get particle velocities for MW
MW_halo = MW.CircularVelocity(1, radii)
MW_disk = MW.CircularVelocity(2, radii)
MW_bulge = MW.CircularVelocity(3, radii)
# Get total mass info
MW_total = MW.CircularVelocityTotal(radii)
# Get total mass of halo
MW_halo_mass = MW.MassEnclosed(1, radii)
# Get Hernquist info
MW_Hern = MW.HernquistVCirc(radii, MW_halo_tot_mass, 60)


# Initialize Plots
fig4, ax4 = plt.subplots()
# Plotting particle masses at given radius
ax4.plot(radii, MW_halo, 'b', linewidth = 2, label = 'Dark Matter Velocity')
ax4.plot(radii, MW_disk, color='magenta', linewidth = 2, label = 'Disk Velocity')
ax4.plot(radii, MW_bulge, 'g', linewidth = 2, label = 'Bulge Velocity')
# Plotting total mass at given radius
ax4.plot(radii, MW_total, color='black', linewidth = 2, label = 'Total Velocity')
# Hernquist fit at given radius
ax4.plot(radii, MW_Hern, color='red', linewidth = 2, linestyle = 'dotted', label = 'Hernquist V Profile (a=60 kpc)')

# Formatting
ax4.grid()
ax4.set(title='Milky Way Velocity Profiles', xlabel='Radius (kpc)', ylabel='Velocity (km/s)')
ax4.legend(loc='best')

'''
M31
'''
# Get particle velocities for M31
M31_halo = M31.CircularVelocity(1, radii)
M31_disk = M31.CircularVelocity(2, radii)
M31_bulge = M31.CircularVelocity(3, radii)
# Get total mass info
M31_total = M31.CircularVelocityTotal(radii)
# Get total mass of halo
M31_halo_mass = M31.MassEnclosed(1, radii)
# Get Hernquist info
M31_Hern = M31.HernquistVCirc(radii, M31_halo_tot_mass, 60)


# Initialize Plots
fig5, ax5 = plt.subplots()
# Plotting particle masses at given radius
ax5.plot(radii, M31_halo, 'b', linewidth = 2, label = 'Dark Matter Velocity')
ax5.plot(radii, M31_disk, color='magenta', linewidth = 2, label = 'Disk Velocity')
ax5.plot(radii, M31_bulge, 'g', linewidth = 2, label = 'Bulge Velocity')
# Plotting total mass at given radius
ax5.plot(radii, M31_total, color='black', linewidth = 2, label = 'Total Velocity')
# Hernquist fit at given radius
ax5.plot(radii, MW_Hern, color='red', linewidth = 2, linestyle = 'dotted', label = 'Hernquist V Profile (a=60 kpc)')

# Formatting
ax5.grid()
ax5.set(title='M31 Velocity Profiles', xlabel='Radius (kpc)', ylabel='Velocity (km/s)')
ax5.legend(loc='best')

'''
M33
'''
# Get particle velocities for M33
M33_halo = M33.CircularVelocity(1, radii)
M33_disk = M33.CircularVelocity(2, radii)
M33_bulge = M33.CircularVelocity(3, radii)
# Get total mass info
M33_total = M33.CircularVelocityTotal(radii)
# Get total mass of halo
M33_halo_mass = M33.MassEnclosed(1, radii)
# Get Hernquist info
M33_Hern = M33.HernquistVCirc(radii, M33_halo_tot_mass, 25)


# Initialize Plots
fig6, ax6 = plt.subplots()
# Plotting particle masses at given radius
ax6.plot(radii, M33_halo, 'b', linewidth = 2, label = 'Dark Matter Velocity')
ax6.plot(radii, M33_disk, color='magenta', linewidth = 2, label = 'Disk Velocity')
# Plotting total mass at given radius
ax6.plot(radii, M33_total, color='black', linewidth = 2, label = 'Total Velocity')
# Hernquist fit at given radius
ax6.plot(radii, M33_Hern, color='red', linewidth = 2, linestyle = 'dotted', label = 'Hernquist V Profile (a=25 kpc)')

# Formatting
ax6.grid()
ax6.set(title='M33 Velocity Profiles', xlabel='Radius (kpc)', ylabel='Velocity (km/s)')
ax6.legend(loc='best')

plt.show()
