#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import numpy as np
from ReadFile import Read
import astropy.units as u
from ReadFile import Read

def ComponentMass(filename, particle_type):
    '''
    Inputs:
        filename: name of the file to extract information from
        particle_type: user input for which type of particle they want to extract info from
            Type 1 is dark matter halo
            Type 2 is disk stars
            Type 3 is bulge stars
    Output:
        total_mass: (astropy unit) total mass of desired galaxy component
    '''
    #Read and store information from 'filename'
    time, particles, data = Read(filename)   
    #Get data for particles of certain type
    index = np.where(data['type'] == particle_type)
    
    #Get masses of particles of that type
    newm = data['m'][index]
    
    #Add up total mass of particles of that type
    total_mass = np.sum(newm)
    
    #Change units from 10^10 to 10^12 solar masses
    total_mass = total_mass/10**2 * u.T * u.Msun
    
    #Rounding
    total_mass = np.around(total_mass, 3)
    
    return total_mass

#test
#print(ComponentMass('MW_000.txt', 1))

stellarmass = .009
dm= .187
print(stellarmass/(stellarmass+dm))