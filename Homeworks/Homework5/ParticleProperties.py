#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import numpy as np
from ReadFile import Read
import astropy.units as u

def ParticleInfo(filename, particletype, particlenum):
    '''
    Extract mass, distance (magnitude), and velocity (magnitude) given particle type and number
    Inputs: 
        filename: name of the file that is being read
        particletype: user input for which type of particle they want to exact info of
            type 1 is dark matter, type 2 is disk stars, type 3 is bulge stars
        particlenum: user-chosen index number of the particle of that type
    Outputs:
        distance: magnitude of distance in kpc of particle from center of mass of milky way
        velocity: magnitude of velocity in km/s of particle with origin at center of mass of milky way
        mass: mass in solar units of particle user has selected
        ...
    '''
    #Get data from file using Read function
    time, particles, data = Read(filename)
    #Get data for particles of certain type
    index = np.where(data['type'] == particletype)
    
    #Get masses, (x,y,z) positions, and (x,y,z) velocities of particles of that type
    newm = data['m'][index]
    
    newx = data['x'][index]
    newy = data['y'][index]
    newz = data['z'][index]
    
    newvx = data['vx'][index]
    newvy = data['vy'][index]
    newvz = data['vz'][index]
    
    #Finding magnitude of displacement and velocity
    distance = np.sqrt(newx**2 + newy**2 + newz**2)
    velocity = np.sqrt(newvx**2 + newvy**2 + newvz**2)
    
    #Units!!
    mass = newm*(10**10) * u.Msun
    distance = distance * u.kpc
    velocity = velocity * (u.km/u.s)
    
    #Rounding
    distance = np.around(distance, 3)
    velocity = np.around(velocity, 3)
    mass = np.around(mass, 3)
    
    #Return distance, velocity, and mass of user-chosen particle
    return distance[particlenum], velocity[particlenum], mass[particlenum]

'''
#Getting info of 100th particle of type 2
test = (ParticleInfo('MW_000.txt', 2, 99))
print(test)

#Converting to lightyears and rounding to 3 decimals
distance = test[0].to(u.lyr)
distance = np.around(distance, 3)
print(distance)
'''