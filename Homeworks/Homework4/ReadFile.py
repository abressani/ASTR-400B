# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
import astropy.units as u

def Read(filename):
    '''
    Get file info
    Inputs: 
        Filename: name of the file that is being read
    Outputs:
        Time: time that the data in the file is representing
        Particles: number of particles in consideration in the file
        Data: #type, mass, x, y, z, vx, vy, vz
            type 1 is dark matter, type 2 is disk stars, type 3 is bulge stars
            mass: in units of 10^10 solar masses
            x, y, z: location from center of mass of galaxy in units of kpc 
            vx, vy, vz: velocity of particles centered milky way center in units of km/s
    '''
    #For time
    file = open(filename,'r')
    line1 = file.readline()
    label, value = line1.split()
    time = float(value)*u.Myr
    #For number of particles
    line2 = file.readline()
    label, value = line2.split()
    particles = float(value)
   
    file.close()
    
    #Storing remaining data of file
    data = np.genfromtxt(filename, dtype=None, names=True, skip_header=3)
    return time, particles, data