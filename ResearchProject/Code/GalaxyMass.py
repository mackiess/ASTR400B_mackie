# From Homework 3
# Feb 6 2020
# S. Mackie

import numpy as np               # imports NumPy module and renames it np in this script for ease of use
import astropy.units as u        # imports AstroPy module and renames it u in this script for ease of use
from ReadFile import Read        # imports Read from ReadFile
#from astropy.table import QTable # imports the Qtable fucntion from astropy.table

# Function to compute the total mass of any desired galaxy component
# total mass of the component = sum of the mass of all particles making up that component

def ComponentMass(filename, parttype):
    # Inputs:
    #    filename = name of file storing the data of the galaxy of interest: so far these are M?(?)_000.txt
    #    parttype = the type of particle making up the component of interest
    #      1 = Halo, 2 = Disk, 3 = Bulge
    # Returns:
    #    Mtot = the total mass of any desired galaxy component in 10^12 solar masses
    
    # receiving the time, total number of particles, and all the data from our target data file
    time, totpart, data = Read(filename)
    
    # makes an index of all particles with the target type
    index = np.where(data['type'] == parttype) 
    
    # storing the masses of the target type
    mindex = data['m'][index]   
    
    # sums all masses in the index of masses, puts the total mass into 10^12 instead of 10^10
    # and rounds to 3 decimal places 
    # actual u.M_Sun units aren't appended in this function because QTable hates that
    Mtot = np.around(np.sum(mindex)/100,3) * u.M_sun
    
    return Mtot 