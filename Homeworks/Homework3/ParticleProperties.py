#!/usr/bin/env python
# coding: utf-8

# In[16]:


# HW2 Part 3 Particle Properties Program

import numpy as np          #imports NumPy module and renames it np in this script for ease of use
import astropy.units as u   #imports AstroPy module and renames it u in this script for ease of use\
from ReadFile import Read   #imports Read from ReadFile

def ParticleInfo(filename, parttype, partnum): #defining new function that will read in 
                                               #a file name, the particle type, and the particle number
    time, totpart, data = Read(filename)       #gets all the info from the file
    
    #makes an index of all particles with the target type
    index = np.where(data['type'] == parttype) 
    
    #storing all the different fields of the target type
    mindex = data['m'][index]   #storing the masses
    xindex = data['x'][index]   #storing x-components of position
    yindex = data['y'][index]   #storing y-components of position
    zindex = data['z'][index]   #storing z-components of position
    vxindex = data['vx'][index] #storing x-components of velocity
    vyindex = data['vy'][index] #storing y-components of velocity
    vzindex = data['vz'][index] #storing z-components of velocity
    
    #DISTANCE
    magdist = np.sqrt(xindex[partnum-1]**2 + yindex[partnum-1]**2 + zindex[partnum-1]**2)
    #the magnitude of a vector is the square root of the sum of its squared components
    #i copied the syntax from the print statement in ReadFile
    
    magdist = np.around(magdist,3)*u.kpc  #rounds the magnitude distance to 3 decimal places
                                          #and putting units
    #VELOCITY
    magv = np.sqrt(vxindex[partnum-1]**2 + vyindex[partnum-1]**2 + vzindex[partnum-1]**2)
    #same math as for the distance
    
    magv = np.around(magv,3)*u.km/(1.0*u.s) #rounds the magnitude velocity to 3 decimal places
                                            #and putting the units in there
    #MASS
    solmass= mindex[partnum-1]*(10**10)*u.M_sun #converting from 10^10 Msun to Msun 
                                              #and putting the units in there
    return magdist, magv, solmass
    #returns the magnitude of the distance, the magnitude of velocity, 
    #and the mass in solar masses


# In[1]:


# testing code
#magdist, magv, solmass = ParticleInfo("../../MW_000.txt", 2, 99) 
# ParticleInfo gets the info from MW_000.txt and returns the 
# magdist, magv, and solmass, which are pointed back to
# new variables in this cell that share the same names

#print(magdist, magv, solmass) #printinggggggg

#lydist = np.around(magdist.to(u.lyr),3) #magdist is converted to ly, rounded to 3 decimal places
                                        #and saved to the new variable lightyear distance 

#print(lydist, magv, solmass) #printinggggggg again


# In[ ]:


#Outputs for line 500104 in MW_000.txt:
#4.245 kpc 312.135 km / s 1000000.0 solMass
#13845.338 lyr 312.135 km / s 1000000.0 solMass
#I did it! :D


# In[ ]:


# 02/02/2020 edit: changed 'partnum' to 'partnum-1' to match Prof. B code
# also so I don't have to remember to subtract one from particle(s) of interest

