#!/usr/bin/env python
# coding: utf-8

# In[1]:


# Homework 3
# Feb 6 2020
# S. Mackie

import numpy as np               # imports NumPy module and renames it np in this script for ease of use
import astropy.units as u        # imports AstroPy module and renames it u in this script for ease of use
from ReadFile import Read        # imports Read from ReadFile
#from astropy.table import QTable # imports the Qtable fucntion from astropy.table


# In[2]:


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


# In[3]:


# testing the function to make sure it works
#Mtot = ComponentMass("../../M31_000.txt", 3)
#print(Mtot)
#print(Mtot*10**12)


# In[4]:


# This cell is for making the table

# This will be the column containing the names of our galaxies and the Local Group
galname = ["Milky Way", "Andromeda (M31)", "Triangulum (M33)", "Local Group"]

# Computing the masses of the halos of the Milky Way, Andromeda, and Triangulum
# LGhalo is the total halo mass of the Local Group
MWhalo = ComponentMass("../../MW_000.txt", 1)
M31halo = ComponentMass("../../M31_000.txt", 1)
M33halo = ComponentMass("../../M33_000.txt", 1)
LGhalo = MWhalo + M31halo + M33halo
# This will be the column containing the masses of the halos
halomass = [MWhalo, M31halo, M33halo, LGhalo]

# Computing the masses of the disks of the Milky Way, Andromeda, and Triangulum
# LGdisk is the total disk mass of the Local Group
MWdisk = ComponentMass("../../MW_000.txt", 2)
M31disk = ComponentMass("../../M31_000.txt", 2)
M33disk = ComponentMass("../../M33_000.txt", 2)
LGdisk = np.around(MWdisk + M31disk + M33disk,3)
# This will be the column containing the masses of the disks
diskmass = [MWdisk, M31disk, M33disk, LGdisk]

# Computing the masses of the disks of the Milky Way and Andromeda
# LGbulge is the total bulge mass of the Local Group
MWbulge = ComponentMass("../../MW_000.txt", 3)
M31bulge = ComponentMass("../../M31_000.txt", 3)
LGbulge = np.around(MWbulge + M31bulge,3)
# This will be the column containing the masses of the bulges of the Milky Way and Andromeda
bulgemass = [MWbulge, M31bulge, 0.0, LGbulge]

# Computing the total masses of the Milky Way, Andromeda, Triangulum, and the Local Group
# LGtot is the total mass of the Local Group
MWtot = MWhalo + MWdisk + MWbulge
M31tot = M31halo + M31disk + M31bulge
M33tot = M33halo + M33disk
LGtot = MWtot + M31tot + M33tot
# This will be the column containing the total masses of the Milky Way, Andromeda, Triangulum, and the Local Group
Mtotgal = [MWtot, M31tot, M33tot, LGtot]

# Computing the baryon fractions of the Milky Way, Andromeda, Triangulum, and the Local Group
MWfbar = np.around((MWdisk + MWbulge) / MWtot,3)
M31fbar = np.around((M31disk + M31bulge) / M31tot,3)
M33fbar = np.around(M33disk / M33tot,3)
LGfbar = np.around((LGdisk + LGbulge) / LGtot,3)
# This will be the column containing the baryon fractions of the Milky Way, Andromeda, Triangulum, and the Local Group
fbar = [MWfbar, M31fbar, M33fbar, LGfbar]

# Printing the table so I can copy the values into LaTeX
#for i in range(len(galname)):
#    print(galname[i], halomass[i], diskmass[i], bulgemass[i], Mtotgal[i], fbar[i])

# I scrapped the below code, but I'll leave it in to look back at later
# assembling the table using Qtable from astropy.table
#LGtable = QTable([galname, halomass, diskmass, bulgemass, Mtotgal, fbar],
#                names = ('Galaxy Name', 'Halo Mass (10^12 solMass)', 'Disk Mass (10^12 solMass)',
#                        'Bulge Mass (10^12 solMass)', 'Total Mass (10^12 solMass)', 'Baryon Fraction'),
#                meta = {'name': 'Table of Sammie'})
                
# This will show the table
# LGtable


# In[ ]:




