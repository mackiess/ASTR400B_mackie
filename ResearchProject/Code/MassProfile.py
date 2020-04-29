#Homework 5
# S Mackie ASTR 400B

# importing modules
import numpy as np
import astropy.units as u
from astropy.constants import G

# getting previously written functions
from ReadFile import Read
from CenterOfMass import CenterOfMass
from GalaxyMass import ComponentMass

# import plotting modules (from Lab 4)
import matplotlib.pyplot as plt
import matplotlib
get_ipython().run_line_magic('matplotlib', 'inline')


# Making the MassProfile class 
class MassProfile:

    # This will construct the filename read in the file, and store all mass and position    
    def __init__(self, galaxy, snap):
        # Inputs:
        #    galaxy = string with galaxy name MW, M31 or M33
        #    snap = snapshot number (ex. 0, 14, 38, etc)
        
        # Making the filename
        # Adds extra 0's to start of the inputted snap number
        ilbl = '000' + str(snap)
        # Gets rid of everything but the last 3 digits
        ilbl = ilbl[-3:]
        # Constructs the filename (ex. MW_000.txt) and saves it to a global property
        self.filename = "%s_"%(galaxy) + ilbl + '.txt'
        
        # Read data in the given file using Read
        self.time, self.total, self.data = Read(self.filename)
        
        # Storing the mass, the position, and the galaxy name as 
        # global properties
        self.m = self.data['m']
        self.x = self.data['x'] * u.kpc
        self.y = self.data['y'] * u.kpc
        self.z = self.data['z'] * u.kpc
        self.gname = galaxy
        
        # Putting G in the units we need to get the circ. vel. into km/s
        # making it a global property so I don't have to do this again
        self.G = G.to(u.kpc*u.km**2/u.s**2/u.Msun)
        
    # This function will calculate the mass enclosed w/in a given radius from the COM of 
    # the chosen galaxy and particle type
    def MassEnclosed(self, parttype, r):
        # Inputs:
        #    self = global properties
        #    parttype = particle type (1 = halo, 2 = disk, 3 = bulge)
        #    r = incremented radii (array) from the COM of the target galaxy (MW, M31, or M33) [kpc]
        # Returns:
        #    An array containing the mass enclosed (of a single particle type) w/in each radius [Msun]
        
        # Makes a CenterOfMass object
        COM = CenterOfMass(self.filename, 2)
        # Calling for the (vector) position of the center of mass [in kpc]
        COMP = COM.COM_P(0.1)
        # Calculates the magnitude of the COM position
        COMPmag = np.sqrt(COMP[0]**2 + COMP[1]**2 + COMP[2]**2)
        
        # calculates how far away each particle is from the galaxy's COM
        xMag = self.x - COMP[0]
        yMag = self.y - COMP[1]
        zMag = self.z - COMP[2]
        rMag = np.sqrt(xMag**2 + yMag**2 + zMag**2)
        
        # This array will store the total mass w/in each radius, where each element 
        # in SumMass corresponds to one radius 
        SumMass = np.zeros(len(r))

        #print(r)
        #print(len(rMag))
        #print(len(self.data['type']))
        
        # Loops through for each radius
        for i in range(len(r)):
            
            # This will hold the index of any particle that is within the specified radius 
            # AND is the particle type we want
            WithinR = np.where((rMag <= r[i]) & \
                               (self.data['type'] == parttype))
            
            # checking the loop step and the radius at that step
            #print("Iteration: ",i," Radius: ",r[i])
            
            # every particle w/in the radius will have its mass added to the sum
            SumMass[i] = sum(self.m[WithinR])
            
            # checking the elements along the way
            #print("SumMass", SumMass[i])
            
        return SumMass * 1e10 * u.Msun
    
    # This function will calculate the total mass enclosed w/in a given radius from the COM of
    # the chosen galaxy
    def TotalMassEnclosed(self, r):
        # Inputs:
        #    r = incremented radii (array) from the COM of the target galaxy (MW, M31, or M33) [kpc]
        # Returns:
        #    An array containing the total mass enclosed w/in each radius [Msun]
    
        # Initializing the arrays that will hold the masses for each type of particle
        MassHalo = np.zeros(len(r))
        MassDisk = np.zeros(len(r))
        MassBulge = np.zeros(len(r))
   
        # getting the mass arrays for the halo and disk
        MassHalo = self.MassEnclosed(1, r)
        MassDisk = self.MassEnclosed(2, r)
    
        # M33 doesn't have a bulge, so this if statement should leave its MassBulge array as zeros 
        if 'M33' not in self.filename:
            MassBulge = self.MassEnclosed(3, r)
        
        MassTot = MassHalo + MassDisk + MassBulge
    
        return MassTot
        
    # This function will calculate the mass enclosed w/in a given radius using the Hernquist 1990 mass profile:
    # M = Mhalo * r**2 / (a + r)**2
    # This function was taken from Lab3
    # In order to get the total halo mass for any particular galaxy, I think I need to use ComponentMass,
    #    but I'm not sure how to do that outside of this cell since ComponentMass needs the file name,
    #    which is put together in this class
    # Also I'm not sure how to use the scale radius without hard coding it in
    # For now I'll just leave the copied function as is for the Milky Way
    def HernquistMass(self, r, a, MHaloTot):
        # Inputs:
        #    r is the distance from the center of the galaxy [kpc]
        #    a is the scale radius [kpc]
        #    MHaloTot is the total dark matter halo mass [1e10 Msun]
        # Return:
        #    The total dark matter mass enclosed within radius r [Msun]
        
        r = r / u.kpc
    
        return np.round(MHaloTot * r**2 / (a + r)**2,2) * 1e10
    
    # This function will calculate the circular velocity at each radius for one particle type:
    # Vcirc = sqrt(G M / r)
    def CircularVelocity(self, parttype, r):
        # Inputs:
        #    parttype = particle type (1 = halo, 2 = disk, 3 = bulge)
        #    r = distance from the center of the galaxy [kpc]
        # Returns: 
        #    an array of circular speeds [km/s]
                
        # Getting the enclosed mass at each radius
        Menc = self.MassEnclosed(parttype,r)
        
        # Calculating Vcirc
        Vcirc = np.sqrt(self.G * Menc / r )
        
        return np.around(Vcirc,2)
    
    # This function will calculate the total circular velocity at each radius:
    # Vcirc = sqrt(G M / r)
    def CircularVelocityTotal(self, r):
        # Input:
        #    r = array of distances from the center of the galaxy [kpc]
        # Returns: 
        #    an array of circular speeds [km/s]
        
        # Getting the total enclosed mass at each radius
        MencTot = self.TotalMassEnclosed(r)
        
        # Calculating Vcirc
        Vcirc = np.sqrt(self.G * MencTot / r )
        
        return np.around(Vcirc,2)
    
    # This function will calculate the circular velocity using the Hernquist mass profile
    def HernquistVCirc(self, r, a, MHaloTot):
        # Inputs:
        #    r is the distance from the center of the galaxy [kpc]
        #    a is the scale radius [kpc]
        #    MHaloTot is the total dark matter halo mass [1e10 Msun]
        # Return:
        #    The circular speed [km/s]
        
        # Getting the enclosed mass at each radius
        Menc = self.HernquistMass(r, a, MHaloTot)
        
        # Need this so it won't yell at me about units
        # astropy units on HernquistG = [kpc]
        HernquistG = self.G * u.Msun * u.s**2 / u.km**2 
        
        # Calculating Vcirc
        Vcirc = np.sqrt(HernquistG * Menc / r )
        
        return np.around(Vcirc,2)
