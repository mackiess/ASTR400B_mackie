#!/usr/bin/env python
# coding: utf-8

# In[10]:


# HW2 Part 2 Read File Program

import numpy as np          #imports NumPy module and renames it np in this script for ease of use
import astropy.units as u   #imports AstroPy module and renames it u in this script for ease of use

def Read(filename):           #defines 'Read' function to take in the name of whatever file we want to read

    file = open(filename,'r') #opens the file of interest to 'r'ead

    #first line
    line1 = file.readline()      #looks at the first line
    label, value = line1.split() #saves the label (Time) and value (0 right now) 
                                 #and they are split apart by spaces
    time = float(value)*u.Myr    #attaches the right units to new variable 'time'
    
    #second line
    line2 = file.readline()      #looks at the second line
    label, value = line2.split() #saves the label (Total (number of particles)) and value (135000)
                                 #and they are split apart by spaces
    totpart = float(value)       #total particles is the number it finds (as a float)
    
    file.close() #closes the MW_000.txt file
    
    data = np.genfromtxt(filename,dtype=None,names=True,skip_header=3)
    #new variable data that gets the column headers (line 4)
    #dtype=None -> white space split; names=True -> arrays named after header elements; 
    #skip_header=3 -> goes to the 4th line (skips first 3)
    
    return time, totpart, data
    #returns time value, total particles, and data header    


# In[13]:


#testing to make sure the function is reading things right

#time, totpart, data = Read("../../MW_000.txt") 
#Read(filename) gets the info from MW_000.txt and returns the 
#time value, total particles, and data header, which are pointed back to
#new variables in this cell that share the same names

#print(data['x'][2])


# In[ ]:




