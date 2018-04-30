
# coding: utf-8

# In[1]:

import sys
#sys.path.append("/Volumes/SmithLab/SmithLab/Modules")
import BraidFuncM as BF
import Functions as Func
import Braid as B
import matplotlib.pyplot as plt
import numpy as np
get_ipython().magic('matplotlib inline')
import math
import cmath
import scipy.io as sio
import itertools


# In[2]:

U_b = 0.06266*60*60/1000 #This is a characteristic velocity measured in Mm/hour (converted from km/s)
r0 = 6.371   #The mean radius of Earth in Mm
T_tot = math.pi*r0/(0.1*U_b) 


# In[ ]:

def getAlg(listW):
    return BF.GeoToAlgBraid(listW)


# In[4]:

def getEnt(filename,period,PON):
    n = period
    time = T_tot*n
    wholeList = BF.openN(filename)
    BraidlistP1Sall = BF.GeoToAlgBraid(filename)
    BraidP1Sall = B.Braid(PON,BraidlistP1Sall)
    enT = BraidP1Sall.DynTopEnt()
    TotalEntro = enT/time
    return TotalEntro


# In[6]:

def findindex(num,total):
    stuff = list(range(total))
    indices = []
    for subset in itertools.combinations(stuff,num):
        indices.append(list(subset))
    return indices


# In[7]:

#use the list  
def getEntro(listW,filename):
    Braidlist = BF.GeoToAlgBraid(filename,listW)
    Braid = B.Braid(len(listW),Braidlist)
    #entriNumber = Braid.DynTopEnt()
    #print("this is entropy: " + entriNumber)
    return Braid.DynTopEnt()


# In[1]:

#get the highest topological entropy from these groups of indices
def getBestValue(indices, listW,filename):
    listEn = []
    for x in range(0,len(indices)):
        listEn.append(getEntro(indices[x],filename))
    return max(listEn)


# In[ ]:

#get the actual combination of the highest ent indices
def getBestIndex(indices, listW, filename):
    listEn = []
    for x in range(0,len(indices)):
        listEn.append(getEntro(indices[x],filename))
    return indices[listEn.index(max(listEn))]

