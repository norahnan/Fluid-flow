
# coding: utf-8

# In[1]:

#imports
import math
# This line configures matplotlib to show figures embedded in the notebook, 
# instead of opening a new window for each figure. More about that later. 
# If you are using an old version of IPython, try using '%pylab inline' instead.
get_ipython().magic('matplotlib inline')
from pylab import *
import matplotlib
import matplotlib.pyplot as plt
import numpy as np

from matplotlib import pyplot as PLT
from matplotlib import cm as CM
from matplotlib import mlab as ML
import numpy as NP


# In[2]:

#open and record data from a file

def openN(fileName):

    #a list of time with evolving x,y coordinates
    wholeList = []

    #open and record
    with open(fileName,"r") as f:
        #for each time
        for line in f:
            #coordinates
            listX = []
            listY = []
            a = line.split(" ")
            #delete the first element which is the time
            words = a[1:]
            m = round((len(words))/2)
            for x in range(0, m):
                listX.append(float(words[2*x]))
                listY.append(float(words[2*x+1]))
            #append xy of each time
            listTime = []
            
            listTime.append(listX)
            listTime.append(listY)
            
            wholeList.append(listTime)
    return wholeList
#print(lister[8])



# In[1]:

#get the distance for each point
def findD(list):
    listD = []
    
    #get start time positions
    listS = list[0]
    #start x
    listS0 = listS[0]
    #start y
    listS1 = listS[1]
    
    #get end
    listE = list[-1]
    #end x
    listE0 = listE[0]
    #end y
    listE1 = listE[1]

    
    #for each pair
    m = len(listS1)
    #print (len(listS))
    for x in range (0,m):
        
        x1 = listS0[x]
        y1 = listS1[x]
        xn = listE0[x]
        yn = listE1[x]
        d = ((xn-x1)**2+(yn-y1)**2)**(1/2)
        listD.append(d)
    return listD



# In[ ]:








# In[5]:


#gridsize=30
#PLT.subplot(111)

# if 'bins=None', then color of each hexagon corresponds directly to its count
# 'C' is optional--it maps values to x-y coordinates; if 'C' is None (default) then 
# the result is a pure 2D histogram 
#x = listxx
#y = listyy
#z = [2,5,3,5,7]
#gridsize=gridsize,
#PLT.hexbin(x, y, C=listDist,  cmap=CM.jet, bins=None)
#PLT.axis([x.min(), x.max(), y.min(), y.max()])

#cb = PLT.colorbar()
#cb.set_label('mean value')
#PLT.show() 


# In[6]:



#cut the list according to the good index
#given a list of indeces, get the according elements from the original list in [[x1,x2...][y1,y2...]]
def cutList(listW,index):
    #hold the new list
    list1 = []
    #time
    for m in range (0,len(listW)):
        listime = []
        listx = []
        listy = []
        #go through the index
        for x in range (0,len(index)):
            #add every pair to the acoording time
            listx.append(listW[m][0][index[x]])
            listy.append(listW[m][1][index[x]])
        #add all the pairs to the time
        listime.append(listx)
        listime.append(listy)
        #add time to the larger list
        list1.append(listime)
    
    #return the trimmed list
    return list1




#get rid of start coordinates of points that does not move much during the whole time
#return the coordinates of moving points
def getR(listW,bound):
    #list of good index of points
    listRest = []
    #loop through all the points
    for y in range(0,len(listW[0][0])):
        
        #start place
        x0 = listW[0][0][y]
        y0 = listW[0][1][y]
        boo = 0
        #append place of all time
        for x in range (0,len(listW)):
            #dis from start through the whole time
            d = ((listW[x][0][y]-x0)**2+(listW[x][1][y]-y0)**2)**(1/2)
            #if the point goes out
            
            if(d>bound):
                #get out of the loop
                x = len(listW)
                #change bool
                boo = 1
        
            
        #if the coordinates do go out of the boundary
        if(boo == 1):
            listRest.append(y)
    return listRest




# In[7]:


#find paths of points with a distance smaller than a certain value
#input the list of all data & distance

def findSmall(listD, bound):
    listPa = []
    for x in range(0,len(listD)):
        if(listD[x]<bound):
            
            listPa.append(x)
    return listPa


    


# In[8]:

def getDis(m1,m2):
    k = ((m1[0]-m2[0])**2+(m1[1]-m2[1])**2)**(1/2)
    return k

#find start coordinates of points in the stirring disk
def findDisk(listW,d1,d2,d3,r):
    listStart = []
    #temp = findStartP(listW)

    for y in range(0,len(listW[0][0])):
        #go through the whole time
        onePoint = []
        #start place
        x0 = listW[0][0][y]
        y0 = listW[0][1][y]
        
        onePoint.append(x0)
        onePoint.append(y0)
        #print("this are the coordinates",x0,y0)
        #if point is not within the disk
        if(getDis(d1,onePoint)>r and getDis(d2,onePoint)>r and getDis(d3,onePoint)>r):
            #print(y,"this is the distance",getDis(d1,onePoint))
            listStart.append(y)
    
    #out put good index
    return listStart


# In[9]:

#get whole track of certain points in [[[T1x1, T1x2],[T1y1,T1y2]][t2]...]
def findTra(listW,index):
    #all points
    listR = []
    
    for x in range (0,len(listW[0][0])):
        
        if (x in index):
            #all time
            listx = []
            listy = []
            for y in range (0,len(listW)):
                listx.append(listW[y][0][x])
                #print("append",listW[y][0][x])
                listy.append(listW[y][1][x])
            onePoint = []
            #print(listx)
            onePoint.append(listx)
            onePoint.append(listy)
            listR.append(onePoint)
        
    return listR


# In[10]:


def findLine(listW):
    num = 0
    for x in range(1,len(listW[0][0])):
        temp = listW[0][0][x]
        #if we have 
        if (temp < listW[0][0][x-1]):
            return x
            
#get the index that is in all three index lists
def findIndex(listW,i1,i2,i3):
    #for every point
    listt = []
    for x in range (0,len(listW[0][0])):
        #print(x in i1 and x in i2 and x in i3)

        if x in i1 and x in i2 and x in i3:
            listt.append(x)
    return listt




#fina local minium
#line is the longer row
def findLM(listD,listW,index):
    #list index of these points
    listM = []
    
    #for each point
    for x in range (0,len(index)):
        m = index[x]
        #group every two rows to make each group have the same number of points
        ingroup = m%(line*2+1)
        groupN = m/(line*2+1)
        #find where in the group is the point
        #if(groupOFtwo<line):
        listAround = [m,m+1,m-1,m+line-1,m+line,m-line,m-line+1]
        listSquare = []
        
        for l in range (0,len(listAround)):
            #print("this is len", len(listD))
            #print(l,listAround[l])
            #if index is legal
            if (listAround[l]<len(listD)):
                 
                listSquare.append(listD[listAround[l]])
            
            
        minIndex0 = listSquare.index(min(listSquare))
        minIndex = listAround[minIndex0]
        if(minIndex == index[x]):
            listM.append(minIndex)
            
            
            #remove same entries
            
        #for ones are not the smallest,remove
        #for k in range (0,len(listAround)):
            #if in the min list and not min
            #if(listAround[k] in listM && listAround[k] != minIndex):
                #listM.remove(listAround[k])
            
    return listM



# In[ ]:




# In[12]:


#print(traPO[0])


# In[13]:




# In[ ]:




# In[ ]:




# In[ ]:



