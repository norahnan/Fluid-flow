
# coding: utf-8

# In[179]:

import copy
import math
import numpy as np
import Braid as B
import sys
get_ipython().magic('matplotlib inline')
import matplotlib.pyplot as plt
import math
import cmath
import scipy.io as sio

#Our Main Function
def GeoToAlgBraidDT(BData,indexlist = None):
    #get the data from file and store in variable BData
    #BData = openN(FileIn)
    #getting the number of strands
    numStrands = NumStrands(BData)
    #getting the time list
    #Time = GetTime(BData)
    #change the data structure to Pos[x or y][strand number][time]
    #define vectors for each center and radius 
    
    if(indexlist == None):
        indexlist = []
        for x in range(0,numStrands):
            indexlist.append(x)
   

    Pos = GetPositions2(BData,indexlist)
    #print(len(Pos))
    #print(len(Pos[0]))
    #print(len(Pos[0][0]))
    #sorted function k:XPos[k][tval]
    #SO = SortOrder(tval,XPos)
    #get the switch times (might need to call one more function inside this one)
    ST = GetSwitchTime(Pos)
    
  
    #get the alg representations without considering the signs (only x value needed)
    #BG = GetSwitchPos(Pos,ST)
    #get the sgin of alg representations considering the Y which returns a list of 1 or -1 length = len(switch time)
    #BS = GetBraidSigma(Pos,BG,ST)
    #the alg representation of braids
    BraidGen = PosToSwi(Pos,ST)
    return BraidGen





#Our Main Function
def GeoToAlgBraid(FileIn,indexlist = None):
    #get the data from file and store in variable BData
    BData = openN(FileIn)
    #getting the number of strands
    numStrands = NumStrands(BData)
    #getting the time list
    #Time = GetTime(BData)
    #change the data structure to Pos[x or y][strand number][time]
    #define vectors for each center and radius 
    
    if(indexlist == None):
        indexlist = []
        for x in range(0,numStrands):
            indexlist.append(x)
   

    
    Pos = GetPositions2(BData,indexlist)
    #sorted function k:XPos[k][tval]
    #SO = SortOrder(tval,XPos)
    #get the switch times (might need to call one more function inside this one)
    ST = GetSwitchTime(Pos)
    #print(len(ST),"lenth of ST")
    
  
    #get the alg representations without considering the signs (only x value needed)
    #BG = GetSwitchPos(Pos,ST)
    #get the sgin of alg representations considering the Y which returns a list of 1 or -1 length = len(switch time)
    #BS = GetBraidSigma(Pos,BG,ST)
    #the alg representation of braids
    BraidGen = PosToSwi(Pos,ST)
    return BraidGen


#tranform list from time xxxxyyy to time xyxyxyxyxy
def tranformL(listW):
    listR = []
    half = (len(listW[0])-1)//2
    #print(half)
    for x in range(0,len(listW)):
        temp = [0]
        for y in range(0,half):
            temp.append(listW[x][y+1])
            temp.append(listW[x][y+1+half])
        #print(len(temp))
        listR.append(temp)
            
    return listR
        
    
    

def openN(fileName):

    #a list of time with evolving x,y coordinates
    wholeList = []
    

    #open and record
    with open(fileName,"r") as f:
        #for each time
        temp = 0
        for line in f:
            #throw away the first line
            if(temp == 0):
                #print(temp)
                temp = 1
            else:
                
                a = line.split(" ")
                #delete the first element which is the space
                words = a[1:]
                #print(len(words))
                words = [0] + words
                #print(words[0][0])
            
                wholeList.append(words)
    
    
    wholeList = tranformL(wholeList)
    #print(len(wholeList[0]))
    #print(len(wholeList[1]))
    return wholeList  
    

def openNNewF(fileName):

    #a list of time with evolving x,y coordinates
    wholeList = []

    #open and record
    with open(fileName,"r") as f:
        #for each time
        temp = 0
        for line in f:
            #throw away the first line
            if(temp == 0):
                #print(temp)
                temp = 1
            else:
                #print(temp)
                #coordinates
                listX = []
                listY = []
                a = line.split(" ")
                #delete the first element which is the space
                words = a[1:]
                #print(len(words))
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
    
    
# abstr data from the datafile
def FileToData(filename):
    with open(filename,'r') as f:
        Data = []
        for line in f:
            Data.append(line.split(" "))
    return Data
    
#time .... list of time values
def GetTime(DataIn):
    Time = []
    for i in range(0,len(DataIn)):
        Time.append(DataIn[i][0])
    return Time
    
    
#number of strands, lots of ways to get this number such as len(pos[0])
def NumStrands(DataIn):
    return (len(DataIn[0])-1)//2




def Pickpoint(DataIn, a, radius):
    biglist = []
    numstrands = NumStrands(DataIn)
    InitialPos = DataIn[0][1:]
    for i in range (0,numstrands):
        #need x and y position for each point 
        ri = [float(InitialPos[2*i]),float(InitialPos[2*i+1])]
        xdiff = ri[0] - a[0]
        ydiff = ri[1] - a[1]
        s = math.sqrt(xdiff**2 + ydiff**2)
        if s < radius:
            #biglist append the index or these points 
            biglist.append(i)
    return biglist

#list4 concludes indexes for all points at time t=0 
#need to exclude three points
def pickPo(DataIn):
    list4=[]
    InitialPos = DataIn[0][1:]
    for i in range(0,len(InitialPos)):
        list4.append(i)
    return list4       



    
#data structure, pos[x or y][strands number][by time steps]
def GetPositions(DataIn):
    X = []
    Y = []
    numstrands = NumStrands(DataIn)
    for x in range (0,numstrands):
        Xtemp = []
        Ytemp = []
        X.append(Xtemp)
        Y.append(Ytemp)
    for i in range (0,numstrands):
        for j in range(0,len(DataIn)):
            X[i].append(float(DataIn[j][2*i+1]))
            #X[i].append(-1*float(DataIn[j][2*i+1]))
            #should be 2*i+1, change x to -x to test
            Y[i].append(float(DataIn[j][2*i+2]))
    Positions = []
    Positions.append(X)
    Positions.append(Y)
    return Positions    
 


#data structure, pos[x or y][strands number][by time steps]
#gets only the data associated with the input indices
def GetPositions2(DataIn,IndexList):
    X = []
    Y = []
    for x in range (0,len(IndexList)):
        Xtemp = []
        Ytemp = []
        X.append(Xtemp)
        Y.append(Ytemp)
    for i in range (0,len(IndexList)):
        for j in range(0,len(DataIn)):
            X[i].append(float(DataIn[j][2*IndexList[i]+1]))
            Y[i].append(float(DataIn[j][2*IndexList[i]+2]))
    Positions = []
    Positions.append(X)
    Positions.append(Y)
    return Positions




#"this returns a list of integers corresponding to the sorting order of Xall at time index tval"
def SortOrder(tval,XPos):
    numstrands = len(XPos)
    #creats an ordered list of integers 
    OrderID = []
    for i in range(0,numstrands):#(0,numstrands) means (0,1,...,numstrands - 1)
        OrderID.append(i)
    #sorted OrderID based on how Xall sorted
    return sorted(OrderID, key=lambda k:XPos[k][tval])
    
#find the swich times (the indices) and find the alg representaion without considering the Ys. 
#Index is the index before the switch
def GetSwitchTime(Pos):
    SwitchTime = []
    a = SortOrder(0,Pos[0])
    for t in range(1,len(Pos[0][0])):
        b = SortOrder(t,Pos[0])
        if (a!=b):
            SwitchTime.append(t-1)
            a = copy.copy(b)   #look here first if this function does not work ... might need a deep copy of b
    return SwitchTime

#briadgen without considering the - or + signs, only need x value at this point to see which pair is crossing
def GetSwitchPos(Pos,ST):
    BraidGen = []
    numstrands = len(Pos[0])
    #already found switch time only need to compare the each element in (sortOrder) before and after switching 
    for t in range(0, len(ST)):
        a = SortOrder(ST[t],Pos[0])
        b = SortOrder(ST[t]+1,Pos[0])
        for s in range(0,numstrands):
            #base case: when there is one switch in the range
            #number of different orders in total
            temp = 0
            #the order number
            listS = []
            if (a[s] != b[s]):
                temp = temp + 1
                listS.append([a[s], b[s]])
        #if two order differences and their places are switched
        if(temp == 2 and lisS[0] == [listS[1][1],listS[1][0]]):
            BraidGen.append(s+1)
                
        
            
            
            
    return BraidGen

#
def PosToSwi(Pos,ST):
    #search for each pair
    #first  and last where is switch
    #time should be inside the 
    #update the ordering after each switch
    
    BraidGen = []
    numstrands = len(Pos[0])
    #already found switch time only need to compare the each element in (sortOrder) before and after switching 
    for t in range(0, len(ST)):
        Ta = ST[t]
        Tb = ST[t]+1
        
        #strabd num list
        a = SortOrder(Ta,Pos[0])
        b = SortOrder(Tb,Pos[0])
        
        #get all the pairs between the first and last switch
        #index of the strand number list
        pairs = findPairList(a,b,numstrands)
        #print((pairs), "length of pairs")
        #intersections [[switch time, strandP, strandQ],[],[],...]
        Isecs = []
        #find the pairs that has a intersection
        #for each pair
        #pos[x or y][strands number][by time steps]
        for x in range(0,len(pairs)):
            #the time
            
            
            #for point P and Q
            #print(pairs[0], Ta)
            Pa = [Pos[0][a[pairs[x][0]]][Ta],Ta]
            Pb = [Pos[0][a[pairs[x][0]]][Tb],Tb]
            Qa = [Pos[0][a[pairs[x][1]]][Ta],Ta]
            Qb = [Pos[0][a[pairs[x][1]]][Tb],Tb]
            
            
            lineP = findKB(Pa, Pb)
            lineQ = findKB(Qa, Qb)
            
            IS = findInter(lineP, lineQ)
            #if there is an intersection and the time is in the range
            if not IS is None:
                if(IS[1] >= Ta and IS[1] <= Tb):
                    #append time, first strand , second strand
                    Isecs.append([IS[1],a[pairs[x][0]],a[pairs[x][1]]])
        
        #sort the list of intersections accourding to switch time
        #print(Isecs)
        sortedIS = sorted(Isecs)
        #temp for transforming
        transTemp = []
        for m in range(0,len(a)):
            transTemp.append(a[m])
        
        #append the very first switch
        #print(sortedIS[0][1],sortedIS[0][2])
        #timeList start int end
        #print("this is time:", ST[t])
        #print("ST[t]",ST[t])
        #print("length of sortedIS",len(sortedIS))
        #print("sortedIS",(sortedIS))

        #get start time int time and end time
        timeList = [Ta,sortedIS[0][0],Tb]
        BraidGen.append(getYSign(timeList,Pos,[sortedIS[0][1],sortedIS[0][2]])*(transTemp.index(sortedIS[0][1])+1))
        #update the tranform list
        pair = [transTemp.index(sortedIS[0][1]),transTemp.index(sortedIS[0][2])]
        transTemp[pair[0]], transTemp[pair[1]] = transTemp[pair[1]], transTemp[pair[0]]
        #for each switch
        for x in range(1,len(sortedIS)):
            #the two strands in their original order
            s1 = sortedIS[x][1]
            s2 = sortedIS[x][2]
            timeListX = [Ta,sortedIS[x][0],Tb]
            BraidGen.append(getYSign(timeListX,Pos,[s1,s2])*(transTemp.index(s1)+1))
            #update the tranform list
            sw1 = transTemp.index(s1)
            sw2 = transTemp.index(s2)
            transTemp[sw1], transTemp[sw2] = transTemp[sw2], transTemp[sw1]
            
            
    return BraidGen

#input the switch time, position list and which two strands switches
def getYSign(t, Pos, pair):
    #use inte y and t
    startT, insT, endT = t
    y1a = Pos[1][pair[0]][startT]
    y1b = Pos[1][pair[0]][endT]
    y2a = Pos[1][pair[1]][startT]
    y2b = Pos[1][pair[1]][endT]
    
    line1 = findKB([startT,y1a],[endT,y1b])
    line2 = findKB([startT,y2a],[endT,y2b])
    
    #print([startT,y1a],[endT,y1b])
    int1 = insT*line1[0] + line1[1]
    int2 = insT*line2[0] + line2[1]
    #print(int1, int2)
    
    #if y at pair[0] is larger
    if(int1>int2):
        return 1
    else:
        return -1
            

            
#return the line function given two points on the line
def findKB(a,b):
    #a = [x1,y1]
    #b = [x2,y2]
    x1,y1 = a
    x2,y2 = b
    #print("a and b")
    #print(a,b)
    k = (y1-y2)/(x1-x2)
    b = y1 -k*x1
    lineKB = [k,b]
    return lineKB



#find intersection of two lines
def findInter(l1,l2):
    #l = [k,b]
    #if the slope is the same
    if(l1[0] == l2[0]):
        print("the slope is too close")
        return None
    point = []
    x = (l2[1]-l1[1])/(l1[0]-l2[0])                
    y = l1[0]*x +l1[1]    
    point = [x,y]  
    return point       

def findPairList(a,b,numstrands):
    
    startS = -10
    endS = -10
    isStart = 0
    isEnd = 0
    
    #check from the beginning
    for s in range (0,numstrands):
        if (a[s] != b[s] and isStart == 0):
            startS = s
            isStart = 1
            
            
        
    
    #check from the end
    for e in range(1,numstrands+1):
        #print(a[numstrands-e],b[numstrands-e])
        if (a[numstrands-e] != b[numstrands-e] and isEnd == 0):
            #print(a[numstrands-e],b[numstrands-e])
            endS = numstrands-e
            isEnd = 1
        
    #print(startS,endS)
    #make into pairsS
    returnL = getPairs(startS, endS)
    # 
    return returnL
            
            
            

                    
def getPairs(start, end):
    #listSingle = list(np.linspace(start, end, num=(end-start+1)))
    #list[pairnumber][first or second]
    listSingle = [start+i for i in range (0,end-start+1)]
    #list = [[a,b],[c,d],[],...]
    returnList = []
    
    for x in range(0,len(listSingle)):
        for y in range(x+1,len(listSingle)):
            returnList.append([listSingle[x], listSingle[y]])
            
    return returnList
        
    
def GetBraidSigma(Pos,BG,ST):
    #return list of -1 or 1 should coresponding to the y position and which strand cross in front of the other
    BraidSigma = []
    numstrands = len(Pos[0])
    XiL = []
    XiR = []
    XfL = []
    XfR = []
    YiL = []
    YiR = []
    YfL = []
    YfR = []
    indL = []
    indR = []
    for x in range (0,len(BG)):
        indL.append(SortOrder(ST[x],Pos[0])[BG[x] - 1])
        indR.append(SortOrder(ST[x],Pos[0])[BG[x]])
        XiL.append(Pos[0][indL[x]][ST[x]])
        XiR.append(Pos[0][indR[x]][ST[x]])
        XfL.append(Pos[0][indR[x]][ST[x] + 1])
        XfR.append(Pos[0][indL[x]][ST[x] + 1])
        YiL.append(Pos[1][indL[x]][ST[x]])
        YiR.append(Pos[1][indR[x]][ST[x]])
        YfL.append(Pos[1][indR[x]][ST[x] + 1])
        YfR.append(Pos[1][indL[x]][ST[x] + 1])
    TRatio = []
    for x in range (0,len(BG)):
        TRatio.append((XiR[x] - XiL[x])/(XiR[x] - XiL[x] + XfR[x] - XfL[x]))
    YBraidL = []
    YBraidR = []
    for x in range (0,len(BG)):
        YBraidL.append((YiL[x]-YfR[x]) * TRatio[x] + YiL[x])
        YBraidR.append((YiR[x]-YfL[x]) * TRatio[x] + YiR[x])
    Braidsigma = []
    for x in range (0,len(BG)):
        if YBraidL[x] < YBraidR[x]:
            Braidsigma.append(1)
        else:
            Braidsigma.append(-1)
    return Braidsigma

def GetBraidGen(BG,BS):
    Braids = []
    for x in range (0,len(BG)):
        Braids.append(BG[x] * BS[x])
    return Braids
          


# In[174]:




# In[180]:




# In[181]:




# In[ ]:




# In[ ]:




# In[ ]:




# In[ ]:




# In[ ]:




# In[ ]:




# In[36]:




# In[31]:




# In[33]:




# In[ ]:




# In[ ]:



