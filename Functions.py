
# coding: utf-8

# In[ ]:

import copy
import math
import numpy as np
import Braid as B


def ExtraPointTopEntropy(FileIn):
    #get the data from file and store in variable BData
    BData = openN(FileIn)
    #change the data structure to Pos[x or y][strand number][time]
    #define vectors for each center and radius 
    a1 = [-1,0]
    a2 = [0,0]
    a3 = [1,0]
    radius = 0.1
    list1 = Pickpoint(BData,a1,radius)
    list2 = Pickpoint(BData,a2,radius)
    list3 = Pickpoint(BData,a3,radius)
    print((list1))
    print(len(list2))
    print(len(list3))
    indexlist = [list1[0],list2[0],list3[0]]
    
    return TopEntropy(BData,indexlist)
    
 



#Our Main Function
def BraidTopEntropy(FileIn):
    #get the data from file and store in variable BData
    BData = openN(FileIn)
    #getting the number of strands
    #NumStrands = numberstrands(BData)
    #getting the time list
    #Time = GetTime(BData)
    #change the data structure to Pos[x or y][strand number][time]
    #define vectors for each center and radius 
    a1 = [-1,0]
    a2 = [0,0]
    a3 = [1,0]
    
    radius = 0.1
    list1 = Pickpoint(BData,a1,radius)
    list2 = Pickpoint(BData,a2,radius)
    list3 = Pickpoint(BData,a3,radius)
    
    indexlist = [list1[0],list2[0],list3[0],list4[i]]
    

    Pos = GetPositions2(BData,indexlist)
    #sorted function k:XPos[k][tval]
    #SO = SortOrder(tval,XPos)
    #get the switch times (might need to call one more function inside this one)
    ST = GetSwitchTime(Pos)
    
  
    #get the alg representations without considering the signs (only x value needed)
    BG = GetSwitchPos(Pos,ST)
    #get the sgin of alg representations considering the Y which returns a list of 1 or -1 length = len(switch time)
    BS = GetBraidSigma(Pos,BG,ST)
    #the alg representation of braids
    BraidGen = GetBraidGen(BG,BS)
    return BraidGen
   
       
    

# abstr data from the datafile
def FileToData(filename):
    with open(filename,'r') as f:
        Data = []
        for line in f:
            Data.append(line.split(" "))
    return Data

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
    
#time .... list of time values
def GetTime(DataIn):
    Time = []
    for i in range(0,len(DataIn)):
        Time.append(DataIn[i][0])
    return Time
    
    
#number of strands, lots of ways to get this number such as len(pos[0])
def NumStrands(DataIn):
    
    return len(DataIn[0][0])
    #return (len(DataIn[0])-1)//2



def Pickpoint(DataIn, a, radius):
    biglist = []
    numstrands = NumStrands(DataIn)
    InitialPos = DataIn[0]
    #for i in range (0,):
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


#the function should return a big list of [X,Y,TopEn] for each point
def TopEntropy(DataIn,fixedlist):
    PosEntropy = []
    X = []
    Y = []
    numstrands = NumStrands(DataIn)
    for i in range (0,numstrands):
        if not i in fixedlist:
            Inlist = [i] + fixedlist
            Pos = GetPositions2(DataIn,Inlist)
            X.append(Pos[0][0][0])
            Y.append(Pos[1][0][0])
            #sorted function k:XPos[k][tval]
            #SO = SortOrder(tval,XPos)
            #get the switch times (might need to call one more function inside this one)
            ST = GetSwitchTime(Pos)
            #get the alg representations without considering the signs (only x value needed)
            BG = GetSwitchPos(Pos,ST)
            #get the sgin of alg representations considering the Y which returns a list of 1 or -1 length = len(switch time)
            BS = GetBraidSigma(Pos,BG,ST)
            #the alg representation of braids
            BraidGen = GetBraidGen(BG,BS)
            Braid = B.Braid(4,BraidGen)
            PosEntropy.append(Braid.DynTopEnt())
    return [X,Y,PosEntropy]
    
    
   




    



     
 


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
    #already founf switch time only need to compare the each element in (sortOrder) before and after switching 
    for t in range(0, len(ST)):
        a = SortOrder(ST[t],Pos[0])
        b = SortOrder(ST[t]+1,Pos[0])
        for s in range(0,numstrands):
            if (a[s] != b[s]):
            #need to do it for every element inside the SO and break it once we found one 
            #"the break commond dose not consider two or more crossing happens at the same time, might need to jump one and then check the next 
                BraidGen.append(s+1)
                break
    return BraidGen
    
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
          

