
# coding: utf-8

# In[1]:

import sys
import numpy as np
get_ipython().magic('matplotlib inline')
import matplotlib.pyplot as plt
import math
import cmath
import scipy.io as sio
from scipy.integrate import odeint

import PeriodicOrbits as PO
from scipy.optimize import minimize, rosen, rosen_der


# In[2]:

#check if the para is correct
def getPeriod():
    return n




# In[4]:

#select 1/10 of the points
def getTenth(listQ):
    #let's make a quick movie of the advected particles
    se = 10  #sample every se time-steps
    #first we need to slice the data for plotting
    listRR = []
    for k in range(0,len(listQ[0]),se):
        listRR.append([[listQ[i][k][0] for i in range(0,len(listQ))],[listQ[i][k][1] for i in range(0,len(listQ))]])
    return listRR

# the transformation will be xnew = (y+R_offset)*cos(x*xmult) and ynew = (y+R_offset)*sin(x*xmult)
def AnnTrans(zin,p):
    x,y = zin
    roffset, xm, r0 = p
    return [(y+roffset)*math.cos(x*xm),(y+roffset)*math.sin(x*xm)]

# the transformation will be xnew = (y+R_offset)*cos(x*xmult) and ynew = (y+R_offset)*sin(x*xmult)
def AnnTransInve(zin,p):
    x,y = zin
    roffset, xm, r0 = p
    if(x<0):
        xp = (1/xm)*(math.atan(y/x)+math.pi)
        return [xp,x/(math.cos(xp*xm))-roffset]
    elif (y<0):
        #
        xp = (1/xm)*(math.atan(y/x)+math.pi*2)
        return [xp,x/(math.cos(xp*xm))-roffset]
    else:
        xp = (1/xm)*math.atan(y/x)
        return [xp,x/(math.cos(xp*xm))-roffset]


#Now we want to use this trajectory data as input into the 2D Etec algorithm.  
#To do this, we first need to project the data onto an annulus, with x as the angular variable, and y as a raduis
#Note that this projection only needs to be topologically accurate (since we are looking at topological quantities)
#We also put the trajectory data in the correct format to be input into Etec


#The stream function, and velocity vector functions
def StreamFunc(z,t,p):
    x,y = z
    U,L,k,c,eps = p
    psi0 = -U*L*math.tanh(y/L)
    psi1_0 = sum([eps[i]*cmath.exp(1j*k[i]*(x-c[i]*t)) for i in range(1,3)])
    psi1 = U*L*math.cosh(y/L)**(-2)*psi1_0.real
    return psi0 + psi1

def VelFunc(z,t,p):
    x,y = z
    U,L,k,c,eps = p
    vx0 = U*math.cosh(y/L)**(-2)
    vx1_0 = sum([eps[i]*cmath.exp(1j*k[i]*(x-c[i]*t)) for i in range(1,3)])
    vx1 = 2*U*math.cosh(y/L)**(-3)*math.sinh(y/L)*vx1_0.real
    vy1_0 = sum([eps[i]*1j*k[i]*cmath.exp(1j*k[i]*(x-c[i]*t)) for i in range(1,3)])
    vy1 = U*L*math.cosh(y/L)**(-2)*vy1_0.real
    return [vx0+vx1,vy1]


# In[5]:

# function
#input initial conditions and number of periods
#out put trajectory
def ICtoTRA(listIC, params):
    
    paraPhy,paraComp,paraTrans,n = params
    
    #x,y = AnnTransInve(z,paraTrans)
    
    #computational parameters
    T_i,T_tot,time =paraComp
    T_num = n*time
    T_f = n*T_tot  #Final time in hours (converted from days)
    dt = (T_f-T_i)/T_num   #The time-step to use
    #times = np.arange(T_i, T_f, dt)  #the set of times used in the ode solver
    times = np.linspace(T_i, T_f, T_num+1)
    #trabsformation parameters
    roffset, xm, r0 = paraTrans
    
    r0 = 6.371
    T_tot = math.pi*r0/(0.1*paraPhy[0])
    T_i = 0
    
    
    IC = listIC
    
    
    
    
    #Now we run each inital condition through the ode solver
    Traj = []
    for i in range(0,len(IC)):
        Traj.append(odeint(VelFunc, IC[i], times, args=(paraPhy,)))
        
    #Now mod the x-pos in each trajectory by pi*r0
    mv = math.pi*r0   #the value to mod by
    for i in range(0,len(Traj)):
        for j in range(0,len(Traj[0])):
            Traj[i][j][0] = Traj[i][j][0]%mv
            
    TrajAnn = []
    for j in range(0,len(Traj[0])):  #iterates first over the times
        TrajAnn.append([AnnTrans([Traj[i][j][0],Traj[i][j][1]],paraTrans) for i in range(0,len(Traj))])
        
    return TrajAnn
    


# In[6]:

#create a grid of IC
def gridIC(gridsize):
    #Here is how the functions are called:
    #StreamFunc([2000+2*math.pi*r0,2000],200000,params)
    #VelFunc([2000+1*math.pi*r0,2000],200000,params)
    #Now we set up the set of initial conditions (500x250 grid used in LCS overview paper) (used (18, 60), and (100,300))
    Nx, Ny = gridsize
    r0 = 6.371
    x_l = 0   #Rectangular Domain, x-range left (in Mm)
    x_r = math.pi*r0   #Rectangular Domain, x-range right (in Mm) ... this works out to be about 20 Mm
    y_b = -3   #Rectangular Domain, y-range bottom (in Mm)
    y_t = 3   #Rectangular Domain, y-range top (in Mm)
    Ntot = Nx*Ny  #the total number of initial points to seed
    dx = (x_r - x_l)/Nx  #the spacing between initial points in the x-direction
    dy = (y_t - y_b)/Ny  #the spacing between initial points in the x-direction
    IC = []
    for i in range(0,Ny):
        for j in range(0,Nx):
            IC.append([x_l+j*dx,y_b+i*dy])
            
    return IC


# In[7]:

#create a grid of IC
def gridICLess():
    #Here is how the functions are called:
    #StreamFunc([2000+2*math.pi*r0,2000],200000,params)
    #VelFunc([2000+1*math.pi*r0,2000],200000,params)
    #Now we set up the set of initial conditions (500x250 grid used in LCS overview paper) (used (18, 60), and (100,300))
    Nx = 20   #the number of initial columns to have along the x-direction
    Ny = 30    #The number of initial rows to have along the y-direction
    Ntot = Nx*Ny  #the total number of initial points to seed
    dx = (x_r - x_l)/Nx  #the spacing between initial points in the x-direction
    dy = (y_t - y_b)/Ny  #the spacing between initial points in the x-direction
    IC = []
    for i in range(0,Ny):
        for j in range(0,Nx):
            IC.append([x_l+j*dx,y_b+i*dy])
            
    return IC


# In[8]:

#plot circle
def plotCirecle(listW, m, n):

    namebase = 'jetData2D'
    extension = '.jpg'
    for i in range(m,n):
    #for i in range(0,len(aTslice)):
        plt.figure(figsize=(14,14))
        axes = plt.gca()
        axes.set_xlim([-10.1,10.1])
        axes.set_ylim([-10.1,10.1])
        axes.set_aspect('equal')
        plt.scatter(listW[i][0],listW[i][1],4.0)
        plt.show()
        filename = namebase + CounterToStr(i) + extension
        #plt.savefig(filename)
        plt.close()


# In[9]:

#len(aTslice = 100
def getPOrbit(listW,bound):
    #start place of points
    startL = listW[0]
    listxx = startL[0]
    listyy = startL[1]

    listDist = []
    listDist = PO.findD(listW)

    #points that travel back
    goodI2 = PO.findSmall(listDist,bound)

    print("this is the W: ",len(listW[0][0]))
    #points that move
    goodI1 = PO.getR(listW,0.3)

    #number of columns   
    line = PO.findLine(listW)  
    #overlapping part
    finalIndex = []
    finalIndex = PO.findIndex(listW, goodI2, goodI1,goodI1)
    #local minium
    #goodI4 = findLM(listDist,listW,finalIndex)


    traPO = []
    traPO = PO.findTra(listW,finalIndex)
    
    return traPO


# In[23]:

#plot
def plotF(traPO,a):
    
    figure(figsize=(9,9))
    xlabel('x')
    ylabel('y')
    string = str(a)
    title(string)
    #listxxx1 = traPO[12][0]
    #listyyy1 = traPO[12][1]
    #print(traPO)
    listxxx1 = []
    listyyy1 = []
    

    for x in range(0,len(traPO)):
        listxxx1.append(traPO[x][0])
        listyyy1.append(traPO[x][1])

    plt.scatter(listxxx1, listyyy1, s=1)
    
        #print(len(traPO[0]))
    listx0 = traPO[0][0]
    listy0 = traPO[0][1]

    plt.scatter(listx0, listy0, color = "C0")
    plt.savefig('PO1_strand'+str(a)+'.png')


# In[11]:

#plot
def plotEachTime(traPO):
    
    figure(figsize=(9,9))
    xlabel('x')
    ylabel('y')
    title('distance')
    #listxxx1 = traPO[12][0]
    #listyyy1 = traPO[12][1]
    #print(traPO)

    namebase = 'jetPO00'
    extension = '.png'
    #for each time step
    for x in range(0,len(traPO[0][0])):
        listxxx1 = []
        listyyy1 = []
        for k in range (0, len(traPO)):
            
            listxxx1.append(traPO[k][0][x])
            listyyy1.append(traPO[k][1][x])

        plt.scatter(listxxx1, listyyy1, s=5)
        plt.show()
        filename = namebase + str(x) + extension
        plt.savefig(filename)


# In[12]:

# get the surrounding points of a point
def getSur(point, bound):
    
    x,y = point
    listREt = []
    listRet.append([x+bound,y])
    listRet.append([x-bound,y])
    listRet.append([x,y+bound])
    listRet.append([x,y-bound])
    
    m = bound/(2^(1/2))
    listRet.append([x+m,y+m])
    listRet.append([x+m,y-m])
    listRet.append([x-m,y+m])
    listRet.append([x-m,y-m])
    
    
    return listRet


# In[13]:

#function to switch the later two indices
def swt23(listW):
    
    #same structure as the whole list [time][x=0/y=1][index]
    qVelo = []
    #start velocity set as 0
    for x in range(0,len(listW)):
        qVelo.append([])
        for m in range(0,len(listW[0][0])):
            qVelo[x].append([])
            for k in range(0,len(listW[0])):
                qVelo[x][m].append(0)
                
                
    #for each time step except the start time 0
    for y in range(0,len(qVelo)):
        
        for k in range(0,len(qVelo[0])):
        
            for z in range(0,len(qVelo[0][0])):
                qVelo[y][k][z] = listW[y][z][k]
                
    return qVelo
                
    
    
#function to switch the first two indices
def swt13(listW):
    
    #same structure as the whole list [time][x=0/y=1][index]
    qVelo = []
    #start velocity set as 0
    for x in range(0,len(listW[0][0])):
        qVelo.append([])
        for m in range(0,len(listW[0])):
            qVelo[x].append([])
            for k in range(0,len(listW)):
                qVelo[x][m].append(0)
                
                
    #for each time step except the start time 0
    for y in range(0,len(qVelo)):
        
        for k in range(0,len(qVelo[0])):
        
            for z in range(0,len(qVelo[0][0])):
                qVelo[y][k][z] = listW[z][k][y]
                
    return qVelo
    
    


# In[14]:

#get ini
def getIni(listW):
    
    listI = []
    for x in range(0,len(listW)):
        listI.append([listW[x][0][0],listW[x][1][0]] )
        
    return listI
    
    
#paramsComp = [T_i,T_tot,time]
def getDif(z, paramsAll):
    #unpack
    
    paraPhy, paraComp, paraTrans, n = paramsAll
    
    x,y = AnnTransInve(z,paraTrans)
    
    #computational parameters
    T_i,T_tot,time =paraComp
    T_num = n*time
    T_f = n*T_tot  #Final time in hours (converted from days)
    dt = (T_f-T_i)/T_num   #The time-step to use
    #times = np.arange(T_i, T_f, dt)  #the set of times used in the ode solver
    times = np.linspace(T_i, T_f, T_num+1)
    #trabsformation parameters
    roffset, xm, r0 = paraTrans
    
    IC = [x,y]
    #Now we run each inital condition through the ode solver
    zf = odeint(VelFunc, IC, times, args=(paraPhy,))[-1]
    #print(zf)
    mv = math.pi*r0   #the value to mod by
    zf[0] = zf[0]%mv
    #Now mod the x-pos in each trajectory by pi*r0
    zf = AnnTrans([zf[0],zf[1]],paraTrans)
    x,y = AnnTrans([x,y],paraTrans)
    
    dif = math.sqrt((zf[0]-x)**2 + (zf[1]-y)**2)
        
    return dif

#paramsComp = [T_i,T_tot,time]
def getEnd(z, paramsAll):
    #unpack
    
    paraPhy, paraComp, paraTrans, n = paramsAll
    
    x,y = AnnTransInve(z,paraTrans)
    
    #computational parameters
    T_i,T_tot,time =paraComp
    T_num = n*time
    T_f = n*T_tot  #Final time in hours (converted from days)
    dt = (T_f-T_i)/T_num   #The time-step to use
    #times = np.arange(T_i, T_f, dt)  #the set of times used in the ode solver
    times = np.linspace(T_i, T_f, T_num+1)
    #trabsformation parameters
    roffset, xm, r0 = paraTrans
    
    IC = [x,y]
    #Now we run each inital condition through the ode solver
    zf = odeint(VelFunc, IC, times, args=(paraPhy,))[-1]
    #print(zf)
    mv = math.pi*r0   #the value to mod by
    zf[0] = zf[0]%mv
    #Now mod the x-pos in each trajectory by pi*r0
    zf = AnnTrans([zf[0],zf[1]],paraTrans)
    x,y = AnnTrans([x,y],paraTrans)
    
    #dif = math.sqrt((zf[0]-x)**2 + (zf[1]-y)**2)
        
    return zf



#paramsComp = [T_i,T_tot,time]
def getNewTra(z, paramsAll):
    #unpack
    
    paraPhy, paraComp, paraTrans, n = paramsAll
    
    x,y = AnnTransInve(z,paraTrans)
    
    #computational parameters
    T_i,T_tot,time =paraComp
    T_num = n*time
    T_f = n*T_tot  #Final time in hours (converted from days)
    dt = (T_f-T_i)/T_num   #The time-step to use
    #times = np.arange(T_i, T_f, dt)  #the set of times used in the ode solver
    times = np.linspace(T_i, T_f, T_num+1)


    #trabsformation parameters
    roffset, xm, r0 = paraTrans
    
     
    IC = [x,y]
    #Now we run each inital condition through the ode solver
    zf = odeint(VelFunc, IC, times, args=(paraPhy,))
    #print(zf)
    mv = math.pi*r0   #the value to mod by
    for j in range(0,len(zf)):
        zf[j][0] = zf[j][0]%mv
    #Now mod the x-pos in each trajectory by pi*r0
    zf1 = []
    for k in range (0,len(zf)):
        zf1.append(AnnTrans([zf[k][0],zf[k][1]],paraTrans))
    
    
        
    return zf1


def plotEachT(listW):
    #print(listofHP) len(listofHP[0][0])
    #for each time step
    for x in range (0,len(listW)):
    #
        figure(figsize=(12,12))
        xlabel('x')
        ylabel('y')
        title('distance')
        listxxx1 = []
        listyyy1 = []
        #print(traPO)

        #for each point
        for m in range(0,len(listW[0][0])):
            listxxx1.append(listW[x][0][m])
            listyyy1.append(listW[x][1][m])

        
        axes = plt.gca()
        axes.set_xlim([-10,10])
        axes.set_ylim([-10,10])
        plt.scatter(listxxx1, listyyy1, s=10)
    
        #plt.show()
    
        plt.savefig('track'+str(x)+'.png')
    
    print("DONE")
    
    
    
#plot the whole track 
def plotAllWP(listW,listPO):
    #print(listofHP) len(listofHP[0][0])
    #for each time step
    for x in range (0,len(listW)):
        figure(figsize=(12,12))
        xlabel('x')
        ylabel('y')
        title('distance')
        listxxx1 = []
        listyyy1 = []
        #print(traPO)
        listynew = []
        listxnew = []

        #for each point
        for m in range(0,len(listW[0][0])):
            listxxx1.append(listW[x][0][m])
            listyyy1.append(listW[x][1][m])

        
        axes = plt.gca()
        axes.set_xlim([-10,10])
        axes.set_ylim([-10,10])
        plt.scatter(listxxx1, listyyy1, s = 7)
        
        
        for k in range(0,len(listPO[0][0])):
            listxnew.append(listPO[x][0][k])
            listynew.append(listPO[x][1][k])

        plt.scatter(listxnew, listynew, color = "C0")
        #print(len(listynew))
        #plt.show()
    
        plt.savefig('track'+str(x)+'.png')
    
    print("DONE")
    


    
def trim(listW):
    #listW[index][x or y]
    #get rid of points starting from the close to each other
    listreturn = []
    
    #for each point in the list
    for x in range(0,len(listW)):
        tempNum = 0
        #for each pair
        for y in range(x,len(listW)):
            #if the distance is too small, we treat them as the same orbit
            dis = math.sqrt((listW[x][0]-listW[y][0])**2 + (listW[x][1]-listW[y][1])**2)
            if(0<dis<0.01):
                tempNum = 1
        if(tempNum == 0):
            listreturn.append(listW[x])
            
    return listreturn


# In[17]:



#test if the two list are the same and return all the different lists
#in form [index][x or y][time]
def dropSame(listW,bound):
    #the list of good orbits
    newList = []
    #input the wholelist
    #get all the pair indices
    for x in range(0,len(listW)):
        isSame = 0
        #start with x+1 to avoid (m,m)
        for y in range(x+1,len(listW)):
            totalDis = 0
            #for each time step
            for t in range(0,len(listW[0][0])):
                #add the distnace
                totalDis = totalDis+ math.sqrt((listW[x][0][t]-listW[y][0][t])**2 + (listW[x][1][t]-listW[y][1][t])**2)
            #if the distance is large enough to distinguish two orbits
            if(totalDis<bound):
                isSame = 1
        #if an orbit is different from all other orbits
        if(isSame == 0):
            newList.append(listW[x])
    
    return newList


# In[ ]:

#func of getting the periodic orbits
def funcPO(gridsize, paramA, bound):
    #set the parameters
    #setPara(period,gridSize = None, returnrad = None)
    paramsPhys,paramsComp,paramstrans,n = paramA
    #starting grid
    initialPoint = gridIC(gridsize)
    #get track
    track = ICtoTRA(initialPoint,paramA)
    #transform to the right structure list[time][x/y][index]
    trackTXI = swt23(track)
    #get potential PO
    potentialPO = getPOrbit(trackTXI, bound)

    #get the initial positions of the potential PO
    ICPPO = getIni(potentialPO)


    #put them into the opt function to find zero function values
    listSt = []
    for x in range (0,len(ICPPO)):
        res = minimize(getDif, ICPPO[x], args=(paramA,), method='nelder-mead',options={'xtol': 1e-7, 'disp': True, 'maxiter': 10000})
        print("progress: ",x, "/", len(ICPPO))
        if(getDif([res.x[0],res.x[1]],paramA)<1e-3):
            listSt.append([res.x[0],res.x[1]])
        
    #trim the list so that close points are not both included
    listStN = trim(listSt)

    #use the new initial points and get the tra
    #for x in range
    listNewTra =[]
    for x in range (0,len(listStN)):
        listNewTra.append(getNewTra(listStN[x],paramA))

    print("before dropsame: ", len(listNewTra))
    listNewTra = dropSame(listNewTra,0.3)
    
    return listNewTra


# In[1]:

#starting grid
#initialPoint = gridIC()

#get track
#trackP1 = ICtoTRA(initialPoint,n,1000,paramsPhys)

#transform to the right structure list[time][x/y][index]
#trackPeriod1 = swt23(trackP1)

#get potential PO
#periodic1 = getPOrbit(trackPeriod1,0.4)

#get the initial positions of the potential PO
#ICP1 = getIni(periodic1)


#put them into the opt function to find zero function values
#listSt = []
#for x in range (0,len(ICP1)):
    #res = minimize(getDif, ICP1[x], args=(paramA,), method='nelder-mead',options={'xtol': 1e-7, 'disp': True, 'maxiter': 1000})
    #if(getDif([res.x[0],res.x[1]],paramA)<1e-3):
        
    
        #listSt.append([res.x[0],res.x[1]])
        
#trim the list so that close points are not both included
#listStN = trim(listSt)

#use the new initial points and get the tra
#for x in range
#listNewTra =[]
#for x in range (0,len(listStN)):
    #listNewTra.append(getNewTra(listStN[x],paramA))

#print("before dropsame: ", len(listNewTra))


# In[5]:

#listNewTra = dropSame(listNewTra,0.3)
#print(len(listNewTra))


# In[ ]:


    
    
#write the track of periodic orbits to the txt file
def writeToPO(listNewTra,filename):
    #write tra data to file
    #save the data to a txt file
    f = open(filename, 'a')
    #data is in the form list[index][x or y][time]  [[[T1x1,T1x2...],[y]],[t2],...]
    #chang it to time xy index
    listNum = swt23(listNewTra)
    listNum = swt13(listNum)
    #f.write()
    knx =''
    for x in range(0,len(listNum)):
        knx = knx+'\n'
        for y in range(0,len(listNum[0])):
        
            for z in range(0, len(listNum[0][0])):
            
                knx = knx+' '
                knx = knx+str(listNum[x][y][z])
                #print(knx)
    f.write(knx)
    f.close()
    print("done")


# In[4]:

#write tra data to file
#save the data to a txt file
#f = open("jetPO1.txt", "w");

#data is in the form list[index][x or y][time]  [[[T1x1,T1x2...],[y]],[t2],...]
#chang it to time xy index
#listNum = swt23(listNewTra)
#listNum = swt13(listNum)
#f.write()
#knx =''
#for x in range(0,len(listNum)):
    #knx = knx+'\n'
    #for y in range(0,len(listNum[0])):
        
        #for z in range(0, len(listNum[0][0])):
            
            #knx = knx+' '
            #knx = knx+str(listNum[x][y][z])
            #print(knx)
#f.write(knx)
#print("done")


# In[ ]:

import math
from math import pi


def getRing(center, r, n):
    
    #print(center)
    ring = [[center[0]+(math.cos(2 * pi / n * x) * r),  
            center[1] + (math.sin(2 * pi / n * x) * r)] for x in range(0, n + 1)]
    #print(ring)
    return ring

def getAng(n, i):
    
    #print(center)
    ring = (2 * pi / n * i) 
    
    #print(ring)
    return ring

#return [L1, L0], the longest and shortest distance from the center
def getL1L0(points,params,center):
    #points = getRing(point, ringR)
    listDiff = getDifls(points,params,center)
    sortedls = sorted(listDiff)
    
    return [sortedls[-1], sortedls[0]]

#get the ending positions 
def getEndls(points,params):
    listEnd = []
    
    
    for x in range (0, len(points)):
        #getDif(z, paramsAll):
        listEnd.append(getEnd(points[x], params))
        
    return listEnd
    
    
#get distance travelled of all points on the ring
def getDifls(points,params,center):
    listDiff = []
    
    #points = getRing(point, ringR)
    
    for x in range (0, len(points)):
        #getDif(z, paramsAll):
        tempEnd = getEnd(points[x], params)
        centerE = getEnd(center, params)
        dif = math.sqrt((tempEnd[0]-centerE[0])**2 + (tempEnd[1]-centerE[1])**2)
        listDiff.append(dif)
    return listDiff
    
#get the lambda for the two axes
#L0/radius = e^(lambda0*T)

def getLambda(L1L0, rad, times):
    L1, L0 = L1L0
    return [math.log(L1/rad)/times,math.log(L0/rad)/times]

#get the start positions for the track
def getStart(PO):
    listSt = []
    #get the start position
    l = list(PO[0])
    #get rid of the first zero
    l.pop(0)
    
    for x in range(0, math.floor(len(PO[5])/2)):
        listSt.append([float(l[x]),float(l[x+1])])
    return listSt
        
        

#return the list of two lambdas of each periodic orbits
def getPOlambda(PO, radius, params, times,n):
    listR = []
    PO = getStart(PO)
    for x in range (0, len(PO)):
        L10 = getL1L0(getRing(PO[x], radius,n),params,PO[x])
    
        listR.append(getLambda(L10, radius, times))
        print(x)
        print(listR[-1])
        
    return listR


# In[2]:



#plot each periodic orbit
#for x in range(0,len(listNewTra)):
    #plotF(listNewTra[x],x)


# In[ ]:

#get the images for running the gif
#def imageGIF():
    


# In[ ]:

#get the right structure
#newP1track = swt13(trackPeriod1)
#newP1track = swt23(newP1track)
#get 1/10 of the time steps
#allTenth1 = getTenth(newP1track)
#plotTenth = getTenth(listNewTra)

#make gifs
#plotEachT(plotTenth)
#plotEachT(allTenth1)
#plotAllWP(allTenth1,plotTenth)




# In[3]:

#print(listStN)


# In[ ]:



