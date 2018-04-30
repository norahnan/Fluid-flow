#This is an attempt to create a braid class in python based on my previous c++ braid class.  I will eventually add extra functionality like braid conjugacy determination. Jan 2017 - Spencer A. Smith

#This class will hold the Artin group representation of a Braid.  It is able to compare braids (solves the word problem) via an order (Dehornoy, or other).
#Operations on the braid will change the word representative, but not the braid (i.e. adhere to the braid relations)
#The braids should be thought of as going from left to right (corresponding to the order of the written generators), and the positive generators will take the
#lesser indexed strand under the greater indexed strand.
#The positive generators will be represented by 1,2,3,4, ... , 0 is the identity, and the inverses are -1,-2,-3,....

#the standard modules to import
import numpy as np
import copy   #needed to have a copy function that makes a deep copy
import math

#This is the braid class 
class Braid:
    #any global shared variables to use?
    def __init__(self, bn = 3, braidword = [0], btype = None, bwsize = 0):
        self.Bnum = bn   #The number of strands in the braid
        if btype == None:
            #the number of words that comprise the braid are given in the lenth of the Bword list
            self.Bword = braidword[:]  #The braid word represented as a list
        elif btype == 'random':
            print("Random braid code to be written")
            self.Bword = []
            #create a random braid with bwsize number of letters
        elif btype == 'delta':
            self.Bword = []
            for i in range(0,self.Bnum-1):
                self.Bword.append(i+1)
        elif btype == 'Delta':
            self.Bword = []
            for i in range(0,self.Bnum-1):
                for j in range(0,self.Bnum-1-i):
                    self.Bword.append(j+1)
        else:
            #all others will be identity
            self.Bword = [0]
            
        self.Sym = []  #This holds the symmetric group member corresponding to the braid
        self.UpdateSym()   #This initializes Sym to the correct symmetric group member

    #Insead of a copy constructor, can use this ... i.e. b = a.copy()
    def copy(self):
        return copy.deepcopy(self)    
    
    #overloaded operators
    
    #Multiplication, understood to be braid composition
    def __mul__(self,other):
        B = []
        if self.Bnum == other.Bnum:
            if self.IsID():
                B = other.Bword[:]
            elif other.IsID():
                B = self.Bword[:]
            else:
                B = self.Bword + other.Bword
        else:
            B = [0]
            print("braids do not have the same number of strands.  Can not compose them")
        return Braid(self.Bnum,B)
    
    # *= overloaded
    def __imul__(self,other):
        self = self*other
        return self
    
    # overload raising to a power
    def __pow__(self,exp):
        exp2 = abs(int(exp))
        result = self.copy()
        for i in range(0,exp2-1):
            result *= self
        return result
            
    #Less Than boolean operator overloaded
    def __lt__(self,other):
        Red = self.Inverse()
        Red2 = Red*other
        less = False
        Red2.FHreduce()
        if Red2.IsSigmaPos():
            less = True
        return less
    
    #Greater than boolean operator overloaded
    def __gt__(self,other):
        Red = other.Inverse()
        Red2 = Red*self
        greater = False
        Red2.FHreduce()
        if Red2.IsSigmaPos():
            greater = True
        return greater
    
    #Equivalence operator overloaded
    def __eq__(self,other):
        Red = self.Inverse()
        Red2 = other*Red
        equal = False
        Red2.FHreduce()
        if (not Red2.IsSigmaPos()) and (not Red2.IsSigmaNeg()):
            equal = True
        return equal
    
    #Returns True if Identity
    def IsID(self):
        IDbool = False
        if len(self.Bword) == 1 and self.Bword[0] == 0:
            IDbool = True
        return IDbool
        
    #Updates Sym to be the symmetric group member corresponding to this braid
    def UpdateSym(self):
        if len(self.Sym) == 0:
            for i in range(0,self.Bnum):
                self.Sym.append(i+1)
        elif len(self.Sym) == len(self.Bword):
            for i in range(0,self.Bnum):
                self.Sym[i] = i+1
        else:
            self.Sym = []
            for i in range(0,self.Bnum):
                self.Sym.append(i+1)
                
        if not self.IsID():
            for i in range(0,len(self.Bword)):
                swid = abs(self.Bword[i])
                temp = self.Sym[swid]
                self.Sym[swid] = self.Sym[swid-1]
                self.Sym[swid-1] = temp      
    
    #define what is returned when print is called
    def __str__(self):
        return ' '.join(str(e) for e in self.Bword)
    
    #Sigma positive determination (for determining order)
    def IsSigmaPos(self):
        sp = False
        MinI = self.MinIndex()
        if self.Contains(MinI) and not self.Contains(-1*MinI):
            sp = True
        return sp
    
    def IsSigmaNeg(self):
        sp = False
        MinI = self.MinIndex()
        if self.Contains(-1*MinI) and not self.Contains(MinI):
            sp = True
        return sp
        
        
    #contains ... see if the braid contains the given generator
    def Contains(self, match):
        contains = False
        if not self.IsID():
            for i in range(0,len(self.Bword)):
                if self.Bword[i] == match:
                    contains = True
                    break
        return contains
    
    #returns the Maximum index for the braid (highest index of a generator)
    def MaxIndex(self):
        MI = 0
        if not self.IsID():
            for i in range(0,len(self.Bword)):
                if abs(self.Bword[i]) > MI:
                    MI = abs(self.Bword[i])
        return MI
    
    #returns Minium index for the braid
    def MinIndex(self):
        MI = self.Bnum - 1
        if not self.IsID():
            for i in range(0,len(self.Bword)):
                if abs(self.Bword[i]) < MI:
                    MI = abs(self.Bword[i])
        return MI
    
    #returns true if the braid is a pure braid (symmetry group rep is id)
    def IsPure(self):
        pure = True
        for i in range(0,self.Bnum):
            if not self.Sym[i] == i+1:
                pure = False
                break
        return pure
    
    #This returns the braid that is the inverse of the current braid (
    #together they are the identity element)
    def Inverse(self):
        Binv = self.copy()
        if not Binv.IsID():
            Binv.Invert()
            Binv.Negative()
        return Binv
    
    #this send sigma(i) to sigma(n-i) .... Warning! This function changes the braid
    def Invert(self):
        if not self.IsID():
            self.Bword.reverse()
            self.UpdateSym()
    
    #this replaces every generator with its negative ... Warning! This function changes the braid
    def Negative(self):
        if not self.IsID():
            self.Bword = [-1*x for x in self.Bword]
    
    #Insert one braid into another
    #convention: new braid will be inserted just before Ipos in the current braid word.
    #size is from from. !!This changes the braid
    def Insert(self, Ipos, Bfrom):
        if self.IsID():
            self = Bfrom.copy()
        else:    
            if not Bfrom.IsID():            
                fws = len(Bfrom.Bword)
                for i in range(0,fws):
                    self.Bword.insert(Ipos+i,Bfrom.Bword[i])
                self.Bnum = max(self.Bnum,Bfrom.Bnum)
                self.UpdateSym()
        
    #Remove a section of the braid
    #convention: the letters Ipos1 through Ipos2 will be removed.
    def Remove(self,Ipos1,Ipos2):
        if Ipos1 <= Ipos2:
            if Ipos1 >= len(self.Bword) or Ipos2 < 0:
                print("Warning, range does not overlap with any Braid index")
            else:
                Imin = max(Ipos1,0)
                Imax = min(Ipos2,len(self.Bword)-1)
                diff = Imax - Imin + 1
                for i in range(0,diff):
                    self.Bword.pop(Imin)
                if len(self.Bword) == 0:
                    #this is the case where everything is removed, and we give the identity
                    self.Bword.append(0)
                self.UpdateSym()
        else:
            print("Warning, not removing any sub-braid, as Ipos1 > Ipos2")
    
    #this returns the sub-braid in this braid with indices from p to q
    def Subword(self,p,q):
        SW = self.copy()
        if not SW.IsID():
            if q >= p:
                SWws = len(SW.Bword)
                if q < SWws - 1:
                    SW.Remove(q+1,SWws-1)
                if p > 0:
                    SW.Remove(0,p-1)
            else:
                #give identity
                SW.Remove(0,len(SW.Bword)-1)
        return SW
    
                    
    #finds the first handle starting from the left
    def FindHandleFromLeft(self,P):
        found = False
        if not self.IsID():
            ws = len(self.Bword)
            MinI = self.MinIndex()
            for i in range(0,ws):
                if abs(self.Bword[i]) == MinI:
                    P[0] = i
                    break
            for i in range(P[0]+1,ws):
                if self.Bword[i] == self.Bword[P[0]]:
                    P[0] = i
                elif self.Bword[i] == -1*self.Bword[P[0]]:
                    P[1] = i
                    found = True
                    break
        return found
    
    #finds the left most handle
    def FindLeftHandle(self,P):
        found = False
        if not self.IsID():
            if self.FindHandleFromLeft(P):
                Bsub = self.Subword(P[0]+1,P[1]-1)
                Pp = [0,0]
                if Bsub.FindLeftHandle(Pp):
                    P[1] = P[0] + Pp[1] + 1
                    P[0] = Pp[0] + 1
                found = True
        return found
    
    #Single Handle Reduction
    def SHreduce(self):
        reduction = False
        if not self.IsID():
            MinI = self.MinIndex()
            if not (self.IsSigmaPos() or self.IsSigmaNeg()):
                #i.e. there is a handle
                P = [0,0]
                if self.FindLeftHandle(P):
                    Ind = abs(self.Bword[P[0]])
                    Isign = -1
                    if self.Bword[P[0]] == Ind:
                        Isign = 1
                    Piece = self.Subword(P[0]+1,P[1]-1)
                    self.Remove(P[0],P[1])
                    negsw = 1
                    for i in range(len(Piece.Bword)-1,-1,-1):
                        if Piece.Bword[i] == Ind + 1:
                            negsw = 1
                        elif Piece.Bword[i] == -1*(Ind + 1):
                            negsw = -1
                        if abs(Piece.Bword[i]) == Ind + 1: #for both cases
                            Bw = []
                            Bw.append(-1*Isign*(Ind+1))
                            Bw.append(negsw*Ind)
                            Bw.append(Isign*(Ind+1))
                            SubPiece = Braid(Piece.Bnum,Bw)
                            Piece.Insert(i,SubPiece)
                            Piece.Remove(i+3,i+3)
                    self.Insert(P[0],Piece)
                else:
                    print("Not sigma positive or negative, not identity, and no handle found!","\n")
                reduction = True
        return reduction
    
    #Single Normal Reduction: single reduction: finds an instance of sigma(i)*sigma(i,-) (if any) and reduces it.  
    #Returns true if a reduction has taken place, and false otherwise
    def SNreduce(self):
        reduction = False
        if not self.IsID():
            for i in range(0,len(self.Bword)-1):
                if self.Bword[i] == -1*self.Bword[i+1]:
                    self.Remove(i,i+1)
                    reduction = True
                    break
        return reduction
    
    #Full normal Reduction
    def FNreduce(self):
        reducible = True
        while reducible:
            reducible = self.SNreduce()
    
    
    #Full Handle Reduction
    def FHreduce(self):
        reduction = False
        if not self.IsID():
            looping = True
            #print("Before FNreduce", self.Bword)
            self.FNreduce()  #First do a full normal reduction
            #print("After FNreduce", self.Bword)
            #then loop a single handle reduction with a subsequent 
            #full normal reduction while there is something to reduce
            while looping:
                looping = self.SHreduce()
                #print("After SHreduce",self.Bword)
                self.FNreduce()
                #print("After FNreduce",self.Bword)
                if looping:
                    reduction = True
        return reduction
    
    #Some Simple functions which return basic info about the braid
    #Returns generator at index i
    def getBword(self,i):
        return self.Bword[i]
    
    #Returns the wordsize
    def getwordsize(self):
        return len(self.Bword)
    
    #Returns the number of strands
    def getBnum(self):
        return self.Bnum
    
    
    #exponent sum is the number of positive generators minus the number of negative generators
    #This integer is a conjugacy invariant
    def ExponentSum(self):
        es = 0
        if not self.IsID():
            for i in range(0,self.getwordsize()):
                if self.Bword[i] < 0:
                    es -= 1
                else:
                    es += 1
        return es

    #This conjugates the braid (a) with the given braid (b) as bab^-1
    def Conjugate(self,Conj):
        result = self.copy()
        if not result.IsID():
            result = Conj*result*(Conj.Inverse())
            result.UpdateSym()
        return result
        
##################################################################
# The Burau representation for finding the topological entropy for 3-stranded braids
    def BurauTopEnt(self):
        if self.Bnum == 3:
            s1 = np.array([[1, 1], [0, 1]])
            sn1 = np.linalg.inv(s1)
            s2 = np.array([[1,0],[-1,1]])
            sn2 = np.linalg.inv(s2)
            Mtot = np.array([[1,0],[0,1]])
            for x in self.Bword:
                if x == 1:
                    Mtot = np.dot(Mtot,s1)
                elif x == -1:
                    Mtot = np.dot(Mtot,sn1)
                elif x == 2:
                    Mtot = np.dot(Mtot,s2)
                elif x == -2:
                    Mtot = np.dot(Mtot,sn2)
                else:
                    print("Incorrect Generator is part of 3-stranded Braid: ",x)
            Eval, Evec = np.linalg.eig(Mtot)
            return math.log(max(abs(Eval[0]),abs(Eval[1])))
                
        else:
            print("This is not a 3-stranded braid, and therefore the Burau estimate of the topological entropy is not appropriate")
            return False
        
        
        
##################################################################
# Now for the Dynnikov coordinate ideas (still part of the braid class)

    #This contains the code to calculate an estimate of the topological entropy of a braid based on the Dynnikov coordinate representation of curves. The basic idea is to pick a random dynnikov coordinate (of appropriate size for the braid in consideration) and have the braid repeatedly act on the coordinates. The number of times the resultant curve intersects a horizontal line going through all points will increase exponentially with period (and therefore needs to be renormalized occasionally). This exponential rate of increase is the topological entropy. It should be relatively insensitive to the particular random dynnikov coordinate that we start with (do few, incase the random curve only interacts with a finite order sub-component of a reducible braid).


    # The update function
    #Note that unlike the convention that Thiffeault has used, I interpret the positive generator to switch pairs of points
    #in a CCW direction (viewed from above), the braid is read/(acts) left to right, and time for the geometric braid
    #flows down.  This has the effect of switching the positive and negative generators as compared to Thiffeault and
    #others (I'm using Birman's conventions).

    def DynUpdate(self,sigma,ui):
        #first calculate the number of strands implicit in the size of ui (2n-4)
        nt = int((len(ui)+4)/2)
        #make a copy of ui to modify
        uout = np.copy(ui)
        #uout = ui[:]
        #cut up ui to help with notation
        a = ui[:nt-2]
        b = ui[nt-2:]
        if sigma < 0:
            nsigma  = -1*sigma
            sigmaind = nsigma - 1
            if nsigma == 1:
                uout[sigmaind] = -1*b[sigmaind] + max((a[sigmaind] + max(b[sigmaind],0)),0)
                uout[sigmaind+(nt-2)] = a[sigmaind] + max(b[sigmaind],0)
            elif nsigma == nt-1:
                uout[sigmaind-1] = -1*b[sigmaind-1] + min((a[sigmaind-1] + min(b[sigmaind-1],0)),0)
                uout[sigmaind-1+(nt-2)] = a[sigmaind-1] + min(b[sigmaind-1],0)
            else:
                c = a[sigmaind-1] - a[sigmaind] - max(b[sigmaind],0) + min(b[sigmaind-1],0)
                uout[sigmaind-1] = a[sigmaind-1] - max(b[sigmaind-1],0) - max((max(b[sigmaind],0)+c),0)
                uout[sigmaind-1+(nt-2)] = b[sigmaind] + min(c,0)
                uout[sigmaind] = a[sigmaind] - min(b[sigmaind],0) - min((min(b[sigmaind-1],0)-c),0)
                uout[sigmaind+(nt-2)] = b[sigmaind-1] - min(c,0)
        elif sigma > 0:
            sigmaind = sigma-1
            if sigma == 1:
                uout[sigmaind] = b[sigmaind] - max((max(b[sigmaind],0)-a[sigmaind]),0)
                uout[sigmaind+(nt-2)] =  max(b[sigmaind],0) - a[sigmaind]
            elif sigma == nt-1:
                uout[sigmaind-1] = b[sigmaind-1] - min((min(b[sigmaind-1],0)-a[sigmaind-1]),0)
                uout[sigmaind-1+(nt-2)] = min(b[sigmaind-1],0) - a[sigmaind-1]
            else:
                d = a[sigmaind-1] - a[sigmaind] + max(b[sigmaind],0) - min(b[sigmaind-1],0)
                uout[sigmaind-1] = a[sigmaind-1] + max(b[sigmaind-1],0) + max((max(b[sigmaind],0)-d),0)
                uout[sigmaind-1+(nt-2)] = b[sigmaind] - max(d,0)
                uout[sigmaind] = a[sigmaind] + min(b[sigmaind],0) + min((min(b[sigmaind-1],0)+d),0)
                uout[sigmaind+(nt-2)] = b[sigmaind-1] + max(d,0)
        else:
            print("Braid Generator is not allowed to be 0")
        return uout


    #this implements the above update on the Dynnikov coordinates for the action of the full braid betai
    def DynFullUpdate(self,uin):
        uout = np.copy(uin)
        if not self.IsID():
            for i in range(0,len(self.Bword)):
                uout = self.DynUpdate(self.Bword[i],uout)
        return uout
    
    
    #This defines the number of intersections the given curve (represented by the Dynnikov coordinates) makes with
    #a horizontal line that goes through each point and terminates on the side boundaries (or at infinity)
    #This measure mimics the exponential increase in total length, and leads to the same topological entropy
    def DynIntL(self,ui):
        nt = int((len(ui)+4)/2) #total number of strands
        #cut up ui to help with notation
        a = ui[:nt-2]
        b = ui[nt-2:]
        L = abs(a[0]) + abs(a[nt-3])
        Ltemp = 0
        for i in range(0,nt-3):
            Ltemp += abs(a[i+1]-a[i])
        L += Ltemp
        b0i = abs(a[0]) + max(b[0],0)
        bacc = 0
        for i in range(1,nt-2):
            bacc += b[i-1]
            b0itemp = abs(a[i]) + max(b[i],0) + bacc
            if b0itemp > b0i:
                b0i = b0itemp
        b0 = -1*b0i
        bftemp = 0
        btemp = 0
        for i in range(0,nt-2):
            bftemp += b[i]
            btemp += abs(b[i])
        bf = -1*b0 - bftemp
        L += btemp + abs(b0) + abs(bf)
        return L
        

    #This will be the main function that take in a braid word (with specified number of strands) and spits out the 
    #approximation of the topological entropy
    def DynTopEnt(self, ErrorTol = 0.000000001, MaxIter = 1000):
        #first create random Dynnikov coordinates of the proper size
        DynNum = 2*self.Bnum - 4
        MaxRand = 20
        uinit = []
        for i in range(0,DynNum):
            uinit.append(np.random.random_sample()*2 - 1)
        u = np.array(uinit)    
        #this is the Large Number that triggers the renormalization procedure
        RenormNum = 281474976710656  #This is 2^48 ... conservative for the 53 bit accuracy
        RenormIter = 0   #the number of times the renormalization has been done
        #Remark ... after renormalization, the elements of u will generically be real insead of integers.  This presents
        #no problem other than the fact that the physical interpretation of the coordinates is no longer true (crossings)
        #MaxIter = 1000  #The maximum number of iterations to try ... should halt earlier due to decreasing error.
        MinIter = 5    #the minimum number of iterations to get rid of transients ... start log slope analysis here
        #ErrorTol = 0.000000001 
        #10^-9 ... if differences in adjacent approximations of h (top ent) (relative to h) are less than this, then halt
        #h = []   #The array holding the sucessively better approximations of the top ent
        h1 = 0
        h2 = 0
        #h.append(0)
        hf = 0
        for i in range(0,MinIter):
            u = self.DynFullUpdate(u)
            maxval = max(max(u),abs(min(u)))
            #check if renorm is needed
            if maxval > RenormNum:
                u = u/RenormNum
                RenormIter += 1
        L0 = self.DynIntL(u)
        hp0 = np.log(L0*RenormNum**RenormIter)
        #print(hp0)
        #tolerance trip flag
        TolFlag = False
        #main loop
        for i in range(0,MaxIter-MinIter):
            u = self.DynFullUpdate(u)
            maxval = max(max(u),abs(min(u)))
            #check if renorm is needed
            if maxval > RenormNum:
                u = u/RenormNum
                RenormIter += 1
            #now find L
            L = self.DynIntL(u)
            #calc h and add to list
            h1 = h2
            h2 = (np.log(L)+RenormIter*np.log(RenormNum)-hp0)/(i+1)
            #print(h2)
            #h.append(h2)
            #now check to see if the error size halting critera has been reached
            if abs(h2) == 0:
                ErrorVal = abs(h2-h1)
            else:
                ErrorVal = abs(h2-h1)/(h2)  #divide by h[i+1] since h[0] is == 0
            hf = h2
            if ErrorVal < ErrorTol:
                TolFlag = True
                print("Tolerance of ", ErrorTol," reached after ", i+MinIter," Iterations","\n")
                break
        if TolFlag == False:
            print("The specified tolerance of ",ErrorTol, " was not reached after the maximum of ", MaxIter, " iterations.","\n")
            print("The Topological Entropy value given is only an upper bound","\n")
        return hf