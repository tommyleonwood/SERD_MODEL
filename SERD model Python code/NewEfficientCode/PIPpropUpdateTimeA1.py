from re import I
import matplotlib.pyplot as plt
import numpy as np
import random


TS=10
St=1
Et=1
propArray=[0,0]
PIPprop=[St]
tauEv=[]
StEv=[]
EtEv=[]
updateTime=0

autoSplit=bool(False)
merged=bool(False)

a=np.ones(St)

print(a)

def updatePIPprop(x,a,St):
    ind1=0
    pos=0
    for i in range(St):
        rand=random.random()
        if rand<1/3:
            x+=a
            St+=1
        elif rand<2/3:
            x-=a
            St-=1

        if i>x[pos]:
            a[pos]=0
            pos+=1
    
    x-=1
    return x

def updatePIPpropShift(t):
    global PIPprop
    global propArray
    global St
    global Et
    global tauEv
    global StEv
    global EtEv
    global updateTime

    shift=0
    pos=0
    ind=bool(True)
    for i in range(St):
        
        rand=random.random()
        if rand<1/3:
            propArray[i+shift]+=1
            propArray.insert(i+shift+1,0)

            PIPprop[pos]+=1
            shift+=1

            if ind==True:
                Et+=propArray[0]
                propArray[0]=0
                ind=False



        elif rand<2/3:
            propArray[i+shift]=propArray[i+shift]+propArray[i+shift+1]-1
            propArray.pop(i+shift+1)

            PIPprop[pos]-=1
            shift-=1

            if ind==True:
                if i+shift>0:
                    Et+=propArray[0]
                    propArray[0]=0
                    ind=False
      


        while i>PIPprop[pos+1]:

            if PIPprop[pos]<1:
                tau=t-updateTime-i
                tauEv.append(tau)
                EtEv.append(Et)
                updateTime+=1
                PIPprop.pop(0)
                pos-=1
            else:
                PIPprop[pos]-=1

            PIPprop[pos+1]+=shift
            pos+=1

    
    St+=shift
    StEv.append(St)


    propArray.pop(0)
    propArray.append(0)

    return PIPprop





def updatePIPpropShiftB(t):
    global PIPprop
    global propArray
    global St
    global Et
    global tauEv
    global StEv
    global EtEv
    global updateTime

    shift=0
    pos=0
    ind=bool(True)
    for i in range(St):
        
        rand=random.random()
        if rand<1/3:
            propArray[i+shift]+=1
            propArray.insert(i+shift+1,0)

            PIPprop[pos]+=1
            shift+=1

        elif rand<2/3:
            propArray[i+shift]=propArray[i+shift]+propArray[i+shift+1]-1
            propArray.pop(i+shift+1)

            PIPprop[pos]-=1
            shift-=1

        while i>PIPprop[pos+1]:

            if PIPprop[pos]<1:
                tau=t-updateTime-i
                tauEv.append(tau)
                EtEv.append(Et)
                updateTime+=1
                PIPprop.pop(0)
                pos-=1
            else:
                PIPprop[pos]-=1

            PIPprop[pos+1]+=shift
            pos+=1


    
    St+=shift
    Et+=propArray[0]

    StEv.append(St)




    propArray.pop(0)
    propArray.append(0)

    return PIPprop




def MergedIE():
    global St
    global merged
    global PIPprop
    global propArray

    if St==0:
        if(autoSplit==True):
            merged=False
            St=1
            PIPprop=[St]
            propArray=[0,0]
        else:
            merged=True

            



def runSingleIEPathEvFromUnionWithSplitA1(TS):
    global updateTime
    global St
    global Et
    global propArray
    global PIPprop
    global tauEv
    global StEv
    global EtEv

    St=1
    Et=1
    propArray=[0,0]
    PIPprop=[St]
    tauEv=[]
    StEv=[]
    EtEv=[]
    updateTime=0


    for t in range(TS):

        a=np.ones(St)

        if t > 0:
            PIPprop.append(St)
            
            updatePIPpropShift(t)

        print(PIPprop)


runSingleIEPathEvFromUnionWithSplitA1(TS)

        




