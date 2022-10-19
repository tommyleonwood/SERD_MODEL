import matplotlib.pyplot as plt
import numpy as np
import random

pd=1/3
pr=1/3
pt=pr+pd

measurements=int(1)

iterations=int(1000)

TS=int(100)

TS1=[]
sumTS=0
countTS=[0]
for i1 in range(measurements+1):
    TS1.append(TS)
    sumTS+=TS1[i1]
    countTS.append(countTS[i1]+TS1[i1])

sampleSize=[]
for i1 in range(measurements+1):
    sampleSize.append(0)

obsStLoc=[]
obsStHalfRange=[]
obsEtLoc=[]
obsEtHalfRange=[]

StEtObservationRatio=float(1)
timeScaleFactor=float(1)
filterWidthFactor=float(0.5)

backToNormalObservations=int(0)

for i1 in range(measurements+1):
    if(i1<measurements+1-backToNormalObservations):
        obsStLoc.append(timeScaleFactor*StEtObservationRatio*countTS[i1+1])
        obsStHalfRange.append(timeScaleFactor*StEtObservationRatio*countTS[i1+1]*filterWidthFactor)
        obsEtLoc.append(timeScaleFactor*countTS[i1+1])
        obsEtHalfRange.append(timeScaleFactor*countTS[i1+1]*filterWidthFactor)
    else:
        obsStLoc.append(timeScaleFactor*countTS[i1+1])
        obsStHalfRange.append(timeScaleFactor*countTS[i1+1]*filterWidthFactor)
        obsEtLoc.append(timeScaleFactor*countTS[i1+1])
        obsEtHalfRange.append(timeScaleFactor*countTS[i1+1]*filterWidthFactor)

print(obsStLoc)
print(obsStHalfRange)
print(obsEtLoc)
print(obsEtHalfRange)



addToRangeLimit=int(1000000)

ObsStEvTot=[]
ObsEtEvTot=[]
ObsTauEvTot=[]
ObsSampleSize=[]

for i1 in range(sumTS):
    ObsStEvTot.append(0)
    ObsEtEvTot.append(0)
    ObsTauEvTot.append(0)
    ObsSampleSize.append(0)

stEvTot=[]
etEvTot=[]
etstEvTotSampleSize=[]
for i1 in range(sumTS):
    stEvTot.append(0)
    etEvTot.append(0)
    etstEvTotSampleSize.append(0)

allowAutoSplits=bool(False)

for i2 in range(iterations):
    
    
    
    merged=False
    propArray=[int(0),int(0)]
    St=int(1)
    Et=int(1)
    PIPprop=[St]
    updateTime=int(0)

    for t1 in range(TS1[0]):
        shift=[int(0)]

        if merged==False:

            for x1 in range(St):
                rand=float(random.random())
                if rand < pd:
                    propArray[x1+shift[x1]]+=1
                    propArray.insert(x1+shift[x1]+1,0)
                    shift.append(shift[x1]+1)
                elif rand < pt:
                    propArray[x1+shift[x1]]=propArray[x1+shift[x1]]+propArray[x1+shift[x1]+1]-1
                    propArray.pop(x1+shift[x1]+1)
                    shift.append(shift[x1]-1)
                else:
                    shift.append(shift[x1])

            St2=St+shift[St]
            Et2=Et+propArray[0]

            if St==0:
                if allowAutoSplits==False:
                    merged=True
                else:
                    St=1
                    Et=1
                    propArray=[int(0),int(0)]
                    shift=[int(0)]
                    PIPprop=[St]
            else:
                ind1=int(0)
                for x2 in range(St):
                    while PIPprop[ind1]==x2:
                        PIPprop[ind1]+=shift[x2]-1
                        ind1+=1
                        if ind1>St-1:
                            break
                    if ind1>St-1:
                        break

                TSupdateTime=int(0)
                for t2 in range(t1-updateTime):
                    if PIPprop[t2]<1:
                        TSupdateTime+=1
                    else:
                        break

                for i2 in range(TSupdateTime):
                    PIPprop.pop(0)
                    ObsStEvTot[t1]+=St2
                    ObsEtEvTot[t1]+=Et
                    ObsTauEvTot[t1]+=t1-updateTime-i2
                    ObsSampleSize[t1]+=1

                updateTime+=TSupdateTime
                St=St2
                PIPprop.append(St)

                Et=Et2
                etEvTot[t1]+=Et
                stEvTot[t1]+=St
                etstEvTotSampleSize[t1]+=1
                propArray.pop(0)
                propArray.append(0)

    if St>obsStLoc[0]-obsStHalfRange[0]:
        if St<=obsStLoc[0]+obsStHalfRange[0]+addToRangeLimit:
            if Et>obsEtLoc[0]-obsEtHalfRange[0]:
                if Et<=obsEtLoc[0]+obsEtHalfRange[0]+addToRangeLimit:
                    inputChosen=bool(True)


                else:
                    inputChosen=bool(False)
            else:
                inputChosen=bool(False)
        else:
            inputChosen=bool(False)
    else:
        inputChosen=bool(False)

    sampleSize[0]+=1

    allowAutoSplits=bool(False)

    if(inputChosen==True):

        for i1 in range(measurements):

            if(inputChosen==True):
                sampleSize[i1+1]+=1

                for t1 in range(TS1[i1+1]):
                    shift=[int(0)]

                    if merged==False:
                        for x1 in range(St):
                            rand=float(random.random())

                            if rand < pd:
                                propArray[x1+shift[x1]]+=1
                                propArray.insert(x1+shift[x1]+1,0)
                                shift.append(shift[x1]+1)

                            elif rand < pt:
                                propArray[x1+shift[x1]]=propArray[x1+shift[x1]]+propArray[x1+shift[x1]+1]-1
                                propArray.pop(x1+shift[x1]+1)
                                shift.append(shift[x1]-1)

                            else:
                                shift.append(shift[x1])

                        St2=St+shift[St]
                        Et2=Et+propArray[0]

                        if St==0:
                            if allowAutoSplits==False:
                                merged=True

                            else:
                                St=1
                                Et=1
                                propArray=[int(0),int(0)]
                                shift=[int(0)]
                                PIPprop=[St]
                        else:
                            ind1=int(0)
                            for x2 in range(St):
                                while PIPprop[ind1]==x2:
                                    PIPprop[ind1]+=shift[x2]-1
                                    ind1+=1
                                    if ind1>St-1:
                                        break

                                if ind1>St-1:
                                    break

                            TSupdateTime=int(0)
                            for t2 in range(TS1[i1+1]*(i1+1)+t1-updateTime):
                                if PIPprop[t2]<1:
                                    TSupdateTime+=1
                                else:
                                    break

                            for i2 in range(TSupdateTime):
                                PIPprop.pop(0)
                                ObsStEvTot[countTS[i1+1]+t1]+=St2
                                ObsEtEvTot[countTS[i1+1]+t1]+=Et
                                ObsTauEvTot[countTS[i1+1]+t1]+=countTS[i1+1]+t1-updateTime-i2
                                ObsSampleSize[countTS[i1+1]+t1]+=1

                            updateTime+=TSupdateTime
                            St=St2
                            PIPprop.append(St)

                            Et=Et2
                            etEvTot[countTS[i1+1]+t1]+=Et
                            stEvTot[countTS[i1+1]+t1]+=St
                            etstEvTotSampleSize[countTS[i1+1]+t1]+=1
                            propArray.pop(0)
                            propArray.append(0)
                            
                if St>obsStLoc[i1+1]-obsStHalfRange[i1+1]:
                    if St<=obsStLoc[i1+1]+obsStHalfRange[i1+1]+addToRangeLimit:
                        if Et>obsEtLoc[i1+1]-obsEtHalfRange[i1+1]:
                            if Et<=obsEtLoc[i1+1]+obsEtHalfRange[i1+1]+addToRangeLimit:
                                inputChosen=bool(True)
                            else:
                                inputChosen=bool(False)
                        else:
                            inputChosen=bool(False)
                    else:
                        inputChosen=bool(False)
                else:
                    inputChosen=bool(False)
            else:
                break

normValObsStEvTot=[]
normValObsEtEvTot=[]
normValObsTauEvTot=[]

normStEvTot=[]
normEtEvTot=[]

for i1 in range(sumTS):
    if ObsSampleSize[i1]>0:
        normValObsStEvTot.append(ObsStEvTot[i1]/ObsSampleSize[i1])
        normValObsEtEvTot.append(ObsEtEvTot[i1]/ObsSampleSize[i1])
        normValObsTauEvTot.append(ObsTauEvTot[i1]/ObsSampleSize[i1])
    else:
        normValObsStEvTot.append(0)
        normValObsEtEvTot.append(0)
        normValObsTauEvTot.append(0)

for i1 in range(sumTS):
    if etstEvTotSampleSize[i1]>0:
        normStEvTot.append(stEvTot[i1]/etstEvTotSampleSize[i1])
        normEtEvTot.append(etEvTot[i1]/etstEvTotSampleSize[i1])
    else:
        normStEvTot.append(0)
        normEtEvTot.append(0)


print(normValObsStEvTot)
print(normValObsEtEvTot)
print(normValObsTauEvTot)

print(sampleSize)

y1 = np.array(normValObsStEvTot)
y2 = np.array(normValObsEtEvTot)
y3 = np.array(normValObsTauEvTot)

x1 = np.array(normStEvTot)
x2 = np.array(normEtEvTot)


plt.plot(y1)
plt.plot(y2)
plt.plot(y3)
plt.plot(x1)
plt.plot(x2)
plt.show()