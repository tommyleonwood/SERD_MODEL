from cmath import exp
from pickle import EXT2
import matplotlib.pyplot as plt
import numpy as np
import random
from numpy import random

pd=1/3
pr=1/3
pt=pr+pd

measurements=int(1)
ObservtionProbabilityFactor=float(0.005)

iterations=int(10000000)

TS=int(500)

obsStLoc=[]
obsStHalfRange=[]
obsEtLoc=[]
obsEtHalfRange=[]

StEtObservationRatio=float(1)
timeScaleFactor=float(1)
filterWidthFactor=float(0.5)

backToNormalObservations=int(0)



stEvTot=[]
etEvTot=[]
etstEvTotSampleSize=[]
for i1 in range(TS):
    stEvTot.append(0)
    etEvTot.append(0)
    etstEvTotSampleSize.append(0)

allowAutoSplits=bool(False)
randomFromObservations=bool(False)

initialPathObservation=[]
initialPathObservationTimes=[]
observationNumber=0

minObs=3
indMinObs=0

while observationNumber<minObs:
    initialPathObservation=[]
    initialPathObservationTimes=[]
    observationNumber=0

    merged=False
    propArray=[int(0),int(0)]
    St=int(1)
    Et=int(1)
    PIPprop=[St]
    updateTime=int(0)

    for t1 in range(TS):
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



            

            Et=Et2
            St=St2
            etEvTot[t1]+=Et
            stEvTot[t1]+=St
            etstEvTotSampleSize[t1]+=1
            propArray.pop(0)
            propArray.append(0)

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
                if randomFromObservations==True:
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

                    updateTime+=TSupdateTime
                    St=St2
                    PIPprop.append(St)

                    for i2 in range(TSupdateTime):
                        PIPprop.pop(0)
                        initialPathObservation.append(float(Et))
                        initialPathObservationTimes.append(int(t1))
                        observationNumber+=1
                else:

                                       
                    rand=float(random.random())
                    if rand<ObservtionProbabilityFactor:
                        initialPathObservation.append(float(Et))
                        initialPathObservationTimes.append(int(t1))
                        observationNumber+=1


                

    


    


print("Hello1")
print(initialPathObservation)
print(initialPathObservationTimes)
print(observationNumber)
print(initialPathObservationTimes[observationNumber-1])


ObsStEvTot=[]
ObsEtEvTot=[]
ObsTauEvTot=[]
ObsSampleSize=[]

for i1 in range(initialPathObservationTimes[observationNumber-1]):
    ObsStEvTot.append(0)
    ObsEtEvTot.append(0)
    ObsTauEvTot.append(0)
    ObsSampleSize.append(0)

samples=[]
for i in range(observationNumber):
    samples.append(0)

inputChosen=bool(True)
probabilityFactor=20

for i2 in range(iterations):
    
    observationInd=0
    
    merged=False
    propArray=[int(0),int(0)]
    St=int(1)
    Et=int(1)
    PIPprop=[St]
    updateTime=int(0)

    for t1 in range(initialPathObservationTimes[observationNumber-1]):
        shift=[int(0)]

        if inputChosen==True:

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
                        ObsEtEvTot[t1]+=Et2
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

                if observationInd<observationNumber:
                
                    if t1==initialPathObservationTimes[observationInd]:
                        Et2=float(Et)          
                        if Et2<initialPathObservationTimes[observationInd]+(probabilityFactor)/2:
                            if Et2>initialPathObservationTimes[observationInd]-(probabilityFactor)/2:
                                inputChosen=True
                                samples[observationInd]+=1
                                observationInd+=1
                            else:
                                inputChosen=False
                                observationInd+=1
                        else:
                            inputChosen=False
                            observationInd+=1

                if inputChosen==False:
                    break

            if inputChosen==False:
                break
    
    inputChosen=True






normValObsStEvTot=[]
normValObsEtEvTot=[]
normValObsTauEvTot=[]

normStEvTot=[]
normEtEvTot=[]

for i1 in range(initialPathObservationTimes[observationNumber-1]):
    if ObsSampleSize[i1]>0:
        normValObsStEvTot.append(ObsStEvTot[i1]/ObsSampleSize[i1])
        normValObsEtEvTot.append(ObsEtEvTot[i1]/ObsSampleSize[i1])
        normValObsTauEvTot.append(ObsTauEvTot[i1]/ObsSampleSize[i1])
    else:
        normValObsStEvTot.append(0)
        normValObsEtEvTot.append(0)
        normValObsTauEvTot.append(0)

for i1 in range(initialPathObservationTimes[observationNumber-1]):
    if etstEvTotSampleSize[i1]>0:
        normStEvTot.append(stEvTot[i1]/etstEvTotSampleSize[i1])
        normEtEvTot.append(etEvTot[i1]/etstEvTotSampleSize[i1])
    else:
        normStEvTot.append(0)
        normEtEvTot.append(0)

y1 = np.array(normValObsStEvTot)
y2 = np.array(normValObsEtEvTot)
y3 = np.array(normValObsTauEvTot)

x1 = np.array(normStEvTot)
x2 = np.array(normEtEvTot)


initialPathObservationTimes

print(observationNumber)
print(initialPathObservation)
print(initialPathObservationTimes)
print(samples)

plt.plot(y1)
plt.plot(y2)
plt.plot(y3)
plt.plot(x1)
plt.plot(x2)
plt.plot(initialPathObservationTimes, initialPathObservation, 'ro')



plt.show()