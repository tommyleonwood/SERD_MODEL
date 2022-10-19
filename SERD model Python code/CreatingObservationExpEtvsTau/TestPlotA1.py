import matplotlib.pyplot as plt
import numpy as np
import random

def FromUnionStateTauVsEtStInputObs(pd,pt,TS1,random,obsStLoc,obsStHalfRange,obsEtLoc,obsEtHalfRange,ObsStEvTot,ObsEtEvTot,ObsTauEvTot,ObsSampleSize,propArray,PIPprop,St,Et,updateTime):

    
    merged=False

    propArray=[int(0),int(0)]
    St=int(1)
    Et=int(1)

    PIPprop=[St]
    updateTime=int(0)

    for t1 in range(TS1):
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
                    ObsTauEvTot[t1]+=TS1+t1-updateTime-i2+1
                    ObsSampleSize[t1]+=1

                updateTime+=TSupdateTime
                St=St2
                PIPprop.append(St)

                Et=Et2
                propArray.pop(0)
                propArray.append(0)




def FromInputStateTauVsEtStInputObs(pd,pt,TS1,TS2,random,obsStLoc,obsStHalfRange,obsEtLoc,obsEtHalfRange,allowAutoSplits,ObsStEvTot,ObsEtEvTot,ObsTauEvTot,ObsSampleSize,propArray,PIPprop,St,Et,updateTime):
    for t1 in range(TS2):
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
                for t2 in range(TS1+t1-updateTime):
                    if PIPprop[t2]<1:
                        TSupdateTime+=1
                    else:
                        break

                for i2 in range(TSupdateTime):
                    PIPprop.pop(0)
                    ObsStEvTot[TS1+t1]+=St2
                    ObsEtEvTot[TS1+t1]+=Et
                    ObsTauEvTot[TS1+t1]+=TS1+t1-updateTime-i2+1
                    ObsSampleSize[TS1+t1]+=1

                updateTime+=TSupdateTime
                St=St2
                PIPprop.append(St)

                Et=Et2
                propArray.pop(0)
                propArray.append(0)




pd=1/3
pr=1/3
pt=pr+pd

measurements=int(3)

iterations=int(100000)

TS=int(50)

TS1=[]
for i1 in range(measurements+1):
    TS1.append(TS)

sampleSize=[]
for i1 in range(measurements+1):
    sampleSize.append(0)

obsStLoc=[]
obsStHalfRange=[]
obsEtLoc=[]
obsEtHalfRange=[]

for i1 in range(measurements+1):
    obsStLoc.append(2*TS*(0.5*(i1+1)))
    obsStHalfRange.append(2*TS/2*(0.5*(i1+1)))
    obsEtLoc.append(TS*(0.5*(i1+1)))
    obsEtHalfRange.append(TS/2*(0.5*(i1+1)))

print(obsStLoc)
print(obsStHalfRange)
print(obsEtLoc)
print(obsEtHalfRange)

ObsStEvTot=[]
ObsEtEvTot=[]
ObsTauEvTot=[]
ObsSampleSize=[]

for i1 in range(TS*(measurements+1)):
    ObsStEvTot.append(0)
    ObsEtEvTot.append(0)
    ObsTauEvTot.append(0)
    ObsSampleSize.append(0)

inputChosen=bool(False)
propArray=[int(0),int(0)]
St=int(1)
Et=int(1)
PIPprop=[St]
updateTime=int(0)


for i2 in range(iterations):

    FromUnionStateTauVsEtStInputObs(pd,pt,TS,random,obsStLoc[0],obsStHalfRange[0],obsEtLoc[0],obsEtHalfRange[0],ObsStEvTot,ObsEtEvTot,ObsTauEvTot,ObsSampleSize,propArray,PIPprop,St,Et,updateTime)
    if St>obsStLoc[0]-obsStHalfRange[0]:
        if St<=obsStLoc[0]+obsStHalfRange[0]:
            if Et>obsEtLoc[0]-obsEtHalfRange[0]:
                if Et<=obsEtLoc[0]+obsEtHalfRange[0]:
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
                FromInputStateTauVsEtStInputObs(pd,pt,(1+i1)*TS,TS,random,obsStLoc[i1+1],obsStHalfRange[i1+1],obsEtLoc[i1+1],obsEtHalfRange[i1+1],allowAutoSplits,ObsStEvTot,ObsEtEvTot,ObsTauEvTot,ObsSampleSize,propArray,PIPprop,St,Et,updateTime)
                if St>obsStLoc[i1+1]-obsStHalfRange[i1+1]:
                    if St<=obsStLoc[i1+1]+obsStHalfRange[i1+1]:
                        if Et>obsEtLoc[i1+1]-obsEtHalfRange[i1+1]:
                            if Et<=obsEtLoc[i1+1]+obsEtHalfRange[i1+1]:
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

for i1 in range(TS*(measurements+1)):
    if ObsSampleSize[i1]>0:
        normValObsStEvTot.append(ObsStEvTot[i1]/ObsSampleSize[i1])
        normValObsEtEvTot.append(ObsEtEvTot[i1]/ObsSampleSize[i1])
        normValObsTauEvTot.append(ObsTauEvTot[i1]/ObsSampleSize[i1])
    else:
        normValObsStEvTot.append(0)
        normValObsEtEvTot.append(0)
        normValObsTauEvTot.append(0)

print(normValObsStEvTot)
print(normValObsEtEvTot)
print(normValObsTauEvTot)

print(sampleSize)

y1 = np.array(normValObsStEvTot)
y2 = np.array(normValObsEtEvTot)
y3 = np.array(normValObsTauEvTot)

plt.plot(y1)
plt.plot(y2)
plt.plot(y3)
plt.show()