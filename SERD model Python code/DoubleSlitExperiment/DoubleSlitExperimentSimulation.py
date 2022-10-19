import matplotlib.pyplot as plt
import numpy as np
import random

pd=1/3
pr=1/3
pt=pr+pd

measurements=int(2)

iterations=int(10000000)

TS=int(1000)

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

obsStLoc2=[]
obsStHalfRange2=[]

obsStMin1=[]
obsStMax1=[]

obsStMin2=[]
obsStMax2=[]


obsEtMin1=[]
obsEtMax1=[]

obsEtMin2=[]
obsEtMax2=[]

StAndEtMeasure=bool(False)

StEtObservationRatio=float(1)
timeScaleFactor=float(0.5)
distanceBetweenSlits=float(20)
slitWidth=float(10)

scaleLimit=sumTS
probabilityDistributionEvA1=[]

for i1 in range(sumTS*scaleLimit):
    probabilityDistributionEvA1.append(int(0));

print(sumTS)
print(scaleLimit)

for i1 in range(measurements+1):
    if(StAndEtMeasure==False):
        obsEtMin1.append(timeScaleFactor*countTS[i1+1]-(distanceBetweenSlits/2)-slitWidth/2)
        obsEtMax1.append(timeScaleFactor*countTS[i1+1]-(distanceBetweenSlits/2)+slitWidth/2)
        obsEtMin2.append(timeScaleFactor*countTS[i1+1]+(distanceBetweenSlits/2)-slitWidth/2)
        obsEtMax2.append(timeScaleFactor*countTS[i1+1]+(distanceBetweenSlits/2)+slitWidth/2)

obsEtMin1[0]=980
obsEtMax1[0]=1020
obsEtMin2[0]=980
obsEtMax2[0]=1020

obsEtMin1[1]=960
obsEtMax1[1]=980

obsEtMin2[1]=1020
obsEtMax2[1]=1040



print(obsEtMin1)
print(obsEtMax1)
print(obsEtMin2)
print(obsEtMax2)

ObsStEvTot=[]
ObsEtEvTot=[]
ObsTauEvTot=[]
ObsSampleSize=[]

for i1 in range(sumTS):
    ObsStEvTot.append(0)
    ObsEtEvTot.append(0)
    ObsTauEvTot.append(0)
    ObsSampleSize.append(0)

sampleSize=[]

for i1 in range(sumTS):
    sampleSize.append(int(0))

inputs=int(0)

etEvTot=[]

for i2 in range(iterations):
    
    
    
    merged=False
    propArray=[int(0),int(0)]
    St=int(1)
    Et=int(1)

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
                merged=True
            else:
                St=St2
                Et=Et2
                propArray.pop(0)
                propArray.append(0)

                if Et<sumTS:
                    probabilityDistributionEvA1[t1*scaleLimit+Et]+=1
                    sampleSize[t1]+=1



    if StAndEtMeasure==False:
        if Et>obsEtMin1[0]:
            if Et<obsEtMax1[0]:
                inputChosen=True
                inputs+=1
            
            else:
                if Et>obsEtMin2[0]:
                    if Et<obsEtMax2[0]:
                        inputChosen=True
                        inputs+=1
                 
                    else:
                        inputChosen=False
                else:
                    inputChosen=False
        else:
            inputChosen=False


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
                        else:




                            St=St2
                            Et=Et2
                            propArray.pop(0)
                            propArray.append(0)

                        if Et<sumTS:
                            probabilityDistributionEvA1[(countTS[i1+1]+t1-1)*scaleLimit+Et]+=1
                            sampleSize[countTS[i1+1]+t1-1]+=1
                            
                if StAndEtMeasure==False:
                        if Et>obsEtMin1[i1+1]:
                            if Et<obsEtMax1[i1+1]:
                                inputChosen=True
                            else:
                                if Et>obsEtMin2[i1+1]:
                                    if Et<obsEtMax2[i1+1]:
                                        inputChosen=True
                                    else:
                                        inputChosen=False
                                else:
                                    inputChosen=False
                        else:
                            inputChosen=False
            else:
                break


probabilityDistributionEvA1str=str(probabilityDistributionEvA1)

g = open("NewTextFile.txt", "w")
g.write(probabilityDistributionEvA1str)
g.close()

print(sumTS)

sampleSizeStr="{"
for i1 in range(sumTS):
    sampleSizeStr=sampleSizeStr+str(sampleSize[i1])+","
    if i1==sumTS-1:
        sampleSizeStr=sampleSizeStr+str(sampleSize[i1])+"}"

g2 = open("sampleSize.txt", "w")
g2.write(sampleSizeStr)
g2.close()
print(inputs)
print(countTS)
print(sampleSize)


