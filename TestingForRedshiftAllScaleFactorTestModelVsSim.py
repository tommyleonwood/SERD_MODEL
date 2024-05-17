import random
import math
import matplotlib.pyplot as plt
import numpy as np

#Probabilities for duplications and reductions
pd=1/3
pr=1/3
pt=pd+pr

#Different choices for TSs and filament scaling to be run through iteratively
timeDiff=10
runs=20
times=[]
for r in range(runs):
    times.append((r+1)*timeDiff)


emmittedGammaValsSim=[]
emmittedGammaValsModel=[]

#Run through different choice of total TSs
for ts in times:
    "TS1 will count size of evolution to measurement and TS2 the size of the evolution after"
    TS1=ts

    "We may need to devise a method for signalling equlilbrium state of static system. This will allow us to scale the initial "
    "state of the system."
    "From there we may be able to implement both a growing and a static run."

    TS2=ts

    #Calculate a suitable number of loops/iterations for the first and second MC runs
    iterations1=400000*math.ceil(ts/10)
    iterations2=100000

    ScaleMin=ts-math.floor(2*ts/timeDiff)
    ScaleMax=ts+math.floor(2*ts/timeDiff)

    #Total SE in filament evolution arrays
    SEnEv=[]
    SEnEvAv=[]

    #Bool for if filament fully reduced and a counter to keep tract of the number of fully reduced filaments in the first tier
    merged=False
    fullyReducedN=0

    #Keeps track of samples outsidfe of scaling range
    notInRange=0

    #Value stores, output of one tier becomes the output of the second tier
    SEnInput=[]
    PIParrayInput=[]
    emissionTimeInput=[]
    SEnEvInput=[]

    #Keep tract of time lag and time difference between observation time and global time
    timeLagStore=[0]
    timeDifference=0

    #Cap for if the time difference is significant
    timeThreshold=50

    #Gap for measuring observed gamma values
    timeEndGap=10

    for i in range(TS1):
        timeLagStore.append(1)
    emmissionTime2=0

    #First MC run for scaling
    for i in range(iterations1):
        SEn=1
        PIParray1=[0]
        for s in range(SEn):
            PIParray1.append(0)
        emissionTime=0
        SEnEv=[]

        merged=False

        #Run through TSs
        for t in range(TS1):
            #Evolution algorithm
            if merged == False:
                shift=0
                for SEin in range(SEn):
                    rand=random.random()
                    if rand<pd:
                        PIParray1.insert(SEin+shift+1,0)
                        shift+=1
                    elif rand<pt:
                        PIParray1[SEin+shift]+=PIParray1[SEin+shift+1]
                        PIParray1.pop(SEin+shift+1)
                        shift-=1
                emissionTime+=PIParray1[0]
                timeLagStore[t]+=emissionTime

                PIParray1.pop(0)
                PIParray1.append(1)
                SEn+=shift
                SEnEv.append(SEn)

                

                if SEn==0:
                    merged=True
                    fullyReducedN+=1
                    break

        timeDifference=t-emissionTime

        #Save outputs of scaled system as sample of inputs for new MC run
        if merged==False:
            if SEn>ScaleMin:
                if SEn<ScaleMax:
                    "Tells us if the system is close to static equlibrium";
                    if timeDifference<timeThreshold:
                        SEnInput.append(SEn)
                        PIParrayInput.append(PIParray1)
                        emissionTimeInput.append(emissionTime)
                        SEnEvInput.append(SEnEv)
                    else:
                        notInRange+=1
                else:
                    notInRange+=1
            else:
                notInRange+=1

    #Caluclate new input sample size
    sampleSize1=iterations1-fullyReducedN-notInRange

    print(sampleSize1)

    PIPupdateDensityEv=[]
    PIPdensityTot=0

    EmmissionScaleTempAv=[]
    EmmissionScaleSample=[]
    for se in range(TS2):
        EmmissionScaleTempAv.append(0)
        EmmissionScaleSample.append(0)

    EOratioScaleEredTot=[]
    for se in range(TS2):
        EOratioScaleEredTot.append(0)

    fullyReducedN2=0

    PIPdensAvStatic=0
    PIPdensAvGrowing=0

    PIPtotStatic=0
    PIPtotGrowing=0

    emissionDiff=[]
    for t in range(TS2):
        emissionDiff.append(0)



    scaleFactorRatio=[]
    gammaObserved=[]
    gammaEmmitted=[]
    gammaEmmittedAv=0
    SampleSize=0


    if sampleSize1>0:
        #Start second MC loop
        for i in range(iterations2):
            randomIn=random.randint(0,sampleSize1-1)
            SEn=SEnInput[randomIn]

            PIParray2=[]
            for s in range(SEn+1):
                PIParray2.append(PIParrayInput[randomIn][s])

            SEnEv=[]
            for t in range(TS1):
                SEnEv.append(SEnEvInput[randomIn][t])

            emissionTime=emissionTimeInput[randomIn]
            emissionTime2=0

            PIPupdateDensityEv=[]
            ScaleAtEmmission=[]

            merged=False
            for t in range(TS2):
                #Evolution algorithm
                if merged == False:
                    shift=0
                    for SEin in range(SEn):
                        rand=random.random()
                        if rand<pd:
                            PIParray2.insert(SEin+shift+1,0)
                            shift+=1
                        elif rand<pt:
                            PIParray2[SEin+shift]+=PIParray2[SEin+shift+1]
                            PIParray2.pop(SEin+shift+1)
                            shift-=1
                    PIPupdateDensityEv.append(PIParray2[0])
                    PIParray2.pop(0)
                    PIParray2.append(1)
                    SEn+=shift
                    SEnEv.append(SEn)

                    if SEn==0:
                        merged=True
                        fullyReducedN2+=1


            emissionTime2=emissionTime

            #Store observed values for gamma, within a given time gap
            if merged==False:
                for t in range(TS2):
                    emissionTime2+=PIPupdateDensityEv[t]
                    if t>TS2-timeEndGap:
                        if PIPupdateDensityEv[t]!=0:
                            gammaObserved.append(PIPupdateDensityEv[t])
                            scaleFactorRatio.append(SEnEv[emissionTime2]/SEn)
                            gammaEmmitted.append(PIPupdateDensityEv[t]*SEn/SEnEv[emissionTime2])

                            #Store gamma emmitted from Sim outputs based on redshift scale factor relation
                            gammaEmmittedAv+=PIPupdateDensityEv[t]*SEn/SEnEv[emissionTime2]
                            SampleSize+=1


        #Calculate the value of gammaEmmitted by extracting the average gammaObserved for the scale factor within a given range r
        #This will correspond to the effective gammaEmmitted as it will be the value of gamma for which the system has not grown
        #the information has propagated through a static space.
        r=0.1
        gammaEmmittedSimAv=0
        sampleSize2=0
        for i in range(len(scaleFactorRatio)):
            if scaleFactorRatio[i]<1+r:
                if scaleFactorRatio[i]>1-r:
                    gammaEmmittedSimAv+=gammaObserved[i]
                    sampleSize2+=1

        gammaEmmittedSimAv/=sampleSize2
        print(sampleSize2)
        gammaRatio=[]
        for i in range(len(gammaObserved)):
            gammaRatio.append(gammaObserved[i]/gammaEmmittedSimAv)

        gammaEmmittedAv/=SampleSize

        print(gammaEmmittedAv)
        print(gammaEmmittedSimAv)
        print("End run with TS = "+str(ts))
        emmittedGammaValsSim.append(gammaEmmittedSimAv)
        emmittedGammaValsModel.append(gammaEmmittedAv)
    else:
        print("sample Size 0")


    
#Plot the simulated emmitted gamma vs the emmitted gamma extracted from using the redshift scale factor relation (PIPupdateDensityEv[t]*SEn/SEnEv[emissionTime2])
plt.plot(emmittedGammaValsSim)
plt.plot(emmittedGammaValsModel)
plt.grid(True)
plt.show()
