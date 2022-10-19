def FromUnionStateTauVsEtStInputObs(pd,pt,TS1,random,obsStLoc,obsStHalfRange,obsEtLoc,obsEtHalfRange,ObsStEvTot,ObsEtEvTot,ObsTauEvTot,ObsSampleSize):

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

                updateTime+=TSupdateTime
                St=St2
                PIPprop.append(St)

                Et=Et2
                propArray.pop(0)
                propArray.append(0)

            
                if St>obsStLoc-obsStHalfRange:
                    if St<=obsStLoc+obsStHalfRange:
                        if Et>obsEtLoc-obsEtHalfRange:
                            if Et<=obsEtLoc+obsEtHalfRange:
                                inputChosen=bool(True)
                                sampleSize+=1






def FromInputStateTauVsEtStInputObs(pd,pt,TS1,TS2,random,obsStLoc,obsStHalfRange,obsEtLoc,obsEtHalfRange,allowAutoSplits,ObsStEvTot,ObsEtEvTot,ObsTauEvTot,ObsSampleSize,normValObsStEvTot,normValObsEtEvTot,normValObsTauEvTot):
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


            for i1 in range(TS2):
                if ObsSampleSize[TS1+i1]>0:
                    normValObsStEvTot.append(ObsStEvTot[i1]/ObsSampleSize[i1])
                    normValObsEtEvTot.append(ObsEtEvTot[i1]/ObsSampleSize[i1])
                    normValObsTauEvTot.append(ObsTauEvTot[i1]/ObsSampleSize[i1])
            else:
                normValObsStEvTot.append(0)
                normValObsEtEvTot.append(0)
                normValObsTauEvTot.append(0)
