pd=float(1/3)
pr=float(1/3)
pt=pd+pr

propArray=[int(0),int(0)]
St=int(1)
Et=int(1)

StEv=[St]
EtEv=[Et]

# Creating an array
StEvIt=[]
EtEvIt=[]

TS=int(25)
iterations=int(15)



allowSplits=bool(True)
merged=bool(False)

for i1 in range(TS+1):
   StEvIt.append(0)
   EtEvIt.append(0)

import random

for i1 in range(iterations):
    propArray=[int(0),int(0)]
    St=int(1)
    Et=int(1)

    PIPprop=[St]
    updateTime=int(0)

    StEv=[St]
    EtEv=[Et]

    StEvIt[0]+=St
    EtEvIt[0]+=Et

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


    ind1=int(0)

    for x2 in range(St+1):
        while PIPprop[ind1]==x2:
            PIPprop[ind1]+=shift[x2]-1
            ind1+=1
            if ind1>St:
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

    St+=shift[St]

    PIPprop.append(St)

    if St==0:
        if allowSplits==False:
            merged=True

        else:
            St=1
            Et=1
            propArray=[int(0),int(0)]
            shift=[int(0)]
            PIPprop=[St]

    Et+=propArray[0]
    propArray.pop(0)
    propArray.append(0)

    StEv.append(St)
    EtEv.append(Et)
  
    StEvIt[t1+1]+=St
    EtEvIt[t1+1]+=Et


print(EtEvIt)
