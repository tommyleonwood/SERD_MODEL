import random

pd=float(1/3)
pr=float(1/3)
pt=pd+pr

m=5;

propArray=[int(0),int(0)]
St=int(1)
Et=int(1)

StIn=[]
EtIn=[]

IEpropIn=[]

StEv=[]
EtEv=[]

IEpropEv=[]

emptyArray=[]
for i1 in range(m):
    emptyArray.append(0)

#In this section we construct the global instananious state matrix St and the information
#propagating tensor IEprop (this has both a left handed and a right handed component).
for i2 in range(m):
    IEproplist1=[]
    for i1 in range(m):
        stList1=[]
        IEproplist2=[]
        if(i1!=i2):
            stList1=[]
            IEproplist3=[]
            stList1.append(1)
            for i3 in range(m):
                IEproplist3.append(0)
            #This makes each entry a list of lists. At each node (IG) in the state graph an array exists. In the constant m case
            #the arrays can superpose on one another and add since they all have the same rank and so the operation of superposition
            #is simply addition. If m was not constant the rank of the arrays would be different and so concatination of arrays of arrays
            #would be the operation. Since all IE only have 1 SE there are only two arrays one for each enclosing IG.
            IEproplist2.append(IEproplist3)
            IEproplist2.append(IEproplist3)
        else:
            stList1.append(0)
            IEproplist2.append(0)
        IEproplist1.append(IEproplist2)
    IEpropIn.append(IEproplist1)

    StIn.append(stList1)

#In this section we construct the initial observed state tensor EtEvIn
for i1 in range(m):
    Etlist1=[]
    for i2 in range(m):
        Etlist2=[]
        for i3 in range(m):
            if i2!=i3:
                Etlist2.append(1)
            else:
                Etlist2.append(0)
        Etlist1.append(Etlist2)
    EtIn.append(Etlist1)

StEv.append(StIn)
print(StEv)
IEpropEv.append(IEpropIn)

StEvIt=[]
EtEvIt=[]

TS=int(25)
iterations=int(15)

#Run through TSs
for i1 in range(TS):
    IEprop=IEpropEv[i1]
    St=StEv[i1]

    #Run through IEs
    for i2 in range(m):
        for i3 in range(m):
            #only need to run through IE in one direction
            if i2<i3:
                actions=[]
                sadd=0
                s=St[i2][i3]

                for i4 in range(s):
                #generate random number to assign element action outcome
                    rand=float(random.random())

                    #assign outcome for each element
                    if rand<pd:
                        actions.append(1)
                        sadd+=1
                    else:
                        if rand<pt:
                            actions.append(-1)
                            sadd-=1
                        else:
                            actions.append(0)

                s2=s+sadd
                St[i2][i3]=s2
                shift=0

                for i4 in range(s):
                    if actions[i4]==1:

                        IEprop[i1][i2][i3][i4+shift][i2]+=1
                        IEprop[i1][i2][i3].insert(i4+shift,emptyArray)

                        IEprop[i1][i3][i2][i4+shift+1][i3]+=1
                        IEprop[i1][i3][i2].insert(i4+shift,emptyArray)

                        shift+=actions[i4]

                    if actions [i4]==-1:

                        IEprop[i1][i2][i3][i4+shift]+=IEprop[i1][i2][i3][i4+shift+1]
                        IEprop[i1][i2][i3][i4+shift][i2]-=1
                        IEprop[i1][i3][i2][i4+shift]+=IEprop[i1][i2][i3][i4+shift+1]
                        IEprop[i1][i3][i2][i4+shift][i3]-=1

                        shift-=actions[i4]

    s=s2
    StEv.append(St)
    IEpropEv.append(IEprop)


                










                















