import random
from re import X

pd=float(1/3)
pr=float(1/3)
pt=pd+pr

m=3;

StIn=[]
EtIn=[]

StEv=[]
EtEv=[]

nonRadIEpropIn=[]
radIEpropIn=[]

#define initial radial IE states
#from each non-central PP
for i1 in range (m-1):
    radIEprop1=[]

    radIEprop1.append([0 for x in range(m)])
    radIEprop1.append([0 for x in range(m)])

    radIEpropIn.append(radIEprop1)



#define initial nonradial IE states
for i1 in range(m-1):
    nonRadIEprop1=[]
    for i2 in range(m-1):
        if i1!=i2:
            nonRadIEprop1.append([0,0])
        else:
            nonRadIEprop1.append(0)
    nonRadIEpropIn.append(nonRadIEprop1)



#define initial Et and St matrices
for i1 in range(m):
    Sta=[]
    Eta=[]
    for i2 in range(m):
        if i1!=i2:

            if i1>i2:
                Sta.append(1)
            else:
                Sta.append(0)
            if i1>0:
                Eta.append(1)
            else:
                Eta.append(0)
        else:
            Sta.append(0)
            Eta.append(0)

    StIn.append(Sta)
    EtIn.append(Eta)

print(EtIn)

StEv=[]
EtEv=[]
StEv.append(StIn)
EtEv.append(EtIn)

St1=StIn
Et1=EtIn



nonRadIEprop2=nonRadIEpropIn
radIEprop2=radIEpropIn

nonRadIEprop=[]
radIEprop=[]

nonRadIEprop.append(nonRadIEpropIn)
radIEprop.append(radIEpropIn)



StEvIt=[]
EtEvIt=[]

TS=int(10)

#run through TSs 
for i1 in range(TS):
    print("TS=")
    print(i1)

    print("a")

    print(St1)
    print(Et1)
    print(nonRadIEprop2)
    print(radIEprop2)

    #run through radial IEs
    for i2 in range(m-1):
        shift=[0]
        Stn=St1[i2+1][0]
        #run through SEs within radial IEs
        for i3 in range(Stn):
            rand=float(random.random())
            #duplicates
            if rand<pd:
                radIEprop2[i2][i3+shift[i3]][0]+=1    
                radIEprop2[i2].insert(i3+shift[i3]+1,[0 for x in range(m)])
                shift.append(shift[i3]+1) 
            #reduces
            elif rand<pt:
                radIEprop2[i2][i3+shift[i3]]=[radIEprop2[i2][i3+shift[i3]][i4]+radIEprop2[i2][i3+shift[i3]+1][i4] for i4 in range(m)]
                radIEprop2[i2][i3+shift[i3]][0]-=1
                radIEprop2[i2].pop(i3+shift[i3]+1)
                shift.append(shift[i3]-1)
            #remains
            else:
                shift.append(shift[i3])

        St1[i2+1][0]+=shift[Stn]

        #auto undo mergers with info conservation
        if St1[i2+1][0]==0:
            St1[i2+1][0]=1

            #info cons

            Et1[i2+1]=[Et1[i2+1][i3]+radIEprop2[i2][0][i3] for i3 in range(m)]

            radIEprop2[i2]=[[0 for x in range(m)],[0 for x in range(m)]]
            Et1[i2+1][0]=1

    print("b")
    print(St1)
    print(Et1)
    print(nonRadIEprop2)
    print(radIEprop2)




    #run through non-radial IE
    for i2 in range(m-1):
        for i3 in range (m-1):
            if i2>i3:
                Stn=St1[i2+1][i3+1]
                shift=[0]
                for i4 in range(Stn):
                    rand=float(random.random())
                    #duplicates
                    if rand<pd:
                        nonRadIEprop2[i2][i3].insert(i4+shift[i4]+1,0)
                        nonRadIEprop2[i2][i3][i4+shift[i4]]+=1

                        nonRadIEprop2[i3][i2][i4+shift[i4]+1]+=1
                        nonRadIEprop2[i3][i2].insert(i4+shift[i4]+1,0)

                        shift.append(shift[i4]+1)
                    else:
                        #reduces
                        if rand<pt:
                            nonRadIEprop2[i2][i3][i4+shift[i4]]+=(nonRadIEprop2[i2][i3][i4+shift[i4]+1]-1)
                            nonRadIEprop2[i2][i3].pop(i4+shift[i4]+1)

                            nonRadIEprop2[i3][i2][i4+shift[i4]+1]+=(nonRadIEprop2[i3][i2][i4+shift[i4]]-1)
                            nonRadIEprop2[i3][i2].pop(i4+shift[i4])

                            shift.append(shift[i4]-1)
                        else:
                            shift.append(shift[i4])

                St1[i2+1][i3+1]+=shift[Stn]

                #auto undo mergers with info conservation
                if St1[i2+1][i3+1]==0:
                    St1[i2+1][i3+1]=1

                    radIEprop2[i2][St1[i2+1][0]][i3+1]+=(nonRadIEprop2[i2][i3][0]+1)
                    radIEprop2[i3][St1[i3+1][0]][i2+1]+=(nonRadIEprop2[i3][i2][0]+1)

                    nonRadIEprop2[i2][i3]=[0,0]
                    nonRadIEprop2[i3][i2]=[0,0]


    print("c")
    print(St1)
    print(Et1)
    print(nonRadIEprop2)
    print(radIEprop2)

    #propagate and transfer information

    #propagate radially and update Et
    for i2 in range(m-1):


        Et1[i2+1]=[Et1[i2+1][i3]+radIEprop2[i2][0][i3] for i3 in range(m)]

        radIEprop2[i2].pop(0)
        radIEprop2[i2].append([0 for x in range(m)])

    print("d")
    print(St1)
    print(Et1)
    print(nonRadIEprop2)
    print(radIEprop2)

    #propagate non-radially, transfer information across PP and update radial IE
    for i2 in range(m-1):
        for i3 in range(m-1):
            if i2>i3:
                radIEprop2[i2][St1[i2+1][0]][i3+1]+=nonRadIEprop2[i2][i3][0]
                nonRadIEprop2[i2][i3].pop(0)
                nonRadIEprop2[i2][i3].append(0)

                radIEprop2[i3][St1[i3+1][0]][i2+1]+=nonRadIEprop2[i3][i2][St1[i2+1][i3+1]]
                nonRadIEprop2[i3][i2].pop(St1[i2+1][i3+1])
                nonRadIEprop2[i3][i2].insert(0,0)

    StEv.append(St1)
    EtEv.append(Et1)

    radIEprop.append(radIEprop2)
    nonRadIEprop.append(nonRadIEprop2)


    print("e")
    print(St1)
    print(Et1)
    print(nonRadIEprop2)
    print(radIEprop2)

    #test for information conservation and correctness
    Ssum=0
    Esum=0
    Asum=0

    for i2 in range(m):
        for i3 in range (m):
            if i2>i3:
                Ssum+=St1[i2][i3]
                if i3>0:
                    Ssum+=St1[i2][i3]
                    
            Esum+=Et1[i2][i3]






    for i2 in range(m-1):
        for i3 in range (m-1):
            if i2>i3:
                for i4 in range (St1[i2+1][i3+1]+1):
                    Asum+=nonRadIEprop2[i2][i3][i4]
                    Asum+=nonRadIEprop2[i3][i2][i4]

    for i2 in range(m-1):
        for i3 in range (St1[i2+1][0]+1):
            for i4 in range(m):
                Asum+=radIEprop2[i2][i3][i4]

    infoTot=Ssum-Asum-Esum
    print(Ssum)
    print(Esum)
    print(Asum)
    print(infoTot)
    if infoTot!=0:
        print("Error")
    else:
        print("success")




x=4
y=0
y+=x 
print(y)


            