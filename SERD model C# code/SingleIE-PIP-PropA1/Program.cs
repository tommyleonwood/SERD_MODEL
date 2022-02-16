using System;
using System.Collections;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace serd
{
    class program
    {


//The following methods give Et as output from unity initial state
        static int[] totEtEvItWASplit(double pd, double pr, int TS, int iterations1)
        {

            int[] totetEv = new int[TS];
            double pt;
            pt=pd+pr;
            int St;
            int et;

            bool merged;
            var random = new Random();
            for(int i3=0; i3 < iterations1; i3++)
            {
            int[] IEProp = new int[2];
            IEProp[0] = 0;
            IEProp[1] = 0;
            merged = false;
            St = 1;
            et = St;
                for(int i1=0; i1<TS; i1++)
                {
                    int shift = 0;
                    int StDiff=0;
                    int[] actionsArray = new int[St];
                    if(merged==false)
                    {
                        for(int i2=0;i2<St;i2++)
                        {
                            double rand = random.NextDouble();
                            if(rand<pd)
                            {
                                actionsArray[i2] = 1;
                                StDiff+=1;
                            }
                            else
                            {
                                if(rand<pt)
                                {
                                    actionsArray[i2] = -1;
                                    StDiff-=1;
                                }
                                else
                                {
                                    actionsArray[i2] = 0;
                                }
                            }
                        }
                        int St2=St+StDiff;
                        if(St2==0){merged=true;}
                        if(merged==false)
                        {
                        int[] IEProp2 = new int[St2+1];
                        for(int i2=0;i2<St;i2++)
                        {
                            if(actionsArray[i2]==1)
                            {
                                IEProp2[i2+shift]=IEProp[i2]+1;
                                IEProp2[i2+shift+1]=0;
                                shift+=1;
                            }
                            else
                            {
                                if(actionsArray[i2]==-1)
                                {
                                    IEProp2[i2+shift]+=IEProp[i2+1]-1;
                                    shift-=1;
                                }
                                else
                                {
                                    IEProp2[i2+shift]=IEProp[i2];
                                }
                            }
                        }
                                
                        et+=IEProp2[0];
                        int[] IEProp3 = new int[St2+1];

                        for(int i2=0;i2<St2;i2++){IEProp3[i2]=IEProp2[i2+1];}
                        IEProp3[St2]=0;
                        IEProp=IEProp3;
                        St=St2;

                        }
                        else
                        {
                            St=1;
                            et=St;
                            int[] IEProp4 = new int[2];
                        IEProp4[0] = 0;
                        IEProp4[1] = 0;
                        IEProp=IEProp4;
                        merged = false;
                        }



                        if(i3==0)
                        {
                        totetEv[i1]=et;
                        }
                        else
                        {
                        totetEv[i1]+=et;
                        }
                        
                    }
                }
            }
        return totetEv;
        }

        static int[] totEtEvItnoSplit(double pd, double pr, int TS, int iterations1)
        {
            int[] totetEv = new int[TS];
            for(int i1=0;i1<TS;i1++)
            {
                totetEv[i1]=0;
            }
            double pt;
            pt=pd+pr;
            int St;
            int et;

            bool merged;
            var random = new Random();
            for(int i3=0; i3 < iterations1; i3++)
            {
            int[] IEProp = new int[2];
            IEProp[0] = 0;
            IEProp[1] = 0;
            merged = false;
            St = 1;
            et = St;
                for(int i1=0; i1<TS; i1++)
                {
                    int shift = 0;
                    int StDiff=0;
                    int[] actionsArray = new int[St];
                    if(merged==false)
                    {
                        for(int i2=0;i2<St;i2++)
                        {
                            double rand = random.NextDouble();
                            if(rand<pd)
                            {
                                actionsArray[i2] = 1;
                                StDiff+=1;
                            }
                            else
                            {
                                if(rand<pt)
                                {
                                    actionsArray[i2] = -1;
                                    StDiff-=1;
                                }
                                else
                                {
                                    actionsArray[i2] = 0;
                                }
                            }
                        }
                        int St2=St+StDiff;
                        if(St2==0){merged=true;}
                        if(merged==false)
                        {
                        int[] IEProp2 = new int[St2+1];
                        for(int i2=0;i2<St;i2++)
                        {
                            if(actionsArray[i2]==1)
                            {
                                IEProp2[i2+shift]=IEProp[i2]+1;
                                IEProp2[i2+shift+1]=0;
                                shift+=1;
                            }
                            else
                            {
                                if(actionsArray[i2]==-1)
                                {
                                    IEProp2[i2+shift]+=IEProp[i2+1]-1;
                                    shift-=1;
                                }
                                else
                                {
                                    IEProp2[i2+shift]=IEProp[i2];
                                }
                            }
                        }
                                
                        et+=IEProp2[0];
                        int[] IEProp3 = new int[St2+1];

                        for(int i2=0;i2<St2;i2++){IEProp3[i2]=IEProp2[i2+1];}
                        IEProp3[St2]=0;
                        IEProp=IEProp3;
                        St=St2;

                        }
                        else
                        {
                            break;
                        }

                        totetEv[i1]+=et;
                        
                    }
                }
            }
        return totetEv;
        }

//The following method gives back an ensemble that is averaged over the size of the ensamble plus the number of decayed IEs.

static int[] totEtEvItnoSplitTayloredEnsembleSize(double pd, double pr, int TS, int iterations1)
        {
            int[] outputData = new int [1+2+2*TS];
            int[] totstEv = new int[TS+1];
            int[] totetEv = new int[TS+1];
            for(int i1=0;i1<TS+1;i1++)
            {
                totetEv[i1]=0;
                totstEv[i1]=0;
            }
            double pt;
            pt=pd+pr;
            int St;
            int et;

            bool merged;
            var random = new Random();
            for(int i3=0; i3 < iterations1; i3++)
            {
            int[] IEProp = new int[2];
            IEProp[0] = 0;
            IEProp[1] = 0;
            merged = false;
            St = 1;
            et = St;

            int[] etEv=new int[TS];
            int[] StEv=new int[TS];


                for(int i1=0; i1<TS; i1++)
                {
                    int shift = 0;
                    int StDiff=0;
                    int[] actionsArray = new int[St];
                    if(merged==false)
                    {
                        for(int i2=0;i2<St;i2++)
                        {
                            double rand = random.NextDouble();
                            if(rand<pd)
                            {
                                actionsArray[i2] = 1;
                                StDiff+=1;
                            }
                            else
                            {
                                if(rand<pt)
                                {
                                    actionsArray[i2] = -1;
                                    StDiff-=1;
                                }
                                else
                                {
                                    actionsArray[i2] = 0;
                                }
                            }
                        }
                        int St2=St+StDiff;
                        if(St2==0){merged=true;}
                        if(merged==false)
                        {
                        int[] IEProp2 = new int[St2+1];
                        for(int i2=0;i2<St;i2++)
                        {
                            if(actionsArray[i2]==1)
                            {
                                IEProp2[i2+shift]=IEProp[i2]+1;
                                IEProp2[i2+shift+1]=0;
                                shift+=1;
                            }
                            else
                            {
                                if(actionsArray[i2]==-1)
                                {
                                    IEProp2[i2+shift]+=IEProp[i2+1]-1;
                                    shift-=1;
                                }
                                else
                                {
                                    IEProp2[i2+shift]=IEProp[i2];
                                }
                            }
                        }
                                
                        et+=IEProp2[0];
                        int[] IEProp3 = new int[St2+1];

                        for(int i2=0;i2<St2;i2++){IEProp3[i2]=IEProp2[i2+1];}
                        IEProp3[St2]=0;
                        IEProp=IEProp3;
                        St=St2;

                        }
                        else
                        {
                            break;
                        }


                        etEv[i1]=et;
                        StEv[i1]=St;


                        
                    }






                }

                if(merged==false)
                {
                    for(int i1=0;i1<TS;i1++){totetEv[i1+1]+=etEv[i1];}
                    for(int i1=0;i1<TS;i1++){totetEv[i1+1]+=etEv[i1];}
                    totetEv[0]+=1;
                    totstEv[0]+=1;
                }

                



            }
            


        return totetEv;
        }

//The following  methods output the results for the paper av Et and av Tau from unitary states

        static int[] totEtandTotTauMethodwithSplit(double pd, double pr, int TS, int iterations1)
        {
            int[] totetEv = new int[TS];
            int[] totStEv = new int[TS];
            int[] totTau = new int[TS];
            int[] sampleSize = new int[TS];

            int[] dataOutPut = new int[4*TS];

            for(int i1=0;i1<TS;i1++)
            {
                totetEv[i1]=0;
                totStEv[i1]=0;
                totTau[i1]=0;
                sampleSize[i1]=0;
            }
            
            double pt;
            pt=pd+pr;
            int St;
            int et;

            bool merged;
            var random = new Random();
            for(int i3=0; i3 < iterations1; i3++)
            {
            int[] IEProp = new int[2];
            int[] PIPprop = new int[TS];

            int updateTime = 0;

            IEProp[0] = 0;
            IEProp[1] = 0;
            merged = false;
            St = 1;
            et = St;


                for(int i1=0; i1<TS; i1++)
                {
                    int shift = 0;
                    int StDiff=0;
                    int[] actionsArray = new int[St];
                    int[] shiftArray = new int[St];
                    shiftArray[0]=0;

                    PIPprop[i1]=St;



                    if(merged==false)
                    {
                        for(int i2=0;i2<St;i2++)
                        {
                            double rand = random.NextDouble();
                            if(rand<pd)
                            {
                                actionsArray[i2] = 1;
                                
                                StDiff+=1;
                            }
                            else
                            {
                                if(rand<pt)
                                {
                                    actionsArray[i2] = -1;
                                    StDiff-=1;
                                }
                                else
                                {
                                    actionsArray[i2] = 0;
                                }
                            }
                        }
                        int St2=St+StDiff;
                        if(St2==0){merged=true;}
                        if(merged==false)
                        {
                        int[] IEProp2 = new int[St2+1];
                        for(int i2=0;i2<St;i2++)
                        {
                            if(actionsArray[i2]==1)
                            {
                                IEProp2[i2+shift]=IEProp[i2]+1;
                                IEProp2[i2+shift+1]=0;
                                shift+=1;
                            }
                            else
                            {
                                if(actionsArray[i2]==-1)
                                {
                                    IEProp2[i2+shift]+=IEProp[i2+1]-1;
                                    shift-=1;
                                }
                                else
                                {
                                    IEProp2[i2+shift]=IEProp[i2];
                                }
                            }
                        }

                        et+=IEProp2[0];
                        int[] IEProp3 = new int[St2+1];

                        for(int i2=0;i2<St2;i2++){IEProp3[i2]=IEProp2[i2+1];}
                        et+=IEProp3[0];
                        IEProp3[0]=0;
                        IEProp3[St2]=0;

                        //The following code runs for shifting the fully traversing PIPs.




                        

                        for(int i2=0;i2<St-1;i2++)
                        {
                            shiftArray[i2+1]=(shiftArray[i2]+actionsArray[i2]);
                        }
                        int ind1=0;
                        int count1=0;
                        for(int i2=0;i2<St;i2++)
                        {
                            if(PIPprop[updateTime+ind1]==i2)
                            {
                                count1=0;
                                while(PIPprop[updateTime+ind1+count1]==i2)
                                {
                                    PIPprop[updateTime+ind1+count1]+=shiftArray[i2];
                                    count1++;
                                    if(ind1+count1>i2-updateTime){break;}
                                    ind1+=count1;
                                }
                            }
                        }

                        int updateTimeInd=0;

                        for(int i2=0;i2<i1-updateTime+1;i2++)
                        {
                            PIPprop[updateTime+i2]-=1;
                            if(PIPprop[updateTime+i2]<1)
                            {
                                updateTimeInd+=1;
                                totetEv[i1]+=et;
                                totStEv[i1]+=St2;
                                totTau[i1]+=(i1-updateTime-i2);
                                sampleSize[i1]+=1;
                            }
                        }

                        updateTime+=updateTimeInd;




                        IEProp=IEProp3;
                        St=St2;

                        }

                        else

                        {

                            St=1;
                            et=St;
                            int[] IEProp4 = new int[2];
                        IEProp4[0] = 0;
                        IEProp4[1] = 0;
                        IEProp=IEProp4;
                        merged = false;

                        }
                    }
                }
            }

        for(int i1=0;i1<4*TS;i1++)
        {
            if(i1<TS)
            {
                dataOutPut[i1]=totetEv[i1];
            }
            else
            {
                if(i1<2*TS)
                {
                    dataOutPut[i1]=totStEv[i1-TS];
                }
                else
                {
                    if(i1<3*TS)
                    {
                    dataOutPut[i1]=totTau[i1-2*TS];
                    }
                    else
                    {
                        dataOutPut[i1]=sampleSize[i1-3*TS];
                    }
                }
            }
        }


        return dataOutPut;
        }

        static int[] totEtandTotTauMethodNoSplit(double pd, double pr, int TS, int iterations1)
        {
            int[] totetEv = new int[TS];
            int[] totStEv = new int[TS];
            int[] totTau = new int[TS];
            int[] sampleSize = new int[TS];

            int[] dataOutPut = new int[4*TS];

            for(int i1=0;i1<TS;i1++)
            {
                totetEv[i1]=0;
                totStEv[i1]=0;
                totTau[i1]=0;
                sampleSize[i1]=0;
            }
            
            double pt;
            pt=pd+pr;
            int St;
            int St2;
            int et;

            bool merged;
            var random = new Random();
            for(int i3=0; i3 < iterations1; i3++)
            {
            int[] IEProp = new int[2];
            int[] PIPprop = new int[TS];

            int updateTime = 0;

            IEProp[0] = 0;
            IEProp[1] = 0;
            merged = false;
            St = 1;
            et = St;


                for(int i1=0; i1<TS; i1++)
                {
                    int shift = 0;
                    int StDiff=0;
                    int[] actionsArray = new int[St];
                    int[] shiftArray = new int[St];
                    shiftArray[0]=0;

                    PIPprop[i1]=St;



                    if(merged==false)
                    {
                        for(int i2=0;i2<St;i2++)
                        {
                            double rand = random.NextDouble();
                            if(rand<pd)
                            {
                                actionsArray[i2] = 1;
                                
                                StDiff+=1;
                            }
                            else
                            {
                                if(rand<pt)
                                {
                                    actionsArray[i2] = -1;
                                    StDiff-=1;
                                }
                                else
                                {
                                    actionsArray[i2] = 0;
                                }
                            }
                        }
                        St2=St+StDiff;
                        if(St2==0){merged=true;}
                        if(merged==false)
                        {
                        int[] IEProp2 = new int[St2+1];
                        for(int i2=0;i2<St2+1;i2++)
                        {
                            IEProp2[i2]=0;
                        }
                        for(int i2=0;i2<St;i2++)
                        {
                            if(actionsArray[i2]==1)
                            {
                                IEProp2[i2+shift]+=1;

                                IEProp2[i2+shift+2]=IEProp[i2+1];
                                shift+=1;
                            }
                            else
                            {
                                if(actionsArray[i2]==-1)
                                {
                                    IEProp2[i2+shift]+=IEProp[i2+1]-1;
                                    shift-=1;
                                }
                                else
                                {
                                    IEProp2[i2+shift+1]+=IEProp[i2+1];
                                }
                            }
                        }

                        et+=IEProp2[0];

                        int[] IEProp3 = new int[St2+1];
                        for(int i2=0;i2<St2;i2++){IEProp3[i2]=IEProp2[i2+1];}
                        et+=IEProp3[0];
                        IEProp3[0]=0;
                        IEProp3[St2]=0;

                        //Test information conservation
                        //int At=0;for(int i4=0;i4<St2+1;i4++){At+=IEProp3[i4];}if(et+At!=St2){Console.WriteLine("Error with St2=" + St2 + " and At=" + At + " and et=" + et);}

                        //The following code runs for shifting the fully traversing PIPs.



                        for(int i2=0;i2<St-1;i2++)
                        {
                            shiftArray[i2+1]=(shiftArray[i2]+actionsArray[i2]);
                        }
                        int ind1=0;
                        int count1=0;
                        for(int i2=0;i2<St;i2++)
                        {
                            if(PIPprop[updateTime+ind1]==i2)
                            {
                                count1=0;
                                while(PIPprop[updateTime+ind1+count1]==i2)
                                {
                                    PIPprop[updateTime+ind1+count1]+=shiftArray[i2];
                                    count1++;
                                    if(ind1+count1>i2-updateTime){break;}
                                    ind1+=count1;
                                }
                            }
                        }

                        int updateTimeInd=0;

                        for(int i2=0;i2<i1-updateTime+1;i2++)
                        {
                            PIPprop[updateTime+i2]-=1;
                            if(PIPprop[updateTime+i2]<1)
                            {
                                updateTimeInd+=1;
                                totetEv[i1]+=et;
                                totStEv[i1]+=St2;
                                totTau[i1]+=(i1-updateTime-i2);
                                sampleSize[i1]+=1;
                            }
                        }

                        updateTime+=updateTimeInd;

                        IEProp=IEProp3;
                        St=St2;

                        }

                        else

                        {
                            break;
                        }
                    }
                }
            }

        for(int i1=0;i1<4*TS;i1++)
        {
            if(i1<TS)
            {
                dataOutPut[i1]=totetEv[i1];
            }
            else
            {
                if(i1<2*TS)
                {
                    dataOutPut[i1]=totStEv[i1-TS];
                }
                else
                {
                    if(i1<3*TS)
                    {
                    dataOutPut[i1]=totTau[i1-2*TS];
                    }
                    else
                    {
                        dataOutPut[i1]=sampleSize[i1-3*TS];
                    }
                }
            }
        }





        return dataOutPut;
        }

//The following give returns an output once an output has been formed from unity initial state. 

        static int[] totEtandTotTauMethodwithSplitOutputNoIt(double pd, double pr, int TS)
        {
            int [] pathOutput = new int[1];
            
            double pt;
            pt=pd+pr;
            int St;
            int St2;
            int et;

            bool outputCreated=false;
            int outputAttempts=0;



            bool merged;
            var random = new Random();

            while(outputCreated==false)
            {

                int[] IEProp = new int[2];
                int[] PIPprop = new int[TS];

                int updateTime = 0;

                IEProp[0] = 0;
                IEProp[1] = 0;
                merged = false;
                St = 1;
                et = St;



                for(int i1=0; i1<TS; i1++)
                {
                int shift = 0;
                int StDiff=0;
                int[] actionsArray = new int[St];
                int[] shiftArray = new int[St];
                shiftArray[0]=0;
                PIPprop[i1]=St;

                    if(merged==false)
                    {
                        for(int i2=0;i2<St;i2++)
                        {
                            double rand = random.NextDouble();
                            if(rand<pd)
                            {
                                actionsArray[i2] = 1;
                                
                                StDiff+=1;
                            }
                            else
                            {
                                if(rand<pt)
                                {
                                    actionsArray[i2] = -1;
                                    StDiff-=1;
                                }
                                else
                                {
                                    actionsArray[i2] = 0;
                                }
                            }
                        }
                        St2=St+StDiff;
                        if(St2==0){merged=true;}

                        if(merged==false)
                        {
                        int[] IEProp2 = new int[St2+1];
                        for(int i2=0;i2<St2+1;i2++)
                        {
                            IEProp2[i2]=0;
                        }
                        for(int i2=0;i2<St;i2++)
                        {
                            if(actionsArray[i2]==1)
                            {
                                IEProp2[i2+shift]+=1;

                                IEProp2[i2+shift+2]=IEProp[i2+1];
                                shift+=1;
                            }
                            else
                            {
                                if(actionsArray[i2]==-1)
                                {
                                    IEProp2[i2+shift]+=IEProp[i2+1]-1;
                                    shift-=1;
                                }
                                else
                                {
                                    IEProp2[i2+shift+1]+=IEProp[i2+1];
                                }
                            }
                        }

                        et+=IEProp2[0];

                        int[] IEProp3 = new int[St2+1];
                        for(int i2=0;i2<St2;i2++){IEProp3[i2]=IEProp2[i2+1];}
                        et+=IEProp3[0];
                        IEProp3[0]=0;
                        IEProp3[St2]=0;

                        //Test information conservation
                        //int At=0;for(int i4=0;i4<St2+1;i4++){At+=IEProp3[i4];}if(et+At!=St2){Console.WriteLine("Error with St2=" + St2 + " and At=" + At + " and et=" + et);}

                        //The following code runs for shifting the fully traversing PIPs.



                        for(int i2=0;i2<St-1;i2++)
                        {
                            shiftArray[i2+1]=(shiftArray[i2]+actionsArray[i2]);
                        }
                        int ind1=0;
                        int count1=0;
                        for(int i2=0;i2<St;i2++)
                        {
                            if(PIPprop[updateTime+ind1]==i2)
                            {
                                count1=0;
                                while(PIPprop[updateTime+ind1+count1]==i2)
                                {
                                    PIPprop[updateTime+ind1+count1]+=shiftArray[i2];
                                    count1++;
                                    if(ind1+count1>i2-updateTime){break;}
                                    ind1+=count1;
                                }
                            }
                        }

                        int updateTimeInd=0;

                        for(int i2=0;i2<i1-updateTime+1;i2++)
                        {
                            PIPprop[updateTime+i2]-=1;
                            if(PIPprop[updateTime+i2]<1){updateTimeInd+=1;}
                        }

                        updateTime+=updateTimeInd;

                        IEProp=IEProp3;
                        St=St2;

                        }

                        else

                        {
                            outputAttempts+=1;
                            break;
                        }
                    }
                }
            
            //Decompose inputs - Order: 1.IEprop 2.PIPprop 3.updateTime 4.st 5.et
            //inputLength1 and inputLength2
            //work out length of PIPprop
            
            //In the following storage allocation we define each individual path output to be separated into 4 sections.
            //1. pathOutput[0] tells us the number of outputs integers;
            //2. pathOutput[1] tells us the number of outputs with length greater than 1;
            //3. pathOutput[2] to pathOutput[3] tells us the length IEprop and PIPprop respectively;
            //4. pathOutput[3] to pathOutput[2+pathOutput[1]] gives IEprop.
            //5. pathOutput[3+pathOutput[1]] to pathOutput
            //5. All outputs with length 1 are put at the end. In this case the order is 1. updateTime 4. St 5. et

            if(merged==false)
            {
                if(St>1.5*et)
                {
            outputCreated=true;
            int pathOutputLength=1+2+(St+1)+(TS-updateTime)+3;
            pathOutput=new int[pathOutputLength];
            pathOutput=new int [pathOutputLength];
            pathOutput[0]=pathOutputLength;
            pathOutput[1]=St+1;
            pathOutput[2]=TS-updateTime;
            for(int i1=0;i1<St+1;i1++){pathOutput[i1+3]=IEProp[i1];}
            for(int i1=0;i1<(TS-updateTime);i1++){pathOutput[i1+pathOutput[1]+3]=PIPprop[updateTime+i1];}
            pathOutput[pathOutputLength-3]=updateTime;
            pathOutput[pathOutputLength-2]=St;
            pathOutput[pathOutputLength-1]=et;
                }
            }
            if(outputCreated==true){break;}


        }

        return pathOutput;

        }

        static int[] totEtandTotTauMethodwithSplitOutputNoItWithConditions(double pd, double pr, int TS, int[] minMaxSt, int[] minMaxEt)
        {
            int [] pathOutput = new int[1];
            
            double pt;
            pt=pd+pr;
            int St;
            int St2;
            int et;

            bool outputCreated=false;
            int outputAttempts=0;



            bool merged;
            var random = new Random();

            while(outputCreated==false)
            {

                int[] IEProp = new int[2];
                int[] PIPprop = new int[TS];

                int updateTime = 0;

                IEProp[0] = 0;
                IEProp[1] = 0;
                merged = false;
                St = 1;
                et = St;



                for(int i1=0; i1<TS; i1++)
                {
                int shift = 0;
                int StDiff=0;
                int[] actionsArray = new int[St];
                int[] shiftArray = new int[St];
                shiftArray[0]=0;
                PIPprop[i1]=St;

                    if(merged==false)
                    {
                        for(int i2=0;i2<St;i2++)
                        {
                            double rand = random.NextDouble();
                            if(rand<pd)
                            {
                                actionsArray[i2] = 1;
                                
                                StDiff+=1;
                            }
                            else
                            {
                                if(rand<pt)
                                {
                                    actionsArray[i2] = -1;
                                    StDiff-=1;
                                }
                                else
                                {
                                    actionsArray[i2] = 0;
                                }
                            }
                        }
                        St2=St+StDiff;
                        if(St2==0){merged=true;}

                        if(merged==false)
                        {
                        int[] IEProp2 = new int[St2+1];
                        for(int i2=0;i2<St2+1;i2++)
                        {
                            IEProp2[i2]=0;
                        }
                        for(int i2=0;i2<St;i2++)
                        {
                            if(actionsArray[i2]==1)
                            {
                                IEProp2[i2+shift]+=1;

                                IEProp2[i2+shift+2]=IEProp[i2+1];
                                shift+=1;
                            }
                            else
                            {
                                if(actionsArray[i2]==-1)
                                {
                                    IEProp2[i2+shift]+=IEProp[i2+1]-1;
                                    shift-=1;
                                }
                                else
                                {
                                    IEProp2[i2+shift+1]+=IEProp[i2+1];
                                }
                            }
                        }

                        et+=IEProp2[0];

                        int[] IEProp3 = new int[St2+1];
                        for(int i2=0;i2<St2;i2++){IEProp3[i2]=IEProp2[i2+1];}
                        et+=IEProp3[0];
                        IEProp3[0]=0;
                        IEProp3[St2]=0;

                        //Test information conservation
                        //int At=0;for(int i4=0;i4<St2+1;i4++){At+=IEProp3[i4];}if(et+At!=St2){Console.WriteLine("Error with St2=" + St2 + " and At=" + At + " and et=" + et);}

                        //The following code runs for shifting the fully traversing PIPs.



                        for(int i2=0;i2<St-1;i2++)
                        {
                            shiftArray[i2+1]=(shiftArray[i2]+actionsArray[i2]);
                        }
                        int ind1=0;
                        int count1=0;
                        for(int i2=0;i2<St;i2++)
                        {
                            if(PIPprop[updateTime+ind1]==i2)
                            {
                                count1=0;
                                while(PIPprop[updateTime+ind1+count1]==i2)
                                {
                                    PIPprop[updateTime+ind1+count1]+=shiftArray[i2];
                                    count1++;
                                    if(ind1+count1>i2-updateTime){break;}
                                    ind1+=count1;
                                }
                            }
                        }

                        int updateTimeInd=0;

                        for(int i2=0;i2<i1-updateTime+1;i2++)
                        {
                            PIPprop[updateTime+i2]-=1;
                            if(PIPprop[updateTime+i2]<1){updateTimeInd+=1;}
                        }

                        updateTime+=updateTimeInd;

                        IEProp=IEProp3;
                        St=St2;

                        }

                        else

                        {
                            outputAttempts+=1;
                            break;
                        }
                    }
                }
            
            //Decompose inputs - Order: 1.IEprop 2.PIPprop 3.updateTime 4.st 5.et
            //inputLength1 and inputLength2
            //work out length of PIPprop
            
            //In the following storage allocation we define each individual path output to be separated into 4 sections.
            //1. pathOutput[0] tells us the number of outputs integers;
            //2. pathOutput[1] tells us the number of outputs with length greater than 1;
            //3. pathOutput[2] to pathOutput[3] tells us the length IEprop and PIPprop respectively;
            //4. pathOutput[3] to pathOutput[2+pathOutput[1]] gives IEprop.
            //5. pathOutput[3+pathOutput[1]] to pathOutput[3+pathOutput[1]+pathOutput[2]]
            //5. All outputs with length 1 are put at the end. In this case the order is 1. updateTime 2. St 3. et

            if(merged==false)
            {
                if(St>minMaxSt[0])
                {
                    if(St<minMaxSt[1])
                    {
                        if(et>minMaxEt[0])
                        {
                            if(et<minMaxEt[1])
                            {
            outputCreated=true;
            int pathOutputLength=1+3+(St+1)+(TS-updateTime)+3;
            pathOutput=new int[pathOutputLength];
            pathOutput[0]=pathOutputLength;
            pathOutput[1]=2;
            pathOutput[2]=St+1;
            pathOutput[3]=TS-updateTime;
            for(int i1=0;i1<St+1;i1++){pathOutput[i1+4]=IEProp[i1];}
            for(int i1=0;i1<(TS-updateTime);i1++){pathOutput[i1+pathOutput[2]+4]=PIPprop[updateTime+i1];}
            pathOutput[pathOutputLength-3]=updateTime;
            pathOutput[pathOutputLength-2]=St;
            pathOutput[pathOutputLength-1]=et;
                            }
                        }
                    }
                }
            }
            if(outputCreated==true){break;}


        }

        return pathOutput;

        }


//The following is to translate the output format into the input format.


//The following are methods that take in input states but give results as outputs.

        static int[] totEtandTotTauMethodnoSplitWithInputs(double pd, double pr, int TS, int iterations1, int[] inputs)
        {




            int[] totetEv = new int[TS];
            int[] totStEv = new int[TS];
            int[] totTau = new int[TS];
            int[] sampleSize = new int[TS];

            int[] dataOutPut = new int[4*TS];

            //Decompose inputs - Order: 1.IEprop 2.PIPprop 3.updateTime 4.merged 5.St 6.et
            int inputLength1=inputs[2];
            int inputLength2=inputs[3];
            int[] IEpropIn = new int[inputLength1];
            for(int i1=0;i1<inputLength1;i1++){IEpropIn[i1]=inputs[4+i1];}
            int[] PIPpropIn = new int[inputLength2];
            for(int i1=0;i1<inputLength2;i1++){PIPpropIn[i1]=inputs[4+inputLength1+i1];}
            int updateTimeIn=inputs[4+inputLength1+inputLength2];
            int StIn=inputs[4+inputLength1+inputLength2 + 1];
            int EtIn=inputs[4+inputLength1+inputLength2 + 2];



            for(int i1=0;i1<TS;i1++)
            {
                totetEv[i1]=0;
                totStEv[i1]=0;
                totTau[i1]=0;
                sampleSize[i1]=0;
            }

            double pt;
            pt=pd+pr;
            int St;
            int St2;
            int et;

            bool merged;
            var random = new Random();
            for(int i3=0;i3<iterations1;i3++)
            {
            int[] IEProp = IEpropIn;
            int TSoverlap=inputLength2+TS;
            int[] PIPprop = new int[TSoverlap];
            for(int i1=0;i1<inputLength2;i1++){PIPprop[i1]=PIPpropIn[i1];}
            

            int updateTime = 0;

            int TS1=updateTimeIn+inputLength2;



            merged = false;
            St = StIn;
            et = EtIn;

                 for(int i1=0; i1<TS; i1++)
                {
                    int shift = 0;
                    int StDiff=0;
                    int[] actionsArray = new int[St];
                    int[] shiftArray = new int[St];
                    shiftArray[0]=0;

                    PIPprop[inputLength2 + i1]=St;



                    if(merged==false)
                    {
                        for(int i2=0;i2<St;i2++)
                        {
                            double rand = random.NextDouble();
                            if(rand<pd)
                            {
                                actionsArray[i2] = 1;
                                
                                StDiff+=1;
                            }
                            else
                            {
                                if(rand<pt)
                                {
                                    actionsArray[i2] = -1;
                                    StDiff-=1;
                                }
                                else
                                {
                                    actionsArray[i2] = 0;
                                }
                            }
                        }
                        St2=St+StDiff;
                        if(St2==0){merged=true;}
                        if(merged==false)
                        {
                        int[] IEProp2 = new int[St2+1];
                        for(int i2=0;i2<St2+1;i2++)
                        {
                            IEProp2[i2]=0;
                        }
                        for(int i2=0;i2<St;i2++)
                        {
                            if(actionsArray[i2]==1)
                            {
                                IEProp2[i2+shift]+=1;

                                IEProp2[i2+shift+2]=IEProp[i2+1];
                                shift+=1;
                            }
                            else
                            {
                                if(actionsArray[i2]==-1)
                                {
                                    IEProp2[i2+shift]+=IEProp[i2+1]-1;
                                    shift-=1;
                                }
                                else
                                {
                                    IEProp2[i2+shift+1]+=IEProp[i2+1];
                                }
                            }
                        }

                        et+=IEProp2[0];

                        int[] IEProp3 = new int[St2+1];
                        for(int i2=0;i2<St2;i2++){IEProp3[i2]=IEProp2[i2+1];}
                        et+=IEProp3[0];
                        IEProp3[0]=0;
                        IEProp3[St2]=0;

                        //Test information conservation
                        //int At=0;for(int i4=0;i4<St2+1;i4++){At+=IEProp3[i4];}if(et+At!=St2){Console.WriteLine("Error with St2=" + St2 + " and At=" + At + " and et=" + et);}

                        //The following code runs for shifting the fully traversing PIPs.



                        for(int i2=0;i2<St-1;i2++)
                        {
                            shiftArray[i2+1]=(shiftArray[i2]+actionsArray[i2]);
                        }
                        int ind1=0;
                        int count1=0;
                        for(int i2=0;i2<St;i2++)
                        {
                            if(PIPprop[updateTime+ind1]==i2)
                            {
                                count1=0;
                                while(PIPprop[updateTime+ind1+count1]==i2)
                                {
                                    PIPprop[updateTime+ind1+count1]+=shiftArray[i2];
                                    count1++;
                                    if(ind1+count1>i2-updateTime){break;}
                                    ind1+=count1;
                                }
                            }
                        }

                        int updateTimeInd=0;

                        for(int i2=0;i2<i1-updateTime+inputLength2;i2++)
                        {
                            PIPprop[updateTime+i2]-=1;
                            if(PIPprop[updateTime+i2]<1)
                            {
                                updateTimeInd+=1;
                                totetEv[i1]+=et;
                                totStEv[i1]+=St2;
                                totTau[i1]+=(i1+TS1)-(updateTimeIn+updateTime+i2);
                                sampleSize[i1]+=1;
                            }
                        }

                        updateTime+=updateTimeInd;

                        IEProp=IEProp3;
                        St=St2;

                        }

                        else

                        {
                            break;
                        }
                    }
                }
            
            }


        for(int i1=0;i1<4*TS;i1++)
        {
            if(i1<TS)
            {
                dataOutPut[i1]=totetEv[i1];
            }
            else
            {
                if(i1<2*TS)
                {
                    dataOutPut[i1]=totStEv[i1-TS];
                }
                else
                {
                    if(i1<3*TS)
                    {
                    dataOutPut[i1]=totTau[i1-2*TS];
                    }
                    else
                    {
                        dataOutPut[i1]=sampleSize[i1-3*TS];
                    }
                }
            }
        }


        return dataOutPut;





        }

//The following take a random input from an ensamble for each loop in the monte carlo loop.

 static int[] totEtandTotTauMethodnoSplitWithInputsChosenAtRandomFromEnsemble(double pd, double pr, int TS, int iterations1, int[] inputsEnsemble)
        {





            int[] totetEv = new int[TS];
            int[] totStEv = new int[TS];
            int[] totTau = new int[TS];
            int[] sampleSize = new int[TS];

            int[] dataOutPut = new int[4*(TS)];

            for(int i1=0;i1<TS;i1++)
            {
                totetEv[i1]=0;
                totStEv[i1]=0;
                totTau[i1]=0;
                sampleSize[i1]=0;
            }




            double pt;
            pt=pd+pr;
            int St;
            int St2;
            int et;

            bool merged;
            var random = new Random();
            for(int i3=0;i3<iterations1;i3++)
            {


            int[] inputs=extractInitialStateFromEnsemble(inputsEnsemble);


            int inputLength1=inputs[2];
            int inputLength2=inputs[3];
            int[] IEpropIn = new int[inputLength1];
            for(int i1=0;i1<inputLength1;i1++){IEpropIn[i1]=inputs[i1+4];}
            int[] PIPpropIn = new int[inputLength2];
            for(int i1=0;i1<inputLength2;i1++){PIPpropIn[i1]=inputs[inputLength1+i1+4];}
            int updateTimeIn=inputs[4+inputLength1+inputLength2];
            int StIn=inputs[inputLength1+inputLength2 + 5];
            int EtIn=inputs[inputLength1+inputLength2 + 6];



            int[] IEProp = IEpropIn;
            int TSoverlap=inputLength2+TS;
            int[] PIPprop = new int[TSoverlap];
            for(int i1=0;i1<inputLength2;i1++){PIPprop[i1]=PIPpropIn[i1];}
            

            int updateTime = 0;

            int TS1=updateTimeIn+inputLength2;

            



            merged = false;
            St = StIn;
            et = EtIn;

                 for(int i1=0; i1<TS; i1++)
                {
                    int shift = 0;
                    int StDiff=0;
                    int[] actionsArray = new int[St];
                    int[] shiftArray = new int[St];
                    shiftArray[0]=0;

                    PIPprop[inputLength2 + i1]=St;



                    if(merged==false)
                    {
                        for(int i2=0;i2<St;i2++)
                        {
                            double rand = random.NextDouble();
                            if(rand<pd)
                            {
                                actionsArray[i2] = 1;
                                
                                StDiff+=1;
                            }
                            else
                            {
                                if(rand<pt)
                                {
                                    actionsArray[i2] = -1;
                                    StDiff-=1;
                                }
                                else
                                {
                                    actionsArray[i2] = 0;
                                }
                            }
                        }
                        St2=St+StDiff;
                        if(St2==0){merged=true;}
                        if(merged==false)
                        {
                        int[] IEProp2 = new int[St2+1];
                        for(int i2=0;i2<St2+1;i2++)
                        {
                            IEProp2[i2]=0;
                        }
                        for(int i2=0;i2<St;i2++)
                        {
                            if(actionsArray[i2]==1)
                            {
                                IEProp2[i2+shift]+=1;

                                IEProp2[i2+shift+2]=IEProp[i2+1];
                                shift+=1;
                            }
                            else
                            {
                                if(actionsArray[i2]==-1)
                                {
                                    IEProp2[i2+shift]+=IEProp[i2+1]-1;
                                    shift-=1;
                                }
                                else
                                {
                                    IEProp2[i2+shift+1]+=IEProp[i2+1];
                                }
                            }
                        }

                        et+=IEProp2[0];

                        int[] IEProp3 = new int[St2+1];
                        for(int i2=0;i2<St2;i2++){IEProp3[i2]=IEProp2[i2+1];}
                        et+=IEProp3[0];
                        IEProp3[0]=0;
                        IEProp3[St2]=0;

                        //Test information conservation
                        //int At=0;for(int i4=0;i4<St2+1;i4++){At+=IEProp3[i4];}if(et+At!=St2){Console.WriteLine("Error with St2=" + St2 + " and At=" + At + " and et=" + et);}

                        //The following code runs for shifting the fully traversing PIPs.



                        for(int i2=0;i2<St-1;i2++)
                        {
                            shiftArray[i2+1]=(shiftArray[i2]+actionsArray[i2]);
                        }
                        int ind1=0;
                        int count1=0;
                        for(int i2=0;i2<St;i2++)
                        {
                            if(PIPprop[updateTime+ind1]==i2)
                            {
                                count1=0;
                                while(PIPprop[updateTime+ind1+count1]==i2)
                                {
                                    PIPprop[updateTime+ind1+count1]+=shiftArray[i2];
                                    count1++;
                                    if(ind1+count1>i2-updateTime){break;}
                                    ind1+=count1;
                                }
                            }
                        }

                        int updateTimeInd=0;

                        for(int i2=0;i2<i1-updateTime+inputLength2;i2++)
                        {
                            PIPprop[updateTime+i2]-=1;
                            if(PIPprop[updateTime+i2]<1)
                            {
                                updateTimeInd+=1;
                                totetEv[i1]+=et;
                                totStEv[i1]+=St2;
                                totTau[i1]+=(i1+TS1)-(updateTimeIn+updateTime+i2);
                                sampleSize[i1]+=1;
                            }
                        }

                        updateTime+=updateTimeInd;

                        IEProp=IEProp3;
                        St=St2;

                        }

                        else

                        {
                            break;
                        }
                    }
                }
            
            }


        for(int i1=0;i1<4*TS;i1++)
        {
            if(i1<TS)
            {
                dataOutPut[i1]=totetEv[i1];
            }
            else
            {
                if(i1<2*TS)
                {
                    dataOutPut[i1]=totStEv[i1-TS];
                }
                else
                {
                    if(i1<3*TS)
                    {
                    dataOutPut[i1]=totTau[i1-2*TS];
                    }
                    else
                    {
                        dataOutPut[i1]=sampleSize[i1-3*TS];
                    }
                }
            }
        }


        return dataOutPut;





        }

//The following method takes in a random initial state generator and creates length marked ensamble integer array as output.

        static int[] createStructuredInitialStateEnsambleFromUnionNoSplits(double pd, double pr, int TS, int sampleSize, int[] minMaxSt, int[] minMaxEt)
        {

            int countTotalSampleIntSize=0;
            int[] sampleSizes = new int[sampleSize];
            int[] sampleEnsableInputs=new int[1];


            for(int i1=0;i1<sampleSize;i1++)
            {
            int[] input1=totEtandTotTauMethodwithSplitOutputNoItWithConditions(pd, pr, TS, minMaxSt, minMaxEt);
            sampleSizes[i1]=+input1[0];
            int countTotalSampleIntSize2=countTotalSampleIntSize+input1[0];

            if(i1==0)
            {
            sampleEnsableInputs=new int[countTotalSampleIntSize2];
                for(int i2=0;i2<countTotalSampleIntSize2;i2++)
                {
                    sampleEnsableInputs[i2]=input1[i2];
                }
            }

            if(i1>0)
            {
            int[] sampleEnsableInputs2=new int[countTotalSampleIntSize2];
            for(int i2=0;i2<countTotalSampleIntSize2;i2++)
            {
                if(i2<countTotalSampleIntSize)
                {
                    sampleEnsableInputs2[i2]=sampleEnsableInputs[i2];
                }
                else
                {
                    sampleEnsableInputs2[i2]=input1[i2-countTotalSampleIntSize];
                }
            }
            sampleEnsableInputs=sampleEnsableInputs2;
            }


            countTotalSampleIntSize+=input1[0];

            }

            int totalArrayLength=2+sampleSize+countTotalSampleIntSize;



//The output ensemble has following form.
//1. [0] = to the total size of the ensemble (array).
//2. [1] = the total size of the sample.
//2. [2] to [n] = the size of each ensamble output state.
//3. [n+1] to [l] = all input data partitioned by the integers in the set from [1] to [n] previous that gives the information of the ensemble of states. 

            int[] outputEnsembleInitialStateArray=new int[totalArrayLength];
            for(int i2=0;i2<totalArrayLength;i2++)
            {
                if(i2==0)
                {
                    outputEnsembleInitialStateArray[0]=totalArrayLength;
                }
                else
                {
                    if(i2==1)
                    {
                        outputEnsembleInitialStateArray[1]=sampleSize;
                    }
                    else
                    {
                        if(i2<(2+sampleSize))
                        {
                            outputEnsembleInitialStateArray[i2]=sampleSizes[i2-2];
                        }
                        else
                        {
                        outputEnsembleInitialStateArray[i2]=sampleEnsableInputs[i2-sampleSize-2];
                        }
                    }
                }

            }







return outputEnsembleInitialStateArray;

        }

//The following method picks a random state out of the initial state ensamble.

        static int[] extractInitialStateFromEnsemble(int[] inputEnsemble)
        {
            int sampleSize=inputEnsemble[1];
            Random rand1 = new Random();
            int chosenInitialStateIndex  = rand1.Next(0, sampleSize);
            int chosenInitialStateAddress = 2+sampleSize;

            if(chosenInitialStateIndex>0)
            {
                for(int i1=0;i1<chosenInitialStateIndex;i1++){chosenInitialStateAddress+=inputEnsemble[chosenInitialStateAddress];}
            }

            int[]  initialState=new int[inputEnsemble[chosenInitialStateAddress]];
            for(int i1=0;i1<inputEnsemble[chosenInitialStateAddress];i1++)
            {
                initialState[i1]=inputEnsemble[chosenInitialStateAddress+i1];
            }
            return initialState;
        }

        static int extractArrayMaxSize(int[] inputEnsemble)
        {
            int maxTime1=0;
            int maxTime2=0;
            int PIPpropArraySizeAddress=0;

            for(int i1=0; i1<inputEnsemble[1]; i1++)
            {
                if(i1==0)
                {
                    PIPpropArraySizeAddress=1+inputEnsemble[1]+3;
                    maxTime1=inputEnsemble[PIPpropArraySizeAddress];
                }
                if(i1>0)
                {
                    PIPpropArraySizeAddress+=inputEnsemble[PIPpropArraySizeAddress-2];
                    maxTime2=inputEnsemble[PIPpropArraySizeAddress];
                }
                if(maxTime2>maxTime1){maxTime1=maxTime2;}

            }
            return maxTime2;
        }

        static void Main(string[] args)
        {
        //int TS=500;
        //int iterations1=1000000;
        //double pd = 0.3333333333333f;
        //double pr = 0.3333333333333f;

            //int[] systEv=totEtandTotTauMethodwithSplit(pd, pr, TS, iterations1);
            //printDataInListMathematicaFormat(systEv,4*TS, "totEtandTotTauMethodwithSplit.txt");

        double pd = 0.3333333333333f;
        double pr = 0.3333333333333f;


        int TS1=200;
        //int iterations1=100000;
        int sampleSize=10;


        
        int TS2=1000;
        int iterations2=10000;
    

int StScale=200;
int StHalfRange=20;
int EtScale=200;
int EtHalfRange=20;

int [] minMaxSt = new int[2];
minMaxSt[0]=StScale-StHalfRange;
minMaxSt[1]=StScale+StHalfRange;
int [] minMaxEt = new int[2];
minMaxEt[0]=EtScale-EtHalfRange;
minMaxEt[1]=EtScale+EtHalfRange;


           //int[] EtEvFromUnion=totEtEvItnoSplitTayloredEnsembleSize(pd, pr, TS1, iterations1);
           //int[] input1=totEtandTotTauMethodwithSplitOutputNoItWithConditions(pd, pr, TS1, minMaxSt, minMaxEt);


        int[] inputEnsemble1=createStructuredInitialStateEnsambleFromUnionNoSplits(pd, pr, TS1, sampleSize, minMaxSt,minMaxEt);

              Console.Beep();
              Console.WriteLine("Initial states found.");

              //int[] input1=extractInitialStateFromEnsemble(inputEnsemble1);


              //int[] systEv=totEtandTotTauMethodnoSplitWithInputs(pd,pr,TS2,iterations2,input1);

              int[] systEv=totEtandTotTauMethodnoSplitWithInputsChosenAtRandomFromEnsemble(pd, pr, TS2, iterations2, inputEnsemble1);
              printDataInListMathematicaFormat(systEv, 4*TS2, "totEtandTotTauMethodnoSplitWithInputsChosenAtRandomFromEnsemble.txt");

              Console.Beep();
        }

        static void printDataInListMathematicaFormat(int[] data,int TS, string filepath)
        {

            string systEvToString = "{";
            for(int i2=0;i2<TS;i2++)
            {  
                if(i2==0)
                {
                systEvToString = systEvToString + data[i2].ToString();       
                }
                else
                {
                if(i2<TS-1)
                {
                systEvToString = systEvToString + "," + data[i2].ToString();
                }
                else
                {
                systEvToString = systEvToString + "," + data[i2].ToString() + "}";
                }
                }
            }
            writeDelDatatoCSV(systEvToString, filepath);
        }

        static void writeDatatoCSV(String etStr, String filepath)
            {
                try
                {
                    using(System.IO.StreamWriter file = new System.IO.StreamWriter(@filepath,true))
                    {
                        file.WriteLine(etStr);
                    }
                }
                catch(Exception ex)
                {
                throw new ApplicationException("Error :", ex);
                }
            }
        static void writeDelDatatoCSV(String etStr, String filepath)
            {
                try
                {
                    using(System.IO.StreamWriter file = new System.IO.StreamWriter(@filepath,false))
                    {
                        file.WriteLine(etStr);
                    }
                }
                catch(Exception ex)
                {
                throw new ApplicationException("Error :", ex);
                }
            }
            

        static void printArrayInput(int[] arrayInput)
        {
            int length=arrayInput[0];
            for(int i1=0;i1<length;i1++)
            {
            Console.WriteLine(arrayInput[i1]);
            }

            Console.WriteLine("end");

        }
    }
}




    
