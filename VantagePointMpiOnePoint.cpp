/*---------------------------------------------
ANTONAKOUDIS SOTIRIOS - IAKOBIDIS IOANNIS
            MPI VANTAGE POINT SEARCH
---------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <math.h>
#include  <time.h>
#include <sys/time.h>
#include <math.h>
#include <stdbool.h>
#include <queue>
#include <limits>
#include <algorithm>
#include <iostream>
#include <fstream>
using namespace std;
MPI_Status Stat;
float tau;
#define COUNT 10
using namespace std;
void partition (float *array, int elements, float pivot, float **arraysmall, float **arraybig, int *endsmall, int *endbig);
float selection(float *array,int number);
void split(float *pointsArray,float* elements,int size,float median,int dimension);
float selectionpPoints(float*pointsArray,int dimension,float *array,int number);
void print2DUtil(struct Node *root, int space,int dimension);
struct Node* createtreelocal(float *pointsArray,int partLength,int dimension,int processId);
struct HeapItem {
        HeapItem( float* p, double dist) :
            p(p), dist(dist) {}
        float* p;
        double dist;
        bool operator<( const HeapItem& o ) const {
            return dist < o.dist;
        }
    };
struct Node {
  float *vpoint;
  float radius;
  struct Node *left;
  struct Node *right;
  bool leaf;
};

float distance(float *pointArray,int intex,float *point,int dimension){
   float tmp;
   float dist=0;

  for(int i=0;i<dimension;i++){
        dist += (pointArray[intex*dimension+i]-point[i])*(pointArray[intex*dimension+i]-point[i]);
        }
  dist=sqrtf(dist);
  return dist;
}


/***Kills processes that have no values left in their arrays****/
void removeElement(int *array, int *size, int element)
{
    int i;
    int flag=0;
    for(i=0;i<*size;i++)
    {
        if(flag==1)
            array[i]=array[i+1];
        if(array[i]==element&& flag==0)
        {
            array[i]=array[i+1];
            flag=1;
        }
    }
    *size=*size-1;
}

/***Calculate Lengths and Send them to the corresponding Node***/
void sendLengths(int size,int noProcesses)
{
    int i,partLength;
    if(size%noProcesses!=0)
    {
        int left=size-(size/noProcesses)*noProcesses;  //Split the size in as close to equal as possible parts
        partLength=(size/noProcesses)+1;
        for(i=1;i<left;i++)      //start from 1 because we create the zero one through the main function
            MPI_Send(&partLength,1,MPI_INT,i,1,MPI_COMM_WORLD);
        partLength-=1;
        for(i=left;i<noProcesses;i++)
            MPI_Send(&partLength,1,MPI_INT,i,1,MPI_COMM_WORLD);
    }
    else
    {
        partLength=size/noProcesses;
        for(i=1;i<noProcesses;i++)
            MPI_Send(&partLength,1,MPI_INT,i,1,MPI_COMM_WORLD);
    }
}

/****Swaps two values in an array****/
void swap_values(float *array,int x,int y)
{
    float temp;
    temp=array[x];
    array[x]=array[y];
    array[y]=temp;
}

void swap_point(float *array,int x,int y,int dimension){

    for(int i=0;i<dimension;i++){
        swap_values(array,x*dimension+i,y*dimension+i);
    }

}

/*****Send random numbers to every Node.*****/
void generateNumbers(float *numberPart,int partLength, int cal)
{
    srand((cal+1)*time(NULL));     //Generate number to fill the array
    int i;
    for(i=0; i<partLength; i++)
        numberPart[i]=((float)rand()/(float)(RAND_MAX)) * 100;

}
void generatePoints(float *pointsArray,int partLength, int cal,int dimension)
{
    srand((cal+1)*time(NULL));     //Generate number to fill the array
    int i;
    for(i=0; i<partLength*dimension; i++)
        pointsArray[i]=((float)rand()/(float)(RAND_MAX)) * 100;

}

/***Validates the stability of the operation****/
void validation(float median,int partLength,int size,float *numberPart,int processId,MPI_Comm group)
{
    MPI_Bcast(&median,1,MPI_FLOAT,0,group);
      	int countMin=0;
    int countMax=0;
    int countEq=0;
    int sumMax,sumMin,sumEq,i;
    for(i=0;i<partLength;i++)
    {
        if(numberPart[i]>median)
            countMax++;
        else if(numberPart[i]<median)
            countMin++;
        else
            countEq++;
    }
    MPI_Reduce(&countMax,&sumMax,1,MPI_INT,MPI_SUM,0,group);
    MPI_Reduce(&countMin,&sumMin,1,MPI_INT,MPI_SUM,0,group);
    MPI_Reduce(&countEq,&sumEq,1,MPI_INT,MPI_SUM,0,group);
    if(processId==0)
    {
        if((sumMax<=size/2)&&(sumMin<=size/2))  //Checks if both the lower and higher values occupy less than 50% of the total array.
            printf("VALIDATION PASSED!\n");
        else
            printf("VALIDATION FAILED!\n");


	      printf("Values greater than median: %d\n",sumMax);
        printf("Values equal to median: %d\n",sumEq);
        printf("Values lower than median: %d\n",sumMin);
    }

}

/***Validates the stability of the operation (Single Threaded)****/
void validationST(float median,int size,float *numberPart)
{
	int countMin=0;
    int countMax=0;
    int countEq=0;
    int i;
    for(i=0;i<size;i++)
    {
        if(numberPart[i]>median)
            countMax++;
        else if(numberPart[i]<median)
            countMin++;
        else
            countEq++;
    }
    if((countMax<=size/2)&&(countMin<=size/2))  //Checks if both the lower and higher values occupy less than 50% of the total array.
        printf("VALIDATION PASSED!\n");
    else
        printf("VALIDATION FAILED!\n");

	printf("Values greater than median: %d\n",countMax);
        printf("Values equal to median: %d\n",countEq);
        printf("Values lower than median: %d\n",countMin);
}

/****Part executed only by the Master Node****/
float masterPart(int noProcesses,int processId,int size,int partLength,float *numberPart,MPI_Comm group) //MASTER NODE CODE
{
    int elements,i,keepBigSet,sumSets,finalize,randomNode,k;
    float median,tempPivot,pivot;
    int endSmall=0;
    int dropoutFlag=0;
    int endBig=0;
    float *arraySmall,*arrayBig,*arrayToUse;
    int* activeNodes;
    int activeSize=noProcesses;
    int stillActive=1;
    int oldSumSets=-1;
    int checkIdentical=0;
    int useNewPivot=0;
    int *pivotArray;
    k=(int)size/2+1; //It is done so in order to find the right median in an even numbered array.
    elements=partLength;
    activeNodes=(int *)malloc(noProcesses*sizeof(int));  //we create the array that contains the active Nodes.
    arrayToUse=numberPart;
    pivotArray=(int*)malloc(noProcesses*sizeof(int));  //Used for special occasions to gather values different than the pivot.
    for(i=0;i<activeSize;i++)
    {
        activeNodes[i]=i;
    }
    int randomCounter=0;
    int randomCounter2=0;
    struct timeval first, second, lapsed;
    struct timezone tzp;
    gettimeofday(&first, &tzp);
    for(;;)   //Begin the infinite loop until the median is found.
    {
        int counter=0;
        useNewPivot=0;
        if(stillActive==1&&checkIdentical!=0)  //If i still have values in my array and the Sumed Big Set is identical to the previous one, check for identical values.
        {
            for(i=0;i<elements;i++)
            {
                if(pivot==arrayToUse[i])
                    counter++;
                else
                {
                    useNewPivot=1;
                    tempPivot=arrayToUse[i];
                    break;
                }
            }
        }
        if(checkIdentical!=0)
        {
            int useNewPivotMax=0;
	        MPI_Reduce(&useNewPivot,&useNewPivotMax,1,MPI_INT,MPI_MAX,0,group); //FIRST(OPTIONAL) REDUCE : MAX useNewPivot
            if(useNewPivotMax!=1)    //That means that the only values left are equal to the pivot!
            {
                median=pivot;
                finalize=1;
                MPI_Bcast(&finalize,1,MPI_INT,0,group); //FIRST(OPTIONAL) BROADCAST : WAIT FOR FINALIZE COMMAND OR NOT

                validation(median,partLength,size,numberPart,processId,group);
                printf("validated result");
                free(pivotArray);
                return median;
            }
            else
            {
                finalize=0;
                int useit=0;
                randomCounter2++;
                MPI_Bcast(&finalize,1,MPI_INT,0,group);
                MPI_Gather(&useNewPivot, 1, MPI_INT, pivotArray, 1, MPI_INT, 0, group); //Gather every value and chose a Node to change the pivot.
                for(i=0;i<activeSize;i++)
                {
                    if(pivotArray[i]==1)
                    {
                        if((randomCounter2>1)&&(randomNode!=activeNodes[i]))  //Check if the same Node has already been used in a similar operation.
                        {
                            useit=1;
                            randomNode=activeNodes[i];
                            randomCounter2=0;
                            break;
                        }
                        else if(randomCounter2<2)
                        {
                            useit=1;
                            randomNode=activeNodes[i];
                            break;
                        }
                    }
                }
                if(useit!=0)
                    useNewPivot=1;
                else
                    useNewPivot=0;
            }
        }
        if(useNewPivot!=0)
            MPI_Bcast(&randomNode,1,MPI_INT,0,group);  //THIRD(OPTIONAL) BROADCAST : BROADCAST THE SPECIAL NODE
        if(useNewPivot==0)  //if we didnt choose a special Node, choose the Node that will pick the pivot in a clockwise manner. Only selects one of the active Nodes.
        {
            if(randomCounter>=activeSize)
                randomCounter=0; //Fail safe
            randomNode=activeNodes[randomCounter];
            randomCounter++;			//Increase the counter
            MPI_Bcast(&randomNode,1,MPI_INT,0,group);   //FIRST BROADCAST : SENDING randomNode, who will chose
        }
        if(randomNode==processId)  //If i am to choose the pivot.....
	    {
            if(useNewPivot==0)
            {
                srand(time(NULL));
                pivot=arrayToUse[rand() % elements];
                MPI_Bcast(&pivot,1,MPI_FLOAT,0,group); //SECOND BROADCAST : SENDING PIVOT   k ton stelnw sto lao
	        }
            else
            {
                MPI_Bcast(&tempPivot,1,MPI_FLOAT,0,group); //SECOND BROADCAST : SENDING PIVOT   k ton stelnw sto lao
                pivot=tempPivot;
            }
        }
        else //If not.. wait for the pivot to be received.
            MPI_Bcast(&pivot,1,MPI_FLOAT,randomNode,group);  // SECOND BROADCAST : RECEIVING PIVOT
        if(stillActive==1)  //If i still have values in my array.. proceed
        {
            partition(arrayToUse,elements,pivot,&arraySmall,&arrayBig,&endSmall,&endBig);  //I partition my array  // endsmall=number of elements in small array, it may be 0
            // endbig=number of elements in big array, it may be 0
            //arraysmall = Points to the position of the small array.NULL if the array is empty
            //Same for arraybig
        }
        else  //If i'm not active endBig/endSmall has zero value.
        {
            endBig=0;
            endSmall=0;
        }
        sumSets=0;
	    //We add the bigSet Values to decide if we keep the small or the big array
	    MPI_Reduce(&endBig,&sumSets,1,MPI_INT,MPI_SUM,0,group);  //FIRST REDUCE : SUM OF BIG
        MPI_Bcast(&sumSets,1,MPI_INT,0,group);
        if(oldSumSets==sumSets)
            checkIdentical=1;
        else
        {
            oldSumSets=sumSets;
            checkIdentical=0;
        }
	    //hmetabliti keepBigSet 0 h 1 einai boolean k me autin enimerwnw ton lao ti na kratisei to bigset h to smallset
	    if(sumSets>k)   //an to sumofbigsets > k tote krataw to big SET
	    {
            keepBigSet=1; //to dilwnw auto gt meta tha to steilw se olous
            if(endBig==0)
                dropoutFlag=1; //wraia.. edw an dw oti to bigset mou einai 0.. alla prepei na kratisw to bigset sikwnw auti ti simaia pou simainei tha ginw inactive ligo pio katw tha to deis
            else
            {
                arrayToUse=arrayBig; //thetw ton neo pinaka na einai o big
                elements=endBig; //thetw arithmo stoixeiwn iso me tou big
            }
	    }
	    else if(sumSets<k) //antistoixa an to sumofbigsets < k tote krataw to small set
	    {
		    keepBigSet=0;
		    k=k-sumSets;
		    if(endSmall==0)
                dropoutFlag=1; //antistoixa koitaw an tha ginw inactive..
		    else
		    {
		    	arrayToUse=arraySmall; //dinw times..
		    	elements=endSmall;
		    }
	    }
	    else  //edw simainei k=sumofbigsetes ara briskw pivot k telos
	    {
		    median=pivot;
		    finalize=1; //dilwnw finalaize =1
		    MPI_Bcast(&finalize,1,MPI_INT,0,group); //to stelnw se olous, oi opoioi an laboun finalize =1 tote kaloun MPI finalize k telos

		    validation(median,partLength,size,numberPart,processId,group);

            free(pivotArray);
            return median;
        }
        finalize=0; //an den exw mpei sta if den exw steilei timi gia finalize.. oi alloi omws perimenoun na laboun kati, stelnw loipon to 0 pou simainei sunexizoume
        MPI_Bcast(&finalize,1,MPI_INT,0,group);	//SECOND BROADCAST : WAIT FOR FINALIZE COMMAND OR NOT
        //edw tous stelnw to keepbigset gia na doun ti tha dialeksoun
	    MPI_Bcast(&keepBigSet,1,MPI_INT,0,group);    //THIRD BROADCAST: SEND keepBigset boolean
        if(dropoutFlag==1 && stillActive==1) //edw sumfwna me to dropoutflag pou orisame prin an einai 1 kalw tin sinartisi pou me petaei apo ton pinaka. episis koitaw na eimai active gt an me exei idi petaksei se proigoumeni epanalispi tote den xreiazetai na me ksanapetaksei
        {
            stillActive=0;
            removeElement(activeNodes, &activeSize, 0);
        }
        int flag;
        //edw perimenw na akousw apo ton kathena an sunexizei active h oxi.. an oxi ton petaw.. an einai idi inactive apo prin stelnei kati allo (oxi 1)k den ton ksanapetaw
        for(i=0;i<activeSize;i++)
        {
            if(activeNodes[i]!=0)
            {
                MPI_Recv(&flag,1,MPI_INT,activeNodes[i],1,group,&Stat);  //FIRST RECEIVE : RECEIVE active or not
                if(flag==1)
                    removeElement(activeNodes, &activeSize, activeNodes[i]);
            }
        }
    }
}

/***Executed only by Slave Nodes!!*****/
void slavePart(int processId,int partLength,float *numberPart,int size,MPI_Comm group)  //code here is for the cheap slaves :P
{
	int dropoutflag,elements,i,sumSets,finalize,keepBigSet,randomNode;
    float pivot,tempPivot;
    int endSmall=0;
    int endBig=0;
    float *arraySmall,*arrayBig,*arrayToUse;
	arrayToUse=numberPart;
	elements=partLength;
	int stillActive=1;
	float *pivotArray;
    int oldSumSets=-1;
    int checkIdentical=0;
    int useNewPivot;
	for(;;)
	{
        finalize=0;
        int counter=0;
        useNewPivot=0;
        if(stillActive==1&&checkIdentical!=0)  //If i still have values in my array..   If the Sumed Big Set is identical to the previous one, check for identical values.
        {
            for(i=0;i<elements;i++)
            {
                if(pivot==arrayToUse[i])
                    counter++;
                else
                {
                    useNewPivot=1;
                    tempPivot=arrayToUse[i];
                    break;
                }
            }
        }
        if(checkIdentical!=0)
        {
            int useNewPivotMax=0;
            MPI_Reduce(&useNewPivot,&useNewPivotMax,1,MPI_INT,MPI_MAX,0,group);
            MPI_Bcast(&finalize,1,MPI_INT,0,group);//an o master apo to keepbigset k apo to count apofasisei oti teleiwsame mou stelnei 1, alliws 0 sunexizoume
            if(finalize==1)
            {
                float median=0;
                validation(median,partLength,size,numberPart,processId,group);

                return ;
            }
            else
            {
                MPI_Gather(&useNewPivot, 1, MPI_INT, pivotArray, 1, MPI_INT, 0, group);
            }
        }
        MPI_Bcast(&randomNode,1,MPI_INT,0,group); //FIRST BROAD CAST : RECEIVING RANDOM NODE, perimenw na dw poios einaito done
        if(randomNode!=processId) //means I am not the one to chose pivot.. so I wait to receive the pivot
            MPI_Bcast(&pivot,1,MPI_FLOAT,randomNode,group);	//SECOND BROADCAST : RECEIVING PIVOT
        else if(randomNode==processId) //I am choosing suckers
        {
            if(useNewPivot==0)
            {
                srand(time(NULL));
                pivot=arrayToUse[rand() % elements];
                MPI_Bcast(&pivot,1,MPI_FLOAT,processId,group); //SECOND BROADCAST : SENDING PIVOT   k ton stelnw sto lao
            }
            else
            {
                MPI_Bcast(&tempPivot,1,MPI_INT,processId,group); //SECOND BROADCAST : SENDING PIVOT   k ton stelnw sto lao
                pivot=tempPivot;
            }
        }
        if(stillActive==1)   //an eksakolouthw na eimai active, trexw tin partition.. k to count kommati to opio eimape kapou exei problima
        {
            partition(arrayToUse,elements,pivot,&arraySmall,&arrayBig,&endSmall,&endBig);
        }
        else
        {
            endBig=0;
            endSmall=0;
        }
        //an eimai inactive stelnw endbig=0 gia to bigset pou den epireazei
        sumSets=0;
        MPI_Reduce(&endBig,&sumSets,1,MPI_INT,MPI_SUM,0,group); //FIRST REDUCE : SUM OF BIG, stelnw ola ta bigset gia na athroistoun sotn master
        MPI_Bcast(&sumSets,1,MPI_INT,0,group);
        if(oldSumSets==sumSets)
            checkIdentical=1;
        else
        {
            oldSumSets=sumSets;
            checkIdentical=0;
        }
        MPI_Bcast(&finalize,1,MPI_INT,0,group);//an o master apo to keepbigset k apo to count apofasisei oti teleiwsame mou stelnei 1, alliws 0 sunexizoume
        if(finalize==1)
        {
            float median=0;
            validation(median,partLength,size,numberPart,processId,group);

            return ;
        }
        MPI_Bcast(&keepBigSet,1,MPI_INT,0,group);//THIRD BROADCAST: Receive keepBigset boolean, edw lambanw an krataw to mikro i megalo set.
            //afou elaba ton keepbigset an eimai active krataw enan apo tous duo pinake small h big.. alliws den kanw tpt
            //edw antistoixa allazw tous pointers, k eksetazw an exw meinei xwris stoixeia tin opoia periptwsi sikwnw to dropoutflag k pio katw tha dilwsw na ginw inactive
        if(stillActive==1)
        {
            if(keepBigSet==1)
            {
                if(endBig==0)
                    dropoutflag=1;
                else
                {
                    arrayToUse=arrayBig;
                    elements=endBig;
                }
            }
            else if(keepBigSet==0)
            {
                if(endSmall==0)
                    dropoutflag=1;
                else
                {
                    arrayToUse=arraySmall;
                    elements=endSmall;
                }
            }
        }
        //edw einai ligo periploka grammeno, isws exei perita mesa alla, an eimai active k thelw na ginw inactive einai i prwti periptwsi, h deuteri einai eimai inactive hdh k i triti einai sunexizw dunamika
        if(dropoutflag==1 && stillActive==1)
        {
            MPI_Send(&dropoutflag,1,MPI_INT,0,1,group); //FIRST SEND : send active or not;
            stillActive=0;
        }
        else if(stillActive==0)
        {
            dropoutflag=-1;
            MPI_Send(&dropoutflag,1,MPI_INT,0,1,group); //FIRST SEND : send active or not;
        }
        else
        {
            dropoutflag=0;
            MPI_Send(&dropoutflag,1,MPI_INT,0,1,group); //FIRST SEND : send active or not;
        }
    }
}
/***                 ***/
void Medians(int processId,int noProcesses,int partLength,float *pointsArray,float *vantagepoints,int size,int depth,int dimension,float* medianArray)
{
  float *vpointDim;
  vpointDim=(float*)malloc(dimension*sizeof(float));

  float *distances;
  distances=(float*)malloc(partLength*sizeof(float));
  float *distances2;
  distances2=(float*)malloc(partLength*sizeof(float));

  int vpoint;
  MPI_Comm group;
  int color = processId / ( noProcesses / (1<<depth) ) ;
  MPI_Comm_split(MPI_COMM_WORLD, color, processId, &group);

  int local_rank, local_size;
  MPI_Comm_rank(group, &local_rank);
  MPI_Comm_size(group, &local_size);
  float median;



  /*---------------SELCECT VP---------------------*/
  if(local_rank==0)
  {
    srand(time(NULL));
    vpoint=(rand() % partLength);
    for(int i=0;i<dimension;i++)
      vpointDim[i]=pointsArray[vpoint+i];
      //printf("selected %f---%f -%d-\n",vpointDim[0],vpointDim[1],processId);
      MPI_Bcast(vpointDim,dimension,MPI_FLOAT,0,group);
  }else{
      MPI_Bcast(vpointDim,dimension,MPI_FLOAT,0,group);
  }
  /*---------------CALCULATE DIST FROM VP---------------------*/
  for(int i=0;i<partLength;i++)
       {
       distances[i]=distance(pointsArray,i,vpointDim,dimension);
       distances2[i]=distances[i];
       }

  /*---------------CALCULATE MEDIAN---------------------*/
  if(local_rank==0)
  {
      median=masterPart(local_size,local_rank,size/(1<<depth),partLength,distances,group);
  }
  else
  {
      slavePart(local_rank,partLength,distances,size/(1<<depth),group);
  }

  /*---------------SEND DATA---------------------*/
  if(local_rank==0)
  {
      if(processId==0)
      {

        medianArray[0]=median;
        for(int i=0;i<dimension;i++)
        vantagepoints[i]=vpointDim[i];

        for(int i=1;i<(1<<depth);i++)
        {
          MPI_Recv(&medianArray[i],1,MPI_FLOAT,i*noProcesses/(1<<depth),0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
          MPI_Recv(&vantagepoints[i*dimension],dimension,MPI_FLOAT,i*noProcesses/(1<<depth),0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        }
      }
      else
      {
        MPI_Send(&median,1,MPI_FLOAT,0,0,MPI_COMM_WORLD);
        MPI_Send(vpointDim,dimension,MPI_FLOAT,0,0,MPI_COMM_WORLD);
      }
  }
  MPI_Bcast(medianArray,(1<<depth),MPI_FLOAT,0,MPI_COMM_WORLD);
  MPI_Bcast(vantagepoints,(1<<depth)*dimension,MPI_FLOAT,0,MPI_COMM_WORLD);

  /*---------------SPLIT DATA---------------------*/

  float* bufferDistance;
  float* bufferPoints;
  if(local_rank==0){
      bufferDistance = (float *)malloc((size/(1<<depth))*sizeof(float));
      bufferPoints = (float *)malloc(( (size/(1<<depth) )*dimension )*sizeof(float));
      }


  MPI_Gather(pointsArray,partLength*dimension,MPI_FLOAT,bufferPoints,partLength*dimension,MPI_FLOAT,0,group);
  MPI_Gather(distances2,partLength,MPI_FLOAT,bufferDistance,partLength,MPI_FLOAT,0,group);


  if(local_rank==0){
      split(bufferPoints,bufferDistance,size/(1<<depth),median,dimension);
      }


  MPI_Scatter(bufferPoints,partLength*dimension,MPI_FLOAT,pointsArray,partLength*dimension,MPI_FLOAT,0,group);
  //MPI_Scatter(bufferDistance,partLength,MPI_FLOAT,distances,partLength,MPI_FLOAT,0,group);

  /*  //VALIDATION
  float dist;
  for(int i=0;i<partLength;i++){
      if(local_rank==1)
           {
             printf("%f--- %f  ",pointsArray[i*2],pointsArray[i*2+1]);
             dist=distance(pointsArray,i,vpointDim,dimension);
             printf("%f---%d\n",dist,processId);
           }
  }
  */

  if(local_rank==0){
  free(bufferDistance);
  free(bufferPoints);

  }
  free(distances);
  free(distances2);

  MPI_Comm_free(&group);


}

void split(float *pointsArray,float* elements,int size,float median,int dimension)
{
  int left=0;
  int right=size-1;
  for(int i=0;i<size;i++)



  while(left<right)
  {
      while(elements[left]<=median)
      {
          left++;
          if(left>=size)
          {
              break;
          }
      }
      while(elements[right]>median)
      {
          right--;
          if(right<0)
          {
              break;
          }
      }
      if(left<right)
      {
          //printf("swap %d %d",left,right);
          swap_values(elements,left,right);
          swap_point(pointsArray,left,right,dimension);
      }
  }
}

void search(struct Node* tree,float* point,int k,int dimension,std::priority_queue<HeapItem>*heap,int pros,int* should){

  int t;
  float dist=0;
  for(int i=0;i<dimension;i++){
      dist=dist+(tree->vpoint[i]-point[i])*(tree->vpoint[i]-point[i]);
    }
  dist=sqrtf(dist);

  if(tree->leaf==false){


    if(dist<=tree->radius){
        t=pros;
        t=t<<1;
      if(tree->left!=NULL){

          search( tree->left, point, k,dimension,heap,t,should);
        }else{

            should[t]=1;
          }
      if ( tau>=tree->radius-dist) {
        t=pros;
        t=t<<1;
        t=t+1;
        if(tree->right!=NULL){

          search( tree->right, point, k,dimension,heap,t,should);
      }else{

        should[t]=1;
      }
    }
    }
    else{
      t=pros;
      t=t<<1;
      t=t+1;
      if(tree->right!=NULL){


      search( tree->right, point, k,dimension,heap,t,should);}
      else{

        should[t]=1;
      }
      if ( tau>=dist-tree->radius) {
        t=pros;
        t=t<<1;
        if(tree->left!=NULL){


        search( tree->left, point, k,dimension,heap,t,should);
      }else{

        should[t]=1;
      }}

    }

  }
  else{

    if ( dist <= tau ) {

            if ( heap->size() == k ) heap->pop();
            heap->push(HeapItem(tree->vpoint,dist));
            if ( heap->size() == k ){ tau = heap->top().dist;

          }
        }

  return ;
}
}


void swap(float *xp, float *yp)
{
    float temp = *xp;
    *xp = *yp;
    *yp = temp;
}



void packHeap(std::priority_queue<HeapItem>*heap,float* arrayPoints,float *arrayRadius,int k,int dimension){
    for(int i=0;i<k;i++){
      arrayRadius[i]=heap->top().dist;

      for(int j=0;j<dimension;j++){
        arrayPoints[i*dimension+j]=heap->top().p[j];
      }
      heap->pop();
    }
  }
void unpackHeap(std::priority_queue<HeapItem>*heap,float* arrayPoints,float *arrayRadius,int k,int dimension){

    for(int i=0;i<k;i++){
    heap->push(HeapItem(arrayPoints+i*dimension,arrayRadius[i]));
  }


}

void sendHeap(priority_queue<HeapItem> *heap,float* heapDist,float *heapPoints,int* should,int*tempshould,
              float* heapDistBuffer,float *heapPointsBuffer,struct Node *tree,float *point,
              int processFlag,int target,int processId,int noProcesses,int k,int dimension,MPI_Comm group)
{

  int sendFlag=processFlag;
  int receiveFlag;
  float *pointBuffer=(float*)malloc(dimension*sizeof(float));

  for(int i=0;i<noProcesses;i++) tempshould[i]=0;
  priority_queue<HeapItem> heapBuffer;

  MPI_Send(&sendFlag,1,MPI_INT ,target,0,group);
  if(sendFlag==1){



    MPI_Send(heapDist,k,MPI_FLOAT,target,0,group);
    MPI_Send(heapPoints,k*dimension,MPI_FLOAT ,target,0,group);
    MPI_Send(point,dimension,MPI_FLOAT ,target,0,group);
  }
  MPI_Recv(&receiveFlag,1,MPI_INT,target,0,group, MPI_STATUS_IGNORE);
  if(receiveFlag==1){
     MPI_Recv(heapDistBuffer,k,MPI_FLOAT,target,0,group, MPI_STATUS_IGNORE);
     MPI_Recv(heapPointsBuffer,k*dimension,MPI_FLOAT,target,0,group, MPI_STATUS_IGNORE);
     MPI_Recv(pointBuffer,dimension,MPI_FLOAT,target,0,group, MPI_STATUS_IGNORE);

     unpackHeap(&heapBuffer,heapPointsBuffer,heapDistBuffer,k,dimension);
     tau=heapBuffer.top().dist;
     printf("received tau:%f heap:%d p: %d from:%d\n",tau,heapBuffer.size() ,processId,target);
     printf("tau: %f\n",tau );
     search(tree,pointBuffer,k,dimension,&heapBuffer,0,tempshould);
     packHeap(&heapBuffer,heapPointsBuffer,heapDistBuffer,k,dimension);
  }
  if(sendFlag==1){

    MPI_Recv(heapDist,k,MPI_FLOAT,target,0,group, MPI_STATUS_IGNORE);
    MPI_Recv(heapPoints,k*dimension,MPI_FLOAT,target,0,group, MPI_STATUS_IGNORE);
    MPI_Recv(should,noProcesses,MPI_FLOAT,target,0,group, MPI_STATUS_IGNORE);

  }
  if(receiveFlag==1){
    MPI_Send(heapDistBuffer,k,MPI_FLOAT,target,0,group);
    MPI_Send(heapPointsBuffer,k*dimension,MPI_FLOAT ,target,0,group);
    MPI_Send(tempshould,noProcesses,MPI_INT,target,0,group);
    printf("send tau:%f heap:%d p: %d from:%d\n",tau,heapBuffer.size() ,processId,target);

  }

}

void receiveHeap(priority_queue<HeapItem> *heap,float* heapDist,float *heapPoints,int *should,int*tempshould,
                 float* heapDistBuffer,float *heapPointsBuffer,struct Node *tree,float *point,
                 int processFlag,int target,int processId,int noProcesses,int k,int dimension,MPI_Comm group)
{
  float *pointBuffer=(float*)malloc(dimension*sizeof(float));
  int sendFlag=processFlag;
  int receiveFlag;

  for(int i=0;i<noProcesses;i++) tempshould[i]=0;
  priority_queue<HeapItem> heapBuffer;

  MPI_Recv(&receiveFlag,1,MPI_INT,target,0,group, MPI_STATUS_IGNORE);

  if(receiveFlag==1){

     MPI_Recv(heapDistBuffer,k,MPI_FLOAT,target,0,group, MPI_STATUS_IGNORE);
     MPI_Recv(heapPointsBuffer,k*dimension,MPI_FLOAT,target,0,group, MPI_STATUS_IGNORE);
     MPI_Recv(pointBuffer,dimension,MPI_FLOAT,target,0,group, MPI_STATUS_IGNORE);

  }
  MPI_Send(&sendFlag,1,MPI_INT ,target,0,group);
  if(sendFlag==1){

    MPI_Send(heapDist,k,MPI_FLOAT,target,0,group);
    MPI_Send(heapPoints,k*dimension,MPI_FLOAT ,target,0,group);
    MPI_Send(point,dimension,MPI_FLOAT ,target,0,group);
  }
  if(receiveFlag==1){
  unpackHeap(&heapBuffer,heapPointsBuffer,heapDistBuffer,k,dimension);
  tau=heapBuffer.top().dist;
  printf("received tau:%f heap:%d p: %d from:%d\n",tau,heapBuffer.size() ,processId,target);
  search(tree,pointBuffer,k,dimension,&heapBuffer,0,tempshould);
  packHeap(&heapBuffer,heapPointsBuffer,heapDistBuffer,k,dimension);

  MPI_Send(heapDistBuffer,k,MPI_FLOAT,target,0,group);
  MPI_Send(heapPointsBuffer,k*dimension,MPI_FLOAT ,target,0,group);
  MPI_Send(tempshould,noProcesses,MPI_INT,target,0,group);
  printf("send tau:%f heap:%d p: %d from:%d\n",tau,heapBuffer.size() ,processId,target);
  }
  if(sendFlag==1){
     MPI_Recv(heapDist,k,MPI_FLOAT,target,0,group, MPI_STATUS_IGNORE);
     MPI_Recv(heapPoints,k*dimension,MPI_FLOAT,target,0,group, MPI_STATUS_IGNORE);
     MPI_Recv(should,noProcesses,MPI_FLOAT,target,0,group, MPI_STATUS_IGNORE);

  }

}

void findNeighbors(priority_queue<HeapItem> *heap,float* heapdist,float *heapPoints,int *should,int*tempshould,
                 float* heapDistBuffer,float *heapPointsBuffer,struct Node *tree,float *point,
                 int processId,int noProcesses,int k,int dimension)
{

  int q=log2(noProcesses);
  int processFlag;
  int target;

    for(int level=q-1;level>=0;level--){

    int groupsize = ( noProcesses / (1<<level) ) ;
    int color = processId / groupsize ;

    MPI_Comm group;
    MPI_Comm_split(MPI_COMM_WORLD, color, processId, &group);
    int local_rank, local_size;
    MPI_Comm_rank(group, &local_rank);
    MPI_Comm_size(group, &local_size);

    int maplocation = groupsize * color ;

    //printf("processId:%d  local_rank:%d  \n",processId,local_rank);


    //number of communications per phase
    for(int i=0;i<groupsize/2;i++){
      if(local_rank<groupsize/2){
      target=groupsize/2+(local_rank+i)%(groupsize/2);
      processFlag=should[maplocation+target];
      sendHeap(heap,heapdist,heapPoints,should,tempshould,heapDistBuffer,heapPointsBuffer,tree,point,
                      processFlag,target,processId,noProcesses,k,dimension,group);


    }else{
      target=(local_rank-i)%(groupsize/2);
      processFlag=should[maplocation+target];
     receiveHeap(heap,heapdist,heapPoints,should,tempshould,heapDistBuffer,heapPointsBuffer,tree,point,
                       processFlag,target,processId,noProcesses,k,dimension,group);
    }

  }

    MPI_Comm_free(&group);
  }

}

/*****MAIN!!!!!!!!!!*****/
int main (int argc, char **argv)
{
    int processId,noProcesses,size,partLength;
    int dimension=2;
    float *pointsArray;
    float *numberPart;
    float median;
    size=(1<<atoi(argv[1]));
    MPI_Init (&argc, &argv);	/* starts MPI */
    MPI_Comm_rank (MPI_COMM_WORLD, &processId);	/* get current process id */
    MPI_Comm_size (MPI_COMM_WORLD, &noProcesses);	/* get number of processes */
    if(processId==0)
    {
        printf("size: %d processes: %d\n",size,noProcesses);
        if(noProcesses>1)
        {
          if(size%noProcesses==0)
              partLength=(size/noProcesses);
          else
              partLength=(size/noProcesses)+1;
            sendLengths(size,noProcesses);
            pointsArray=(float*)malloc(partLength*dimension*sizeof(float));
            generatePoints(pointsArray,partLength,processId,dimension);

        }
        else
        {
          pointsArray=(float*)malloc(partLength*dimension*sizeof(float));
          generatePoints(pointsArray,partLength,processId,dimension);



                    struct Node *tree;

                  //  printf("----Creating trees----\n");

                    tree=createtreelocal(pointsArray, partLength, dimension,processId);
                    //SEARCH YOUR OWN

                    int* should;
                    int k=4;
                    should=(int*)malloc(sizeof(int)*noProcesses*partLength);
                    priority_queue<HeapItem> heap ;
                    float* heapdist=(float *)malloc(sizeof(float)*k);
                    float* heapPoints=(float *)malloc(sizeof(float)*dimension*k);
                    int point=0;
                       *should=0;
                        tau = numeric_limits<float>::max();
                        search(tree,pointsArray+dimension*point,k,dimension,&heap,0,should);
                        packHeap(&heap,heapPoints,heapdist,k,dimension);


            MPI_Finalize();
            return 0;
        }
    }
    else
    {
        MPI_Recv(&partLength,1,MPI_INT,0,1,MPI_COMM_WORLD,&Stat);
        pointsArray=(float*)malloc(partLength*dimension*sizeof(float));
        generatePoints(pointsArray,partLength,processId,dimension);

    }




    float* vantagepoints;
    float* medianArray;
    int q = log2(noProcesses);
    struct Node tree;
    struct Node *lvl;
    struct Node *dad;
    int depth;
    depth=q;
    medianArray=(float*)malloc((1<<depth)*sizeof(float));
    vantagepoints=(float*)malloc((1<<depth)*dimension*sizeof(float));
    struct Node* localroot;
    struct Node** temp;
    dad=&tree;

    for(depth=0;depth<q;depth++)
    {

      Medians(processId,noProcesses,partLength,pointsArray,vantagepoints,size,depth,dimension,medianArray);



          lvl=(struct Node*)malloc((1<<depth)*2*sizeof(struct Node));
          for(int i=0;i<(1<<depth);i++){
            dad[i].vpoint=(float *)malloc(dimension*sizeof(float));
            for(int j=0;j<dimension;j++)
              {dad[i].vpoint[j]=vantagepoints[i*dimension+j];}
            dad[i].leaf=false;
            dad[i].radius=medianArray[i];
            dad[i].left=&lvl[i*2];
            dad[i].right=&lvl[i*2+1];
            if(depth==q-1){
              dad[i].right=NULL;
              dad[i].left=NULL;
              if(i==floor(processId/2)){
              if(processId%2==0){
              temp=&dad[i].left;
              }else{
              temp=&dad[i].right;
            }}
    }
    }
    dad=lvl;
    MPI_Barrier(MPI_COMM_WORLD);
  //}

}


    int k=16;


    localroot=createtreelocal(pointsArray, partLength, dimension,processId);
    *temp=localroot;

    //printf("---------------Tree------------ \n");
    //print2DUtil(&tree,0,2);
    //printf("\n" );
    //printf("---------------Tree LOCAL------------\n \n");
    //print2DUtil(localroot,0,2);


    float* point;
    point=(float *)malloc(sizeof(float)*dimension);
    for(int j=0;j<dimension;j++)
      {point[j]=pointsArray[5*dimension+j];      }


    int* should;
    should=(int*)malloc(sizeof(int)*noProcesses);

    for(int j=0;j<noProcesses;j++)
      {
        should[j]=0;
      }

      priority_queue<HeapItem> heap ;

      tau = numeric_limits<float>::max();

      search(&tree,point,k,dimension,&heap,0,should);


      if(processId==0){
      printf("--PROCESS:%d--\n",processId);
      for(int j=0;j<noProcesses;j++)
         {
         printf("%-d-",should[j]);
          }
      printf("\n");
      }




       float* d= (float*)malloc(size*sizeof(float));
       for(int i=0;i<partLength;i++)
            {
            d[i]=distance(pointsArray,i,point,dimension);
            }

       for (int i = 0; i < partLength-1; i++){
          // Last i elements are already in place
          for (int j = 0; j < partLength-i-1; j++)
              if (d[j] > d[j+1])
                 swap(&d[j], &d[j+1]);
        }


        float* heapdist=(float *)malloc(sizeof(float)*k);
        float* heapPoints=(float *)malloc(sizeof(float)*dimension*k);

        float* heapDistBuffer=(float *)malloc(sizeof(float)*k);
        float* heapPointsBuffer=(float *)malloc(sizeof(float)*dimension*k);
        int *tempshould=(int*)malloc(noProcesses*sizeof(int));

        packHeap(&heap,heapPoints,heapdist,k,dimension);

        if(processId==0){
        for(int j=0;j<k;j++)
           {
             printf(" %f--%f \n",d[k-1-j],heapdist[j]);
            }
        }



        findNeighbors(&heap,heapdist,heapPoints,should,tempshould,heapDistBuffer,heapPointsBuffer,&tree,point
                      ,processId,noProcesses,k,dimension);




       printf("process: %d FINISH \n", processId );



        //tau = numeric_limits<float>::max();
      //Search(localroot,point,k,dimension,&heap);

     float *bufferPoints;
     //unpackHeap(&heap,heapPoints,heapdist,k,dimension);

      if(processId==0){
          bufferPoints = (float *)malloc(dimension*size*sizeof(float));
          }

      MPI_Gather(pointsArray,partLength*dimension,MPI_FLOAT,bufferPoints,partLength*dimension,MPI_FLOAT,0,MPI_COMM_WORLD);

      if(processId==0){
      for(int i=0;i<size;i++)
           {
           d[i]=distance(bufferPoints,i,point,dimension);
           }

      for (int i = 0; i < size-1; i++){
         // Last i elements are already in place
         for (int j = 0; j < size-i-1; j++)
             if (d[j] > d[j+1])
                swap(&d[j], &d[j+1]);
       }
       for(int j=0;j<k;j++)
          {
         printf(" %f--%f \n",d[k-1-j],heapdist[j]);
           }

         free(bufferPoints);
     }









    MPI_Finalize();



    return 0;


}


void print2DUtil(struct Node *root, int space,int dimension)
{
    // Base case
    if (root->leaf == true){
       printf("       " );
    for(int i=0;i<dimension;i++){
       printf("leaf-%d->%f",i,root->vpoint[i]);}
    printf("\n");

        return;
}
    // Increase distance between levels
    space += COUNT;

    // Process right child first
    print2DUtil(root->right, space,dimension);

    // Print current node after space
    // count
    printf("\n");
    for (int i = COUNT; i < space; i++)
        {printf(" ");}
    printf("R %f -", root->radius);
    for(int i=0;i<dimension;i++){
       printf("P-%d->%f",i,root->vpoint[i]);}


    // Process left child
    print2DUtil(root->left, space,dimension);
}
/*
 Given a binary search tree, print out
 its data elements in increasing
 sorted order.
*/
void printTree(struct Node* node,int dimension) {
  if (node->leaf == true) return;

  printTree(node->left,dimension);
  printf("Radius %f ", node->radius);
  for(int i=0;i<dimension;i++){
    printf("Point %d ->%f",i,node->vpoint[i]);}
  printf("\n" );
  printTree(node->right,dimension);
}


struct Node* createtreelocal(float *pointsArray,int partLength,int dimension,int processId){
  //printf("============================NODE============================\n");
  float* distances;
  distances=(float*)malloc(partLength*sizeof(float));

  float *vpointDim;
  vpointDim=(float *)malloc(sizeof(float)*dimension);
  float median;
  /*---------------RANDOM VP---------------------*/
  srand(time(NULL));
  int vpoint=(rand() % partLength);
  for(int i=0;i<dimension;i++){
    vpointDim[i]=pointsArray[vpoint*dimension+i];
  }
  /*---------------CALCULATE DIST FROM VP---------------------*/
  for(int i=0;i<partLength;i++)
       {
       distances[i]=distance(pointsArray,i,vpointDim,dimension);
       }
  /*---------------CALCULATE MEDIAN---------------------*/
  if(partLength==2){
  median=(distances[0]+distances[1])/2;
  }else{
  median=selection(distances,partLength);
  }

  for(int i=0;i<partLength;i++)
       {
       distances[i]=distance(pointsArray,i,vpointDim,dimension);
       }


  split(pointsArray,distances,partLength,median,dimension);

  struct Node *Nodev;
  Nodev=(struct Node *)malloc(sizeof(struct Node));

  Nodev->radius=median;

  Nodev->vpoint=(float *)malloc(sizeof(float)*dimension);
  for(int i=0;i<dimension;i++){
  Nodev->vpoint[i]=vpointDim[i];}
  Nodev->leaf=true;

  if(partLength>1){
  Nodev->leaf=false;
   int offset=(partLength/2)*dimension;
   int k;
   if(partLength%2==0){
     k=partLength/2-1;
   }else{k=partLength/2;}
   //printf("============================LEFT============================\n");
  Nodev->left=createtreelocal(pointsArray,partLength/2, dimension,processId);
  //printf("========================================================\n");
  Nodev->right=createtreelocal(pointsArray+offset,partLength/2,dimension,processId);




}else{
//printf("--createtreelocal-- %f,%f \n",Nodev->vpoint[0],Nodev->vpoint[1]);

  }



  //printf("====>%d----%d\n", processId,partLength);
  return Nodev;
}


/****Partitions the Array into larger and smaller than the pivot values****/
void partition (float *array, int elements, float pivot, float **arraysmall, float **arraybig, int *endsmall, int *endbig)
{
    int right=elements-1;
    int left=0;
    int pos;
    if(elements==1)
    {
        if(pivot>array[0])
        {
            *endsmall=1;  //One value in the small part
            *endbig=0;   //Zero on the big one
            *arraysmall=array;   //There is no big array therefore NULL value
            *arraybig=NULL;
        }
        else if(pivot<=array[0])
        {
            *endsmall=0;    //The exact opposite of the above actions.
            *endbig=1;
            *arraysmall=NULL;
            *arraybig=array;
        }
    }
    else if(elements>1)
    {
        while(left<right)
        {
            while(array[left]<pivot)
            {
                left++;
                if(left>=elements)
                {
                    break;
                }
            }
            while(array[right]>=pivot)
            {
                right--;
                if(right<0)
                {
                    break;
                }
            }
            if(left<right)
            {
                swap_values(array,left,right);
            }
        }
        pos=right;
        if(pos<0)                   //Arrange the arrays so that they are split into two smaller ones.
        {                               //One containing the small ones. And one the big ones.
            *arraysmall=NULL;           //However these arrays are virtual meaning that we only save the pointer values of the beging and end
        }                               //of the "real" one.
        else
        {
            *arraysmall=array;
        }
        *endsmall=pos+1;
        *arraybig=&array[pos+1];
        *endbig=elements-pos-1;
    }
}


/***==============================================***/
/***==============================================***/
/***=============SERIAL SELECTION==============***/
/***==============================================***/
/***==============================================***/

float selection(float *array,int number)
{
    float *arraybig;
    float *arraysmall;
    int endsmall=0;
    int endbig=0;
    float *arraytobeused;
    int i;
    int counter=0;
    int k;
    float pivot;
    float median;
    k=(int)number/2+1;
    arraytobeused=array;
    for(;;)
    {
        pivot=arraytobeused[rand() % number];
        partition(arraytobeused,number,pivot,&arraysmall,&arraybig,&endsmall,&endbig);
        if(endbig>k)
        {
            number=endbig;
            arraytobeused=arraybig;
            for(i=0;i<endbig;i++)
            {
                if(pivot==arraybig[i])
                    counter++;
                else
                    break;
            }
            if(counter==endbig)
            {
                median=arraybig[0];
                break;
            }
            else
                counter=0;
            //end of count equals
        }
        else if(endbig<k)
        {
            number=endsmall;
            arraytobeused=arraysmall;
            k=k-endbig;
        }
        else
        {
            median=pivot;
            break;
        }
    }
    return median;
}
void partitionPoints (float*pointsArray,int dimension,float *array, int elements, float pivot, float **arraysmall, float **arraybig, int *endsmall, int *endbig)
{
    int right=elements-1;
    int left=0;
    int pos;
    if(elements==1)
    {
        if(pivot>array[0])
        {
            *endsmall=1;  //One value in the small part
            *endbig=0;   //Zero on the big one
            *arraysmall=array;   //There is no big array therefore NULL value
            *arraybig=NULL;
        }
        else if(pivot<=array[0])
        {
            *endsmall=0;    //The exact opposite of the above actions.
            *endbig=1;
            *arraysmall=NULL;
            *arraybig=array;
        }
    }
    else if(elements>1)
    {
        while(left<right)
        {
            while(array[left]<pivot)
            {
                left++;
                if(left>=elements)
                {
                    break;
                }
            }
            while(array[right]>=pivot)
            {
                right--;
                if(right<0)
                {
                    break;
                }
            }
            if(left<right)
            {
                swap_values(array,left,right);
                swap_point(pointsArray,left,right,dimension);
            }
        }
        pos=right;
        if(pos<0)                   //Arrange the arrays so that they are split into two smaller ones.
        {                               //One containing the small ones. And one the big ones.
            *arraysmall=NULL;           //However these arrays are virtual meaning that we only save the pointer values of the beging and end
        }                               //of the "real" one.
        else
        {
            *arraysmall=array;
        }
        *endsmall=pos+1;
        *arraybig=&array[pos+1];
        *endbig=elements-pos-1;
    }
}
float selectionpPoints(float*pointsArray,int dimension,float *array,int number)
{
    float *arraybig;
    float *arraysmall;
    int endsmall=0;
    int endbig=0;
    float *arraytobeused;
    int i;
    int counter=0;
    int k;
    float pivot;
    float median;
    k=(int)number/2+1;
    arraytobeused=array;
    for(;;)
    {
        pivot=arraytobeused[rand() % number];
        partitionPoints(pointsArray,dimension,arraytobeused,number,pivot,&arraysmall,&arraybig,&endsmall,&endbig);
        if(endbig>k)
        {
            number=endbig;
            arraytobeused=arraybig;
            for(i=0;i<endbig;i++)
            {
                if(pivot==arraybig[i])
                    counter++;
                else
                    break;
            }
            if(counter==endbig)
            {
                median=arraybig[0];
                break;
            }
            else
                counter=0;
            //end of count equals
        }
        else if(endbig<k)
        {
            number=endsmall;
            arraytobeused=arraysmall;
            k=k-endbig;
        }
        else
        {
            median=pivot;
            break;
        }
    }

    return median;
}
