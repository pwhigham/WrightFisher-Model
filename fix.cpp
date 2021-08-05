/*
Network Model
Hill & Robertson Paper (1966)
Space implicit. 
9nd Sept., 2018
Author: P A Whigham

Updated: 22/2/2020
This is a single locus model - otherwise same structure.
Used for theoretical models and simplified behaviour.
This is NOT Hill & Robertson.


Updated: 10/10/2018

Changed to have all inputs defined (but still fixed space)
This will allow both Fig 1 and other measures to be done

We currently have defined: "fixed" networks:
 
Ring, SM , SF 
 
Setup that you read in the space - and can just define the folder where the 1 or more
spaces are loaded.  This allows many SM & SF spaces if required, and a single one for 
the ring.

C compiler:  bcc32c <file>
With space it is approx x3 in terms of time, but still acceptable.
*/
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <ctype.h>


/* 0 == AB
   1 == aB
   2 == Ab
   3 == ab
*/
// Here we define the parameters of the model that are "fixed"
// Allows us to use main memory vs stack

#define N 200 // Number of individuals

#define SMALLFITNESS 0.0000001
#define FALSE 0
#define TRUE 1

static int pop[N];  //  Population with 1 locus
static int nextpop[N];  // Next generation
static float fitness[N];  // Fitness for N individuals
static float cumsum[N];   // Holds calculation for cumulative sum/sum
static int fixedResult[2]; // <fix1/0 gens>
                           // where fix* = TRUE/FALSE, 
						   // gens is number generations to fixation
static float mutantfitness; // Set by Nalpha as starting parameter

/*************************************************************************************/
/*   NETWORK STRUCTURE 
***************************************************************************************/

static float network[N][N];  
					   
/************************** UTILITIES ***************************/

/*****************************************************************
********* trimstr - trim spaces at the end of a string
*****************************************************************/ 
char *trimstr(char *s)
{
    if(s) { /* Don't forget to check for NULL! */
        while(*s && isspace(*s))
            ++s;
        if(*s) {
            register char *p = s;
            while(*p)
                ++p;
            do {
                --p;
            } while((p != s) && isspace(*p));
            *(p + 1) = '\0';
        }
    }
    return(s);
}
/******************************************************************
*********** intToString - convert an integer to the string
******************************************************************/
char* intToString(int x)
{
	int length = snprintf( NULL, 0, "%d",x);
	char* str = (char *) malloc( length + 1 );
	snprintf( str, length + 1, "%d", x );
	return(str);
}
/*
randrange(lower,upper)
Generate random integers in range lower to upper inclusive.
i.e. randrange(0,5)  generates random samples of 0,1,2,3,4,5

*/
int randrange(int lower, int upper) 
{ 
	return (rand() % 
           (upper - lower + 1)) + lower; 
} 

/**************************** NETWORK FUNCTIONS ******************/
  
/*   loadNetwork(fname) Load a network into network array from file
********************************************************************/

void loadNetwork(char* fname)
{
	FILE *net;
    net = fopen(fname, "r");
	
	if (net==NULL)
	{
		printf("FAILED TO OPEN NETWORK: %s\n",fname);
		exit(0);
	}
	for(int i=0;i<N;++i)
		for(int j=0;j<N;++j)
		{
			if (fscanf(net, "%f", &network[i][j])!=1) printf("FAILED\n");
		}
	fclose(net);
}
/*   printNetwork() - print current network to stdoutj
***************************************************************************/
void printNetwork()
{
	for(int i=0;i<N;++i)
	{
		for(int j=0;j<N;++j)
			printf("%f",network[i][j]);
		printf("\n");
	}
}

/*
NOTE: This is only required when setting p or q to zero
Since we don't do this in the model runs, it should make no difference
to the result.  
*/
double rand01()
{
	double res=(((double)rand())/((double)(RAND_MAX)));
	return((res==(double)0.0)?(double)0.0000000001:res); 
}
/*
Initialisation based on p
IF p is 0, just put one mutate (1) randomly
in the population

pop - a 1 represents the mutate, 0 is the normal (fitness 1)

*/

void initpop(double p)
{
	if (p==0.0)
	{
		for (int i=0;i<N;++i)
			pop[i] = 0;
		pop[randrange(0,N-1)] = 1;
	}
	else
	{
		for (int i=0;i<N;++i)
		{
			pop[i] = (rand01() > p)?0:1;
		}
	}
}

void calc_fitness()
{
//	for (int i=0;i<N;++i)
//		fitness[i] = 1;   /* DRIFT */
	
	
	for (int i=0;i<N;++i)
		fitness[i] = pop[i]==0?1:mutantfitness;
}
/*
*  fixed
   Is population fixed? All 0 or 1's	
*/
int fixed()
{
	int first=pop[0];
	
	for(int i=1;i<N;++i)
	{
		if (pop[i]!=first) return(FALSE);
	}
	return(TRUE);
}
void print_current_pop()
{
	for (int i=0;i<N;++i)
		printf("<IND %d> <L %d>: %f\n",i,pop[i],fitness[i]);
}
/*
Set the normalised cumulative sum of fitness for location <locn> 
Uses the network[locn][N] - if 1, then include in cumsum running total
*/
void calc_cumsum(int locn)
{
	register float total=0.0;
	
	for(register int i=0;i<N;++i)
	{
		if (network[locn][i] > 0.0)
		{
			total += fitness[i];   // /* ORIGINAL MODEL */
		}
		cumsum[i] = total; /* Update or set to what was total before */
	}
	
	for (int i=0;i<N;++i)
		cumsum[i]/=total;  // Normalise by the sum of fitness...
}
/*
Select parent using the cumsum that has been set
*/
int select_parent()
{
	float rnd = ((float)rand()/(float)(RAND_MAX));
	float sum=0.0;
	
	for(int i=0;i<N;++i)
	{
		if (cumsum[i] >= rnd) return(i);
	}
	printf("Failed select parent \n");
	print_current_pop();
	exit(0);
}
/*
Create next generation
Uses nextpop to hold, and then copy back to pop at the end.
OUT: Changes pop
*/
void next_generation()
{
	int p1;  // parent
		
	for (int i=0;i<N;++i)  // for each location/individual
	{
		calc_cumsum(i); // Calculate cumulative fitness for location <i>
		p1 = select_parent(); // Select parents

		// and create the child....		
		nextpop[i] = pop[p1];
	}
	
	// and finally copy the next pop back to pop....
	for(int i=0;i<N;++i)
	{
		pop[i] = nextpop[i];
	}
}
/*
Run until population fixes
Return data on final population - #gens, final locus values

*/
void popfixed(float p)
{
	int gens=0;
	int fix=FALSE;  
	
	initpop((double)p);  

	while(TRUE)
	{	
		if (fixed()) break;
		
		calc_fitness();  // calculate fitness for current population
						// Fitness table has been previously created.
						
		next_generation();
		gens+=1;  
	}
	fixedResult[0]=pop[1];  
	fixedResult[1]=gens;
}

/*********************************************************************
*** Network rewiring model.
***********************************************************************/

int main(int argc, char *argv[])
{
	int numruns,
		numspaces;  // number of spaces in the network folder
					// these are numbered 1 .. numspaces
	int writeheader=1;
	float p;
    float s;
	
	char netfname[512]; // pathname for folder holding fixed spaces
	char *networkpath;

	srand((unsigned int)time(NULL));	// initialise random number generator
	
	if (argc<7)
	{
		printf("<p float> <s float> <numruns int> \n <network folder path> <numspaces int> <writeheader int> \n");
		printf("N = %d\n", N);
		printf("\nSingle Locus Model - Nialpha gives selection added value (i.e. 1+s)\n");
		return(1);
	}
	sscanf (argv[1],"%f",&p);
	sscanf (argv[2],"%f",&s);
	sscanf(argv[3],"%d",&numruns);
	networkpath = argv[4];
	sscanf (argv[5],"%d",&numspaces);
	sscanf(argv[6],"%d",&writeheader);
	
	srand((unsigned int)time(NULL));	// initialise random number generator
		
	mutantfitness = 1.0 + s; // Set mutant fitness
	
	if (writeheader) printf("N,s,p,fix,gens\n");

	if (numspaces==1) /* Can just load space once - assume 1.txt */
	{
		sprintf(netfname,"%s\\\\1.txt\n",networkpath);
		loadNetwork(trimstr(netfname));		
	}
	for (int runs=0;runs < numruns; ++runs)
	{			
		char* strexp = intToString((runs%numspaces)+1);
			
		if (numspaces > 1)
		{
			sprintf(netfname,"%s\\\\%s.txt\n",networkpath,strexp);
			loadNetwork(trimstr(netfname));
		}
		popfixed(p);
		printf("%d,%.2f,%.2f,%d,%d\n",
				N,(float)s, (float)p,fixedResult[0],fixedResult[1]);
	}
}
