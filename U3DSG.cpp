//my Navier Stokes solution


#include<stdio.h>
#include<stdlib.h>
#include<assert.h>
#include<time.h>
#include<math.h>
#include <algorithm>
#include <cstring>
#include"files/variables.h"
#include"files/U3DSG.h"
#include"files/preProcessing.h"
#include"files/postProcessing.h"
#include"files/simple.h"
#include"files/simpler.h"
#include"files/schemes.h"
#include"files/deferredCorrection.h"



int main()
{

	defaultInput(); 	   	//Inputting the values required

	dynamicAllocate();		//dynamically allocates the arrays

	zeroTimeValues();   	//Gives initial values to all the variables and fields at time 0

	writeMesh();	      	//writes internal mesh to mesh.csv

	stepSize();

	startCPUClock();



	while(simulationTime<=endTime)         	//time loop
 	{
	  initialize();

	  boundaryConditions();     			//updates the values of the boundary cells

	  ghostCells();

	  simple(); 		    				//use simple or simpler algorithm

	  makeVelocityVectors();    			//Interpolation of the velocity vectors back from staggerd cell centers to the pressure cell centers

	  if(t % writeOutputSteps == 0)
	  writeOutput();           		 		//Outputs values to files

	  shift();		    					//shifts the values to l=0

	  courantNumber();	    				//calculates courant number

	  incrementTime();	    				//increments time variables

	  finalize();

	}


}
