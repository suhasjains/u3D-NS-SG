//Declaration of all functions


#ifndef NAVIERSTOKESSOLVER_H
#define NAVIERSTOKESSOLVER_H

    void defaultInput(); 	   	//Inputting the values required

	void dynamicAllocate();		//dynamically allocates the arrays

	void zeroTimeValues();   	//Gives initial values to all the variables and fields at time 0

	void writeMesh();	      	//writes internal mesh to mesh.csv

	void stepSize();

	void startCPUClock();

	void initialize();

    void boundaryConditions();     //updates the values of the boundary cells

	void ghostCells();

    void simple(); 		    //simple algorithm
    
    void simpler(); 		    //simpler algorithm
    
    void makeVelocityVectors();    //Interpolation of the velocity vectors back from staggerd cell centers to the pressure cell centers

    void writeOutput();            //Outputs values to files

    void shift();		    //shifts the values to l=0

    void courantNumber();	    //calculates courant number

    void incrementTime();	    //increments time variables

    void finalize();
    
    void Thomas();
    
    void pseudoVelocity();
    
    void pEqn();

    void uEqnCoeff();

    void vEqnCoeff();

    void wEqnCoeff();

    void pCorrEqnCoeff();

    void pEqnCoeff();
    
    void upwindU();
	
	void upwindV();
	
	void upwindW();
    
    void hybridU();
	
	void hybridV();
	
	void hybridW();
	
	void powerLawU();
	
	void powerLawV();
	
	void powerLawW();
	
	void quickU();
	 
	void quickV();
	
	void quickW();
	
    void massResidual();



#endif

