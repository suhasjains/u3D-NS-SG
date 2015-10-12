//************************************************************
//SOURCE CODE FUNCTIONS
//************************************************************

//SIMPLER algorithm loop
void simpler()
{

	printf("Starting SIMPLER loop\n");
	printf("Initial mass Residual: %e\n",maxMassResidual);


	for(simplerVar=0;simplerVar<nSimplerLoops;simplerVar++)
 {

	nLoops();

	for(nonLinear=0;nonLinear<nVelocityLoops;nonLinear++)
	{

		pseudoVelocity();	//finding pseudo velocities

		pEqn();			//finding absolute pressures

		uEqn();
		vEqn();			//finding initial velocities
		wEqn();


		massResidual();		//finding intermediate max mass residual
		
		nMassResidual=maxMassResidual*10e6;
		if(nMassResidual!=pMassResidual)	nVelocityLoops = nVelocityLoops + 1;
		pMassResidual=maxMassResidual*10e6;
	}

	printf("Intermediate Mass residual: %e\n",maxMassResidual);

	for(nonLinear=0;nonLinear<nPressureLoops;nonLinear++)
	{

	pCorrEqn();		//finding pressure correction based on initial velocities
	updateVelocity();	//updating velocties
	massResidual();		//finding final max mass residual

	}

    printf("Final Mass residual: %e\n",maxMassResidual);
	if(maxMassResidual>10e-6||maxMassResidual<-10e-6)		nSimplerLoops = nSimplerLoops + 1;;

 }
    
	printf("No of pressure loops: %d\n",nPressureLoops);
	printf("No of Velocity loops: %d\n",nVelocityLoops);
	printf("No of SIMPLER loops: %d\n",nSimplerLoops);
}

//Calculates the coefficients of the pressure correction equation using hybrid scheme
void pEqnCoeff()
{
    uEqnCoeff();    //updates the value of Ae
	vEqnCoeff();	//updates the value of An
	wEqnCoeff();	//updates the value of At

    for  (k=2; k<=nz+1; k++)

	 for (j=2; j<=ny+1; j++)
   {
	  for(i=2; i<=nx+1; i++)
	{

	aE[i][j][k] = rho * Ae[i][j][k]   * dy*dz; //rho at e
	aW[i][j][k] = rho * Ae[i-1][j][k] * dy*dz; //rho at w
	aN[i][j][k] = rho * An[i][j][k]   * dx*dz; //rho at n
	aS[i][j][k] = rho * An[i][j-1][k] * dx*dz; //rho at s
	aT[i][j][k] = rho * At[i][j][k]   * dx*dy; //rho at t
	aB[i][j][k] = rho * At[i][j][k-1] * dx*dy; //rho at b
	aP[i][j][k] = aE[i][j][k] + aW[i][j][k] + aN[i][j][k] + aS[i][j][k] + aT[i][j][k] + aB[i][j][k];

	if(i==nx+1) aE[i][j][k]=0;
	if(i==2)    aW[i][j][k]=0;
	if(j==ny+1) aN[i][j][k]=0;
	if(j==2)    aS[i][j][k]=0;
	if(k==nz+1) aT[i][j][k]=0;
	if(k==2)    aB[i][j][k]=0;


	//residual mass
	Sm[i][j][k] = (rho*pseudoU[i][j][k][l]-rho*pseudoU[i-1][j][k][l])*dy*dz + (rho*pseudoV[i][j][k][l]-rho*pseudoV[i][j-1][k][l])*dx*dz + (rho*pseudoW[i][j][k][l]-rho*pseudoW[i][j][k-1][l])*dx*dy;  // rho at respective places
	Sc[i][j][k]=-Sm[i][j][k];

	}}

}

//calculates the actual value of pressure
void pEqn()
{
    pEqnCoeff();  //updating pressure correction coefficients

	for  (k=2; k<=nz+1; k++)

	 for (j=2; j<=ny+1; j++)
   {
	  for(i=2; i<=nx+1; i++)
	{

		lower[i]=-aW[i][j][k]/aP[i][j][k];  //lower diagonal of A Matrix

		upper[i]=-aE[i][j][k]/aP[i][j][k];  //upper diagonal of A Matrix

		//Matrix B
		if(i==2)
         c[i] = aN[i][j][k]/aP[i][j][k]*p[i][j+1][k][l] + aS[i][j][k]/aP[i][j][k]*p[i][j-1][k][l] + aT[i][j][k]/aP[i][j][k]*p[i][j][k+1][l] + aB[i][j][k]/aP[i][j][k]*p[i][j][k-1][l] + Sc[i][j][k]/aP[i][j][k] + aW[i][j][k]/aP[i][j][k]*p[i-1][j][k][l]; //check this

		else if(i==nx+1)
		 c[i] = aN[i][j][k]/aP[i][j][k]*p[i][j+1][k][l] + aS[i][j][k]/aP[i][j][k]*p[i][j-1][k][l] + aT[i][j][k]/aP[i][j][k]*p[i][j][k+1][l] + aB[i][j][k]/aP[i][j][k]*p[i][j][k-1][l] + Sc[i][j][k]/aP[i][j][k] + aE[i][j][k]/aP[i][j][k]*p[i+1][j][k][l];

		else
		 c[i] = aN[i][j][k]/aP[i][j][k]*p[i][j+1][k][l] + aS[i][j][k]/aP[i][j][k]*p[i][j-1][k][l] + aT[i][j][k]/aP[i][j][k]*p[i][j][k+1][l] + aB[i][j][k]/aP[i][j][k]*p[i][j][k-1][l] + Sc[i][j][k]/aP[i][j][k];
	}

	n=nx;  //size of Matrix A and B

	Thomas();   //Thomas algorithm solution


	//Equating solved value back to pressure correction
	for(thomasI=2;thomasI<=n+1;thomasI++)
	  p[thomasI][j][k][l]=x[thomasI];


   }
	boundaryConditions();		//updates boundary conditions
	ghostCells();  			//updates the ghost cells

}

//Calculates the values of pseudo velocity 
void pseudoVelocity()
{
    uEqnCoeff();
    vEqnCoeff();
	wEqnCoeff();


	for  (k=2; k<=nz+1; k++)
	 for (i=2; i<=nx;   i++)
	  for(j=2; j<=ny+1; j++)
        pseudoU[i][j][k][l] = (aee[i][j][k]*u[i+1][j][k][l] + aw[i][j][k]*u[i-1][j][k][l] + aTe[i][j][k]*u[i][j][k+1][l] + aBe[i][j][k]*u[i][j][k-1][l] + aNe[i][j][k]*u[i][j+1][k][l] + aSe[i][j][k]*u[i][j-1][k][l])/ae[i][j][k] + B[i][j][k];


	for  (i=2; i<=nx+1; i++)
	 for (j=2; j<=ny;   j++)
	  for(k=2; k<=nz+1; k++)
        pseudoV[i][j][k][l] = (anE[i][j][k]*v[i+1][j][k][l] + anW[i][j][k]*v[i-1][j][k][l] + ann[i][j][k]*v[i][j+1][k][l] + as[i][j][k]*v[i][j-1][k][l] + anB[i][j][k]*v[i][j][k-1][l] + anT[i][j][k]*v[i][j][k+1][l])/an[i][j][k] + BB[i][j][k];


	for  (j=2; j<=ny+1; j++)
	 for (k=2; k<=nz;   k++)
	  for(i=2; i<=nx+1; i++)
        pseudoW[i][j][k][l] = (att[i][j][k]*w[i][j][k+1][l] + ab[i][j][k]*w[i][j][k-1][l] + aNt[i][j][k]*w[i][j+1][k][l] + aSt[i][j][k]*w[i][j-1][k][l] + atW[i][j][k]*w[i-1][j][k][l] + atE[i][j][k]*w[i+1][j][k][l])/at[i][j][k] + BBB[i][j][k];


}






