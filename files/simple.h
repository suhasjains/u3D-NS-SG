//************************************************************
//SOURCE CODE FUNCTIONS
//************************************************************

//calculates the value of initial u velocity
void uEqn()
{

    uEqnCoeff();  //updating coefficients

	for  (k=2; k<=nz+1; k++)
 {
	 for (i=2; i<=nx;   i++)
   {
	  for(j=2; j<=ny+1; j++)
	{

		lower[j]=-aSe[i][j][k]/ae[i][j][k];  //lower diagonal of A Matrix

		upper[j]=-aNe[i][j][k]/ae[i][j][k];  //upper diagonal of A Matrix

		//Matrix B
		if(j==2)
         c[j]=aee[i][j][k]/ae[i][j][k]*u[i+1][j][k][l] + aw[i][j][k]/ae[i][j][k]*u[i-1][j][k][l] + aTe[i][j][k]/ae[i][j][k]*u[i][j][k+1][l] + aBe[i][j][k]/ae[i][j][k]*u[i][j][k-1][l] + (p[i][j][k][l] - p[i+1][j][k][l])*Ae[i][j][k] + B[i][j][k] + uDC[i][j][k]/ae[i][j][k] + aSe[i][j][k]/ae[i][j][k]*u[i][j-1][k][l];

		else if(j==ny+1)
		 c[j]=aee[i][j][k]/ae[i][j][k]*u[i+1][j][k][l] + aw[i][j][k]/ae[i][j][k]*u[i-1][j][k][l] + aTe[i][j][k]/ae[i][j][k]*u[i][j][k+1][l] + aBe[i][j][k]/ae[i][j][k]*u[i][j][k-1][l] + (p[i][j][k][l] - p[i+1][j][k][l])*Ae[i][j][k] + B[i][j][k] + uDC[i][j][k]/ae[i][j][k] + aNe[i][j][k]/ae[i][j][k]*u[i][j+1][k][l];

		else
		 c[j]=aee[i][j][k]/ae[i][j][k]*u[i+1][j][k][l] + aw[i][j][k]/ae[i][j][k]*u[i-1][j][k][l] + aTe[i][j][k]/ae[i][j][k]*u[i][j][k+1][l] + aBe[i][j][k]/ae[i][j][k]*u[i][j][k-1][l] + (p[i][j][k][l] - p[i+1][j][k][l])*Ae[i][j][k] + B[i][j][k] + uDC[i][j][k]/ae[i][j][k];	
	}

	n=ny;  //size of Matrix A and B

	Thomas();   //Thomas algorithm solution


	//Equating solved value back to velocity
	for(thomasI=2;thomasI<=n+1;thomasI++)
	  u[i][thomasI][k][l]=x[thomasI]*vUnderRelaxCoeff+u[i][thomasI][k][l]*(1-vUnderRelaxCoeff);
   }
 }

	boundaryConditions();
	ghostCells();  //updates the ghost cells

}

//Updates u Equation Coefficients using hybrid scheme
void uEqnCoeff()
{

	for  (k=2; k<=nz+1; k++)
 {
	 for (i=2; i<=nx;   i++)
   {
	  for(j=2; j<=ny+1; j++)
	{


	//Convective Fluxes
	 	Fte = rho * ((1 - weight)*std::max(w[i][j][k][l],w[i+1][j][k][l])     + weight*std::min(w[i][j][k][l],w[i+1][j][k][l])    ) *dx*dy ;  //weighted mean of wt and wtE

		Fbe = rho * ((1 - weight)*std::max(w[i][j][k-1][l],w[i+1][j][k-1][l]) + weight*std::min(w[i][j][k-1][l],w[i+1][j][k-1][l])) *dx*dy ;  //weighted mean of wb and wbE

		FE  = rho * ((1 - weight)*std::max(u[i][j][k][l],u[i+1][j][k][l])     + weight*std::min(u[i][j][k][l],u[i+1][j][k][l])    ) *dy*dz ;  //weighted mean of ue and uee

		FP  = rho * ((1 - weight)*std::max(u[i][j][k][l],u[i-1][j][k][l])     + weight*std::min(u[i][j][k][l],u[i-1][j][k][l])    ) *dy*dz ;  //weighted mean of ue and uw

		Fne = rho * ((1 - weight)*std::max(v[i][j][k][l],v[i+1][j][k][l])     + weight*std::min(v[i][j][k][l],v[i+1][j][k][l])    ) *dx*dz ;  //weighted mean of vn and vnE

		Fse = rho * ((1 - weight)*std::max(v[i][j-1][k][l],v[i+1][j-1][k][l]) + weight*std::min(v[i][j-1][k][l],v[i+1][j-1][k][l])) *dx*dz ;  //weighted mean of vs and vsE



	//constant
 	b[i][j][k] = mu*(v[i+1][j][k][l]-v[i][j][k][l])*dz - mu*(v[i+1][j-1][k][l]-v[i][j-1][k][l])*dz + mu*(w[i+1][j][k][l]-w[i][j][k][l])*dy - mu*(w[i+1][j][k-1][l]-w[i][j][k-1][l])*dy + rho*xGrav*dx*dy*dz + rho*u[i][j][k][l-1]*dx*dy*dz/dt; //rho at e

	//Diffusion Fluxes
	DE  = 2 * mu * dy * dz / dx; //mu at E
	DP  = 2 * mu * dy * dz / dx; //mu at P
	Dne =     mu * dx * dz / dy; //mu at ne
	Dse =     mu * dx * dz / dy; //mu at se
	Dte =     mu * dx * dy / dz; //mu at te
	Dbe =     mu * dx * dy / dz; //mu at be

	if((i==2)||(i==(nx))||(j==2)||(j==(ny))) 			DE = DP = 0;
    
    powerLawU();
	
	ae[i][j][k]  = aee[i][j][k] + aw[i][j][k] + aNe[i][j][k] + aSe[i][j][k] + aTe[i][j][k] + aBe[i][j][k] + rho*dx*dy*dz/dt; //rho at e

	quickU();
	
	Ae[i][j][k] = dy*dz/ae[i][j][k];

	B[i][j][k]  = b[i][j][k]/ ae[i][j][k];
	
	}}}

}

//calculates the value of initial v velocity
void vEqn()
{

    vEqnCoeff();  //updating coefficients

	for  (i=2; i<=nx+1; i++)
 {
	 for (j=2; j<=ny;   j++)
   {
	  for(k=2; k<=nz+1; k++)
	{

		lower[k]=-anB[i][j][k]/an[i][j][k];  //lower diagonal of A Matrix

		upper[k]=-anT[i][j][k]/an[i][j][k];  //upper diagonal of A Matrix

		//Matrix B
		if(k==2)
         c[k]=anE[i][j][k]/an[i][j][k]*v[i+1][j][k][l] + anW[i][j][k]/an[i][j][k]*v[i-1][j][k][l] + ann[i][j][k]/an[i][j][k]*v[i][j+1][k][l] + as[i][j][k]/an[i][j][k]*v[i][j-1][k][l] + (p[i][j][k][l] - p[i][j+1][k][l])*An[i][j][k] + BB[i][j][k] + vDC[i][j][k]/an[i][j][k] + anB[i][j][k]/an[i][j][k]*v[i][j][k-1][l];

		else if(k==nz+1)
		 c[k]=anE[i][j][k]/an[i][j][k]*v[i+1][j][k][l] + anW[i][j][k]/an[i][j][k]*v[i-1][j][k][l] + ann[i][j][k]/an[i][j][k]*v[i][j+1][k][l] + as[i][j][k]/an[i][j][k]*v[i][j-1][k][l] + (p[i][j][k][l] - p[i][j+1][k][l])*An[i][j][k] + BB[i][j][k] + vDC[i][j][k]/an[i][j][k] + anT[i][j][k]/an[i][j][k]*v[i][j][k+1][l];

		else
		 c[k]=anE[i][j][k]/an[i][j][k]*v[i+1][j][k][l] + anW[i][j][k]/an[i][j][k]*v[i-1][j][k][l] + ann[i][j][k]/an[i][j][k]*v[i][j+1][k][l] + as[i][j][k]/an[i][j][k]*v[i][j-1][k][l] + (p[i][j][k][l] - p[i][j+1][k][l])*An[i][j][k] + BB[i][j][k] + vDC[i][j][k]/an[i][j][k];  //check this
	}

	n=nz;  //size of Matrix A and B

	Thomas();   //Thomas algorithm solution

	//Equating solved value back to velocity
	for(thomasI=2;thomasI<=n+1;thomasI++)
	  v[i][j][thomasI][l]=x[thomasI]*vUnderRelaxCoeff+v[i][j][thomasI][l]*(1-vUnderRelaxCoeff);
  }
 }

	boundaryConditions();
	ghostCells();  //updates the ghost cells

}

//Updates v Equation Coefficients using hybrid scheme
void vEqnCoeff()
{

	for  (i=2; i<=nx+1; i++)
 {
	 for (j=2; j<=ny;   j++)
   {
	  for(k=2; k<=nz+1; k++)
	{


	//Convective Fluxes
        Fned = rho * ((1 - weight)*std::max(u[i][j+1][k][l],u[i][j][k][l])     + weight*std::min(u[i][j+1][k][l],u[i][j][k][l])     )*dy*dz; //weighted mean of ue and uNe

		Fnw  = rho * ((1 - weight)*std::max(u[i-1][j+1][k][l],u[i-1][j][k][l]) + weight*std::min(u[i-1][j+1][k][l],u[i-1][j][k][l]) )*dy*dz; //weighted mean of uw and uNw

		FN   = rho * ((1 - weight)*std::max(v[i][j+1][k][l],v[i][j][k][l])     + weight*std::min(v[i][j+1][k][l],v[i][j][k][l])     )*dx*dz; //weighted mean of vn and vnn

		FPd  = rho * ((1 - weight)*std::max(v[i][j-1][k][l],v[i][j][k][l])     + weight*std::min(v[i][j-1][k][l],v[i][j][k][l])     )*dx*dz; //weighted mean of vn and vs

		Fnt  = rho * ((1 - weight)*std::max(w[i][j+1][k][l],w[i][j][k][l])     + weight*std::min(w[i][j+1][k][l],w[i][j][k][l])     )*dx*dy; //weighted mean of wt and wNt

		Fnb  = rho * ((1 - weight)*std::max(w[i][j+1][k-1][l],w[i][j][k-1][l]) + weight*std::min(w[i][j+1][k-1][l],w[i][j][k-1][l]) )*dx*dy; //weighted mean of wb and wNb

	//constant
	bb[i][j][k] = mu*(u[i][j+1][k][l]-u[i][j][k][l])*dz - mu*(u[i-1][j+1][k][l]-u[i-1][j][k][l])*dz + mu*(w[i][j+1][k][l]-w[i][j][k][l])*dx - mu*(w[i][j+1][k-1][l]-w[i][j][k-1][l])*dx + rho*yGrav*dx*dy*dz + rho*v[i][j][k][l-1]*dx*dy*dz/dt;  //rho at n

	//Diffusion Fluxes
	Dned=     mu * dy * dz/dx; //mu at ne
	Dnw =     mu * dy * dz/dx; //mu at nw
	DN  = 2 * mu * dx * dz/dy; //mu at N
	DPd = 2 * mu * dx * dz/dy; //mu at P
	Dnt =     mu * dx * dy/dz; //mu at nt
	Dnb =     mu * dx * dy/dz; //mu at nb

	if((i==2)||(i==(nx))||(j==2)||(j==(ny))) 			DN = DPd = 0;
	
	powerLawV();
	
	an[i][j][k]  = anE[i][j][k] + anW[i][j][k] + ann[i][j][k] + as[i][j][k] + anT[i][j][k] + anB[i][j][k] + rho*dx*dy*dz/dt; //rho at n

	quickV();
	
	An[i][j][k] = dx*dz/an[i][j][k];

	BB[i][j][k] = bb[i][j][k]/an[i][j][k];

	}}}

}

//calculates the value of initial w velocity
void wEqn()
{
    wEqnCoeff();  //updating coefficients


	for  (j=2; j<=ny+1; j++)
 {
	 for (k=2; k<=nz;   k++)
   {
	  for(i=2; i<=nx+1; i++)
	{

		lower[i]=-atW[i][j][k]/at[i][j][k];  //lower diagonal of A Matrix

		upper[i]=-atE[i][j][k]/at[i][j][k];  //upper diagonal of A Matrix

		//Matrix B
		if(i==2)
         c[i]=att[i][j][k]/at[i][j][k]*w[i][j][k+1][l] + ab[i][j][k]/at[i][j][k]*w[i][j][k-1][l] + aNt[i][j][k]/at[i][j][k]*w[i][j+1][k][l] + aSt[i][j][k]/at[i][j][k]*w[i][j-1][k][l] + (p[i][j][k][l] - p[i][j][k+1][l])*At[i][j][k] + BBB[i][j][k] + wDC[i][j][k]/at[i][j][k] + atW[i][j][k]/at[i][j][k]*w[i-1][j][k][l];

		else if(i==nx+1)
		 c[i]=att[i][j][k]/at[i][j][k]*w[i][j][k+1][l] + ab[i][j][k]/at[i][j][k]*w[i][j][k-1][l] + aNt[i][j][k]/at[i][j][k]*w[i][j+1][k][l] + aSt[i][j][k]/at[i][j][k]*w[i][j-1][k][l] + (p[i][j][k][l] - p[i][j][k+1][l])*At[i][j][k] + BBB[i][j][k] + wDC[i][j][k]/at[i][j][k] + atE[i][j][k]/at[i][j][k]*w[i+1][j][k][l];

		else
		 c[i]=att[i][j][k]/at[i][j][k]*w[i][j][k+1][l] + ab[i][j][k]/at[i][j][k]*w[i][j][k-1][l] + aNt[i][j][k]/at[i][j][k]*w[i][j+1][k][l] + aSt[i][j][k]/at[i][j][k]*w[i][j-1][k][l] + (p[i][j][k][l] - p[i][j][k+1][l])*At[i][j][k] + BBB[i][j][k] + wDC[i][j][k]/at[i][j][k];  
	}

	n=nx;  //size of Matrix A and B

	Thomas();   //Thomas algorithm solution

	//Equating solved value back to velocity
	for(thomasI=2;thomasI<=n+1;thomasI++)
	  w[thomasI][j][k][l]=x[thomasI]*vUnderRelaxCoeff+w[thomasI][j][k][l]*(1-vUnderRelaxCoeff);


   }
 }

	boundaryConditions();
	ghostCells();  //updates the ghost cells

}


//Updates w Equation Coefficients using hybrid scheme
void wEqnCoeff()
{

   for  (j=2; j<=ny+1; j++)
   {
	 for (k=2; k<=nz;   k++)
    {
	  for(i=2; i<=nx+1; i++)
	 {
	 //Convective Fluxes
		Fted = rho * ((1 - weight)*std::max(u[i][j][k+1][l],u[i][j][k][l])     + weight*std::min(u[i][j][k+1][l],u[i][j][k][l])       )*dy*dz; //weighted mean of ue and uTe

		Ftw  = rho * ((1 - weight)*std::max(u[i-1][j][k+1][l],u[i-1][j][k][l]) + weight*std::min(u[i-1][j][k+1][l],u[i-1][j][k][l])   )*dy*dz; //weighted mean of uw and uTw

		Fntd = rho * ((1 - weight)*std::max(v[i][j][k+1][l],v[i][j][k][l])     + weight*std::min(v[i][j][k+1][l],v[i][j][k][l])       )*dx*dz; //weighted mean of vn and vnT

		Fst  = rho * ((1 - weight)*std::max(v[i][j-1][k+1][l],v[i][j-1][k][l]) + weight*std::min(v[i][j-1][k+1][l],v[i][j-1][k][l])   )*dx*dz; //weighted mean of vs and vsT

		FT   = rho * ((1 - weight)*std::max(w[i][j][k+1][l],w[i][j][k][l])     + weight*std::min(w[i][j][k+1][l],w[i][j][k][l])       )*dx*dy; //weighted mean of wt and wtt

		FPdd = rho * ((1 - weight)*std::max(w[i][j][k-1][l],w[i][j][k][l])     + weight*std::min(w[i][j][k-1][l],w[i][j][k][l])	      ) *dx*dy; //weighted mean of wt and wb


	 //constant
	 bbb[i][j][k] = mu*(u[i][j][k+1][l]-u[i][j][k][l])*dy - mu*(u[i-1][j][k+1][l]-u[i-1][j][k][l])*dy + mu*(v[i][j][k+1][l]-v[i][j][k][l])*dx - mu*(v[i][j-1][k+1][l]-v[i][j-1][k][l])*dx + rho*zGrav*dx*dy*dz + rho*w[i][j][k][l-1]*dx*dy*dz/dt;  //rho at t


	 //Diffusion Fluxes
	 DT   = 2 * mu * dx*dy/dz; //mu at T
	 DPdd = 2 * mu * dx*dy/dz; //mu at P
	 Dted =     mu * dy*dz/dx; //mu at te
	 Dtw  =     mu * dy*dz/dx; //mu at tw
	 Dntd =     mu * dx*dz/dy; //mu at nt
	 Dst  =     mu * dx*dz/dy; //mu at st


	 if((i==2)||(i==(nx))||(j==2)||(j==(ny))) 			DT = DPdd = 0;
	 
	 
	 powerLawW();
	 
	 at[i][j][k]  = att[i][j][k] + ab[i][j][k] + atE[i][j][k] + atW[i][j][k] + aNt[i][j][k] + aSt[i][j][k] + rho*dx*dy*dz/dt; //rho at t
           
	 quickW();
	 
	 At[i][j][k]  = dx*dy/at[i][j][k];

	 BBB[i][j][k] = bbb[i][j][k]/ at[i][j][k];

	 }}}
}

//Interpolation of the velocity vectors back from staggerd cell centers to the pressure cell centers
void makeVelocityVectors()
{

	for  (k=2; k<=nz+1; k++)
	 for (i=2; i<=nx+1; i++)
	  for(j=2; j<=ny+1; j++)
	{
	 finalU[i][j][k][l] = (u[i][j][k][l] + u[i-1][j][k][l])/(2);
	 finalV[i][j][k][l] = (v[i][j][k][l] + v[i][j-1][k][l])/(2);
	 finalW[i][j][k][l] = (w[i][j][k][l] + w[i][j][k-1][l])/(2);
	}
}

//Updates the values of the velocity
void momentumResidual()
{

	uMaxRes=0;
	vMaxRes=0;
	wMaxRes=0;

	for  (k=2; k<=nz+1; k++)
	 for (i=2; i<=nx;   i++)
	  for(j=2; j<=ny+1; j++)
	{
	  uRes[i][j][k]=(ae[i][j][k]*u[i][j][k][l]-aee[i][j][k]*u[i+1][j][k][l] - aw[i][j][k]*u[i-1][j][k][l] - aTe[i][j][k]*u[i][j][k+1][l] - aBe[i][j][k]*u[i][j][k-1][l] - aNe[i][j][k]*u[i][j+1][k][l] - aSe[i][j][k]*u[i][j-1][k][l] - B[i][j][k]*ae[i][j][k] - uDC[i][j][k] - Ae[i][j][k] * (p[i][j][k][l] - p[i+1][j][k][l]) );

	  if((uRes[i][j][k]>uMaxRes)||(uRes[i][j][k]<-uMaxRes))	uMaxRes=uRes[i][j][k];
	}


	for  (i=2; i<=nx+1; i++)
	 for (j=2; j<=ny;   j++)
	  for(k=2; k<=nz+1; k++)
	{
	  vRes[i][j][k]=(an[i][j][k]*v[i][j][k][l] - anE[i][j][k]*v[i+1][j][k][l] - anW[i][j][k]*v[i-1][j][k][l] - ann[i][j][k]*v[i][j+1][k][l] - as[i][j][k]*v[i][j-1][k][l] - anB[i][j][k]*v[i][j][k-1][l] - anT[i][j][k]*v[i][j][k+1][l] - BB[i][j][k]*an[i][j][k] - vDC[i][j][k] - An[i][j][k] * (p[i][j][k][l] - p[i][j+1][k][l]));

	  if((vRes[i][j][k]>vMaxRes)||(vRes[i][j][k]<-vMaxRes))	vMaxRes=vRes[i][j][k];
	}


	for  (j=2; j<=ny+1; j++)
	 for (k=2; k<=nz;   k++)
	  for(i=2; i<=nx+1; i++)
	{
	 wRes[i][j][k]=(at[i][j][k]*w[i][j][k][l] - att[i][j][k]*w[i][j][k+1][l] - ab[i][j][k]*w[i][j][k-1][l] - aNt[i][j][k]*w[i][j+1][k][l] - aSt[i][j][k]*w[i][j-1][k][l] - atW[i][j][k]*w[i-1][j][k][l] - atE[i][j][k]*w[i+1][j][k][l] - BBB[i][j][k]*at[i][j][k] - wDC[i][j][k] - At[i][j][k] * (p[i][j][k][l] - p[i][j][k+1][l]));

	 if((wRes[i][j][k]>wMaxRes)||(wRes[i][j][k]<-wMaxRes))	wMaxRes=wRes[i][j][k];
	}

}

//Updates the values of the velocity
void updateVelocity()
{
    //updates u velocity
	for  (k=2; k<=nz+1; k++)
	 for (i=2; i<=nx;   i++)
	  for(j=2; j<=ny+1; j++)
        u[i][j][k][l] = u[i][j][k][l] + Ae[i][j][k] * (pCorr[i][j][k][l] - pCorr[i+1][j][k][l]);



	//updates v velocity
	for  (i=2; i<=nx+1; i++)
	 for (j=2; j<=ny;   j++)
	  for(k=2; k<=nz+1; k++)
        v[i][j][k][l] = v[i][j][k][l] + An[i][j][k] * (pCorr[i][j][k][l] - pCorr[i][j+1][k][l]);



	//updates w velocity
	for  (j=2; j<=ny+1; j++)
	 for (k=2; k<=nz;   k++)
	  for(i=2; i<=nx+1; i++)
        w[i][j][k][l] = w[i][j][k][l] + At[i][j][k] * (pCorr[i][j][k][l] - pCorr[i][j][k+1][l]);


}

//calculates the value of pressure correction
void pCorrEqn()
{
    pCorrEqnCoeff();  //updating pressure correction coefficients

	for  (k=2; k<=nz+1; k++)
  {
	 for (j=2; j<=ny+1; j++)
   {
	  for(i=2; i<=nx+1; i++)
	{
		lower[i]=-aW[i][j][k]/aP[i][j][k];  //lower diagonal of A Matrix

		upper[i]=-aE[i][j][k]/aP[i][j][k];  //upper diagonal of A Matrix

		//Matrix B
		if(i==2)
         c[i] = aN[i][j][k]/aP[i][j][k]*pCorr[i][j+1][k][l] + aS[i][j][k]/aP[i][j][k]*pCorr[i][j-1][k][l] + aT[i][j][k]/aP[i][j][k]*pCorr[i][j][k+1][l] + aB[i][j][k]/aP[i][j][k]*pCorr[i][j][k-1][l] + Sc[i][j][k]/aP[i][j][k] + aW[i][j][k]/aP[i][j][k]*pCorr[i-1][j][k][l]; //check this

		else if(i==nx+1)
		 c[i] = aN[i][j][k]/aP[i][j][k]*pCorr[i][j+1][k][l] + aS[i][j][k]/aP[i][j][k]*pCorr[i][j-1][k][l] + aT[i][j][k]/aP[i][j][k]*pCorr[i][j][k+1][l] + aB[i][j][k]/aP[i][j][k]*pCorr[i][j][k-1][l] + Sc[i][j][k]/aP[i][j][k] + aE[i][j][k]/aP[i][j][k]*pCorr[i+1][j][k][l];

		else
		 c[i] = aN[i][j][k]/aP[i][j][k]*pCorr[i][j+1][k][l] + aS[i][j][k]/aP[i][j][k]*pCorr[i][j-1][k][l] + aT[i][j][k]/aP[i][j][k]*pCorr[i][j][k+1][l] + aB[i][j][k]/aP[i][j][k]*pCorr[i][j][k-1][l] + Sc[i][j][k]/aP[i][j][k];
	}

	n=nx;  //size of Matrix A and B

	Thomas();   //Thomas algorithm solution


	//Equating solved value back to pressure correction
	for(thomasI=2;thomasI<=n+1;thomasI++)
	  pCorr[thomasI][j][k][l]=x[thomasI];


   }}
	boundaryConditions();
	ghostCells();  //updates the ghost cells

}

//Calculates the coefficients of the pressure correction equation using hybrid scheme
void pCorrEqnCoeff()
{
    uEqnCoeff();  //updates the value of Ae
	vEqnCoeff();	//updates the value of An
	wEqnCoeff();	//updates the value of At

    for  (k=2; k<=nz+1; k++)
  {
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
	Sm[i][j][k] = (rho*u[i][j][k][l]-rho*u[i-1][j][k][l])*dy*dz + (rho*v[i][j][k][l]-rho*v[i][j-1][k][l])*dx*dz + (rho*w[i][j][k][l]-rho*w[i][j][k-1][l])*dx*dy;  // rho at respective places
	Sc[i][j][k]=-Sm[i][j][k];

	}}}

}

//Updates the values of the pressure
void updatePressure()
{

	for  (k=2; k<=nz+1; k++)
	 for (j=2; j<=ny+1; j++)
	  for(i=2; i<=nx+1; i++)
	p[i][j][k][l] = p[i][j][k][l] + pUnderRelaxCoeff * pCorr[i][j][k][l];

}

//SIMPLE algorithm loop
void simple()
{

	printf("Starting SIMPLE loop\n");
	printf("Initial mass Residual: %e\n",maxMassResidual);

	for(simplerVar=0;simplerVar<nSimpleLoops;simplerVar++)
 {

	nLoops();

	for(nonLinear=0;nonLinear<nVelocityLoops;nonLinear++)
	{

		uEqn();
		vEqn();			//finding initial velocities
		wEqn();


		massResidual();		//finding intermediate max mass residual



		nMassResidual=maxMassResidual*10e5;
		if((nMassResidual!=pMassResidual))	nVelocityLoops = nVelocityLoops + 1;
		pMassResidual=maxMassResidual*10e5;


	}

	printf("Intermediate Mass residual: %e\n",maxMassResidual);


	for(nonLinear=0;nonLinear<nPressureLoops;nonLinear++)
	{
	pCorrEqn();		//finding pressure correction based on initial velocities
	updatePressure();       //updating pressure
	updateVelocity();	//updating velocties
      
	momentumResidual();

	massResidual();		//finding final max mass residual
	}

	printf("Final Mass residual: %e\n",maxMassResidual);

	if(maxMassResidual>10e-6||maxMassResidual<-10e-6)		nSimpleLoops = nSimpleLoops + 1;


 }

	printf("No of SIMPLE loops: %d\n",nSimpleLoops);
	printf("No of pressure loops: %d\n",nPressureLoops);
	printf("No of Velocity loops: %d\n",nVelocityLoops);
}

//calculates the max courant number
void courantNumber()
{
	uMax=0;
	vMax=0;
	wMax=0;

	for  (k=2; k<=nz+1; k++)
	 for (j=2; j<=ny+1; j++)
	  for(i=2; i<=nx+1; i++)
	{
		if(finalU[i][j][k][l]>uMax)	uMax=finalU[i][j][k][l];
		if(finalV[i][j][k][l]>vMax)	vMax=finalV[i][j][k][l];			//finding out max velocity
		if(finalW[i][j][k][l]>wMax)	wMax=finalW[i][j][k][l];
	}

	maxCourantNumber = (uMax/dx + vMax/dy + wMax/dz)*dt;

	printf("Max courant number: %lf\n",maxCourantNumber);


}

//Increments the time variables
void incrementTime()
{

	if(maxCourantNumber>allowableCourantNumber)
	{
		dt = allowableCourantNumber/(4*(uMax/dx + vMax/dy + wMax/dz));
		courantNumberControl=0;
	}

	if(courantNumberControl==10)						dt = deltaT;


	simulationTime = simulationTime + dt;	//time increment

	t=t+1;			//output file name parameter
	courantNumberControl = courantNumberControl+1;

}

//Calculates max mass residual
void massResidual()
{

	maxMassResidual = 0;

	for  (k=2; k<=nz+1; k++)
	 for (j=2; j<=ny+1; j++)
	  for(i=2; i<=nx+1; i++)
	{
       	  Sm[i][j][k] = (rho*u[i][j][k][l]-rho*u[i-1][j][k][l])*dy*dz + (rho*v[i][j][k][l]-rho*v[i][j-1][k][l])*dx*dz + (rho*w[i][j][k][l]-rho*w[i][j][k-1][l])*dx*dy;
	  if((Sm[i][j][k]>maxMassResidual)||(Sm[i][j][k]<-maxMassResidual))		maxMassResidual = Sm[i][j][k];
	}


}

//shifts the values to l=0
void shift()
{

	for  (k=0; k<(zNumberOfCells+5); k++)
	 for (j=0; j<(yNumberOfCells+5); j++)
	  for(i=0; i<(xNumberOfCells+5); i++)
	{
	  u[i][j][k][0]=u[i][j][k][1];
	  v[i][j][k][0]=v[i][j][k][1];
	  w[i][j][k][0]=w[i][j][k][1];
	  p[i][j][k][0]=p[i][j][k][1];			//copying values to l=0
	  pCorr[i][j][k][0]=pCorr[i][j][k][1];


	  pCorr[i][j][k][1]=0;	                          //equating values in l=1 to zero

	}

}

//Starts CPU clock
void startCPUClock()
{

	t=1;			//output file starts from 1st time step

	simulationTime=dt;		//time counting starts from dt

	assert((start = clock())!=-1);			//starting real time

}

//Calculation of the spacial and temporal step sizes
void stepSize()
{

	dx = xDomainLength/xNumberOfCells;
	dy = yDomainLength/yNumberOfCells;
	dz = zDomainLength/zNumberOfCells;
	deltaT=(endTime-startTime)/nTimeSteps;  //initial time step

	dt=deltaT;  				//equating initial time step to time step

	nx=xNumberOfCells;
	ny=yNumberOfCells;
	nz=zNumberOfCells;

}

//Solving AX=B using Thomas algorithm
void Thomas()
{

	m[2]=1;  //changed principal diagonal of A matrix
	cc[2]=c[2];  //changed B matrix

	for(thomasI=3;thomasI<=n+1;thomasI++)
   {

  		m[thomasI]  = 1 - (lower[thomasI] * upper[thomasI-1] / m[thomasI-1]);
  		cc[thomasI] = c[thomasI] - (cc[thomasI-1] * lower[thomasI] / m[thomasI-1]);
   }


	x[n+1] = cc[n+1] / m[n+1];


	for(thomasI=n;thomasI>=2;thomasI--)
	{
  	 x[thomasI] = (cc[thomasI] - upper[thomasI] * x[thomasI+1]) / m[thomasI];
	}

}



