//*******************************************************
//Semi upwind schemes
//******************************************************

//Coefficients of U Equation for semi upwind scheme
void semiUpwindU()
{
	aee[i][j][k] = std::max(-FE,0.0) + DE;  //upwind Scheme
	aw[i][j][k]  = std::max(FP,0.0)  + DP;  //upwind Scheme
	aNe[i][j][k] = Dne - Fne/2;				//central difference
	aSe[i][j][k] = Dse + Fse/2;				//central difference
	aTe[i][j][k] = Dte - Fte/2;				//central difference
	aBe[i][j][k] = Dbe + Fbe/2;				//central difference
	
	//Deferred Correction
	uUDS[i][j][k] = u[i][j][k][l] * std::max(FE,0.0) - u[i+1][j][k][l] * std::max(-FE,0.0) - u[i-1][j][k][l] * std::max(FP,0.0) + u[i][j][k][l] * std::max(-FP,0.0);

	uCDS[i][j][k] = FE*(u[i][j][k][l] + u[i+1][j][k][l])/2 - FP*(u[i][j][k][l] + u[i-1][j][k][l])/2;

	uDC[i][j][k] = DCvalue*( uUDS[i][j][k] - uCDS[i][j][k]);
}

//Coefficients for V Equation for semi upwind scheme
void semiUpwindV()
{
	anE[i][j][k] = Dned- Fned/ 2;				//central differencing
	anW[i][j][k] = Dnw + Fnw / 2;				//central differencing
	ann[i][j][k] = DN  + std::max(-FN,0.0);  	//upwind Scheme
	as[i][j][k]  = DPd + std::max(FPd,0.0);  	//upwind Scheme
	anT[i][j][k] = Dnt - Fnt / 2;				//central differencing
	anB[i][j][k] = Dnb + Fnb / 2;				//central differencing
	
	//Deferred Correction
	vUDS[i][j][k] = v[i][j][k][l] * std::max(FN,0.0) - v[i][j+1][k][l] * std::max(-FN,0.0) - v[i][j-1][k][l] * std::max(FPd,0.0) + v[i][j][k][l] * std::max(-FPd,0.0);

	vCDS[i][j][k] = FN*(v[i][j][k][l] + v[i][j+1][k][l])/2 - FPd*(v[i][j][k][l] + v[i][j-1][k][l])/2;

	vDC[i][j][k] = DCvalue*( vUDS[i][j][k] - vCDS[i][j][k]);
}

//Coefficients for W Equation for semi upwind scheme
void semiUpwindW()
{
	 att[i][j][k] = DT   + std::max(-FT,0.0);	//upwind
	 ab[i][j][k]  = DPdd + std::max(FPdd,0.0);	//upwind
	 atE[i][j][k] = Dted - Fted/2;				//central differencing
	 atW[i][j][k] = Dtw  + Ftw/2;				//central differencing
	 aNt[i][j][k] = Dntd - Fntd/2;				//central differencing
	 aSt[i][j][k] = Dst  + Fst/2;				//central differencing
	 
	 //Deferred Correction
	 wUDS[i][j][k] = w[i][j][k][l] * std::max(FT,0.0) - w[i][j][k+1][l] * std::max(-FT,0.0) - w[i][j][k-1][l] * std::max(FPdd,0.0) + w[i][j][k][l] * std::max(-FPdd,0.0);

	 wCDS[i][j][k] = FT*(w[i][j][k][l] + w[i][j][k+1][l])/2 - FPdd*(w[i][j][k][l] + w[i][j][k-1][l])/2;

	 wDC[i][j][k] = DCvalue*( wUDS[i][j][k] - wCDS[i][j][k]);
}


//*******************************************************
//Upwind schemes
//*******************************************************

//Coeffcients for U Equation for upwind scheme
void upwindU()
{
	aee[i][j][k] = std::max(-FE,0.0) + DE;  //upwind Scheme
	aw[i][j][k]  = std::max(FP,0.0)  + DP;  //upwind Scheme
	aNe[i][j][k] = std::max(-Fne,0.0)+ Dne; //upwind scheme
	aSe[i][j][k] = std::max(Fse,0.0) + Dse; //upwind scheme
	aTe[i][j][k] = std::max(-Fte,0.0)+ Dte; //upwind scheme
	aBe[i][j][k] = std::max(Fbe,0.0) + Dbe; //upwind scheme
	
	//Deferred Correction
    //uUDS[i][j][k] = u[i][j][k][l] * (std::max(FE,0.0) + std::max(-FP,0.0) + std::max(Fne,0.0) + std::max(-Fse,0.0) + std::max(Fte,0.0) + std::max(-Fbe,0.0)) - u[i+1][j][k][l] * std::max(-FE,0.0) - u[i-1][j][k][l] * std::max(FP,0.0) - u[i][j+1][k][l] * std::max(-Fne,0.0) - u[i][j-1][k][l] * std::max(Fse,0.0) - u[i][j][k+1][l] * std::max(-Fte,0.0) - u[i][j][k-1][l] * std::max(Fbe,0.0);                     ;

	//uCDS[i][j][k] = FE*(u[i][j][k][l] + u[i+1][j][k][l])/2 - FP*(u[i][j][k][l] + u[i-1][j][k][l])/2 + Fne*(u[i][j][k][l] + u[i][j+1][k][l])/2 - Fse*(u[i][j][k][l] + u[i][j-1][k][l])/2 + Fte*(u[i][j][k][l] + u[i][j][k+1][l])/2 - Fbe*(u[i][j][k][l] + u[i][j][k-1][l])/2;

	//uDC[i][j][k] = DCvalue*( uUDS[i][j][k] - uCDS[i][j][k]);
}

//Coeffcients for V Equation for upwind scheme
void upwindV()
{
	anE[i][j][k] = Dned+ std::max(-Fned,0.0);//upwind scheme
	anW[i][j][k] = Dnw + std::max(Fnw,0.0);  //upwind scheme
	ann[i][j][k] = DN  + std::max(-FN,0.0);  //upwind Scheme
	as[i][j][k]  = DPd + std::max(FPd,0.0);  //upwind Scheme
	anT[i][j][k] = Dnt + std::max(-Fnt,0.0); //upwind scheme
	anB[i][j][k] = Dnb + std::max(Fnb,0.0);  //upwind scheme
	
	//Deferred Correction
	//vUDS[i][j][k] = v[i][j][k][l] * (std::max(FN,0.0) + std::max(-FPd,0.0) + std::max(Fnt,0.0) + std::max(-Fnb,0.0) + std::max(Fned,0.0) + std::max(-Fnw,0.0)) - v[i][j+1][k][l] * std::max(-FN,0.0) - v[i][j-1][k][l] * std::max(FPd,0.0) - v[i][j][k+1][l] * std::max(-Fnt,0.0) - v[i][j][k-1][l] * std::max(Fnb,0.0) - v[i+1][j][k][l] * std::max(-Fned,0.0) - v[i-1][j][k][l] * std::max(Fnw,0.0);

	//vCDS[i][j][k] = FN*(v[i][j][k][l] + v[i][j+1][k][l])/2 - FPd*(v[i][j][k][l] + v[i][j-1][k][l])/2 + Fnt*(v[i][j][k][l] + v[i][j][k+1][l])/2 - Fnb*(v[i][j][k][l] + v[i][j][k-1][l])/2 + Fned*(v[i][j][k][l] + v[i+1][j][k][l])/2 - Fnw*(v[i][j][k][l] + v[i-1][j][k][l])/2;

	//vDC[i][j][k] = DCvalue*( vUDS[i][j][k] - vCDS[i][j][k]);

}

//Coeffcients for W Equation for upwind scheme
void upwindW()
{
	 att[i][j][k] = DT   + std::max(-FT,0.0);		//upwind 
	 ab[i][j][k]  = DPdd + std::max(FPdd,0.0);		//upwind
	 atE[i][j][k] = Dted + std::max(-Fted,0.0);		//upwind
	 atW[i][j][k] = Dtw  + std::max(Ftw,0.0);		//upwind
	 aNt[i][j][k] = Dntd + std::max(-Fntd,0.0);		//upwind
	 aSt[i][j][k] = Dst  + std::max(Fst,0.0);		//upwind
	 
	 //Deferred Correction
	//wUDS[i][j][k] = w[i][j][k][l] * (std::max(FT,0.0) + std::max(-FPdd,0.0) + std::max(Fted,0.0) + std::max(-Ftw,0.0) + std::max(Fntd,0.0) + std::max(-Fst,0.0)) - w[i][j][k+1][l] * std::max(-FT,0.0) - w[i][j][k-1][l] * std::max(FPdd,0.0) - w[i+1][j][k][l] * std::max(-Fted,0.0) - w[i-1][j][k][l] * std::max(Ftw,0.0) - w[i][j+1][k][l] * std::max(-Fntd,0.0) - w[i][j-1][k][l] * std::max(Fst,0.0);

    //	 wCDS[i][j][k] = FT*(w[i][j][k][l] + w[i][j][k+1][l])/2 - FPdd*(w[i][j][k][l] + w[i][j][k-1][l])/2 + Fted*(w[i][j][k][l] + w[i+1][j][k][l])/2 - Ftw*(w[i][j][k][l] + w[i-1][j][k][l])/2 + Fntd*(w[i][j][k][l] + w[i][j+1][k][l])/2 - Fst*(w[i][j][k][l] + w[i][j-1][k][l])/2;

	 //wDC[i][j][k] = DCvalue*( wUDS[i][j][k] - wCDS[i][j][k]);
}



//********************************************************
//Hybrid schemes
//*******************************************************


//Coefficients of  U Equation for hybrid scheme
void hybridU()
{
	//Coeff of the Equation
    if(DE==0)
     aee[i][j][k] = std::max(-FE,0.0);
    else
	{
		PE = FE / DE;
	    aee[i][j][k] = std::max(-FE,0.0) + DE*std::max(1 - fabs(PE)/2,0.0);   //hybrid Scheme	
	}
	 
    if(DP==0)
     aw[i][j][k]  = std::max(FP,0.0);
    else
    {
    	PP = FP / DP;
		aw[i][j][k]  = std::max(FP,0.0)  + DP*std::max(1 - fabs(PP)/2,0.0);   //hybrid Scheme	
    }
    
	if(Dne==0)
     aNe[i][j][k] = std::max(-Fne,0.0);
    else
    {
    	Pne= Fne/ Dne;
    	aNe[i][j][k] = std::max(-Fne,0.0)+ Dne*std::max(1 - fabs(Pne)/2,0.0); //hybrid Scheme
    }
	
	if(Dse==0)
     aSe[i][j][k] = std::max(Fse,0.0);
    else
    {
    	Pse= Fse/ Dse;
    	aSe[i][j][k] = std::max(Fse,0.0) + Dse*std::max(1 - fabs(Pse)/2,0.0); //hybrid Scheme
    }
    
	if(Dte==0)
     aTe[i][j][k] = std::max(-Fte,0.0);
    else
    {
    	Pte= Fte/ Dte;
    	aTe[i][j][k] = std::max(-Fte,0.0)+ Dte*std::max(1 - fabs(Pte)/2,0.0); //hybrid Scheme
    }
    
	if(Dbe==0)
     aBe[i][j][k] = std::max(Fbe,0.0);
    else
    {
    	Pbe= Fbe/ Dbe;
    	aBe[i][j][k] = std::max(Fbe,0.0) + Dbe*std::max(1 - fabs(Pbe)/2,0.0); //hybrid Scheme
    }
    
		
}

//Coefficients of  V Equation for hybrid scheme
void hybridV()
{

	if(Dned==0)
	 anE[i][j][k] = std::max(-Fned,0.0);
	else
	{
		Pned= Fned/ Dned;
		anE[i][j][k] = Dned*std::max(1 - fabs(Pned)/2,0.0) + std::max(-Fned,0.0); //hybrid scheme
	}
	 
	if(Dnw==0)
	 anW[i][j][k] = std::max(Fnw,0.0);
	else
	{
		Pnw = Fnw/ Dnw;
		anW[i][j][k] = Dnw *std::max(1 - fabs(Pnw)/2,0.0)  + std::max(Fnw,0.0);   //hybrid scheme
	}
	
	if(DN==0)
	 ann[i][j][k] = std::max(-FN,0.0);
	else
	{
		PN  = FN  / DN;
		ann[i][j][k] = DN  *std::max(1 - fabs(PN)/2,0.0)   + std::max(-FN,0.0);   //hybrid Scheme
	} 
	
	if(DPd==0)
	 as[i][j][k]  = std::max(FPd,0.0);
	else
	{
		PPd = FPd / DPd;
		as[i][j][k]  = DPd *std::max(1 - fabs(PPd)/2,0.0)  + std::max(FPd,0.0);   //hybrid Scheme
	} 
	
	if(Dnt==0)
	 anT[i][j][k] = std::max(-Fnt,0.0);
	else
	{
		Pnt = Fnt / Dnt;
		anT[i][j][k] = Dnt *std::max(1 - fabs(Pnt)/2,0.0)  + std::max(-Fnt,0.0);  //hybrid scheme	
	} 
	
	if(Dnb==0)
	 anB[i][j][k] = std::max(Fnb,0.0);
	else
	{
		Pnb = Fnb / Dnb;
		anB[i][j][k] = Dnb *std::max(1 - fabs(Pnb)/2,0.0)  + std::max(Fnb,0.0);   //hybrid scheme
	} 
	
	
}

////Coefficients of  W Equation for hybrid scheme
void hybridW()
{
	 
	 if(DT==0)
	  att[i][j][k] = std::max(-FT,0.0);    //to eliminate zero in denominator
	 else
	 {
	 	PT   = FT  / DT;
	 	att[i][j][k] = DT  *std::max(1 - 0.5*fabs(PT),0.0) + std::max(-FT,0.0);		//hybrid scheme
     }
	 
	 if(DPdd==0)
	  ab[i][j][k]  = std::max(FPdd,0.0);
	 else
	 {
	    PPdd = FPdd/ DPdd;
	    ab[i][j][k]  = DPdd*std::max(1 - 0.5*fabs(PPdd),0.0) + std::max(FPdd,0.0); 	//hybrid scheme		  	
	 } 
	 
	 if(Dted==0)
	  atE[i][j][k] = std::max(-Fted,0.0);
	 else
	 {
	 	Pted = Fted/ Dted;
	 	atE[i][j][k] = Dted*std::max(1 - 0.5*fabs(Pted),0.0) + std::max(-Fted,0.0);	//hybrid scheme
	 }
	 	 
	 if(Dtw==0)
      atW[i][j][k] = std::max(Ftw,0.0);
	 else
	 {
	 	Ptw  = Ftw / Dtw;
	 	atW[i][j][k] = Dtw *std::max(1 - 0.5*fabs(Ptw),0.0) + std::max(Ftw,0.0);
	 } 
	 
	 if(Dntd==0)
	  aNt[i][j][k] = std::max(-Fntd,0.0);
	 else
	 {
	 	Pntd = Fntd/ Dntd;
		aNt[i][j][k] = Dntd*std::max(1 - 0.5*fabs(Pntd),0.0) + std::max(-Fntd,0.0);
	 } 
	 
	 if(Dst==0)	
	  aSt[i][j][k] = std::max(Fst,0.0);
	 else
	 {
	 	 Pst  = Fst / Dst;
	 	 aSt[i][j][k] = Dst *std::max(1 - 0.5*fabs(Pst),0.0) + std::max(Fst,0.0);
	 } 
	 
 
}	


//***************************************************************
//Power Law Schemes
//***************************************************************

//Coefficients of  U Equation for hybrid scheme
void powerLawU()
{
	if(DE==0)
     aee[i][j][k] = std::max(-FE,0.0);
    else
	{
		PE = FE / DE;
	    aee[i][j][k] = std::max(-FE,0.0) + DE*std::max(pow(1 - 0.1*fabs(PE),5),0.0);   //hybrid Scheme	
	}
	 
    if(DP==0)
     aw[i][j][k]  = std::max(FP,0.0);
    else
    {
    	PP = FP / DP;
		aw[i][j][k]  = std::max(FP,0.0)  + DP*std::max(pow(1 - 0.1*fabs(PP),5),0.0);   //hybrid Scheme	
    }
    
	if(Dne==0)
     aNe[i][j][k] = std::max(-Fne,0.0);
    else
    {
    	Pne= Fne/ Dne;
    	aNe[i][j][k] = std::max(-Fne,0.0)+ Dne*std::max(pow(1 - 0.1*fabs(Pne),5),0.0); //hybrid Scheme
    }
	
	if(Dse==0)
     aSe[i][j][k] = std::max(Fse,0.0);
    else
    {
    	Pse= Fse/ Dse;
    	aSe[i][j][k] = std::max(Fse,0.0) + Dse*std::max(pow(1 - 0.1*fabs(Pse),5),0.0); //hybrid Scheme
    }
    
	if(Dte==0)
     aTe[i][j][k] = std::max(-Fte,0.0);
    else
    {
    	Pte= Fte/ Dte;
    	aTe[i][j][k] = std::max(-Fte,0.0)+ Dte*std::max(pow(1 - 0.1*fabs(Pte),5),0.0); //hybrid Scheme
    }
    
	if(Dbe==0)
     aBe[i][j][k] = std::max(Fbe,0.0);
    else
    {
    	Pbe= Fbe/ Dbe;
    	aBe[i][j][k] = std::max(Fbe,0.0) + Dbe*std::max(pow(1 - 0.1*fabs(Pbe),5),0.0); //hybrid Scheme
    }
    
		
}

//Coefficients of  V Equation for hybrid scheme
void powerLawV()
{

	if(Dned==0)
	 anE[i][j][k] = std::max(-Fned,0.0);
	else
	{
		Pned= Fned/ Dned;
		anE[i][j][k] = Dned*std::max(pow(1 - 0.1*fabs(Pned),5),0.0) + std::max(-Fned,0.0); //hybrid scheme
	}
	 
	if(Dnw==0)
	 anW[i][j][k] = std::max(Fnw,0.0);
	else
	{
		Pnw = Fnw/ Dnw;
		anW[i][j][k] = Dnw *std::max(pow(1 - 0.1*fabs(Pnw),5),0.0)  + std::max(Fnw,0.0);   //hybrid scheme
	}
	
	if(DN==0)
	 ann[i][j][k] = std::max(-FN,0.0);
	else
	{
		PN  = FN  / DN;
		ann[i][j][k] = DN  *std::max(pow(1 - 0.1*fabs(PN),5),0.0)   + std::max(-FN,0.0);   //hybrid Scheme
	} 
	
	if(DPd==0)
	 as[i][j][k]  = std::max(FPd,0.0);
	else
	{
		PPd = FPd / DPd;
		as[i][j][k]  = DPd *std::max(pow(1 - 0.1*fabs(PPd),5),0.0)  + std::max(FPd,0.0);   //hybrid Scheme
	} 
	
	if(Dnt==0)
	 anT[i][j][k] = std::max(-Fnt,0.0);
	else
	{
		Pnt = Fnt / Dnt;
		anT[i][j][k] = Dnt *std::max(pow(1 - 0.1*fabs(Pnt),5),0.0)  + std::max(-Fnt,0.0);  //hybrid scheme	
	} 
	
	if(Dnb==0)
	 anB[i][j][k] = std::max(Fnb,0.0);
	else
	{
		Pnb = Fnb / Dnb;
		anB[i][j][k] = Dnb *std::max(pow(1 - 0.1*fabs(Pnb),5),0.0)  + std::max(Fnb,0.0);   //hybrid scheme
	} 
	
	
}

////Coefficients of  W Equation for hybrid scheme
void powerLawW()
{
	 
	 if(DT==0)
	  att[i][j][k] = std::max(-FT,0.0);    //to eliminate zero in denominator
	 else
	 {
	 	PT   = FT  / DT;
	 	att[i][j][k] = DT  *std::max(pow(1 - 0.1*fabs(PT),5),0.0) + std::max(-FT,0.0);		//hybrid scheme
     }
	 
	 if(DPdd==0)
	  ab[i][j][k]  = std::max(FPdd,0.0);
	 else
	 {
	    PPdd = FPdd/ DPdd;
	    ab[i][j][k]  = DPdd*std::max(pow(1 - 0.1*fabs(PPdd),5),0.0) + std::max(FPdd,0.0); 	//hybrid scheme		  	
	 } 
	 
	 if(Dted==0)
	  atE[i][j][k] = std::max(-Fted,0.0);
	 else
	 {
	 	Pted = Fted/ Dted;
	 	atE[i][j][k] = Dted*std::max(pow(1 - 0.1*fabs(Pted),5),0.0) + std::max(-Fted,0.0);	//hybrid scheme
	 }
	 	 
	 if(Dtw==0)
      atW[i][j][k] = std::max(Ftw,0.0);
	 else
	 {
	 	Ptw  = Ftw / Dtw;
	 	atW[i][j][k] = Dtw *std::max(pow(1 - 0.1*fabs(Ptw),5),0.0) + std::max(Ftw,0.0);
	 } 
	 
	 if(Dntd==0)
	  aNt[i][j][k] = std::max(-Fntd,0.0);
	 else
	 {
	 	Pntd = Fntd/ Dntd;
		aNt[i][j][k] = Dntd*std::max(pow(1 - 0.1*fabs(Pntd),5),0.0) + std::max(-Fntd,0.0);
	 } 
	 
	 if(Dst==0)	
	  aSt[i][j][k] = std::max(Fst,0.0);
	 else
	 {
	 	 Pst  = Fst / Dst;
	 	 aSt[i][j][k] = Dst *std::max(pow(1 - 0.1*fabs(Pst),5),0.0) + std::max(Fst,0.0);
	 } 
	 
}	

	


	
