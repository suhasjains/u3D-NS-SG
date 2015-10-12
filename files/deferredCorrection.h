void quickU()
{
	if(FP>0) 	alphaP=1;
	else	 	alphaP=0;
	if(FE>0)	alphaE=1;
	else		alphaE=0;
	if(Fse>0) 	alphase=1;
	else	 	alphase=0;
	if(Fne>0) 	alphane=1;
	else	 	alphane=0;
	if(Fbe>0) 	alphabe=1;
	else	 	alphabe=0;
	if(Fte>0) 	alphate=1;
	else	 	alphate=0;
	
	Qaee[i][j][k]	= DE - 3*alphaE*FE/8 - 6*(1-alphaE)*FE/8 - (1-alphaP)*FP/8;
	Qaw[i][j][k] 	= DP + 6*alphaP*FP/8 + alphaE*FE/8 + 3*(1-alphaP)*FP/8;
	Qaww[i][j][k]	= -alphaP*FP/8;
	Qaeee[i][j][k]	= (1-alphaE)*FE/8;
	
	QaNe[i][j][k]	= Dne - 3*alphane*Fne/8 - 6*(1-alphane)*Fne/8 - (1-alphase)*Fse/8;
	QaSe[i][j][k] 	= Dse + 6*alphase*Fse/8 + alphane*Fne/8 + 3*(1-alphase)*Fse/8;
	QaSSe[i][j][k]	= -alphase*Fse/8;
	QaNNe[i][j][k]	= (1-alphane)*Fne/8;
	
	QaTe[i][j][k]	= Dte - 3*alphate*Fte/8 - 6*(1-alphate)*Fte/8 - (1-alphabe)*Fbe/8;
	QaBe[i][j][k] 	= Dbe + 6*alphabe*Fbe/8 + alphate*Fte/8 + 3*(1-alphabe)*Fbe/8;
	QaBBe[i][j][k]	= -alphabe*Fbe/8;
	QaTTe[i][j][k]	= (1-alphate)*Fte/8;
	
	Qae[i][j][k]	= Qaee[i][j][k] + Qaw[i][j][k] + Qaww[i][j][k] + Qaeee[i][j][k] + QaNe[i][j][k] + QaSe[i][j][k] + QaSSe[i][j][k] + QaNNe[i][j][k] + QaTe[i][j][k] + QaBe[i][j][k] + QaBBe[i][j][k] + QaTTe[i][j][k];

    uLower[i][j][k]	= (ae[i][j][k] - rho*dx*dy*dz/dt)*u[i][j][k][l] - (aee[i][j][k]*u[i+1][j][k][l] + aw[i][j][k]*u[i-1][j][k][l] + aNe[i][j][k]*u[i][j+1][k][l] + aSe[i][j][k]*u[i][j-1][k][l] + aTe[i][j][k]*u[i][j][k+1][l] + aBe[i][j][k]*u[i][j][k-1][l]);
    
	uHigher[i][j][k]= Qae[i][j][k]*u[i][j][k][l] - (Qaee[i][j][k]*u[i+1][j][k][l] + Qaw[i][j][k]*u[i-1][j][k][l] + Qaww[i][j][k]*u[i-2][j][k][l] + Qaeee[i][j][k]*u[i+2][j][k][l] + QaNe[i][j][k]*u[i][j+1][k][l] + QaSe[i][j][k]*u[i][j-1][k][l] + QaSSe[i][j][k]*u[i][j-2][k][l] + QaNNe[i][j][k]*u[i][j+2][k][l] + QaTe[i][j][k]*u[i][j][k+1][l] + QaBe[i][j][k]*u[i][j][k-1][l] + QaBBe[i][j][k]*u[i][j][k-2][l] + QaTTe[i][j][k]*u[i][j][k+2][l]);		 
	
	uDC[i][j][k]	= DCvalue*(uLower[i][j][k] - uHigher[i][j][k]); 
	 	
}

void quickV()
{
	if(Fnw>0) 	alphanw=1;
	else	 	alphanw=0;
	if(Fned>0)	alphaned=1;
	else		alphaned=0;
	if(FPd>0) 	alphaPd=1;
	else	 	alphaPd=0;
	if(FN>0) 	alphaN=1;
	else	 	alphaN=0;
	if(Fnb>0) 	alphanb=1;
	else	 	alphanb=0;
	if(Fnt>0) 	alphant=1;
	else	 	alphant=0;
	
	QanE[i][j][k]	= Dned - 3*alphaned*Fned/8 - 6*(1-alphaned)*Fned/8 - (1-alphanw)*Fnw/8;
	QanW[i][j][k] 	= Dnw + 6*alphanw*Fnw/8 + alphaned*Fned/8 + 3*(1-alphanw)*Fnw/8;
	QanWW[i][j][k]	= -alphanw*Fnw/8;
	QanEE[i][j][k]	= (1-alphaned)*Fned/8;
	
	Qann[i][j][k]	= DN - 3*alphaN*FN/8 - 6*(1-alphaN)*FN/8 - (1-alphaPd)*FPd/8;
	Qas[i][j][k] 	= DPd + 6*alphaPd*FPd/8 + alphaN*FN/8 + 3*(1-alphaPd)*FPd/8;
	Qass[i][j][k]	= -alphaPd*FPd/8;
	Qannn[i][j][k]	= (1-alphaN)*FN/8;
	
	QanT[i][j][k]	= Dnt - 3*alphant*Fnt/8 - 6*(1-alphant)*Fnt/8 - (1-alphanb)*Fnb/8;
	QanB[i][j][k] 	= Dnb + 6*alphanb*Fnb/8 + alphant*Fnt/8 + 3*(1-alphanb)*Fnb/8;
	QanBB[i][j][k]	= -alphanb*Fnb/8;
	QanTT[i][j][k]	= (1-alphant)*Fnt/8;
	
	Qan[i][j][k]	= QanE[i][j][k] + QanW[i][j][k] + QanEE[i][j][k] + QanWW[i][j][k] + Qann[i][j][k] + Qas[i][j][k] + Qannn[i][j][k] + Qass[i][j][k] + QanT[i][j][k] + QanB[i][j][k] + QanTT[i][j][k] + QanBB[i][j][k];

    vLower[i][j][k]	= (an[i][j][k] - rho*dx*dy*dz/dt)*v[i][j][k][l] - (anE[i][j][k]*v[i+1][j][k][l] + anW[i][j][k]*v[i-1][j][k][l] + ann[i][j][k]*v[i][j+1][k][l] + as[i][j][k]*v[i][j-1][k][l] + anT[i][j][k]*v[i][j][k+1][l] + anB[i][j][k]*v[i][j][k-1][l]); 
    
	vHigher[i][j][k]= Qan[i][j][k]*v[i][j][k][l] - (QanE[i][j][k]*v[i+1][j][k][l] + QanW[i][j][k]*v[i-1][j][k][l] + QanWW[i][j][k]*v[i-2][j][k][l] + QanEE[i][j][k]*v[i+2][j][k][l] + Qann[i][j][k]*v[i][j+1][k][l] + Qas[i][j][k]*v[i][j-1][k][l] + Qass[i][j][k]*v[i][j-2][k][l] + Qannn[i][j][k]*v[i][j+2][k][l] + QanT[i][j][k]*v[i][j][k+1][l] + QanB[i][j][k]*v[i][j][k-1][l] + QanBB[i][j][k]*v[i][j][k-2][l] + QanTT[i][j][k]*v[i][j][k+2][l]);		 
	
	vDC[i][j][k]	= DCvalue*(vLower[i][j][k] - vHigher[i][j][k]); 
	 
}


void quickW()
{
	if(FPdd>0) 	alphaPdd=1;
	else	 	alphaPdd=0;
	if(FT>0)	alphaT=1;
	else		alphaT=0;
	if(Fted>0) 	alphated=1;
	else	 	alphated=0;
	if(Ftw>0) 	alphatw=1;
	else	 	alphatw=0;
	if(Fntd>0) 	alphantd=1;
	else	 	alphantd=0;
	if(Fst>0) 	alphast=1;
	else	 	alphast=0;
	
	Qatt[i][j][k]	= DT - 3*alphaT*FT/8 - 6*(1-alphaT)*FT/8 - (1-alphaPdd)*FPdd/8;
	Qab[i][j][k] 	= DPdd + 6*alphaPdd*FPdd/8 + alphaT*FT/8 + 3*(1-alphaPdd)*FPdd/8;
	Qabb[i][j][k]	= -alphanw*FPdd/8;
	Qattt[i][j][k]	= (1-alphaT)*FT/8;
	
	QatE[i][j][k]	= Dted - 3*alphated*Fted/8 - 6*(1-alphated)*Fted/8 - (1-alphatw)*Ftw/8;
	QatW[i][j][k] 	= Dtw + 6*alphatw*Ftw/8 + alphated*Fted/8 + 3*(1-alphatw)*Ftw/8;
	QatWW[i][j][k]	= -alphatw*Ftw/8;
	QatEE[i][j][k]	= (1-alphated)*Fted/8;
	
	QaNt[i][j][k]	= Dntd - 3*alphantd*Fntd/8 - 6*(1-alphantd)*Fntd/8 - (1-alphast)*Fst/8;
	QaSt[i][j][k] 	= Dst + 6*alphast*Fst/8 + alphantd*Fntd/8 + 3*(1-alphast)*Fst/8;
	QaSSt[i][j][k]	= -alphast*Fst/8;
	QaNNt[i][j][k]	= (1-alphantd)*Fntd/8;
	
	Qat[i][j][k]	= Qatt[i][j][k] + Qab[i][j][k] + Qabb[i][j][k] + Qattt[i][j][k] + QatE[i][j][k] + QatW[i][j][k] + QatWW[i][j][k] + QatEE[i][j][k] + QaNt[i][j][k] + QaSt[i][j][k] + QaNNt[i][j][k] + QaSSt[i][j][k];

    wLower[i][j][k]	= (at[i][j][k] - rho*dx*dy*dz/dt)*w[i][j][k][l] - (att[i][j][k]*w[i][j][k+1][l] + ab[i][j][k]*w[i][j][k-1][l] + atE[i][j][k]*w[i+1][j][k][l] + atW[i][j][k]*w[i-1][j][k][l] + aNt[i][j][k]*w[i][j+1][k][l] + aSt[i][j][k]*w[i][j-1][k][l]); 
    
	wHigher[i][j][k]= Qat[i][j][k]*w[i][j][k][l] - (Qatt[i][j][k]*w[i][j][k+1][l] + Qab[i][j][k]*w[i][j][k-1][l] + Qabb[i][j][k]*w[i][j][k-2][l] + Qattt[i][j][k]*w[i][j][k+2][l] + QatE[i][j][k]*w[i+1][j][k][l] + QatW[i][j][k]*w[i-1][j][k][l] + QatWW[i][j][k]*w[i-2][j][k][l] + QatEE[i][j][k]*w[i+2][j][k][l] + QaNt[i][j][k]*w[i][j+1][k][l] + QaSt[i][j][k]*w[i][j-1][k][l] + QaNNt[i][j][k]*w[i][j+2][k][l] + QaSSt[i][j][k]*w[i][j-2][k][l]);
	
	wDC[i][j][k]	= DCvalue*(wLower[i][j][k] - wHigher[i][j][k]); 
	 
}
