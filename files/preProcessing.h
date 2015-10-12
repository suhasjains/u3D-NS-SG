//*********************************************************
//PRE-PROCESSING FUNCTIONS
//**********************************************************

//Inputs the default value that is required for simulation
void defaultInput()
{

		xDomainLength=1;
		yDomainLength=1;
		zDomainLength=1;

		xNumberOfCells=20;
		yNumberOfCells=20;
		zNumberOfCells=20;

		startTime=0;

		endTime=100;

		nTimeSteps=1000;

		writeOutputSteps=10;

		allowableCourantNumber=10;

		rho=0.001;
		mu=0.000001;

		uTopWallVel=0;
		vTopWallVel=0;
		wTopWallVel=0;

		uBottomWallVel=0;
		vBottomWallVel=0;
		wBottomWallVel=0;

		uLeftWallVel=0;
		vLeftWallVel=0;
		wLeftWallVel=0;

		uRightWallVel=0;
		vRightWallVel=0;
		wRightWallVel=0;

		uFrontWallVel=0;
		vFrontWallVel=0;
		wFrontWallVel=0;

		uBackWallVel=1;
		vBackWallVel=0;
		wBackWallVel=0;

		pUnderRelaxCoeff=0.5;
		vUnderRelaxCoeff=0.5;

		xGrav=0;
		yGrav=0;
		zGrav=0;

		weight=0.5;
		DCvalue=1;



}

//dynamically allocates the size of array
void dynamicAllocate()
{

   	u	 = new double***[xNumberOfCells+5];
	v	 = new double***[xNumberOfCells+5];
	w	 = new double***[xNumberOfCells+5];
	finalU	 = new double***[xNumberOfCells+5];
	finalV	 = new double***[xNumberOfCells+5];
	finalW	 = new double***[xNumberOfCells+5];
	pseudoU	 = new double***[xNumberOfCells+5];
	pseudoV	 = new double***[xNumberOfCells+5];
	pseudoW	 = new double***[xNumberOfCells+5];
	p	 = new double***[xNumberOfCells+5];
	pCorr	 = new double***[xNumberOfCells+5];
	xCoord	 = new double**[xNumberOfCells+5];
	yCoord	 = new double**[xNumberOfCells+5];
	zCoord	 = new double**[xNumberOfCells+5];
	aee	 = new double**[xNumberOfCells+5];
	aw	 = new double**[xNumberOfCells+5];
	aNe	 = new double**[xNumberOfCells+5];
	aSe	 = new double**[xNumberOfCells+5];
	aTe	 = new double**[xNumberOfCells+5];
	aBe	 = new double**[xNumberOfCells+5];
	ae	 = new double**[xNumberOfCells+5];
	anE	 = new double**[xNumberOfCells+5];
	anW	 = new double**[xNumberOfCells+5];
	ann	 = new double**[xNumberOfCells+5];
	as	 = new double**[xNumberOfCells+5];
	anT	 = new double**[xNumberOfCells+5];
	anB	 = new double**[xNumberOfCells+5];
	an	 = new double**[xNumberOfCells+5];
	att	 = new double**[xNumberOfCells+5];
	ab	 = new double**[xNumberOfCells+5];
	atE	 = new double**[xNumberOfCells+5];
	atW	 = new double**[xNumberOfCells+5];
	aNt	 = new double**[xNumberOfCells+5];
	aSt	 = new double**[xNumberOfCells+5];
	at	 = new double**[xNumberOfCells+5];
	Qaee = new double**[xNumberOfCells+5];
	Qaw	 = new double**[xNumberOfCells+5];
	QaNe = new double**[xNumberOfCells+5];
	QaSe = new double**[xNumberOfCells+5];
	QaTe = new double**[xNumberOfCells+5];
	QaBe = new double**[xNumberOfCells+5];
	Qae	 = new double**[xNumberOfCells+5];
	QanE = new double**[xNumberOfCells+5];
	QanW = new double**[xNumberOfCells+5];
	Qann = new double**[xNumberOfCells+5];
	Qas	 = new double**[xNumberOfCells+5];
	QanT = new double**[xNumberOfCells+5];
	QanB = new double**[xNumberOfCells+5];
	Qan	 = new double**[xNumberOfCells+5];
	Qatt = new double**[xNumberOfCells+5];
	Qab	 = new double**[xNumberOfCells+5];
	QatE = new double**[xNumberOfCells+5];
	QatW = new double**[xNumberOfCells+5];
	QaNt = new double**[xNumberOfCells+5];
	QaSt = new double**[xNumberOfCells+5];
	Qat	 = new double**[xNumberOfCells+5];
	Qaeee= new double**[xNumberOfCells+5];
	Qaww = new double**[xNumberOfCells+5];
	QaNNe= new double**[xNumberOfCells+5];
	QaSSe= new double**[xNumberOfCells+5];
	QaTTe= new double**[xNumberOfCells+5];
	QaBBe= new double**[xNumberOfCells+5];
	QanEE= new double**[xNumberOfCells+5];
	QanWW= new double**[xNumberOfCells+5];
	Qannn= new double**[xNumberOfCells+5];
	Qass = new double**[xNumberOfCells+5];
	QanTT= new double**[xNumberOfCells+5];
	QanBB= new double**[xNumberOfCells+5];
	Qattt= new double**[xNumberOfCells+5];
	Qabb = new double**[xNumberOfCells+5];
	QatEE= new double**[xNumberOfCells+5];
	QatWW= new double**[xNumberOfCells+5];
	QaNNt= new double**[xNumberOfCells+5];
	QaSSt= new double**[xNumberOfCells+5];
	At	 = new double**[xNumberOfCells+5];
	An	 = new double**[xNumberOfCells+5];
	As	 = new double**[xNumberOfCells+5];
	Ae	 = new double**[xNumberOfCells+5];
	Aw	 = new double**[xNumberOfCells+5];
	Ab	 = new double**[xNumberOfCells+5];
	B	 = new double**[xNumberOfCells+5];
	BB	 = new double**[xNumberOfCells+5];
	BBB	 = new double**[xNumberOfCells+5];
	b	 = new double**[xNumberOfCells+5];
	bb	 = new double**[xNumberOfCells+5];
	bbb	 = new double**[xNumberOfCells+5];
	aE	 = new double**[xNumberOfCells+5];
	aW	 = new double**[xNumberOfCells+5];
	aN	 = new double**[xNumberOfCells+5];
	aS	 = new double**[xNumberOfCells+5];
	aT	 = new double**[xNumberOfCells+5];
	aB	 = new double**[xNumberOfCells+5];
	aP	 = new double**[xNumberOfCells+5];
	uUDS = new double**[xNumberOfCells+5];
	vUDS = new double**[xNumberOfCells+5];
	wUDS = new double**[xNumberOfCells+5];
	uCDS = new double**[xNumberOfCells+5];
	vCDS = new double**[xNumberOfCells+5];
	wCDS = new double**[xNumberOfCells+5];
	uDC  = new double**[xNumberOfCells+5];
	vDC  = new double**[xNumberOfCells+5];
	wDC  = new double**[xNumberOfCells+5];
	uRes = new double**[xNumberOfCells+5];
	vRes = new double**[xNumberOfCells+5];
	wRes = new double**[xNumberOfCells+5];
	Sm = new double**[xNumberOfCells+5];
	Sc = new double**[xNumberOfCells+5];
	uLower  = new double**[xNumberOfCells+5];
	vHigher = new double**[xNumberOfCells+5];
	wLower  = new double**[xNumberOfCells+5];
	uHigher = new double**[xNumberOfCells+5];
	vLower  = new double**[xNumberOfCells+5];
	wHigher = new double**[xNumberOfCells+5];




	for (int i=0; i<(xNumberOfCells+5); ++i)
	{
	  u[i]		 = new double**[yNumberOfCells+5];
	  v[i]		 = new double**[yNumberOfCells+5];
	  w[i]		 = new double**[yNumberOfCells+5];
	  finalU[i]	 = new double**[yNumberOfCells+5];
	  finalV[i]	 = new double**[yNumberOfCells+5];
	  finalW[i]	 = new double**[yNumberOfCells+5];
	  pseudoU[i]	 = new double**[yNumberOfCells+5];
	  pseudoV[i]	 = new double**[yNumberOfCells+5];
	  pseudoW[i]	 = new double**[yNumberOfCells+5];
	  p[i]		 = new double**[yNumberOfCells+5];
	  pCorr[i]	 = new double**[yNumberOfCells+5];
	  xCoord[i]	 = new double*[yNumberOfCells+5];
	  yCoord[i]	 = new double*[yNumberOfCells+5];
	  zCoord[i]	 = new double*[yNumberOfCells+5];
	  aee[i]	 = new double*[yNumberOfCells+5];
      aw[i]	     = new double*[yNumberOfCells+5];
	  aNe[i]	 = new double*[yNumberOfCells+5];
	  aSe[i]	 = new double*[yNumberOfCells+5];
	  aTe[i]	 = new double*[yNumberOfCells+5];
	  aBe[i]	 = new double*[yNumberOfCells+5];
	  ae[i]	     = new double*[yNumberOfCells+5];
	  anE[i]	 = new double*[yNumberOfCells+5];
	  anW[i]	 = new double*[yNumberOfCells+5];
	  ann[i]	 = new double*[yNumberOfCells+5];
	  as[i]	     = new double*[yNumberOfCells+5];
	  anT[i]	 = new double*[yNumberOfCells+5];
	  anB[i]	 = new double*[yNumberOfCells+5];
	  an[i]	     = new double*[yNumberOfCells+5];
	  att[i]	 = new double*[yNumberOfCells+5];
	  ab[i]	     = new double*[yNumberOfCells+5];
	  atE[i]	 = new double*[yNumberOfCells+5];
	  atW[i]	 = new double*[yNumberOfCells+5];
	  aNt[i]	 = new double*[yNumberOfCells+5];
	  aSt[i]	 = new double*[yNumberOfCells+5];
	  at[i]	     = new double*[yNumberOfCells+5];
	  Qaee[i]	 = new double*[yNumberOfCells+5];
      Qaw[i]	 = new double*[yNumberOfCells+5];
	  QaNe[i]	 = new double*[yNumberOfCells+5];
	  QaSe[i]	 = new double*[yNumberOfCells+5];
	  QaTe[i]	 = new double*[yNumberOfCells+5];
	  QaBe[i]	 = new double*[yNumberOfCells+5];
	  Qae[i]     = new double*[yNumberOfCells+5];
	  QanE[i]	 = new double*[yNumberOfCells+5];
	  QanW[i]	 = new double*[yNumberOfCells+5];
	  Qann[i]	 = new double*[yNumberOfCells+5];
	  Qas[i]     = new double*[yNumberOfCells+5];
	  QanT[i]	 = new double*[yNumberOfCells+5];
	  QanB[i]	 = new double*[yNumberOfCells+5];
	  Qan[i]     = new double*[yNumberOfCells+5];
	  Qatt[i]	 = new double*[yNumberOfCells+5];
	  Qab[i]     = new double*[yNumberOfCells+5];
	  QatE[i]	 = new double*[yNumberOfCells+5];
	  QatW[i]	 = new double*[yNumberOfCells+5];
	  QaNt[i]	 = new double*[yNumberOfCells+5];
	  QaSt[i]	 = new double*[yNumberOfCells+5];
	  Qat[i]     = new double*[yNumberOfCells+5];
	  Qaeee[i]	 = new double*[yNumberOfCells+5];
      Qaww[i]	 = new double*[yNumberOfCells+5];
	  QaNNe[i]	 = new double*[yNumberOfCells+5];
	  QaSSe[i]	 = new double*[yNumberOfCells+5];
	  QaTTe[i]	 = new double*[yNumberOfCells+5];
	  QaBBe[i]	 = new double*[yNumberOfCells+5];
	  QanEE[i]	 = new double*[yNumberOfCells+5];
	  QanWW[i]	 = new double*[yNumberOfCells+5];
	  Qannn[i]	 = new double*[yNumberOfCells+5];
	  Qass[i]     = new double*[yNumberOfCells+5];
	  QanTT[i]	 = new double*[yNumberOfCells+5];
	  QanBB[i]	 = new double*[yNumberOfCells+5];
	  Qattt[i]	 = new double*[yNumberOfCells+5];
	  Qabb[i]     = new double*[yNumberOfCells+5];
	  QatEE[i]	 = new double*[yNumberOfCells+5];
	  QatWW[i]	 = new double*[yNumberOfCells+5];
	  QaNNt[i]	 = new double*[yNumberOfCells+5];
	  QaSSt[i]	 = new double*[yNumberOfCells+5];
	  At[i]	     = new double*[yNumberOfCells+5];
	  An[i]	     = new double*[yNumberOfCells+5];
	  As[i]	     = new double*[yNumberOfCells+5];
	  Ae[i]	     = new double*[yNumberOfCells+5];
	  Aw[i]	     = new double*[yNumberOfCells+5];
	  Ab[i]	     = new double*[yNumberOfCells+5];
	  B[i]	     = new double*[yNumberOfCells+5];
	  BB[i]	     = new double*[yNumberOfCells+5];
	  BBB[i]	 = new double*[yNumberOfCells+5];
	  b[i]   	 = new double*[yNumberOfCells+5];
	  bb[i]	     = new double*[yNumberOfCells+5];
	  bbb[i]	 = new double*[yNumberOfCells+5];
	  aE[i]	     = new double*[yNumberOfCells+5];
	  aW[i]	     = new double*[yNumberOfCells+5];
	  aN[i]	     = new double*[yNumberOfCells+5];
	  aS[i]	     = new double*[yNumberOfCells+5];
	  aT[i]	     = new double*[yNumberOfCells+5];
	  aB[i]	     = new double*[yNumberOfCells+5];
	  aP[i]      = new double*[yNumberOfCells+5];
	  uUDS[i]    = new double*[yNumberOfCells+5];
	  vUDS[i]    = new double*[yNumberOfCells+5];
	  wUDS[i]    = new double*[yNumberOfCells+5];
	  uCDS[i]    = new double*[yNumberOfCells+5];
	  vCDS[i]    = new double*[yNumberOfCells+5];
	  wCDS[i]    = new double*[yNumberOfCells+5];
	  uDC[i]     = new double*[yNumberOfCells+5];
	  vDC[i]     = new double*[yNumberOfCells+5];
	  wDC[i]     = new double*[yNumberOfCells+5];
      uRes[i]    = new double*[yNumberOfCells+5];
	  vRes[i]    = new double*[yNumberOfCells+5];
	  wRes[i]    = new double*[yNumberOfCells+5];
	  Sc[i]    = new double*[yNumberOfCells+5];
	  Sm[i]    = new double*[yNumberOfCells+5];
	  uLower[i]    = new double*[yNumberOfCells+5];
	  vHigher[i]   = new double*[yNumberOfCells+5];
	  wLower[i]    = new double*[yNumberOfCells+5];
	  uHigher[i]   = new double*[yNumberOfCells+5];
	  vLower[i]    = new double*[yNumberOfCells+5];
	  wHigher[i]   = new double*[yNumberOfCells+5];

  	  for (int j=0; j<(yNumberOfCells+5); ++j)
	  {
		u[i][j]		 = new double*[zNumberOfCells+5];
		v[i][j]		 = new double*[zNumberOfCells+5];
		w[i][j]		 = new double*[zNumberOfCells+5];
	  	finalU[i][j]	 = new double*[zNumberOfCells+5];
	  	finalV[i][j]	 = new double*[zNumberOfCells+5];
	  	finalW[i][j]	 = new double*[zNumberOfCells+5];
		pseudoU[i][j]	 = new double*[zNumberOfCells+5];
	  	pseudoV[i][j]	 = new double*[zNumberOfCells+5];
	  	pseudoW[i][j]	 = new double*[zNumberOfCells+5];
	  	p[i][j]		 = new double*[zNumberOfCells+5];
	  	pCorr[i][j]	 = new double*[zNumberOfCells+5];
	  	xCoord[i][j]	 = new double[zNumberOfCells+5];
	  	yCoord[i][j]	 = new double[zNumberOfCells+5];
	  	zCoord[i][j]	 = new double[zNumberOfCells+5];
	  	aee[i][j]	 = new double[zNumberOfCells+5];
        aw[i][j]	 = new double[zNumberOfCells+5];
        aNe[i][j]	 = new double[zNumberOfCells+5];
        aSe[i][j]	 = new double[zNumberOfCells+5];
        aTe[i][j]	 = new double[zNumberOfCells+5];
        aBe[i][j]	 = new double[zNumberOfCells+5];
        ae[i][j]	 = new double[zNumberOfCells+5];
        anE[i][j]	 = new double[zNumberOfCells+5];
        anW[i][j]	 = new double[zNumberOfCells+5];
        ann[i][j]	 = new double[zNumberOfCells+5];
        as[i][j]	 = new double[zNumberOfCells+5];
        anT[i][j]	 = new double[zNumberOfCells+5];
        anB[i][j]	 = new double[zNumberOfCells+5];
        an[i][j]	 = new double[zNumberOfCells+5];
        att[i][j]	 = new double[zNumberOfCells+5];
        ab[i][j]	 = new double[zNumberOfCells+5];
        atE[i][j]	 = new double[zNumberOfCells+5];
        atW[i][j]	 = new double[zNumberOfCells+5];
        aNt[i][j]	 = new double[zNumberOfCells+5];
        aSt[i][j]	 = new double[zNumberOfCells+5];
        at[i][j]	 = new double[zNumberOfCells+5];
        Qaee[i][j]	 = new double[zNumberOfCells+5];
        Qaw[i][j]	 = new double[zNumberOfCells+5];
        QaNe[i][j]	 = new double[zNumberOfCells+5];
        QaSe[i][j]	 = new double[zNumberOfCells+5];
        QaTe[i][j]	 = new double[zNumberOfCells+5];
        QaBe[i][j]	 = new double[zNumberOfCells+5];
        Qae[i][j]	 = new double[zNumberOfCells+5];
        QanE[i][j]	 = new double[zNumberOfCells+5];
        QanW[i][j]	 = new double[zNumberOfCells+5];
        Qann[i][j]	 = new double[zNumberOfCells+5];
        Qas[i][j]	 = new double[zNumberOfCells+5];
        QanT[i][j]	 = new double[zNumberOfCells+5];
        QanB[i][j]	 = new double[zNumberOfCells+5];
        Qan[i][j]	 = new double[zNumberOfCells+5];
        Qatt[i][j]	 = new double[zNumberOfCells+5];
        Qab[i][j]	 = new double[zNumberOfCells+5];
        QatE[i][j]	 = new double[zNumberOfCells+5];
        QatW[i][j]	 = new double[zNumberOfCells+5];
        QaNt[i][j]	 = new double[zNumberOfCells+5];
        QaSt[i][j]	 = new double[zNumberOfCells+5];
        Qat[i][j]	 = new double[zNumberOfCells+5];
        Qaeee[i][j]	 = new double[zNumberOfCells+5];
        Qaww[i][j]	 = new double[zNumberOfCells+5];
        QaNNe[i][j]	 = new double[zNumberOfCells+5];
        QaSSe[i][j]	 = new double[zNumberOfCells+5];
        QaTTe[i][j]	 = new double[zNumberOfCells+5];
        QaBBe[i][j]	 = new double[zNumberOfCells+5];
        QanEE[i][j]	 = new double[zNumberOfCells+5];
        QanWW[i][j]	 = new double[zNumberOfCells+5];
        Qannn[i][j]	 = new double[zNumberOfCells+5];
        Qass[i][j]	 = new double[zNumberOfCells+5];
        QanTT[i][j]	 = new double[zNumberOfCells+5];
        QanBB[i][j]	 = new double[zNumberOfCells+5];
        Qattt[i][j]	 = new double[zNumberOfCells+5];
        Qabb[i][j]	 = new double[zNumberOfCells+5];
        QatEE[i][j]	 = new double[zNumberOfCells+5];
        QatWW[i][j]	 = new double[zNumberOfCells+5];
        QaNNt[i][j]	 = new double[zNumberOfCells+5];
        QaSSt[i][j]	 = new double[zNumberOfCells+5];
        At[i][j]	 = new double[zNumberOfCells+5];
        An[i][j]	 = new double[zNumberOfCells+5];
        As[i][j]	 = new double[zNumberOfCells+5];
        Ae[i][j]	 = new double[zNumberOfCells+5];
        Aw[i][j]	 = new double[zNumberOfCells+5];
        Ab[i][j]	 = new double[zNumberOfCells+5];
        B[i][j]	     = new double[zNumberOfCells+5];
        BB[i][j]	 = new double[zNumberOfCells+5];
        BBB[i][j]	 = new double[zNumberOfCells+5];
        b[i][j]      = new double[zNumberOfCells+5];
        bb[i][j]	 = new double[zNumberOfCells+5];
        bbb[i][j]	 = new double[zNumberOfCells+5];
        aE[i][j]	 = new double[zNumberOfCells+5];
        aW[i][j]	 = new double[zNumberOfCells+5];
        aN[i][j]	 = new double[zNumberOfCells+5];
        aS[i][j]	 = new double[zNumberOfCells+5];
        aT[i][j]	 = new double[zNumberOfCells+5];
        aB[i][j]	 = new double[zNumberOfCells+5];
        aP[i][j]	 = new double[zNumberOfCells+5];
        uUDS[i][j]   = new double[zNumberOfCells+5];
	    vUDS[i][j]   = new double[zNumberOfCells+5];
	    wUDS[i][j]   = new double[zNumberOfCells+5];
	    uCDS[i][j]   = new double[zNumberOfCells+5];
	    vCDS[i][j]   = new double[zNumberOfCells+5];
	    wCDS[i][j]   = new double[zNumberOfCells+5];
	    uDC[i][j]    = new double[zNumberOfCells+5];
	    vDC[i][j]    = new double[zNumberOfCells+5];
	    wDC[i][j]    = new double[zNumberOfCells+5];
	    uRes[i][j]   = new double[zNumberOfCells+5];
	    vRes[i][j]   = new double[zNumberOfCells+5];
	    wRes[i][j]   = new double[zNumberOfCells+5];
	    Sm[i][j]   = new double[zNumberOfCells+5];
	    Sc[i][j]   = new double[zNumberOfCells+5];
	    uLower[i][j]   = new double[zNumberOfCells+5];
	    vHigher[i][j]   = new double[zNumberOfCells+5];
	    wLower[i][j]   = new double[zNumberOfCells+5];
	    uHigher[i][j]   = new double[zNumberOfCells+5];
	    vLower[i][j]   = new double[zNumberOfCells+5];
	    wHigher[i][j]   = new double[zNumberOfCells+5];

		 for (int k=0; k<(zNumberOfCells+5); ++k)
		 {
		  u[i][j][k] 		 = new double[2];
		  v[i][j][k]		 = new double[2];
	 	  w[i][j][k]		 = new double[2];
	  	  finalU[i][j][k]	 = new double[2];
	  	  finalV[i][j][k]	 = new double[2];
	  	  finalW[i][j][k]	 = new double[2];
		  pseudoU[i][j][k]	 = new double[2];
	  	  pseudoV[i][j][k]	 = new double[2];
	  	  pseudoW[i][j][k]	 = new double[2];
	  	  p[i][j][k]		 = new double[2];
	  	  pCorr[i][j][k]	 = new double[2];
	 	 }

	  }
	}

}

//Finalizes things inside simple loop
void finalize()
{

	  stop = clock();
	  realTime = (double) (stop-start)/CLOCKS_PER_SEC;
          printf("Run time: %2.2f seconds \n", realTime);				//printing real time

}

//Initializes things inside simple loop
void initialize()
{

	l = 1;   //making the simulation run at l=1 so that previous time step values are stored in l=0

	nSimpleLoops = 1;  //making initial no of simple loops = 3
	nSimplerLoops = 1; //making initial no of simpler loops = 1

	printf("\n\n\nsimulationTime :%2.10lf seconds\n",simulationTime);
}

//Initializing number of loops
void nLoops()
{


	nVelocityLoops = 1;
	nPressureLoops = 1;

}

//Gives initial values to all the variables and fields at time 0
void zeroTimeValues()
{

	//INPUT OUTPUT VARIABLES

	nx = ny = nz = n = 0;
	deltaT = simulationTime = 0;
	adjustTimeStep = 0;
	topBoundary = bottomBoundary = leftBoundary = rightBoundary = frontBoundary = backBoundary = 0;


	//COEFFICIENTS OF VELOCITY EQUATION
	dx = dy = dz = dt = 0;


	Fte = Fbe = FE = FP = Fne = Fse = 0;
	Fned = Fnw = FN = FPd = Fnt = Fnb = 0;
	Fted = Ftw = Fntd = Fst = FT = FPdd = 0;

	DE = DP = Dne = Dse = Dte = Dbe = 0;
	Dned = Dnw = DN = DPd = Dnt = Dnb = 0;
	DT = DPdd = Dted = Dtw = Dntd = Dst = 0;

	//COEFFICIENTS OF PRESSURE CORRECTION EQUATION


	maxMassResidual=0;
	uMaxRes = vMaxRes = wMaxRes = 0;




	//VELOCITIES and PRESSURE
	for   (l=0; l<2; l++)
	 for  (k=0; k<(zNumberOfCells+5); k++)
	  for (j=0; j<(yNumberOfCells+5); j++)
	   for(i=0; i<(xNumberOfCells+5); i++)
	{
		u[i][j][k][l]	  = 0;
		v[i][j][k][l]	  = 0;
		w[i][j][k][l] 	  = 0;

		finalU[i][j][k][l]	  = 0;
		finalV[i][j][k][l]	  = 0;
		finalW[i][j][k][l] 	  = 0;

		pseudoU[i][j][k][l]	  = 0;
		pseudoV[i][j][k][l]	  = 0;
		pseudoW[i][j][k][l] 	  = 0;

		p[i][j][k][l]     = 0;
		pCorr[i][j][k][l] = 0;

		xCoord[i][j][k] = 0;
		yCoord[i][j][k] = 0;
		zCoord[i][j][k] = 0;

		aee[i][j][k] = aw[i][j][k] = aNe[i][j][k] = aSe[i][j][k] = aTe[i][j][k] = aBe[i][j][k] = ae[i][j][k] = 0;
	    anE[i][j][k] = anW[i][j][k] = ann[i][j][k] = as[i][j][k] = anT[i][j][k] = anB[i][j][k] = an[i][j][k] = 0;
        att[i][j][k] = ab[i][j][k] = atE[i][j][k] = atW[i][j][k] = aNt[i][j][k] = aSt[i][j][k] = at[i][j][k] = 0;

        Ae[i][j][k] = An[i][j][k] = At[i][j][k] = Aw[i][j][k] = As[i][j][k] = Ab[i][j][k] = 0;
        B[i][j][k] = BB[i][j][k] = BBB[i][j][k] = 0;

        aE[i][j][k] = aW[i][j][k] = aN[i][j][k] = aS[i][j][k] = aT[i][j][k] = aB[i][j][k] = aP[i][j][k] = 0;
        b[i][j][k] = bb[i][j][k] = bbb[i][j][k] = 0;

        Sc[i][j][k] = 0;
	    Sm[i][j][k] = 0;
	    uRes[i][j][k] = vRes[i][j][k] = wRes[i][j][k] = 0;

	}


	for(i=0; i<(xNumberOfCells+5); i++)
	{
		lower[i] =  0;
		upper[i] = 0;
		c[i] = 0;
		x[i] = 0;
		m[i] = 0;
		cc[i] = 0;
	}

	for(i=0; i<2; i++)
	{
		uMaxResidual[i]=0;
		vMaxResidual[i]=0;
		wMaxResidual[i]=0;
	}

	i = j = k = l = t = thomasI = simpleVar = nonLinear = 0;   //indices

	nx = ny = nz = n = 0;                                  //Number of Cells

	uMax = vMax = wMax = maxCourantNumber = 0;

	check = f0 = f1 = f2 = f3 = f4 = f5 = f6 = f7 = 0;

}

//Updates the values of the boundary cells
void boundaryConditions()
{

	//velocity boundary conditions
    for  (k=1; k<=nz+2; k++)
	 for (i=1; i<=nx+2; i++)
	  for(j=1; j<=ny+2; j++)
	{

	//u[1][j][k][l]    = u[2][j][k][l];					//u zero gradient along x on left
	u[1][j][k][l]    = uLeftWallVel;							//u[1][j][k][l]    = uLeftWallVel;
	 //u[nx+1][j][k][l] = u[nx][j][k][l];					//u zero gradient along x on right
	u[nx+1][j][k][l] = uRightWallVel;

	 v[i][1][k][l]    = vFrontWallVel;					//v[i][1][k][l] = v[i][2][k[l];
	 v[i][ny+1][k][l] = vBackWallVel;					//v[i][ny+1][k][l] = v[i][ny][k][l];

	//w[i][j][1][l]    = wBottomWallVel;
	w[i][j][1][l]    = w[i][j][2][l];								//
	w[i][j][nz+1][l] = w[i][j][nz][l];
	//w[i][j][nz+1][l] = wTopWallVel;					//								//

	}

/*


	for  (k=1; k<=nz+2; k++)
	 for (i=1; i<=nx+2; i++)
	  for(j=((ny+2)/2); j<=ny+2; j++)                                    //backward Facing Step Velocity
	{

	 u[1][j][k][l]    = 3 * ((j/(ny+2)) - 1) * ((j/(ny+2))- 0.5);



	}
*/

}

//Updates the values of the pressure and velocity ghost cells
void ghostCells()
{
	//velocity ghost cells
	for  (k=1; k<=nz+2; k++)
	 for (i=1; i<=nx+2; i++)
	  for(j=1; j<=ny+2; j++)
	{

	 //u[i][1][k][l]    = u[i][2][k][l];
	u[i][1][k][l]    = 2 * uFrontWallVel  - u[i][2][k][l];
	 //u[i][ny+2][k][l] = u[i][ny+1][k][l];
	u[i][ny+2][k][l] = 2 * uBackWallVel   - u[i][ny+1][k][l] ;
	u[i][j][1][l]    = u[i][j][2][l];
	//u[i][j][1][l]    = 2 * uBottomWallVel - u[i][j][2][l];
	u[i][j][nz+2][l] = u[i][j][nz+1][l];
	//u[i][j][nz+2][l] = 2 * uTopWallVel    - u[i][j][nz+1][l];

	//v[1][j][k][l]    = v[2][j][k][l];					//v zero gradient along x on left
	v[1][j][k][l]    = 2 * vLeftWallVel   - v[2][j][k][l];
	//v[nx+2][j][k][l] = v[nx+1][j][k][l];					//v zero gradient along x on right
	v[nx+2][j][k][l] = 2 * vRightWallVel  - v[nx+1][j][k][l];
	v[i][j][1][l]    = v[i][j][2][l];
	//v[i][j][1][l]    = 2 * vBottomWallVel - v[i][j][2][l];
	v[i][j][nz+2][l] = v[i][j][nz+1][l];
	//v[i][j][nz+2][l] = 2 * vTopWallVel    - v[i][j][nz+1][l];

	//w[1][j][k][l]    = w[2][j][k][l];					//w zero gradient along x on left
	w[1][j][k][l]    = 2 * wLeftWallVel   - w[2][j][k][l];
	//w[nx+2][j][k][l] = w[nx+1][j][k][l];					//w zero gradient along x on right
 	w[nx+2][j][k][l] = 2 * wRightWallVel  - w[nx+1][j][k][l];
	//w[i][1][k][l]    = w[i][2][k][l];
	w[i][1][k][l]    = 2 * wFrontWallVel  - w[i][2][k][l];
	//w[i][ny+2][k][l] = w[i][ny+1][k][l];
	w[i][ny+2][k][l] = 2 * wBackWallVel   - w[i][ny+1][k][l];

	}


	//pressure ghost cells
	for  (k=1; k<=nz+2; k++)
	 for (i=1; i<=nx+2; i++)
	  for(j=1; j<=ny+2; j++)
	{
	 //pCorr[1][j][k][l]    =0;
	 pCorr[1][j][k][l]    = pCorr[2][j][k][l];
	 pCorr[i][1][k][l]    = pCorr[i][2][k][l];
	 pCorr[i][j][1][l]    = pCorr[i][j][2][l];
	 //pCorr[nx+2][j][k][l] = 0;
	 pCorr[nx+2][j][k][l] = pCorr[nx+1][j][k][l];
	 pCorr[i][ny+2][k][l] = pCorr[i][ny+1][k][l];
	 pCorr[i][j][nz+2][l] = pCorr[i][j][nz+1][l];

	 //p[1][j][k][l]	  = 0;
	 //p[1][j][k][l]    = p[2][j][k][l];
	 //p[i][1][k][l]    = p[i][2][k][l];
	 //p[i][j][1][l]    = p[i][j][2][l];
	 //p[nx+2][j][k][l] = 0;
	 //p[nx+2][j][k][l] = p[nx+1][j][k][l];
	 //p[i][ny+2][k][l] = p[i][ny+1][k][l];
	 //p[i][j][nz+2][l] = p[i][j][nz+1][l];

	}



}


