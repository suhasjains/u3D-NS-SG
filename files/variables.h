//Declaration of all variables


#ifndef VARIABLES_H
#define VARIABLES_H


	//INPUT OUTPUT VARIABLES
	double xDomainLength,yDomainLength,zDomainLength;
	int   xNumberOfCells,yNumberOfCells,zNumberOfCells,nx,ny,nz,n;   //number of cells
	int   nTimeSteps;
	double startTime,endTime,deltaT,simulationTime,realTime,start,stop;
	char  adjustTimeStep;
	char  topBoundary,bottomBoundary,leftBoundary,rightBoundary,frontBoundary,backBoundary;
	double uTopWallVel,uBottomWallVel,uLeftWallVel,uRightWallVel,uFrontWallVel,uBackWallVel;
	double vTopWallVel,vBottomWallVel,vLeftWallVel,vRightWallVel,vFrontWallVel,vBackWallVel;
	double wTopWallVel,wBottomWallVel,wLeftWallVel,wRightWallVel,wFrontWallVel,wBackWallVel;


	//COEFFICIENTS OF VELOCITY EQUATION
	double dx,dy,dz,dt;
	double ***b,***bb,***bbb;
	double xGrav,yGrav,zGrav;

	double Fte,Fbe,FE,FP,Fne,Fse;
	double Fned,Fnw,FN,FPd,Fnt,Fnb;
	double Fted,Ftw,Fntd,Fst,FT,FPdd;

	double DE,DP,Dne,Dse,Dte,Dbe;
	double Dned,Dnw,DN,DPd,Dnt,Dnb;
	double DT,DPdd,Dted,Dtw,Dntd,Dst;
	
	double PE,PP,Pne,Pse,Pte,Pbe;
	double Pned,Pnw,PN,PPd,Pnt,Pnb;
	double PT,PPdd,Pted,Ptw,Pntd,Pst;

	double ***aee,***aw,***aNe,***aSe,***aTe,***aBe,***ae;
	double ***anE,***anW,***ann,***as,***anT,***anB,***an;
	double ***att,***ab,***atE,***atW,***aNt,***aSt,***at;

	double ***Ae,***An,***At,***Aw,***As,***Ab;
	double ***B,***BB,***BBB;


	//COEFFICIENTS OF PRESSURE CORRECTION EQUATION
	double ***aE,***aW,***aN,***aS,***aT,***aB,***aP;
	double ***Sm,***Sc,maxMassResidual;
	double uMaxRes,vMaxRes,wMaxRes,***uRes,***vRes,***wRes,uMaxResidual[2],vMaxResidual[2],wMaxResidual[2];
	int pMassResidual,nMassResidual,nsMassResidual,psMassResidual,nuMaxRes,nvMaxRes,nwMaxRes,puMaxRes,pvMaxRes,pwMaxRes;

	//VELOCITIES and PRESSURE
	double ****u,****v,****w;
	double ****finalU,****finalV,****finalW;
	double ****p,****pCorr;
	double ***xCoord,***yCoord,***zCoord;
	double ****pseudoW,****pseudoV,****pseudoU;

	double uMax,vMax,wMax;										//max velocities
	double maxCourantNumber,allowableCourantNumber;									//max courant number
	int nSimpleLoops,nVelocityLoops,nPressureLoops,nSimplerLoops;

	int i,j,k,l,t,courantNumberControl,thomasI,simpleVar,simplerVar,nonLinear;         		//indices
	double lower[250],upper[250],c[250],x[250],m[250],cc[250],middle[250];       //MATRIX A,X AND B

	double pUnderRelaxCoeff,vUnderRelaxCoeff;                               	// under relaxation coefficient for pressure

	double rho,mu; 						//density and viscosity


	FILE *fp; 						//file pointer
	char f0,f1,f2,f3,f4,f5,f6,f7;
	char str[100],check;
	double accuracy;

	int writeOutputSteps;
	double weight;


	//Deferred correction
	double ***uUDS, ***vUDS, ***wUDS, ***uCDS, ***vCDS, ***wCDS, ***uDC, ***vDC, ***wDC, DCvalue;
	double ***uLower,***uHigher,***vLower,***vHigher,***wLower,***wHigher;
	
	double alphaP,alphaE,alphase,alphane,alphabe,alphate;
	double alphanw,alphaned,alphaPd,alphaN,alphanb,alphant;
	double alphaPdd,alphaT,alphated,alphatw,alphantd,alphast;
	
	double ***Qaee,***Qaw,***QaNe,***QaSe,***QaTe,***QaBe,***Qae,***Qaeee,***Qaww,***QaNNe,***QaSSe,***QaTTe,***QaBBe;
	double ***QanE,***QanW,***Qann,***Qas,***QanT,***QanB,***Qan,***QanEE,***QanWW,***Qannn,***Qass,***QanTT,***QanBB;
	double ***Qatt,***Qab,***QatE,***QatW,***QaNt,***QaSt,***Qat,***Qattt,***Qabb,***QatEE,***QatWW,***QaNNt,***QaSSt;
	
	

	double any;



#endif

