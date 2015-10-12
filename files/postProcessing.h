//writes internal mesh to mesh.csv which can be opened using paraview
void writeMesh()
{

	stepSize();

	fp=fopen("/home/suhas/Desktop/myNavierStokesSimplerSolver/postProcessing/mesh/mesh.csv","w");  //opening the mesh file

	fprintf(fp,"x Coordinate,y Coordinate,z Coordinate\n");

	for(k=0; k<=nz; k++)
	 for(j=0; j<=ny; j++)
	   for(i=0; i<=nx; i++)
	{
	  xCoord[i][j][k] = dx * i;
	  yCoord[i][j][k] = dy * j;
	  zCoord[i][j][k] = dz * k;

	  fprintf(fp,"%e,%e,%e\n",xCoord[i][j][k],yCoord[i][j][k],zCoord[i][j][k]);
	}

	fclose(fp);

	printf("Creating mesh\n");

}

//Outputs values to files
void writeOutput()
{


	//output file names
	if(t<100)
	{
		f0 = (char)('0'+ (t % 10));
		f1 = (char)('0'+ (t / 10));
		f2 = '0';
		f3 = '0';
		f4 = '0';
		f5 = '0';
		f6 = '0';
		f7 = '0';
	}

	if((t>=100)&&(t<1000))
	{
		f0 = (char)('0'+ (t % 10));
		f1 = (char)('0'+ (t / 10 % 10));
		f2 = (char)('0'+ (t / 100));
		f3 = '0';
		f4 = '0';
		f5 = '0';
		f6 = '0';
		f7 = '0';
	}

	if((t>=1000)&&(t<10000))
	{
		f0 = (char)('0'+ (t % 10));
		f1 = (char)('0'+ (t % 100 / 10));
		f2 = (char)('0'+ (t % 1000 / 100));
		f3 = (char)('0'+ (t / 1000));
		f4 = '0';
		f5 = '0';
		f6 = '0';
		f7 = '0';
	}

	if((t>=10000)&&(t<100000))
	{
		f0 = (char)('0'+ (t % 10));
		f1 = (char)('0'+ (t % 100 / 10));
		f2 = (char)('0'+ (t % 1000 / 100));
		f3 = (char)('0'+ (t % 10000 / 1000));
		f4 = (char)('0'+ (t / 10000));
		f5 = '0';
		f6 = '0';
		f7 = '0';
	}

	if((t>=100000)&&(t<1000000))
	{
		f0 = (char)('0'+ (t % 10));
		f1 = (char)('0'+ (t % 100 / 10));
		f2 = (char)('0'+ (t % 1000 / 100));
		f3 = (char)('0'+ (t % 10000 / 1000));
		f4 = (char)('0'+ (t % 100000 / 100000));
		f5 = (char)('0'+ (t / 100000));
		f6 = '0';
		f7 = '0';
	}

	if((t>=1000000)&&(t<10000000))
	{
		f0 = (char)('0'+ (t % 10));
		f1 = (char)('0'+ (t % 100 / 10));
		f2 = (char)('0'+ (t % 1000 / 100));
		f3 = (char)('0'+ (t % 10000 / 1000));
		f4 = (char)('0'+ (t % 100000 / 10000));
		f5 = (char)('0'+ (t % 1000000 / 100000));
		f6 = (char)('0'+ (t / 1000000));
		f7 = '0';
	}

	if((t>=10000000)&&(t<100000000))
	{
		f0 = (char)('0'+ (t % 10));
		f1 = (char)('0'+ (t % 100 / 10));
		f2 = (char)('0'+ (t % 1000 / 100));
		f3 = (char)('0'+ (t % 10000 / 1000));
		f4 = (char)('0'+ (t % 100000 / 10000));
		f5 = (char)('0'+ (t % 1000000 / 100000));
		f6 = (char)('0'+ (t % 10000000 / 1000000));
		f7 = (char)('0'+ (t / 10000000));
	}



	strcpy (str,"output/output");
	str[13]=f7;
	str[14]=f6;
	str[15]=f5;
	str[16]=f4;
	str[17]=f3;
	str[18]=f2;
	str[19]=f1;
	str[20]='.';
	str[21]='c';
	str[22]='s';
	str[23]='v';


	fp=fopen(str,"w");  //opening the output file

	fprintf(fp,"x Coordinate,y Coordinate,z Coordinate,u Velocity,v Velocity,w Velocity,pressure\n");

	for  (k=1; k<=nz+2; k++)
	 for (j=1; j<=ny+2; j++)
	  for(i=1; i<=nx+2; i++)
	{

		xCoord[i][j][k] = dx * i ;
	  	yCoord[i][j][k] = dy * j ;
	  	zCoord[i][j][k] = dz * k ;


		if((i==1)&&(j==1)&&(k==1))
		fprintf(fp,"%e,%e,%e,%e,%e,%e,%e\n",xCoord[i][j][k]+dx,yCoord[i][j][k]+dy,zCoord[i][j][k]+dz,0.0,0.0,0.0,p[i][j][k][l]);

		else if((i==(nx+2))&&(j==(ny+2))&&(k==(nz+2)))
		fprintf(fp,"%e,%e,%e,%e,%e,%e,%e\n",xCoord[i][j][k],yCoord[i][j][k],zCoord[i][j][k],0.0,0.0,0.0,p[i][j][k][l]);

		else if((i==1)&&(j==(ny+2))&&(k==(nz+2)))
		fprintf(fp,"%e,%e,%e,%e,%e,%e,%e\n",xCoord[i][j][k]+dx,yCoord[i][j][k],zCoord[i][j][k],0.0,0.0,0.0,p[i][j][k][l]);

		else if((i==(nx+2))&&(j==1)&&(k==(nz+2)))
		fprintf(fp,"%e,%e,%e,%e,%e,%e,%e\n",xCoord[i][j][k],yCoord[i][j][k]+dy,zCoord[i][j][k],0.0,0.0,0.0,p[i][j][k][l]);

		else if((i==(nx+2))&&(j==(ny+2))&&(k==1))
		fprintf(fp,"%e,%e,%e,%e,%e,%e,%e\n",xCoord[i][j][k],yCoord[i][j][k],zCoord[i][j][k]+dz,0.0,0.0,0.0,p[i][j][k][l]);

		else if((i==1)&&(j==1)&&(k==(nz+2)))
		fprintf(fp,"%e,%e,%e,%e,%e,%e,%e\n",xCoord[i][j][k]+dx,yCoord[i][j][k]+dy,zCoord[i][j][k],0.0,0.0,0.0,p[i][j][k][l]);

		else if((i==(nx+2))&&(j==1)&&(k==1))
		fprintf(fp,"%e,%e,%e,%e,%e,%e,%e\n",xCoord[i][j][k],yCoord[i][j][k]+dy,zCoord[i][j][k]+dz,0.0,0.0,0.0,p[i][j][k][l]);

		else if((i==1)&&(j==(ny+2))&&(k==1))
		fprintf(fp,"%e,%e,%e,%e,%e,%e,%e\n",xCoord[i][j][k]+dx,yCoord[i][j][k],zCoord[i][j][k]+dz,0.0,0.0,0.0,p[i][j][k][l]);

		else if((i==1)&&(j==1))
		fprintf(fp,"%e,%e,%e,%e,%e,%e,%e\n",xCoord[i][j][k]+dx,yCoord[i][j][k]+dy,zCoord[i][j][k]+dz/2,0.0,0.0,0.0,p[i][j][k][l]);

		else if((k==1)&&(j==1))
		fprintf(fp,"%e,%e,%e,%e,%e,%e,%e\n",xCoord[i][j][k]+dx/2,yCoord[i][j][k]+dy,zCoord[i][j][k]+dz,0.0,0.0,0.0,p[i][j][k][l]);

		else if((i==1)&&(k==1))
		fprintf(fp,"%e,%e,%e,%e,%e,%e,%e\n",xCoord[i][j][k]+dx,yCoord[i][j][k]+dy/2,zCoord[i][j][k]+dz,0.0,0.0,0.0,p[i][j][k][l]);

		else if((i==(nx+2))&&(j==(ny+2)))
		fprintf(fp,"%e,%e,%e,%e,%e,%e,%e\n",xCoord[i][j][k],yCoord[i][j][k],zCoord[i][j][k]+dz/2,0.0,0.0,0.0,p[i][j][k][l]);

		else if((k==(nz+2))&&(j==(ny+2)))
		fprintf(fp,"%e,%e,%e,%e,%e,%e,%e\n",xCoord[i][j][k]+dx/2,yCoord[i][j][k],zCoord[i][j][k],0.0,0.0,0.0,p[i][j][k][l]);

		else if((i==(nx+2))&&(k==(nz+2)))
	        fprintf(fp,"%e,%e,%e,%e,%e,%e,%e\n",xCoord[i][j][k],yCoord[i][j][k]+dy/2,zCoord[i][j][k],0.0,0.0,0.0,p[i][j][k][l]);

		else if((i==1)&&(j==(ny+2)))
		fprintf(fp,"%e,%e,%e,%e,%e,%e,%e\n",xCoord[i][j][k]+dx,yCoord[i][j][k],zCoord[i][j][k]+dz/2,0.0,0.0,0.0,p[i][j][k][l]);

		else if((k==1)&&(j==(ny+2)))
		fprintf(fp,"%e,%e,%e,%e,%e,%e,%e\n",xCoord[i][j][k]+dx/2,yCoord[i][j][k],zCoord[i][j][k]+dz,0.0,0.0,0.0,p[i][j][k][l]);

		else if((i==(nx+2))&&(k==1))
		fprintf(fp,"%e,%e,%e,%e,%e,%e,%e\n",xCoord[i][j][k],yCoord[i][j][k]+dy/2,zCoord[i][j][k]+dz,0.0,0.0,0.0,p[i][j][k][l]);

		else if((i==(nx+2))&&(j==1))
		fprintf(fp,"%e,%e,%e,%e,%e,%e,%e\n",xCoord[i][j][k],yCoord[i][j][k]+dy,zCoord[i][j][k]+dz/2,0.0,0.0,0.0,p[i][j][k][l]);

		else if((k==(nz+2))&&(j==1))
		fprintf(fp,"%e,%e,%e,%e,%e,%e,%e\n",xCoord[i][j][k]+dx/2,yCoord[i][j][k]+dy,zCoord[i][j][k],0.0,0.0,0.0,p[i][j][k][l]);

		else if((i==1)&&(k==(nz+2)))
	        fprintf(fp,"%e,%e,%e,%e,%e,%e,%e\n",xCoord[i][j][k]+dx,yCoord[i][j][k]+dy/2,zCoord[i][j][k],0.0,0.0,0.0,p[i][j][k][l]);

		else if(i==1)
		fprintf(fp,"%e,%e,%e,%e,%e,%e,%e\n",xCoord[i][j][k]+dx,yCoord[i][j][k]+dy/2,zCoord[i][j][k]+dz/2,u[i][j][k][l],(v[i][j][k][l]+v[i+1][j][k][l])/2,(w[i][j][k][l]+w[i+1][j][k][l])/2,p[i][j][k][l]);

		else if(j==1)
		fprintf(fp,"%e,%e,%e,%e,%e,%e,%e\n",xCoord[i][j][k]+dx/2,yCoord[i][j][k]+dy,zCoord[i][j][k]+dz/2,(u[i][j][k][l]+u[i][j+1][k][l])/2,v[i][j][k][l],(w[i][j][k][l]+w[i][j+1][k][l])/2,p[i][j][k][l]);

		else if(k==1)
		fprintf(fp,"%e,%e,%e,%e,%e,%e,%e\n",xCoord[i][j][k]+dx/2,yCoord[i][j][k]+dy/2,zCoord[i][j][k]+dz,(u[i][j][k][l]+u[i][j][k+1][l])/2,(v[i][j][k][l]+v[i][j][k+1][l])/2,w[i][j][k][l],p[i][j][k][l]);

		else if(i==(nx+2))
		fprintf(fp,"%e,%e,%e,%e,%e,%e,%e\n",xCoord[i][j][k],yCoord[i][j][k]+dy/2,zCoord[i][j][k]+dz/2,u[i-1][j][k][l],(v[i][j][k][l]+v[i-1][j][k][l])/2,(w[i][j][k][l]+w[i-1][j][k][l])/2,p[i][j][k][l]);

		else if(j==(ny+2))
		fprintf(fp,"%e,%e,%e,%e,%e,%e,%e\n",xCoord[i][j][k]+dx/2,yCoord[i][j][k],zCoord[i][j][k]+dz/2,(u[i][j][k][l]+u[i][j-1][k][l])/2,v[i][j-1][k][l],(w[i][j][k][l]+w[i][j-1][k][l])/2,p[i][j][k][l]);

		else if(k==(nz+2))
		fprintf(fp,"%e,%e,%e,%e,%e,%e,%e\n",xCoord[i][j][k]+dx/2,yCoord[i][j][k]+dy/2,zCoord[i][j][k],(u[i][j][k][l]+u[i][j][k-1][l])/2,(v[i][j][k][l]+v[i][j][k-1][l])/2,w[i][j][k-1][l],p[i][j][k][l]);

		else
	     {
		xCoord[i][j][k] = dx * i + dx/2;
	  	yCoord[i][j][k] = dy * j + dy/2;
	  	zCoord[i][j][k] = dz * k + dz/2;

		fprintf(fp,"%e,%e,%e,%e,%e,%e,%e\n",xCoord[i][j][k],yCoord[i][j][k],zCoord[i][j][k],finalU[i][j][k][l],finalV[i][j][k][l],finalW[i][j][k][l],p[i][j][k][l]);
	     }
	}


	fclose(fp);

}


