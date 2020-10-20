#include "petiga.h"
#include <stdio.h>
#include <dirent.h>
#include <math.h>
#include <time.h>

#if PETSC_VERSION_LT(3,5,0)
#define KSPSetOperators(ksp,A,B) KSPSetOperators(ksp,A,B,SAME_NONZERO_PATTERN)
#endif

#define ConstPi 3.14159265358979323846

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc, char *argv[]) {
//Creation of variables
	PetscErrorCode  ierr;
	ierr = PetscInitialize(&argc,&argv,0,0);CHKERRQ(ierr);										//Always initialize PETSc
	PetscPrintf(PETSC_COMM_WORLD,"Start of generateInputS \n");

	PetscInt commsize,rank;
	ierr = MPI_Comm_size(PETSC_COMM_WORLD,&commsize);CHKERRQ(ierr);
	ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank);CHKERRQ(ierr);
//

time_t T=time(NULL);
struct tm tm=*localtime(&T);

PetscPrintf(PETSC_COMM_WORLD,"Current time is %02d:%02d:%02d \n",tm.tm_hour,tm.tm_min,tm.tm_sec);

//Generate Mesh and copy cpp to result folder to save for reproduction (check that parameters are the same on file to run)

	PetscInt b=581;					//Parmeter to choose size of cores, must always be odd, core will be of size 1 unit, rest of the body will be of size b-1 units in each direction
	PetscReal Lx=80.0;
	PetscReal Ly=80.0;
	PetscInt  nx=b;
	PetscInt  ny=b;

	PetscReal h=fmin(Lx/nx,Ly/ny);

	PetscPrintf(PETSC_COMM_WORLD,"h= %f\n",h);

	char      nombreMalla[12]="geometry";
	char 	  comando[512];
	sprintf(comando,"exec python rectangle.py %s %f %f %d %d\n",nombreMalla,Lx,Ly,nx,ny);
	ierr=system(comando);

	strcpy(nombreMalla,"geometry2");
	sprintf(comando,"exec python rectangle2.py %s %f %f %d %d\n",nombreMalla,Lx,Ly,nx,ny);
	ierr=system(comando);

	strcpy(nombreMalla,"geometry3");
	sprintf(comando,"exec python rectangle3.py %s %f %f %d %d\n",nombreMalla,Lx,Ly,nx,ny);
	ierr=system(comando);

//

//Delete files in result folder
	char direct1[]="./results";
	char direct2[]="../Results";
	char *direct;

	DIR *pDir1=opendir(direct1);

	if (pDir1 != NULL)
	{
		closedir(pDir1);
		direct=direct1;
	}

	DIR *pDir2=opendir(direct2);

	if (pDir2 != NULL)
	{
		closedir(pDir2);
		direct=direct2;
	}

	DIR *folder = opendir(direct);
	struct dirent *next_file;
	char filepath[512];
	while ( (next_file = readdir(folder)) != NULL )
	{
		sprintf(filepath, "%s/%s", direct, next_file->d_name);
		remove(filepath);
	}
	closedir(folder);

	FILE *source, *dest;
	char buffer[8192];
	size_t bytes;
	source = fopen("./generateInputS.c","r");
	dest   = fopen("../Results/generateInputS.c","w");

	while (0 < (bytes = fread(buffer, 1, sizeof(buffer), source)))
    	fwrite(buffer, 1, bytes, dest);

	fclose(source);
	fclose(dest);
//

//General parameters
	PetscInt N,M, numEls; 	
	N=0;	M=0;	numEls=1;
//

//Creation of Initialization of S
	PetscPrintf(PETSC_COMM_WORLD,"System for initial state of S started \n");
	IGA igaS;
	ierr = IGACreate(PETSC_COMM_WORLD,&igaS);CHKERRQ(ierr);
	ierr = IGASetDim(igaS,2);CHKERRQ(ierr);														//Spatial dimension of the problem
	ierr = IGASetDof(igaS,8);CHKERRQ(ierr);														//Number of degrees of freedom, per node
	ierr = IGASetOrder(igaS,1);CHKERRQ(ierr);													//Number of spatial derivatives to calculate
	ierr = IGASetFromOptions(igaS);CHKERRQ(ierr);												//Note: The order (or degree) of the shape functions is given by the mesh!
	ierr = IGARead(igaS,"./geometry.dat");CHKERRQ(ierr);
	ierr = IGASetUp(igaS);CHKERRQ(ierr);

	Vec s0;
	ierr = IGACreateVec(igaS,&s0);CHKERRQ(ierr);

	PetscInt nf,nc,numF,cord,counter, *pointsS, dist;
	PetscReal *valoresS,c,t,l,g[8];

	//nf=(b-1)/2; 		//For odd values of b 
	nf=(b-2)/2;			//For even values of b 
	//nc=(b+1)/2;			//Half the length of the body, for a disclination on the center
	nc=b+1;				//Full length of the body, for something like a through twin
	numF=8;			//Number of node rows to assign, rows of elements will be one less
	dist=80;				//Distance from center to eigenwall (Distance center to center is 2*dist [elements])
	
	ierr = PetscCalloc1(1048576,&pointsS);CHKERRQ(ierr);
	ierr = PetscCalloc1(1048576,&valoresS);CHKERRQ(ierr);

	ierr = VecAssemblyBegin(s0);CHKERRQ(ierr);
	ierr = VecAssemblyEnd  (s0);CHKERRQ(ierr);

	c=Lx/nx;
	t=Ly/ny;
	l=t*(numF-1);

	//fullS[0][0][0]=S[0]; fullS[0][0][1]=S[1]; fullS[0][1][0]=S[2]; fullS[0][1][1]=S[3];		//Expand S to full 3rd order form, only non-zero elements
	//fullS[1][0][0]=S[4]; fullS[1][0][1]=S[5]; fullS[1][1][0]=S[6]; fullS[1][1][1]=S[7];

	/*
	PetscReal alpha=0.5;																	//Alpha linearly interpolates between a pure twin (alpha=0.0) and a pure rotation grain boundary (alpha=1.0)
	
	g[0]=0.0;
	g[1]=alpha*(cos(85.0/180.0*ConstPi)-cos(90.0/180.0*ConstPi));								//S_112
	g[2]=0.0;
	g[3]=alpha*(-sin(85.0/180.0*ConstPi)+sin(90.0/180.0*ConstPi))+(1.0-alpha)*1.0;				//S_122
	g[4]=0.0;
	g[5]=alpha*(sin(85.0/180.0*ConstPi)-sin(90.0/180.0*ConstPi));								//S_212
	g[6]=0.0;
	g[7]=alpha*(cos(85.0/180.0*ConstPi)-cos(90.0/180.0*ConstPi));								//S_222

	//nf=nf-(numF-2)/2;		//For odd values of b
	nf=nf-(numF-3)/2;		//For even values of b
	counter=0;
	for (int i=nf; i<nf+numF; i++)
	{
		cord=(nx+1)*i*8;
		for (int j=0; j<nc; j++)
		{
			for (int k=0; k<8; k++)
			{
				pointsS[counter]=cord-8*(nx+1)*dist;
				//pointsS[counter+1]=cord-8*(nx+1)*dist;
				valoresS[counter]=g[k]/l;
				//valoresS[counter+1]=g[k]/l;
				cord=cord+1;
				counter=counter+1;
			}
		}
	}
	*/

	//Uncomment for 2 eigenwalls
	/*
	for (int i=nf; i<nf+numF; i++)
	{
		cord=(nx+1)*i*8;
		for (int j=0; j<nc; j++)
		{
			for (int k=0; k<8; k++)
			{
				//pointsS[counter]=cord+8*(nx+1)*dist;
				pointsS[counter]=cord-8*(nx+1)*dist;
				//valoresS[counter]=g[k];
				valoresS[counter]=g[k];
				cord=cord+1;
				counter=counter+1;
			}
		}
	}
	*/

	//This is for coherency defect (see Amit's email on 11/08/2020, subject: "some thoughts")
	nf=(b-2)/2;			//For even values of b 
	//nc=(b+1)/2;		//Half the length of the body, for a disclination on the center
	nc=b+1;				//Full length of the body, for something like a through twin
	numF=48;				//Number of node rows to assign, rows of elements will be one less
	dist=22;				//Distance from center to eigenwall (Distance center to center is 2*dist [elements])
	counter=0;
	l=t*(numF+1);

	/*
	g[0]=0.0;								//S_111 (these indices in 1 based numbering)
	g[1]=0.0;								//S_112
	g[2]=0.0;								//S_121
	g[3]=1.0;								//S_122
	g[4]=0.0;								//S_211
	g[5]=0.0;								//S_212
	g[6]=0.0;								//S_221
	g[7]=0.0;								//S_222
	*/

	g[0]=0.0;
	g[1]=0.0;																					//S_112
	g[2]=1.0;
	g[3]=0.0;																					//S_122
	g[4]=-1.0;
	g[5]=0.0;																					//S_212
	g[6]=0.0;
	g[7]=0.0;																					//S_222

	PetscReal factor=1.0/(64.0*c*t*numF/8.0);

	for (int i=nf; i<nf+numF; i++)
	{
		cord=(nx+1)*i*8;
		for (int j=0; j<nc; j++)
		{
			for (int k=0; k<8; k++)
			{ 
				if(j>=295)
				//if(j>=154)
				{
					pointsS[counter]=cord-8*(nx+1)*dist;
					valoresS[counter]=0.0;
					cord=cord+1;
					counter=counter+1;	
				}
				else if(j>=287)
				//else if(j>=147)
				{
					pointsS[counter]=cord-8*(nx+1)*dist;
					valoresS[counter]=factor*g[k];
					//valoresS[counter]=factor*(g[k]/l)* (double)(j-147.0)/(7.0);
					cord=cord+1;
					counter=counter+1;
				}
				else
				{
					pointsS[counter]=cord-8*(nx+1)*dist;
					valoresS[counter]=0.0;
					cord=cord+1;
					counter=counter+1;
				}
			}
		}
	}

	ierr = VecSetValues(s0,8*nc*numF*4,pointsS,valoresS,ADD_VALUES);	
	ierr = VecAssemblyBegin(s0);CHKERRQ(ierr);
	ierr = VecAssemblyEnd  (s0);CHKERRQ(ierr);

	char nameS[]="/Input-S-2d-0.dat";
	char pathS[512];
	sprintf(pathS,"%s%s",direct,nameS);
	ierr = IGAWriteVec(igaS,s0,pathS);CHKERRQ(ierr);
//

//Creation of Initialization of Alfa0
	PetscPrintf(PETSC_COMM_WORLD,"System for initial state of Alpha starting \n");
	IGA igaAl;
	ierr = IGACreate(PETSC_COMM_WORLD,&igaAl);CHKERRQ(ierr);
	ierr = IGASetDim(igaAl,2);CHKERRQ(ierr);														//Spatial dimension of the problem
	ierr = IGASetDof(igaAl,2);CHKERRQ(ierr);														//Number of degrees of freedom, per node
	ierr = IGASetOrder(igaAl,1);CHKERRQ(ierr);														//Number of spatial derivatives to calculate
	ierr = IGASetFromOptions(igaAl);CHKERRQ(ierr);													//Note: The order (or degree) of the shape functions is given by the mesh!
	ierr = IGARead(igaAl,"./geometry.dat");CHKERRQ(ierr);
	ierr = IGASetUp(igaAl);CHKERRQ(ierr);
	
	Vec alp0;
	ierr = IGACreateVec(igaAl,&alp0);CHKERRQ(ierr);

	PetscReal a=0.0;	c=Lx/nx;	t=Ly/ny;
	PetscInt center;

	PetscReal e1=(0.0+0.0*a)/(c*t);
	PetscReal e2=0.0*(a-0.100211)/(c*t);
	//PetscReal e2=(2.0)/(c*t);
	//PetscReal lx=2.0*N*Lx/nx;

	PetscInt *pointsAl;
	PetscReal *valoresAl;

	ierr = PetscCalloc1(8192,&pointsAl);CHKERRQ(ierr);
	ierr = PetscCalloc1(8192,&valoresAl);CHKERRQ(ierr);

	ierr = VecAssemblyBegin(alp0);CHKERRQ(ierr);
	ierr = VecAssemblyEnd  (alp0);CHKERRQ(ierr);

	//Borrar esto despues
	N=0; M=0;
	//N=N-1; M=M-1;
	//N=-N+1; M=-M+1;
	//hasta aqui

	counter=0;

	if(b%2==0)		
	{
		PetscInt realNumEls=pow(2,numEls-1);
		center=(nx+1)*(ny/2-(numEls+1))+nx/2-(numEls+1);
		for (int i=0; i<pow(2,numEls-1)+1; i++)
		{
			for (int j=0; j<pow(2,numEls-1)+1; j++)
			{
				//pointsAl[counter]=  2*center-N*2*(nx+1)-M*2+1-(numEls-i-1)*2*(nx+1)-(numEls-j-1)*2-1;				valoresAl[counter]=  e1/((realNumEls+1.0)*(realNumEls+1.0));
				//pointsAl[counter+1]=2*center-N*2*(nx+1)-M*2+1-(numEls-i-1)*2*(nx+1)-(numEls-j-1)*2;							valoresAl[counter+1]=e2/((realNumEls+1.0)*(realNumEls+1.0));
									
				pointsAl[counter]  =2*center+(N+1)*2*(nx+1)+(M+1)*2+(i)*2*(nx+1)+(j)*2;									valoresAl[counter]  =e1/((realNumEls+1.0)*(realNumEls+1.0));
				pointsAl[counter+1]=2*center+(N+1)*2*(nx+1)+(M+1)*2+(i)*2*(nx+1)+(j)*2+1;								valoresAl[counter+1]=e2/((realNumEls+1.0)*(realNumEls+1.0));

				counter=counter+2;					//Modify accodringly to hoy many gdl you are setting
			}
		}

	}
	else
	{
		center=(nx+1)*ny/2-1;
		for (int i=0; i<2*numEls; i++)
		{
			for (int j=0; j<2*numEls; j++)
			{
				//pointsAl[counter]=  2*center+N*2*(nx+1)+M*2+1-(numEls-i-1)*2*(nx+1)-(numEls-j-1)*2-1;				valoresAl[counter]=  e1/((8*numEls-4)*0.5+4*0.25+(2*numEls-1)*(2*numEls-1));
				//pointsAl[counter+1]=2*center+N*2*(nx+1)+M*2+1-(numEls-i-1)*2*(nx+1)-(numEls-j-1)*2;					valoresAl[counter+1]=e2/((8*numEls-4)*0.5+4*0.25+(2*numEls-1)*(2*numEls-1));

				//pointsAl[counter+2]=2*center-N*2*(nx+1)-M*2+1-(numEls-i-1)*2*(nx+1)-(numEls-j-1)*2-1;				valoresAl[counter+2]=e/((8*numEls-4)*0.5+4*0.25+(2*numEls-1)*(2*numEls-1));
				//pointsAl[counter+3]=2*center-N*2*(nx+1)-M*2+1-(numEls-i-1)*2*(nx+1)-(numEls-j-1)*2;					valoresAl[counter+3]=e2/((8*numEls-4)*0.5+4*0.25+(2*numEls-1)*(2*numEls-1));

				//points[counter  ]=4*center+N*4*(nx+1)+M*4+1-(numEls-i-1)*4*(nx+1)-(numEls-j-1)*4;
				pointsAl[counter  ]=2*center+N*2*(nx+1)+M*2  -(numEls-i-1)*2*(nx+1)-(numEls-j-1)*2;			valoresAl[counter  ]=e1/((8*numEls-4)*0.5+4*0.25+(2*numEls-1)*(2*numEls-1));
				pointsAl[counter+1]=2*center+N*2*(nx+1)+M*2+1-(numEls-i-1)*2*(nx+1)-(numEls-j-1)*2;			valoresAl[counter+1]=e2/((8*numEls-4)*0.5+4*0.25+(2*numEls-1)*(2*numEls-1));

				counter=counter+2;					//Modify accodringly to hoy many gdl you are setting
			}
		}
	}

	ierr = VecSetValues(alp0,counter+1,pointsAl,valoresAl,ADD_VALUES);

	ierr = VecAssemblyBegin(alp0);CHKERRQ(ierr);
	ierr = VecAssemblyEnd  (alp0);CHKERRQ(ierr);
	
	char nameAlp[]="/Input-Al-2d-0.dat";
	char pathAlp[512];
	sprintf(pathAlp,"%s%s",direct,nameAlp);
	ierr = IGAWriteVec(igaAl,alp0,pathAlp);CHKERRQ(ierr);
//

//Destroy all objects not needed anymore (Better to do it here in case different codes call the same IGA, move if memory is a problem)
	ierr = IGADestroy(&igaAl);CHKERRQ(ierr);
	ierr = IGADestroy(&igaS);CHKERRQ(ierr);
//

ierr = PetscFinalize();CHKERRQ(ierr);

return 0;
}