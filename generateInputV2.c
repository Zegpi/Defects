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
	PetscPrintf(PETSC_COMM_WORLD,"Start of generateInput \n");

	PetscInt commsize,rank;
	ierr = MPI_Comm_size(PETSC_COMM_WORLD,&commsize);CHKERRQ(ierr);
	ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank);CHKERRQ(ierr);
//

time_t T=time(NULL);
struct tm tm=*localtime(&T);

PetscPrintf(PETSC_COMM_WORLD,"Current time is %02d:%02d:%02d \n",tm.tm_hour,tm.tm_min,tm.tm_sec);

//Generate Mesh and copy cpp to result folder to save for reproduction (check that parameters are the same on file to run)

	PetscInt b=200;					//Parmeter to choose size of cores, must always be odd, core will be of size 1 unit, rest of the body will be of size b-1 units in each direction
	PetscReal Lx=20.0;
	PetscReal Ly=20.0;
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
	char direct2[]="../../sharedResults";
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

	ierr=system("cp ~/shared/CodigosPetIGA/generateInput.c ~/sharedResults/");
//

//General parameters
	PetscInt N,M, numEls; 	
	N=2;	M=1;	numEls=1;

	//Para el caso de ahora, usar (112,6,6,2)
//

//Creation of Initialization of Pi
	PetscPrintf(PETSC_COMM_WORLD,"System for initial state of Pi starting \n");
	IGA igaPi;
	ierr = IGACreate(PETSC_COMM_WORLD,&igaPi);CHKERRQ(ierr);
	ierr = IGASetDim(igaPi,2);CHKERRQ(ierr);														//Spatial dimension of the problem
	ierr = IGASetDof(igaPi,4);CHKERRQ(ierr);														//Number of degrees of freedom, per node
	ierr = IGASetOrder(igaPi,1);CHKERRQ(ierr);														//Number of spatial derivatives to calculate
	ierr = IGASetFromOptions(igaPi);CHKERRQ(ierr);													//Note: The order (or degree) of the shape functions is given by the mesh!
	ierr = IGARead(igaPi,"./geometry.dat");CHKERRQ(ierr);
	ierr = IGASetUp(igaPi);CHKERRQ(ierr);

	Vec pi0;
	ierr = IGACreateVec(igaPi,&pi0);CHKERRQ(ierr);  

	PetscInt *points;
	PetscReal *valores;

	points=(PetscInt*)calloc(8192,sizeof(PetscInt));
	valores=(PetscReal*)calloc(8192,sizeof(PetscReal));

	PetscReal c,t;	c=Lx/nx;	t=Ly/ny;
	PetscInt center;

	PetscReal a=3.585284301;
	PetscPrintf(PETSC_COMM_WORLD,"a= %f \n",a);
	PetscReal beta=-0.100211;
	PetscReal gamma=-beta*a/(2.0*a*a)/(c*t);
	PetscPrintf(PETSC_COMM_WORLD,"gamma= %f \n",gamma*c*t);
	PetscReal tanTh=tan(45.0/180.0*ConstPi)/(c*t);

	PetscInt counter=0;

	if(b%2==0)		
	{
		PetscInt realNumEls=pow(2,numEls-1);
		center=(nx+1)*ny/2+nx/2;
		if(M<2)
		{
			//M=2;
		}
		if(N<2)
		{
			//N=2;
		}
		for (int i=0; i<pow(2,numEls-1)+1; i++)
		{
			for (int j=0; j<pow(2,numEls-1)+1; j++)
			{
				//points[counter]  =4*center+(N-1)*4*(nx+1)+(M-1)*4-(numEls-i-1)*4*(nx+1)-(numEls-j-1)*4+1;		valores[counter]=(-tanTh)/((realNumEls+1)*(realNumEls+1));
				//points[counter+1]=4*center+(N-1)*4*(nx+1)+(M-1)*4-(numEls-i-1)*4*(nx+1)-(numEls-j-1)*4+2;		valores[counter+1]=(tanTh-gamma)/((realNumEls+1)*(realNumEls+1));
				//points[counter+2]=4*center+(N-1)*4*(nx+1)+(M-1)*4-(numEls-i-1)*4*(nx+1)-(numEls-j-1)*4+3;		valores[counter+2]=(-gamma)/((realNumEls+1)*(realNumEls+1));

				points[counter]  =4*center+N*4*(nx+1)+M*4-(i)*4*(nx+1)-(j)*4+1;								valores[counter]=(-tanTh)/((realNumEls+1)*(realNumEls+1));
				points[counter+1]=4*center+N*4*(nx+1)+M*4-(i)*4*(nx+1)-(j)*4+2;								valores[counter+1]=(tanTh-gamma)/((realNumEls+1)*(realNumEls+1));
				points[counter+2]=4*center+N*4*(nx+1)+M*4-(i)*4*(nx+1)-(j)*4+3;								valores[counter+2]=(-gamma)/((realNumEls+1)*(realNumEls+1));
				
				//points[counter+3]=4*center-N*4*(nx+1)-M*4-(numEls-i-1)*4*(nx+1)-(numEls-j-1)*4+1;				valores[counter+3]=(tanTh)/((realNumEls+1)*(realNumEls+1));
				//points[counter+4]=4*center-N*4*(nx+1)-M*4-(numEls-i-1)*4*(nx+1)-(numEls-j-1)*4+2;				valores[counter+4]=((-tanTh+gamma))/((realNumEls+1)*(realNumEls+1));
				//points[counter+5]=4*center-N*4*(nx+1)-M*4-(numEls-i-1)*4*(nx+1)-(numEls-j-1)*4+3;				valores[counter+5]=(gamma)/((realNumEls+1)*(realNumEls+1));

				points[counter+3]=4*center-N*4*(nx+1)-M*4-(i)*4*(nx+1)-(j)*4+1;								valores[counter+3]=(tanTh)/((realNumEls+1)*(realNumEls+1));
				points[counter+4]=4*center-N*4*(nx+1)-M*4-(i)*4*(nx+1)-(j)*4+2;								valores[counter+4]=((-tanTh+gamma))/((realNumEls+1)*(realNumEls+1));
				points[counter+5]=4*center-N*4*(nx+1)-M*4-(i)*4*(nx+1)-(j)*4+3;								valores[counter+5]=(gamma)/((realNumEls+1)*(realNumEls+1));

				counter=counter+6;
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
				points[counter  ]=4*center+N*4*(nx+1)+M*4+1-(numEls-i-1)*4*(nx+1)-(numEls-j-1)*4;			valores[counter  ]=(-tanTh)/((8*numEls-4)*0.5+4*0.25+(2*numEls-1)*(2*numEls-1));
				points[counter+1]=4*center+N*4*(nx+1)+M*4+2-(numEls-i-1)*4*(nx+1)-(numEls-j-1)*4;			valores[counter+1]=(tanTh-gamma)/((8*numEls-4)*0.5+4*0.25+(2*numEls-1)*(2*numEls-1));
				points[counter+2]=4*center+N*4*(nx+1)+M*4+3-(numEls-i-1)*4*(nx+1)-(numEls-j-1)*4;			valores[counter+2]=(-gamma)/((8*numEls-4)*0.5+4*0.25+(2*numEls-1)*(2*numEls-1));
				
				points[counter+3]=4*center-N*4*(nx+1)-M*4+1-(numEls-i-1)*4*(nx+1)-(numEls-j-1)*4;			valores[counter+3]=(tanTh)/((8*numEls-4)*0.5+4*0.25+(2*numEls-1)*(2*numEls-1));
				points[counter+4]=4*center-N*4*(nx+1)-M*4+2-(numEls-i-1)*4*(nx+1)-(numEls-j-1)*4;			valores[counter+4]=(-tanTh+gamma)/((8*numEls-4)*0.5+4*0.25+(2*numEls-1)*(2*numEls-1));
				points[counter+5]=4*center-N*4*(nx+1)-M*4+3-(numEls-i-1)*4*(nx+1)-(numEls-j-1)*4;			valores[counter+5]=(gamma)/((8*numEls-4)*0.5+4*0.25+(2*numEls-1)*(2*numEls-1));

				counter=counter+6;
			}
		}
	}
	
	ierr = VecAssemblyBegin(pi0);CHKERRQ(ierr);
	ierr = VecAssemblyEnd  (pi0);CHKERRQ(ierr);

	ierr = VecSetValues(pi0,8192,points,valores,ADD_VALUES);

	ierr = VecAssemblyBegin(pi0);CHKERRQ(ierr);
	ierr = VecAssemblyEnd  (pi0);CHKERRQ(ierr);

	free(points);
	free(valores);

	char namePi[]="/Input-Pi-2d-0.dat";
	char pathPi[512];
	sprintf(pathPi,"%s%s",direct,namePi);
	ierr = IGAWriteVec(igaPi,pi0,pathPi);CHKERRQ(ierr);
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

	PetscReal e1=-0.0*a/(c*t);
	PetscReal e2=0.0*(a-0.100211)/(c*t);
	//PetscReal e2=(2.0)/(c*t);
	//PetscReal lx=2.0*N*Lx/nx;

	PetscInt *pointsAl;
	PetscReal *valoresAl;

	pointsAl=(PetscInt*)calloc(8192,sizeof(PetscInt));
	valoresAl=(PetscReal*)calloc(8192,sizeof(PetscReal));

	ierr = VecAssemblyBegin(alp0);CHKERRQ(ierr);
	ierr = VecAssemblyEnd  (alp0);CHKERRQ(ierr);

	PetscPrintf(PETSC_COMM_WORLD,"Commsize es %d \n",commsize);

	//Borrar esto despues
	//N=0; M=0;
	//N=N-1; M=M-1;
	//N=-N+1; M=-M+1;
	//hasta aqui

	counter=0;

	if(b%2==0)		
	{
		PetscInt realNumEls=pow(2,numEls-1);
		center=(nx+1)*ny/2+nx/2;
		if(M<2)
		{
			//M=2;
		}
		if(N<2)
		{
			//N=2;
		}
		for (int i=0; i<pow(2,numEls-1)+1; i++)
		{
			for (int j=0; j<pow(2,numEls-1)+1; j++)
			{
				//pointsAl[counter]=  2*center-N*2*(nx+1)-M*2+1-(numEls-i-1)*2*(nx+1)-(numEls-j-1)*2-1;				valoresAl[counter]=  e1/((realNumEls+1.0)*(realNumEls+1.0));
				//pointsAl[counter+1]=2*center-N*2*(nx+1)-M*2+1-(numEls-i-1)*2*(nx+1)-(numEls-j-1)*2;							valoresAl[counter+1]=e2/((realNumEls+1.0)*(realNumEls+1.0));
									
				pointsAl[counter]  =2*center+(N+1)*2*(nx+1)+(M+1)*2-(i)*2*(nx+1)-(j)*2;									valoresAl[counter]  =e1/((realNumEls+1.0)*(realNumEls+1.0));
				pointsAl[counter+1]=2*center+(N+1)*2*(nx+1)+(M+1)*2-(i)*2*(nx+1)-(j)*2+1;								valoresAl[counter+1]=e2/((realNumEls+1.0)*(realNumEls+1.0));

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


	/*
	for (int i=0; i<commsize; i++)
	{
		if (rank==i)
		{
			pointsAl[0]=2*center+2*(nx+1)*N+2*M;						valoresAl[0]=coef/4.0;
			pointsAl[1]=2*center+2*(nx+1)*N+2*M+2;						valoresAl[1]=coef/4.0;
			pointsAl[2]=2*center+2*(nx+1)*(N+1)+2*M;					valoresAl[2]=coef/4.0;
			pointsAl[3]=2*center+2*(nx+1)*(N+1)+2*M+2;					valoresAl[3]=coef/4.0;

			pointsAl[4]=2*center+2*(nx+1)*N+2*M+1;						valoresAl[4]=(coef+e)/4.0;
			pointsAl[5]=2*center+2*(nx+1)*N+2*M+2+1;					valoresAl[5]=(coef+e)/4.0;
			pointsAl[6]=2*center+2*(nx+1)*(N+1)+2*M+1;					valoresAl[6]=(coef+e)/4.0;
			pointsAl[7]=2*center+2*(nx+1)*(N+1)+2*M+2+1;				valoresAl[7]=(coef+e)/4.0;

			for (int j=0; j<8; j++)
			{
				PetscPrintf(PETSC_COMM_SELF," rank = %d points[%d]=%d \n ", rank, j, pointsAl[j]);
				PetscPrintf(PETSC_COMM_SELF," rank = %d valores[%d]=%f \n ", rank, j, valoresAl[j]);
			}	
		}
	}
	*/

	//ierr = VecSetValues(pi0,24+(Ndisl-0)*48,points,valores,INSERT_VALUES);
	ierr = VecSetValues(alp0,8192,pointsAl,valoresAl,ADD_VALUES);

	ierr = VecAssemblyBegin(alp0);CHKERRQ(ierr);
	ierr = VecAssemblyEnd  (alp0);CHKERRQ(ierr);

	
	char nameAlp[]="/Input-Al-2d-0.dat";
	char pathAlp[512];
	sprintf(pathAlp,"%s%s",direct,nameAlp);
	ierr = IGAWriteVec(igaAl,alp0,pathAlp);CHKERRQ(ierr);
//

//Destroy all objects not needed anymore (Better to do it here in case different codes call the same IGA, move if memory is a problem)
	ierr = IGADestroy(&igaAl);CHKERRQ(ierr);
	ierr = IGADestroy(&igaPi);CHKERRQ(ierr);
//

ierr = PetscFinalize();CHKERRQ(ierr);

return 0;
}