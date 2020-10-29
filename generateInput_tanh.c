#include "petiga.h"
#include <stdio.h>
#include <dirent.h>
#include <math.h>
#include <time.h>

#if PETSC_VERSION_LT(3,5,0)
#define KSPSetOperators(ksp,A,B) KSPSetOperators(ksp,A,B,SAME_NONZERO_PATTERN)
#endif

#define ConstPi 3.14159265358979323846

//Definition of AppCtxL2
typedef struct 
{
	PetscReal Lx;
	PetscReal Ly;
	PetscReal nx;
	PetscReal ny;
}	AppCtxL2;

//L2 projection of S(0)
	#undef  __FUNCT__
	#define __FUNCT__ "L2ProjectionS"
	PetscErrorCode L2ProjectionS(IGAPoint p,PetscReal *K,PetscReal *F,void *ctx)
	{
		if (p->atboundary)
		{
			return 0;																		//Zero Neumann boundary condition
		}

		AppCtxL2    *user  = (AppCtxL2 *)ctx;
		PetscReal Lx     = user->Lx;
		PetscReal Ly     = user->Ly;
		PetscReal nx     = user->nx;
		PetscReal ny     = user->ny;

		PetscInt a,b,i;
		PetscInt nen = p->nen;																//Number of shape functions
		PetscInt dim = p->dim;																//Spatial dimensions of the problem
		PetscInt dof = p->dof;

		PetscReal x[dim];																	//Vector of reals, size equal to problem's dimension
		IGAPointFormGeomMap(p,x);															//Fills x with the coordinates of p, Gauss's point

		//g is the function to L2 project
		PetscReal dx=Lx/nx;
		PetscReal dy=Ly/ny;

		//S has 8 components, in order S(1,1,1), S(1,1,2), S(1,2,1), S(1,2,2), S(2,1,1), S(2,1,2), S(2,2,1), S(2,2,2)

		PetscReal a1,a2,a3,a4,a5,a6,b0,b1;
		a1=-145.0*dy;
		a2=-138.0*dy;
		a3=  60.0*dx;
		a4=  68.0*dx;
		b0=1.0/5.0;
		b1=1.0/4.5;
		b1=1.0/2.0;

		PetscReal tanh_y=0.5*(tanh((x[1]-a1)/b0)-tanh((x[1]-a2)/b0));
		PetscReal tanh_x=0.5*(tanh((x[0]-a3)/b1)-tanh((x[0]-a4)/b1));
		PetscReal half_tanh_x=0.5*(1.0-tanh((x[0]-a3)/b1));
		PetscReal half_tanh_y=0.5*(tanh((x[1]-a1+5.0*dy)/b1)+1.0);;

		PetscReal g[dof];

		g[0]=0.0;
		g[1]=0.0;
		g[2]= tan(5.0*ConstPi/180.0)*tanh_x*half_tanh_y*0.92;
		g[3]=-tan(5.0*ConstPi/180.0)*tanh_y*half_tanh_x;
		g[4]=-tan(5.0*ConstPi/180.0)*tanh_x*half_tanh_y*0.92;
		g[5]= tan(5.0*ConstPi/180.0)*tanh_y*half_tanh_x;
		g[6]=0.0;
		g[7]=0.0;

		//Consider changing all this to just assigning S directly

		const PetscReal (*N) = (typeof(N)) p->shape[0];
		PetscReal (*FF)[dof] = (PetscReal (*)[dof])F;
		PetscReal (*KK)[dof][nen][dof] = (PetscReal (*)[dof][nen][dof])K;
		for(a=0; a<nen; a++)
		{
			for(i=0; i<dof; i++) 
			{
				for(b=0; b<nen; b++)
				{
		  			KK[a][i][b][i] = N[a]*N[b];
		  		}
		  		FF[a][i] = N[a]*g[i];
			}
		}
		return 0;
	}
//

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc, char *argv[]) {
//Creation of variables
	PetscErrorCode  ierr;
	PetscInt dir,side;
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

	AppCtxL2 userL2;
	userL2.Lx     = Lx;
	userL2.Ly     = Ly;
	userL2.nx     = nx;
	userL2.ny     = ny;

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

//Creation of types and systems for the L2 projection of S0
	PetscPrintf(PETSC_COMM_WORLD,"\n System for L2 projection for S starting \n\n");
	IGA igaS;
	ierr = IGACreate(PETSC_COMM_WORLD,&igaS);CHKERRQ(ierr);
	ierr = IGASetDim(igaS,2);CHKERRQ(ierr);														//Spatial dimension of the problem
	ierr = IGASetDof(igaS,8);CHKERRQ(ierr);														//Number of degrees of freedom, per node
	ierr = IGASetOrder(igaS,1);CHKERRQ(ierr);													//Number of spatial derivatives to calculate
	ierr = IGASetFromOptions(igaS);CHKERRQ(ierr);												//Note: The order (or degree) of the shape functions is given by the mesh!
	ierr = IGARead(igaS,"./geometry.dat");CHKERRQ(ierr);
	
	for (dir=0; dir<2; dir++)
	{
		ierr = IGASetRuleType(igaS,dir,IGA_RULE_LEGENDRE);CHKERRQ(ierr);
		ierr = IGASetRuleSize(igaS,dir,6);CHKERRQ(ierr);
	}
	ierr = IGASetUp(igaS);CHKERRQ(ierr);
	
	for (dir=0; dir<2; dir++) 
	{
		for (side=0; side<2; side++) 
		{
			//ierr = IGASetBoundaryValue(iga,dir,side,dof,0.0);CHKERRQ(ierr);    				// Dirichlet boundary conditions
			//ierr = IGASetBoundaryForm(igaS,dir,side,PETSC_TRUE);CHKERRQ(ierr);  				// Neumann boundary conditions
		}
	}

	Vec S0;
	Mat K_L2_8GDL;
	Vec F_L2_S;
	ierr = IGACreateVec(igaS,&S0);CHKERRQ(ierr);  
	ierr = IGACreateMat(igaS,&K_L2_8GDL);CHKERRQ(ierr);
	ierr = IGACreateVec(igaS,&F_L2_S);CHKERRQ(ierr);
	ierr = IGASetFormSystem(igaS,L2ProjectionS,&userL2);CHKERRQ(ierr);
	ierr = IGAComputeSystem(igaS,K_L2_8GDL,F_L2_S);CHKERRQ(ierr);

	//This parts set and calls KSP to solve the linear system
	KSP ksp_L2_S;
	ierr = IGACreateKSP(igaS,&ksp_L2_S);CHKERRQ(ierr);										
	ierr = KSPSetOperators(ksp_L2_S,K_L2_8GDL,K_L2_8GDL);CHKERRQ(ierr); 									//This function creates the matrix for the system on the second parameter and uses the 3rd parameter as a preconditioner
	ierr = KSPSetType(ksp_L2_S,KSPCG);CHKERRQ(ierr);											//Using KSPCG (conjugated gradient) because the matrix is symmetric
	ierr = KSPSetFromOptions(ksp_L2_S);CHKERRQ(ierr);
	ierr = KSPSetTolerances(ksp_L2_S,1.0e-35,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);CHKERRQ(ierr);
	ierr = KSPSolve(ksp_L2_S,F_L2_S,S0);CHKERRQ(ierr);											//This is a simple system, so it can be solved with just this command

	ierr = KSPDestroy(&ksp_L2_S);CHKERRQ(ierr);
	ierr = MatDestroy(&K_L2_8GDL);CHKERRQ(ierr);
	ierr = VecDestroy(&F_L2_S);CHKERRQ(ierr);
	char nameS[]="/Input-S-2d-0.dat";
	char pathS[512];
	sprintf(pathS,"%s%s",direct,nameS);
	ierr = IGAWriteVec(igaS,S0,pathS);CHKERRQ(ierr);
//

//Creation of Initialization of Alfa0
	//General parameters
	PetscInt N,M, numEls, counter; 	
	N=0;	M=0;	numEls=1;
	//

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

	PetscReal a=0.0;	
	PetscReal c=Lx/nx;
	PetscReal t=Ly/ny;
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