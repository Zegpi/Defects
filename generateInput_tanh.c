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

		PetscInt i,j,k,l;
		PetscInt nen = p->nen;																//Number of shape functions
		PetscInt dim = p->dim;																//Spatial dimensions of the problem
		PetscInt dof = p->dof;

		PetscReal x[dim];																	//Vector of reals, size equal to problem's dimension
		IGAPointFormGeomMap(p,x);															//Fills x with the coordinates of p, Gauss's point

		//g is the function to L2 project
		PetscReal dx=Lx/nx;
		PetscReal dy=Ly/ny;

		//S has 8 components, in order S(1,1,1), S(1,1,2), S(1,2,1), S(1,2,2), S(2,1,1), S(2,1,2), S(2,2,1), S(2,2,2)

		const PetscReal e[3][3][3]=
		{
			{{0.0,0.0,0.0},{0.0,0.0,1.0},{0.0,-1.0,0.0}},
			{{0.0,0.0,-1.0},{0.0,0.0,0.0},{1.0,0.0,0.0}},
			{{0.0,1.0,0.0},{-1.0,0.0,0.0},{0.0,0.0,0.0}}
		};

		PetscReal a,b,c,s,x_ast1, x_ast2, y_ast1, y_ast2;;
		a= 15.0;
		b=-12.0;
		c= 9.0;
		s= 1.0/3.0;

		//a= 15.0;
		//b=-12.0;
		//c=-12.0;
		//s= 1.0/3.0;

		PetscReal factor=tan(5.0/180.0*ConstPi);
		PetscReal S[3][3][3]={0};

		PetscReal g[dof];
		

		
		/*
		//														     ______
		//This code is for the terraces that look like this:         |
		//													   ______|
		PetscReal V[3][3];
		PetscReal dfx,dfy;	

		dfx=-factor*0.25*( tanh( (x[1]-b)/s )+1.0 )*( 1.0/s *1.0/( cosh((a-x[0])/s)*cosh((a-x[0])/s ) ) )
			+factor*0.25*( tanh( (x[1]-c)/s )+1.0 )*( 1.0/s *1.0/( cosh((a-x[0])/s)*cosh((a-x[0])/s ) ) );

		dfy= factor*0.25*( 1.0/s *1.0/( cosh((b-x[1])/s )*cosh((b-x[1])/s ) ) )*( 1.0-tanh( (x[0]-a)/s ) )
			+factor*0.25*( 1.0/s *1.0/( cosh((c-x[1])/s )*cosh((c-x[1])/s ) ) )*( tanh( (x[0]-a)/s )+1.0 );    

		V[0][0]=0.0; V[0][1]=0.0; V[0][2]=0.0;
		V[1][0]=0.0; V[1][1]=0.0; V[1][2]=0.0;
		V[2][0]=dfx; V[2][1]=dfy; V[2][2]=0.0;

		for(i=0;i<3;i++)
		{
			for(j=0;j<3;j++)
			{
				for(k=0;k<3;k++)
				{
					for(l=0;l<3;l++)
					{
						S[i][j][k]=S[i][j][k]+e[i][j][l]*V[l][k];
					}
				}
			}
		}
		*/
			
		//																	\
		//This code is for the terminated eigenwall that looks like this:    \
		//																	  \
		//
		//Note that this does NOT go to the boundaries, i.e. \Pi is NOT zero

		PetscReal func,func2,inclin,normN;

		a= 0.1;
		b= 10.0;
		s= 1.0/3.0;

		inclin=0.22;					//Factor that changes the inclination of the eigenwall, inclin=1 is a 45° line, values less than one make it steeper, with 0.0 making it vertical.

		func=factor	*(0.5*tanh((x[0]+inclin*x[1]+a)/s)-0.5*tanh((x[0]+inclin*x[1]-a)/s))
					*(0.5*tanh((inclin*x[0]-x[1]+b)/s)-0.5*tanh((inclin*x[0]-x[1]-b)/s));

		PetscReal n[3],W[3][3]={0};

		normN=sqrt(1.0+inclin*inclin);
		n[0]=1.0/normN; n[1]=inclin/normN; n[2]=0.0; 

		W[0][0]=-func;
		W[0][1]=-func;

		for(i=0;i<3;i++)
		{
			for(j=0;j<3;j++)
			{
				for(k=0;k<3;k++)
				{
					S[i][j][k]=S[i][j][k]+W[i][j]*n[k];
				}
			}
		}

		//Recall that S has 8 components, in order S(1,1,1), S(1,1,2), S(1,2,1), S(1,2,2), S(2,1,1), S(2,1,2), S(2,2,1), S(2,2,2)

		func=func/2.2419902995516;			//This factor depends on inclin and n (and probably other things). Must be recalculated for every run so the integral is equal to the factors
											//being set up in the following section.

		x_ast1=(-a+-inclin*b)/(inclin*inclin+1.0);
		y_ast1=inclin*x_ast1+b;

		x_ast2=(a+inclin*b)/(inclin*inclin+1.0);
		y_ast2=inclin*x_ast2-b;

		func2 = (0.5*tanh((x[1]-y_ast1+a)/s)-0.5*tanh((x[1]-y_ast1-a)/s))*(1.0-tanh((x[0]+inclin*x[1]+a)/s))/2.0
			  + (0.5*tanh((x[1]-y_ast2+a)/s)-0.5*tanh((x[1]-y_ast2-a)/s))*(tanh((x[0]+inclin*x[1]-a)/s)+1.0)/2.0;

		func2=factor*func2/2.2419902995516;			//This factor depends on inclin and n (and probably other things). Must be recalculated for every run so the integral is equal to the factors
												//being set up in the following section.			  

		func=1.35316947369629*func;
		func2=1.35316947369629*func2;

		g[0]=-0.5*func;
		g[1]=-0.5*func-func2;
		g[2]= 0.5*func;
		g[3]= 0.5*func;
		g[4]= 0.0;
		g[5]= 0.0;
		g[6]= 0.0;
		g[7]= 0.0;

		//Consider changing all this to just assigning S directly

		const PetscReal (*N) = (typeof(N)) p->shape[0];
		PetscReal (*FF)[dof] = (PetscReal (*)[dof])F;
		PetscReal (*KK)[dof][nen][dof] = (PetscReal (*)[dof][nen][dof])K;
		for(i=0; i<nen; i++)
		{
			for(j=0; j<dof; j++) 
			{
				for(k=0; k<nen; k++)
				{
		  			KK[i][j][k][j] = N[i]*N[k];
		  		}
		  		FF[i][j] = N[i]*g[j];
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

	PetscInt b=601;//1743;//581;					//Parmeter to choose size of cores, must always be odd, core will be of size 1 unit, rest of the body will be of size b-1 units in each direction
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

	strcpy(nombreMalla,"geometry4");
	sprintf(comando,"exec python rectangle4.py %s %f %f %d %d\n",nombreMalla,Lx,Ly,nx,ny);
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
	source = fopen("./generateInput_tanh.c","r");
	dest   = fopen("../Results/generateInput_tanh.c","w");

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
	ierr = KSPSetOperators(ksp_L2_S,K_L2_8GDL,K_L2_8GDL);CHKERRQ(ierr); 						//This function creates the matrix for the system on the second parameter and uses the 3rd parameter as a preconditioner
	PC pcL2;
	ierr = KSPGetPC(ksp_L2_S,&pcL2);CHKERRQ(ierr);
	ierr = PCSetType(pcL2,PCSOR);CHKERRQ(ierr);
	ierr = PCSORSetIterations(pcL2,2,1);CHKERRQ(ierr);
	ierr = KSPSetFromOptions(ksp_L2_S);CHKERRQ(ierr);
	ierr = KSPSetType(ksp_L2_S,KSPFCG);CHKERRQ(ierr);

	ierr = KSPSetTolerances(ksp_L2_S,1.0e-11,1.0e-15,PETSC_DEFAULT,PETSC_DEFAULT);CHKERRQ(ierr);
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