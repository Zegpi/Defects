#include "petiga.h"
#include <stdio.h>
#include <dirent.h>
#include <math.h>

#if PETSC_VERSION_LT(3,5,0)
#define KSPSetOperators(ksp,A,B) KSPSetOperators(ksp,A,B,SAME_NONZERO_PATTERN)
#endif

//Definition of AppCtx
typedef struct 
{
	PetscReal Lx;
	PetscReal Ly;
	PetscReal nx;
	PetscReal ny;
}	AppCtxL2;

typedef struct 
{
	PetscReal dt;
}	AppCtx;

PETSC_STATIC_INLINE
PetscBool IGAElementNextFormSystem(IGAElement element,IGAFormSystem *sys,void **ctx)
{
	IGAForm form = element->parent->form;
	if (!IGAElementNextForm(element,form->visit)) return PETSC_FALSE;
	*sys = form->ops->System;
	*ctx = form->ops->SysCtx;
	return PETSC_TRUE;
}


//L2 projection of Al(0)
	#undef  __FUNCT__
	#define __FUNCT__ "L2ProjectionAl"
	PetscErrorCode L2ProjectionAl(IGAPoint p,PetscReal *K,PetscReal *F,void *ctx)
	{
		if (p->atboundary)
		{
			return 0;														//Condición de borde Neumann nula
		}

		AppCtxL2    *user  = (AppCtxL2 *)ctx;
		PetscReal Lx     = user->Lx;
		PetscReal Ly     = user->Ly;
		PetscReal nx     = user->nx;
		PetscReal ny     = user->ny;

		PetscInt a,b,i;
		PetscInt nen = p->nen;																//Numero de funciones de forma, en este caso es 9
		PetscInt dim = p->dim;																//Dimension del problema
		PetscInt dof = p->dof;

		PetscReal x[dim];																	//Vector de reales, tamaño igual a dimensión del problema
		IGAPointFormGeomMap(p,x);															//llena x con las coordenadas de p, punto de Gauss

		//g es la función a proyectar
		PetscReal dx=Lx/nx;
		PetscReal dy=Ly/ny;
		PetscReal eps=0.005;
		PetscReal g[dof];
		for (i=0; i<dof; i++)
		{
			if (i==0)
			{
				g[i]=1.917545*(0.5*(tanh((x[0]+dx)/eps)+1.0)-0.5*(tanh((x[0]-dx)/eps)+1.0))*(0.5*(tanh((x[1]+dy)/eps)+1.0)-0.5*(tanh((x[1]-dy)/eps)+1.0));
			}
			else if (i==1)
			{
				g[i]=1.917545*(0.5*(tanh((x[0]+10*dx)/eps)+1.0)-0.5*(tanh((x[0]+8*dx)/eps)+1.0))*(0.5*(tanh((x[1]+dy)/eps)+1.0)-0.5*(tanh((x[1]-dy)/eps)+1.0));;
			}
		}

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

//L2 projection of test function
	#undef  __FUNCT__
	#define __FUNCT__ "L2ProjectionF"
	PetscErrorCode L2ProjectionF(IGAPoint p,PetscReal *K,PetscReal *F,void *ctx)
	{
		if (p->atboundary)
		{
			return 0;														//Condición de borde Neumann nula
		}

		PetscInt a,b,i;
		PetscInt nen = p->nen;																//Numero de funciones de forma, en este caso es 9
		PetscInt dim = p->dim;																//Dimension del problema
		PetscInt dof = p->dof;

		PetscReal x[dim];																	//Vector de reales, tamaño igual a dimensión del problema
		IGAPointFormGeomMap(p,x);															//llena x con las coordenadas de p, punto de Gauss

		//f es la función a proyectar
		PetscReal f[dof];
		for (i=0; i<dof; i++)
		{
			if (i==0)
			{
				f[i]=256.0*pow(x[1],4)-32.0*pow(x[1],2)+1;
			}
			else if (i==1)
			{
				f[i]=256.0*pow(x[0],4)-32.0*pow(x[0],2)+1;
			}
		}

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
		  		FF[a][i] = N[a]*f[i];
			}
		}
		return 0;
	}
//

//System for div(grad(z\dot))=curl(alpha x V)
	#undef  __FUNCT__
	#define __FUNCT__ "divGradZ"
	PetscErrorCode divGradZ(IGAPoint p,IGAPoint pAl,PetscReal *K,PetscReal *F, PetscReal *U,void *ctx)
	{
		//AppCtx    *user  = (AppCtx *)ctx;
		//PetscReal dt     = user->dt;

		//Definition of alternating tensor
		const PetscReal e[3][3][3]=
		{
			{{0.0,0.0,0.0},{0.0,0.0,1.0},{0.0,-1.0,0.0}},
			{{0.0,0.0,-1.0},{0.0,0.0,0.0},{1.0,0.0,0.0}},
			{{0.0,1.0,0.0},{-1.0,0.0,0.0},{0.0,0.0,0.0}}
		};
		
		PetscReal alfa[3][3]={0};
		PetscReal dAlfa[3][3][3]={0};
		PetscReal N[3]={0};

		const PetscReal *N0,(*N1)[2];
		IGAPointGetShapeFuns(p,0,(const PetscReal**)&N0);
		IGAPointGetShapeFuns(p,1,(const PetscReal**)&N1);
		PetscInt a,b,i,j,k,l,m,n,nen=p->nen, dof=p->dof;

		PetscReal Al0[2];																		//Create array to recieve Alfa
		PetscReal dAl0[dof][2];																	//Create array to recieve dAlfa
		IGAPointFormValue(pAl,U,&Al0[0]);														//This fills the values
		IGAPointFormGrad (pAl,U,&dAl0[0][0]);													//This fills the values of the derivatives

		PetscReal x[2];																			//Vector de reales, tamaño igual a dimensión del problema
		IGAPointFormGeomMap(p,x);																//llena x con las coordenadas de p, punto de Gauss

		PetscReal V[3]={0};		V[0]=x[0]; V[1]=x[1];
		PetscReal dV[3][3]={0}; dV[0][0]=1.0; dV[0][1]=1.0;

		//if(x[0]>=-0.25/101.0 && x[0]<=0.25/101.0 && x[1]>=-0.25/101.0 && x[1]<=0.25/101.0)
		//{
		//	Al0[0]=1.0; Al0[1]=1.0;
		//}
		//else
		//{
		//	Al0[0]=0.0; Al0[1]=0.0;
		//}

		alfa[0][2]=Al0[0];	alfa[1][2]=Al0[1];												//Assign Al0 to the full alfa expression
		dAlfa[0][2][0]=dAl0[0][0]; dAlfa[0][2][1]=dAl0[0][1];
		dAlfa[1][2][0]=dAl0[1][0]; dAlfa[1][2][1]=dAl0[1][1];

		PetscReal (*Kz)[dof][nen][dof] = (typeof(Kz)) K;
		PetscReal (*Fz)[dof] = (PetscReal (*)[dof])F;

		if (p->atboundary)
		{
			PetscInt dir  = p->boundary_id / 2;
			PetscInt side = p->boundary_id % 2;

			if(dir==0)
			{
				if(side==0)
					N[0]=-1.0;	
				else
					N[0]=1.0;
			}
			else
			{
				if(side==0)
					N[1]=-1.0;	
				else
					N[1]=1.0;
			}

			for (a=0;a<nen;a++)
			{
				PetscReal Na   = N0[a];
				for(i=0;i<dof;i++)
				{
					PetscReal v[3]={0};
					v[i]=Na;

					for(k=0;k<3;k++)
					{
						for(l=0;l<3;l++)
						{
							for(m=0;m<3;m++)
							{
								for(n=0;n<3;n++)
								{
									Fz[a][i]+=e[l][n][m]*alfa[k][n]*V[m]*N[l]*v[k];
								}
							}
						}
					}
				}
			}

			return 0;
		}
		else
		{
			for (a=0; a<nen; a++)
			{
				PetscReal Na   = N0[a];
				PetscReal Na_x = N1[a][0];
				PetscReal Na_y = N1[a][1];

				for (i=0; i<dof; i++)
				{
					PetscReal v[3]={0};
					v[i]=Na;

					PetscReal dv[3][3]={0};
					dv[i][0]=Na_x; dv[i][1]=Na_y;

					for (b=0; b<nen; b++)
					{
						PetscReal Nb_x = N1[b][0];
						PetscReal Nb_y = N1[b][1];

						for(j=0; j<dof; j++)
						{
							PetscReal dz[3][3]={0};
							dz[j][0]=Nb_x; dz[j][1]=Nb_y;

							for(k=0;k<3;k++)
							{	
								for(l=0;l<3;l++)
								{
									Kz[a][i][b][j]+=dz[k][l]*dv[k][l];
								}
							}
						}
					}

					for(k=0;k<3;k++)
					{
						for(l=0;l<3;l++)
						{
							for(m=0;m<3;m++)
							{
								for(n=0;n<3;n++)
								{
									Fz[a][i]+=e[l][n][m]*dAlfa[k][n][l]*V[m]*v[k]+e[l][n][m]*alfa[k][n]*dV[m][l]*v[k];
								}
							}
						}
					}
				}
			}
		}
		return 0;
	}
//

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc, char *argv[]) {
	PetscErrorCode  ierr;

//Delete files in result folder
	DIR *folder = opendir("./results");
	struct dirent *next_file;
	char filepath[256];

	while ( (next_file = readdir(folder)) != NULL )
	{
		sprintf(filepath, "%s/%s", "./results", next_file->d_name);
		remove(filepath);
	}
	closedir(folder);
//

//Generate Mesh
	PetscReal Lx=0.5;
	PetscReal Ly=0.5;
	PetscInt  nx=101;
	PetscInt  ny=101;
	PetscReal h=fmin(Lx/nx,Ly/ny);

	printf("h= %f\n",h);

	char      nombreMalla[12]="geometry2";
	char 	  comando[512];
	sprintf(comando,"exec python rectangle2.py %s %f %f %d %d\n",nombreMalla,Lx,Ly,nx,ny);
	ierr=system(comando);

	//Time parameters
	PetscReal T=0.1;
	PetscReal fact=0.5;
	PetscReal dt=fact*(h/2.0);																	//dt must be less than h/2
	printf("dt<= %f",dt/fact);	printf("     dt= %f\n",dt);
	PetscInt numStep=ceil(T/dt);
	printf("Nt= %d\n",numStep);

//

//Creation of solution systems
	AppCtxL2 userL2;
	userL2.Lx     = Lx;
	userL2.Ly     = Ly;
	userL2.nx     = nx;
	userL2.ny     = ny;
	//void* user;
	ierr = PetscInitialize(&argc,&argv,0,0);CHKERRQ(ierr);										//Always initialize PETSc

	KSP kspl2;

	//System for Alfa
	IGA igaAl;
	ierr = IGACreate(PETSC_COMM_WORLD,&igaAl);CHKERRQ(ierr);
	ierr = IGASetDim(igaAl,2);CHKERRQ(ierr);													//Dimnsion espacial del problema
	ierr = IGASetDof(igaAl,2);CHKERRQ(ierr);													//Number of degrees of freedom, per node
	ierr = IGASetOrder(igaAl,1);CHKERRQ(ierr);													//Number of spatial derivatives to calculate
	ierr = IGASetFromOptions(igaAl);CHKERRQ(ierr);												//Nota: El orden (o grado) de las funciones de forma viene dado por la malla!
	ierr = IGARead(igaAl,"./geometry2.dat");CHKERRQ(ierr);
	ierr = IGASetUp(igaAl);CHKERRQ(ierr);

	PetscInt dir,side;
	for (dir=0; dir<2; dir++) 
	{
		for (side=0; side<2; side++) 
		{
			//ierr = IGASetBoundaryValue(iga,dir,side,dof,0.0);CHKERRQ(ierr);    				// Dirichlet boundary conditions
			ierr = IGASetBoundaryForm(igaAl,dir,side,PETSC_TRUE);CHKERRQ(ierr);  				// Neumann boundary conditions
		}
	}
	
	//System for z\dot
	IGA igaZ;
	ierr = IGACreate(PETSC_COMM_WORLD,&igaZ);CHKERRQ(ierr);
	ierr = IGASetDim(igaZ,2);CHKERRQ(ierr);														//Spatial dimension of the problem
	ierr = IGASetDof(igaZ,2);CHKERRQ(ierr);														//Number of degrees of freedom, per node
	ierr = IGASetOrder(igaZ,1);CHKERRQ(ierr);													//Numero máximo de derivadas a calcular
	ierr = IGASetFromOptions(igaZ);CHKERRQ(ierr);												//Nota: El orden (o grado) de las funciones de forma viene dado por la malla!
	ierr = IGARead(igaZ,"./geometry2.dat");CHKERRQ(ierr);
	ierr = IGASetUp(igaZ);CHKERRQ(ierr);
	//PetscInt dir,side;
	for (dir=0; dir<2; dir++) 
	{
		for (side=0; side<2; side++) 
		{
			//ierr = IGASetBoundaryValue(iga,dir,side,dof,0.0);CHKERRQ(ierr);    				// Dirichlet boundary conditions
			ierr = IGASetBoundaryForm(igaZ,dir,side,PETSC_TRUE);CHKERRQ(ierr);  				// Neumann boundary conditions
		}
	}
//

//Creation of types and systems for the L2 projection of Alfa0
	Vec al0;
	Mat Kl2Al;
	Vec Fl2Al;
	ierr = IGACreateVec(igaAl,&al0);CHKERRQ(ierr);  
	ierr = IGACreateMat(igaAl,&Kl2Al);CHKERRQ(ierr);
	ierr = IGACreateVec(igaAl,&Fl2Al);CHKERRQ(ierr);
	ierr = IGASetFormSystem(igaAl,L2ProjectionAl,&userL2);CHKERRQ(ierr);
	ierr = IGAComputeSystem(igaAl,Kl2Al,Fl2Al);CHKERRQ(ierr);

	//Esta parte llama a KSP para resolver el sistema lineal
	ierr = IGACreateKSP(igaAl,&kspl2);CHKERRQ(ierr);
	ierr = KSPSetOperators(kspl2,Kl2Al,Kl2Al);CHKERRQ(ierr); 								//Esta función crea la matriz a resolver en el segundo parámetro y usa el 3er parámetro como preacondicionador
	ierr = KSPSetType(kspl2,KSPCG);CHKERRQ(ierr);										//Usando KSPCG (gradiente conjugado) porque la matriz es simétrica
	ierr = KSPSetOptionsPrefix(kspl2,"l2pAl_");CHKERRQ(ierr);
	ierr = KSPSetFromOptions(kspl2);CHKERRQ(ierr);
	//ierr = KSPSetTolerances(kspl2,1e-10,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);CHKERRQ(ierr);
	ierr = KSPSolve(kspl2,Fl2Al,al0);CHKERRQ(ierr);										//This is a simple system, so it can be solved with just this command

	ierr = KSPDestroy(&kspl2);CHKERRQ(ierr);
	ierr = MatDestroy(&Kl2Al);CHKERRQ(ierr);
	ierr = VecDestroy(&Fl2Al);CHKERRQ(ierr);
	ierr = IGAWriteVec(igaAl,al0,"./results/Al-2d-0.dat");CHKERRQ(ierr);
//

//Creation of types and systems for the L2 projection of F
	Vec Ff;
	Mat Kl2F;
	Vec Fl2F;
	ierr = IGACreateVec(igaAl,&Ff);CHKERRQ(ierr);  
	ierr = IGACreateMat(igaAl,&Kl2F);CHKERRQ(ierr);
	ierr = IGACreateVec(igaAl,&Fl2F);CHKERRQ(ierr);
	ierr = IGASetFormSystem(igaAl,L2ProjectionF,&userL2);CHKERRQ(ierr);
	ierr = IGAComputeSystem(igaAl,Kl2F,Fl2F);CHKERRQ(ierr);

	//Esta parte llama a KSP para resolver el sistema lineal
	ierr = IGACreateKSP(igaAl,&kspl2);CHKERRQ(ierr);
	ierr = KSPSetOperators(kspl2,Kl2F,Kl2F);CHKERRQ(ierr); 								//Esta función crea la matriz a resolver en el segundo parámetro y usa el 3er parámetro como preacondicionador
	ierr = KSPSetType(kspl2,KSPCG);CHKERRQ(ierr);										//Usando KSPCG (gradiente conjugado) porque la matriz es simétrica
	ierr = KSPSetOptionsPrefix(kspl2,"l2pF_");CHKERRQ(ierr);
	ierr = KSPSetFromOptions(kspl2);CHKERRQ(ierr);
	//ierr = KSPSetTolerances(kspl2,1e-10,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);CHKERRQ(ierr);
	ierr = KSPSolve(kspl2,Fl2F,Ff);CHKERRQ(ierr);										//This is a simple system, so it can be solved with just this command

	ierr = KSPDestroy(&kspl2);CHKERRQ(ierr);
	ierr = MatDestroy(&Kl2F);CHKERRQ(ierr);
	ierr = VecDestroy(&Fl2F);CHKERRQ(ierr);
	ierr = IGAWriteVec(igaAl,Ff,"./results/Ff-2d-0.dat");CHKERRQ(ierr);
//

//Creation of types and systems for z
	Mat Kz;
	Vec z0,Fz;

	ierr = IGACreateMat(igaZ,&Kz);CHKERRQ(ierr);
	ierr = IGACreateVec(igaZ,&z0);CHKERRQ(ierr);
	ierr = IGACreateVec(igaZ,&Fz);CHKERRQ(ierr);

	IGAPoint        pointZ,pointAl;
	IGAElement      elemZ,elemAl;			//element
	PetscReal       *KlocZ,*FlocZ;			//AA y BB
	PetscReal       *KpointZ,*FpointZ;		//KKK y FFF
	const PetscReal *arrayAl0Z;				//arrayU
	Vec  			localAl0Z;				//localU
	PetscReal       *Al0Z;					//U

  	IGAFormSystem  wtfZ;
 	void           *wtf2Z;

 	KSP kspZ;
	ierr = IGACreateKSP(igaZ,&kspZ);CHKERRQ(ierr);

	PetscInt indices[nx+2], indices2[nx+2], indices3[nx+2];

	for(int k=0;k<nx+2;k++)
	{
		indices[k]=2*k;
		indices2[k]=2*(nx+2)+2*k;
		indices3[k]=4*(nx+2)+2*k;
	}

	PetscInt rows=0;
	PetscReal vals; 

	PetscInt i;
	char nombreZ[24];

	i=0;
	//for(i=0;i<4*(nx+1);i=i+2)
	{
		ierr = MatZeroEntries(Kz);CHKERRQ(ierr);
		ierr = VecZeroEntries(Fz);CHKERRQ(ierr);
		ierr = IGAGetLocalVecArray(igaAl,al0,&localAl0Z,&arrayAl0Z);CHKERRQ(ierr);
		//Element loop
		ierr = IGABeginElement(igaZ,&elemZ);CHKERRQ(ierr);
		ierr = IGABeginElement(igaAl,&elemAl);CHKERRQ(ierr);

		while (IGANextElement(igaZ,elemZ))
		{
			IGANextElement(igaAl,elemAl);
			ierr = IGAElementGetWorkMat(elemZ,&KlocZ);CHKERRQ(ierr);
			ierr = IGAElementGetWorkVec(elemZ,&FlocZ);CHKERRQ(ierr);
			ierr = IGAElementGetValues(elemAl,arrayAl0Z,&Al0Z);CHKERRQ(ierr);

			//FormSystem loop
			while (IGAElementNextFormSystem(elemZ,&wtfZ,&wtf2Z)) 
			{
				//Quadrature loop
				ierr = IGAElementBeginPoint(elemZ,&pointZ);CHKERRQ(ierr);
				ierr = IGAElementBeginPoint(elemAl,&pointAl);CHKERRQ(ierr);
				if(pointZ->atboundary==0 && pointAl->atboundary==0)
				{
					while (IGAElementNextPoint(elemZ,pointZ)) 
					{
						IGAElementNextPoint(elemAl,pointAl);
						ierr = IGAPointGetWorkMat(pointZ,&KpointZ);CHKERRQ(ierr);
						ierr = IGAPointGetWorkVec(pointZ,&FpointZ);CHKERRQ(ierr);
						ierr = divGradZ(pointZ,pointAl,KpointZ,FpointZ,Al0Z,NULL);CHKERRQ(ierr);
						ierr = IGAPointAddMat(pointZ,KpointZ,KlocZ);CHKERRQ(ierr);
						ierr = IGAPointAddVec(pointZ,FpointZ,FlocZ);CHKERRQ(ierr);
					}
					IGAElementNextPoint(elemAl,pointAl);
					ierr = IGAElementEndPoint(elemZ,&pointZ);CHKERRQ(ierr);
					ierr = IGAElementEndPoint(elemAl,&pointAl);CHKERRQ(ierr);
				}
			}
			ierr = IGAElementFixSystem(elemZ,KlocZ,FlocZ);CHKERRQ(ierr);					//Esto pone condicion de Dirichlet ¿?
			ierr = IGAElementAssembleMat(elemZ,KlocZ,Kz);CHKERRQ(ierr);
			ierr = IGAElementAssembleVec(elemZ,FlocZ,Fz);CHKERRQ(ierr);

		}
		IGANextElement(igaAl,elemAl);
		ierr = IGAEndElement(igaZ,&elemZ);CHKERRQ(ierr);
		ierr = IGAEndElement(igaAl,&elemAl);CHKERRQ(ierr);

		// Restore local vectors Al0 and arrays
		ierr = IGARestoreLocalVecArray(igaAl,al0,&localAl0Z,&arrayAl0Z);CHKERRQ(ierr);

		ierr = MatAssemblyBegin(Kz,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
		ierr = MatAssemblyEnd  (Kz,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
		//Aquí meto lo que me mandó Hugo
		PetscInt mz,nz;
		ierr = MatGetSize(Kz,&mz,&nz);CHKERRQ(ierr);

		/*
		rows=0;
		vals=1.0;
		ierr = MatZeroRows(Kz,1,&rows,1.0e16,0,0);
		rows=1;
		vals=1.0;
		ierr = MatZeroRows(Kz,1,&rows,1.0e16,0,0);

		ierr = MatAssemblyBegin(Kz,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
		ierr = MatAssemblyEnd  (Kz,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

		//Hasta aquí
		ierr = VecAssemblyBegin(Fz);CHKERRQ(ierr);
		ierr = VecAssemblyEnd  (Fz);CHKERRQ(ierr);
		//Aquí hay que meter lo que mandó hugo

		rows=0;
		vals=0.0;
		ierr = VecSetValue(Fz,rows,vals,INSERT_VALUES);CHKERRQ(ierr);
		ierr = VecAssemblyBegin(Fz);CHKERRQ(ierr);
		ierr = VecAssemblyEnd  (Fz);CHKERRQ(ierr);

		//Hasta aquí
		ierr = KSPSetOperators(kspZ,Kz,Kz);CHKERRQ(ierr);
		ierr = KSPSetFromOptions(kspZ);CHKERRQ(ierr);
		ierr = KSPSolve(kspZ,Fz,z0);CHKERRQ(ierr);
		*/

		
		//rows=2*i;
		//ierr = MatZeroRows(Kz,1,&rows,1.0e16,0,0);
		
		ierr = MatZeroRows(Kz,nx+2,&indices,1.0e16,0,0);
		//ierr = MatZeroRows(Kz,nx+2,&indices2,1.0e16,0,0);
		//ierr = MatZeroRows(Kz,nx+2,&indices3,1.0e16,0,0);
		
		//rows=308;
		//ierr = MatZeroRows(Kz,1,&rows,1.0e16,0,0);
		
		//rows=514;
		//ierr = MatZeroRows(Kz,1,&rows,1.0e16,0,0);
		
		rows=1;
		ierr = MatZeroRows(Kz,1,&rows,1.0e16,0,0);

		ierr = MatAssemblyBegin(Kz,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
		ierr = MatAssemblyEnd  (Kz,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

		ierr = VecAssemblyBegin(Fz);CHKERRQ(ierr);
		ierr = VecAssemblyEnd  (Fz);CHKERRQ(ierr);
		rows=2*i+1;
		vals=0.0;
		ierr = VecSetValue(Fz,rows,vals,INSERT_VALUES);CHKERRQ(ierr);

		for(i=0;i<2*(nx+2);i=i+2)
		{
			ierr = VecSetValue(Fz,i,0,INSERT_VALUES);CHKERRQ(ierr);
		}

		//for(i=2*(nx+2);i<4*(nx+2);i=i+2)
		//{
		//	ierr = VecSetValue(Fz,i,0,INSERT_VALUES);CHKERRQ(ierr);
		//}
		//for(i=4*(nx+2);i<6*(nx+2);i=i+2)
		//{
		//	ierr = VecSetValue(Fz,i,0,INSERT_VALUES);CHKERRQ(ierr);
		//}


		i=0;

		//rows=308;
		//vals=1.0e16;
		//ierr = VecSetValue(Fz,rows,vals,INSERT_VALUES);CHKERRQ(ierr);
		//rows=514;
		//vals=1.e16;
		//ierr = VecSetValue(Fz,rows,vals,INSERT_VALUES);CHKERRQ(ierr);

		ierr = VecAssemblyBegin(Fz);CHKERRQ(ierr);
		ierr = VecAssemblyEnd  (Fz);CHKERRQ(ierr);

		ierr = KSPSetOperators(kspZ,Kz,Kz);CHKERRQ(ierr);
		ierr = KSPSetFromOptions(kspZ);CHKERRQ(ierr);
		ierr = KSPSolve(kspZ,Fz,z0);CHKERRQ(ierr);
		sprintf(nombreZ, "./results/Z-2d-%d.dat", (i+1));	
		ierr = IGAWriteVec(igaZ,z0,nombreZ);CHKERRQ(ierr);
	}

	ierr = KSPDestroy(&kspZ);CHKERRQ(ierr);
	ierr = MatDestroy(&Kz);CHKERRQ(ierr);
	ierr = VecDestroy(&Fz);CHKERRQ(ierr);
	ierr = IGADestroy(&igaZ);CHKERRQ(ierr);
//

return 0;
}