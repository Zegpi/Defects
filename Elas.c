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

PetscReal delta(PetscInt i, PetscInt j)
{
	if (i==j)
		return 1.0;
	else
		return 0.0;
}


//System for equilibrium
	#undef  __FUNCT__
	#define __FUNCT__ "Eq"
	PetscErrorCode Eq(IGAPoint p,IGAPoint pAl,PetscReal *K,PetscReal *F, PetscReal *UAl,void *ctx)		//When other fields like Pi or S (etc.) come into this, declare a IGAPoint pPi or pS and a new PescReal *UPi or *US for each
	{
		PetscReal C[3][3][3][3]={0};														//Definition of elastic tensor
		const PetscReal *N0,(*N1)[2];
		IGAPointGetShapeFuns(p,0,(const PetscReal**)&N0);
		IGAPointGetShapeFuns(p,1,(const PetscReal**)&N1);
		PetscInt a,b,i,j,k,l,u,w,nen=p->nen, dof=p->dof;

		PetscReal Al0[dof];																	//Assign Al(t) to a vector
		PetscReal d_Al0[dof][2];															//Same for partial derivatives
		IGAPointFormValue(pAl,UAl,&Al0[0]);
		IGAPointFormGrad (pAl,UAl,&d_Al0[0][0]);

		PetscReal (*Keq)[dof][nen][dof] = (typeof(Keq)) K;
		PetscReal (*Feq)[dof] = (PetscReal (*)[dof])F;

		PetscInt dim = p->dim;																//Dimension of the problem
		PetscReal x[dim];
		IGAPointFormGeomMap(p,x);

		//E and nu should come from AppCtx in the future
		const PetscReal E=2100000*9.81*10000;
		const PetscReal nu=0.33;
		const PetscReal lambda=(E*nu)/((1+nu)*(1-2*nu));
		const PetscReal mu=E/(2*(1+nu));

		//////////////Delete this parte later, loads should come from appCtx or be 0
		PetscReal f[2]={E , 0.0};
		PetscReal g[2]={0.0 , 0.0};
		/////////////////////////////////////


		for (i=0; i<3; i++)
		{
			for (j=0; j<3; j++)
			{
				for (k=0; k<3; k++)
				{
					for (l=0; l<3; l++)
					{
						C[i][j][k][l]=lambda*delta(i,j)*delta(k,l)+mu*(delta(i,k)*delta(j,l)+delta(i,l)*delta(j,k));
					}
				}
			}
		}

		PetscReal v[2]={0};
		PetscReal eu[3][3]={0};
		PetscReal ev[3][3]={0};

		if (p->atboundary)
		{
			/*
			PetscInt dir  = p->boundary_id / 2;
			PetscInt side = p->boundary_id % 2;

			for (a=0; a<nen; a++)
			{
				PetscReal Na   = N0[a];

				for (i=0; i<dof; i++)
				{
					if (dir==0)
					{
						if(side==1)
						{
							Feq[a][i]+=Na*g[i];
						}
					}
				}
			}
			*/
		}
		else
		{
			for (a=0; a<nen; a++) 
			{
				PetscReal Na_x = N1[a][0];
				PetscReal Na_y = N1[a][1];

				for (b=0; b<nen; b++) 
				{
					PetscReal Nb_x = N1[b][0];
					PetscReal Nb_y = N1[b][1];

					for (i=0; i<dof; i++)
					{
						if (i==0)
						{
							ev[0][0]=Na_x; ev[0][1]=0.5*Na_y; ev[1][0]=0.5*Na_y; ev[1][1]=0.0; ev[2][0]=0.0; ev[2][1]=0.0; ev[2][2]=0.0; ev[0][2]=0.0; ev[1][2]=0.0;
						}
						else if(i==1)
						{
							ev[0][0]=0.0; ev[0][1]=0.5*Na_x; ev[1][0]=0.5*Na_x; ev[1][1]=Na_y; ev[2][0]=0.0; ev[2][1]=0.0; ev[2][2]=0.0; ev[0][2]=0.0; ev[1][2]=0.0;
						}

						for (j=0; j<dof; j++)
						{
							if (j==0)
							{
								eu[0][0]=Nb_x; eu[0][1]=0.5*Nb_y; eu[1][0]=0.5*Nb_y; eu[1][1]=0.0; eu[2][0]=0.0; eu[2][1]=0.0; eu[2][2]=0.0; eu[0][2]=0.0; eu[1][2]=0.0;
							}
							else if(j==1)
							{
								eu[0][0]=0.0; eu[0][1]=0.5*Nb_x; eu[1][0]=0.5*Nb_x; eu[1][1]=Nb_y; eu[2][0]=0.0; eu[2][1]=0.0; eu[2][2]=0.0; eu[0][2]=0.0; eu[1][2]=0.0;
							}

							/*
							Keq[a][i][b][j]=lambda*N1[a][i]*N1[b][j];
							for (k=0; k<3; k++)
							{
								for (l=0; l<3; l++)
								{
									Keq[a][i][b][j]+=2*mu*eu[k][l]*ev[k][l];
								}
							}
							*/

							
							Keq[a][i][b][j]=0.0;
							for (k=0; k<3; k++)
							{
								for (l=0; l<3; l++)
								{
									for (u=0; u<3; u++)
									{
										for (w=0; w<3; w++)
										{
											Keq[a][i][b][j]+=C[k][l][u][w]*eu[u][w]*ev[k][l];
										}
									}
								}
							}
						}
					}
				}
			}

			for(a=0 ;a<nen; a++)
			{
				for (i=0; i<dof; i++)
				{
					v[0]=0; v[1]=0;
					v[i]=N0[a];

					Feq[a][i] = v[i]*f[i];
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
//

//Generate Mesh
	PetscReal Lx=5.0;
	PetscReal Ly=5.0;
	PetscInt  nx=101;
	PetscInt  ny=101;
	PetscReal h=fmin(Lx/nx,Ly/ny);

	printf("h= %f\n",h);

	char      nombreMalla[12]="geometry";
	char 	  comando[512];
	sprintf(comando,"exec python rectangle.py %s %f %f %d %d\n",nombreMalla,Lx,Ly,nx,ny);
	ierr=system(comando);

	//Time parameters
	PetscReal T=0.05;
	PetscReal fact=0.5;
	PetscReal dt=fact*(h/2.0);																	//dt must be less than h/2
	printf("dt<= %f     dt= %f\n",dt/fact,dt);
	PetscInt numStep=ceil(T/dt);
	printf("Nt= %d\n",numStep);
//

//Creation of solution systems
	ierr = PetscInitialize(&argc,&argv,0,0);CHKERRQ(ierr);										//Always initialize PETSc
	PetscInt dir,side;
//

//App context creation
	AppCtxL2 userL2;
	userL2.Lx     = Lx;
	userL2.Ly     = Ly;
	userL2.nx     = nx;
	userL2.ny     = ny;
	//void* user;
//

	//System for equilibrium
	IGA igaEq;
	ierr = IGACreate(PETSC_COMM_WORLD,&igaEq);CHKERRQ(ierr);
	ierr = IGASetDim(igaEq,2);CHKERRQ(ierr);													//Spatial dimension of the problem
	ierr = IGASetDof(igaEq,2);CHKERRQ(ierr);													//Number of degrees of freedom, per node
	ierr = IGASetOrder(igaEq,1);CHKERRQ(ierr);													//Number of spatial derivatives to calculate
	ierr = IGASetFromOptions(igaEq);CHKERRQ(ierr);												//Note: The order (or degree) of the shape functions is given by the mesh!
	ierr = IGARead(igaEq,"./geometry.dat");CHKERRQ(ierr);
	ierr = IGASetUp(igaEq);CHKERRQ(ierr);

	Mat KEq;
	Vec u0,FEq;

	ierr = IGACreateMat(igaEq,&KEq);CHKERRQ(ierr);
	ierr = IGACreateVec(igaEq,&u0);CHKERRQ(ierr);
	ierr = IGACreateVec(igaEq,&FEq);CHKERRQ(ierr);

	IGAPoint		pointEq;					//point
	IGAElement		elemEq;						//element
	PetscReal		*KlocEq,*FlocEq;			//AA y BB
	PetscReal		*KpointEq,*FpointEq;		//KKK y FFF
	const PetscReal	*arrayAl0Eq;				//arrayU
	Vec				localAl0Eq;					//localU
	PetscReal		*Al0Eq;						//U

  	IGAFormSystem	wtfEq;
 	void			*wtf2Eq;

 	KSP kspEq;
	ierr = IGACreateKSP(igaEq,&kspEq);CHKERRQ(ierr);

	//ierr = IGASetBoundaryValue(iga,dir,side,dof,value);CHKERRQ(ierr);							//Syntax
	//ierr = IGASetBoundaryValue(igaEq,1 ,0 ,0 ,0.0);CHKERRQ(ierr);    							// Dirichlet boundary conditions x-dir, left side, u(0)=0
	//ierr = IGASetBoundaryValue(igaEq,1, 0, 1, 0.0);CHKERRQ(ierr);    							// Dirichlet boundary conditions x-dir, left side, u(1)=0 //This means clamped in left side

	//ierr = IGASetBoundaryForm(iga,dir,side,PETSC_TRUE);CHKERRQ(ierr);								// Syntax	
	//ierr = IGASetBoundaryForm(igaEq,0,1,PETSC_TRUE);CHKERRQ(ierr);								// Neumann boundary conditions right side
	//ierr = IGASetBoundaryForm(igaEq,1,0,PETSC_TRUE);CHKERRQ(ierr);								// Neumann boundary conditions lower side
	//ierr = IGASetBoundaryForm(igaEq,1,1,PETSC_TRUE);CHKERRQ(ierr);								// Neumann boundary conditions upper side


	// Get local vectors Al0 and arrays
	ierr = IGAGetLocalVecArray(igaAl,al0,&localAl0Eq,&arrayAl0Eq);CHKERRQ(ierr);

	// Element loop
	ierr = IGABeginElement(igaEq,&elemEq);CHKERRQ(ierr);
	ierr = IGABeginElement(igaAl,&elemAl);CHKERRQ(ierr);
	while (IGANextElement(igaEq,elemEq)) 
	{
		IGANextElement(igaAl,elemAl);
		ierr = IGAElementGetWorkMat(elemEq,&KlocEq);CHKERRQ(ierr);
		ierr = IGAElementGetWorkVec(elemEq,&FlocEq);CHKERRQ(ierr);
		ierr = IGAElementGetValues(elemAl,arrayAl0Eq,&Al0Eq);CHKERRQ(ierr);

		// FormSystem loop
		while (IGAElementNextFormSystem(elemEq,&wtfEq,&wtf2Eq)) 
		{
		// Quadrature loop
			ierr = IGAElementBeginPoint(elemEq,&pointEq);CHKERRQ(ierr);
			ierr = IGAElementBeginPoint(elemAl,&pointAl);CHKERRQ(ierr);
			if(pointEq->atboundary==0 && pointAl->atboundary==0)
			{
				while (IGAElementNextPoint(elemEq,pointEq)) 
				{
					IGAElementNextPoint(elemAl,pointAl);
					ierr = IGAPointGetWorkMat(pointEq,&KpointEq);CHKERRQ(ierr);
					ierr = IGAPointGetWorkVec(pointEq,&FpointEq);CHKERRQ(ierr);
					ierr = Eq(pointEq,pointAl,KpointEq,FpointEq,Al0Eq,NULL);CHKERRQ(ierr);
					ierr = IGAPointAddMat(pointEq,KpointEq,KlocEq);CHKERRQ(ierr);
					ierr = IGAPointAddVec(pointEq,FpointEq,FlocEq);CHKERRQ(ierr);
				}
				IGAElementNextPoint(elemAl,pointAl);
				ierr = IGAElementEndPoint(elemEq,&pointEq);CHKERRQ(ierr);
				ierr = IGAElementEndPoint(elemAl,&pointAl);CHKERRQ(ierr);
			}
		}
		ierr = IGAElementFixSystem(elemEq,KlocEq,FlocEq);CHKERRQ(ierr);					//This sets Dirichlet condition ¿?
		ierr = IGAElementAssembleMat(elemEq,KlocEq,KEq);CHKERRQ(ierr);
		ierr = IGAElementAssembleVec(elemEq,FlocEq,FEq);CHKERRQ(ierr);

	}
	IGANextElement(igaAl,elemAl);
	ierr = IGAEndElement(igaEq,&elemEq);CHKERRQ(ierr);
	ierr = IGAEndElement(igaAl,&elemAl);CHKERRQ(ierr);

	// Restore local vectors Al0 and arrays
	ierr = IGARestoreLocalVecArray(igaAl,al0,&localAl0Eq,&arrayAl0Eq);CHKERRQ(ierr);

	ierr = MatAssemblyBegin(KEq,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyEnd  (KEq,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	

	//Here we set values to the matrix directly (to impose Dirichlet condition in a single point)
	ierr = MatGetSize(KEq,&mz,&nz);CHKERRQ(ierr);

	PetscInt rowsEq=0;

	/*
	for (int ind=0; ind<204; ind++)
	{
	rowsEq=ind;
	ierr = MatZeroRows(KEq,1,&rowsEq,1.0e30,0,0);
	}
	*/

	ierr = MatAssemblyBegin(KEq,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyEnd  (KEq,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);


	ierr = VecAssemblyBegin(FEq);CHKERRQ(ierr);
	ierr = VecAssemblyEnd  (FEq);CHKERRQ(ierr);

	ierr = KSPSetOperators(kspEq,KEq,KEq);CHKERRQ(ierr);
	ierr = KSPSetFromOptions(kspEq);CHKERRQ(ierr);
	ierr = KSPSolve(kspEq,FEq,u0);CHKERRQ(ierr);

	//ierr = MatGetSize(KEq,&mz,&nz);CHKERRQ(ierr);
	//printf("El tamaño de la matriz es %d X %d \n",mz,nz);

	//ierr = VecGetSize(FEq,&mz);CHKERRQ(ierr);
	//printf("El tamaño del vector es %d \n",mz);

	/*
	PetscReal *array;
	VecGetArray(FEq, &array);
	for (int index=0; index<mz; index++)
	{
		printf("%.8f \n",array[index]);
	}
	*/

	ierr = KSPDestroy(&kspEq);CHKERRQ(ierr);
	ierr = MatDestroy(&KEq);CHKERRQ(ierr);
	ierr = VecDestroy(&FEq);CHKERRQ(ierr);
	//ierr = IGAWriteVec(igaEq,u0,"./results/u-2d-0.dat");CHKERRQ(ierr);
	char nameEq[]="/u-2d-0.dat";
	char pathEq[512];
	sprintf(pathEq,"%s%s",direct,nameEq);
	ierr = IGAWriteVec(igaEq,u0,pathEq);CHKERRQ(ierr);
	ierr = IGADestroy(&igaEq);CHKERRQ(ierr);
	ierr = IGADestroy(&igaAl);CHKERRQ(ierr);
//

ierr = PetscFinalize();CHKERRQ(ierr);

return 0;
}