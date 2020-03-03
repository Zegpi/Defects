#include "petiga.h"
#include <stdio.h>
#include <dirent.h>
#include <math.h>

#if PETSC_VERSION_LT(3,5,0)
#define KSPSetOperators(ksp,A,B) KSPSetOperators(ksp,A,B,SAME_NONZERO_PATTERN)
#endif

#define ConstPi 3.14159265358979323846

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

//System for Energy
	#undef  __FUNCT__
	#define __FUNCT__ "Energy"
	PetscErrorCode Energy(IGAPoint p,IGAPoint pChi,IGAPoint pZu,PetscReal *F, PetscReal *Chi,PetscReal *Zu,void *ctx)
	{
		PetscInt a,b,i,j,k,l,m,n,nen=p->nen, dof=p->dof;

		PetscReal x[2];																	//Vector of reals, size equal to problem's dimension
		IGAPointFormGeomMap(p,x);															//Fills x with the coordinates of p, Gauss's point

		PetscReal chi0[4], d_chi0[4][2];													//Array to contain the vector chi(0)
		IGAPointFormValue(pChi,Chi,&chi0[0]);												//Assign chi to its container
		IGAPointFormGrad (pChi,Chi,&d_chi0[0][0]);											//Same for the gradient

		PetscReal d_Z0[2][2], d2_Z0[2][2][2];												//Array to contain the gradient of z(0)
		IGAPointFormGrad (pZu,Zu,&d_Z0[0][0]);												//Assign grad(z0) to array
		IGAPointFormHess (pZu,Zu,&d2_Z0[0][0][0]);											//Same for the hessian (second derivatives)

		//const PetscReal E=1.0; //2100000.0*9.81*10000.0;					//[Pa]
		//const PetscReal nu=0.33;
		//const PetscReal lambda=(E*nu)/((1.0+nu)*(1.0-2.0*nu));
		//const PetscReal mu=E/(2.0*(1.0+nu));
		//const PetscReal eps=0.0*E/1000.0;							//Choose later based on whatever Amit says :)

		//Change to consider G=1
		const PetscReal nu=0.33;
		const PetscReal mu=1.0;
		const PetscReal lambda=2.0*mu*nu/(1.0-2.0*nu);
		const PetscReal eps=0.0*mu/100.0;

		PetscReal C[3][3][3][3]={0};														//Initialization of elastic tensor
		for (i=0; i<3; i++)
		{
			for (j=0; j<3; j++)
			{
				for (k=0; k<3; k++)
				{
					for (l=0; l<3; l++)
					{
						C[i][j][k][l]=lambda*delta(i,j)*delta(k,l)+mu*(delta(i,k)*delta(j,l)+delta(i,l)*delta(j,k));		//Definition of elastic tensor
					}
				}
			}
		}

		//The four non-zero components of Chi are stored as a vector, restore them to an array with the correct indexing for value and derivative
		PetscReal fullChi[3][3]={0};
		fullChi[0][0]=chi0[0]; 	fullChi[0][1]=chi0[1];
		fullChi[1][0]=chi0[2]; 	fullChi[1][1]=chi0[3];

		PetscReal fulld_Chi[3][3][3]={0};
		fulld_Chi[0][0][0]=d_chi0[0][0]; fulld_Chi[0][0][1]=d_chi0[0][1];
		fulld_Chi[0][1][0]=d_chi0[1][0]; fulld_Chi[0][1][1]=d_chi0[1][1];
		fulld_Chi[1][0][0]=d_chi0[2][0]; fulld_Chi[1][0][1]=d_chi0[2][1];
		fulld_Chi[1][1][0]=d_chi0[3][0]; fulld_Chi[1][1][1]=d_chi0[3][1];


		//Expanding z (and derivatives) to 3 components, more convenient for sums in for loops
		PetscReal fulld_z[3][3]={0};
		fulld_z[0][0]=d_Z0[0][0]; fulld_z[0][1]=d_Z0[0][1];
		fulld_z[1][0]=d_Z0[1][0]; fulld_z[1][1]=d_Z0[1][1];

		PetscReal fulld2_z[3][3][3]={0};
		fulld2_z[0][0][0]=d2_Z0[0][0][0]; fulld2_z[0][0][1]=d2_Z0[0][0][1];
		fulld2_z[0][1][0]=d2_Z0[0][1][0]; fulld2_z[0][1][1]=d2_Z0[0][1][1];
		fulld2_z[1][0][0]=d2_Z0[1][0][0]; fulld2_z[1][0][1]=d2_Z0[1][0][1];
		fulld2_z[1][1][0]=d2_Z0[1][1][0]; fulld2_z[1][1][1]=d2_Z0[1][1][1];

		PetscReal (*Fstress)[dof] = (PetscReal (*)[dof])F;

		if (p->atboundary)
		{
			return 0;
		}
		else
		{
			//if (x[0]>-3.0 && x[0]<3.0 && x[1]>-3.0 && x[1]<3.0)
			//{
			//	Fstress[0][0]=0.0;
			//}
			//else
			//{
				Fstress[0][0]=0.0;
				for(int i=0; i<3;i++)
				{
					for(int j=0;j<3;j++)
					{
						for(int k=0;k<3;k++)
						{
							for(int l=0;l<3;l++)
							{
								Fstress[0][0]+=0.5*(fulld_z[i][j]-fullChi[i][j])*C[i][j][k][l]*(fulld_z[k][l]-fullChi[k][l]);
							}
						}
					}
				}

				for(int i=0; i<3;i++)
				{
					for(int j=0;j<3;j++)
					{
						for(int k=0;k<3;k++)
						{
							Fstress[0][0]+=0.5*eps*(fulld2_z[i][j][k]-fulld_Chi[i][j][k])*(fulld2_z[i][j][k]-fulld_Chi[i][j][k]);
						}
					}
				}
			//}

		}
		return 0;
	}
//

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc, char *argv[]) {
//Creation of solution systems
	PetscErrorCode  ierr;
	ierr = PetscInitialize(&argc,&argv,0,0);CHKERRQ(ierr);										//Always initialize PETSc
	PetscInt dir,side;
	PetscPrintf(PETSC_COMM_WORLD,"Start of Energy \n");
//

//Check for folders for files
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
//

//System for Energy
	PetscPrintf(PETSC_COMM_WORLD,"System for Energy starting \n");
	IGA igaEnergy, igaZ0, igachiUp;

	ierr = IGACreate(PETSC_COMM_WORLD,&igaZ0);CHKERRQ(ierr);
	ierr = IGASetDim(igaZ0,2);CHKERRQ(ierr);													//Spatial dimension of the problem
	ierr = IGASetDof(igaZ0,2);CHKERRQ(ierr);													//Number of degrees of freedom, per node
	ierr = IGASetOrder(igaZ0,3);CHKERRQ(ierr);													//Number of spatial derivatives to calculate
	ierr = IGASetFromOptions(igaZ0);CHKERRQ(ierr);												//Note: The order (or degree) of the shape functions is given by the mesh!
	ierr = IGARead(igaZ0,"./geometry3.dat");CHKERRQ(ierr);
	ierr = IGASetUp(igaZ0);CHKERRQ(ierr);

	ierr = IGACreate(PETSC_COMM_WORLD,&igachiUp);CHKERRQ(ierr);
	ierr = IGASetDim(igachiUp,2);CHKERRQ(ierr);													//Spatial dimension of the problem
	ierr = IGASetDof(igachiUp,4);CHKERRQ(ierr);													//Number of degrees of freedom, per node
	ierr = IGASetOrder(igachiUp,2);CHKERRQ(ierr);												//Number of spatial derivatives to calculate
	ierr = IGASetFromOptions(igachiUp);CHKERRQ(ierr);											//Note: The order (or degree) of the shape functions is given by the mesh!
	ierr = IGARead(igachiUp,"./geometry2.dat");CHKERRQ(ierr);
	ierr = IGASetUp(igachiUp);CHKERRQ(ierr);

	ierr = IGACreate(PETSC_COMM_WORLD,&igaEnergy);CHKERRQ(ierr);
	ierr = IGASetDim(igaEnergy,2);CHKERRQ(ierr);												//Spatial dimension of the problem
	ierr = IGASetDof(igaEnergy,1);CHKERRQ(ierr);												//Number of degrees of freedom, per node
	ierr = IGASetOrder(igaEnergy,2);CHKERRQ(ierr);												//Number of spatial derivatives to calculate
	ierr = IGASetFromOptions(igaEnergy);CHKERRQ(ierr);											//Note: The order (or degree) of the shape functions is given by the mesh!
	ierr = IGARead(igaEnergy,"./geometry3.dat");CHKERRQ(ierr);
	ierr = IGASetUp(igaEnergy);CHKERRQ(ierr);

	Vec FEnergy,chiUp0,Z0;

	ierr = IGACreateVec(igaEnergy,&FEnergy);CHKERRQ(ierr);
	ierr = IGACreateVec(igachiUp,&chiUp0);CHKERRQ(ierr);
	ierr = IGACreateVec(igaZ0,&Z0);CHKERRQ(ierr);

	IGAPoint		pointZ0, pointchiUp, pointEnergy;								//point
	IGAElement		elemZ0, elemchiUp, elemEnergy;									//element
	PetscReal		*FlocEnergy;																//AA y BB
	PetscReal		*FpointEnergy;																//KKK y FFF
	const PetscReal	*arrayZ0Energy,*arrayChi0Energy;							//arrayU
	Vec				localZ0Energy,localChi0Energy;								//localU
	PetscReal		*Chi0Energy,*Z0Energy;											//U

  	IGAFormSystem	wtfEnergy;
 	void			*wtf2Energy;

	// Get local vectors Z0 and Chi0 and arrays
 	//PetscErrorCode IGAReadVec(IGA iga,Vec vec,const char filename[])
 	char pathChi[512];
 	char nameChi[]="/ChiUp-2d-0.dat";
 	sprintf(pathChi,"%s%s",direct,nameChi);
 	ierr = IGAReadVec(igachiUp,chiUp0,pathChi);CHKERRQ(ierr);
	ierr = IGAGetLocalVecArray(igachiUp,chiUp0,&localChi0Energy,&arrayChi0Energy);CHKERRQ(ierr);

	char pathZ[512];
 	char nameZ[]="/Z0-2d-0.dat";
 	sprintf(pathZ,"%s%s",direct,nameZ);
 	ierr = IGAReadVec(igaZ0,Z0,pathZ);CHKERRQ(ierr);
	ierr = IGAGetLocalVecArray(igaZ0,Z0,&localZ0Energy,&arrayZ0Energy);CHKERRQ(ierr);

	// Element loop
	ierr = IGABeginElement(igaEnergy,&elemEnergy);CHKERRQ(ierr);
	ierr = IGABeginElement(igaZ0,&elemZ0);CHKERRQ(ierr);
	ierr = IGABeginElement(igachiUp,&elemchiUp);CHKERRQ(ierr);

	while (IGANextElement(igaEnergy,elemEnergy)) 
	{
		IGANextElement(igaZ0,elemZ0);
		IGANextElement(igachiUp,elemchiUp);

		ierr = IGAElementGetWorkVec(elemEnergy,&FlocEnergy);CHKERRQ(ierr);
		ierr = IGAElementGetValues(elemZ0,arrayZ0Energy,&Z0Energy);CHKERRQ(ierr);
		ierr = IGAElementGetValues(elemchiUp,arrayChi0Energy,&Chi0Energy);CHKERRQ(ierr);

		// FormSystem loop
		while (IGAElementNextFormSystem(elemEnergy,&wtfEnergy,&wtf2Energy)) 
		{
		// Quadrature loop
			ierr = IGAElementBeginPoint(elemEnergy,&pointEnergy);CHKERRQ(ierr);
			ierr = IGAElementBeginPoint(elemZ0,&pointZ0);CHKERRQ(ierr);
			ierr = IGAElementBeginPoint(elemchiUp,&pointchiUp);CHKERRQ(ierr);

			while (IGAElementNextPoint(elemEnergy,pointEnergy))
			{
				if(pointEnergy->atboundary==1)
				{
					
				}

				if(pointEnergy->atboundary==0 && pointZ0->atboundary==0 && pointchiUp->atboundary==0)
				{
					IGAElementNextPoint(elemZ0,pointZ0);
					IGAElementNextPoint(elemchiUp,pointchiUp);

					ierr = IGAPointGetWorkVec(pointEnergy,&FpointEnergy);CHKERRQ(ierr);
					//PetscErrorCode Energy(IGAPoint p,IGAPoint pChi,IGAPoint pZu,IGAPoint pStress,PetscReal *F, PetscReal *Chi,PetscReal *Zu,PetscReal *S,void *ctx)
					//ierr = Energy(pointEnergy,pointchiUp,pointZ0,pointStress,FpointEnergy,Chi0Energy,Z0Energy,SigmaEnergy,NULL);CHKERRQ(ierr);
					ierr = Energy(pointEnergy,pointchiUp,pointZ0,FpointEnergy,Chi0Energy,Z0Energy,NULL);CHKERRQ(ierr);
					ierr = IGAPointAddVec(pointEnergy,FpointEnergy,FlocEnergy);CHKERRQ(ierr);
				}
			}
			while (pointZ0->index != -1)
			{
				IGAElementNextPoint(elemZ0,pointZ0);
			}
			while (pointchiUp->index != -1)
			{
				IGAElementNextPoint(elemchiUp,pointchiUp);
			}

			ierr = IGAElementEndPoint(elemEnergy,&pointEnergy);CHKERRQ(ierr);
			ierr = IGAElementEndPoint(elemZ0,&pointZ0);CHKERRQ(ierr);
			ierr = IGAElementEndPoint(elemchiUp,&pointchiUp);CHKERRQ(ierr);
		}

		//ierr = IGAElementFixSystem(elemStress,KlocStress,FlocStress);CHKERRQ(ierr);					//This sets Dirichlet condition ¿? (Yes, this applies the conditions from IGASetBoundaryValue)
		ierr = IGAElementAssembleVec(elemEnergy,FlocEnergy,FEnergy);CHKERRQ(ierr);
	}
	IGANextElement(igaZ0,elemZ0);
	IGANextElement(igachiUp,elemchiUp);

	ierr = IGAEndElement(igaEnergy,&elemEnergy); CHKERRQ(ierr);
	ierr = IGAEndElement(igaZ0,&elemZ0);CHKERRQ(ierr);
	ierr = IGAEndElement(igachiUp,&elemchiUp);CHKERRQ(ierr);

	// Restore local vectors u, Z0, Chi0 and arrays
	ierr = IGARestoreLocalVecArray(igaZ0,Z0,&localZ0Energy,&arrayZ0Energy);CHKERRQ(ierr);
	ierr = IGARestoreLocalVecArray(igachiUp,chiUp0,&localChi0Energy,&arrayChi0Energy);CHKERRQ(ierr);

	ierr = VecAssemblyBegin(FEnergy);CHKERRQ(ierr);
	ierr = VecAssemblyEnd  (FEnergy);CHKERRQ(ierr);

	PetscReal result=0.0;
	ierr = VecSum(FEnergy,&result);

	PetscPrintf(PETSC_COMM_WORLD,"\n %f \n \n",result);

	ierr = VecDestroy(&FEnergy);CHKERRQ(ierr);

//

//Destroy all objects not needed anymore (Better to do it here in case different codes call the same IGA, move if memory is a problem)
	ierr = IGADestroy(&igaEnergy);CHKERRQ(ierr);
	ierr = IGADestroy(&igaZ0);CHKERRQ(ierr);
	ierr = IGADestroy(&igachiUp);CHKERRQ(ierr);
//

ierr = PetscFinalize();CHKERRQ(ierr);

return 0;
}