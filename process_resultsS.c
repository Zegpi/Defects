#include "petiga.h"
#include <stdio.h>
#include <dirent.h>
#include <math.h>
#include <time.h>

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

//System for L2 projection of stress (sym part)
	#undef  __FUNCT__
	#define __FUNCT__ "Stress"
	//PetscErrorCode Stress(IGAPoint p,IGAPoint pU, IGAPoint pHs,IGAPoint pChi,IGAPoint pZu,PetscReal *K,PetscReal *F,PetscReal *U,PetscReal *HS, PetscReal *Chi,PetscReal *Zu,void *ctx)		//When other fields like Pi or S (etc.) come into this, declare a IGAPoint pPi or pS and a new PescReal *UPi or *US for each
	PetscErrorCode Stress(IGAPoint p,IGAPoint pChi,IGAPoint pZu,IGAPoint pS, IGAPoint pGradZ, PetscReal *K,PetscReal *F, PetscReal *Chi,PetscReal *Zu,PetscReal *S, PetscReal *UgradZ, void *ctx)		//When other fields like Pi or S (etc.) come into this, declare a IGAPoint pPi or pS and a new PestcReal *UPi or *US for each
	{
		const PetscReal *N0;
		IGAPointGetShapeFuns(p,0,(const PetscReal**)&N0);									//Value of the shape functions
		PetscInt a,b,i,j,k,l,m,n,nen=p->nen, dof=p->dof;

		//Change to consider G=1
		const PetscReal nu=0.33;
		const PetscReal mu=1.0;
		const PetscReal lambda=2.0*mu*nu/(1.0-2.0*nu);
		const PetscReal eps=mu/100.0;															//Choose later based on whatever Amit says :)

		PetscReal chi0[4];																	//Array to contain the vector chi(0)
		PetscReal d2_Chi0[4][2][2];															//Same for its Hessian
		IGAPointFormValue(pChi,Chi,&chi0[0]);												//Assign chi to its container
		IGAPointFormHess (pChi,Chi,&d2_Chi0[0][0][0]);										//This should be the 3-rd order tensor Chi_{i,jk} (remember that we are storing Chi_{mn} as a column vector Chi_{i})

		//The four non-zero components of Chi are stored as a vector, restore them to an array with the correct indexing for value and derivative
		PetscReal fullChi[3][3]={0};
		fullChi[0][0]=chi0[0]; 	fullChi[0][1]=chi0[1];
		fullChi[1][0]=chi0[2]; 	fullChi[1][1]=chi0[3];

		PetscReal fulld2_Chi[3][3][3][3]={0};
		fulld2_Chi[0][0][0][0]=d2_Chi0[0][0][0]; fulld2_Chi[0][0][0][1]=d2_Chi0[0][0][1]; fulld2_Chi[0][0][1][0]=d2_Chi0[0][1][0]; fulld2_Chi[0][0][1][1]=d2_Chi0[0][1][1];
		fulld2_Chi[0][1][0][0]=d2_Chi0[1][0][0]; fulld2_Chi[0][1][0][1]=d2_Chi0[1][0][1]; fulld2_Chi[0][1][1][0]=d2_Chi0[1][1][0]; fulld2_Chi[0][1][1][1]=d2_Chi0[1][1][1];
		fulld2_Chi[1][0][0][0]=d2_Chi0[2][0][0]; fulld2_Chi[1][0][0][1]=d2_Chi0[2][0][1]; fulld2_Chi[1][0][1][0]=d2_Chi0[2][1][0]; fulld2_Chi[1][0][1][1]=d2_Chi0[2][1][1];
		fulld2_Chi[1][1][0][0]=d2_Chi0[3][0][0]; fulld2_Chi[1][1][0][1]=d2_Chi0[3][0][1]; fulld2_Chi[1][1][1][0]=d2_Chi0[3][1][0]; fulld2_Chi[1][1][1][1]=d2_Chi0[3][1][1];

		PetscReal d_Z0[2][2];																//Same for its gradient
		PetscReal d3_Z0[2][2][2][2];														//Same for its 3rd order partial derivatives
		IGAPointFormGrad (pZu,Zu,&d_Z0[0][0]);												//Same for the gradient
		IGAPointFormDer3 (pZu,Zu,&d3_Z0[0][0][0][0]);										//Same for the thir derivatives

		//Expanding z (and derivatives) to 3 components, more convenient for sums in for loops
		PetscReal fulld_z[3][3]={0};
		fulld_z[0][0]=d_Z0[0][0]; fulld_z[0][1]=d_Z0[0][1];
		fulld_z[1][0]=d_Z0[1][0]; fulld_z[1][1]=d_Z0[1][1];

		PetscReal fulld3_z[3][3][3][3]={0};
		fulld3_z[0][0][0][0]=d3_Z0[0][0][0][0]; fulld3_z[0][0][0][1]=d3_Z0[0][0][0][1]; fulld3_z[0][0][1][0]=d3_Z0[0][0][1][0]; fulld3_z[0][0][1][1]=d3_Z0[0][0][1][1];
		fulld3_z[0][1][0][0]=d3_Z0[0][1][0][0]; fulld3_z[0][1][0][1]=d3_Z0[0][1][0][1]; fulld3_z[0][1][1][0]=d3_Z0[0][1][1][0]; fulld3_z[0][1][1][1]=d3_Z0[0][1][1][1];
		fulld3_z[1][0][0][0]=d3_Z0[1][0][0][0]; fulld3_z[1][0][0][1]=d3_Z0[1][0][0][1]; fulld3_z[1][0][1][0]=d3_Z0[1][0][1][0]; fulld3_z[1][0][1][1]=d3_Z0[1][0][1][1]; 
		fulld3_z[1][1][0][0]=d3_Z0[1][1][0][0]; fulld3_z[1][1][0][1]=d3_Z0[1][1][0][1]; fulld3_z[1][1][1][0]=d3_Z0[1][1][1][0]; fulld3_z[1][1][1][1]=d3_Z0[1][1][1][1];

		PetscReal dS[8][2];																
		IGAPointFormGrad (pS,S,&dS[0][0]);													//Same for the gradient

		PetscReal fulld_S[3][3][3][3]={0};													//Expand grad(S) to full tensor order form, only non-zero elements
		fulld_S[0][0][0][0]=dS[0][0]; fulld_S[0][0][0][1]=dS[0][1]; 
		fulld_S[0][0][1][0]=dS[1][0]; fulld_S[0][0][1][1]=dS[1][1]; 
		fulld_S[0][1][0][0]=dS[2][0]; fulld_S[0][1][0][1]=dS[2][1]; 
		fulld_S[0][1][1][0]=dS[3][0]; fulld_S[0][1][1][1]=dS[3][1];
		fulld_S[1][0][0][0]=dS[4][0]; fulld_S[1][0][0][1]=dS[4][1]; 
		fulld_S[1][0][1][0]=dS[5][0]; fulld_S[1][0][1][1]=dS[5][1]; 
		fulld_S[1][1][0][0]=dS[6][0]; fulld_S[1][1][0][1]=dS[6][1]; 
		fulld_S[1][1][1][0]=dS[7][0]; fulld_S[1][1][1][1]=dS[7][1];

		PetscReal d2_GradZ[4][2][2];
		IGAPointFormHess (pGradZ,UgradZ,&d2_GradZ[0][0][0]);								//This is the 3-rd order tensor Z_{i,kl} (remember that we are storing Z_{mn} as a column vector Z_{i})

		PetscReal fulld_GradZ[3][3][3][3]={0};
		fulld_GradZ[0][0][0][0]=d2_GradZ[0][0][0]; fulld_GradZ[0][0][0][1]=d2_GradZ[0][0][1]; fulld_GradZ[0][0][1][0]=d2_GradZ[0][1][0]; fulld_GradZ[0][0][1][1]=d2_GradZ[0][1][1];
		fulld_GradZ[0][1][0][0]=d2_GradZ[1][0][0]; fulld_GradZ[0][1][0][1]=d2_GradZ[1][0][1]; fulld_GradZ[0][1][1][0]=d2_GradZ[1][1][0]; fulld_GradZ[0][1][1][1]=d2_GradZ[1][1][1];
		fulld_GradZ[1][0][0][0]=d2_GradZ[2][0][0]; fulld_GradZ[1][0][0][1]=d2_GradZ[2][0][1]; fulld_GradZ[1][0][1][0]=d2_GradZ[2][1][0]; fulld_GradZ[1][0][1][1]=d2_GradZ[2][1][1];
		fulld_GradZ[1][1][0][0]=d2_GradZ[3][0][0]; fulld_GradZ[1][1][0][1]=d2_GradZ[3][0][1]; fulld_GradZ[1][1][1][0]=d2_GradZ[3][1][0]; fulld_GradZ[1][1][1][1]=d2_GradZ[3][1][1];

		PetscReal (*Kstress)[dof][nen][dof] = (typeof(Kstress)) K;
		PetscReal (*Fstress)[dof] = (PetscReal (*)[dof])F;

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

		PetscReal v[3][3]={0};

		if (p->atboundary)
		{
			return 0;
		}
		else
		{
			for (a=0; a<nen; a++) 
			{
				for (b=0; b<nen; b++) 
				{
					for (i=0; i<dof; i++)
					{
						Kstress[a][i][b][i]=N0[a]*N0[b];
					}
				}
			}

			for(a=0 ;a<nen; a++)
			{
				for (i=0; i<dof; i++)
				{
					if (i==0)
					{
						v[0][0]=N0[a]; v[0][1]=0.0; v[1][0]=0.0; v[1][1]=0.0;
						m=0; n=0;
					}
					else if (i==1)
					{
						v[0][0]=0.0; v[0][1]=N0[a]; v[1][0]=0.0; v[1][1]=0.0;
						m=0; n=1;
					}
					else if (i==2)
					{
						v[0][0]=0.0; v[0][1]=0.0; v[1][0]=N0[a]; v[1][1]=0.0;
						m=1; n=0;
					}
					else if (i==3)
					{
						v[0][0]=0.0; v[0][1]=0.0; v[1][0]=0.0; v[1][1]=N0[a];
						m=1; n=1;
					}
					
					Fstress[a][i]=0.0;
					for (k=0;k<3;k++)
					{
						for(l=0;l<3;l++)
						{
							Fstress[a][i]+=0.5*(C[m][n][k][l]*(-fulld_z[k][l]-fullChi[k][l])+C[n][m][k][l]*(-fulld_z[k][l]-fullChi[k][l]))*v[m][n];
						}

						Fstress[a][i]+= -0.25*eps*(-fulld3_z[m][n][k][k]-fulld2_Chi[m][n][k][k]-fulld3_z[m][k][n][k]-fulld2_Chi[m][k][n][k]
												   -fulld3_z[n][m][k][k]-fulld2_Chi[n][m][k][k]-fulld3_z[n][k][m][k]-fulld2_Chi[n][k][m][k])*v[m][n];
						//Next part is the contribution from having S in the energy function
						Fstress[a][i]+= -0.25*eps*(-fulld_S[m][n][k][k]-fulld_S[m][k][n][k]-fulld_S[n][m][k][k]-fulld_S[n][k][m][k])*v[m][n];
						//Next part is the contribution of adding grad(Z) from the energy function, to turn J_hat into J
						Fstress[a][i]+= -0.25*eps*(+fulld_GradZ[m][n][k][k]+fulld_GradZ[m][k][n][k]+fulld_GradZ[n][m][k][k]+fulld_GradZ [n][k][m][k])*v[m][n];

					}
				}
			}
		}
		return 0;
	}
//

//System for L2 projection of classic stress (C*Ue)
	#undef  __FUNCT__
	#define __FUNCT__ "ClassicStress"
	//PetscErrorCode Stress(IGAPoint p,IGAPoint pU, IGAPoint pHs,IGAPoint pChi,IGAPoint pZu,PetscReal *K,PetscReal *F,PetscReal *U,PetscReal *HS, PetscReal *Chi,PetscReal *Zu,void *ctx)		//When other fields like Pi or S (etc.) come into this, declare a IGAPoint pPi or pS and a new PescReal *UPi or *US for each
	PetscErrorCode ClassicStress(IGAPoint p,IGAPoint pChi,IGAPoint pZu,PetscReal *K,PetscReal *F, PetscReal *Chi,PetscReal *Zu,void *ctx)		//When other fields like Pi or S (etc.) come into this, declare a IGAPoint pPi or pS and a new PetscReal *UPi or *US for each
	{
		const PetscReal *N0;
		IGAPointGetShapeFuns(p,0,(const PetscReal**)&N0);									//Value of the shape functions
		PetscInt a,b,i,j,k,l,m,n,nen=p->nen, dof=p->dof;

		//Change to consider G=1
		const PetscReal nu=0.33;
		const PetscReal mu=1.0;
		const PetscReal lambda=2.0*mu*nu/(1.0-2.0*nu);

		PetscReal chi0[4];																	//Array to contain the vector chi(0)
		IGAPointFormValue(pChi,Chi,&chi0[0]);												//Assign chi to its container

		PetscReal d_Z0[2][2];																//Same for its gradient
		IGAPointFormGrad (pZu,Zu,&d_Z0[0][0]);												//Same for the gradient

		//The four non-zero components of Chi are stored as a vector, restore them to an array with the correct indexing for value and derivative
		PetscReal fullChi[3][3]={0};
		fullChi[0][0]=chi0[0]; 	fullChi[0][1]=chi0[1];
		fullChi[1][0]=chi0[2]; 	fullChi[1][1]=chi0[3];

		//Expanding z (and derivatives) to 3 components, more convenient for sums in for loops
		PetscReal fulld_z[3][3]={0};
		fulld_z[0][0]=d_Z0[0][0]; fulld_z[0][1]=d_Z0[0][1];
		fulld_z[1][0]=d_Z0[1][0]; fulld_z[1][1]=d_Z0[1][1];

		PetscReal (*Kstress)[dof][nen][dof] = (typeof(Kstress)) K;
		PetscReal (*Fstress)[dof] = (PetscReal (*)[dof])F;

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

		PetscReal v[3][3]={0};

		if (p->atboundary)
		{
			return 0;
		}
		else
		{
			for (a=0; a<nen; a++) 
			{
				for (b=0; b<nen; b++) 
				{
					for (i=0; i<dof; i++)
					{
						Kstress[a][i][b][i]=N0[a]*N0[b];
					}
				}
			}

			for(a=0 ;a<nen; a++)
			{
				for (i=0; i<dof; i++)
				{
					if (i==0)
					{
						v[0][0]=N0[a]; v[0][1]=0.0; v[1][0]=0.0; v[1][1]=0.0;
						m=0; n=0;
					}
					else if (i==1)
					{
						v[0][0]=0.0; v[0][1]=N0[a]; v[1][0]=0.0; v[1][1]=0.0;
						m=0; n=1;
					}
					else if (i==2)
					{
						v[0][0]=0.0; v[0][1]=0.0; v[1][0]=N0[a]; v[1][1]=0.0;
						m=1; n=0;
					}
					else if (i==3)
					{
						v[0][0]=0.0; v[0][1]=0.0; v[1][0]=0.0; v[1][1]=N0[a];
						m=1; n=1;
					}
					
					Fstress[a][i]=0.0;
					for (k=0;k<3;k++)
					{
						for(l=0;l<3;l++)
						{
							Fstress[a][i]+=0.5*( C[m][n][k][l]*(-fulld_z[k][l]-fullChi[k][l]) + C[n][m][k][l]*(-fulld_z[k][l]-fullChi[k][l]) )*v[m][n];
						}
					}
				}
			}
		}
		return 0;
	}
//

//System for L2 projection of couple-stress
	#undef  __FUNCT__
	#define __FUNCT__ "CoupleStress"
	//			   CoupleStress(pointCS,pointchiUp,pointZ0,pointZS,pointS,KpointCS,FpointCS,Chi0CS,Z0CS,ZSCS,S0CS,NULL)			//When other fields like Pi or S (etc.) come into this, declare a IGAPoint pPi or pS and a new PescReal *UPi or *US for each
	PetscErrorCode CoupleStress(IGAPoint p,IGAPoint pChi,IGAPoint pZu,IGAPoint pZS,IGAPoint pS, PetscReal *K,PetscReal *F, PetscReal *Chi,PetscReal *Zu,PetscReal *ZS,PetscReal *S,void *ctx)		//When other fields like Pi or S (etc.) come into this, declare a IGAPoint pPi or pS and a new PetscReal *UPi or *US for each
	{
		const PetscReal *N0;
		IGAPointGetShapeFuns(p,0,(const PetscReal**)&N0);									//Value of the shape functions
		PetscInt a,b,i,k,l,m,n,nen=p->nen, dof=p->dof;

		//Change to consider G=1
		//const PetscReal nu=0.33;
		const PetscReal mu=1.0;
		//const PetscReal lambda=2.0*mu*nu/(1.0-2.0*nu);
		const PetscReal eps=mu/100.0;															//Choose later based on whatever Amit says :)

		//Definition of alternating tensor
		const PetscReal e[3][3][3]=
		{
			{{0.0,0.0,0.0},{0.0,0.0,1.0},{0.0,-1.0,0.0}},
			{{0.0,0.0,-1.0},{0.0,0.0,0.0},{1.0,0.0,0.0}},
			{{0.0,1.0,0.0},{-1.0,0.0,0.0},{0.0,0.0,0.0}}
		};

		PetscReal d_Chi0[4][2];																//Same for its gradient
		IGAPointFormGrad(pChi,Chi,&d_Chi0[0][0]);											//Assign grad chi to its container
		
		PetscReal d2_Z0[2][2][2];															//Same for its 2nd order partial derivatives
		IGAPointFormHess (pZu,Zu,&d2_Z0[0][0][0]);											//Same for the hessian (second derivatives)

		PetscReal d_ZS[4][2];
		IGAPointFormGrad (pZS,ZS,&d_ZS[0][0]);												//This is the 3-rd order tensor Z_{i,kl} (remember that we are storing Z_{mn} as a column vector Z_{i})

		PetscReal S0[8];																	//Assign S to a vector
		IGAPointFormValue(pS,S,&S0[0]);

		//The four non-zero components of Chi are stored as a vector, restore them to an array with the correct indexing for value and derivative
		PetscReal fulld_Chi[3][3][3]={0};
		fulld_Chi[0][0][0]=d_Chi0[0][0]; fulld_Chi[0][0][1]=d_Chi0[0][1];
		fulld_Chi[0][1][0]=d_Chi0[1][0]; fulld_Chi[0][1][1]=d_Chi0[1][1];
		fulld_Chi[1][0][0]=d_Chi0[2][0]; fulld_Chi[1][0][1]=d_Chi0[2][1];
		fulld_Chi[1][1][0]=d_Chi0[3][0]; fulld_Chi[1][1][1]=d_Chi0[3][1];

		//Expanding z (and derivatives) to 3 components, more convenient for sums in for loops
		PetscReal fulld2_z[3][3][3]={0};
		fulld2_z[0][0][0]=d2_Z0[0][0][0]; fulld2_z[0][0][1]=d2_Z0[0][0][1];
		fulld2_z[0][1][0]=d2_Z0[0][1][0]; fulld2_z[0][1][1]=d2_Z0[0][1][1];
		fulld2_z[1][0][0]=d2_Z0[1][0][0]; fulld2_z[1][0][1]=d2_Z0[1][0][1];
		fulld2_z[1][1][0]=d2_Z0[1][1][0]; fulld2_z[1][1][1]=d2_Z0[1][1][1];

		PetscReal fulld_ZS[3][3][3]={0};
		fulld_ZS[0][0][0]=d_ZS[0][0]; fulld_ZS[0][0][1]=d_ZS[0][1];
		fulld_ZS[0][1][0]=d_ZS[1][0]; fulld_ZS[0][1][1]=d_ZS[1][1];
		fulld_ZS[1][0][0]=d_ZS[2][0]; fulld_ZS[1][0][1]=d_ZS[2][1];
		fulld_ZS[1][1][0]=d_ZS[3][0]; fulld_ZS[1][1][1]=d_ZS[3][1];

		//Inflate stored vectors to full tensor form
		PetscReal fullS[3][3][3]={0};
		fullS[0][0][0]=S0[0]; fullS[0][0][1]=S0[1];											//Expand S to full 3rd order form, only non-zero elements
		fullS[0][1][0]=S0[2]; fullS[0][1][1]=S0[3];
		fullS[1][0][0]=S0[4]; fullS[1][0][1]=S0[5];
		fullS[1][1][0]=S0[6]; fullS[1][1][1]=S0[7];

		PetscReal (*KCS)[dof][nen][dof] = (typeof(KCS)) K;
		PetscReal (*FCS)[dof] = (PetscReal (*)[dof])F;

		PetscReal v[3][3]={0};

		if (p->atboundary)
		{
			return 0;
		}
		else
		{
			for (a=0; a<nen; a++) 
			{
				for (b=0; b<nen; b++) 
				{
					for (i=0; i<dof; i++)
					{
						KCS[a][i][b][i]=N0[a]*N0[b];
					}
				}
			}

			for(a=0 ;a<nen; a++)
			{
				for (i=0; i<dof; i++)
				{
					if (i==0)
					{
						v[2][0]=N0[a]; v[2][1]=0.0; v[2][2]=0.0;
						m=2; n=0;
					}
					else if (i==1)
					{
						v[2][0]=0.0; v[2][1]=N0[a]; v[2][2]=0.0;
						m=2; n=1;
					}
					
					FCS[a][i]=0.0;
					for (k=0;k<3;k++)
					{
						for(l=0;l<3;l++)
						{
							FCS[a][i]+=0.5*eps*e[m][k][l]*(fullS[k][l][n]+fulld2_z[k][l][n]+fulld_Chi[k][l][n]+fulld_ZS[k][l][n]
														  +fullS[k][n][l]+fulld2_z[k][n][l]+fulld_Chi[k][n][l]+fulld_ZS[k][n][l])*v[m][n];
						}
					}
				}
			}
		}
		return 0;
	}
//

//System for L2 projection of energy density  //Corrected for J instead of J_hat
	#undef  __FUNCT__
	#define __FUNCT__ "EnergyDensity"
	//PetscErrorCode Stress(IGAPoint p,IGAPoint pU, IGAPoint pHs,IGAPoint pChi,IGAPoint pZu,PetscReal *K,PetscReal *F,PetscReal *U,PetscReal *HS, PetscReal *Chi,PetscReal *Zu,void *ctx)		//When other fields like Pi or S (etc.) come into this, declare a IGAPoint pPi or pS and a new PescReal *UPi or *US for each
	PetscErrorCode EnergyDensity(IGAPoint p, IGAPoint pChi, IGAPoint pZu, IGAPoint pS, IGAPoint pZS, PetscReal *K, PetscReal *F, PetscReal *Chi, PetscReal *Zu, PetscReal *S, PetscReal *ZS, void *ctx)		//When other fields like Pi or S (etc.) come into this, declare a IGAPoint pPi or pS and a new PestcReal *UPi or *US for each
	{
		const PetscReal *N0;
		IGAPointGetShapeFuns(p,0,(const PetscReal**)&N0);									//Value of the shape functions
		PetscInt a,b,i,j,k,l,m,n,nen=p->nen, dof=p->dof;

		PetscReal x[2];																		//Vector of reals, size equal to problem's dimension
		IGAPointFormGeomMap(p,x);

		//Change to consider G=1
		const PetscReal nu=0.33;
		const PetscReal mu=1.0;
		const PetscReal lambda=2.0*mu*nu/(1.0-2.0*nu);
		const PetscReal eps=mu/100.0;															//Choose later based on whatever Amit says :)

		PetscReal C[3][3][3][3]={0};
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

		PetscReal chi0[4];																	//Array to contain the vector chi(0)
		PetscReal d_Chi0[4][2];																//Same for its gradient
		IGAPointFormValue(pChi,Chi,&chi0[0]);												//Assign chi to its container
		IGAPointFormGrad(pChi,Chi,&d_Chi0[0][0]);											//Assign grad chi to its container

		PetscReal d_Z0[2][2];																//Same for its gradient
		PetscReal d2_Z0[2][2][2];															//Same for its 2nd order partial derivatives
		IGAPointFormGrad (pZu,Zu,&d_Z0[0][0]);												//Same for the gradient
		IGAPointFormHess (pZu,Zu,&d2_Z0[0][0][0]);											//Same for the hessian (second derivatives)

		PetscReal S0[8];																	//Array to contain the vector S0
		IGAPointFormValue(pS,S,&S0[0]);														//Assign chi to its container

		PetscReal d_ZS[4][2];
		IGAPointFormGrad (pZS,ZS,&d_ZS[0][0]);												//This is the 3-rd order tensor Z_{i,kl} (remember that we are storing Z_{mn} as a column vector Z_{i})

		//The four non-zero components of Chi are stored as a vector, restore them to an array with the correct indexing for value and derivative
		PetscReal fullChi[3][3]={0};
		fullChi[0][0]=chi0[0]; 	fullChi[0][1]=chi0[1];
		fullChi[1][0]=chi0[2]; 	fullChi[1][1]=chi0[3];

		PetscReal fulld_Chi[3][3][3]={0};
		fulld_Chi[0][0][0]=d_Chi0[0][0]; fulld_Chi[0][0][1]=d_Chi0[0][1];
		fulld_Chi[0][1][0]=d_Chi0[1][0]; fulld_Chi[0][1][1]=d_Chi0[1][1];
		fulld_Chi[1][0][0]=d_Chi0[2][0]; fulld_Chi[1][0][1]=d_Chi0[2][1];
		fulld_Chi[1][1][0]=d_Chi0[3][0]; fulld_Chi[1][1][1]=d_Chi0[3][1];

		//Expanding z (and derivatives) to 3 components, more convenient for sums in for loops
		PetscReal fulld_z[3][3]={0};
		fulld_z[0][0]=d_Z0[0][0]; fulld_z[0][1]=d_Z0[0][1];
		fulld_z[1][0]=d_Z0[1][0]; fulld_z[1][1]=d_Z0[1][1];

		PetscReal fulld2_z[3][3][3]={0};
		fulld2_z[0][0][0]=d2_Z0[0][0][0]; fulld2_z[0][0][1]=d2_Z0[0][0][1];
		fulld2_z[0][1][0]=d2_Z0[0][1][0]; fulld2_z[0][1][1]=d2_Z0[0][1][1];
		fulld2_z[1][0][0]=d2_Z0[1][0][0]; fulld2_z[1][0][1]=d2_Z0[1][0][1];
		fulld2_z[1][1][0]=d2_Z0[1][1][0]; fulld2_z[1][1][1]=d2_Z0[1][1][1];

		PetscReal fullS[3][3][3]={0};
		fullS[0][0][0]=S0[0]; fullS[0][0][1]=S0[1];											//Expand S to full 3rd order form, only non-zero elements
		fullS[0][1][0]=S0[2]; fullS[0][1][1]=S0[3];
		fullS[1][0][0]=S0[4]; fullS[1][0][1]=S0[5];
		fullS[1][1][0]=S0[6]; fullS[1][1][1]=S0[7];

		PetscReal fulld_ZS[3][3][3]={0};
		fulld_ZS[0][0][0]=d_ZS[0][0]; fulld_ZS[0][0][1]=d_ZS[0][1];
		fulld_ZS[0][1][0]=d_ZS[1][0]; fulld_ZS[0][1][1]=d_ZS[1][1];
		fulld_ZS[1][0][0]=d_ZS[2][0]; fulld_ZS[1][0][1]=d_ZS[2][1];
		fulld_ZS[1][1][0]=d_ZS[3][0]; fulld_ZS[1][1][1]=d_ZS[3][1];

		PetscReal (*KED)[dof][nen][dof] = (typeof(KED)) K;
		PetscReal (*FED)[dof] = (PetscReal (*)[dof])F;

		if (p->atboundary)
		{
			return 0;
		}
		else
		{
			for (a=0; a<nen; a++) 
			{
				for (b=0; b<nen; b++) 
				{
					for (i=0; i<dof; i++)
					{
						KED[a][i][b][i]=N0[a]*N0[b];
					}
				}
			}

			for(a=0 ;a<nen; a++)
			{
				for (i=0; i<dof; i++)
				{
					FED[a][i]=0.0;
					//if(x[0]>-3.0 && x[0]<3.0 && x[1]>-3.0 && x[1]<3.0)
					//{
					//	FED[a][i]=0.0;
					//}
					//else
					//{
						for (k=0;k<3;k++)
					{
						for(l=0;l<3;l++)
						{
							for (m=0;m<3;m++)
							{
								for (n=0;n<3;n++)
								{
									FED[a][i]+=0.5*((-fulld_z[k][l]-fullChi[k][l])*C[k][l][m][n]*(-fulld_z[m][n]-fullChi[m][n]))*N0[a];
								}
								FED[a][i]+=0.5*eps*(fullS[k][l][m]+fulld2_z[k][l][m]+fulld_Chi[k][l][m]+fulld_ZS[k][l][m])*(fullS[k][l][m]+fulld2_z[k][l][m]+fulld_Chi[k][l][m]+fulld_ZS[k][l][m])*N0[a];
							}
						}
					}
					//}
				}
			}
		}
		return 0;
	}
//

//System for L2 projection of full stress
	#undef  __FUNCT__
	#define __FUNCT__ "FullS"
	PetscErrorCode FullS(IGAPoint p,IGAPoint pChi,IGAPoint pZu, PetscReal *K,PetscReal *F, PetscReal *Chi,PetscReal *Zu, void *ctx)		//When other fields like Pi or S (etc.) come into this, declare a IGAPoint pPi or pS and a new PestcReal *UPi or *US for each
	{
		const PetscReal *N0;
		IGAPointGetShapeFuns(p,0,(const PetscReal**)&N0);									//Value of the shape functions
		PetscInt a,b,i,j,k,l,m,n,nen=p->nen, dof=p->dof;

		//Change to consider G=1
		const PetscReal nu=0.33;
		const PetscReal mu=1.0;
		const PetscReal lambda=2.0*mu*nu/(1.0-2.0*nu);
		const PetscReal eps=mu/100.0;															//Choose later based on whatever Amit says :)

		PetscReal chi0[4];																	//Array to contain the vector chi(0)
		PetscReal d2_Chi0[4][2][2];															//Same for its Hessian
		IGAPointFormValue(pChi,Chi,&chi0[0]);												//Assign chi to its container
		IGAPointFormHess (pChi,Chi,&d2_Chi0[0][0][0]);										//This should be the 3-rd order tensor Chi_{i,jk} (remember that we are storing Chi_{kl} as a column vector)

		PetscReal d_Z0[2][2];																//Same for its gradient
		PetscReal d3_Z0[2][2][2][2];														//Same for its 3rd order partial derivatives
		IGAPointFormGrad (pZu,Zu,&d_Z0[0][0]);												//Same for the gradient
		IGAPointFormDer3 (pZu,Zu,&d3_Z0[0][0][0][0]);										//Same for the thir derivatives

		//The four non-zero components of Chi are stored as a vector, restore them to an array with the correct indexing for value and derivative
		PetscReal fullChi[3][3]={0};
		fullChi[0][0]=chi0[0]; 	fullChi[0][1]=chi0[1];
		fullChi[1][0]=chi0[2]; 	fullChi[1][1]=chi0[3];

		PetscReal fulld2_Chi[3][3][3][3]={0};
		fulld2_Chi[0][0][0][0]=d2_Chi0[0][0][0]; fulld2_Chi[0][0][0][1]=d2_Chi0[0][0][1]; fulld2_Chi[0][0][1][0]=d2_Chi0[0][1][0]; fulld2_Chi[0][0][1][1]=d2_Chi0[0][1][1];
		fulld2_Chi[0][1][0][0]=d2_Chi0[1][0][0]; fulld2_Chi[0][1][0][1]=d2_Chi0[1][0][1]; fulld2_Chi[0][1][1][0]=d2_Chi0[1][1][0]; fulld2_Chi[0][1][1][1]=d2_Chi0[1][1][1];
		fulld2_Chi[1][0][0][0]=d2_Chi0[2][0][0]; fulld2_Chi[1][0][0][1]=d2_Chi0[2][0][1]; fulld2_Chi[1][0][1][0]=d2_Chi0[2][1][0]; fulld2_Chi[1][0][1][1]=d2_Chi0[2][1][1];
		fulld2_Chi[1][1][0][0]=d2_Chi0[3][0][0]; fulld2_Chi[1][1][0][1]=d2_Chi0[3][0][1]; fulld2_Chi[1][1][1][0]=d2_Chi0[3][1][0]; fulld2_Chi[1][1][1][1]=d2_Chi0[3][1][1];

		//Expanding z (and derivatives) to 3 components, more convenient for sums in for loops
		PetscReal fulld_z[3][3]={0};
		fulld_z[0][0]=d_Z0[0][0]; fulld_z[0][1]=d_Z0[0][1];
		fulld_z[1][0]=d_Z0[1][0]; fulld_z[1][1]=d_Z0[1][1];

		PetscReal fulld3_z[3][3][3][3]={0};
		fulld3_z[0][0][0][0]=d3_Z0[0][0][0][0]; fulld3_z[0][0][0][1]=d3_Z0[0][0][0][1]; fulld3_z[0][0][1][0]=d3_Z0[0][0][1][0]; fulld3_z[0][0][1][1]=d3_Z0[0][0][1][1];
		fulld3_z[0][1][0][0]=d3_Z0[0][1][0][0]; fulld3_z[0][1][0][1]=d3_Z0[0][1][0][1]; fulld3_z[0][1][1][0]=d3_Z0[0][1][1][0]; fulld3_z[0][1][1][1]=d3_Z0[0][1][1][1];
		
		fulld3_z[1][0][0][0]=d3_Z0[1][0][0][0]; fulld3_z[1][0][0][1]=d3_Z0[1][0][0][1]; fulld3_z[1][0][1][0]=d3_Z0[1][0][1][0]; fulld3_z[1][0][1][1]=d3_Z0[1][0][1][1]; 
		fulld3_z[1][1][0][0]=d3_Z0[1][1][0][0]; fulld3_z[1][1][0][1]=d3_Z0[1][1][0][1]; fulld3_z[1][1][1][0]=d3_Z0[1][1][1][0]; fulld3_z[1][1][1][1]=d3_Z0[1][1][1][1];

		PetscReal (*Kstress)[dof][nen][dof] = (typeof(Kstress)) K;
		PetscReal (*Fstress)[dof] = (PetscReal (*)[dof])F;

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

		PetscReal v[3][3]={0};

		if (p->atboundary)
		{
			return 0;
		}
		else
		{
			for (a=0; a<nen; a++) 
			{
				for (b=0; b<nen; b++) 
				{
					for (i=0; i<dof; i++)
					{
						Kstress[a][i][b][i]=N0[a]*N0[b];
					}
				}
			}

			for(a=0 ;a<nen; a++)
			{
				for (i=0; i<dof; i++)
				{
					if (i==0)
					{
						v[0][0]=N0[a]; v[0][1]=0.0;
						v[1][0]=0.0;   v[1][1]=0.0;
						m=0; n=0;
					}
					else if (i==1)
					{
						v[0][0]=0.0; v[0][1]=N0[a];
						v[1][0]=0.0; v[1][1]=0.0;
						m=0; n=1;
					}
					else if (i==2)
					{
						v[0][0]=0.0;   v[0][1]=0.0;
						v[1][0]=N0[a]; v[1][1]=0.0;
						m=1; n=0;
					}
					else if (i==3)
					{
						v[0][0]=0.0; v[0][1]=0.0;
						v[1][0]=0.0; v[1][1]=N0[a];
						m=1; n=1;
					}
					
					Fstress[a][i]=0.0;
					for (k=0;k<3;k++)
					{
						for(l=0;l<3;l++)
						{
							Fstress[a][i]+=C[m][n][k][l]*(-fulld_z[k][l]-fullChi[k][l])*v[m][n];
						}

						Fstress[a][i]+= -0.5*eps*(-fulld3_z[m][n][k][k]-fulld2_Chi[m][n][k][k]-fulld3_z[m][k][n][k]-fulld2_Chi[m][k][n][k])*v[m][n];
					}
				}
			}
		}
		return 0;
	}
//

//System for L2 projection of skew stress
	#undef  __FUNCT__
	#define __FUNCT__ "SkewS"
	PetscErrorCode SkewS(IGAPoint p,IGAPoint pChi,IGAPoint pZu,IGAPoint pZS, IGAPoint pS, PetscReal *K,PetscReal *F, PetscReal *Chi,PetscReal *Zu,PetscReal *ZS,PetscReal *S, void *ctx)		//When other fields like Pi or S (etc.) come into this, declare a IGAPoint pPi or pS and a new PestcReal *UPi or *US for each
	{
		const PetscReal *N0;
		IGAPointGetShapeFuns(p,0,(const PetscReal**)&N0);									//Value of the shape functions
		PetscInt a,b,i,k,m,n,nen=p->nen, dof=p->dof;

		//Change to consider G=1
		const PetscReal mu=1.0;
		const PetscReal eps=mu/100.0;														//Choose later based on whatever Amit says :)

		PetscReal d2_Chi0[4][2][2];															//Same for its Hessian
		IGAPointFormHess (pChi,Chi,&d2_Chi0[0][0][0]);										//This should be the 3-rd order tensor Chi_{i,jk} (remember that we are storing Chi_{kl} as a column vector)

		PetscReal d3_Z0[2][2][2][2];														//Same for its 3rd order partial derivatives
		IGAPointFormDer3 (pZu,Zu,&d3_Z0[0][0][0][0]);										//Same for the thir derivatives

		PetscReal dS[8][2];																
		IGAPointFormGrad (pS,S,&dS[0][0]);													//Same for the gradient

		PetscReal d2_ZS[4][2][2];
		IGAPointFormHess (pZS,ZS,&d2_ZS[0][0][0]);											//This is the 3-rd order tensor Z_{i,kl} (remember that we are storing Z_{mn} as a column vector Z_{i})

		//The four non-zero components of Chi are stored as a vector, restore them to an array with the correct indexing for value and derivative
		PetscReal fulld2_Chi[3][3][3][3]={0};
		fulld2_Chi[0][0][0][0]=d2_Chi0[0][0][0]; fulld2_Chi[0][0][0][1]=d2_Chi0[0][0][1]; fulld2_Chi[0][0][1][0]=d2_Chi0[0][1][0]; fulld2_Chi[0][0][1][1]=d2_Chi0[0][1][1];
		fulld2_Chi[0][1][0][0]=d2_Chi0[1][0][0]; fulld2_Chi[0][1][0][1]=d2_Chi0[1][0][1]; fulld2_Chi[0][1][1][0]=d2_Chi0[1][1][0]; fulld2_Chi[0][1][1][1]=d2_Chi0[1][1][1];
		fulld2_Chi[1][0][0][0]=d2_Chi0[2][0][0]; fulld2_Chi[1][0][0][1]=d2_Chi0[2][0][1]; fulld2_Chi[1][0][1][0]=d2_Chi0[2][1][0]; fulld2_Chi[1][0][1][1]=d2_Chi0[2][1][1];
		fulld2_Chi[1][1][0][0]=d2_Chi0[3][0][0]; fulld2_Chi[1][1][0][1]=d2_Chi0[3][0][1]; fulld2_Chi[1][1][1][0]=d2_Chi0[3][1][0]; fulld2_Chi[1][1][1][1]=d2_Chi0[3][1][1];

		//Expanding z (and derivatives) to 3 components, more convenient for sums in for loops
		PetscReal fulld3_z[3][3][3][3]={0};
		fulld3_z[0][0][0][0]=d3_Z0[0][0][0][0]; fulld3_z[0][0][0][1]=d3_Z0[0][0][0][1]; fulld3_z[0][0][1][0]=d3_Z0[0][0][1][0]; fulld3_z[0][0][1][1]=d3_Z0[0][0][1][1];
		fulld3_z[0][1][0][0]=d3_Z0[0][1][0][0]; fulld3_z[0][1][0][1]=d3_Z0[0][1][0][1]; fulld3_z[0][1][1][0]=d3_Z0[0][1][1][0]; fulld3_z[0][1][1][1]=d3_Z0[0][1][1][1];
		
		fulld3_z[1][0][0][0]=d3_Z0[1][0][0][0]; fulld3_z[1][0][0][1]=d3_Z0[1][0][0][1]; fulld3_z[1][0][1][0]=d3_Z0[1][0][1][0]; fulld3_z[1][0][1][1]=d3_Z0[1][0][1][1]; 
		fulld3_z[1][1][0][0]=d3_Z0[1][1][0][0]; fulld3_z[1][1][0][1]=d3_Z0[1][1][0][1]; fulld3_z[1][1][1][0]=d3_Z0[1][1][1][0]; fulld3_z[1][1][1][1]=d3_Z0[1][1][1][1];

		//Expand grad(S) to full tensor order form, only non-zero elements
		PetscReal fulld_S[3][3][3][3]={0};
		fulld_S[0][0][0][0]=dS[0][0]; fulld_S[0][0][0][1]=dS[0][1]; 
		fulld_S[0][0][1][0]=dS[1][0]; fulld_S[0][0][1][1]=dS[1][1]; 
		fulld_S[0][1][0][0]=dS[2][0]; fulld_S[0][1][0][1]=dS[2][1]; 
		fulld_S[0][1][1][0]=dS[3][0]; fulld_S[0][1][1][1]=dS[3][1];
		fulld_S[1][0][0][0]=dS[4][0]; fulld_S[1][0][0][1]=dS[4][1]; 
		fulld_S[1][0][1][0]=dS[5][0]; fulld_S[1][0][1][1]=dS[5][1]; 
		fulld_S[1][1][0][0]=dS[6][0]; fulld_S[1][1][0][1]=dS[6][1]; 
		fulld_S[1][1][1][0]=dS[7][0]; fulld_S[1][1][1][1]=dS[7][1];

		PetscReal fulld_ZS[3][3][3][3]={0};
		fulld_ZS[0][0][0][0]=d2_ZS[0][0][0]; fulld_ZS[0][0][0][1]=d2_ZS[0][0][1]; fulld_ZS[0][0][1][0]=d2_ZS[0][1][0]; fulld_ZS[0][0][1][1]=d2_ZS[0][1][1];
		fulld_ZS[0][1][0][0]=d2_ZS[1][0][0]; fulld_ZS[0][1][0][1]=d2_ZS[1][0][1]; fulld_ZS[0][1][1][0]=d2_ZS[1][1][0]; fulld_ZS[0][1][1][1]=d2_ZS[1][1][1];
		fulld_ZS[1][0][0][0]=d2_ZS[2][0][0]; fulld_ZS[1][0][0][1]=d2_ZS[2][0][1]; fulld_ZS[1][0][1][0]=d2_ZS[2][1][0]; fulld_ZS[1][0][1][1]=d2_ZS[2][1][1];
		fulld_ZS[1][1][0][0]=d2_ZS[3][0][0]; fulld_ZS[1][1][0][1]=d2_ZS[3][0][1]; fulld_ZS[1][1][1][0]=d2_ZS[3][1][0]; fulld_ZS[1][1][1][1]=d2_ZS[3][1][1];

		PetscReal (*Kskewstress)[dof][nen][dof] = (typeof(Kskewstress)) K;
		PetscReal (*Fskewstress)[dof] = (PetscReal (*)[dof])F;

		PetscReal v[3][3]={0};

		if (p->atboundary)
		{
			return 0;
		}
		else
		{
			for (a=0; a<nen; a++) 
			{
				for (b=0; b<nen; b++) 
				{
					for (i=0; i<dof; i++)
					{
						Kskewstress[a][i][b][i]=N0[a]*N0[b];
					}
				}
			}

			for(a=0 ;a<nen; a++)
			{
				for (i=0; i<dof; i++)
				{
					if (i==0)
					{
						v[0][0]=N0[a]; v[0][1]=0.0;
						v[1][0]=0.0;   v[1][1]=0.0;
						m=0; n=0;
					}
					else if (i==1)
					{
						v[0][0]=0.0; v[0][1]=N0[a];
						v[1][0]=0.0; v[1][1]=0.0;
						m=0; n=1;
					}
					else if (i==2)
					{
						v[0][0]=0.0;   v[0][1]=0.0;
						v[1][0]=N0[a]; v[1][1]=0.0;
						m=1; n=0;
					}
					else if (i==3)
					{
						v[0][0]=0.0; v[0][1]=0.0;
						v[1][0]=0.0; v[1][1]=N0[a];
						m=1; n=1;
					}
					
					Fskewstress[a][i]=0.0;
					for (k=0;k<3;k++)
					{
						Fskewstress[a][i]+= 0.5*eps*(fulld_S[m][n][k][k]+fulld3_z[m][n][k][k]+fulld2_Chi[m][n][k][k]+fulld_ZS[m][n][k][k]
													+fulld_S[m][k][n][k]+fulld3_z[m][k][n][k]+fulld2_Chi[m][k][n][k]+fulld_ZS[m][k][n][k]
													-fulld_S[n][m][k][k]+fulld3_z[n][m][k][k]+fulld2_Chi[n][m][k][k]+fulld_ZS[n][m][k][k]
													-fulld_S[n][k][m][k]+fulld3_z[n][k][m][k]+fulld2_Chi[n][k][m][k]+fulld_ZS[n][k][m][k])*v[m][n];
					}
				}
			}
		}
		return 0;
	}
//

//System for L2 projection of V^{alpha}
	#undef  __FUNCT__
	#define __FUNCT__ "Valpha"
	//PetscErrorCode Stress(IGAPoint p,IGAPoint pU, IGAPoint pHs,IGAPoint pChi,IGAPoint pZu,PetscReal *K,PetscReal *F,PetscReal *U,PetscReal *HS, PetscReal *Chi,PetscReal *Zu,void *ctx)		//When other fields like Pi or S (etc.) come into this, declare a IGAPoint pPi or pS and a new PescReal *UPi or *US for each
	PetscErrorCode Valpha(IGAPoint p,IGAPoint pChi,IGAPoint pZu,IGAPoint pAl,PetscReal *K,PetscReal *F, PetscReal *Chi,PetscReal *Zu,PetscReal *U,void *ctx)		//When other fields like Pi or S (etc.) come into this, declare a IGAPoint pPi or pS and a new PetscReal *UPi or *US for each
	{
		const PetscReal *N0;
		IGAPointGetShapeFuns(p,0,(const PetscReal**)&N0);									//Value of the shape functions
		PetscInt a,b,c,d,i,j,k,l,m,nen=p->nen, dof=p->dof;

		//Change to consider G=1
		const PetscReal nu=0.33;
		const PetscReal mu=1.0;
		const PetscReal lambda=2.0*mu*nu/(1.0-2.0*nu);
		const PetscReal eps=mu/100.0;														//Choose later based on whatever Amit says :)

		PetscReal x[2];																		//Vector of reals, size equal to problem's dimension
		IGAPointFormGeomMap(p,x);															//Fills x with the coordinates of p, Gauss's point

		//Definition of alternating tensor
		const PetscReal e[3][3][3]=
		{
			{{0.0,0.0,0.0},{0.0,0.0,1.0},{0.0,-1.0,0.0}},
			{{0.0,0.0,-1.0},{0.0,0.0,0.0},{1.0,0.0,0.0}},
			{{0.0,1.0,0.0},{-1.0,0.0,0.0},{0.0,0.0,0.0}}
		};

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

		PetscReal alfa[2];																	//Create array to recieve Alfa
		IGAPointFormValue(pAl,U,&alfa[0]);													//This fills the values

		PetscReal chi0[4];																	//Array to contain the vector chi(0)
		PetscReal d2_Chi0[4][2][2];															//Same for its Hessian
		IGAPointFormValue(pChi,Chi,&chi0[0]);												//Assign chi to its container
		IGAPointFormHess (pChi,Chi,&d2_Chi0[0][0][0]);										//This should be the 3-rd order tensor Chi_{i,jk} (remember that we are storing Chi_{kl} as a column vector)

		PetscReal d_Z0[2][2];																//Same for its gradient
		PetscReal d3_Z0[2][2][2][2];														//Same for its 3rd order partial derivatives
		IGAPointFormGrad (pZu,Zu,&d_Z0[0][0]);												//Same for the gradient
		IGAPointFormDer3 (pZu,Zu,&d3_Z0[0][0][0][0]);										//Same for the thir derivatives

		//Inflate stored vectors to full tensor form
		PetscReal fullAlfa[3][3]={0};
		fullAlfa[0][2]=alfa[0]; fullAlfa[1][2]=alfa[1];										//Expand Alfa to full 2nd order form, only non-zero elements

		//The four non-zero components of Chi are stored as a vector, restore them to an array with the correct indexing for value and derivative
		PetscReal fullChi[3][3]={0};
		fullChi[0][0]=chi0[0]; 	fullChi[0][1]=chi0[1];
		fullChi[1][0]=chi0[2]; 	fullChi[1][1]=chi0[3];

		PetscReal fulld2_Chi[3][3][3][3]={0};
		fulld2_Chi[0][0][0][0]=d2_Chi0[0][0][0]; fulld2_Chi[0][0][0][1]=d2_Chi0[0][0][1]; fulld2_Chi[0][0][1][0]=d2_Chi0[0][1][0]; fulld2_Chi[0][0][1][1]=d2_Chi0[0][1][1];
		fulld2_Chi[0][1][0][0]=d2_Chi0[1][0][0]; fulld2_Chi[0][1][0][1]=d2_Chi0[1][0][1]; fulld2_Chi[0][1][1][0]=d2_Chi0[1][1][0]; fulld2_Chi[0][1][1][1]=d2_Chi0[1][1][1];
		fulld2_Chi[1][0][0][0]=d2_Chi0[2][0][0]; fulld2_Chi[1][0][0][1]=d2_Chi0[2][0][1]; fulld2_Chi[1][0][1][0]=d2_Chi0[2][1][0]; fulld2_Chi[1][0][1][1]=d2_Chi0[2][1][1];
		fulld2_Chi[1][1][0][0]=d2_Chi0[3][0][0]; fulld2_Chi[1][1][0][1]=d2_Chi0[3][0][1]; fulld2_Chi[1][1][1][0]=d2_Chi0[3][1][0]; fulld2_Chi[1][1][1][1]=d2_Chi0[3][1][1];

		//Expanding z (and derivatives) to 3 components, more convenient for sums in for loops
		PetscReal fulld_z[3][3]={0};
		fulld_z[0][0]=d_Z0[0][0]; fulld_z[0][1]=d_Z0[0][1];
		fulld_z[1][0]=d_Z0[1][0]; fulld_z[1][1]=d_Z0[1][1];

		PetscReal fulld3_z[3][3][3][3]={0};
		fulld3_z[0][0][0][0]=d3_Z0[0][0][0][0]; fulld3_z[0][0][0][1]=d3_Z0[0][0][0][1]; fulld3_z[0][0][1][0]=d3_Z0[0][0][1][0]; fulld3_z[0][0][1][1]=d3_Z0[0][0][1][1];
		fulld3_z[0][1][0][0]=d3_Z0[0][1][0][0]; fulld3_z[0][1][0][1]=d3_Z0[0][1][0][1]; fulld3_z[0][1][1][0]=d3_Z0[0][1][1][0]; fulld3_z[0][1][1][1]=d3_Z0[0][1][1][1];
		fulld3_z[1][0][0][0]=d3_Z0[1][0][0][0]; fulld3_z[1][0][0][1]=d3_Z0[1][0][0][1]; fulld3_z[1][0][1][0]=d3_Z0[1][0][1][0]; fulld3_z[1][0][1][1]=d3_Z0[1][0][1][1]; 
		fulld3_z[1][1][0][0]=d3_Z0[1][1][0][0]; fulld3_z[1][1][0][1]=d3_Z0[1][1][0][1]; fulld3_z[1][1][1][0]=d3_Z0[1][1][1][0]; fulld3_z[1][1][1][1]=d3_Z0[1][1][1][1];

		PetscReal (*KCS)[dof][nen][dof] = (typeof(KCS)) K;
		PetscReal (*FCS)[dof] = (PetscReal (*)[dof])F;

		PetscReal v[3]={0};

		if (p->atboundary)
		{
			return 0;
		}
		else
		{
			for (a=0; a<nen; a++) 
			{
				for (b=0; b<nen; b++) 
				{
					for (i=0; i<dof; i++)
					{
						KCS[a][i][b][i]=N0[a]*N0[b];
					}
				}
			}

			for(a=0 ;a<nen; a++)
			{
				for (i=0; i<dof; i++)
				{
					if (i==0)
					{
						v[0]=N0[a]; v[1]=0.0; v[2]=0.0;
					}
					else if (i==1)
					{
						v[0]=0.0; v[1]=N0[a]; v[2]=0.0;
					}
					
					FCS[a][i]=0.0;

					for (j=0;j<3;j++)
					{
						for (k=0;k<3;k++)
						{
							for (l=0;l<3;l++)
							{
								for (m=0;m<3;m++)
								{
									for (c=0;c<3;c++)
									{
										for (d=0;d<3;d++)
										{
											FCS[a][i]+=C[j][k][l][m]*(-fulld_z[l][m]-fullChi[l][m])*e[k][c][d]*fullAlfa[j][c]*v[d];
										}
									}
								}
							}
						}
					}
					for (j=0;j<3;j++)
					{
						for (k=0;k<3;k++)
						{
							for (l=0;l<3;l++)
							{
								for (m=0;m<3;m++)
								{
									for (c=0;c<3;c++)
									{
										FCS[a][i]+=-eps*(-fulld3_z[j][k][l][l]-fulld2_Chi[j][k][l][l])*e[k][m][c]*fullAlfa[j][m]*v[c];
									}
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

//System for L2 projection of V^{alpha} smooth (Integrates Va multiplied with xi(alpha), where xi(alpha)=1 if any component of alpha is different than 0, and 0 in other cases)
	#undef  __FUNCT__
	#define __FUNCT__ "IntValpha"
	//IntValpha(pointVa,pointAlp,FpointVaInt,PointInt1,PointInt2,Va0Values,Al1Va,NULL);CHKERRQ(ierr);		//When other fields like Pi or S (etc.) come into this, declare a IGAPoint pPi or pS and a new PescReal *UPi or *US for each
	PetscErrorCode IntValpha(IGAPoint pV,IGAPoint pAl1,PetscReal *FInt1a,PetscReal *FInt2a,PetscReal *FInt1b,PetscReal *FInt2b,PetscReal *Valpha,PetscReal *Al1,PetscReal *Al2,void *ctx)		//When other fields like Pi or S (etc.) come into this, declare a IGAPoint pPi or pS and a new PetscReal *UPi or *US for each
	{
		PetscInt dof=pV->dof;

		PetscReal Va[2];																	//Create array to recieve Alfa
		IGAPointFormValue(pV,Valpha,&Va[0]);												//This fills the values

		PetscReal alfa1[2];																	//Create array to recieve Alfa
		IGAPointFormValue(pAl1,Al1,&alfa1[0]);												//This fills the values

		PetscReal alfa2[2];																	//Create array to recieve Alfa
		IGAPointFormValue(pAl1,Al2,&alfa2[0]);												//This fills the values

		PetscReal xi1=0.0,xi2=0.0;
		if (alfa1[0]>0 || alfa1[1]>0 || alfa1[0]<0 || alfa1[1]<0)
		{
			xi1=1.0;
		}

		if (alfa2[0]>0 || alfa2[1]>0 || alfa2[0]<0 || alfa2[1]<0)
		{
			xi2=1.0;
		}

		PetscReal (*FI1a)[dof] = (PetscReal (*)[dof])FInt1a;
		PetscReal (*FI2a)[dof] = (PetscReal (*)[dof])FInt2a;
		PetscReal (*FI1b)[dof] = (PetscReal (*)[dof])FInt1b;
		PetscReal (*FI2b)[dof] = (PetscReal (*)[dof])FInt2b;

		if (pV->atboundary)
		{
			return 0;
		}
		else
		{
			FI1a[0][0]+=xi1*Va[0];
			FI2a[0][1]+=xi1*Va[1];
			FI1b[0][0]+=xi2*Va[0];
			FI2b[0][1]+=xi2*Va[1];

		}
		return 0;
	}

	#undef  __FUNCT__
	#define __FUNCT__ "ProjValpha"
	PetscErrorCode ProjValpha(IGAPoint pV,IGAPoint pAl1,PetscReal *FVa1,PetscReal *FVa2,PetscReal *Al1,PetscReal *Al2,void *ctx)		//When other fields like Pi or S (etc.) come into this, declare a IGAPoint pPi or pS and a new PetscReal *UPi or *US for each
	{
		const PetscReal *N0;
		IGAPointGetShapeFuns(pV,0,(const PetscReal**)&N0);									//Value of the shape functions
		PetscInt a,nen=pV->nen,dof=pV->dof;

		PetscReal alfa1[2];																	//Create array to recieve Alfa
		IGAPointFormValue(pAl1,Al1,&alfa1[0]);												//This fills the values

		PetscReal alfa2[2];																	//Create array to recieve Alfa
		IGAPointFormValue(pAl1,Al2,&alfa2[0]);												//This fills the values

		PetscReal xi1=0.0,xi2=0.0;
		if (alfa1[0]>0 || alfa1[1]>0 || alfa1[0]<0 || alfa1[1]<0)
		{
			xi1=1.0;
		}

		if (alfa2[0]>0 || alfa2[1]>0 || alfa2[0]<0 || alfa2[1]<0)
		{
			xi2=1.0;
		}

		PetscReal (*FV1)[dof] = (PetscReal (*)[dof])FVa1;
		PetscReal (*FV2)[dof] = (PetscReal (*)[dof])FVa2;

		if (pV->atboundary)
		{
			return 0;
		}
		else
		{
			for(a=0;a<nen;a++)
			{
				FV1[a][0]+=xi1*N0[a];
				FV1[a][1]+=xi1*N0[a];

				FV2[a][0]+=xi2*N0[a];
				FV2[a][1]+=xi2*N0[a];
			}
		}
		return 0;
	}
//

//System for L2 projection of V^{S}
	#undef  __FUNCT__
	#define __FUNCT__ "VS"
				// VS(pointVs,pointchiUp,pointZ0,pointS,pointZS,KpointVs,FpointVs,Chi0Vs,Z0Vs,S0Vs,ZSVs,NULL);CHKERRQ(ierr);
	PetscErrorCode VS(IGAPoint p,IGAPoint pChi,IGAPoint pZu,IGAPoint pS,IGAPoint pZS,PetscReal *K,PetscReal *F, PetscReal *Chi,PetscReal *Zu,PetscReal *S,PetscReal *ZS,void *ctx)		//When other fields like Pi or S (etc.) come into this, declare a IGAPoint pPi or pS and a new PestcReal *UPi or *US for each
	{
		const PetscReal *N0;
		IGAPointGetShapeFuns(p,0,(const PetscReal**)&N0);									//Value of the shape functions
		PetscInt a,b,i,j,k,l,m,n,nen=p->nen, dof=p->dof;

		//Change to consider G=1
		const PetscReal nu=0.33;
		const PetscReal mu=1.0;
		const PetscReal lambda=2.0*mu*nu/(1.0-2.0*nu);
		const PetscReal eps=mu/100.0;														//Choose later based on whatever Amit says :)

		PetscReal C[3][3][3][3]={0};
		//Creation of elasticity tensor
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

		PetscReal x[2];																		//Vector of reals, size equal to problem's dimension
		IGAPointFormGeomMap(p,x);															//Fills x with the coordinates of p, Gauss's point

		PetscReal S0[8];																	//Assign S to a vector
		PetscReal dS[8][2];																	//And its derivative
		IGAPointFormValue(pS,S,&S0[0]);																
		IGAPointFormGrad (pS,S,&dS[0][0]);													//Same for the gradient

		PetscReal chi0[4];																	//Array to contain the vector chi(0)
		PetscReal d2_Chi0[4][2][2];															//Same for its Hessian
		IGAPointFormValue(pChi,Chi,&chi0[0]);												//Assign chi to its container
		IGAPointFormHess (pChi,Chi,&d2_Chi0[0][0][0]);										//This should be the 3-rd order tensor Chi_{i,jk} (remember that we are storing Chi_{kl} as a column vector)

		PetscReal d_z0[2][2];																//Same for its gradient
		PetscReal d3_z0[2][2][2][2];														//Same for its 3rd order partial derivatives
		IGAPointFormGrad (pZu,Zu,&d_z0[0][0]);												//Same for the gradient
		IGAPointFormDer3 (pZu,Zu,&d3_z0[0][0][0][0]);										//Same for the thir derivatives

		PetscReal d2_Z[4][2][2];																	//gradZ, grad part of S
		IGAPointFormHess (pZS,ZS,&d2_Z[0][0][0]);

		//The four non-zero components of Chi are stored as a vector, restore them to an array with the correct indexing for value and derivative
		PetscReal fullChi[3][3]={0};
		fullChi[0][0]=chi0[0]; 	fullChi[0][1]=chi0[1];
		fullChi[1][0]=chi0[2]; 	fullChi[1][1]=chi0[3];

		PetscReal fulld2_Chi[3][3][3][3]={0};
		fulld2_Chi[0][0][0][0]=d2_Chi0[0][0][0]; fulld2_Chi[0][0][0][1]=d2_Chi0[0][0][1]; fulld2_Chi[0][0][1][0]=d2_Chi0[0][1][0]; fulld2_Chi[0][0][1][1]=d2_Chi0[0][1][1];
		fulld2_Chi[0][1][0][0]=d2_Chi0[1][0][0]; fulld2_Chi[0][1][0][1]=d2_Chi0[1][0][1]; fulld2_Chi[0][1][1][0]=d2_Chi0[1][1][0]; fulld2_Chi[0][1][1][1]=d2_Chi0[1][1][1];
		fulld2_Chi[1][0][0][0]=d2_Chi0[2][0][0]; fulld2_Chi[1][0][0][1]=d2_Chi0[2][0][1]; fulld2_Chi[1][0][1][0]=d2_Chi0[2][1][0]; fulld2_Chi[1][0][1][1]=d2_Chi0[2][1][1];
		fulld2_Chi[1][1][0][0]=d2_Chi0[3][0][0]; fulld2_Chi[1][1][0][1]=d2_Chi0[3][0][1]; fulld2_Chi[1][1][1][0]=d2_Chi0[3][1][0]; fulld2_Chi[1][1][1][1]=d2_Chi0[3][1][1];

		//Expanding z (and derivatives) to 3 components, more convenient for sums in for loops
		PetscReal fulld_z[3][3]={0};
		fulld_z[0][0]=d_z0[0][0]; fulld_z[0][1]=d_z0[0][1];
		fulld_z[1][0]=d_z0[1][0]; fulld_z[1][1]=d_z0[1][1];

		PetscReal fulld3_z[3][3][3][3]={0};
		fulld3_z[0][0][0][0]=d3_z0[0][0][0][0]; fulld3_z[0][0][0][1]=d3_z0[0][0][0][1]; fulld3_z[0][0][1][0]=d3_z0[0][0][1][0]; fulld3_z[0][0][1][1]=d3_z0[0][0][1][1];
		fulld3_z[0][1][0][0]=d3_z0[0][1][0][0]; fulld3_z[0][1][0][1]=d3_z0[0][1][0][1]; fulld3_z[0][1][1][0]=d3_z0[0][1][1][0]; fulld3_z[0][1][1][1]=d3_z0[0][1][1][1];
		fulld3_z[1][0][0][0]=d3_z0[1][0][0][0]; fulld3_z[1][0][0][1]=d3_z0[1][0][0][1]; fulld3_z[1][0][1][0]=d3_z0[1][0][1][0]; fulld3_z[1][0][1][1]=d3_z0[1][0][1][1]; 
		fulld3_z[1][1][0][0]=d3_z0[1][1][0][0]; fulld3_z[1][1][0][1]=d3_z0[1][1][0][1]; fulld3_z[1][1][1][0]=d3_z0[1][1][1][0]; fulld3_z[1][1][1][1]=d3_z0[1][1][1][1];

		//Inflate stored vectors to full tensor form
		PetscReal fullS[3][3][3]={0};
		fullS[0][0][0]=S0[0]; fullS[0][0][1]=S0[1];											//Expand S to full 3rd order form, only non-zero elements
		fullS[0][1][0]=S0[2]; fullS[0][1][1]=S0[3];
		fullS[1][0][0]=S0[4]; fullS[1][0][1]=S0[5];
		fullS[1][1][0]=S0[6]; fullS[1][1][1]=S0[7];

		PetscReal fulld_S[3][3][3][3]={0};													//Expand grad(S) to full tensor order form, only non-zero elements
		fulld_S[0][0][0][0]=dS[0][0]; fulld_S[0][0][0][1]=dS[0][1]; 
		fulld_S[0][0][1][0]=dS[1][0]; fulld_S[0][0][1][1]=dS[1][1]; 
		fulld_S[0][1][0][0]=dS[2][0]; fulld_S[0][1][0][1]=dS[2][1]; 
		fulld_S[0][1][1][0]=dS[3][0]; fulld_S[0][1][1][1]=dS[3][1];
		fulld_S[1][0][0][0]=dS[4][0]; fulld_S[1][0][0][1]=dS[4][1]; 
		fulld_S[1][0][1][0]=dS[5][0]; fulld_S[1][0][1][1]=dS[5][1]; 
		fulld_S[1][1][0][0]=dS[6][0]; fulld_S[1][1][0][1]=dS[6][1]; 
		fulld_S[1][1][1][0]=dS[7][0]; fulld_S[1][1][1][1]=dS[7][1];

		PetscReal fulld2_Z[3][3][3][3]={0};
		fulld2_Z[0][0][0][0]=d2_Z[0][0][0];	fulld2_Z[0][0][0][1]=d2_Z[0][0][1];	fulld2_Z[0][0][1][0]=d2_Z[0][1][0];	fulld2_Z[0][0][1][1]=d2_Z[0][1][1];
		fulld2_Z[0][1][0][0]=d2_Z[1][0][0];	fulld2_Z[0][1][0][1]=d2_Z[1][0][1];	fulld2_Z[0][1][1][0]=d2_Z[1][1][0];	fulld2_Z[0][1][1][1]=d2_Z[1][1][1];
		fulld2_Z[1][0][0][0]=d2_Z[2][0][0];	fulld2_Z[1][0][0][1]=d2_Z[2][0][1];	fulld2_Z[1][0][1][0]=d2_Z[2][1][0];	fulld2_Z[1][0][1][1]=d2_Z[2][1][1];
		fulld2_Z[1][1][0][0]=d2_Z[3][0][0];	fulld2_Z[1][1][0][1]=d2_Z[3][0][1];	fulld2_Z[1][1][1][0]=d2_Z[3][1][0];	fulld2_Z[1][1][1][1]=d2_Z[3][1][1];

		PetscReal (*KVS)[dof][nen][dof] = (typeof(KVS)) K;
		PetscReal (*FVS)[dof] = (PetscReal (*)[dof])F;

		PetscReal v[3]={0.0};

		if (p->atboundary)
		{
			return 0;
		}
		else
		{
			for (a=0; a<nen; a++) 
			{
				for (b=0; b<nen; b++) 
				{
					for (i=0; i<dof; i++)
					{
						KVS[a][i][b][i]=N0[a]*N0[b];
					}
				}
			}

			for(a=0 ;a<nen; a++)
			{
				for (i=0; i<dof; i++)
				{
					if (i==0)
					{
						v[0]=N0[a]; v[1]=0.0; v[2]=0.0;
					}
					else if (i==1)
					{
						v[0]=0.0; v[1]=N0[a]; v[2]=0.0;
					}
					else if (i==2)
					{
						v[0]=0.0; v[1]=0.0; v[2]=N0[a];
					}

					FVS[a][i]=0.0;
					//Part for div(dPsi/dS)
					for (j=0;j<3;j++)
					{
						for(k=0;k<3;k++)
						{
							for(l=0;l<3;l++)
							{
								for(m=0;m<3;m++)
								{
									//FVS[a][i]+=eps*(fulld_S[j][k][l][l]+fulld3_z[j][k][l][l]+fulld2_Chi[j][k][l][l]-fulld2_Z[j][k][l][l])*fullS[j][k][m]*v[m];
									//This part is when considering dPsi/dS symmetrized on the last 2 indices
									FVS[a][i]+=eps*0.5*((fulld_S[j][k][l][l]+fulld3_z[j][k][l][l]+fulld2_Chi[j][k][l][l]-fulld2_Z[j][k][l][l])
													   +(fulld_S[j][l][k][l]+fulld3_z[j][l][k][l]+fulld2_Chi[j][l][k][l]-fulld2_Z[j][l][k][l]))*fullS[j][k][m]*v[m];
								}
							}
						}
					}
					//Part for div(grad(p))
					for (j=0;j<3;j++)
					{
						for(k=0;k<3;k++)
						{
							for(l=0;l<3;l++)
							{
								for(m=0;m<3;m++)
								{
									for(n=0;n<3;n++)
									{
										FVS[a][i]+=C[j][k][l][m]*(-fulld_z[l][m]-fullChi[l][m])*fullS[j][k][n]*v[n];
									}
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

//System for L2 projection of (\partial_{S}\psi)_{ijk,k}
	#undef  __FUNCT__
	#define __FUNCT__ "ProjectDiv"
				// VS(pointVs,pointchiUp,pointZ0,pointS,pointZS,KpointVs,FpointVs,Chi0Vs,Z0Vs,S0Vs,ZSVs,NULL);CHKERRQ(ierr);
	PetscErrorCode ProjectDiv(IGAPoint p,IGAPoint pChi,IGAPoint pZu,IGAPoint pS,IGAPoint pZS,PetscReal *K,PetscReal *F, PetscReal *Chi,PetscReal *Zu,PetscReal *S,PetscReal *ZS,void *ctx)		//When other fields like Pi or S (etc.) come into this, declare a IGAPoint pPi or pS and a new PestcReal *UPi or *US for each
	{
		const PetscReal *N0;
		IGAPointGetShapeFuns(p,0,(const PetscReal**)&N0);									//Value of the shape functions
		PetscInt a,b,i,l,m,n,nen=p->nen, dof=p->dof;

		const PetscReal mu=1.0;
		const PetscReal eps=mu/100.0;														//Choose later based on whatever Amit says :)

		PetscReal x[2];																		//Vector of reals, size equal to problem's dimension
		IGAPointFormGeomMap(p,x);															//Fills x with the coordinates of p, Gauss's point

		PetscReal dS[8][2];																	//And its derivative														
		IGAPointFormGrad (pS,S,&dS[0][0]);													//Same for the gradient

		PetscReal d2_Chi0[4][2][2];															//Same for its Hessian
		IGAPointFormHess (pChi,Chi,&d2_Chi0[0][0][0]);										//This should be the 3-rd order tensor Chi_{i,jk} (remember that we are storing Chi_{kl} as a column vector)

		PetscReal d3_z0[2][2][2][2];														//Same for z 3rd order partial derivatives
		IGAPointFormDer3 (pZu,Zu,&d3_z0[0][0][0][0]);										//Put in container

		PetscReal d2_Z[4][2][2];															//gradZ, grad part of S=grad(Z)+S^\perp
		IGAPointFormHess (pZS,ZS,&d2_Z[0][0][0]);

		//The four non-zero components of Chi are stored as a vector, restore them to an array with the correct indexing for value and derivative

		PetscReal fulld2_Chi[3][3][3][3]={0};
		fulld2_Chi[0][0][0][0]=d2_Chi0[0][0][0]; fulld2_Chi[0][0][0][1]=d2_Chi0[0][0][1]; fulld2_Chi[0][0][1][0]=d2_Chi0[0][1][0]; fulld2_Chi[0][0][1][1]=d2_Chi0[0][1][1];
		fulld2_Chi[0][1][0][0]=d2_Chi0[1][0][0]; fulld2_Chi[0][1][0][1]=d2_Chi0[1][0][1]; fulld2_Chi[0][1][1][0]=d2_Chi0[1][1][0]; fulld2_Chi[0][1][1][1]=d2_Chi0[1][1][1];
		fulld2_Chi[1][0][0][0]=d2_Chi0[2][0][0]; fulld2_Chi[1][0][0][1]=d2_Chi0[2][0][1]; fulld2_Chi[1][0][1][0]=d2_Chi0[2][1][0]; fulld2_Chi[1][0][1][1]=d2_Chi0[2][1][1];
		fulld2_Chi[1][1][0][0]=d2_Chi0[3][0][0]; fulld2_Chi[1][1][0][1]=d2_Chi0[3][0][1]; fulld2_Chi[1][1][1][0]=d2_Chi0[3][1][0]; fulld2_Chi[1][1][1][1]=d2_Chi0[3][1][1];

		PetscReal fulld3_z[3][3][3][3]={0};
		fulld3_z[0][0][0][0]=d3_z0[0][0][0][0]; fulld3_z[0][0][0][1]=d3_z0[0][0][0][1]; fulld3_z[0][0][1][0]=d3_z0[0][0][1][0]; fulld3_z[0][0][1][1]=d3_z0[0][0][1][1];
		fulld3_z[0][1][0][0]=d3_z0[0][1][0][0]; fulld3_z[0][1][0][1]=d3_z0[0][1][0][1]; fulld3_z[0][1][1][0]=d3_z0[0][1][1][0]; fulld3_z[0][1][1][1]=d3_z0[0][1][1][1];
		fulld3_z[1][0][0][0]=d3_z0[1][0][0][0]; fulld3_z[1][0][0][1]=d3_z0[1][0][0][1]; fulld3_z[1][0][1][0]=d3_z0[1][0][1][0]; fulld3_z[1][0][1][1]=d3_z0[1][0][1][1]; 
		fulld3_z[1][1][0][0]=d3_z0[1][1][0][0]; fulld3_z[1][1][0][1]=d3_z0[1][1][0][1]; fulld3_z[1][1][1][0]=d3_z0[1][1][1][0]; fulld3_z[1][1][1][1]=d3_z0[1][1][1][1];

		PetscReal fulld_S[3][3][3][3]={0};													//Expand grad(S) to full tensor order form, only non-zero elements
		fulld_S[0][0][0][0]=dS[0][0]; fulld_S[0][0][0][1]=dS[0][1]; 
		fulld_S[0][0][1][0]=dS[1][0]; fulld_S[0][0][1][1]=dS[1][1]; 
		fulld_S[0][1][0][0]=dS[2][0]; fulld_S[0][1][0][1]=dS[2][1]; 
		fulld_S[0][1][1][0]=dS[3][0]; fulld_S[0][1][1][1]=dS[3][1];
		fulld_S[1][0][0][0]=dS[4][0]; fulld_S[1][0][0][1]=dS[4][1]; 
		fulld_S[1][0][1][0]=dS[5][0]; fulld_S[1][0][1][1]=dS[5][1]; 
		fulld_S[1][1][0][0]=dS[6][0]; fulld_S[1][1][0][1]=dS[6][1]; 
		fulld_S[1][1][1][0]=dS[7][0]; fulld_S[1][1][1][1]=dS[7][1];

		PetscReal fulld2_Z[3][3][3][3]={0};
		fulld2_Z[0][0][0][0]=d2_Z[0][0][0];	fulld2_Z[0][0][0][1]=d2_Z[0][0][1];	fulld2_Z[0][0][1][0]=d2_Z[0][1][0];	fulld2_Z[0][0][1][1]=d2_Z[0][1][1];
		fulld2_Z[0][1][0][0]=d2_Z[1][0][0];	fulld2_Z[0][1][0][1]=d2_Z[1][0][1];	fulld2_Z[0][1][1][0]=d2_Z[1][1][0];	fulld2_Z[0][1][1][1]=d2_Z[1][1][1];
		fulld2_Z[1][0][0][0]=d2_Z[2][0][0];	fulld2_Z[1][0][0][1]=d2_Z[2][0][1];	fulld2_Z[1][0][1][0]=d2_Z[2][1][0];	fulld2_Z[1][0][1][1]=d2_Z[2][1][1];
		fulld2_Z[1][1][0][0]=d2_Z[3][0][0];	fulld2_Z[1][1][0][1]=d2_Z[3][0][1];	fulld2_Z[1][1][1][0]=d2_Z[3][1][0];	fulld2_Z[1][1][1][1]=d2_Z[3][1][1];

		PetscReal (*KProjDiv)[dof][nen][dof] = (typeof(KProjDiv)) K;
		PetscReal (*FProjDiv)[dof] = (PetscReal (*)[dof])F;

		PetscReal v[3][3]={0};

		if (p->atboundary)
		{
			return 0;
		}
		else
		{
			for (a=0; a<nen; a++) 
			{
				for (b=0; b<nen; b++) 
				{
					for (i=0; i<dof; i++)
					{
						KProjDiv[a][i][b][i]=N0[a]*N0[b];
					}
				}
			}

			for(a=0 ;a<nen; a++)
			{
				for (i=0; i<dof; i++)
				{
					if (i==0)
					{
						v[0][0]=N0[a]; v[0][1]=0.0; v[1][0]=0.0; v[1][1]=0.0;
						m=0; n=0;
					}
					else if (i==1)
					{
						v[0][0]=0.0; v[0][1]=N0[a]; v[1][0]=0.0; v[1][1]=0.0;
						m=0; n=1;
					}
					else if (i==2)
					{
						v[0][0]=0.0; v[0][1]=0.0; v[1][0]=N0[a]; v[1][1]=0.0;
						m=1; n=0;
					}
					else if (i==3)
					{
						v[0][0]=0.0; v[0][1]=0.0; v[1][0]=0.0; v[1][1]=N0[a];
						m=1; n=1;
					}

					FProjDiv[a][i]=0.0;
					//Part for div(dPsi/dS)
					for(l=0;l<3;l++)
					{
						FProjDiv[a][i]+=eps*(fulld_S[m][n][l][l]+fulld3_z[m][n][l][l]+fulld2_Chi[m][n][l][l]-fulld2_Z[m][n][l][l])*v[m][n];
					}
				}
			}
		}
		return 0;
	}
//

//System for L2 projection of \hat{Ue} and skew part of \hat{Ue}
	#undef  __FUNCT__
	#define __FUNCT__ "skwUe"
	//PetscErrorCode Stress(IGAPoint p,IGAPoint pU, IGAPoint pHs,IGAPoint pChi,IGAPoint pZu,PetscReal *K,PetscReal *F,PetscReal *U,PetscReal *HS, PetscReal *Chi,PetscReal *Zu,void *ctx)		//When other fields like Pi or S (etc.) come into this, declare a IGAPoint pPi or pS and a new PescReal *UPi or *US for each
	PetscErrorCode skwUe(IGAPoint p,IGAPoint pChi,IGAPoint pZu,PetscReal *K,PetscReal *F, PetscReal *Chi,PetscReal *Zu,PetscInt flag,void *ctx)		//When other fields like Pi or S (etc.) come into this, declare a IGAPoint pPi or pS and a new PestcReal *UPi or *US for each
	{
		const PetscReal *N0;
		IGAPointGetShapeFuns(p,0,(const PetscReal**)&N0);									//Value of the shape functions
		PetscInt a,b,i,k,l,nen=p->nen, dof=p->dof;

		PetscReal chi0[4];																	//Array to contain the vector chi(0)
		IGAPointFormValue(pChi,Chi,&chi0[0]);												//Assign chi to its container

		PetscReal d_Z0[2][2];																//Same for its gradient
		IGAPointFormGrad (pZu,Zu,&d_Z0[0][0]);												//Same for the gradient

		//The four non-zero components of Chi are stored as a vector, restore them to an array with the correct indexing for value and derivative
		PetscReal fullChi[3][3]={0};
		fullChi[0][0]=chi0[0]; 	fullChi[0][1]=chi0[1];
		fullChi[1][0]=chi0[2]; 	fullChi[1][1]=chi0[3];

		//Expanding z (and derivatives) to 3 components, more convenient for sums in for loops
		PetscReal fulld_z[3][3]={0};
		fulld_z[0][0]=d_Z0[0][0]; fulld_z[0][1]=d_Z0[0][1];
		fulld_z[1][0]=d_Z0[1][0]; fulld_z[1][1]=d_Z0[1][1];

		PetscReal (*Kstress)[dof][nen][dof] = (typeof(Kstress)) K;
		PetscReal (*Fstress)[dof] = (PetscReal (*)[dof])F;

		PetscReal v[3][3]={0};

		if (p->atboundary)
		{
			return 0;
		}
		else
		{
			for (a=0; a<nen; a++) 
			{
				for (b=0; b<nen; b++) 
				{
					for (i=0; i<dof; i++)
					{
						Kstress[a][i][b][i]=N0[a]*N0[b];
					}
				}
			}

			for(a=0 ;a<nen; a++)
			{
				for (i=0; i<dof; i++)
				{
					if (i==0)
					{
						v[0][0]=N0[a]; v[0][1]=0.0; v[1][0]=0.0; v[1][1]=0.0;
					}
					else if (i==1)
					{
						v[0][0]=0.0; v[0][1]=N0[a]; v[1][0]=0.0; v[1][1]=0.0;
					}
					else if (i==2)
					{
						v[0][0]=0.0; v[0][1]=0.0; v[1][0]=N0[a]; v[1][1]=0.0;
					}
					else if (i==3)
					{
						v[0][0]=0.0; v[0][1]=0.0; v[1][0]=0.0; v[1][1]=N0[a];
					}
					if(flag==0)
					{
						Fstress[a][i]=0.0;
						for (k=0;k<3;k++)
						{
							for(l=0;l<3;l++)
							{
								Fstress[a][i]+=(-fulld_z[k][l]-fullChi[k][l])*v[k][l];
							}
						}
					}
					else if(flag==1)
					{
						Fstress[a][i]=0.0;
						for (k=0;k<3;k++)
						{
							for(l=0;l<3;l++)
							{
								Fstress[a][i]+=0.5*((-fulld_z[k][l]-fullChi[k][l])-(-fulld_z[l][k]-fullChi[l][k]))*v[k][l];
							}
						}
					}
				}
			}
		}
		return 0;
	}
//

//System for L2 projection of norm of skew part of \hat{Ue}
	#undef  __FUNCT__
	#define __FUNCT__ "NormUeSkwSys"
	PetscErrorCode NormUeSkwSys(IGAPoint p,IGAPoint pChi,IGAPoint pZu,PetscReal *K,PetscReal *F, PetscReal *Chi,PetscReal *Zu,void *ctx)
	{
		const PetscReal *N0;
		IGAPointGetShapeFuns(p,0,(const PetscReal**)&N0);									//Value of the shape functions
		PetscInt a,b,i,k,l,nen=p->nen, dof=p->dof;

		PetscReal chi0[4];																	//Array to contain the vector chi(0)
		IGAPointFormValue(pChi,Chi,&chi0[0]);												//Assign chi to its container

		PetscReal d_Z0[2][2];																//Same for its gradient
		IGAPointFormGrad (pZu,Zu,&d_Z0[0][0]);												//Same for the gradient

		//The four non-zero components of Chi are stored as a vector, restore them to an array with the correct indexing for value and derivative
		PetscReal fullChi[3][3]={0};
		fullChi[0][0]=chi0[0]; 	fullChi[0][1]=chi0[1];
		fullChi[1][0]=chi0[2]; 	fullChi[1][1]=chi0[3];

		//Expanding z (and derivatives) to 3 components, more convenient for sums in for loops
		PetscReal fulld_z[3][3]={0};
		fulld_z[0][0]=d_Z0[0][0]; fulld_z[0][1]=d_Z0[0][1];
		fulld_z[1][0]=d_Z0[1][0]; fulld_z[1][1]=d_Z0[1][1];

		PetscReal (*Kstress)[dof][nen][dof] = (typeof(Kstress)) K;
		PetscReal (*Fstress)[dof] = (PetscReal (*)[dof])F;

		PetscReal val;

		if (p->atboundary)
		{
			return 0;
		}
		else
		{
			for (a=0; a<nen; a++) 
			{
				for (b=0; b<nen; b++) 
				{
					for (i=0; i<dof; i++)
					{
						Kstress[a][i][b][i]=N0[a]*N0[b];
					}
				}
			}

			for(a=0 ;a<nen; a++)
			{
				for (i=0; i<dof; i++)
				{	
					Fstress[a][i]=0.0;
					val=0.0;
					for (k=0;k<3;k++)
					{
						for(l=0;l<3;l++)
						{
							val=val+0.5*(-fulld_z[k][l]-fullChi[k][l]+fulld_z[l][k]+fullChi[l][k])*0.5*(-fulld_z[k][l]-fullChi[k][l]+fulld_z[l][k]+fullChi[l][k]);
						}
					}
					Fstress[a][i]+=sqrt(val)*N0[a];
				}
			}
		}
		return 0;
	}
//

//System for L2 projection of exact stress
	#undef  __FUNCT__
	#define __FUNCT__ "L2ProjectionExactStress"
	PetscErrorCode L2ProjectionExactStress(IGAPoint p,PetscReal *K,PetscReal *F,void *ctx)
	{
		if (p->atboundary)
		{
			return 0;																		//Zero Neumann boundary condition
		}

		PetscInt a,b,i;
		PetscInt nen = p->nen;																//Number of shape functions
		PetscInt dim = p->dim;																//Spatial dimensions of the problem
		PetscInt dof = p->dof;

		PetscReal x[dim];																	//Vector of reals, size equal to problem's dimension
		IGAPointFormGeomMap(p,x);															//Fills x with the coordinates of p, Gauss's point

		//g is the function to L2 project
		//Stress has 4 components, in order S(1,1), S(1,2), S(2,1), S(2,2)
		PetscReal g[dof];

		//const PetscReal E=2100000.0*9.81*10000.0;
		//const PetscReal nu=0.33;
		//const PetscReal mu=E/(2.0*(1.0+nu));
		
		const PetscReal nu=0.33;
		const PetscReal mu=1.0;

		//const PetscReal burgers[2]={1.0,0.0};
		//PetscReal Omega=tan(5.0/180.0*ConstPi);
		//PetscReal rho=sqrt(x[0]*x[0]+x[1]*x[1]);
		PetscReal burgers=1.0;

		//This is for a single disclination
		//g[0]=mu*Omega/(2.0*ConstPi*(1.0-nu))*(log(rho)+(x[1]*x[1])/(rho*rho)+nu/(1.0-2.0*nu));
		//g[1]=-mu*Omega/(2.0*ConstPi*(1.0-nu))*x[0]*x[1]/(rho*rho);
		//g[2]=-mu*Omega/(2.0*ConstPi*(1.0-nu))*x[0]*x[1]/(rho*rho);
		//g[3]=mu*Omega/(2.0*ConstPi*(1.0-nu))*(log(rho)+(x[0]*x[0])/(rho*rho)+nu/(1.0-2.0*nu));

		//This for dislocation with burgers vector in x axis
		g[0]=-mu*burgers/(ConstPi*(1.0-nu))/2.0*x[1]*(x[1]*x[1]+3*x[0]*x[0])/((x[0]*x[0]+x[1]*x[1])*(x[0]*x[0]+x[1]*x[1]));
		g[1]=-mu*burgers/(ConstPi*(1.0-nu))/2.0*x[0]*(x[1]*x[1]-x[0]*x[0])/((x[0]*x[0]+x[1]*x[1])*(x[0]*x[0]+x[1]*x[1]));
		g[2]=-mu*burgers/(ConstPi*(1.0-nu))/2.0*x[0]*(x[1]*x[1]-x[0]*x[0])/((x[0]*x[0]+x[1]*x[1])*(x[0]*x[0]+x[1]*x[1]));
		g[3]=-mu*burgers/(ConstPi*(1.0-nu))/2.0*x[1]*(x[1]*x[1]-x[0]*x[0])/((x[0]*x[0]+x[1]*x[1])*(x[0]*x[0]+x[1]*x[1]));

		//This for when burgers vector points in y axis
		//g[0]=-mu*burgers/(ConstPi*(1.0-nu))/2.0*-x[0]*(x[0]*x[0]+3.0*x[1]*x[1])/((x[0]*x[0]+x[1]*x[1])*(x[0]*x[0]+x[1]*x[1]));
		//g[1]=-mu*burgers/(ConstPi*(1.0-nu))/2.0*x[1]*(x[0]*x[0]-x[1]*x[1])/((x[0]*x[0]+x[1]*x[1])*(x[0]*x[0]+x[1]*x[1]));
		//g[2]=-mu*burgers/(ConstPi*(1.0-nu))/2.0*x[1]*(x[0]*x[0]-x[1]*x[1])/((x[0]*x[0]+x[1]*x[1])*(x[0]*x[0]+x[1]*x[1]));
		//g[3]=-mu*burgers/(ConstPi*(1.0-nu))/2.0*-x[0]*(x[0]*x[0]-x[1]*x[1])/((x[0]*x[0]+x[1]*x[1])*(x[0]*x[0]+x[1]*x[1]));

		//This is for a general dislocation
		//g[0]=-mu*burgers[0]/(ConstPi*(1.0-nu))/2.0*x[1]*(x[1]*x[1]+3*x[0]*x[0])/((x[0]*x[0]+x[1]*x[1])*(x[0]*x[0]+x[1]*x[1]))
		//	 -mu*burgers[1]/(ConstPi*(1.0-nu))/2.0*-x[0]*(x[0]*x[0]+3.0*x[1]*x[1])/((x[0]*x[0]+x[1]*x[1])*(x[0]*x[0]+x[1]*x[1]));
		//g[1]=-mu*burgers[0]/(ConstPi*(1.0-nu))/2.0*x[0]*(x[1]*x[1]-x[0]*x[0])/((x[0]*x[0]+x[1]*x[1])*(x[0]*x[0]+x[1]*x[1]))
		//	 -mu*burgers[1]/(ConstPi*(1.0-nu))/2.0*x[1]*(x[0]*x[0]-x[1]*x[1])/((x[0]*x[0]+x[1]*x[1])*(x[0]*x[0]+x[1]*x[1]));
		//g[2]=-mu*burgers[0]/(ConstPi*(1.0-nu))/2.0*x[0]*(x[1]*x[1]-x[0]*x[0])/((x[0]*x[0]+x[1]*x[1])*(x[0]*x[0]+x[1]*x[1]))
		//	 -mu*burgers[1]/(ConstPi*(1.0-nu))/2.0*x[1]*(x[0]*x[0]-x[1]*x[1])/((x[0]*x[0]+x[1]*x[1])*(x[0]*x[0]+x[1]*x[1]));
		//g[3]=-mu*burgers[0]/(ConstPi*(1.0-nu))/2.0*x[1]*(x[1]*x[1]-x[0]*x[0])/((x[0]*x[0]+x[1]*x[1])*(x[0]*x[0]+x[1]*x[1]))
		//	 -mu*burgers[1]/(ConstPi*(1.0-nu))/2.0*-x[0]*(x[0]*x[0]-x[1]*x[1])/((x[0]*x[0]+x[1]*x[1])*(x[0]*x[0]+x[1]*x[1]));

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

//System for L2 projection of Peierls stress
	#undef  __FUNCT__
	#define __FUNCT__ "L2ProjectionPeierlsStress"
	PetscErrorCode L2ProjectionPeierlsStress(IGAPoint p,PetscReal *K,PetscReal *F,void *ctx)
	{
		if (p->atboundary)
		{
			return 0;																		//Zero Neumann boundary condition
		}

		PetscInt a,b,i;
		PetscInt nen = p->nen;																//Number of shape functions
		PetscInt dim = p->dim;																//Spatial dimensions of the problem
		PetscInt dof = p->dof;

		PetscReal x[dim];																	//Vector of reals, size equal to problem's dimension
		IGAPointFormGeomMap(p,x);															//Fills x with the coordinates of p, Gauss's point

		//g is the function to L2 project
		//Stress has 4 components, in order S(1,1), S(1,2), S(2,1), S(2,2)
		PetscReal g[dof];
		
		const PetscReal nu=0.33;
		const PetscReal mu=1.0;

		PetscReal burgers=1.0;

		PetscReal zeta=burgers/(2.0*(1.0-nu));

		//This for dislocation with burgers vector in x axis
		g[0]=0.0;
		g[1]=mu*burgers/(2*ConstPi*(1.0-nu))*x[0]/(x[0]*x[0]+zeta*zeta);
		g[2]=mu*burgers/(2*ConstPi*(1.0-nu))*x[0]/(x[0]*x[0]+zeta*zeta);
		g[3]=0.0;

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

//System for L2 projection of Grad(z0)
	#undef  __FUNCT__
	#define __FUNCT__ "L2ProjectionGradZ"
	PetscErrorCode L2ProjectionGradZ(IGAPoint p,IGAPoint pZ0, PetscReal *K,PetscReal *F,PetscReal *U, void *ctx)
	{
		if (p->atboundary)
		{
			return 0;																		//Zero Neumann boundary condition
		}

		PetscInt a,b,i;
		PetscInt nen = p->nen;																//Number of shape functions
		PetscInt dim = p->dim;																//Spatial dimensions of the problem
		PetscInt dof = p->dof;

		PetscReal x[dim];																	//Vector of reals, size equal to problem's dimension
		IGAPointFormGeomMap(p,x);															//Fills x with the coordinates of p, Gauss's point

		PetscReal dZ0[2][2];																//Same for partial derivatives
		IGAPointFormGrad (pZ0,U,&dZ0[0][0]);

		PetscReal full_dZ0[3][3]={0};
		full_dZ0[0][0]=dZ0[0][0]; full_dZ0[0][1]=dZ0[0][1];
		full_dZ0[1][0]=dZ0[1][0]; full_dZ0[1][1]=dZ0[1][1];

		//g is the function to L2 project
		//div(Chi) has 2 components, div(Chi)(1) and div(Chi)(2)
		PetscReal g[dof];

		g[0]=full_dZ0[0][0];
		g[1]=full_dZ0[0][1];
		g[2]=full_dZ0[1][0];
		g[3]=full_dZ0[1][1];

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

//Creation of solution systems
	PetscErrorCode  ierr;
	ierr = PetscInitialize(&argc,&argv,0,0);CHKERRQ(ierr);										//Always initialize PETSc
	PetscInt dir,side;
	PetscPrintf(PETSC_COMM_WORLD,"Start of process_results \n");

	PetscInt commsize,rank;
	ierr = MPI_Comm_size(PETSC_COMM_WORLD,&commsize);CHKERRQ(ierr);
	ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank);CHKERRQ(ierr);
//

//App context creation and some data
	//Mesh parameters (to fix specific points in z0 system)
	PetscInt b=581;				//Parmeter to choose size of cores, must always be odd, core will be of size 1 unit, rest of the body will be of size b-1 units in each direction
	PetscReal Lx=80.0;
	PetscReal Ly=80.0;
	PetscInt  nx=b;
	PetscInt  ny=b;
	
	AppCtxL2 userL2;
	userL2.Lx     = Lx;
	userL2.Ly     = Ly;
	userL2.nx     = nx;
	userL2.ny     = ny;
	//void* user;

	//Directory to write files to
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

	time_t T=time(NULL);
	struct tm tm=*localtime(&T);
	PetscPrintf(PETSC_COMM_WORLD,"Current time is %02d:%02d:%02d \n",tm.tm_hour,tm.tm_min,tm.tm_sec);

	if(rank==0)
	{
		FILE *source, *dest;
		char buffer[8192];
		size_t bytes;
		source = fopen("./process_resultsS.c","r");
		dest   = fopen("../Results/process_resultsS.c","w");				//Change this to use directory variable "direct"

		while (0 < (bytes = fread(buffer, 1, sizeof(buffer), source)))
			fwrite(buffer, 1, bytes, dest);

		fclose(source);
		fclose(dest);
	}
//

//Read vectors (In case I get the solution for z0 but need to recreate the post processing results)
	IGA igachiUp;
	ierr = IGACreate(PETSC_COMM_WORLD,&igachiUp);CHKERRQ(ierr);
	ierr = IGASetDim(igachiUp,2);CHKERRQ(ierr);														//Spatial dimension of the problem
	ierr = IGASetDof(igachiUp,4);CHKERRQ(ierr);														//Number of degrees of freedom, per node
	ierr = IGASetOrder(igachiUp,2);CHKERRQ(ierr);													//Number of spatial derivatives to calculate
	ierr = IGASetFromOptions(igachiUp);CHKERRQ(ierr);												//Note: The order (or degree) of the shape functions is given by the mesh!
	ierr = IGARead(igachiUp,"./geometry2.dat");CHKERRQ(ierr);
	
	for (dir=0; dir<2; dir++)
	{
		ierr = IGASetRuleType(igachiUp,dir,IGA_RULE_LEGENDRE);CHKERRQ(ierr);
		ierr = IGASetRuleSize(igachiUp,dir,6);CHKERRQ(ierr);
	}
	ierr = IGASetUp(igachiUp);CHKERRQ(ierr);

	Vec chiUp0;
	ierr = IGACreateVec(igachiUp,&chiUp0);CHKERRQ(ierr);  

	char nameChi[]="/ChiUp-2d-0.dat";
	char pathChi[512];
	sprintf(pathChi,"%s%s",direct,nameChi);
	ierr = IGAReadVec(igachiUp,chiUp0,pathChi); CHKERRQ(ierr);

	IGA igaZ0;
	ierr = IGACreate(PETSC_COMM_WORLD,&igaZ0);CHKERRQ(ierr);
	ierr = IGASetDim(igaZ0,2);CHKERRQ(ierr);													//Spatial dimension of the problem
	ierr = IGASetDof(igaZ0,2);CHKERRQ(ierr);													//Number of degrees of freedom, per node
	ierr = IGASetOrder(igaZ0,3);CHKERRQ(ierr);													//Number of spatial derivatives to calculate
	ierr = IGASetFromOptions(igaZ0);CHKERRQ(ierr);												//Note: The order (or degree) of the shape functions is given by the mesh!
	ierr = IGARead(igaZ0,"./geometry3.dat");CHKERRQ(ierr);
	
	for (dir=0; dir<2; dir++)
	{
		ierr = IGASetRuleType(igaZ0,dir,IGA_RULE_LEGENDRE);CHKERRQ(ierr);
		ierr = IGASetRuleSize(igaZ0,dir,6);CHKERRQ(ierr);
	}
	ierr = IGASetUp(igaZ0);CHKERRQ(ierr);

	Vec Z0;
	ierr = IGACreateVec(igaZ0,&Z0);CHKERRQ(ierr);

	char nameZ0[]="/Z0-2d-0.dat";
	char pathZ0[512];
	sprintf(pathZ0,"%s%s",direct,nameZ0);
	ierr = IGAReadVec(igaZ0,Z0,pathZ0); CHKERRQ(ierr);

	IGA igaAl;
	ierr = IGACreate(PETSC_COMM_WORLD,&igaAl);CHKERRQ(ierr);
	ierr = IGASetDim(igaAl,2);CHKERRQ(ierr);														//Spatial dimension of the problem
	ierr = IGASetDof(igaAl,2);CHKERRQ(ierr);														//Number of degrees of freedom, per node
	ierr = IGASetOrder(igaAl,1);CHKERRQ(ierr);													//Number of spatial derivatives to calculate
	ierr = IGASetFromOptions(igaAl);CHKERRQ(ierr);												//Note: The order (or degree) of the shape functions is given by the mesh!
	ierr = IGARead(igaAl,"./geometry.dat");CHKERRQ(ierr);
	
	for (dir=0; dir<2; dir++)
	{
		ierr = IGASetRuleType(igaAl,dir,IGA_RULE_LEGENDRE);CHKERRQ(ierr);
		ierr = IGASetRuleSize(igaAl,dir,6);CHKERRQ(ierr);
	}
	ierr = IGASetUp(igaAl);CHKERRQ(ierr);

	Vec alp0,alInput,alInput1;
	ierr = IGACreateVec(igaAl,&alp0);CHKERRQ(ierr);
	ierr = IGACreateVec(igaAl,&alInput);CHKERRQ(ierr);
	ierr = IGACreateVec(igaAl,&alInput1);CHKERRQ(ierr);

	char nameAlInput1[]="/Input-Al-2d-0.dat";
	char pathAlInput1[512];
	sprintf(pathAlInput1,"%s%s",direct,nameAlInput1);
	ierr = IGAReadVec(igaAl,alInput1,pathAlInput1); CHKERRQ(ierr);

	ierr = VecAXPY(alInput,1.0,alInput1); CHKERRQ(ierr);

	char nameAlp[]="/Alp-2d-0.dat";
	char pathAlp[512];
	sprintf(pathAlp,"%s%s",direct,nameAlp);
	ierr = IGAReadVec(igaAl,alp0,pathAlp); CHKERRQ(ierr);

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

	Vec s0;
	ierr = IGACreateVec(igaS,&s0);CHKERRQ(ierr);  

	char nameS[]="/Input-S-2d-0.dat";
	char pathS[512];
	sprintf(pathS,"%s%s",direct,nameS);
	ierr = IGAReadVec(igaS,s0,pathS); CHKERRQ(ierr);

	IGA igaZS;
	ierr = IGACreate(PETSC_COMM_WORLD,&igaZS);CHKERRQ(ierr);
	ierr = IGASetDim(igaZS,2);CHKERRQ(ierr);														//Spatial dimension of the problem
	ierr = IGASetDof(igaZS,4);CHKERRQ(ierr);														//Number of degrees of freedom, per node
	ierr = IGASetOrder(igaZS,2);CHKERRQ(ierr);														//Number of spatial derivatives to calculate
	ierr = IGASetFromOptions(igaZS);CHKERRQ(ierr);													//Note: The order (or degree) of the shape functions is given by the mesh!
	ierr = IGARead(igaZS,"./geometry.dat");CHKERRQ(ierr);
	
	for (dir=0; dir<2; dir++)
	{
		ierr = IGASetRuleType(igaZS,dir,IGA_RULE_LEGENDRE);CHKERRQ(ierr);
		ierr = IGASetRuleSize(igaZS,dir,6);CHKERRQ(ierr);
	}

	ierr = IGASetUp(igaZS);CHKERRQ(ierr);

	Vec ZS0;
	ierr = IGACreateVec(igaZS,&ZS0);CHKERRQ(ierr);
	char nameZS[]="/ZS-2d-0.dat";
	char pathZS[512];
	sprintf(pathZS,"%s%s",direct,nameZS);
	ierr = IGAReadVec(igaZS,ZS0,pathZS);CHKERRQ(ierr);
//

//Create points and elements
	IGAPoint		pointAlp,pointchiUp,pointZ0,pointS,pointZS;						//point
	IGAElement		elemAlp,elemchiUp,elemZ0,elemS,elemZS;							//element
//

//System for L2 projection of stress (sym part)
	PetscPrintf(PETSC_COMM_WORLD,"\nSystem for Stress starting \n\n");
	T=time(NULL);
	tm=*localtime(&T);
	PetscPrintf(PETSC_COMM_WORLD,"Current time is %02d:%02d:%02d \n",tm.tm_hour,tm.tm_min,tm.tm_sec);
	IGA igaStress;
	ierr = IGACreate(PETSC_COMM_WORLD,&igaStress);CHKERRQ(ierr);
	ierr = IGASetDim(igaStress,2);CHKERRQ(ierr);													//Spatial dimension of the problem
	ierr = IGASetDof(igaStress,4);CHKERRQ(ierr);													//Number of degrees of freedom, per node
	ierr = IGASetOrder(igaStress,2);CHKERRQ(ierr);													//Number of spatial derivatives to calculate
	ierr = IGASetFromOptions(igaStress);CHKERRQ(ierr);												//Note: The order (or degree) of the shape functions is given by the mesh!
	ierr = IGARead(igaStress,"./geometry3.dat");CHKERRQ(ierr);
	
	for (dir=0; dir<2; dir++)
	{
		ierr = IGASetRuleType(igaStress,dir,IGA_RULE_LEGENDRE);CHKERRQ(ierr);
		ierr = IGASetRuleSize(igaStress,dir,6);CHKERRQ(ierr);
	}
	ierr = IGASetUp(igaStress);CHKERRQ(ierr);

	for (dir=0; dir<2; dir++) 
	{
		for (side=0; side<2; side++) 
		{
			//ierr = IGASetBoundaryValue(iga,dir,side,dof,0.0);CHKERRQ(ierr);					// Dirichlet boundary conditions
			ierr = IGASetBoundaryForm(igaStress,dir,side,PETSC_TRUE);CHKERRQ(ierr);  				// Neumann boundary conditions
		}
	}

	Mat KStress;
	Vec sigma0,FStress;

	ierr = IGACreateMat(igaStress,&KStress);CHKERRQ(ierr);
	ierr = IGACreateVec(igaStress,&sigma0);CHKERRQ(ierr);
	ierr = IGACreateVec(igaStress,&FStress);CHKERRQ(ierr);

	IGAPoint		pointStress;															//point
	IGAElement		elemStress;																//element
	PetscReal		*KlocStress,*FlocStress;												//AA y BB
	PetscReal		*KpointStress,*FpointStress;											//KKK y FFF
	const PetscReal	*arrayZ0Stress,*arrayChi0Stress, *arraySStress, *arrayGradZStress;		//arrayU
	Vec				localZ0Stress,localChi0Stress, localSStress, localGradZStress;			//localU
	PetscReal		*Chi0Stress,*Z0Stress, *SStress, *GradZStress;							//U

  	IGAFormSystem	wtfStress;
 	void			*wtf2Stress;

 	KSP kspStress;
	ierr = IGACreateKSP(igaStress,&kspStress);CHKERRQ(ierr);

	// Get local vectors Z0 and Chi0 and arrays
	ierr = IGAGetLocalVecArray(igachiUp,chiUp0,&localChi0Stress,&arrayChi0Stress);CHKERRQ(ierr);
	ierr = IGAGetLocalVecArray(igaZ0,Z0,&localZ0Stress,&arrayZ0Stress);CHKERRQ(ierr);
	ierr = IGAGetLocalVecArray(igaS,s0,&localSStress,&arraySStress);CHKERRQ(ierr);
	ierr = IGAGetLocalVecArray(igaZS,ZS0,&localGradZStress,&arrayGradZStress);CHKERRQ(ierr);

	// Element loop
	ierr = IGABeginElement(igaStress,&elemStress);CHKERRQ(ierr);
	ierr = IGABeginElement(igaZ0,&elemZ0);CHKERRQ(ierr);
	ierr = IGABeginElement(igachiUp,&elemchiUp);CHKERRQ(ierr);
	ierr = IGABeginElement(igaS,&elemS);CHKERRQ(ierr);
	ierr = IGABeginElement(igaZS,&elemZS);CHKERRQ(ierr);

	while (IGANextElement(igaStress,elemStress)) 
	{
		IGANextElement(igaZ0,elemZ0);
		IGANextElement(igachiUp,elemchiUp);
		IGANextElement(igaS,elemS);
		IGANextElement(igaZS,elemZS);

		ierr = IGAElementGetWorkMat(elemStress,&KlocStress);CHKERRQ(ierr);
		ierr = IGAElementGetWorkVec(elemStress,&FlocStress);CHKERRQ(ierr);
		ierr = IGAElementGetValues(elemZ0,arrayZ0Stress,&Z0Stress);CHKERRQ(ierr);
		ierr = IGAElementGetValues(elemchiUp,arrayChi0Stress,&Chi0Stress);CHKERRQ(ierr);
		ierr = IGAElementGetValues(elemS,arraySStress,&SStress);CHKERRQ(ierr);
		ierr = IGAElementGetValues(elemZS,arrayGradZStress,&GradZStress);CHKERRQ(ierr);

		// FormSystem loop
		while (IGAElementNextFormSystem(elemStress,&wtfStress,&wtf2Stress)) 
		{
		// Quadrature loop
			ierr = IGAElementBeginPoint(elemStress,&pointStress);CHKERRQ(ierr);
			ierr = IGAElementBeginPoint(elemZ0,&pointZ0);CHKERRQ(ierr);
			ierr = IGAElementBeginPoint(elemchiUp,&pointchiUp);CHKERRQ(ierr);
			ierr = IGAElementBeginPoint(elemS,&pointS);CHKERRQ(ierr);
			ierr = IGAElementBeginPoint(elemZS,&pointZS);CHKERRQ(ierr);

			while (IGAElementNextPoint(elemStress,pointStress))
			{
				if(pointStress->atboundary==1)
				{
					//IGAElementNextPoint(elemZ0,pointZ0);
					//IGAElementNextPoint(elemchiUp,pointchiUp);
				}

				if(pointZ0->atboundary==1)
				{
					PetscPrintf(PETSC_COMM_WORLD,"Hola1 en stress \n");
				}

				if(pointStress->atboundary==0 && pointZ0->atboundary==0 && pointchiUp->atboundary==0)
				{
					IGAElementNextPoint(elemZ0,pointZ0);
					IGAElementNextPoint(elemchiUp,pointchiUp);
					IGAElementNextPoint(elemS,pointS);
					IGAElementNextPoint(elemZS,pointZS);

					ierr = IGAPointGetWorkMat(pointStress,&KpointStress);CHKERRQ(ierr);
					ierr = IGAPointGetWorkVec(pointStress,&FpointStress);CHKERRQ(ierr);
					//	   Stress(IGAPoint p,IGAPoint pChi,IGAPoint pZu,IGAPoint pS, IGAPoint pGradZ, PetscReal *K,PetscReal *F, PetscReal *Chi,PetscReal *Zu,PetscReal *S, PetscReal *UgradZ, void *ctx)
					ierr = Stress(pointStress,pointchiUp,pointZ0,pointS,pointZS,KpointStress,FpointStress,Chi0Stress,Z0Stress,SStress,GradZStress,NULL);CHKERRQ(ierr);
					ierr = IGAPointAddMat(pointStress,KpointStress,KlocStress);CHKERRQ(ierr);
					ierr = IGAPointAddVec(pointStress,FpointStress,FlocStress);CHKERRQ(ierr);
				}
			}
			if (pointZ0->index != -1)
			{
				IGAElementNextPoint(elemZ0,pointZ0);
			}
			while (pointchiUp->index != -1)
			{
				IGAElementNextPoint(elemchiUp,pointchiUp);
			}
			while (pointS->index != -1)
			{
				IGAElementNextPoint(elemS,pointS);
			}
			while (pointZS->index != -1)
			{
				IGAElementNextPoint(elemZS,pointZS);
			}

			ierr = IGAElementEndPoint(elemStress,&pointStress);CHKERRQ(ierr);
			ierr = IGAElementEndPoint(elemZ0,&pointZ0);CHKERRQ(ierr);
			ierr = IGAElementEndPoint(elemchiUp,&pointchiUp);CHKERRQ(ierr);
			ierr = IGAElementEndPoint(elemS,&pointS);CHKERRQ(ierr);
			ierr = IGAElementEndPoint(elemZS,&pointZS);CHKERRQ(ierr);
		}

		//ierr = IGAElementFixSystem(elemStress,KlocStress,FlocStress);CHKERRQ(ierr);					//This sets Dirichlet condition ¿? (Yes, this applies the conditions from IGASetBoundaryValue)
		ierr = IGAElementAssembleMat(elemStress,KlocStress,KStress);CHKERRQ(ierr);
		ierr = IGAElementAssembleVec(elemStress,FlocStress,FStress);CHKERRQ(ierr);
	}
	IGANextElement(igaZ0,elemZ0);
	IGANextElement(igachiUp,elemchiUp);
	IGANextElement(igaS,elemS);
	IGANextElement(igaZS,elemZS);

	ierr = IGAEndElement(igaStress,&elemStress);CHKERRQ(ierr);
	ierr = IGAEndElement(igaZ0,&elemZ0);CHKERRQ(ierr);
	ierr = IGAEndElement(igachiUp,&elemchiUp);CHKERRQ(ierr);
	ierr = IGAEndElement(igaS,&elemS);CHKERRQ(ierr);
	ierr = IGAEndElement(igaZS,&elemZS);CHKERRQ(ierr);

	// Restore local vectors u, Z0, Chi0 and arrays
	ierr = IGARestoreLocalVecArray(igaZ0,Z0,&localZ0Stress,&arrayZ0Stress);CHKERRQ(ierr);
	ierr = IGARestoreLocalVecArray(igachiUp,chiUp0,&localChi0Stress,&arrayChi0Stress);CHKERRQ(ierr);
	ierr = IGARestoreLocalVecArray(igaS,s0,&localSStress,&arraySStress);CHKERRQ(ierr);
	ierr = IGARestoreLocalVecArray(igaZS,ZS0,&localGradZStress,&arrayGradZStress);CHKERRQ(ierr);

	ierr = MatAssemblyBegin(KStress,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyEnd  (KStress,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

	ierr = VecAssemblyBegin(FStress);CHKERRQ(ierr);
	ierr = VecAssemblyEnd  (FStress);CHKERRQ(ierr);

	ierr = KSPSetOperators(kspStress,KStress,KStress);CHKERRQ(ierr);
	PC pcStress;
	ierr = KSPGetPC(kspStress,&pcStress); CHKERRQ(ierr);
	ierr = PCSetType(pcStress,PCSOR); CHKERRQ(ierr);
	ierr = KSPSetFromOptions(kspStress);CHKERRQ(ierr);
	ierr = KSPSetType(kspStress,KSPFCG);
	ierr = KSPSetTolerances(kspStress,1e-20,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);CHKERRQ(ierr);
	ierr = KSPSolve(kspStress,FStress,sigma0);CHKERRQ(ierr);

	ierr = KSPDestroy(&kspStress);CHKERRQ(ierr);
	ierr = MatDestroy(&KStress);CHKERRQ(ierr);
	ierr = VecDestroy(&FStress);CHKERRQ(ierr);

	char nameStress[]="/sigma-2d-0.dat";
	char pathStress[512];
	sprintf(pathStress,"%s%s",direct,nameStress);
	ierr = IGAWriteVec(igaStress,sigma0,pathStress);CHKERRQ(ierr);	
//

//System for L2 projection of V^{S}
	PetscPrintf(PETSC_COMM_WORLD,"\nSystem for V-S starting \n\n");
	T=time(NULL);
	tm=*localtime(&T);
	PetscPrintf(PETSC_COMM_WORLD,"Current time is %02d:%02d:%02d \n",tm.tm_hour,tm.tm_min,tm.tm_sec);
	IGA igaVs;
	ierr = IGACreate(PETSC_COMM_WORLD,&igaVs);CHKERRQ(ierr);
	ierr = IGASetDim(igaVs,2);CHKERRQ(ierr);													//Spatial dimension of the problem
	ierr = IGASetDof(igaVs,2);CHKERRQ(ierr);													//Number of degrees of freedom, per node
	ierr = IGASetOrder(igaVs,1);CHKERRQ(ierr);													//Number of spatial derivatives to calculate
	ierr = IGASetFromOptions(igaVs);CHKERRQ(ierr);												//Note: The order (or degree) of the shape functions is given by the mesh!
	ierr = IGARead(igaVs,"./geometry3.dat");CHKERRQ(ierr);
	
	for (dir=0; dir<2; dir++)
	{
		ierr = IGASetRuleType(igaVs,dir,IGA_RULE_LEGENDRE);CHKERRQ(ierr);
		ierr = IGASetRuleSize(igaVs,dir,6);CHKERRQ(ierr);
	}
	ierr = IGASetUp(igaVs);CHKERRQ(ierr);

	for (dir=0; dir<2; dir++) 
	{
		for (side=0; side<2; side++) 
		{
			//ierr = IGASetBoundaryValue(iga,dir,side,dof,0.0);CHKERRQ(ierr);					// Dirichlet boundary conditions
			ierr = IGASetBoundaryForm(igaVs,dir,side,PETSC_TRUE);CHKERRQ(ierr);  				// Neumann boundary conditions
		}
	}

	Mat KVs;
	Vec Vs0,FVs;

	ierr = IGACreateMat(igaVs,&KVs);CHKERRQ(ierr);
	ierr = IGACreateVec(igaVs,&Vs0);CHKERRQ(ierr);
	ierr = IGACreateVec(igaVs,&FVs);CHKERRQ(ierr);

	IGAPoint		pointVs;															//point
	IGAElement		elemVs;																//element
	PetscReal		*KlocVs,*FlocVs;													//AA y BB
	PetscReal		*KpointVs,*FpointVs;												//KKK y FFF
	const PetscReal	*arrayZ0Vs,*arrayChi0Vs,*arrayS0Vs,*arrayZSVs;						//arrayU
	Vec				localZ0Vs,localChi0Vs,localS0Vs,localZSVs;							//localU
	PetscReal		*Chi0Vs,*Z0Vs,*S0Vs,*ZSVs;											//U

  	IGAFormSystem	wtfVs;
 	void			*wtf2Vs;

 	KSP kspVs;
	ierr = IGACreateKSP(igaVs,&kspVs);CHKERRQ(ierr);

	// Get local vectors Z0 and Chi0 and arrays
	ierr = IGAGetLocalVecArray(igachiUp,chiUp0,&localChi0Vs,&arrayChi0Vs);CHKERRQ(ierr);
	ierr = IGAGetLocalVecArray(igaZ0,Z0,&localZ0Vs,&arrayZ0Vs);CHKERRQ(ierr);
	ierr = IGAGetLocalVecArray(igaS,s0,&localS0Vs,&arrayS0Vs);CHKERRQ(ierr);
	ierr = IGAGetLocalVecArray(igaZS,ZS0,&localZSVs,&arrayZSVs);CHKERRQ(ierr);

	// Element loop
	ierr = IGABeginElement(igaVs,&elemVs);CHKERRQ(ierr);
	ierr = IGABeginElement(igaZ0,&elemZ0);CHKERRQ(ierr);
	ierr = IGABeginElement(igachiUp,&elemchiUp);CHKERRQ(ierr);
	ierr = IGABeginElement(igaS,&elemS);CHKERRQ(ierr);
	ierr = IGABeginElement(igaZS,&elemZS);CHKERRQ(ierr);

	while (IGANextElement(igaVs,elemVs)) 
	{
		IGANextElement(igaZ0,elemZ0);
		IGANextElement(igachiUp,elemchiUp);
		IGANextElement(igaS,elemS);
		IGANextElement(igaZS,elemZS);

		ierr = IGAElementGetWorkMat(elemVs,&KlocVs);CHKERRQ(ierr);
		ierr = IGAElementGetWorkVec(elemVs,&FlocVs);CHKERRQ(ierr);
		ierr = IGAElementGetValues(elemZ0,arrayZ0Vs,&Z0Vs);CHKERRQ(ierr);
		ierr = IGAElementGetValues(elemchiUp,arrayChi0Vs,&Chi0Vs);CHKERRQ(ierr);
		ierr = IGAElementGetValues(elemS,arrayS0Vs,&S0Vs);CHKERRQ(ierr);
		ierr = IGAElementGetValues(elemZS,arrayZSVs,&ZSVs);CHKERRQ(ierr);

		// FormSystem loop
		while (IGAElementNextFormSystem(elemVs,&wtfVs,&wtf2Vs)) 
		{
		// Quadrature loop
			ierr = IGAElementBeginPoint(elemVs,&pointVs);CHKERRQ(ierr);
			ierr = IGAElementBeginPoint(elemZ0,&pointZ0);CHKERRQ(ierr);
			ierr = IGAElementBeginPoint(elemchiUp,&pointchiUp);CHKERRQ(ierr);
			ierr = IGAElementBeginPoint(elemS,&pointS);CHKERRQ(ierr);
			ierr = IGAElementBeginPoint(elemZS,&pointZS);CHKERRQ(ierr);

			while (IGAElementNextPoint(elemVs,pointVs))
			{
				if(pointVs->atboundary==1)
				{
					//IGAElementNextPoint(elemZ0,pointZ0);
					//IGAElementNextPoint(elemchiUp,pointchiUp);
				}

				if(pointVs->atboundary==0 && pointZ0->atboundary==0 && pointchiUp->atboundary==0 && pointS->atboundary==0)
				{
					IGAElementNextPoint(elemZ0,pointZ0);
					IGAElementNextPoint(elemchiUp,pointchiUp);
					IGAElementNextPoint(elemS,pointS);
					IGAElementNextPoint(elemZS,pointZS);

					ierr = IGAPointGetWorkMat(pointVs,&KpointVs);CHKERRQ(ierr);
					ierr = IGAPointGetWorkVec(pointVs,&FpointVs);CHKERRQ(ierr);
					
					ierr = VS(pointVs,pointchiUp,pointZ0,pointS,pointZS,KpointVs,FpointVs,Chi0Vs,Z0Vs,S0Vs,ZSVs,NULL);CHKERRQ(ierr);
					ierr = IGAPointAddMat(pointVs,KpointVs,KlocVs);CHKERRQ(ierr);
					ierr = IGAPointAddVec(pointVs,FpointVs,FlocVs);CHKERRQ(ierr);
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
			while (pointS->index != -1)
			{
				IGAElementNextPoint(elemS,pointS);
			}
			while (pointZS->index != -1)
			{
				IGAElementNextPoint(elemZS,pointZS);
			}
			ierr = IGAElementEndPoint(elemVs,&pointVs);CHKERRQ(ierr);
			ierr = IGAElementEndPoint(elemZ0,&pointZ0);CHKERRQ(ierr);
			ierr = IGAElementEndPoint(elemchiUp,&pointchiUp);CHKERRQ(ierr);
			ierr = IGAElementEndPoint(elemS,&pointS);CHKERRQ(ierr);
			ierr = IGAElementEndPoint(elemZS,&pointZS);CHKERRQ(ierr);
		}

		//ierr = IGAElementFixSystem(elemStress,KlocStress,FlocStress);CHKERRQ(ierr);					//This sets Dirichlet condition ¿? (Yes, this applies the conditions from IGASetBoundaryValue)
		ierr = IGAElementAssembleMat(elemVs,KlocVs,KVs);CHKERRQ(ierr);
		ierr = IGAElementAssembleVec(elemVs,FlocVs,FVs);CHKERRQ(ierr);
	}
	IGANextElement(igaZ0,elemZ0);
	IGANextElement(igachiUp,elemchiUp);
	IGANextElement(igaS,elemS);
	IGANextElement(igaZS,elemZS);

	ierr = IGAEndElement(igaVs,&elemVs);CHKERRQ(ierr);
	ierr = IGAEndElement(igaZ0,&elemZ0);CHKERRQ(ierr);
	ierr = IGAEndElement(igachiUp,&elemchiUp);CHKERRQ(ierr);
	ierr = IGAEndElement(igaS,&elemS);CHKERRQ(ierr);
	ierr = IGAEndElement(igaZS,&elemZS);CHKERRQ(ierr);

	// Restore local vectors u, Z0, Chi0 and arrays
	ierr = IGARestoreLocalVecArray(igaZ0,Z0,&localZ0Vs,&arrayZ0Vs);CHKERRQ(ierr);
	ierr = IGARestoreLocalVecArray(igachiUp,chiUp0,&localChi0Vs,&arrayChi0Vs);CHKERRQ(ierr);
	ierr = IGARestoreLocalVecArray(igaS,s0,&localS0Vs,&arrayS0Vs);CHKERRQ(ierr);
	ierr = IGARestoreLocalVecArray(igaZS,ZS0,&localZSVs,&arrayZSVs);CHKERRQ(ierr);

	ierr = MatAssemblyBegin(KVs,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyEnd  (KVs,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

	ierr = VecAssemblyBegin(FVs);CHKERRQ(ierr);
	ierr = VecAssemblyEnd  (FVs);CHKERRQ(ierr);

	ierr = KSPSetOperators(kspVs,KVs,KVs);CHKERRQ(ierr);
	//PC pcVs;
	//ierr = KSPGetPC(kspVs,&pcVs); CHKERRQ(ierr);
	//ierr = PCSetType(pcVs,PCLU); CHKERRQ(ierr);
	//ierr = PCFactorSetMatSolverType(pcVs,MATSOLVERMUMPS); CHKERRQ(ierr);
	ierr = KSPSetFromOptions(kspVs);CHKERRQ(ierr);
	ierr = KSPSetType(kspVs,KSPFCG);
	ierr = KSPSetTolerances(kspVs,1e-20,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);CHKERRQ(ierr);
	ierr = KSPSolve(kspVs,FVs,Vs0);CHKERRQ(ierr);

	ierr = KSPDestroy(&kspVs);CHKERRQ(ierr);
	ierr = MatDestroy(&KVs);CHKERRQ(ierr);
	ierr = VecDestroy(&FVs);CHKERRQ(ierr);

	char nameVs[]="/Vs-2d-0.dat";
	char pathVs[512];
	sprintf(pathVs,"%s%s",direct,nameVs);
	ierr = IGAWriteVec(igaVs,Vs0,pathVs);CHKERRQ(ierr);	
//

//System for L2 projection of div(\partial_S\psi)
	PetscPrintf(PETSC_COMM_WORLD,"\nSystem for Projection of div (\\partial_{S}\\psi) \n\n");
	T=time(NULL);
	tm=*localtime(&T);
	PetscPrintf(PETSC_COMM_WORLD,"Current time is %02d:%02d:%02d \n",tm.tm_hour,tm.tm_min,tm.tm_sec);
	IGA iga_div_proj;
	ierr = IGACreate(PETSC_COMM_WORLD,&iga_div_proj);CHKERRQ(ierr);
	ierr = IGASetDim(iga_div_proj,2);CHKERRQ(ierr);													//Spatial dimension of the problem
	ierr = IGASetDof(iga_div_proj,4);CHKERRQ(ierr);													//Number of degrees of freedom, per node
	ierr = IGASetOrder(iga_div_proj,1);CHKERRQ(ierr);													//Number of spatial derivatives to calculate
	ierr = IGASetFromOptions(iga_div_proj);CHKERRQ(ierr);												//Note: The order (or degree) of the shape functions is given by the mesh!
	ierr = IGARead(iga_div_proj,"./geometry3.dat");CHKERRQ(ierr);
	
	for (dir=0; dir<2; dir++)
	{
		ierr = IGASetRuleType(iga_div_proj,dir,IGA_RULE_LEGENDRE);CHKERRQ(ierr);
		ierr = IGASetRuleSize(iga_div_proj,dir,6);CHKERRQ(ierr);
	}
	ierr = IGASetUp(iga_div_proj);CHKERRQ(ierr);

	for (dir=0; dir<2; dir++) 
	{
		for (side=0; side<2; side++) 
		{
			//ierr = IGASetBoundaryValue(iga,dir,side,dof,0.0);CHKERRQ(ierr);					// Dirichlet boundary conditions
			ierr = IGASetBoundaryForm(iga_div_proj,dir,side,PETSC_TRUE);CHKERRQ(ierr);  				// Neumann boundary conditions
		}
	}

	Mat K_div_proj;
	Vec div_proj,F_div_proj;

	ierr = IGACreateMat(iga_div_proj,&K_div_proj);CHKERRQ(ierr);
	ierr = IGACreateVec(iga_div_proj,&div_proj);CHKERRQ(ierr);
	ierr = IGACreateVec(iga_div_proj,&F_div_proj);CHKERRQ(ierr);

	IGAPoint		point_div_proj;																	//point
	IGAElement		elem_div_proj;																	//element
	PetscReal		*Kloc_div_proj,*Floc_div_proj;													//AA y BB
	PetscReal		*Kpoint_div_proj,*Fpoint_div_proj;												//KKK y FFF
	const PetscReal	*arrayZ0_div_proj,*arrayChi0_div_proj,*arrayS0_div_proj,*arrayZS_div_proj;		//arrayU
	Vec				localZ0_div_proj,localChi0_div_proj,localS0_div_proj,localZS_div_proj;			//localU
	PetscReal		*Chi0_div_proj,*Z0_div_proj,*S0_div_proj,*ZS_div_proj;							//U

  	IGAFormSystem	wtf_div_proj;
 	void			*wtf2_div_proj;

 	KSP ksp_div_proj;
	ierr = IGACreateKSP(iga_div_proj,&ksp_div_proj);CHKERRQ(ierr);

	// Get local vectors Z0 and Chi0 and arrays
	ierr = IGAGetLocalVecArray(igachiUp,chiUp0,&localChi0_div_proj,&arrayChi0_div_proj);CHKERRQ(ierr);
	ierr = IGAGetLocalVecArray(igaZ0,Z0,&localZ0_div_proj,&arrayZ0_div_proj);CHKERRQ(ierr);
	ierr = IGAGetLocalVecArray(igaS,s0,&localS0_div_proj,&arrayS0_div_proj);CHKERRQ(ierr);
	ierr = IGAGetLocalVecArray(igaZS,ZS0,&localZS_div_proj,&arrayZS_div_proj);CHKERRQ(ierr);

	// Element loop
	ierr = IGABeginElement(iga_div_proj,&elem_div_proj);CHKERRQ(ierr);
	ierr = IGABeginElement(igaZ0,&elemZ0);CHKERRQ(ierr);
	ierr = IGABeginElement(igachiUp,&elemchiUp);CHKERRQ(ierr);
	ierr = IGABeginElement(igaS,&elemS);CHKERRQ(ierr);
	ierr = IGABeginElement(igaZS,&elemZS);CHKERRQ(ierr);

	ierr = PetscPrintf(PETSC_COMM_WORLD,"Hola1 \n");

	while (IGANextElement(iga_div_proj,elem_div_proj)) 
	{
		IGANextElement(igaZ0,elemZ0);
		IGANextElement(igachiUp,elemchiUp);
		IGANextElement(igaS,elemS);
		IGANextElement(igaZS,elemZS);

		ierr = IGAElementGetWorkMat(elem_div_proj,&Kloc_div_proj);CHKERRQ(ierr);
		ierr = IGAElementGetWorkVec(elem_div_proj,&Floc_div_proj);CHKERRQ(ierr);
		ierr = IGAElementGetValues(elemZ0,arrayZ0_div_proj,&Z0_div_proj);CHKERRQ(ierr);
		ierr = IGAElementGetValues(elemchiUp,arrayChi0_div_proj,&Chi0_div_proj);CHKERRQ(ierr);
		ierr = IGAElementGetValues(elemS,arrayS0_div_proj,&S0_div_proj);CHKERRQ(ierr);
		ierr = IGAElementGetValues(elemZS,arrayZS_div_proj,&ZS_div_proj);CHKERRQ(ierr);

		// FormSystem loop
		while (IGAElementNextFormSystem(elem_div_proj,&wtf_div_proj,&wtf2_div_proj)) 
		{
		// Quadrature loop
			ierr = IGAElementBeginPoint(elem_div_proj,&point_div_proj);CHKERRQ(ierr);
			ierr = IGAElementBeginPoint(elemZ0,&pointZ0);CHKERRQ(ierr);
			ierr = IGAElementBeginPoint(elemchiUp,&pointchiUp);CHKERRQ(ierr);
			ierr = IGAElementBeginPoint(elemS,&pointS);CHKERRQ(ierr);
			ierr = IGAElementBeginPoint(elemZS,&pointZS);CHKERRQ(ierr);

			while (IGAElementNextPoint(elem_div_proj,point_div_proj))
			{
				if(point_div_proj->atboundary==1)
				{
					//IGAElementNextPoint(elemZ0,pointZ0);
					//IGAElementNextPoint(elemchiUp,pointchiUp);
				}

				if(point_div_proj->atboundary==0 && pointZ0->atboundary==0 && pointchiUp->atboundary==0 && pointS->atboundary==0)
				{
					IGAElementNextPoint(elemZ0,pointZ0);
					IGAElementNextPoint(elemchiUp,pointchiUp);
					IGAElementNextPoint(elemS,pointS);
					IGAElementNextPoint(elemZS,pointZS);

					ierr = IGAPointGetWorkMat(point_div_proj,&Kpoint_div_proj);CHKERRQ(ierr);
					ierr = IGAPointGetWorkVec(point_div_proj,&Fpoint_div_proj);CHKERRQ(ierr);
					
					ierr = ProjectDiv(point_div_proj,pointchiUp,pointZ0,pointS,pointZS,Kpoint_div_proj,Fpoint_div_proj,Chi0_div_proj,Z0_div_proj,S0_div_proj,ZS_div_proj,NULL);CHKERRQ(ierr);
					ierr = IGAPointAddMat(point_div_proj,Kpoint_div_proj,Kloc_div_proj);CHKERRQ(ierr);
					ierr = IGAPointAddVec(point_div_proj,Fpoint_div_proj,Floc_div_proj);CHKERRQ(ierr);
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
			while (pointS->index != -1)
			{
				IGAElementNextPoint(elemS,pointS);
			}
			while (pointZS->index != -1)
			{
				IGAElementNextPoint(elemZS,pointZS);
			}
			ierr = IGAElementEndPoint(elem_div_proj,&point_div_proj);CHKERRQ(ierr);
			ierr = IGAElementEndPoint(elemZ0,&pointZ0);CHKERRQ(ierr);
			ierr = IGAElementEndPoint(elemchiUp,&pointchiUp);CHKERRQ(ierr);
			ierr = IGAElementEndPoint(elemS,&pointS);CHKERRQ(ierr);
			ierr = IGAElementEndPoint(elemZS,&pointZS);CHKERRQ(ierr);
		}

		//ierr = IGAElementFixSystem(elemStress,KlocStress,FlocStress);CHKERRQ(ierr);					//This sets Dirichlet condition ¿? (Yes, this applies the conditions from IGASetBoundaryValue)
		ierr = IGAElementAssembleMat(elem_div_proj,Kloc_div_proj,K_div_proj);CHKERRQ(ierr);
		ierr = IGAElementAssembleVec(elem_div_proj,Floc_div_proj,F_div_proj);CHKERRQ(ierr);
	}
	ierr = PetscPrintf(PETSC_COMM_WORLD,"Hola2 \n");
	IGANextElement(igaZ0,elemZ0);
	IGANextElement(igachiUp,elemchiUp);
	IGANextElement(igaS,elemS);
	IGANextElement(igaZS,elemZS);

	ierr = PetscPrintf(PETSC_COMM_WORLD,"Hola3 \n");

	ierr = IGAEndElement(iga_div_proj,&elem_div_proj);CHKERRQ(ierr);
	ierr = IGAEndElement(igaZ0,&elemZ0);CHKERRQ(ierr);
	ierr = IGAEndElement(igachiUp,&elemchiUp);CHKERRQ(ierr);
	ierr = IGAEndElement(igaS,&elemS);CHKERRQ(ierr);
	ierr = IGAEndElement(igaZS,&elemZS);CHKERRQ(ierr);

	ierr = PetscPrintf(PETSC_COMM_WORLD,"Hola4 \n");

	// Restore local vectors u, Z0, Chi0 and arrays
	ierr = IGARestoreLocalVecArray(igaZ0,Z0,&localZ0_div_proj,&arrayZ0_div_proj);CHKERRQ(ierr);
	ierr = IGARestoreLocalVecArray(igachiUp,chiUp0,&localChi0_div_proj,&arrayChi0_div_proj);CHKERRQ(ierr);
	ierr = IGARestoreLocalVecArray(igaS,s0,&localS0_div_proj,&arrayS0_div_proj);CHKERRQ(ierr);
	ierr = IGARestoreLocalVecArray(igaZS,ZS0,&localZS_div_proj,&arrayZS_div_proj);CHKERRQ(ierr);

	ierr = MatAssemblyBegin(K_div_proj,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyEnd  (K_div_proj,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

	ierr = VecAssemblyBegin(F_div_proj);CHKERRQ(ierr);
	ierr = VecAssemblyEnd  (F_div_proj);CHKERRQ(ierr);

	ierr = KSPSetOperators(ksp_div_proj,K_div_proj,K_div_proj);CHKERRQ(ierr);
	ierr = KSPSetFromOptions(ksp_div_proj);CHKERRQ(ierr);
	PC pc_div_proj;
	ierr = KSPGetPC(ksp_div_proj,&pc_div_proj); CHKERRQ(ierr);
	ierr = PCSetType(pc_div_proj,PCSOR); CHKERRQ(ierr);
	ierr = PCSORSetIterations(pc_div_proj,1,1);CHKERRQ(ierr);
	ierr = KSPSetType(ksp_div_proj,KSPFCG);
	ierr = KSPSetTolerances(ksp_div_proj,1e-20,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);CHKERRQ(ierr);
	ierr = KSPSolve(ksp_div_proj,F_div_proj,div_proj);CHKERRQ(ierr);

	ierr = KSPDestroy(&ksp_div_proj);CHKERRQ(ierr);
	ierr = MatDestroy(&K_div_proj);CHKERRQ(ierr);
	ierr = VecDestroy(&F_div_proj);CHKERRQ(ierr);

	char name_div_proj[]="/DivPsi_S-2d-0.dat";
	char path_div_proj[512];
	sprintf(path_div_proj,"%s%s",direct,name_div_proj);
	ierr = IGAWriteVec(iga_div_proj,div_proj,path_div_proj);CHKERRQ(ierr);	
//

/* Not adapted yet
//System for L2 projection of classic stress (C*Ue)
	PetscPrintf(PETSC_COMM_WORLD,"\nSystem for Classic Stress starting \n\n");
	T=time(NULL);
	tm=*localtime(&T);
	PetscPrintf(PETSC_COMM_WORLD,"Current time is %02d:%02d:%02d \n",tm.tm_hour,tm.tm_min,tm.tm_sec);
	IGA igaClassicStress;
	ierr = IGACreate(PETSC_COMM_WORLD,&igaClassicStress);CHKERRQ(ierr);
	ierr = IGASetDim(igaClassicStress,2);CHKERRQ(ierr);													//Spatial dimension of the problem
	ierr = IGASetDof(igaClassicStress,4);CHKERRQ(ierr);													//Number of degrees of freedom, per node
	ierr = IGASetOrder(igaClassicStress,2);CHKERRQ(ierr);													//Number of spatial derivatives to calculate
	ierr = IGASetFromOptions(igaClassicStress);CHKERRQ(ierr);												//Note: The order (or degree) of the shape functions is given by the mesh!
	ierr = IGARead(igaClassicStress,"./geometry3.dat");CHKERRQ(ierr);
	
	for (dir=0; dir<2; dir++)
	{
		ierr = IGASetRuleType(igaClassicStress,dir,IGA_RULE_LEGENDRE);CHKERRQ(ierr);
		ierr = IGASetRuleSize(igaClassicStress,dir,6);CHKERRQ(ierr);
	}
	ierr = IGASetUp(igaClassicStress);CHKERRQ(ierr);

	for (dir=0; dir<2; dir++) 
	{
		for (side=0; side<2; side++) 
		{
			ierr = IGASetBoundaryForm(igaClassicStress,dir,side,PETSC_TRUE);CHKERRQ(ierr);  				// Neumann boundary conditions
		}
	}

	Mat KClassicStress;
	Vec classicSigma0,FClassicStress;

	ierr = IGACreateMat(igaClassicStress,&KClassicStress);CHKERRQ(ierr);
	ierr = IGACreateVec(igaClassicStress,&classicSigma0);CHKERRQ(ierr);
	ierr = IGACreateVec(igaClassicStress,&FClassicStress);CHKERRQ(ierr);

	IGAPoint		pointClassicStress;															//point
	IGAElement		elemClassicStress;																//element
	PetscReal		*KlocClassicStress,*FlocClassicStress;												//AA y BB
	PetscReal		*KpointClassicStress,*FpointClassicStress;											//KKK y FFF
	const PetscReal	*arrayZ0ClassicStress,*arrayChi0ClassicStress;			//arrayU
	Vec				localZ0ClassicStress,localChi0ClassicStress;				//localU
	PetscReal		*Chi0ClassicStress,*Z0ClassicStress;											//U

  	IGAFormSystem	wtfClassicStress;
 	void			*wtf2ClassicStress;

 	KSP kspClassicStress;
	ierr = IGACreateKSP(igaClassicStress,&kspClassicStress);CHKERRQ(ierr);

	// Get local vectors Z0 and Chi0 and arrays
	ierr = IGAGetLocalVecArray(igachiUp,chiUp0,&localChi0ClassicStress,&arrayChi0ClassicStress);CHKERRQ(ierr);
	ierr = IGAGetLocalVecArray(igaZ0,Z0,&localZ0ClassicStress,&arrayZ0ClassicStress);CHKERRQ(ierr);

	// Element loop
	ierr = IGABeginElement(igaClassicStress,&elemClassicStress);CHKERRQ(ierr);
	ierr = IGABeginElement(igaZ0,&elemZ0);CHKERRQ(ierr);
	ierr = IGABeginElement(igachiUp,&elemchiUp);CHKERRQ(ierr);

	while (IGANextElement(igaClassicStress,elemClassicStress)) 
	{
		IGANextElement(igaZ0,elemZ0);
		IGANextElement(igachiUp,elemchiUp);

		ierr = IGAElementGetWorkMat(elemClassicStress,&KlocClassicStress);CHKERRQ(ierr);
		ierr = IGAElementGetWorkVec(elemClassicStress,&FlocClassicStress);CHKERRQ(ierr);
		ierr = IGAElementGetValues(elemZ0,arrayZ0ClassicStress,&Z0ClassicStress);CHKERRQ(ierr);
		ierr = IGAElementGetValues(elemchiUp,arrayChi0ClassicStress,&Chi0ClassicStress);CHKERRQ(ierr);

		// FormSystem loop
		while (IGAElementNextFormSystem(elemClassicStress,&wtfClassicStress,&wtf2ClassicStress)) 
		{
		// Quadrature loop
			ierr = IGAElementBeginPoint(elemClassicStress,&pointClassicStress);CHKERRQ(ierr);
			ierr = IGAElementBeginPoint(elemZ0,&pointZ0);CHKERRQ(ierr);
			ierr = IGAElementBeginPoint(elemchiUp,&pointchiUp);CHKERRQ(ierr);

			while (IGAElementNextPoint(elemClassicStress,pointClassicStress))
			{
				if(pointClassicStress->atboundary==1)
				{
					//IGAElementNextPoint(elemZ0,pointZ0);
					//IGAElementNextPoint(elemchiUp,pointchiUp);
				}

				if(pointZ0->atboundary==1)
				{
					PetscPrintf(PETSC_COMM_WORLD,"Hola1 en ClassicStress \n");
				}

				if(pointClassicStress->atboundary==0 && pointZ0->atboundary==0 && pointchiUp->atboundary==0)
				{
					IGAElementNextPoint(elemZ0,pointZ0);
					IGAElementNextPoint(elemchiUp,pointchiUp);

					ierr = IGAPointGetWorkMat(pointClassicStress,&KpointClassicStress);CHKERRQ(ierr);
					ierr = IGAPointGetWorkVec(pointClassicStress,&FpointClassicStress);CHKERRQ(ierr);
					//PetscErrorCode Stress(IGAPoint p,IGAPoint pChi,IGAPoint pZu,PetscReal *K,PetscReal *F, PetscReal *Chi,PetscReal *Zu,void *ctx)
					ierr = ClassicStress(pointClassicStress,pointchiUp,pointZ0,KpointClassicStress,FpointClassicStress,Chi0ClassicStress,Z0ClassicStress,NULL);CHKERRQ(ierr);
					ierr = IGAPointAddMat(pointClassicStress,KpointClassicStress,KlocClassicStress);CHKERRQ(ierr);
					ierr = IGAPointAddVec(pointClassicStress,FpointClassicStress,FlocClassicStress);CHKERRQ(ierr);
				}
			}
			if (pointZ0->index != -1)
			{
				IGAElementNextPoint(elemZ0,pointZ0);
			}
			while (pointchiUp->index != -1)
			{
				IGAElementNextPoint(elemchiUp,pointchiUp);
			}
			ierr = IGAElementEndPoint(elemClassicStress,&pointClassicStress);CHKERRQ(ierr);
			ierr = IGAElementEndPoint(elemZ0,&pointZ0);CHKERRQ(ierr);
			ierr = IGAElementEndPoint(elemchiUp,&pointchiUp);CHKERRQ(ierr);
		}

		//ierr = IGAElementFixSystem(elemStress,KlocStress,FlocStress);CHKERRQ(ierr);					//This sets Dirichlet condition ¿? (Yes, this applies the conditions from IGASetBoundaryValue)
		ierr = IGAElementAssembleMat(elemClassicStress,KlocClassicStress,KClassicStress);CHKERRQ(ierr);
		ierr = IGAElementAssembleVec(elemClassicStress,FlocClassicStress,FClassicStress);CHKERRQ(ierr);
	}
	IGANextElement(igaZ0,elemZ0);
	IGANextElement(igachiUp,elemchiUp);

	ierr = IGAEndElement(igaClassicStress,&elemClassicStress);CHKERRQ(ierr);
	ierr = IGAEndElement(igaZ0,&elemZ0);CHKERRQ(ierr);
	ierr = IGAEndElement(igachiUp,&elemchiUp);CHKERRQ(ierr);

	// Restore local vectors u, Z0, Chi0 and arrays
	ierr = IGARestoreLocalVecArray(igaZ0,Z0,&localZ0ClassicStress,&arrayZ0ClassicStress);CHKERRQ(ierr);
	ierr = IGARestoreLocalVecArray(igachiUp,chiUp0,&localChi0ClassicStress,&arrayChi0ClassicStress);CHKERRQ(ierr);
	ierr = MatAssemblyBegin(KClassicStress,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyEnd  (KClassicStress,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

	ierr = VecAssemblyBegin(FClassicStress);CHKERRQ(ierr);
	ierr = VecAssemblyEnd  (FClassicStress);CHKERRQ(ierr);

	ierr = KSPSetOperators(kspClassicStress,KClassicStress,KClassicStress);CHKERRQ(ierr);
	PC pcClassicStress;
	ierr = KSPGetPC(kspClassicStress,&pcClassicStress); CHKERRQ(ierr);
	ierr = PCSetType(pcClassicStress,PCSOR); CHKERRQ(ierr);
	//ierr = PCFactorSetMatSolverType(pcStress,MATSOLVERMUMPS); CHKERRQ(ierr);
	ierr = KSPSetFromOptions(kspClassicStress);CHKERRQ(ierr);
	ierr = KSPSetTolerances(kspClassicStress,1.0e-12,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);CHKERRQ(ierr);
	ierr = KSPSolve(kspClassicStress,FClassicStress,classicSigma0);CHKERRQ(ierr);

	ierr = KSPDestroy(&kspClassicStress);CHKERRQ(ierr);
	ierr = MatDestroy(&KClassicStress);CHKERRQ(ierr);
	ierr = VecDestroy(&FClassicStress);CHKERRQ(ierr);

	char nameClassicStress[]="/classicSigma-2d-0.dat";
	char pathClassicStress[512];
	sprintf(pathClassicStress,"%s%s",direct,nameClassicStress);
	ierr = IGAWriteVec(igaClassicStress,classicSigma0,pathClassicStress);CHKERRQ(ierr);

	ierr = VecDestroy(&classicSigma0);CHKERRQ(ierr);
	ierr = IGADestroy(&igaClassicStress);CHKERRQ(ierr);
//
*/

// Adapted, but not debugged
//System for L2 projection of couple stress
	PetscPrintf(PETSC_COMM_WORLD,"\nSystem for CoupleStress starting \n\n");
	T=time(NULL);
	tm=*localtime(&T);
	PetscPrintf(PETSC_COMM_WORLD,"Current time is %02d:%02d:%02d \n",tm.tm_hour,tm.tm_min,tm.tm_sec);
	IGA igaCS;
	ierr = IGACreate(PETSC_COMM_WORLD,&igaCS);CHKERRQ(ierr);
	ierr = IGASetDim(igaCS,2);CHKERRQ(ierr);													//Spatial dimension of the problem
	ierr = IGASetDof(igaCS,2);CHKERRQ(ierr);													//Number of degrees of freedom, per node
	ierr = IGASetOrder(igaCS,2);CHKERRQ(ierr);													//Number of spatial derivatives to calculate
	ierr = IGASetFromOptions(igaCS);CHKERRQ(ierr);												//Note: The order (or degree) of the shape functions is given by the mesh!
	ierr = IGARead(igaCS,"./geometry3.dat");CHKERRQ(ierr);
	
	for (dir=0; dir<2; dir++)
	{
		ierr = IGASetRuleType(igaCS,dir,IGA_RULE_LEGENDRE);CHKERRQ(ierr);
		ierr = IGASetRuleSize(igaCS,dir,6);CHKERRQ(ierr);
	}
	ierr = IGASetUp(igaCS);CHKERRQ(ierr);

	for (dir=0; dir<2; dir++) 
	{
		for (side=0; side<2; side++) 
		{
			//ierr = IGASetBoundaryValue(iga,dir,side,dof,0.0);CHKERRQ(ierr);					// Dirichlet boundary conditions
			ierr = IGASetBoundaryForm(igaCS,dir,side,PETSC_TRUE);CHKERRQ(ierr);  				// Neumann boundary conditions
		}
	}

	Mat KCStress;
	Vec lambda0,FCStress;

	ierr = IGACreateMat(igaCS,&KCStress);CHKERRQ(ierr);
	ierr = IGACreateVec(igaCS,&lambda0);CHKERRQ(ierr);
	ierr = IGACreateVec(igaCS,&FCStress);CHKERRQ(ierr);

	IGAPoint		pointCS;														//point
	IGAElement		elemCS;															//element
	PetscReal		*KlocCS,*FlocCS;												//AA y BB
	PetscReal		*KpointCS,*FpointCS;											//KKK y FFF
	const PetscReal	*arrayZ0CS,*arrayChi0CS,*arrayZSCS,*arrayS0CS;					//arrayU
	Vec				localZ0CS,localChi0CS,localZSCS,localS0CS;						//localU
	PetscReal		*Chi0CS,*Z0CS,*ZSCS,*S0CS;										//U

  	IGAFormSystem	wtfCS;
 	void			*wtf2CS;

 	KSP kspCS;
	ierr = IGACreateKSP(igaCS,&kspCS);CHKERRQ(ierr);

	// Get local vectors Z0 and Chi0 and arrays
	ierr = IGAGetLocalVecArray(igachiUp,chiUp0,&localChi0CS,&arrayChi0CS);CHKERRQ(ierr);
	ierr = IGAGetLocalVecArray(igaZ0,Z0,&localZ0CS,&arrayZ0CS);CHKERRQ(ierr);
	ierr = IGAGetLocalVecArray(igaZS,ZS0,&localZSCS,&arrayZSCS);CHKERRQ(ierr);
	ierr = IGAGetLocalVecArray(igaS,s0,&localS0CS,&arrayS0CS);CHKERRQ(ierr);

	// Element loop
	ierr = IGABeginElement(igaCS,&elemCS);CHKERRQ(ierr);
	ierr = IGABeginElement(igaZ0,&elemZ0);CHKERRQ(ierr);
	ierr = IGABeginElement(igachiUp,&elemchiUp);CHKERRQ(ierr);
	ierr = IGABeginElement(igaZS,&elemZS);CHKERRQ(ierr);
	ierr = IGABeginElement(igaS,&elemS);CHKERRQ(ierr);

	while (IGANextElement(igaCS,elemCS)) 				
	{
		IGANextElement(igaZ0,elemZ0);
		IGANextElement(igachiUp,elemchiUp);
		IGANextElement(igaZS,elemZS);
		IGANextElement(igaS,elemS);

		ierr = IGAElementGetWorkMat(elemCS,&KlocCS);CHKERRQ(ierr);
		ierr = IGAElementGetWorkVec(elemCS,&FlocCS);CHKERRQ(ierr);
		ierr = IGAElementGetValues(elemZ0,arrayZ0CS,&Z0CS);CHKERRQ(ierr);
		ierr = IGAElementGetValues(elemchiUp,arrayChi0CS,&Chi0CS);CHKERRQ(ierr);
		ierr = IGAElementGetValues(elemZS,arrayZSCS,&ZSCS);CHKERRQ(ierr);
		ierr = IGAElementGetValues(elemS,arrayS0CS,&S0CS);CHKERRQ(ierr);

		// FormSystem loop
		while (IGAElementNextFormSystem(elemCS,&wtfCS,&wtf2CS))
		{
		// Quadrature loop
			ierr = IGAElementBeginPoint(elemCS,&pointCS);CHKERRQ(ierr);
			ierr = IGAElementBeginPoint(elemZ0,&pointZ0);CHKERRQ(ierr);
			ierr = IGAElementBeginPoint(elemchiUp,&pointchiUp);CHKERRQ(ierr);
			ierr = IGAElementBeginPoint(elemZS,&pointZS);CHKERRQ(ierr);
			ierr = IGAElementBeginPoint(elemS,&pointS);CHKERRQ(ierr);

			while (IGAElementNextPoint(elemCS,pointCS))
			{
				if(pointCS->atboundary==1)
				{
					//delete this block??
				}

				if(pointZ0->atboundary==1)
				{
					PetscPrintf(PETSC_COMM_WORLD,"Hola1 en CoupleStress \n");
				}

				if(pointCS->atboundary==0 && pointZ0->atboundary==0 && pointchiUp->atboundary==0 && pointZS->atboundary==0 && pointS->atboundary==0)
				{
					IGAElementNextPoint(elemZ0,pointZ0);
					IGAElementNextPoint(elemchiUp,pointchiUp);
					IGAElementNextPoint(elemZS,pointZS);
					IGAElementNextPoint(elemS,pointS);

					ierr = IGAPointGetWorkMat(pointCS,&KpointCS);CHKERRQ(ierr);
					ierr = IGAPointGetWorkVec(pointCS,&FpointCS);CHKERRQ(ierr);
					//	   CoupleStress(IGAPoint p,IGAPoint pChi,IGAPoint pZu,IGAPoint pZS,IGAPoint pS, PetscReal *K,PetscReal *F, PetscReal *Chi,PetscReal *Zu,PetscReal *ZS,PetscReal *S,void *ctx)
					ierr = CoupleStress(pointCS,pointchiUp,pointZ0,pointZS,pointS,KpointCS,FpointCS,Chi0CS,Z0CS,ZSCS,S0CS,NULL);CHKERRQ(ierr);
					ierr = IGAPointAddMat(pointCS,KpointCS,KlocCS);CHKERRQ(ierr);
					ierr = IGAPointAddVec(pointCS,FpointCS,FlocCS);CHKERRQ(ierr);
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
			while (pointZS->index != -1)
			{
				IGAElementNextPoint(elemZS,pointZS);
			}
			while (pointS->index != -1)
			{
				IGAElementNextPoint(elemS,pointS);
			}

			ierr = IGAElementEndPoint(elemCS,&pointCS);CHKERRQ(ierr);
			ierr = IGAElementEndPoint(elemZ0,&pointZ0);CHKERRQ(ierr);
			ierr = IGAElementEndPoint(elemchiUp,&pointchiUp);CHKERRQ(ierr);
			ierr = IGAElementEndPoint(elemZS,&pointZS);CHKERRQ(ierr);
			ierr = IGAElementEndPoint(elemS,&pointS);CHKERRQ(ierr);
		}

		//ierr = IGAElementFixSystem(elemStress,KlocStress,FlocStress);CHKERRQ(ierr);					//This sets Dirichlet condition ¿? (Yes, this applies the conditions from IGASetBoundaryValue)
		ierr = IGAElementAssembleMat(elemCS,KlocCS,KCStress);CHKERRQ(ierr);
		ierr = IGAElementAssembleVec(elemCS,FlocCS,FCStress);CHKERRQ(ierr);
	}
	IGANextElement(igaZ0,elemZ0);
	IGANextElement(igachiUp,elemchiUp);
	IGANextElement(igaZS,elemZS);
	IGANextElement(igaS,elemS);

	ierr = IGAEndElement(igaCS,&elemCS);CHKERRQ(ierr);
	ierr = IGAEndElement(igaZ0,&elemZ0);CHKERRQ(ierr);
	ierr = IGAEndElement(igachiUp,&elemchiUp);CHKERRQ(ierr);
	ierr = IGAEndElement(igaZS,&elemZS);CHKERRQ(ierr);
	ierr = IGAEndElement(igaS,&elemS);CHKERRQ(ierr);

	// Restore local vectors u, Z0, Chi0 and arrays
	ierr = IGARestoreLocalVecArray(igaZ0,Z0,&localZ0CS,&arrayZ0CS);CHKERRQ(ierr);
	ierr = IGARestoreLocalVecArray(igachiUp,chiUp0,&localChi0CS,&arrayChi0CS);CHKERRQ(ierr);
	ierr = IGARestoreLocalVecArray(igaZS,ZS0,&localZSCS,&arrayZSCS);CHKERRQ(ierr);
	ierr = IGARestoreLocalVecArray(igaS,s0,&localS0CS,&arrayS0CS);CHKERRQ(ierr);

	ierr = MatAssemblyBegin(KCStress,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyEnd  (KCStress,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

	ierr = VecAssemblyBegin(FCStress);CHKERRQ(ierr);
	ierr = VecAssemblyEnd  (FCStress);CHKERRQ(ierr);

	ierr = KSPSetOperators(kspCS,KCStress,KCStress);CHKERRQ(ierr);
	PC pcCS;
	ierr = KSPGetPC(kspCS,&pcCS); CHKERRQ(ierr);
	ierr = PCSetType(pcCS,PCSOR); CHKERRQ(ierr);
	ierr = KSPSetFromOptions(kspCS);CHKERRQ(ierr);
	ierr = KSPSetType(kspCS,KSPFCG);
	ierr = KSPSetTolerances(kspCS,1.0e-16,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);CHKERRQ(ierr);
	ierr = KSPSolve(kspCS,FCStress,lambda0);CHKERRQ(ierr);

	ierr = KSPDestroy(&kspCS);CHKERRQ(ierr);
	ierr = MatDestroy(&KCStress);CHKERRQ(ierr);
	ierr = VecDestroy(&FCStress);CHKERRQ(ierr);

	char nameCStress[]="/lambda-2d-0.dat";
	char pathCStress[512];
	sprintf(pathCStress,"%s%s",direct,nameCStress);
	ierr = IGAWriteVec(igaCS,lambda0,pathCStress);CHKERRQ(ierr);

	ierr = VecDestroy(&lambda0);CHKERRQ(ierr);
	ierr = IGADestroy(&igaCS);CHKERRQ(ierr);
//

// Adapted but not debugged
//System for L2 projection of Energy Density
	PetscPrintf(PETSC_COMM_WORLD,"\nSystem for Energy Density starting \n\n");
	T=time(NULL);
	tm=*localtime(&T);
	PetscPrintf(PETSC_COMM_WORLD,"Current time is %02d:%02d:%02d \n",tm.tm_hour,tm.tm_min,tm.tm_sec);
	IGA igaED;
	ierr = IGACreate(PETSC_COMM_WORLD,&igaED);CHKERRQ(ierr);
	ierr = IGASetDim(igaED,2);CHKERRQ(ierr);													//Spatial dimension of the problem
	ierr = IGASetDof(igaED,1);CHKERRQ(ierr);													//Number of degrees of freedom, per node
	ierr = IGASetOrder(igaED,2);CHKERRQ(ierr);												//Number of spatial derivatives to calculate
	ierr = IGASetFromOptions(igaED);CHKERRQ(ierr);											//Note: The order (or degree) of the shape functions is given by the mesh!
	ierr = IGARead(igaED,"./geometry3.dat");CHKERRQ(ierr);
	
	for (dir=0; dir<2; dir++)
	{
		ierr = IGASetRuleType(igaED,dir,IGA_RULE_LEGENDRE);CHKERRQ(ierr);
		ierr = IGASetRuleSize(igaED,dir,6);CHKERRQ(ierr);
	}
	ierr = IGASetUp(igaED);CHKERRQ(ierr);

	for (dir=0; dir<2; dir++) 
	{
		for (side=0; side<2; side++) 
		{
			//ierr = IGASetBoundaryValue(iga,dir,side,dof,0.0);CHKERRQ(ierr);					// Dirichlet boundary conditions
			ierr = IGASetBoundaryForm(igaED,dir,side,PETSC_TRUE);CHKERRQ(ierr);  				// Neumann boundary conditions
		}
	}

	Mat KED;
	Vec ed0,FED;

	ierr = IGACreateMat(igaED,&KED);CHKERRQ(ierr);
	ierr = IGACreateVec(igaED,&ed0);CHKERRQ(ierr);
	ierr = IGACreateVec(igaED,&FED);CHKERRQ(ierr);

	IGAPoint		pointED;															//point
	IGAElement		elemED;																//element
	PetscReal		*KlocED,*FlocED;													//AA y BB
	PetscReal		*KpointED,*FpointED;												//KKK y FFF
	const PetscReal	*arrayChi0ED,*arrayZ0ED,*arrayS0ED,*arrayZSED;						//arrayU
	Vec				localChi0ED,localZ0ED,localS0ED,localZSED;							//localU
	PetscReal		*chi0ED,*Z0ED,*S0ED,*ZSED;											//U

  	IGAFormSystem	wtfED;
 	void			*wtf2ED;

 	KSP kspED;
	ierr = IGACreateKSP(igaED,&kspED);CHKERRQ(ierr);

	// Get local vectors Z0 and Chi0 and arrays
	ierr = IGAGetLocalVecArray(igachiUp,chiUp0,&localChi0ED,&arrayChi0ED);CHKERRQ(ierr);
	ierr = IGAGetLocalVecArray(igaZ0,Z0,&localZ0ED,&arrayZ0ED);CHKERRQ(ierr);
	ierr = IGAGetLocalVecArray(igaS,s0,&localS0ED,&arrayS0ED);CHKERRQ(ierr);
	ierr = IGAGetLocalVecArray(igaZS,ZS0,&localZSED,&arrayZSED);CHKERRQ(ierr);

	// Element loop
	ierr = IGABeginElement(igaED,&elemED);CHKERRQ(ierr);
	ierr = IGABeginElement(igaZ0,&elemZ0);CHKERRQ(ierr);
	ierr = IGABeginElement(igachiUp,&elemchiUp);CHKERRQ(ierr);
	ierr = IGABeginElement(igaS,&elemS);CHKERRQ(ierr);
	ierr = IGABeginElement(igaZS,&elemZS);CHKERRQ(ierr);

	while (IGANextElement(igaED,elemED)) 				
	{
		IGANextElement(igachiUp,elemchiUp);
		IGANextElement(igaZ0,elemZ0);
		IGANextElement(igaS,elemS);
		IGANextElement(igaZS,elemZS);

		ierr = IGAElementGetWorkMat(elemED,&KlocED);CHKERRQ(ierr);
		ierr = IGAElementGetWorkVec(elemED,&FlocED);CHKERRQ(ierr);
		ierr = IGAElementGetValues(elemchiUp,arrayChi0ED,&chi0ED);CHKERRQ(ierr);
		ierr = IGAElementGetValues(elemZ0,arrayZ0ED,&Z0ED);CHKERRQ(ierr);
		ierr = IGAElementGetValues(elemS,arrayS0ED,&S0ED);CHKERRQ(ierr);
		ierr = IGAElementGetValues(elemZS,arrayZSED,&ZSED);CHKERRQ(ierr);

		// FormSystem loop
		while (IGAElementNextFormSystem(elemED,&wtfED,&wtf2ED))
		{
		// Quadrature loop
			ierr = IGAElementBeginPoint(elemED,&pointED);CHKERRQ(ierr);
			ierr = IGAElementBeginPoint(elemZ0,&pointZ0);CHKERRQ(ierr);
			ierr = IGAElementBeginPoint(elemchiUp,&pointchiUp);CHKERRQ(ierr);
			ierr = IGAElementBeginPoint(elemS,&pointS);CHKERRQ(ierr);
			ierr = IGAElementBeginPoint(elemZS,&pointZS);CHKERRQ(ierr);

			while (IGAElementNextPoint(elemED,pointED))
			{
				if(pointED->atboundary==0 && pointchiUp->atboundary==0 && pointZ0->atboundary==0)
				{
					IGAElementNextPoint(elemZ0,pointZ0);
					IGAElementNextPoint(elemchiUp,pointchiUp);
					IGAElementNextPoint(elemS,pointS);
					IGAElementNextPoint(elemZS,pointZS);

					ierr = IGAPointGetWorkMat(pointED,&KpointED);CHKERRQ(ierr);
					ierr = IGAPointGetWorkVec(pointED,&FpointED);CHKERRQ(ierr);
					//	   EnergyDensity(IGAPoint p, IGAPoint pChi, IGAPoint pZu, IGAPoint pS, PetscReal *K, PetscReal *F, PetscReal *Chi, PetscReal *Zu, PetscReal *S0, void *ctx)	
					ierr = EnergyDensity(pointED,pointchiUp,pointZ0,pointS,pointZS,KpointED,FpointED,chi0ED,Z0ED,S0ED,ZSED,NULL);CHKERRQ(ierr);
					ierr = IGAPointAddMat(pointED,KpointED,KlocED);CHKERRQ(ierr);
					ierr = IGAPointAddVec(pointED,FpointED,FlocED);CHKERRQ(ierr);
				}
			}
			while (pointchiUp->index != -1)
			{
				IGAElementNextPoint(elemchiUp,pointchiUp);
			}
			while (pointZ0->index != -1)
			{
				IGAElementNextPoint(elemZ0,pointZ0);
			}
			while (pointS->index != -1)
			{
				IGAElementNextPoint(elemS,pointS);
			}
			while (pointZS->index != -1)
			{
				IGAElementNextPoint(elemZS,pointZS);
			}

			ierr = IGAElementEndPoint(elemED,&pointED);CHKERRQ(ierr);
			ierr = IGAElementEndPoint(elemchiUp,&pointchiUp);CHKERRQ(ierr);
			ierr = IGAElementEndPoint(elemZ0,&pointZ0);CHKERRQ(ierr);
			ierr = IGAElementEndPoint(elemS,&pointS);CHKERRQ(ierr);
			ierr = IGAElementEndPoint(elemZS,&pointZS);CHKERRQ(ierr);
		}

		//ierr = IGAElementFixSystem(elemStress,KlocStress,FlocStress);CHKERRQ(ierr);					//This sets Dirichlet condition ¿? (Yes, this applies the conditions from IGASetBoundaryValue)
		ierr = IGAElementAssembleMat(elemED,KlocED,KED);CHKERRQ(ierr);
		ierr = IGAElementAssembleVec(elemED,FlocED,FED);CHKERRQ(ierr);
	}
	IGANextElement(igachiUp,elemchiUp);
	IGANextElement(igaZ0,elemZ0);
	IGANextElement(igaS,elemS);
	IGANextElement(igaZS,elemZS);

	ierr = IGAEndElement(igaED,&elemED);CHKERRQ(ierr);
	ierr = IGAEndElement(igaZ0,&elemZ0);CHKERRQ(ierr);
	ierr = IGAEndElement(igachiUp,&elemchiUp);CHKERRQ(ierr);
	ierr = IGAEndElement(igaS,&elemS);CHKERRQ(ierr);
	ierr = IGAEndElement(igaZS,&elemZS);CHKERRQ(ierr);

	// Restore local vectors u, Z0, Chi0 and arrays
	ierr = IGARestoreLocalVecArray(igachiUp,chiUp0,&localChi0ED,&arrayChi0ED);CHKERRQ(ierr);
	ierr = IGARestoreLocalVecArray(igaZ0,Z0,&localZ0ED,&arrayZ0ED);CHKERRQ(ierr);
	ierr = IGARestoreLocalVecArray(igaS,s0,&localS0ED,&arrayS0ED);CHKERRQ(ierr);
	ierr = IGARestoreLocalVecArray(igaZS,ZS0,&localZSED,&arrayZSED);CHKERRQ(ierr);

	ierr = MatAssemblyBegin(KED,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyEnd  (KED,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

	ierr = VecAssemblyBegin(FED);CHKERRQ(ierr);
	ierr = VecAssemblyEnd  (FED);CHKERRQ(ierr);

	ierr = KSPSetOperators(kspED,KED,KED);CHKERRQ(ierr);
	//PC pcED;
	//ierr = KSPGetPC(kspED,&pcED); CHKERRQ(ierr);
	//ierr = PCSetType(pcED,PCLU); CHKERRQ(ierr);
	//ierr = PCFactorSetMatSolverType(pcED,MATSOLVERMUMPS); CHKERRQ(ierr);
	ierr = KSPSetFromOptions(kspED);CHKERRQ(ierr);
	ierr = KSPSetTolerances(kspED,1.0e-27,1.0e-28,PETSC_DEFAULT,PETSC_DEFAULT);CHKERRQ(ierr);
	ierr = KSPSolve(kspED,FED,ed0);CHKERRQ(ierr);

	ierr = KSPDestroy(&kspED);CHKERRQ(ierr);
	ierr = MatDestroy(&KED);CHKERRQ(ierr);
	ierr = VecDestroy(&FED);CHKERRQ(ierr);

	char nameED[]="/DensityEnergy.dat";
	char pathED[512];
	sprintf(pathED,"%s%s",direct,nameED);
	ierr = IGAWriteVec(igaED,ed0,pathED);CHKERRQ(ierr);	
//

/* Not adapted yet
//System for L2 projection of full stress
	PetscPrintf(PETSC_COMM_WORLD,"\nSystem for full Stress starting \n\n");
	T=time(NULL);
	tm=*localtime(&T);
	PetscPrintf(PETSC_COMM_WORLD,"Current time is %02d:%02d:%02d \n",tm.tm_hour,tm.tm_min,tm.tm_sec);
	IGA igaFullS;
	ierr = IGACreate(PETSC_COMM_WORLD,&igaFullS);CHKERRQ(ierr);
	ierr = IGASetDim(igaFullS,2);CHKERRQ(ierr);													//Spatial dimension of the problem
	ierr = IGASetDof(igaFullS,4);CHKERRQ(ierr);													//Number of degrees of freedom, per node
	ierr = IGASetOrder(igaFullS,2);CHKERRQ(ierr);													//Number of spatial derivatives to calculate
	ierr = IGASetFromOptions(igaFullS);CHKERRQ(ierr);												//Note: The order (or degree) of the shape functions is given by the mesh!
	ierr = IGARead(igaFullS,"./geometry3.dat");CHKERRQ(ierr);
	
	for (dir=0; dir<2; dir++)
	{
		ierr = IGASetRuleType(igaFullS,dir,IGA_RULE_LEGENDRE);CHKERRQ(ierr);
		ierr = IGASetRuleSize(igaFullS,dir,6);CHKERRQ(ierr);
	}
	ierr = IGASetUp(igaFullS);CHKERRQ(ierr);

	for (dir=0; dir<2; dir++) 
	{
		for (side=0; side<2; side++) 
		{
			//ierr = IGASetBoundaryValue(iga,dir,side,dof,0.0);CHKERRQ(ierr);					// Dirichlet boundary conditions
			ierr = IGASetBoundaryForm(igaFullS,dir,side,PETSC_FALSE);CHKERRQ(ierr);  				// Neumann boundary conditions
		}
	}

	Mat KFullS;
	Vec Fs0,FFullS;

	ierr = IGACreateMat(igaFullS,&KFullS);CHKERRQ(ierr);
	ierr = IGACreateVec(igaFullS,&Fs0);CHKERRQ(ierr);
	ierr = IGACreateVec(igaFullS,&FFullS);CHKERRQ(ierr);

	IGAPoint		pointFullS;													//point
	IGAElement		elemFullS;													//element
	PetscReal		*KlocFullS,*FlocFullS;										//AA y BB
	PetscReal		*KpointFullS,*FpointFullS;									//KKK y FFF
	const PetscReal	*arrayZ0FullS,*arrayChi0FullS;								//arrayU
	Vec				localZ0FullS,localChi0FullS;								//localU
	PetscReal		*Chi0FullS,*Z0FullS;										//U

  	IGAFormSystem	wtfFullS;
 	void			*wtf2FullS;

 	KSP kspFullS;
	ierr = IGACreateKSP(igaFullS,&kspFullS);CHKERRQ(ierr);

	// Get local vectors Z0 and Chi0 and arrays
	ierr = IGAGetLocalVecArray(igachiUp,chiUp0,&localChi0FullS,&arrayChi0FullS);CHKERRQ(ierr);
	ierr = IGAGetLocalVecArray(igaZ0,Z0,&localZ0FullS,&arrayZ0FullS);CHKERRQ(ierr);

	// Element loop
	ierr = IGABeginElement(igaFullS,&elemFullS);CHKERRQ(ierr);
	ierr = IGABeginElement(igaZ0,&elemZ0);CHKERRQ(ierr);
	ierr = IGABeginElement(igachiUp,&elemchiUp);CHKERRQ(ierr);

	ierr = PetscPrintf(PETSC_COMM_WORLD,"Antes del ciclo, todos los objetos creados \n");CHKERRQ(ierr);
	while (IGANextElement(igaFullS,elemFullS)) 
	{
		IGANextElement(igaZ0,elemZ0);
		IGANextElement(igachiUp,elemchiUp);

		ierr = IGAElementGetWorkMat(elemFullS,&KlocFullS);CHKERRQ(ierr);
		ierr = IGAElementGetWorkVec(elemFullS,&FlocFullS);CHKERRQ(ierr);
		ierr = IGAElementGetValues(elemZ0,arrayZ0FullS,&Z0FullS);CHKERRQ(ierr);
		ierr = IGAElementGetValues(elemchiUp,arrayChi0FullS,&Chi0FullS);CHKERRQ(ierr);

		// FormSystem loop
		while (IGAElementNextFormSystem(elemFullS,&wtfFullS,&wtf2FullS)) 
		{
		// Quadrature loop
			ierr = IGAElementBeginPoint(elemFullS,&pointFullS);CHKERRQ(ierr);
			ierr = IGAElementBeginPoint(elemZ0,&pointZ0);CHKERRQ(ierr);
			ierr = IGAElementBeginPoint(elemchiUp,&pointchiUp);CHKERRQ(ierr);

			while (IGAElementNextPoint(elemFullS,pointFullS))
			{
				if(pointFullS->atboundary==1)
				{
					//IGAElementNextPoint(elemZ0,pointZ0);
					//IGAElementNextPoint(elemchiUp,pointchiUp);
				}

				if(pointFullS->atboundary==0 && pointZ0->atboundary==0 && pointchiUp->atboundary==0)
				{
					IGAElementNextPoint(elemZ0,pointZ0);
					IGAElementNextPoint(elemchiUp,pointchiUp);

					ierr = IGAPointGetWorkMat(pointFullS,&KpointFullS);CHKERRQ(ierr);
					ierr = IGAPointGetWorkVec(pointFullS,&FpointFullS);CHKERRQ(ierr);
					//	   FullS(IGAPoint p,IGAPoint pChi,IGAPoint pZu, PetscReal *K,PetscReal *F, PetscReal *Chi,PetscReal *Zu, void *ctx)
					ierr = FullS(pointFullS,pointchiUp,pointZ0,KpointFullS,FpointFullS,Chi0FullS,Z0FullS,NULL);CHKERRQ(ierr);
					ierr = IGAPointAddMat(pointFullS,KpointFullS,KlocFullS);CHKERRQ(ierr);
					ierr = IGAPointAddVec(pointFullS,FpointFullS,FlocFullS);CHKERRQ(ierr);
				}
			}
			if (pointZ0->index != -1)
			{
				IGAElementNextPoint(elemZ0,pointZ0);
			}
			while (pointchiUp->index != -1)
			{
				IGAElementNextPoint(elemchiUp,pointchiUp);
			}

			ierr = IGAElementEndPoint(elemFullS,&pointFullS);CHKERRQ(ierr);
			ierr = IGAElementEndPoint(elemZ0,&pointZ0);CHKERRQ(ierr);
			ierr = IGAElementEndPoint(elemchiUp,&pointchiUp);CHKERRQ(ierr);
		}

		//ierr = IGAElementFixSystem(elemFullS,KlocFullS,FlocFullS);CHKERRQ(ierr);					//This sets Dirichlet condition ¿? (Yes, this applies the conditions from IGASetBoundaryValue)
		ierr = IGAElementAssembleMat(elemFullS,KlocFullS,KFullS);CHKERRQ(ierr);
		ierr = IGAElementAssembleVec(elemFullS,FlocFullS,FFullS);CHKERRQ(ierr);
	}
	IGANextElement(igaZ0,elemZ0);
	IGANextElement(igachiUp,elemchiUp);

	ierr = IGAEndElement(igaFullS,&elemFullS);CHKERRQ(ierr);
	ierr = IGAEndElement(igaZ0,&elemZ0);CHKERRQ(ierr);
	ierr = IGAEndElement(igachiUp,&elemchiUp);CHKERRQ(ierr);

	// Restore local vectors u, Z0, Chi0 and arrays
	ierr = IGARestoreLocalVecArray(igaZ0,Z0,&localZ0FullS,&arrayZ0FullS);CHKERRQ(ierr);
	ierr = IGARestoreLocalVecArray(igachiUp,chiUp0,&localChi0FullS,&arrayChi0FullS);CHKERRQ(ierr);

	ierr = PetscPrintf(PETSC_COMM_WORLD,"Despues del ciclo, antes de armar la matriz y el vector \n");CHKERRQ(ierr);

	ierr = MatAssemblyBegin(KFullS,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyEnd  (KFullS,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

	ierr = VecAssemblyBegin(FFullS);CHKERRQ(ierr);
	ierr = VecAssemblyEnd  (FFullS);CHKERRQ(ierr);

	ierr = PetscPrintf(PETSC_COMM_WORLD,"Matriz y vector armados \n");CHKERRQ(ierr);

	ierr = KSPSetOperators(kspFullS,KFullS,KFullS);CHKERRQ(ierr);
	PC pcFullS;
	ierr = KSPGetPC(kspFullS,&pcFullS); CHKERRQ(ierr);
	ierr = PCSetType(pcFullS,PCSOR); CHKERRQ(ierr);
	ierr = KSPSetFromOptions(kspFullS);CHKERRQ(ierr);
	ierr = KSPSetType(kspFullS,KSPFCG);
	ierr = KSPSetTolerances(kspFullS,1.0e-16,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);CHKERRQ(ierr);
	ierr = PetscPrintf(PETSC_COMM_WORLD,"Antes del solve \n");CHKERRQ(ierr);
	ierr = KSPSolve(kspFullS,FFullS,Fs0);CHKERRQ(ierr);

	ierr = KSPDestroy(&kspFullS);CHKERRQ(ierr);
	ierr = MatDestroy(&KFullS);CHKERRQ(ierr);
	ierr = VecDestroy(&FFullS);CHKERRQ(ierr);

	char nameFullS[]="/FullSigma-2d-0.dat";
	char pathFullS[512];
	sprintf(pathFullS,"%s%s",direct,nameFullS);
	ierr = IGAWriteVec(igaFullS,Fs0,pathFullS);CHKERRQ(ierr);	

	ierr = VecDestroy(&Fs0);CHKERRQ(ierr);
	ierr = IGADestroy(&igaFullS);CHKERRQ(ierr);
//
*/

// Not adapted yet
//System for L2 projection of skew part of stress
	PetscPrintf(PETSC_COMM_WORLD,"\nSystem for Skew Stress starting \n\n");
	T=time(NULL);
	tm=*localtime(&T);
	PetscPrintf(PETSC_COMM_WORLD,"Current time is %02d:%02d:%02d \n",tm.tm_hour,tm.tm_min,tm.tm_sec);
	IGA igaSkewS;
	ierr = IGACreate(PETSC_COMM_WORLD,&igaSkewS);CHKERRQ(ierr);
	ierr = IGASetDim(igaSkewS,2);CHKERRQ(ierr);													//Spatial dimension of the problem
	ierr = IGASetDof(igaSkewS,4);CHKERRQ(ierr);													//Number of degrees of freedom, per node
	ierr = IGASetOrder(igaSkewS,2);CHKERRQ(ierr);													//Number of spatial derivatives to calculate
	ierr = IGASetFromOptions(igaSkewS);CHKERRQ(ierr);												//Note: The order (or degree) of the shape functions is given by the mesh!
	ierr = IGARead(igaSkewS,"./geometry3.dat");CHKERRQ(ierr);
	
	for (dir=0; dir<2; dir++)
	{
		ierr = IGASetRuleType(igaSkewS,dir,IGA_RULE_LEGENDRE);CHKERRQ(ierr);
		ierr = IGASetRuleSize(igaSkewS,dir,6);CHKERRQ(ierr);
	}
	ierr = IGASetUp(igaSkewS);CHKERRQ(ierr);

	for (dir=0; dir<2; dir++) 
	{
		for (side=0; side<2; side++) 
		{
			//ierr = IGASetBoundaryValue(iga,dir,side,dof,0.0);CHKERRQ(ierr);					// Dirichlet boundary conditions
			ierr = IGASetBoundaryForm(igaSkewS,dir,side,PETSC_FALSE);CHKERRQ(ierr);  				// Neumann boundary conditions
		}
	}

	Mat KSkewS;
	Vec Fss0,FSkewS;

	ierr = IGACreateMat(igaSkewS,&KSkewS);CHKERRQ(ierr);
	ierr = IGACreateVec(igaSkewS,&Fss0);CHKERRQ(ierr);
	ierr = IGACreateVec(igaSkewS,&FSkewS);CHKERRQ(ierr);

	IGAPoint		pointSkewS;													//point
	IGAElement		elemSkewS;													//element
	PetscReal		*KlocSkewS,*FlocSkewS;										//AA y BB
	PetscReal		*KpointSkewS,*FpointSkewS;									//KKK y FFF
	const PetscReal	*arrayZ0SkewS,*arrayChi0SkewS,*arrayZSSkewS,*arrayS0SkewS;								//arrayU
	Vec				localZ0SkewS,localChi0SkewS,localZSSkewS,localS0SkewS;								//localU
	PetscReal		*Chi0SkewS,*Z0SkewS,*ZSSkewS,*S0SkewS;										//U

  	IGAFormSystem	wtfSkewS;
 	void			*wtf2SkewS;

 	KSP kspSkewS;
	ierr = IGACreateKSP(igaSkewS,&kspSkewS);CHKERRQ(ierr);

	// Get local vectors Z0 and Chi0 and arrays
	ierr = IGAGetLocalVecArray(igachiUp,chiUp0,&localChi0SkewS,&arrayChi0SkewS);CHKERRQ(ierr);
	ierr = IGAGetLocalVecArray(igaZ0,Z0,&localZ0SkewS,&arrayZ0SkewS);CHKERRQ(ierr);
	ierr = IGAGetLocalVecArray(igaZS,ZS0,&localZSSkewS,&arrayZSSkewS);CHKERRQ(ierr);
	ierr = IGAGetLocalVecArray(igaS,s0,&localS0SkewS,&arrayS0SkewS);CHKERRQ(ierr);

	// Element loop
	ierr = IGABeginElement(igaSkewS,&elemSkewS);CHKERRQ(ierr);
	ierr = IGABeginElement(igaZ0,&elemZ0);CHKERRQ(ierr);
	ierr = IGABeginElement(igachiUp,&elemchiUp);CHKERRQ(ierr);
	ierr = IGABeginElement(igaZS,&elemZS);CHKERRQ(ierr);
	ierr = IGABeginElement(igaS,&elemS);CHKERRQ(ierr);

	ierr = PetscPrintf(PETSC_COMM_WORLD,"Antes del ciclo, todos los objetos creados \n");CHKERRQ(ierr);
	while (IGANextElement(igaSkewS,elemSkewS)) 
	{
		IGANextElement(igaZ0,elemZ0);
		IGANextElement(igachiUp,elemchiUp);
		IGANextElement(igaZS,elemZS);
		IGANextElement(igaS,elemS);

		ierr = IGAElementGetWorkMat(elemSkewS,&KlocSkewS);CHKERRQ(ierr);
		ierr = IGAElementGetWorkVec(elemSkewS,&FlocSkewS);CHKERRQ(ierr);
		ierr = IGAElementGetValues(elemZ0,arrayZ0SkewS,&Z0SkewS);CHKERRQ(ierr);
		ierr = IGAElementGetValues(elemchiUp,arrayChi0SkewS,&Chi0SkewS);CHKERRQ(ierr);
		ierr = IGAElementGetValues(elemZS,arrayZSSkewS,&ZSSkewS);CHKERRQ(ierr);
		ierr = IGAElementGetValues(elemS,arrayS0SkewS,&S0SkewS);CHKERRQ(ierr);

		// FormSystem loop
		while (IGAElementNextFormSystem(elemSkewS,&wtfSkewS,&wtf2SkewS)) 
		{
		// Quadrature loop
			ierr = IGAElementBeginPoint(elemSkewS,&pointSkewS);CHKERRQ(ierr);
			ierr = IGAElementBeginPoint(elemZ0,&pointZ0);CHKERRQ(ierr);
			ierr = IGAElementBeginPoint(elemchiUp,&pointchiUp);CHKERRQ(ierr);
			ierr = IGAElementBeginPoint(elemZS,&pointZS);CHKERRQ(ierr);
			ierr = IGAElementBeginPoint(elemS,&pointS);CHKERRQ(ierr);

			while (IGAElementNextPoint(elemSkewS,pointSkewS))
			{
				if(pointSkewS->atboundary==1)
				{
					//IGAElementNextPoint(elemZ0,pointZ0);
					//IGAElementNextPoint(elemchiUp,pointchiUp);
				}

				if(pointSkewS->atboundary==0 && pointZ0->atboundary==0 && pointchiUp->atboundary==0 && pointZS->atboundary==0 && pointS->atboundary==0)
				{
					IGAElementNextPoint(elemZ0,pointZ0);
					IGAElementNextPoint(elemchiUp,pointchiUp);

					ierr = IGAPointGetWorkMat(pointSkewS,&KpointSkewS);CHKERRQ(ierr);
					ierr = IGAPointGetWorkVec(pointSkewS,&FpointSkewS);CHKERRQ(ierr);
					//	   SkewS(IGAPoint p,IGAPoint pChi,IGAPoint pZu, PetscReal *K,PetscReal *F, PetscReal *Chi,PetscReal *Zu, void *ctx)
					ierr = SkewS(pointSkewS,pointchiUp,pointZ0,pointZS,pointS,KpointSkewS,FpointSkewS,Chi0SkewS,Z0SkewS,ZSSkewS,S0SkewS,NULL);CHKERRQ(ierr);
					ierr = IGAPointAddMat(pointSkewS,KpointSkewS,KlocSkewS);CHKERRQ(ierr);
					ierr = IGAPointAddVec(pointSkewS,FpointSkewS,FlocSkewS);CHKERRQ(ierr);
				}
			}
			if (pointZ0->index != -1)
			{
				IGAElementNextPoint(elemZ0,pointZ0);
			}
			while (pointchiUp->index != -1)
			{
				IGAElementNextPoint(elemchiUp,pointchiUp);
			}
			while (pointZS->index != -1)
			{
				IGAElementNextPoint(elemZS,pointZS);
			}
			while (pointS->index != -1)
			{
				IGAElementNextPoint(elemS,pointS);
			}

			ierr = IGAElementEndPoint(elemSkewS,&pointSkewS);CHKERRQ(ierr);
			ierr = IGAElementEndPoint(elemZ0,&pointZ0);CHKERRQ(ierr);
			ierr = IGAElementEndPoint(elemchiUp,&pointchiUp);CHKERRQ(ierr);
			ierr = IGAElementEndPoint(elemZS,&pointZS);CHKERRQ(ierr);
			ierr = IGAElementEndPoint(elemS,&pointS);CHKERRQ(ierr);
		}

		//ierr = IGAElementFixSystem(elemSkewS,KlocSkewS,FlocSkewS);CHKERRQ(ierr);					//This sets Dirichlet condition ¿? (Yes, this applies the conditions from IGASetBoundaryValue)
		ierr = IGAElementAssembleMat(elemSkewS,KlocSkewS,KSkewS);CHKERRQ(ierr);
		ierr = IGAElementAssembleVec(elemSkewS,FlocSkewS,FSkewS);CHKERRQ(ierr);
	}
	IGANextElement(igaZ0,elemZ0);
	IGANextElement(igachiUp,elemchiUp);
	IGANextElement(igaZS,elemZS);
	IGANextElement(igaS,elemS);

	ierr = IGAEndElement(igaSkewS,&elemSkewS);CHKERRQ(ierr);
	ierr = IGAEndElement(igaZ0,&elemZ0);CHKERRQ(ierr);
	ierr = IGAEndElement(igachiUp,&elemchiUp);CHKERRQ(ierr);
	ierr = IGAEndElement(igaZS,&elemZS);CHKERRQ(ierr);
	ierr = IGAEndElement(igaS,&elemS);CHKERRQ(ierr);

	// Restore local vectors u, Z0, Chi0 and arrays
	ierr = IGARestoreLocalVecArray(igaZ0,Z0,&localZ0SkewS,&arrayZ0SkewS);CHKERRQ(ierr);
	ierr = IGARestoreLocalVecArray(igachiUp,chiUp0,&localChi0SkewS,&arrayChi0SkewS);CHKERRQ(ierr);
	ierr = IGARestoreLocalVecArray(igaZS,ZS0,&localZSSkewS,&arrayZSSkewS);CHKERRQ(ierr);
	ierr = IGARestoreLocalVecArray(igaS,s0,&localS0SkewS,&arrayS0SkewS);CHKERRQ(ierr);

	ierr = PetscPrintf(PETSC_COMM_WORLD,"Despues del ciclo, antes de armar la matriz y el vector \n");CHKERRQ(ierr);

	ierr = MatAssemblyBegin(KSkewS,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyEnd  (KSkewS,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

	ierr = VecAssemblyBegin(FSkewS);CHKERRQ(ierr);
	ierr = VecAssemblyEnd  (FSkewS);CHKERRQ(ierr);

	ierr = PetscPrintf(PETSC_COMM_WORLD,"Matriz y vector armados \n");CHKERRQ(ierr);

	ierr = KSPSetOperators(kspSkewS,KSkewS,KSkewS);CHKERRQ(ierr);
	PC pcSkewS;
	ierr = KSPGetPC(kspSkewS,&pcSkewS); CHKERRQ(ierr);
	ierr = PCSetType(pcSkewS,PCSOR); CHKERRQ(ierr);
	ierr = KSPSetFromOptions(kspSkewS);CHKERRQ(ierr);
	ierr = KSPSetType(kspSkewS,KSPFCG);
	ierr = KSPSetTolerances(kspSkewS,1.0e-16,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);CHKERRQ(ierr);
	ierr = PetscPrintf(PETSC_COMM_WORLD,"Antes del solve \n");CHKERRQ(ierr);
	ierr = KSPSolve(kspSkewS,FSkewS,Fss0);CHKERRQ(ierr);

	ierr = KSPDestroy(&kspSkewS);CHKERRQ(ierr);
	ierr = MatDestroy(&KSkewS);CHKERRQ(ierr);
	ierr = VecDestroy(&FSkewS);CHKERRQ(ierr);

	char nameSkewS[]="/SkewSigma-2d-0.dat";
	char pathSkewS[512];
	sprintf(pathSkewS,"%s%s",direct,nameSkewS);
	ierr = IGAWriteVec(igaSkewS,Fss0,pathSkewS);CHKERRQ(ierr);	

	ierr = VecDestroy(&Fss0);CHKERRQ(ierr);
	ierr = IGADestroy(&igaSkewS);CHKERRQ(ierr);
//

/*
//System for L2 projection of V^{alpha}
	PetscPrintf(PETSC_COMM_WORLD,"\nSystem for V-alpha starting \n\n");
	T=time(NULL);
	tm=*localtime(&T);
	PetscPrintf(PETSC_COMM_WORLD,"Current time is %02d:%02d:%02d \n",tm.tm_hour,tm.tm_min,tm.tm_sec);
	IGA igaVa;
	ierr = IGACreate(PETSC_COMM_WORLD,&igaVa);CHKERRQ(ierr);
	ierr = IGASetDim(igaVa,2);CHKERRQ(ierr);													//Spatial dimension of the problem
	ierr = IGASetDof(igaVa,2);CHKERRQ(ierr);													//Number of degrees of freedom, per node
	ierr = IGASetOrder(igaVa,1);CHKERRQ(ierr);													//Number of spatial derivatives to calculate
	ierr = IGASetFromOptions(igaVa);CHKERRQ(ierr);												//Note: The order (or degree) of the shape functions is given by the mesh!
	ierr = IGARead(igaVa,"./geometry.dat");CHKERRQ(ierr);
	
	for (dir=0; dir<2; dir++)
	{
		ierr = IGASetRuleType(igaVa,dir,IGA_RULE_LEGENDRE);CHKERRQ(ierr);
		ierr = IGASetRuleSize(igaVa,dir,6);CHKERRQ(ierr);
	}
	ierr = IGASetUp(igaVa);CHKERRQ(ierr);

	for (dir=0; dir<2; dir++) 
	{
		for (side=0; side<2; side++) 
		{
			//ierr = IGASetBoundaryValue(iga,dir,side,dof,0.0);CHKERRQ(ierr);					// Dirichlet boundary conditions
			ierr = IGASetBoundaryForm(igaVa,dir,side,PETSC_TRUE);CHKERRQ(ierr);  				// Neumann boundary conditions
		}
	}

	Mat KVa;
	Vec Va0,FVa;

	ierr = IGACreateMat(igaVa,&KVa);CHKERRQ(ierr);
	ierr = IGACreateVec(igaVa,&Va0);CHKERRQ(ierr);
	ierr = IGACreateVec(igaVa,&FVa);CHKERRQ(ierr);

	IGAPoint		pointVa;															//point
	IGAElement		elemVa;																//element
	PetscReal		*KlocVa,*FlocVa;												//AA y BB
	PetscReal		*KpointVa,*FpointVa;											//KKK y FFF
	const PetscReal	*arrayZ0Va,*arrayChi0Va,*arrayAl0Va;			//arrayU
	Vec				localZ0Va,localChi0Va,localAl0Va;				//localU
	PetscReal		*Chi0Va,*Z0Va,*Al0Va;											//U

	IGAFormSystem	wtfVa;
	void			*wtf2Va;

 	KSP kspVa;
	ierr = IGACreateKSP(igaVa,&kspVa);CHKERRQ(ierr);

	// Get local vectors Z0 and Chi0 and arrays
	ierr = IGAGetLocalVecArray(igachiUp,chiUp0,&localChi0Va,&arrayChi0Va);CHKERRQ(ierr);
	ierr = IGAGetLocalVecArray(igaZ0,Z0,&localZ0Va,&arrayZ0Va);CHKERRQ(ierr);
	ierr = IGAGetLocalVecArray(igaAl,alInput,&localAl0Va,&arrayAl0Va);CHKERRQ(ierr);

	// Element loop
	ierr = IGABeginElement(igaVa,&elemVa);CHKERRQ(ierr);
	ierr = IGABeginElement(igaZ0,&elemZ0);CHKERRQ(ierr);
	ierr = IGABeginElement(igachiUp,&elemchiUp);CHKERRQ(ierr);
	ierr = IGABeginElement(igaAl,&elemAlp);CHKERRQ(ierr);

	while (IGANextElement(igaVa,elemVa)) 
	{
		IGANextElement(igaZ0,elemZ0);
		IGANextElement(igachiUp,elemchiUp);
		IGANextElement(igaAl,elemAlp);

		ierr = IGAElementGetWorkMat(elemVa,&KlocVa);CHKERRQ(ierr);
		ierr = IGAElementGetWorkVec(elemVa,&FlocVa);CHKERRQ(ierr);
		ierr = IGAElementGetValues(elemZ0,arrayZ0Va,&Z0Va);CHKERRQ(ierr);
		ierr = IGAElementGetValues(elemchiUp,arrayChi0Va,&Chi0Va);CHKERRQ(ierr);
		ierr = IGAElementGetValues(elemAlp,arrayAl0Va,&Al0Va);CHKERRQ(ierr);

		// FormSystem loop
		while (IGAElementNextFormSystem(elemVa,&wtfVa,&wtf2Va)) 
		{
		// Quadrature loop
			ierr = IGAElementBeginPoint(elemVa,&pointVa);CHKERRQ(ierr);
			ierr = IGAElementBeginPoint(elemZ0,&pointZ0);CHKERRQ(ierr);
			ierr = IGAElementBeginPoint(elemchiUp,&pointchiUp);CHKERRQ(ierr);
			ierr = IGAElementBeginPoint(elemAlp,&pointAlp);CHKERRQ(ierr);

			while (IGAElementNextPoint(elemVa,pointVa))
			{
				if(pointVa->atboundary==1)
				{
					//IGAElementNextPoint(elemZ0,pointZ0);
					//IGAElementNextPoint(elemchiUp,pointchiUp);
				}

				if(pointVa->atboundary==0 && pointZ0->atboundary==0 && pointchiUp->atboundary==0 && pointAlp->atboundary==0)
				{
					IGAElementNextPoint(elemZ0,pointZ0);
					IGAElementNextPoint(elemchiUp,pointchiUp);
					IGAElementNextPoint(elemAlp,pointAlp);

					ierr = IGAPointGetWorkMat(pointVa,&KpointVa);CHKERRQ(ierr);
					ierr = IGAPointGetWorkVec(pointVa,&FpointVa);CHKERRQ(ierr);
					//Valpha(IGAPoint p,IGAPoint pChi,IGAPoint pZu,IGAPoint pAl,PetscReal *K,PetscReal *F, PetscReal *Chi,PetscReal *Zu,PetscReal *Alu,void *ctx)
					ierr = Valpha(pointVa,pointchiUp,pointZ0,pointAlp,KpointVa,FpointVa,Chi0Va,Z0Va,Al0Va,NULL);CHKERRQ(ierr);
					ierr = IGAPointAddMat(pointVa,KpointVa,KlocVa);CHKERRQ(ierr);
					ierr = IGAPointAddVec(pointVa,FpointVa,FlocVa);CHKERRQ(ierr);
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
			while (pointAlp->index != -1)
			{
				IGAElementNextPoint(elemAlp,pointAlp);
			}
			ierr = IGAElementEndPoint(elemVa,&pointVa);CHKERRQ(ierr);
			ierr = IGAElementEndPoint(elemZ0,&pointZ0);CHKERRQ(ierr);
			ierr = IGAElementEndPoint(elemchiUp,&pointchiUp);CHKERRQ(ierr);
			ierr = IGAElementEndPoint(elemAlp,&pointAlp);CHKERRQ(ierr);
		}

		//ierr = IGAElementFixSystem(elemStress,KlocStress,FlocStress);CHKERRQ(ierr);					//This sets Dirichlet condition ¿? (Yes, this applies the conditions from IGASetBoundaryValue)
		ierr = IGAElementAssembleMat(elemVa,KlocVa,KVa);CHKERRQ(ierr);
		ierr = IGAElementAssembleVec(elemVa,FlocVa,FVa);CHKERRQ(ierr);
	}
	IGANextElement(igaZ0,elemZ0);
	IGANextElement(igachiUp,elemchiUp);
	IGANextElement(igaAl,elemAlp);

	ierr = IGAEndElement(igaVa,&elemVa);CHKERRQ(ierr);
	ierr = IGAEndElement(igaZ0,&elemZ0);CHKERRQ(ierr);
	ierr = IGAEndElement(igachiUp,&elemchiUp);CHKERRQ(ierr);
	ierr = IGAEndElement(igaAl,&elemAlp);CHKERRQ(ierr);

	// Restore local vectors u, Z0, Chi0 and arrays
	ierr = IGARestoreLocalVecArray(igaZ0,Z0,&localZ0Va,&arrayZ0Va);CHKERRQ(ierr);
	ierr = IGARestoreLocalVecArray(igachiUp,chiUp0,&localChi0Va,&arrayChi0Va);CHKERRQ(ierr);
	ierr = IGARestoreLocalVecArray(igaAl,alInput,&localAl0Va,&arrayAl0Va);CHKERRQ(ierr);
	ierr = MatAssemblyBegin(KVa,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyEnd  (KVa,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

	ierr = VecAssemblyBegin(FVa);CHKERRQ(ierr);
	ierr = VecAssemblyEnd  (FVa);CHKERRQ(ierr);

	ierr = KSPSetOperators(kspVa,KVa,KVa);CHKERRQ(ierr);
	PC pcVa;
	ierr = KSPGetPC(kspVa,&pcVa); CHKERRQ(ierr);
	ierr = PCSetType(pcVa,PCSOR); CHKERRQ(ierr);
	//ierr = PCFactorSetMatSolverType(pcStress,MATSOLVERMUMPS); CHKERRQ(ierr);
	ierr = KSPSetFromOptions(kspVa);CHKERRQ(ierr);
	ierr = KSPSetTolerances(kspVa,1.0e-12,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);CHKERRQ(ierr);
	ierr = KSPSolve(kspVa,FVa,Va0);CHKERRQ(ierr);

	//ierr = KSPDestroy(&kspVa);CHKERRQ(ierr);
	//ierr = MatDestroy(&KVa);CHKERRQ(ierr);
	ierr = VecDestroy(&FVa);CHKERRQ(ierr);

	char nameVa[]="/Va-2d-0.dat";
	char pathVa[512];
	sprintf(pathVa,"%s%s",direct,nameVa);
	ierr = IGAWriteVec(igaVa,Va0,pathVa);CHKERRQ(ierr);
//

//System for L2 proyection of smoothed V^{alpha}
	PetscPrintf(PETSC_COMM_WORLD,"\nSystem for Int(V-alpha) starting \n\n");
	T=time(NULL);
	tm=*localtime(&T);
	PetscPrintf(PETSC_COMM_WORLD,"Current time is %02d:%02d:%02d \n",tm.tm_hour,tm.tm_min,tm.tm_sec);

	//First integrate velocity for both defects
	Vec FVaInt, FVaXi, Int1aVec, Int2aVec, Int1bVec, Int2bVec;
	ierr = IGACreateVec(igaVa,&FVaInt);CHKERRQ(ierr);
	ierr = IGACreateVec(igaVa,&FVaXi);CHKERRQ(ierr);
	ierr = IGACreateVec(igaVa,&Int1aVec);CHKERRQ(ierr);
	ierr = IGACreateVec(igaVa,&Int2aVec);CHKERRQ(ierr);
	ierr = IGACreateVec(igaVa,&Int1bVec);CHKERRQ(ierr);
	ierr = IGACreateVec(igaVa,&Int2bVec);CHKERRQ(ierr);
																																		//Fix these comments
	PetscReal		*PointInt1a,*PointInt2a,*PointInt1b,*PointInt2b,*locInt1a,*locInt2a,*locInt1b,*locInt2b;	//AA y BB
	const PetscReal	*arrayVa0,*arrayAl1aVa,*arrayAl1bVa;																								//arrayU
	Vec				localVa0,localAl1aVa,localAl1bVa;																								//localU
	PetscReal 		*Al1aVa,*Al1bVa,*Va0Values,Int1a=0.0,Int2a=0.0,Int1b=0.0,Int2b=0.0;

	IGAFormSystem	wtfVaInt;
	void			*wtf2VaInt;

	// Get local vectors V0 and InputAl and arrays
	ierr = IGAGetLocalVecArray(igaVa,Va0,&localVa0,&arrayVa0);CHKERRQ(ierr);
	ierr = IGAGetLocalVecArray(igaAl,alInput1,&localAl1aVa,&arrayAl1aVa);CHKERRQ(ierr);
	ierr = IGAGetLocalVecArray(igaAl,alInput2,&localAl1bVa,&arrayAl1bVa);CHKERRQ(ierr);

	// Element loop
	ierr = IGABeginElement(igaVa,&elemVa);CHKERRQ(ierr);
	ierr = IGABeginElement(igaAl,&elemAlp);CHKERRQ(ierr);

	while (IGANextElement(igaVa,elemVa))
	{
		IGANextElement(igaAl,elemAlp);

		ierr = IGAElementGetWorkVec(elemVa,&locInt1a);CHKERRQ(ierr);
		ierr = IGAElementGetWorkVec(elemVa,&locInt2a);CHKERRQ(ierr);
		ierr = IGAElementGetWorkVec(elemVa,&locInt1b);CHKERRQ(ierr);
		ierr = IGAElementGetWorkVec(elemVa,&locInt2b);CHKERRQ(ierr);
		ierr = IGAElementGetValues(elemVa,arrayVa0,&Va0Values);CHKERRQ(ierr);
		ierr = IGAElementGetValues(elemAlp,arrayAl1aVa,&Al1aVa);CHKERRQ(ierr);
		ierr = IGAElementGetValues(elemAlp,arrayAl1bVa,&Al1bVa);CHKERRQ(ierr);

		// FormSystem loop
		while (IGAElementNextFormSystem(elemVa,&wtfVaInt,&wtf2VaInt))
		{
		// Quadrature loop
			ierr = IGAElementBeginPoint(elemVa,&pointVa);CHKERRQ(ierr);
			ierr = IGAElementBeginPoint(elemAlp,&pointAlp);CHKERRQ(ierr);
			
			while (IGAElementNextPoint(elemVa,pointVa))
			{
				if(pointVa->atboundary==1)
				{

				}

				if(pointVa->atboundary==0 && pointAlp->atboundary==0)
				{
					IGAElementNextPoint(elemAlp,pointAlp);
					ierr = IGAPointGetWorkVec(pointVa,&PointInt1a);CHKERRQ(ierr);
					ierr = IGAPointGetWorkVec(pointVa,&PointInt2a);CHKERRQ(ierr);
					ierr = IGAPointGetWorkVec(pointVa,&PointInt1b);CHKERRQ(ierr);
					ierr = IGAPointGetWorkVec(pointVa,&PointInt2b);CHKERRQ(ierr);
					//IntValpha(IGAPoint pV,IGAPoint pAl1,PetscReal *FInt1a,PetscReal *FInt2a,PetscReal *FInt1b,PetscReal *FInt2b,PetscReal *Valpha,PetscReal *Al1,PetscReal *Al2,void *ctx)
					ierr = IntValpha(pointVa,pointAlp,PointInt1a,PointInt2a,PointInt1b,PointInt2b,Va0Values,Al1aVa,Al1bVa,NULL);CHKERRQ(ierr);
					ierr = IGAPointAddVec(pointVa,PointInt1a,locInt1a);CHKERRQ(ierr);
					ierr = IGAPointAddVec(pointVa,PointInt2a,locInt2a);CHKERRQ(ierr);
					ierr = IGAPointAddVec(pointVa,PointInt1b,locInt1b);CHKERRQ(ierr);
					ierr = IGAPointAddVec(pointVa,PointInt2b,locInt2b);CHKERRQ(ierr);
				}
			}
			while (pointAlp->index != -1)
			{
				IGAElementNextPoint(elemAlp,pointAlp);
			}
			ierr = IGAElementEndPoint(elemVa,&pointVa);CHKERRQ(ierr);
			ierr = IGAElementEndPoint(elemAlp,&pointAlp);CHKERRQ(ierr);
		}
		ierr = IGAElementAssembleVec(elemVa,locInt1a,Int1aVec);CHKERRQ(ierr);
		ierr = IGAElementAssembleVec(elemVa,locInt2a,Int2aVec);CHKERRQ(ierr);
		ierr = IGAElementAssembleVec(elemVa,locInt1b,Int1bVec);CHKERRQ(ierr);
		ierr = IGAElementAssembleVec(elemVa,locInt2b,Int2bVec);CHKERRQ(ierr);
	}
	IGANextElement(igaAl,elemAlp);

	ierr = IGAEndElement(igaVa,&elemVa);CHKERRQ(ierr);
	ierr = IGAEndElement(igaAl,&elemAlp);CHKERRQ(ierr);

	ierr = IGARestoreLocalVecArray(igaAl,alInput1,&localAl1aVa,&arrayAl1aVa);CHKERRQ(ierr);
	ierr = IGARestoreLocalVecArray(igaAl,alInput2,&localAl1bVa,&arrayAl1bVa);CHKERRQ(ierr);
	ierr = IGARestoreLocalVecArray(igaVa,Va0,&localVa0,&arrayVa0);CHKERRQ(ierr);

	ierr = VecAssemblyBegin(Int1aVec);CHKERRQ(ierr);
	ierr = VecAssemblyEnd  (Int1aVec);CHKERRQ(ierr);
	ierr = VecAssemblyBegin(Int2aVec);CHKERRQ(ierr);
	ierr = VecAssemblyEnd  (Int2aVec);CHKERRQ(ierr);

	ierr = VecAssemblyBegin(Int1bVec);CHKERRQ(ierr);
	ierr = VecAssemblyEnd  (Int1bVec);CHKERRQ(ierr);
	ierr = VecAssemblyBegin(Int2bVec);CHKERRQ(ierr);
	ierr = VecAssemblyEnd  (Int2bVec);CHKERRQ(ierr);

	//Here add all values at Gauss points
	ierr = VecSum(Int1aVec,&Int1a);CHKERRQ(ierr);
	ierr = VecSum(Int2aVec,&Int2a);CHKERRQ(ierr);
	ierr = VecSum(Int1bVec,&Int1b);CHKERRQ(ierr);
	ierr = VecSum(Int2bVec,&Int2b);CHKERRQ(ierr);

	//Now put those velocities on each defect location
	Vec FVaSmooth1,FVaSmooth2,Va1Smooth,Va2Smooth;
	ierr = IGACreateVec(igaVa,&FVaSmooth1);CHKERRQ(ierr);
	ierr = IGACreateVec(igaVa,&FVaSmooth2);CHKERRQ(ierr);
	ierr = IGACreateVec(igaVa,&Va1Smooth);CHKERRQ(ierr);
	ierr = IGACreateVec(igaVa,&Va2Smooth);CHKERRQ(ierr);

	PetscReal		*VaSmooth1Point,*VaSmooth2Point,*locSmooth1,*locSmooth2;				//AA y BB
	//const PetscReal	*arrayVa0,*arrayAl1aVa,*arrayAl1bVa;									//arrayU 	//Probably this line will go, variables already declared
	//Vec				localVa0,localAl1aVa,localAl1bVa;										//localU 	//Probably this line will go, variables already declared
	//PetscReal 		*Al1aVa,*Al1bVa,*Va0Values;															//Probably this line will go, variables already declared

	//IGAFormSystem	wtfVaInt;																			//Probably this line will go, variables already declared
	//void			*wtf2VaInt;																			//Probably this line will go, variables already declared

	// Get local vectors V0 and InputAl and arrays
	ierr = IGAGetLocalVecArray(igaVa,Va0,&localVa0,&arrayVa0);CHKERRQ(ierr);
	ierr = IGAGetLocalVecArray(igaAl,alInput1,&localAl1aVa,&arrayAl1aVa);CHKERRQ(ierr);
	ierr = IGAGetLocalVecArray(igaAl,alInput2,&localAl1bVa,&arrayAl1bVa);CHKERRQ(ierr);

	// Element loop
	ierr = IGABeginElement(igaVa,&elemVa);CHKERRQ(ierr);
	ierr = IGABeginElement(igaAl,&elemAlp);CHKERRQ(ierr);

	while (IGANextElement(igaVa,elemVa))
	{
		IGANextElement(igaAl,elemAlp);

		ierr = IGAElementGetWorkVec(elemVa,&locSmooth1);CHKERRQ(ierr);
		ierr = IGAElementGetWorkVec(elemVa,&locSmooth2);CHKERRQ(ierr);
		ierr = IGAElementGetValues(elemVa,arrayVa0,&Va0Values);CHKERRQ(ierr);
		ierr = IGAElementGetValues(elemAlp,arrayAl1aVa,&Al1aVa);CHKERRQ(ierr);
		ierr = IGAElementGetValues(elemAlp,arrayAl1bVa,&Al1bVa);CHKERRQ(ierr);

		// FormSystem loop
		while (IGAElementNextFormSystem(elemVa,&wtfVaInt,&wtf2VaInt))
		{
		// Quadrature loop
			ierr = IGAElementBeginPoint(elemVa,&pointVa);CHKERRQ(ierr);
			ierr = IGAElementBeginPoint(elemAlp,&pointAlp);CHKERRQ(ierr);
			
			while (IGAElementNextPoint(elemVa,pointVa))
			{
				if(pointVa->atboundary==1)
				{

				}

				if(pointVa->atboundary==0 && pointAlp->atboundary==0)
				{
					IGAElementNextPoint(elemAlp,pointAlp);
					ierr = IGAPointGetWorkVec(pointVa,&VaSmooth1Point);CHKERRQ(ierr);
					ierr = IGAPointGetWorkVec(pointVa,&VaSmooth2Point);CHKERRQ(ierr);
					//	   ProjValpha(IGAPoint pV,IGAPoint pAl1,PetscReal *FVa1,PetscReal *FVa2,PetscReal Valpha1a,PetscReal Valpha2a,PetscReal Valpha1b,PetscReal Valpha2b,PetscReal *Al1,PetscReal *Al2,void *ctx)
					ierr = ProjValpha(pointVa,pointAlp,VaSmooth1Point,VaSmooth2Point,Al1aVa,Al1bVa,NULL);CHKERRQ(ierr);
					ierr = IGAPointAddVec(pointVa,VaSmooth1Point,locSmooth1);CHKERRQ(ierr);
					ierr = IGAPointAddVec(pointVa,VaSmooth2Point,locSmooth2);CHKERRQ(ierr);
				}
			}
			while (pointAlp->index != -1)
			{
				IGAElementNextPoint(elemAlp,pointAlp);
			}
			ierr = IGAElementEndPoint(elemVa,&pointVa);CHKERRQ(ierr);
			ierr = IGAElementEndPoint(elemAlp,&pointAlp);CHKERRQ(ierr);
		}
		ierr = IGAElementAssembleVec(elemVa,locSmooth1,FVaSmooth1);CHKERRQ(ierr);
		ierr = IGAElementAssembleVec(elemVa,locSmooth2,FVaSmooth2);CHKERRQ(ierr);
	}
	IGANextElement(igaAl,elemAlp);

	ierr = IGAEndElement(igaVa,&elemVa);CHKERRQ(ierr);
	ierr = IGAEndElement(igaAl,&elemAlp);CHKERRQ(ierr);

	ierr = IGARestoreLocalVecArray(igaAl,alInput1,&localAl1aVa,&arrayAl1aVa);CHKERRQ(ierr);
	ierr = IGARestoreLocalVecArray(igaAl,alInput2,&localAl1bVa,&arrayAl1bVa);CHKERRQ(ierr);
	ierr = IGARestoreLocalVecArray(igaVa,Va0,&localVa0,&arrayVa0);CHKERRQ(ierr);

	//Here we have a vector with an L2 Projection of the integrated velocity.
	ierr = VecAssemblyBegin(FVaSmooth1);CHKERRQ(ierr);
	ierr = VecAssemblyEnd  (FVaSmooth1);CHKERRQ(ierr);
	ierr = VecAssemblyBegin(FVaSmooth2);CHKERRQ(ierr);
	ierr = VecAssemblyEnd  (FVaSmooth2);CHKERRQ(ierr);

	//Here we solved the L2 projection system
	ierr = KSPSolve(kspVa,FVaSmooth1,Va1Smooth);CHKERRQ(ierr);
	ierr = KSPSolve(kspVa,FVaSmooth2,Va2Smooth);CHKERRQ(ierr);

	//Solutions are in Va1Smooth (for core #1) and Va2Smooth (for core #2)
	PetscInt indLow1,indHigh1,numInd1,indLow2,indHigh2,numInd2;
	ierr = VecGetOwnershipRange(Va1Smooth,&indLow1,&indHigh1);CHKERRQ(ierr);
	ierr = VecGetOwnershipRange(Va2Smooth,&indLow2,&indHigh2);CHKERRQ(ierr);
	numInd1=indHigh1-indLow1;
	numInd2=indHigh2-indLow2;

	PetscInt *indices1,*indices2,i,count;
	indices1=(PetscInt*)calloc(numInd1,sizeof(PetscInt));
	indices2=(PetscInt*)calloc(numInd2,sizeof(PetscInt));

	count=0;
	while(count<numInd1)
	{
		indices1[count]=indLow1+count;
		count++;
	}

	count=0;
	while(count<numInd2)
	{
		indices2[count]=indLow2+count;
		count++;
	}

	PetscReal *valores1,*valores2;
	valores1=(PetscReal*)calloc(numInd1,sizeof(PetscReal));
	valores2=(PetscReal*)calloc(numInd2,sizeof(PetscReal));

	ierr = VecGetValues(Va1Smooth,numInd1,indices1,valores1);
	ierr = VecGetValues(Va2Smooth,numInd2,indices2,valores2);

	for (i=0;i<numInd1;i=i+2)
	{
		valores1[i]=valores1[i]*valores1[i]*Int1a;
		valores1[i+1]=valores1[i+1]*valores1[i+1]*Int2a;
	}
	for (i=0;i<numInd2;i=i+2)
	{
		valores2[i]=valores2[i]*valores2[i]*Int1b;
		valores2[i+1]=valores2[i+1]*valores2[i+1]*Int2b;
	}

	ierr = VecSetValues(Va1Smooth,numInd1,indices1,valores1,INSERT_VALUES);CHKERRQ(ierr);
	ierr = VecSetValues(Va2Smooth,numInd2,indices2,valores2,INSERT_VALUES);CHKERRQ(ierr);

	ierr = VecAssemblyBegin(Va1Smooth);CHKERRQ(ierr);
	ierr = VecAssemblyEnd  (Va1Smooth);CHKERRQ(ierr);
	ierr = VecAssemblyBegin(Va2Smooth);CHKERRQ(ierr);
	ierr = VecAssemblyEnd  (Va2Smooth);CHKERRQ(ierr);

	char nameVaS1[]="/VaSmooth1-2d-0.dat";
	char pathVaS1[512];
	sprintf(pathVaS1,"%s%s",direct,nameVaS1);
	ierr = IGAWriteVec(igaVa,Va1Smooth,pathVaS1);CHKERRQ(ierr);

	char nameVaS2[]="/VaSmooth2-2d-0.dat";
	char pathVaS2[512];
	sprintf(pathVaS2,"%s%s",direct,nameVaS2);
	ierr = IGAWriteVec(igaVa,Va2Smooth,pathVaS2);CHKERRQ(ierr);

	ierr = KSPDestroy(&kspVa);CHKERRQ(ierr);
	ierr = MatDestroy(&KVa);CHKERRQ(ierr);
	ierr = VecDestroy(&FVaSmooth1);CHKERRQ(ierr);
	ierr = VecDestroy(&FVaSmooth2);CHKERRQ(ierr);
	ierr = VecDestroy(&Va1Smooth);CHKERRQ(ierr);
	ierr = VecDestroy(&Va2Smooth);CHKERRQ(ierr);
	ierr = VecDestroy(&Va0);CHKERRQ(ierr);
	ierr = IGADestroy(&igaVa);CHKERRQ(ierr);
//

//System for L2 projection of \hat{Ue}
	PetscPrintf(PETSC_COMM_WORLD,"\nSystem for Ue starting \n\n");
	T=time(NULL);
	tm=*localtime(&T);
	PetscPrintf(PETSC_COMM_WORLD,"Current time is %02d:%02d:%02d \n",tm.tm_hour,tm.tm_min,tm.tm_sec);
	IGA igaUe;
	ierr = IGACreate(PETSC_COMM_WORLD,&igaUe);CHKERRQ(ierr);
	ierr = IGASetDim(igaUe,2);CHKERRQ(ierr);													//Spatial dimension of the problem
	ierr = IGASetDof(igaUe,4);CHKERRQ(ierr);													//Number of degrees of freedom, per node
	ierr = IGASetOrder(igaUe,2);CHKERRQ(ierr);													//Number of spatial derivatives to calculate
	ierr = IGASetFromOptions(igaUe);CHKERRQ(ierr);												//Note: The order (or degree) of the shape functions is given by the mesh!
	ierr = IGARead(igaUe,"./geometry3.dat");CHKERRQ(ierr);
	
	for (dir=0; dir<2; dir++)
	{
		ierr = IGASetRuleType(igaUe,dir,IGA_RULE_LEGENDRE);CHKERRQ(ierr);
		ierr = IGASetRuleSize(igaUe,dir,6);CHKERRQ(ierr);
	}
	ierr = IGASetUp(igaUe);CHKERRQ(ierr);

	for (dir=0; dir<2; dir++) 
	{
		for (side=0; side<2; side++) 
		{
			ierr = IGASetBoundaryForm(igaUe,dir,side,PETSC_TRUE);CHKERRQ(ierr);  				// Neumann boundary conditions
		}
	}

	Mat KUe;
	Vec Ue,UeSkw,FUe,FUeSkw;

	ierr = IGACreateMat(igaUe,&KUe);CHKERRQ(ierr);
	ierr = IGACreateVec(igaUe,&Ue);CHKERRQ(ierr);
	ierr = IGACreateVec(igaUe,&UeSkw);CHKERRQ(ierr);
	ierr = IGACreateVec(igaUe,&FUe);CHKERRQ(ierr);
	ierr = IGACreateVec(igaUe,&FUeSkw);CHKERRQ(ierr);

	IGAPoint		pointUe;										//point
	IGAElement		elemUe;											//element
	PetscReal		*KlocUe,*FlocUe,*FlocUeSkw;						//AA y BB
	PetscReal		*KpointUe,*FpointUe,*FpointUeSkw;				//KKK y FFF
	const PetscReal	*arrayZ0Ue,*arrayChi0Ue;						//arrayU
	Vec				localZ0Ue,localChi0Ue;							//localU
	PetscReal		*Chi0Ue,*Z0Ue;									//U

  	IGAFormSystem	wtfUe;
 	void			*wtf2Ue;

 	KSP kspUe;
	ierr = IGACreateKSP(igaUe,&kspUe);CHKERRQ(ierr);

	// Get local vectors Z0 and Chi0 and arrays
	ierr = IGAGetLocalVecArray(igachiUp,chiUp0,&localChi0Ue,&arrayChi0Ue);CHKERRQ(ierr);
	ierr = IGAGetLocalVecArray(igaZ0,Z0,&localZ0Ue,&arrayZ0Ue);CHKERRQ(ierr);

	// Element loop
	ierr = IGABeginElement(igaUe,&elemUe);CHKERRQ(ierr);
	ierr = IGABeginElement(igaZ0,&elemZ0);CHKERRQ(ierr);
	ierr = IGABeginElement(igachiUp,&elemchiUp);CHKERRQ(ierr);

	while (IGANextElement(igaUe,elemUe)) 
	{
		IGANextElement(igaZ0,elemZ0);
		IGANextElement(igachiUp,elemchiUp);

		ierr = IGAElementGetWorkMat(elemUe,&KlocUe);CHKERRQ(ierr);
		ierr = IGAElementGetWorkVec(elemUe,&FlocUe);CHKERRQ(ierr);
		ierr = IGAElementGetWorkVec(elemUe,&FlocUeSkw);CHKERRQ(ierr);
		ierr = IGAElementGetValues(elemZ0,arrayZ0Ue,&Z0Ue);CHKERRQ(ierr);
		ierr = IGAElementGetValues(elemchiUp,arrayChi0Ue,&Chi0Ue);CHKERRQ(ierr);

		// FormSystem loop
		while (IGAElementNextFormSystem(elemUe,&wtfUe,&wtf2Ue)) 
		{
		// Quadrature loop
			ierr = IGAElementBeginPoint(elemUe,&pointUe);CHKERRQ(ierr);
			ierr = IGAElementBeginPoint(elemZ0,&pointZ0);CHKERRQ(ierr);
			ierr = IGAElementBeginPoint(elemchiUp,&pointchiUp);CHKERRQ(ierr);

			while (IGAElementNextPoint(elemUe,pointUe))
			{
				if(pointUe->atboundary==1)
				{
					//IGAElementNextPoint(elemZ0,pointZ0);
					//IGAElementNextPoint(elemchiUp,pointchiUp);
				}

				if(pointZ0->atboundary==1)
				{
					PetscPrintf(PETSC_COMM_WORLD,"Hola1 en Ue \n");
				}

				if(pointUe->atboundary==0 && pointZ0->atboundary==0 && pointchiUp->atboundary==0)
				{
					IGAElementNextPoint(elemZ0,pointZ0);
					IGAElementNextPoint(elemchiUp,pointchiUp);

					ierr = IGAPointGetWorkMat(pointUe,&KpointUe);CHKERRQ(ierr);
					ierr = IGAPointGetWorkVec(pointUe,&FpointUe);CHKERRQ(ierr);
					ierr = IGAPointGetWorkVec(pointUe,&FpointUeSkw);CHKERRQ(ierr);
					//PetscErrorCode Stress(IGAPoint p,IGAPoint pChi,IGAPoint pZu,PetscReal *K,PetscReal *F, PetscReal *Chi,PetscReal *Zu,PetscInt flag,void *ctx)
					ierr = skwUe(pointUe,pointchiUp,pointZ0,KpointUe,FpointUe,Chi0Ue,Z0Ue,0,NULL);CHKERRQ(ierr);
					ierr = skwUe(pointUe,pointchiUp,pointZ0,KpointUe,FpointUeSkw,Chi0Ue,Z0Ue,1,NULL);CHKERRQ(ierr);
					ierr = IGAPointAddMat(pointUe,KpointUe,KlocUe);CHKERRQ(ierr);
					ierr = IGAPointAddVec(pointUe,FpointUe,FlocUe);CHKERRQ(ierr);
					ierr = IGAPointAddVec(pointUe,FpointUeSkw,FlocUeSkw);CHKERRQ(ierr);
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
			ierr = IGAElementEndPoint(elemUe,&pointUe);CHKERRQ(ierr);
			ierr = IGAElementEndPoint(elemZ0,&pointZ0);CHKERRQ(ierr);
			ierr = IGAElementEndPoint(elemchiUp,&pointchiUp);CHKERRQ(ierr);
		}

		//ierr = IGAElementFixSystem(elemStress,KlocStress,FlocStress);CHKERRQ(ierr);					//This sets Dirichlet condition ¿? (Yes, this applies the conditions from IGASetBoundaryValue)
		ierr = IGAElementAssembleMat(elemUe,KlocUe,KUe);CHKERRQ(ierr);
		ierr = IGAElementAssembleVec(elemUe,FlocUe,FUe);CHKERRQ(ierr);
		ierr = IGAElementAssembleVec(elemUe,FlocUeSkw,FUeSkw);CHKERRQ(ierr);
	}
	IGANextElement(igaZ0,elemZ0);
	IGANextElement(igachiUp,elemchiUp);

	ierr = IGAEndElement(igaUe,&elemUe);CHKERRQ(ierr);
	ierr = IGAEndElement(igaZ0,&elemZ0);CHKERRQ(ierr);
	ierr = IGAEndElement(igachiUp,&elemchiUp);CHKERRQ(ierr);

	// Restore local vectors u, Z0, Chi0 and arrays
	ierr = IGARestoreLocalVecArray(igaZ0,Z0,&localZ0Ue,&arrayZ0Ue);CHKERRQ(ierr);
	ierr = IGARestoreLocalVecArray(igachiUp,chiUp0,&localChi0Ue,&arrayChi0Ue);CHKERRQ(ierr);
	ierr = MatAssemblyBegin(KUe,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyEnd  (KUe,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

	ierr = VecAssemblyBegin(FUe);CHKERRQ(ierr);
	ierr = VecAssemblyEnd  (FUe);CHKERRQ(ierr);

	ierr = VecAssemblyBegin(FUeSkw);CHKERRQ(ierr);
	ierr = VecAssemblyEnd  (FUeSkw);CHKERRQ(ierr);

	ierr = KSPSetOperators(kspUe,KUe,KUe);CHKERRQ(ierr);
	PC pcUe;
	ierr = KSPGetPC(kspUe,&pcUe); CHKERRQ(ierr);
	ierr = PCSetType(pcUe,PCSOR); CHKERRQ(ierr);
	//ierr = PCFactorSetMatSolverType(pcUe,MATSOLVERMUMPS); CHKERRQ(ierr);
	ierr = KSPSetFromOptions(kspUe);CHKERRQ(ierr);
	ierr = KSPSetTolerances(kspUe,1.0e-12,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);CHKERRQ(ierr);
	ierr = KSPSolve(kspUe,FUe,Ue);CHKERRQ(ierr);
	ierr = KSPSolve(kspUe,FUeSkw,UeSkw);CHKERRQ(ierr);

	ierr = KSPDestroy(&kspUe);CHKERRQ(ierr);
	ierr = MatDestroy(&KUe);CHKERRQ(ierr);
	ierr = VecDestroy(&FUe);CHKERRQ(ierr);
	ierr = VecDestroy(&FUeSkw);CHKERRQ(ierr);

	char nameUe[]="/Ue-2d-0.dat";
	char pathUe[512];
	sprintf(pathUe,"%s%s",direct,nameUe);
	ierr = IGAWriteVec(igaUe,Ue,pathUe);CHKERRQ(ierr);

	char nameUeSkw[]="/UeSkw-2d-0.dat";
	char pathUeSkw[512];
	sprintf(pathUeSkw,"%s%s",direct,nameUeSkw);
	ierr = IGAWriteVec(igaUe,UeSkw,pathUeSkw);CHKERRQ(ierr);

	ierr = VecDestroy(&UeSkw);CHKERRQ(ierr);
	ierr = IGADestroy(&igaUe);CHKERRQ(ierr);
//

//System for L2 projection of Norm \hat{UeSkw}
	//System for L2 projection of Ue
	PetscPrintf(PETSC_COMM_WORLD,"\n System for Elastic Energy Density starting \n\n");
	IGA igaNormUeSkw;
	ierr = IGACreate(PETSC_COMM_WORLD,&igaNormUeSkw);CHKERRQ(ierr);
	ierr = IGASetDim(igaNormUeSkw,2);CHKERRQ(ierr);													//Spatial dimension of the problem
	ierr = IGASetDof(igaNormUeSkw,1);CHKERRQ(ierr);													//Number of degrees of freedom, per node
	ierr = IGASetOrder(igaNormUeSkw,2);CHKERRQ(ierr);												//Number of spatial derivatives to calculate
	ierr = IGASetFromOptions(igaNormUeSkw);CHKERRQ(ierr);											//Note: The order (or degree) of the shape functions is given by the mesh!
	ierr = IGARead(igaNormUeSkw,"./geometry3.dat");CHKERRQ(ierr);
	ierr = IGASetUp(igaNormUeSkw);CHKERRQ(ierr);

	for (dir=0; dir<2; dir++) 
	{
		for (side=0; side<2; side++) 
		{
			//ierr = IGASetBoundaryValue(iga,dir,side,dof,0.0);CHKERRQ(ierr);					// Dirichlet boundary conditions
			ierr = IGASetBoundaryForm(igaNormUeSkw,dir,side,PETSC_TRUE);CHKERRQ(ierr);  				// Neumann boundary conditions
		}
	}

	Mat KNormUeSkw;
	Vec NormUeSkw,FNormUeSkw;

	ierr = IGACreateMat(igaNormUeSkw,&KNormUeSkw);CHKERRQ(ierr);
	ierr = IGACreateVec(igaNormUeSkw,&NormUeSkw);CHKERRQ(ierr);
	ierr = IGACreateVec(igaNormUeSkw,&FNormUeSkw);CHKERRQ(ierr);

	IGAPoint		pointNormUeSkw;															//point
	IGAElement		elemNormUeSkw;																//element
	PetscReal		*KlocNormUeSkw,*FlocNormUeSkw;													//AA y BB
	PetscReal		*KpointNormUeSkw,*FpointNormUeSkw;												//KKK y FFF
	const PetscReal	*arrayChi0NormUeSkw,*arrayZ0NormUeSkw;									//arrayU
	Vec				localChi0NormUeSkw,localZ0NormUeSkw;									//localU
	PetscReal		*chi0NormUeSkw,*Z0NormUeSkw;												//U

  	IGAFormSystem	wtfNormUeSkw;
 	void			*wtf2NormUeSkw;

 	KSP kspNormUeSkw;
	ierr = IGACreateKSP(igaNormUeSkw,&kspNormUeSkw);CHKERRQ(ierr);

	// Get local vectors Z0 and Chi0 and arrays
	ierr = IGAGetLocalVecArray(igachiUp,chiUp0,&localChi0NormUeSkw,&arrayChi0NormUeSkw);CHKERRQ(ierr);
	ierr = IGAGetLocalVecArray(igaZ0,Z0,&localZ0NormUeSkw,&arrayZ0NormUeSkw);CHKERRQ(ierr);

	// Element loop
	ierr = IGABeginElement(igaNormUeSkw,&elemNormUeSkw);CHKERRQ(ierr);
	ierr = IGABeginElement(igaZ0,&elemZ0);CHKERRQ(ierr);
	ierr = IGABeginElement(igachiUp,&elemchiUp);CHKERRQ(ierr);

	while (IGANextElement(igaNormUeSkw,elemNormUeSkw)) 				
	{
		IGANextElement(igachiUp,elemchiUp);
		IGANextElement(igaZ0,elemZ0);

		ierr = IGAElementGetWorkMat(elemNormUeSkw,&KlocNormUeSkw);CHKERRQ(ierr);
		ierr = IGAElementGetWorkVec(elemNormUeSkw,&FlocNormUeSkw);CHKERRQ(ierr);
		ierr = IGAElementGetValues(elemchiUp,arrayChi0NormUeSkw,&chi0NormUeSkw);CHKERRQ(ierr);
		ierr = IGAElementGetValues(elemZ0,arrayZ0NormUeSkw,&Z0NormUeSkw);CHKERRQ(ierr);

		// FormSystem loop
		while (IGAElementNextFormSystem(elemNormUeSkw,&wtfNormUeSkw,&wtf2NormUeSkw))
		{
		// Quadrature loop
			ierr = IGAElementBeginPoint(elemNormUeSkw,&pointNormUeSkw);CHKERRQ(ierr);
			ierr = IGAElementBeginPoint(elemZ0,&pointZ0);CHKERRQ(ierr);
			ierr = IGAElementBeginPoint(elemchiUp,&pointchiUp);CHKERRQ(ierr);

			while (IGAElementNextPoint(elemNormUeSkw,pointNormUeSkw))
			{
				if(pointNormUeSkw->atboundary==0 && pointchiUp->atboundary==0 && pointZ0->atboundary==0)
				{
					IGAElementNextPoint(elemZ0,pointZ0);
					IGAElementNextPoint(elemchiUp,pointchiUp);

					ierr = IGAPointGetWorkMat(pointNormUeSkw,&KpointNormUeSkw);CHKERRQ(ierr);
					ierr = IGAPointGetWorkVec(pointNormUeSkw,&FpointNormUeSkw);CHKERRQ(ierr);
						 //NormUeSkwSys(IGAPoint p,IGAPoint pChi,IGAPoint pZu,PetscReal *K,PetscReal *F, PetscReal *Chi,PetscReal *Zu,void *ctx)
					ierr = NormUeSkwSys(pointNormUeSkw,pointchiUp,pointZ0,KpointNormUeSkw,FpointNormUeSkw,chi0NormUeSkw,Z0NormUeSkw,NULL);CHKERRQ(ierr);
					ierr = IGAPointAddMat(pointNormUeSkw,KpointNormUeSkw,KlocNormUeSkw);CHKERRQ(ierr);
					ierr = IGAPointAddVec(pointNormUeSkw,FpointNormUeSkw,FlocNormUeSkw);CHKERRQ(ierr);
				}
			}
			while (pointchiUp->index != -1)
			{
				IGAElementNextPoint(elemchiUp,pointchiUp);
			}
			while (pointZ0->index != -1)
			{
				IGAElementNextPoint(elemZ0,pointZ0);
			}
			ierr = IGAElementEndPoint(elemNormUeSkw,&pointNormUeSkw);CHKERRQ(ierr);
			ierr = IGAElementEndPoint(elemchiUp,&pointchiUp);CHKERRQ(ierr);
			ierr = IGAElementEndPoint(elemZ0,&pointZ0);CHKERRQ(ierr);
		}

		ierr = IGAElementAssembleMat(elemNormUeSkw,KlocNormUeSkw,KNormUeSkw);CHKERRQ(ierr);
		ierr = IGAElementAssembleVec(elemNormUeSkw,FlocNormUeSkw,FNormUeSkw);CHKERRQ(ierr);
	}
	IGANextElement(igachiUp,elemchiUp);
	IGANextElement(igaZ0,elemZ0);

	ierr = IGAEndElement(igaNormUeSkw,&elemNormUeSkw);CHKERRQ(ierr);
	ierr = IGAEndElement(igaZ0,&elemZ0);CHKERRQ(ierr);
	ierr = IGAEndElement(igachiUp,&elemchiUp);CHKERRQ(ierr);

	// Restore local vectors u, Z0, Chi0 and arrays
	ierr = IGARestoreLocalVecArray(igachiUp,chiUp0,&localChi0NormUeSkw,&arrayChi0NormUeSkw);CHKERRQ(ierr);
	ierr = IGARestoreLocalVecArray(igaZ0,Z0,&localZ0NormUeSkw,&arrayZ0NormUeSkw);CHKERRQ(ierr);

	ierr = MatAssemblyBegin(KNormUeSkw,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyEnd  (KNormUeSkw,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

	ierr = VecAssemblyBegin(FNormUeSkw);CHKERRQ(ierr);
	ierr = VecAssemblyEnd  (FNormUeSkw);CHKERRQ(ierr);

	ierr = KSPSetOperators(kspNormUeSkw,KNormUeSkw,KNormUeSkw);CHKERRQ(ierr);
	PC pcNormUeSkw;
	ierr = KSPGetPC(kspNormUeSkw,&pcNormUeSkw); CHKERRQ(ierr);
	ierr = PCSetType(pcNormUeSkw,PCSOR); CHKERRQ(ierr);
	//ierr = PCFactorSetMatSolverType(pcNormUeSkw,MATSOLVERMUMPS); CHKERRQ(ierr);
	ierr = KSPSetFromOptions(kspNormUeSkw);CHKERRQ(ierr);
	ierr = KSPSetTolerances(kspNormUeSkw,1.0e-12,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);CHKERRQ(ierr);
	ierr = KSPSolve(kspNormUeSkw,FNormUeSkw,NormUeSkw);CHKERRQ(ierr);

	ierr = KSPDestroy(&kspNormUeSkw);CHKERRQ(ierr);
	ierr = MatDestroy(&KNormUeSkw);CHKERRQ(ierr);
	ierr = VecDestroy(&FNormUeSkw);CHKERRQ(ierr);

	char nameNormUeSkw[]="/NormUeSkw.dat";
	char pathNormUeSkw[512];
	sprintf(pathNormUeSkw,"%s%s",direct,nameNormUeSkw);
	ierr = IGAWriteVec(igaNormUeSkw,NormUeSkw,pathNormUeSkw);CHKERRQ(ierr);

	ierr = VecDestroy(&NormUeSkw);CHKERRQ(ierr);
	ierr = IGADestroy(&igaNormUeSkw);CHKERRQ(ierr);
//
*/

/*
//System for L2 projection of exact stress
	PetscPrintf(PETSC_COMM_WORLD,"\nSystem for L2 projection for exact stress starting \n\n");
	T=time(NULL);
	tm=*localtime(&T);
	PetscPrintf(PETSC_COMM_WORLD,"Current time is %02d:%02d:%02d \n",tm.tm_hour,tm.tm_min,tm.tm_sec);
	ierr = PetscPrintf(PETSC_COMM_WORLD,"Elastic exact system, step 0 \n");CHKERRQ(ierr);
	IGA igaExact;
	ierr = IGACreate(PETSC_COMM_WORLD,&igaExact);CHKERRQ(ierr);
	ierr = IGASetDim(igaExact,2);CHKERRQ(ierr);														//Spatial dimension of the problem
	ierr = IGASetDof(igaExact,4);CHKERRQ(ierr);														//Number of degrees of freedom, per node
	ierr = IGASetOrder(igaExact,1);CHKERRQ(ierr);													//Number of spatial derivatives to calculate
	ierr = PetscPrintf(PETSC_COMM_WORLD,"Elastic exact system, step 0.5 \n");CHKERRQ(ierr);
	ierr = IGASetFromOptions(igaExact);CHKERRQ(ierr);												//Note: The order (or degree) of the shape functions is given by the mesh!
	ierr = IGARead(igaExact,"./geometry3.dat");CHKERRQ(ierr);
	
	ierr = PetscPrintf(PETSC_COMM_WORLD,"Elastic exact system, step 1 \n");CHKERRQ(ierr);

	for (dir=0; dir<2; dir++)
	{
		ierr = IGASetRuleType(igaExact,dir,IGA_RULE_LEGENDRE);CHKERRQ(ierr);
		ierr = IGASetRuleSize(igaExact,dir,4);CHKERRQ(ierr);
	}
	ierr = PetscPrintf(PETSC_COMM_WORLD,"Elastic exact system, step 2 \n");CHKERRQ(ierr);

	ierr = IGASetUp(igaExact);CHKERRQ(ierr);

	ierr = PetscPrintf(PETSC_COMM_WORLD,"Elastic exact system, step 3 \n");CHKERRQ(ierr);
	Vec e0;
	Mat Kl2e;
	Vec Fl2e;
	ierr = IGACreateVec(igaExact,&e0);CHKERRQ(ierr);  
	ierr = IGACreateMat(igaExact,&Kl2e);CHKERRQ(ierr);
	ierr = IGACreateVec(igaExact,&Fl2e);CHKERRQ(ierr);
	ierr = IGASetFormSystem(igaExact,L2ProjectionExactStress,&userL2);CHKERRQ(ierr);
	ierr = IGAComputeSystem(igaExact,Kl2e,Fl2e);CHKERRQ(ierr);

	ierr = PetscPrintf(PETSC_COMM_WORLD,"Elastic exact system, step 4 \n");CHKERRQ(ierr);

	//This parts set and calls KSP to solve the linear system
	KSP kspl2e;
	ierr = IGACreateKSP(igaExact,&kspl2e);CHKERRQ(ierr);										
	ierr = KSPSetOperators(kspl2e,Kl2e,Kl2e);CHKERRQ(ierr); 								//This function creates the matrix for the system on the second parameter and uses the 3rd parameter as a preconditioner
	PC pcl2e;
	ierr = KSPGetPC(kspl2e,&pcl2e); CHKERRQ(ierr);
	ierr = PCSetType(pcl2e,PCSOR); CHKERRQ(ierr);
	ierr = KSPSetFromOptions(kspl2e);CHKERRQ(ierr);
	ierr = KSPSetType(kspl2e,KSPFCG);CHKERRQ(ierr);
	ierr = KSPSetTolerances(kspl2e,1.0e-16,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);CHKERRQ(ierr);
	ierr = KSPSolve(kspl2e,Fl2e,e0);CHKERRQ(ierr);										//This is a simple system, so it can be solved with just this command

	ierr = KSPDestroy(&kspl2e);CHKERRQ(ierr);
	ierr = MatDestroy(&Kl2e);CHKERRQ(ierr);
	ierr = VecDestroy(&Fl2e);CHKERRQ(ierr);
	char nameE[]="/ExactStress-2d-0.dat";
	char pathE[512];
	sprintf(pathE,"%s%s",direct,nameE);
	ierr = IGAWriteVec(igaExact,e0,pathE);CHKERRQ(ierr);
//
*/

/*
//System for L2 projection of Peierls stress
	PetscPrintf(PETSC_COMM_WORLD,"\nSystem for L2 projection for Peierls stress starting \n\n");
	T=time(NULL);
	tm=*localtime(&T);
	PetscPrintf(PETSC_COMM_WORLD,"Current time is %02d:%02d:%02d \n",tm.tm_hour,tm.tm_min,tm.tm_sec);
	IGA igaPeierls;
	ierr = IGACreate(PETSC_COMM_WORLD,&igaPeierls);CHKERRQ(ierr);
	ierr = IGASetDim(igaPeierls,2);CHKERRQ(ierr);														//Spatial dimension of the problem
	ierr = IGASetDof(igaPeierls,4);CHKERRQ(ierr);														//Number of degrees of freedom, per node
	ierr = IGASetOrder(igaPeierls,2);CHKERRQ(ierr);													//Number of spatial derivatives to calculate
	ierr = IGASetFromOptions(igaPeierls);CHKERRQ(ierr);												//Note: The order (or degree) of the shape functions is given by the mesh!
	ierr = IGARead(igaPeierls,"./geometry3.dat");CHKERRQ(ierr);
	
	for (dir=0; dir<2; dir++)
	{
		ierr = IGASetRuleType(igaPeierls,dir,IGA_RULE_LEGENDRE);CHKERRQ(ierr);
		ierr = IGASetRuleSize(igaPeierls,dir,6);CHKERRQ(ierr);
	}
	ierr = IGASetUp(igaPeierls);CHKERRQ(ierr);
	
	for (dir=0; dir<2; dir++) 
	{
		for (side=0; side<2; side++) 
		{
			//ierr = IGASetBoundaryValue(iga,dir,side,dof,0.0);CHKERRQ(ierr);    				// Dirichlet boundary conditions
			//ierr = IGASetBoundaryForm(igaPeierls,dir,side,PETSC_TRUE);CHKERRQ(ierr);  				// Neumann boundary conditions
		}
	}

	Vec stress_peierls;
	Mat Kl2p;
	Vec Fl2p;
	ierr = IGACreateVec(igaPeierls,&stress_peierls);CHKERRQ(ierr);  
	ierr = IGACreateMat(igaPeierls,&Kl2p);CHKERRQ(ierr);
	ierr = IGACreateVec(igaPeierls,&Fl2p);CHKERRQ(ierr);
	ierr = IGASetFormSystem(igaPeierls,L2ProjectionPeierlsStress,&userL2);CHKERRQ(ierr);
	ierr = IGAComputeSystem(igaPeierls,Kl2p,Fl2p);CHKERRQ(ierr);

	//This parts set and calls KSP to solve the linear system
	KSP kspl2p;
	ierr = IGACreateKSP(igaPeierls,&kspl2p);CHKERRQ(ierr);										
	ierr = KSPSetOperators(kspl2p,Kl2p,Kl2p);CHKERRQ(ierr); 								//This function creates the matrix for the system on the second parameter and uses the 3rd parameter as a preconditioner
	ierr = KSPSetType(kspl2p,KSPCG);CHKERRQ(ierr);											//Using KSPCG (conjugated gradient) because the matrix is symmetric
	//ierr = KSPSetOptionsPrefix(kspl2S,"l2pS_");CHKERRQ(ierr);
	ierr = KSPSetFromOptions(kspl2p);CHKERRQ(ierr);
	ierr = KSPSetTolerances(kspl2p,1.0e-12,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);CHKERRQ(ierr);
	ierr = KSPSolve(kspl2p,Fl2p,stress_peierls);CHKERRQ(ierr);										//This is a simple system, so it can be solved with just this command

	ierr = KSPDestroy(&kspl2p);CHKERRQ(ierr);
	ierr = MatDestroy(&Kl2p);CHKERRQ(ierr);
	ierr = VecDestroy(&Fl2p);CHKERRQ(ierr);
	char nameP[]="/PeierlsStress-2d-0.dat";
	char pathP[512];
	sprintf(pathP,"%s%s",direct,nameP);
	ierr = IGAWriteVec(igaPeierls,stress_peierls,pathP);CHKERRQ(ierr);
//
*/

//Destroy all objects not needed anymore (Better to do it here in case different codes call the same IGA, move if memory is a problem)
	T=time(NULL);
	tm=*localtime(&T);
	PetscPrintf(PETSC_COMM_WORLD,"Current time is %02d:%02d:%02d \n",tm.tm_hour,tm.tm_min,tm.tm_sec);
	//ierr = IGADestroy(&igaAl);CHKERRQ(ierr);
	//ierr = IGADestroy(&igachiS);CHKERRQ(ierr);
	//ierr = IGADestroy(&igachiUp);CHKERRQ(ierr);
	//ierr = IGADestroy(&igaExact);CHKERRQ(ierr);
//

ierr = PetscFinalize();CHKERRQ(ierr);

return 0;
}
