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
	PetscReal dt;
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

//Helmholtz decomposition of S, curl part for S based input
	#undef  __FUNCT__
	#define __FUNCT__ "curlChiSM"
	PetscErrorCode curlChiSM(IGAPoint p,PetscReal *K,void *ctx)
	{
		const PetscReal (*N1)[2];
		IGAPointGetShapeFuns(p,1,(const PetscReal**)&N1);

		PetscInt a,b,u,w,i,j,k,m,nen=p->nen, dof=p->dof;

		PetscReal (*KchiS)[dof][nen][dof] = (typeof(KchiS)) K;

		if (p->atboundary)
		{
			return 0;
		}
		else
		{
			for (a=0; a<nen; a++)
			{
				PetscReal Na_x = N1[a][0];
				PetscReal Na_y = N1[a][1];
				for (u=0; u<dof; u++) 
				{
					PetscReal dv[3][3][3][3]={0};

					if(u==0){
						dv[0][0][0][0]=Na_x; dv[0][0][0][1]=Na_y;}
					if(u==1){
						dv[0][0][1][0]=Na_x; dv[0][0][1][1]=Na_y;}
					if(u==2){
						dv[0][1][0][0]=Na_x; dv[0][1][0][1]=Na_y;}
					if(u==3){
						dv[0][1][1][0]=Na_x; dv[0][1][1][1]=Na_y;}
					if(u==4){
						dv[1][0][0][0]=Na_x; dv[1][0][0][1]=Na_y;}
					if(u==5){
						dv[1][0][1][0]=Na_x; dv[1][0][1][1]=Na_y;}
					if(u==6){
						dv[1][1][0][0]=Na_x; dv[1][1][0][1]=Na_y;}
					if(u==7){
						dv[1][1][1][0]=Na_x; dv[1][1][1][1]=Na_y;}

					for (b=0; b<nen; b++)
					{
						PetscReal Nb_x = N1[b][0];
						PetscReal Nb_y = N1[b][1];
						for (w=0; w<dof; w++)
						{
							PetscReal dchiS[3][3][3][3]={0};

							if(w==0){
								dchiS[0][0][0][0]=Nb_x; dchiS[0][0][0][1]=Nb_y;}
							if(w==1){
								dchiS[0][0][1][0]=Nb_x; dchiS[0][0][1][1]=Nb_y;}
							if(w==2){
								dchiS[0][1][0][0]=Nb_x; dchiS[0][1][0][1]=Nb_y;}
							if(w==3){
								dchiS[0][1][1][0]=Nb_x; dchiS[0][1][1][1]=Nb_y;}
							if(w==4){
								dchiS[1][0][0][0]=Nb_x; dchiS[1][0][0][1]=Nb_y;}
							if(w==5){
								dchiS[1][0][1][0]=Nb_x; dchiS[1][0][1][1]=Nb_y;}
							if(w==6){
								dchiS[1][1][0][0]=Nb_x; dchiS[1][1][0][1]=Nb_y;}
							if(w==7){
								dchiS[1][1][1][0]=Nb_x; dchiS[1][1][1][1]=Nb_y;}

							KchiS[a][u][b][w]=0.0;
							for(i=0;i<3;i++)
							{
								for(j=0;j<3;j++)
								{
									for(k=0;k<3;k++)
									{
										for (m=0;m<3;m++)
										{
											//First two terms are curl(chi):curl(v), it's an order of magnitude faster to expand the product of e_{kmn}e_{kab}
											KchiS[a][u][b][w]+=dchiS[i][j][k][m]*(dv[i][j][k][m]-dv[i][j][m][k])+dchiS[i][j][k][k]*dv[i][j][m][m];
										}
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

	#undef  __FUNCT__
	#define __FUNCT__ "curlChiSF"
	PetscErrorCode curlChiSF(IGAPoint p,IGAPoint pS,PetscReal *F,PetscReal *US,void *ctx)
	{
		const PetscReal (*N1)[2];
		IGAPointGetShapeFuns(p,1,(const PetscReal**)&N1);

		PetscInt a,u,i,j,k,m,nen=p->nen, dof=p->dof;

		PetscReal dS[8][2];																		//Create array to receive dS
		IGAPointFormGrad (pS,US,&dS[0][0]);														//This fills the values of the derivatives
		
		PetscReal fulld_S[3][3][3][3]={0};
		//Expanding grad(S) to full 4th orden tensor form, easier to work in for loops
		fulld_S[0][0][0][0]=dS[0][0]; fulld_S[0][0][0][1]=dS[0][1]; 
		fulld_S[0][0][1][0]=dS[1][0]; fulld_S[0][0][1][1]=dS[1][1]; 
		fulld_S[0][1][0][0]=dS[2][0]; fulld_S[0][1][0][1]=dS[2][1]; 
		fulld_S[0][1][1][0]=dS[3][0]; fulld_S[0][1][1][1]=dS[3][1];
		fulld_S[1][0][0][0]=dS[4][0]; fulld_S[1][0][0][1]=dS[4][1]; 
		fulld_S[1][0][1][0]=dS[5][0]; fulld_S[1][0][1][1]=dS[5][1]; 
		fulld_S[1][1][0][0]=dS[6][0]; fulld_S[1][1][0][1]=dS[6][1]; 
		fulld_S[1][1][1][0]=dS[7][0]; fulld_S[1][1][1][1]=dS[7][1];

		PetscReal (*FSp)[dof] = (PetscReal (*)[dof]) F;

		if (p->atboundary)
		{
			return 0;
		}
		else
		{
			for (a=0; a<nen; a++)
			{
				PetscReal Na_x = N1[a][0];
				PetscReal Na_y = N1[a][1];
				for (u=0; u<dof; u++) 
				{
					PetscReal dv[3][3][3][3]={0};
					if(u==0){
						dv[0][0][0][0]=Na_x; dv[0][0][0][1]=Na_y;}
					if(u==1){
						dv[0][0][1][0]=Na_x; dv[0][0][1][1]=Na_y;}
					if(u==2){
						dv[0][1][0][0]=Na_x; dv[0][1][0][1]=Na_y;}
					if(u==3){
						dv[0][1][1][0]=Na_x; dv[0][1][1][1]=Na_y;}
					if(u==4){
						dv[1][0][0][0]=Na_x; dv[1][0][0][1]=Na_y;}
					if(u==5){
						dv[1][0][1][0]=Na_x; dv[1][0][1][1]=Na_y;}
					if(u==6){
						dv[1][1][0][0]=Na_x; dv[1][1][0][1]=Na_y;}
					if(u==7){
						dv[1][1][1][0]=Na_x; dv[1][1][1][1]=Na_y;}

					FSp[a][u]=0.0;
					for(i=0;i<3;i++)
					{
						for(j=0;j<3;j++)
						{
							for(k=0;k<3;k++)
							{
								for(m=0;m<3;m++)
								{
									FSp[a][u]+=fulld_S[i][j][k][m]*(dv[i][j][k][m]-dv[i][j][m][k]);
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

//Helmholtz decomposition of S, grad part
	#undef  __FUNCT__
	#define __FUNCT__ "gradZSM"
	PetscErrorCode gradZSM(IGAPoint p,PetscReal *K,void *ctx)
	{

		const PetscReal (*N1)[2];
		IGAPointGetShapeFuns(p,1,(const PetscReal**)&N1);

		PetscInt a,b,u,w,i,j,k,nen=p->nen, dof=p->dof;

		PetscReal (*K_ZS)[dof][nen][dof] = (typeof(K_ZS)) K;

		if (p->atboundary)
		{
			return 0;
		}
		else
		{
			for (a=0; a<nen; a++)
			{
				PetscReal Na_x = N1[a][0];
				PetscReal Na_y = N1[a][1];
				for (u=0; u<dof; u++) 
				{
					PetscReal dv[3][3][3]={0};

					if(u==0)
					{
						dv[0][0][0]=Na_x; dv[0][0][1]=Na_y;
					}
					if(u==1)
					{
						dv[0][1][0]=Na_x; dv[0][1][1]=Na_y;
					}
					if(u==2)
					{
						dv[1][0][0]=Na_x; dv[1][0][1]=Na_y;
					}
					if(u==3)
					{
						dv[1][1][0]=Na_x; dv[1][1][1]=Na_y;
					}

					for (b=0; b<nen; b++)
					{
						PetscReal Nb_x = N1[b][0];
						PetscReal Nb_y = N1[b][1];
						for (w=0; w<dof; w++)
						{
							PetscReal dZS[3][3][3]={0};

							if(w==0)
							{
								dZS[0][0][0]=Nb_x; dZS[0][0][1]=Nb_y;
							}
							if(w==1)
							{
								dZS[0][1][0]=Nb_x; dZS[0][1][1]=Nb_y;
							}
							if(w==2)
							{
								dZS[1][0][0]=Nb_x; dZS[1][0][1]=Nb_y;
							}
							if(w==3)
							{
								dZS[1][1][0]=Nb_x; dZS[1][1][1]=Nb_y;
							}

							K_ZS[a][u][b][w]=0.0;
							for(i=0;i<3;i++)
							{
								for(j=0;j<3;j++)
								{
									for(k=0;k<3;k++)
									{
										K_ZS[a][u][b][w]+=dZS[i][j][k]*dv[i][j][k];
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

	#undef  __FUNCT__
	#define __FUNCT__ "gradZSF"
	PetscErrorCode gradZSF(IGAPoint p,IGAPoint pS,PetscReal *F,PetscReal *US,void *ctx)
	{

		const PetscReal (*N1)[2];
		IGAPointGetShapeFuns(p,1,(const PetscReal**)&N1);

		PetscInt a,u,i,j,k,nen=p->nen, dof=p->dof;

		PetscReal S[8];																		//Create array to receive Alfa
		IGAPointFormValue(pS,US,&S[0]);															//This fills the values

		PetscReal fullS[3][3][3]={0};
		//Expand S to full 3rd order form, only non-zero elements
		fullS[0][0][0]=S[0]; fullS[0][0][1]=S[1];
		fullS[0][1][0]=S[2]; fullS[0][1][1]=S[3];		
		fullS[1][0][0]=S[4]; fullS[1][0][1]=S[5]; 
		fullS[1][1][0]=S[6]; fullS[1][1][1]=S[7];

		PetscReal (*F_ZS)[dof] = (PetscReal (*)[dof]) F;

		if (p->atboundary)
		{
			return 0;
		}
		else
		{
			for (a=0; a<nen; a++)
			{
				PetscReal Na_x = N1[a][0];
				PetscReal Na_y = N1[a][1];
				for (u=0; u<dof; u++) 
				{
					PetscReal dv[3][3][3]={0};

					if(u==0)
					{
						dv[0][0][0]=Na_x; dv[0][0][1]=Na_y;
					}
					if(u==1)
					{
						dv[0][1][0]=Na_x; dv[0][1][1]=Na_y;
					}
					if(u==2)
					{
						dv[1][0][0]=Na_x; dv[1][0][1]=Na_y;
					}
					if(u==3)
					{
						dv[1][1][0]=Na_x; dv[1][1][1]=Na_y;
					}

					for(i=0;i<3;i++)
					{
						for(j=0;j<3;j++)
						{
							for(k=0;k<3;k++)
							{
								F_ZS[a][u]+=fullS[i][j][k]*dv[i][j][k];
							}
						}
					}
				}
			}
		}

		return 0;
	}
//

//L2 projection of alpha_tilde and alpha_hat
	#undef  __FUNCT__
	#define __FUNCT__ "L2Projection2DOF"
	PetscErrorCode L2Projection2DOF(IGAPoint p,PetscReal *K,void *ctx)
	{
		//This function creates an L2 projection matrix for a 2 DOF system
		if (p->atboundary)
		{
			return 0;																		//No BCs required for L2 Projection
		}

		PetscInt a,b,i;
		PetscInt nen = p->nen;																//Number of shape functions per element
		PetscInt dof = p->dof;																//Number of degrees of freedom per node

		const PetscReal (*N) = (typeof(N)) p->shape[0];
		PetscReal (*KK)[dof][nen][dof] = (PetscReal (*)[dof][nen][dof])K;
		for(a=0; a<nen; a++)
		{
			for(i=0; i<dof; i++) 
			{
				for(b=0; b<nen; b++)
				{
					KK[a][i][b][i] = N[a]*N[b];
				}
			}
		}
		return 0;
	}

	#undef  __FUNCT__
	#define __FUNCT__ "S_X"
	PetscErrorCode S_X(IGAPoint p,IGAPoint pS,PetscReal *F,PetscReal *US,void *ctx)
	{
		//This function creates the right hand side vector for an L2 projection of -S:X, alpha is read from file and added to it on the function call routine, as it's an input already stored on the correct mesh
		if (p->atboundary)
		{
			return 0.0;																		//No BCs for L2 Proj
		}

		//Definition of alternating tensor
		const PetscReal e[3][3][3]=
		{
			{{0.0,0.0,0.0},{0.0,0.0,1.0},{0.0,-1.0,0.0}},
			{{0.0,0.0,-1.0},{0.0,0.0,0.0},{1.0,0.0,0.0}},
			{{0.0,1.0,0.0},{-1.0,0.0,0.0},{0.0,0.0,0.0}}
		};

		PetscReal S[8];
		IGAPointFormValue(pS,US,&S[0]);														//This fills the values of S
		
		PetscReal fullS[3][3][3]={0};
		//Expand S^perp to full 3rd order form, only non-zero elements
		//The rounding kills spurious gradients due to numerics, could be removed, but makes S^perp look nicer
		fullS[0][0][0]=round(1.0e12*S[0])/1.0e12; fullS[0][0][1]=round(1.0e12*S[1])/1.0e12; 
		fullS[0][1][0]=round(1.0e12*S[2])/1.0e12; fullS[0][1][1]=round(1.0e12*S[3])/1.0e12;
		fullS[1][0][0]=round(1.0e12*S[4])/1.0e12; fullS[1][0][1]=round(1.0e12*S[5])/1.0e12; 
		fullS[1][1][0]=round(1.0e12*S[6])/1.0e12; fullS[1][1][1]=round(1.0e12*S[7])/1.0e12;

		PetscInt a,i,j,k,l;
		PetscInt nen = p->nen;																//Number of shape functions per element
		PetscInt dof = p->dof;																//Number of degrees of freedom per node

		PetscReal SX[3][3]={0};																//Stores -S:X

		for(i=0; i<3;i++)
		{
			for (j=0; j<3; j++)
			{
				for (k=0; k<3; k++)
				{
					for (l=0; l<3; l++)
					{
						SX[i][j]=SX[i][j]-fullS[i][k][l]*e[k][l][j];
					}
				}
			}
		}

		PetscReal sx[2]={0};
		sx[0]=SX[0][2]; sx[1]=SX[1][2];														//Due to S being 0 if i,j,k=3 and the properties of X, these are the only 2 non-zero components

		const PetscReal *N0;
		IGAPointGetShapeFuns(p,0,(const PetscReal**)&N0);
		PetscReal (*FF)[dof] = (PetscReal (*)[dof])F;
		for(a=0; a<nen; a++)
		{
			for(i=0; i<dof; i++) 
			{
				FF[a][i] = N0[a]*sx[i];
			}
		}
		return 0;
	}

	#undef  __FUNCT__
	#define __FUNCT__ "Sp_X"
	PetscErrorCode Sp_X(IGAPoint p,IGAPoint pSp,PetscReal *F,PetscReal *USp,void *ctx)
	{
		//This function creates the right hand side vector for an L2 projection of -S^perp:X, alpha is read from file and added to it on the function call routine, as it's an input already stored on the correct mesh
		if (p->atboundary)
		{
			return 0.0;																		//No BCs for L2 Proj
		}

		//Definition of alternating tensor
		const PetscReal e[3][3][3]=
		{
			{{0.0,0.0,0.0},{0.0,0.0,1.0},{0.0,-1.0,0.0}},
			{{0.0,0.0,-1.0},{0.0,0.0,0.0},{1.0,0.0,0.0}},
			{{0.0,1.0,0.0},{-1.0,0.0,0.0},{0.0,0.0,0.0}}
		};

		PetscReal Sp[8];
		IGAPointFormValue(pSp,USp,&Sp[0]);														//This fills the values of S
		
		PetscReal fullSp[3][3][3]={0};
		//Expand S^perp to full 3rd order form, only non-zero elements
		//The rounding kills spurious gradients due to numerics, could be removed, but makes S^perp look nicer
		fullSp[0][0][0]=round(1.0e12*Sp[0])/1.0e12; fullSp[0][0][1]=round(1.0e12*Sp[1])/1.0e12; 
		fullSp[0][1][0]=round(1.0e12*Sp[2])/1.0e12; fullSp[0][1][1]=round(1.0e12*Sp[3])/1.0e12;
		fullSp[1][0][0]=round(1.0e12*Sp[4])/1.0e12; fullSp[1][0][1]=round(1.0e12*Sp[5])/1.0e12; 
		fullSp[1][1][0]=round(1.0e12*Sp[6])/1.0e12; fullSp[1][1][1]=round(1.0e12*Sp[7])/1.0e12;

		PetscInt a,i,j,k,l;
		PetscInt nen = p->nen;																//Number of shape functions per element
		PetscInt dof = p->dof;																//Number of degrees of freedom per node

		PetscReal SpX[3][3]={0};															//Stores -Sp:X

		for(i=0; i<3;i++)
		{
			for (j=0; j<3; j++)
			{
				for (k=0; k<3; k++)
				{
					for (l=0; l<3; l++)
					{
						SpX[i][j]=SpX[i][j]-fullSp[i][k][l]*e[k][l][j];
					}
				}
			}
		}

		PetscReal sp[2]={0};
		sp[0]=SpX[0][2]; sp[1]=SpX[1][2];													//Due to Sp being 0 if i,j,k=3 and the properties of X, these are the only 2 non-zero components

		const PetscReal *N0;
		IGAPointGetShapeFuns(p,0,(const PetscReal**)&N0);
		PetscReal (*FF)[dof] = (PetscReal (*)[dof])F;
		for(a=0; a<nen; a++)
		{
			for(i=0; i<dof; i++) 
			{
				FF[a][i] = N0[a]*sp[i];
			}
		}
		return 0;
	}
//

//Helmholtz decomposition of -Ûp (same as Ûe), curl part
	#undef  __FUNCT__
	#define __FUNCT__ "curlChiUM"
	PetscErrorCode curlChiUM(IGAPoint p,PetscReal *K,void *ctx)
	{
		const PetscReal (*N1)[2];
		IGAPointGetShapeFuns(p,1,(const PetscReal**)&N1);

		PetscInt a,b,u,w,i,j,k,nen=p->nen, dof=p->dof;

		PetscReal (*KchiU)[dof][nen][dof] = (typeof(KchiU)) K;

		PetscReal dChiU[3][3][3]={0};
		PetscReal dv[3][3][3]={0};

		if (p->atboundary)
		{
			return 0;
		}
		else
		{
			for (a=0; a<nen; a++)
			{
				PetscReal Na_x = N1[a][0];
				PetscReal Na_y = N1[a][1];
				for (u=0; u<dof; u++) 
				{
					if(u==0)
					{
						dv[0][0][0]=Na_x; dv[0][0][1]=Na_y;
						dv[0][1][0]=0.0;  dv[0][1][1]=0.0;
						dv[1][0][0]=0.0;  dv[1][0][1]=0.0;
						dv[1][1][0]=0.0;  dv[1][1][1]=0.0;
					}
					if(u==1)
					{
						dv[0][0][0]=0.0;  dv[0][0][1]=0.0;
						dv[0][1][0]=Na_x; dv[0][1][1]=Na_y;
						dv[1][0][0]=0.0;  dv[1][0][1]=0.0;
						dv[1][1][0]=0.0;  dv[1][1][1]=0.0;
					}
					if(u==2)
					{
						dv[0][0][0]=0.0;  dv[0][0][1]=0.0;
						dv[0][1][0]=0.0;  dv[0][1][1]=0.0;
						dv[1][0][0]=Na_x; dv[1][0][1]=Na_y;
						dv[1][1][0]=0.0;  dv[1][1][1]=0.0;
					}
					if(u==3)
					{
						dv[0][0][0]=0.0;  dv[0][0][1]=0.0;
						dv[0][1][0]=0.0;  dv[0][1][1]=0.0;
						dv[1][0][0]=0.0;  dv[1][0][1]=0.0;
						dv[1][1][0]=Na_x; dv[1][1][1]=Na_y;
					}

					for (b=0; b<nen; b++)
					{
						PetscReal Nb_x = N1[b][0];
						PetscReal Nb_y = N1[b][1];
						for (w=0; w<dof; w++)
						{
							if(w==0)
							{
								dChiU[0][0][0]=Nb_x; dChiU[0][0][1]=Nb_y;
								dChiU[0][1][0]=0.0;  dChiU[0][1][1]=0.0;
								dChiU[1][0][0]=0.0;  dChiU[1][0][1]=0.0;
								dChiU[1][1][0]=0.0;  dChiU[1][1][1]=0.0;
							}
							if(w==1)
							{
								dChiU[0][0][0]=0.0;  dChiU[0][0][1]=0.0;
								dChiU[0][1][0]=Nb_x; dChiU[0][1][1]=Nb_y;
								dChiU[1][0][0]=0.0;  dChiU[1][0][1]=0.0;
								dChiU[1][1][0]=0.0;  dChiU[1][1][1]=0.0;
							}
							if(w==2)
							{
								dChiU[0][0][0]=0.0;  dChiU[0][0][1]=0.0;
								dChiU[0][1][0]=0.0;  dChiU[0][1][1]=0.0;
								dChiU[1][0][0]=Nb_x; dChiU[1][0][1]=Nb_y;
								dChiU[1][1][0]=0.0;  dChiU[1][1][1]=0.0;
							}
							if(w==3)
							{
								dChiU[0][0][0]=0.0;  dChiU[0][0][1]=0.0;
								dChiU[0][1][0]=0.0;  dChiU[0][1][1]=0.0;
								dChiU[1][0][0]=0.0;  dChiU[1][0][1]=0.0;
								dChiU[1][1][0]=Nb_x; dChiU[1][1][1]=Nb_y;
							}

							KchiU[a][u][b][w]=0.0;
							for(i=0;i<3;i++)
							{
								for(j=0;j<3;j++)
								{
									for(k=0;k<3;k++)
									{
										//Here the product of the two alternating tensors is expanded, this is about 2 orders of magnitude faster
										KchiU[a][u][b][w]+=dChiU[i][j][k]*dv[i][j][k]-dChiU[i][j][k]*dv[i][k][j]+dChiU[i][j][j]*dv[i][k][k];
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

	#undef  __FUNCT__
	#define __FUNCT__ "curlChiUF"
	PetscErrorCode curlChiUF(IGAPoint p,IGAPoint pAl,PetscReal *F,PetscReal *U, void *ctx)
	{
		//Generates right hand side only. Faster for iterations, as K is constant 

		const PetscReal (*N1)[2];
		IGAPointGetShapeFuns(p,1,(const PetscReal**)&N1);

		PetscInt a,u,i,j,k,l,nen=p->nen, dof=p->dof;

		//Definition of alternating tensor
		const PetscReal e[3][3][3]=
		{
			{{0.0,0.0,0.0},{0.0,0.0,1.0},{0.0,-1.0,0.0}},
			{{0.0,0.0,-1.0},{0.0,0.0,0.0},{1.0,0.0,0.0}},
			{{0.0,1.0,0.0},{-1.0,0.0,0.0},{0.0,0.0,0.0}}
		};

		PetscReal alfa[2];																		//Create array to receive Alfa
		IGAPointFormValue(pAl,U,&alfa[0]);														//This fills the values

		PetscReal fullAlfa[3][3]={0};
		fullAlfa[0][2]=alfa[0]; fullAlfa[1][2]=alfa[1];												//Expand Alfa to full 2nd order form, only non-zero elements

		PetscReal (*FchiU)[dof] = (PetscReal (*)[dof]) F;

		PetscReal dv[3][3][3]={0};

		if (p->atboundary)
		{
			return 0;
		}
		else
		{
			for (a=0; a<nen; a++)
			{
				PetscReal Na_x = N1[a][0];
				PetscReal Na_y = N1[a][1];
				for (u=0; u<dof; u++) 
				{
					if(u==0)
					{
						dv[0][0][0]=Na_x;	dv[0][0][1]=Na_y;
						dv[0][1][0]=0.0;	dv[0][1][1]=0.0;
						dv[1][0][0]=0.0;	dv[1][0][1]=0.0;
						dv[1][1][0]=0.0;	dv[1][1][1]=0.0;
					}
					if(u==1)
					{
						dv[0][0][0]=0.0;	dv[0][0][1]=0.0;
						dv[0][1][0]=Na_x;	dv[0][1][1]=Na_y;
						dv[1][0][0]=0.0;	dv[1][0][1]=0.0;
						dv[1][1][0]=0.0;	dv[1][1][1]=0.0;
					}
					if(u==2)
					{
						dv[0][0][0]=0.0;	dv[0][0][1]=0.0;
						dv[0][1][0]=0.0;	dv[0][1][1]=0.0;
						dv[1][0][0]=Na_x;	dv[1][0][1]=Na_y;
						dv[1][1][0]=0.0;	dv[1][1][1]=0.0;
					}
					if(u==3)
					{
						dv[0][0][0]=0.0;	dv[0][0][1]=0.0;
						dv[0][1][0]=0.0;	dv[0][1][1]=0.0;
						dv[1][0][0]=0.0;	dv[1][0][1]=0.0;
						dv[1][1][0]=Na_x;	dv[1][1][1]=Na_y;
					}

					for(i=0;i<3;i++)
					{
						for(j=0;j<3;j++)
						{
							for(k=0;k<3;k++)
							{
								for (l=0;l<3;l++)
								{
									FchiU[a][u]+=-fullAlfa[i][j]*e[j][k][l]*dv[i][l][k];
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

//System for z(0)  (I'm not separating this into M and F because it's only done one)
	#undef  __FUNCT__
	#define __FUNCT__ "Z0sys"
	PetscErrorCode Z0sys(IGAPoint p, IGAPoint pChi, PetscReal *K, PetscReal *F, PetscReal *UChi, void *ctx)
	{
		const PetscReal *N0,(*N1)[2];
		IGAPointGetShapeFuns(p,0,(const PetscReal**)&N0);									//Value of the shape functions
		IGAPointGetShapeFuns(p,1,(const PetscReal**)&N1);									//Derivatives of the shape functions
		PetscInt a,b,i,j,k,l,u,w,nen=p->nen, dof=p->dof;

		PetscReal x[2];																		//Vector of reals, size equal to problem's dimension
		IGAPointFormGeomMap(p,x);															//Fills x with the coordinates of p, Gauss's point

		//E and nu should come from AppCtx in the future
		//const PetscReal E=1.0; //2100000.0*9.81*10000.0;					//[Pa]
		//const PetscReal nu=0.33;
		//const PetscReal lambda=(E*nu)/((1.0+nu)*(1.0-2.0*nu));
		//const PetscReal mu=E/(2.0*(1.0+nu));

		//Change for G=1
		const PetscReal nu=0.33;
		const PetscReal mu=1.0;
		const PetscReal lambda=2.0*mu*nu/(1.0-2.0*nu);

		PetscReal Chi0[4];																	//Assign chi to a vector
		IGAPointFormValue(pChi,UChi,&Chi0[0]);

		PetscReal fullChi[3][3]={0};
		fullChi[0][0]=Chi0[0]; 	fullChi[0][1]=Chi0[1];
		fullChi[1][0]=Chi0[2]; 	fullChi[1][1]=Chi0[3];


		PetscReal (*Keq)[dof][nen][dof] = (typeof(Keq)) K;
		PetscReal (*Feq)[dof] = (PetscReal (*)[dof])F;

		//Creation of elasticity tensor
		PetscReal C[3][3][3][3]={0};														//Initialization of elastic tensor
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

		PetscReal dv[3][3]={0};
		PetscReal dz[3][3]={0};

		if (p->atboundary)
		{
			return 0;
		}
		else
		{
			for (a=0; a<nen; a++) 
			{
				PetscReal Na_x = N1[a][0];			PetscReal Na_y = N1[a][1];

				for (b=0; b<nen; b++) 
				{
					PetscReal Nb_x = N1[b][0];		PetscReal Nb_y = N1[b][1];

					for (i=0; i<dof; i++)
					{
						if (i==0)
						{
							dv[0][0]=Na_x; dv[0][1]=Na_y;
							dv[1][0]=0.0;  dv[1][1]=0.0;
						}
						else if(i==1)
						{
							dv[0][0]=0.0;  dv[0][1]=0.0;
							dv[1][0]=Na_x; dv[1][1]=Na_y;
						}

						for (j=0; j<dof; j++)
						{
							if (j==0)
							{
								dz[0][0]=Nb_x; dz[0][1]=Nb_y;
								dz[1][0]=0.0;  dz[1][1]=0.0;
							}
							else if(j==1)
							{
								dz[0][0]=0.0;  dz[0][1]=0.0;
								dz[1][0]=Nb_x; dz[1][1]=Nb_y;
							}

							Keq[a][i][b][j]=0.0;
							for (k=0; k<3; k++)
							{
								for (l=0; l<3; l++)
								{
									for (u=0; u<3; u++)
									{
										for (w=0; w<3; w++)
										{
											Keq[a][i][b][j]+=C[k][l][u][w]*dz[u][w]*dv[k][l];
										}
									}
								}
							}
						}
					}
				}
			}

			for(a=0;a<nen;a++)
			{
				PetscReal Na_x=N1[a][0];		PetscReal Na_y=N1[a][1];

				for (i=0; i<dof; i++)
				{
					if (i==0)
					{
						dv[0][0]=Na_x; dv[0][1]=Na_y;
						dv[1][0]=0.0;  dv[1][1]=0.0;
					}
					else if (i==1)
					{
						dv[0][0]=0.0;  dv[0][1]=0.0;
						dv[1][0]=Na_x; dv[1][1]=Na_y;
					}

					Feq[a][i] = 0.0;
					for (k=0;k<3;k++)
					{
						for (l=0;l<3;l++)
						{
							for (u=0;u<3;u++)
							{
								for (w=0;w<3;w++)
								{
									Feq[a][i]+=C[k][l][u][w]*(-fullChi[u][w])*dv[k][l];
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

//System for u (debug, everything seems ok, but a second look won't hurt)
	#undef  __FUNCT__
	#define __FUNCT__ "UsysM"
				 //Usys(point_u,point_chi,point_z,Kpoint_u,Fpoint_u,Chi_u,z_u,NULL);CHKERRQ(ierr);
	PetscErrorCode UsysM(IGAPoint p,PetscReal *K,void *ctx)
	{
		const PetscReal (*N1)[2];
		IGAPointGetShapeFuns(p,1,(const PetscReal**)&N1);

		PetscInt a,b,i,j,u,w,k,l,nen=p->nen, dof=p->dof;

		//Change for G=1
		const PetscReal nu=0.33;
		const PetscReal mu=1.0;
		const PetscReal lambda=2.0*mu*nu/(1.0-2.0*nu);

		//Creation of elasticity tensor
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

		PetscReal dv[3][3]={0};
		PetscReal du[3][3]={0};

		PetscReal (*Keq)[dof][nen][dof] = (typeof(Keq)) K;
		if (p->atboundary)
		{
			return 0;
		}
		else
		{
			for (a=0;a<nen;a++) 
			{
				PetscReal Na_x = N1[a][0];			PetscReal Na_y = N1[a][1];

				for (b=0;b<nen;b++) 
				{
					PetscReal Nb_x = N1[b][0];		PetscReal Nb_y = N1[b][1];

					for (i=0;i<dof;i++)
					{
						if (i==0)
						{
							dv[0][0]=Na_x; dv[0][1]=Na_y;
							dv[1][0]=0.0;  dv[1][1]=0.0;
						}
						else if(i==1)
						{
							dv[0][0]=0.0;  dv[0][1]=0.0;
							dv[1][0]=Na_x; dv[1][1]=Na_y;
						}

						for (j=0;j<dof;j++)
						{
							if (j==0)
							{
								du[0][0]=Nb_x; du[0][1]=Nb_y;
								du[1][0]=0.0;  du[1][1]=0.0;
							}
							else if(j==1)
							{
								du[0][0]=0.0;  du[0][1]=0.0;
								du[1][0]=Nb_x; du[1][1]=Nb_y;
							}

							Keq[a][i][b][j]=0.0;
							for (k=0;k<3;k++)
							{
								for (l=0;l<3;l++)
								{
									for (u=0;u<3;u++)
									{
										for(w=0;w<3;w++)
										{
											Keq[a][i][b][j]+=C[k][l][u][w]*du[u][w]*dv[k][l];
										}
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

	#undef  __FUNCT__
	#define __FUNCT__ "UsysF"
				 //Usys(point_u,point_chi,point_z,Kpoint_u,Fpoint_u,Chi_u,z_u,NULL);CHKERRQ(ierr);
	PetscErrorCode UsysF(IGAPoint p,IGAPoint pChi,IGAPoint pZ,PetscReal *F,PetscReal *UChi,PetscReal *UZ,void *ctx)
	{
		const PetscReal *N0,(*N1)[2];
		IGAPointGetShapeFuns(p,0,(const PetscReal**)&N0);
		IGAPointGetShapeFuns(p,1,(const PetscReal**)&N1);

		PetscInt a,i,j,u,w,k,l,nen=p->nen, dof=p->dof;

		//Change for G=1
		const PetscReal nu=0.33;
		const PetscReal mu=1.0;
		const PetscReal lambda=2.0*mu*nu/(1.0-2.0*nu);

		PetscReal Chi0[4];																	//Assign chi to a vector
		IGAPointFormValue(pChi,UChi,&Chi0[0]);

		PetscReal full_Chi[3][3]={0};
		full_Chi[0][0]=Chi0[0]; 	full_Chi[0][1]=Chi0[1];
		full_Chi[1][0]=Chi0[2]; 	full_Chi[1][1]=Chi0[3];
	
		PetscReal dz0[2][2];																//Same for partial derivatives of z
		IGAPointFormGrad (pZ,UZ,&dz0[0][0]);

		PetscReal full_dz[3][3]={0};														//Inflate to 3 indices, simplifies sums in loops
		full_dz[0][0]=dz0[0][0]; 	full_dz[0][1]=dz0[0][1];
		full_dz[1][0]=dz0[1][0]; 	full_dz[1][1]=dz0[1][1];

		//Creation of elasticity tensor
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

		PetscReal v[3]={0};
		PetscReal dv[3][3]={0};
		PetscReal n[3]={0};

		PetscReal (*Feq)[dof] = (PetscReal (*)[dof])F;

		if (p->atboundary)
		{
			PetscReal Sborde[3][3]={0};

			//Stress in boundary
			Sborde[0][0]=0.0;
			Sborde[0][1]=mu/1000.0;
			Sborde[1][0]=mu/1000.0;
			Sborde[1][1]=0.0;
			Sborde[0][2]=0.0; Sborde[1][2]=0.0; Sborde[2][0]=0.0; Sborde[2][1]=0.0; Sborde[2][2]=0.0;

			PetscInt dir  = p->boundary_id / 2;
			PetscInt side = p->boundary_id % 2;

			for (a=0;a<nen;a++)
			{
				PetscReal Na = N0[a];

				for (i=0;i<dof;i++)
				{
					if(i==0)
					{
						v[0]=Na; v[1]=0.0; v[2]=0.0;
					}
					if(i==1)
					{
						v[0]=0.0; v[1]=Na; v[2]=0.0;
					}

					Feq[a][i] = 0.0;
					if(dir==0)
					{
						if(side==0)
						{
							n[0]=-1.0; n[1]=0.0; n[2]=0.0;
							for (j=0;j<3;j++)
							{
								for(k=0;k<3;k++)
								{
									Feq[a][i]+=Sborde[j][k]*n[k]*v[j];
								}
							}
						}
						if(side==1)
						{
							n[0]=1.0; n[1]=0.0; n[2]=0.0;
							for (j=0;j<3;j++)
							{
								for(k=0;k<3;k++)
								{
									Feq[a][i]+=Sborde[j][k]*n[k]*v[j];
								}
							}
						}
					}
					if(dir==1)
					{
						if(side==0)
						{
							n[0]=0.0; n[1]=-1.0; n[2]=0.0;
							for (j=0;j<3;j++)
							{
								for(k=0;k<3;k++)
								{
									Feq[a][i]+=Sborde[j][k]*n[k]*v[j];
								}
							}
						}
						if(side==1)
						{
							n[0]=0.0; n[1]=1.0; n[2]=0.0;
							for (j=0;j<3;j++)
							{
								for(k=0;k<3;k++)
								{
									Feq[a][i]+=Sborde[j][k]*n[k]*v[j];
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
			for(a=0;a<nen;a++)
			{
				PetscReal Na_x = N1[a][0];			PetscReal Na_y = N1[a][1];

				for (i=0;i<dof;i++)
				{
					if (i==0)
					{
						dv[0][0]=Na_x; dv[0][1]=Na_y;
						dv[1][0]=0.0;  dv[1][1]=0.0;
					}
					else if (i==1)
					{
						dv[0][0]=0.0;  dv[0][1]=0.0;
						dv[1][0]=Na_x; dv[1][1]=Na_y;
					}

					Feq[a][i] = 0.0;
					for (k=0;k<3;k++)
					{
						for (l=0;l<3;l++)
						{
							for (u=0;u<3;u++)
							{
								for (w=0;w<3;w++)
								{
									Feq[a][i]+=C[k][l][u][w]*(full_dz[u][w]+full_Chi[u][w])*dv[k][l];
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

//System for L2 projection of V^{S} (debug both functions)
	#undef  __FUNCT__
	#define __FUNCT__ "VS"
	     		 //VS(pointVs,point_chi,point_z,point_u,point_S,KpointVs,FpointVs,Chi0Vs,Z0Vs,UVs,S0Vs,NULL);        
	PetscErrorCode VS(IGAPoint p,IGAPoint pChi,IGAPoint pZu,IGAPoint pU,IGAPoint pS,PetscReal *K,PetscReal *F,PetscReal *UChi,PetscReal *UZu,PetscReal *Uu,PetscReal *US,void *ctx)
	{
		//This functions does two things
		//a) Generates the L2 projection matrix
		//b) RHS is V^S with no u, as at t=0 there is only z and chi		
		const PetscReal *N0;
		IGAPointGetShapeFuns(p,0,(const PetscReal**)&N0);									//Value of the shape functions
		PetscInt a,b,i,j,k,l,m,w,nen=p->nen, dof=p->dof;

		//Change to consider G=1
		const PetscReal nu=0.33;
		const PetscReal mu=1.0;
		const PetscReal lambda=2.0*mu*nu/(1.0-2.0*nu);

		PetscReal mob=50.0;																	//Mobility coefficient such that V^S has units of velocity, for now chosen so V^S has a
																							//value that permits simulation in a reasonable number of timesteps.
		PetscReal S0[8];																	//Assign S to a vector
		IGAPointFormValue(pS,US,&S0[0]);																

		//The four non-zero components of Chi are stored as a vector, restore them to an array with the correct indexing for value and derivative
		PetscReal Chi0[4];																	//Array to contain the vector chi(0)
		IGAPointFormValue(pChi,UChi,&Chi0[0]);												//Assign chi to its container

		PetscReal dZ0[2][2];																//Array to store 3rd order partial derivative of z
		IGAPointFormGrad (pZu,UZu,&dZ0[0][0]);												//Calculate and store in array

		PetscReal fullChi[3][3]={0};
		fullChi[0][0]=Chi0[0]; 	fullChi[0][1]=Chi0[1];
		fullChi[1][0]=Chi0[2]; 	fullChi[1][1]=Chi0[3];

		PetscReal du[2][2];
		IGAPointFormGrad (pU,Uu,&du[0][0]);
		PetscReal full_du[3][3]={0};
		full_du[0][0]=du[0][0];
		full_du[0][1]=du[0][1];
		full_du[1][0]=du[1][0];
		full_du[1][1]=du[1][1];

		//Expanding z (and derivatives) to 3 components, more convenient for sums in for loops
		PetscReal full_dz[3][3]={0};
		full_dz[0][0]=dZ0[0][0];
		full_dz[0][1]=dZ0[0][1];
		full_dz[1][0]=dZ0[1][0];
		full_dz[1][1]=dZ0[1][1];

		//Inflate stored vectors to full tensor form
		PetscReal fullS[3][3][3]={0};
		fullS[0][0][0]=S0[0]; fullS[0][0][1]=S0[1];											//Expand S to full 3rd order form, only non-zero elements
		fullS[0][1][0]=S0[2]; fullS[0][1][1]=S0[3];
		fullS[1][0][0]=S0[4]; fullS[1][0][1]=S0[5];
		fullS[1][1][0]=S0[6]; fullS[1][1][1]=S0[7];

		PetscReal C[3][3][3][3]={0};														//Initialization of elastic tensor
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
					for (j=0;j<3;j++)
					{
						for(k=0;k<3;k++)
						{
							for(l=0;l<3;l++)
							{
								for(m=0;m<3;m++)
								{
									for(w=0;w<3;w++)
									{
										FVS[a][i]+=mob*C[j][k][l][m]*(full_du[l][m]-full_dz[l][m]-fullChi[l][m])*fullS[j][k][w]*v[w];
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

	#undef  __FUNCT__
	#define __FUNCT__ "VSF"
	PetscErrorCode VSF(IGAPoint p,IGAPoint pChi,IGAPoint pZu,IGAPoint pu,IGAPoint pS,PetscReal *F, PetscReal *UChi,PetscReal *UZu,PetscReal *Uu,PetscReal *US,void *ctx)		//When other fields like Pi or S (etc.) come into this, declare a IGAPoint pPi or pS and a new PetscReal *UPi or *US for each
	{
		//This function only generates the RHS, but considers sigma=C(grad(u)-grad(z)-chi)
		const PetscReal *N0;
		IGAPointGetShapeFuns(p,0,(const PetscReal**)&N0);									//Value of the shape functions
		PetscInt a,i,j,k,l,m,w,nen=p->nen, dof=p->dof;

		//Change to consider G=1
		const PetscReal nu=0.33;
		const PetscReal mu=1.0;
		const PetscReal lambda=2.0*mu*nu/(1.0-2.0*nu);

		PetscReal mob=50.0;																	//Mobility coefficient such that V^S has units of velocity, for now chosen so V^S has a
																							//value that permits simulation in a reasonable number of timesteps.
		PetscReal S0[8];																	//Assign S to a vector
		IGAPointFormValue(pS,US,&S0[0]);																

		//The four non-zero components of Chi are stored as a vector, restore them to an array with the correct indexing for value and derivative
		PetscReal Chi0[4];																	//Array to contain the vector chi(0)
		IGAPointFormValue(pChi,UChi,&Chi0[0]);												//Assign chi to its container

		PetscReal dZ0[2][2];																//Array to store 3rd order partial derivative of z
		IGAPointFormGrad (pZu,UZu,&dZ0[0][0]);												//Calculate and store in array

		PetscReal du[2][2];
		IGAPointFormGrad (pu,Uu,&du[0][0]);

		PetscReal fullChi[3][3]={0};
		fullChi[0][0]=Chi0[0]; 	fullChi[0][1]=Chi0[1];
		fullChi[1][0]=Chi0[2]; 	fullChi[1][1]=Chi0[3];

		//Expanding z (and derivatives) to 3 components, more convenient for sums in for loops
		PetscReal full_dz[3][3]={0};
		full_dz[0][0]=dZ0[0][0];
		full_dz[0][1]=dZ0[0][1];
		full_dz[1][0]=dZ0[1][0];
		full_dz[1][1]=dZ0[1][1];

		//Expanding u (and derivatives) to 3 components, more convenient for sums in for loops
		PetscReal full_du[3][3]={0};
		full_du[0][0]=du[0][0];
		full_du[0][1]=du[0][1];
		full_du[1][0]=du[1][0];
		full_du[1][1]=du[1][1];

		//Inflate stored vectors to full tensor form
		PetscReal fullS[3][3][3]={0};
		fullS[0][0][0]=S0[0]; fullS[0][0][1]=S0[1];											//Expand S to full 3rd order form, only non-zero elements
		fullS[0][1][0]=S0[2]; fullS[0][1][1]=S0[3];
		fullS[1][0][0]=S0[4]; fullS[1][0][1]=S0[5];
		fullS[1][1][0]=S0[6]; fullS[1][1][1]=S0[7];

		PetscReal C[3][3][3][3]={0};														//Initialization of elastic tensor
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

		PetscReal (*FVS)[dof] = (PetscReal (*)[dof])F;

		PetscReal v[3]={0.0};

		if (p->atboundary)
		{
			return 0;
		}
		else
		{
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
					for (j=0;j<3;j++)
					{
						for(k=0;k<3;k++)
						{
							for(l=0;l<3;l++)
							{
								for(m=0;m<3;m++)
								{
									for(w=0;w<3;w++)
									{
										FVS[a][i]+=mob*C[j][k][l][m]*(full_du[l][m]-full_dz[l][m]-fullChi[l][m])*fullS[j][k][w]*v[w];
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

//System for L2 projection of V^{S} smooth (Integrates Vs multiplied with xi(S), where xi(S)=1 if any component of S is different than 0, and 0 in other case)
	#undef  __FUNCT__
	#define __FUNCT__ "Int_Xi_Vs"
				 //Int_Xi_Vs(pointVs,point_S,PointInt1a,PointInt2a,Vs0Values,SVs,NULL);CHKERRQ(ierr);
	PetscErrorCode Int_Xi_Vs(IGAPoint pV,IGAPoint pS,PetscReal *FInt1a,PetscReal *FInt2a,PetscReal *VS,PetscReal *US,void *ctx)
	{
		PetscInt dof=pV->dof;
		PetscInt i;

		PetscReal Vs[2];																	//Create array to receive Alfa
		IGAPointFormValue(pV,VS,&Vs[0]);													//This fills the values

		PetscReal S[8];																		//Create array to receive Alfa
		IGAPointFormValue(pS,US,&S[0]);														//This fills the values

		PetscReal xi=0.0;
		for (i=0;i<8;i++)
		{
			if (S[i]>0)
			{
				xi=1.0;
			}	
		}
		
		if (S[0]<0 || S[1]<0 || S[2]<0 || S[3]<0 || S[4]<0 || S[5]<0 || S[6]<0 || S[7]<0)
		{
			xi=1.0;
		}

		PetscReal (*FI1a)[dof] = (PetscReal (*)[dof])FInt1a;		//This vector will have just VS[0]
		PetscReal (*FI2a)[dof] = (PetscReal (*)[dof])FInt2a;		//This vector will have just VS[1]

		if (pV->atboundary)
		{
			return 0;
		}
		else
		{
			FI1a[0][0]+=xi*Vs[0];
			FI2a[0][1]+=xi*Vs[1];
		}
		return 0;
	}

	#undef  __FUNCT__
	#define __FUNCT__ "ProjVS"		//Don't use this one anymore? Use the direct integral value?
	//ierr = ProjVS(pointVs,point_S,VsSmoothPoint,SVs,Int1a,Int2a,IntS,absSmax,NULL);CHKERRQ(ierr);
	PetscErrorCode ProjVS(IGAPoint pV,IGAPoint pS,PetscReal *FS1,PetscReal *US,PetscReal Vs1,PetscReal Vs2,PetscReal IntS,PetscReal Smax, void *ctx)		//When other fields like Pi or S (etc.) come into this, declare a IGAPoint pPi or pS and a new PetscReal *UPi or *US for each
	{
		const PetscReal *N0;
		IGAPointGetShapeFuns(pV,0,(const PetscReal**)&N0);									//Value of the shape functions
		PetscInt a,i,nen=pV->nen,dof=pV->dof;

		PetscReal S[8];																		//Create array to receive Alfa
		IGAPointFormValue(pS,US,&S[0]);														//This fills the values

		PetscReal xi=0.0;
		PetscInt ind=0;
		for (i=0;i<8;i++)
		{
			if (S[i]>0.15 || S[i]<-0.15)
			{
				ind=i,
				xi=1.0;
			}	
		}

		PetscReal (*FVs1)[dof] = (PetscReal (*)[dof])FS1;

		if (pV->atboundary)
		{
			return 0;
		}
		else
		{
			for(a=0;a<nen;a++)
			{
				FVs1[a][0]+=xi*(Vs1/IntS)*fabs(S[ind])/Smax*N0[a];
				FVs1[a][1]+=xi*(Vs2/IntS)*fabs(S[ind])/Smax*N0[a];
			}
		}
		return 0;
	}
//

//System for V^{alpha}
	#undef  __FUNCT__
	#define __FUNCT__ "Valpha"
				 //Valpha(pointVa,point_chi,point_z,point_u,point_Al_hat,point_S,KpointVa,FpointVa,Chi0Va,Z0Va,UVa,SVa,Al0Va,NULL);CHKERRQ(ierr);
	PetscErrorCode Valpha(IGAPoint p,IGAPoint pChi,IGAPoint pZu,IGAPoint pU,IGAPoint pAl,IGAPoint pS,PetscReal *K,PetscReal *F,PetscReal *Chi,PetscReal *Zu,PetscReal *Uu,PetscReal *UAl,PetscReal *US,void *ctx)
	{
		//This functions does two things
		//a) Generates the L2 projection matrix
		//b) RHS is V^alpha with no u, as at t=0 there is only z and chi		
		const PetscReal *N0;
		IGAPointGetShapeFuns(p,0,(const PetscReal**)&N0);									//Value of the shape functions
		PetscInt a,b,c,d,i,j,k,l,m,nen=p->nen, dof=p->dof;

		//Change to consider G=1
		const PetscReal nu=0.33;
		const PetscReal mu=1.0;
		const PetscReal lambda=2.0*mu*nu/(1.0-2.0*nu);

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

		PetscReal alfa[2];																	//Create array to receive Alfa
		IGAPointFormValue(pAl,UAl,&alfa[0]);												//This fills the values

		PetscReal chi[4];																	//Array to contain the vector chi(0)
		IGAPointFormValue(pChi,Chi,&chi[0]);												//Assign chi to its container

		PetscReal d_Z0[2][2];
		IGAPointFormGrad (pZu,Zu,&d_Z0[0][0]);

		PetscReal du[2][2];
		IGAPointFormGrad (pU,Uu,&du[0][0]);
		PetscReal full_du[3][3]={0};
		full_du[0][0]=du[0][0];
		full_du[0][1]=du[0][1];
		full_du[1][0]=du[1][0];
		full_du[1][1]=du[1][1];

		//Inflate stored vectors to full tensor form
		PetscReal fullAlfa[3][3]={0};
		fullAlfa[0][2]=alfa[0]; fullAlfa[1][2]=alfa[1];										//Expand Alfa to full 2nd order form, only non-zero elements

		//The four non-zero components of Chi are stored as a vector, restore them to an array with the correct indexing for value and derivative
		PetscReal fullChi[3][3]={0};
		fullChi[0][0]=chi[0]; 	fullChi[0][1]=chi[1];
		fullChi[1][0]=chi[2]; 	fullChi[1][1]=chi[3];

		//Expanding z (and derivatives) to 3 components, more convenient for sums in for loops
		PetscReal fulld_z[3][3]={0};
		fulld_z[0][0]=d_Z0[0][0]; fulld_z[0][1]=d_Z0[0][1];
		fulld_z[1][0]=d_Z0[1][0]; fulld_z[1][1]=d_Z0[1][1];

		PetscReal S[8];																		//Create array to receive S
		IGAPointFormValue(pS,US,&S[0]);														//This fills the values of S (remember that S has 8 non zero components in 2D)

		PetscReal fullS[3][3][3]={0};
		fullS[0][0][0]=S[0]; fullS[0][0][1]=S[1];											//Expand S to full 3rd order form, only non-zero elements
		fullS[0][1][0]=S[2]; fullS[0][1][1]=S[3];
		fullS[1][0][0]=S[4]; fullS[1][0][1]=S[5];
		fullS[1][1][0]=S[6]; fullS[1][1][1]=S[7];

		//This transforms \alfa into \tilde{\alpha}=\alpha-S:X
		for (j=0;j<3;j++)
		{
			for (k=0;k<3;k++)
			{
				for (l=0;l<3;l++)
				{
					for (m=0;m<3;m++)
					{
						fullAlfa[j][k]=fullAlfa[j][k]-fullS[j][l][m]*e[k][l][m];
					}
				}
			}
		}

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
											FCS[a][i]+=C[j][k][l][m]*(full_du[l][m]-fulld_z[l][m]-fullChi[l][m])*e[k][c][d]*fullAlfa[j][c]*v[d];
										}
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

	#undef  __FUNCT__
	#define __FUNCT__ "ValphaF"
	//PetscErrorCode Stress(IGAPoint p,IGAPoint pU, IGAPoint pHs,IGAPoint pChi,IGAPoint pZu,PetscReal *K,PetscReal *F,PetscReal *U,PetscReal *HS, PetscReal *Chi,PetscReal *Zu,void *ctx)		//When other fields like Pi or S (etc.) come into this, declare a IGAPoint pPi or pS and a new PescReal *UPi or *US for each
	PetscErrorCode ValphaF(IGAPoint p,IGAPoint pChi,IGAPoint pu,IGAPoint pZu,IGAPoint pAl,IGAPoint pS,PetscReal *F,PetscReal *Chi,PetscReal *Zu,PetscReal *US,PetscReal *Uu,PetscReal *UAl,void *ctx)		//When other fields like Pi or S (etc.) come into this, declare a IGAPoint pPi or pS and a new PetscReal *UPi or *US for each
	{
		//This functions generates only right hand side and considers u. The stiffness matrix is reused from the previous case
		const PetscReal *N0;
		IGAPointGetShapeFuns(p,0,(const PetscReal**)&N0);									//Value of the shape functions
		PetscInt a,c,d,i,j,k,l,m,nen=p->nen, dof=p->dof;

		//Change to consider G=1
		const PetscReal nu=0.33;
		const PetscReal mu=1.0;
		const PetscReal lambda=2.0*mu*nu/(1.0-2.0*nu);

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

		PetscReal alfa[2];																	//Create array to receive Alfa
		IGAPointFormValue(pAl,UAl,&alfa[0]);													//This fills the values

		PetscReal chi[4];																	//Array to contain the vector chi(0)
		IGAPointFormValue(pChi,Chi,&chi[0]);												//Assign chi to its container

		PetscReal d_Z0[2][2];																//Same for its gradient
		IGAPointFormGrad (pZu,Zu,&d_Z0[0][0]);												//Same for the gradient

		PetscReal d_u[2][2];
		IGAPointFormGrad (pu,Uu,&d_u[0][0]);												//Same for the gradient
		
		PetscReal S[8];																		//Create array to receive S
		IGAPointFormValue(pS,US,&S[0]);														//This fills the values of S (remember that S has 8 non zero components in 2D)

		//Inflate stored vectors to full tensor form
		PetscReal fullAlfa[3][3]={0};
		fullAlfa[0][2]=alfa[0]; fullAlfa[1][2]=alfa[1];										//Expand Alfa to full 2nd order form, only non-zero elements

		//The four non-zero components of Chi are stored as a vector, restore them to an array with the correct indexing for value and derivative
		PetscReal fullChi[3][3]={0};
		fullChi[0][0]=chi[0]; 	fullChi[0][1]=chi[1];
		fullChi[1][0]=chi[2]; 	fullChi[1][1]=chi[3];

		//Expanding z (and derivatives) to 3 components, more convenient for sums in for loops
		PetscReal fulld_z[3][3]={0};
		fulld_z[0][0]=d_Z0[0][0]; fulld_z[0][1]=d_Z0[0][1];
		fulld_z[1][0]=d_Z0[1][0]; fulld_z[1][1]=d_Z0[1][1];

		//Expanding u (and derivatives) to 3 components, more convenient for sums in for loops
		PetscReal full_du[3][3]={0};
		full_du[0][0]=d_u[0][0]; full_du[0][1]=d_u[0][1];
		full_du[1][0]=d_u[1][0]; full_du[1][1]=d_u[1][1];

		PetscReal fullS[3][3][3]={0};
		fullS[0][0][0]=S[0]; fullS[0][0][1]=S[1];											//Expand S to full 3rd order form, only non-zero elements
		fullS[0][1][0]=S[2]; fullS[0][1][1]=S[3];
		fullS[1][0][0]=S[4]; fullS[1][0][1]=S[5];
		fullS[1][1][0]=S[6]; fullS[1][1][1]=S[7];

		for (j=0;j<3;j++)
		{
			for (k=0;k<3;k++)
			{
				for (l=0;l<3;l++)
				{
					for (m=0;m<3;m++)
					{
						fullAlfa[j][k]=fullAlfa[j][k]-fullS[j][l][m]*e[k][l][m];
					}
				}
			}
		}

		PetscReal (*FCS)[dof] = (PetscReal (*)[dof])F;

		PetscReal v[3]={0};

		if (p->atboundary)
		{
			return 0;
		}
		else
		{
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
											FCS[a][i]+=C[j][k][l][m]*(full_du[l][m]-fulld_z[l][m]-fullChi[l][m])*e[k][c][d]*fullAlfa[j][c]*v[d];
										}
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

//System for updating S using equation for Sdot (debug)
	//To do list
	//This function does not take Pi as an input yet, still to implement that

	#undef  __FUNCT__
	#define __FUNCT__ "SdotFunc"
	//SdotFunc(point_Sdot,point_S,Kpoint_Sdot,Fpoint_Sdot,Int1a,Int2a,IntS,S0Sdot,&user);CHKERRQ(ierr);
	PetscErrorCode SdotFunc(IGAPoint p, IGAPoint pS, PetscReal *K, PetscReal *F, PetscReal Int1a, PetscReal Int2a, PetscReal normS, PetscReal *US, void *ctx)
	{
		PetscInt a,b,i,j,k,l,u,w,dof=p->dof,nen=p->nen;	

		const PetscReal *N0,(*N1)[2];
		IGAPointGetShapeFuns(p,0,(const PetscReal**)&N0);									//Value of the shape functions
		IGAPointGetShapeFuns(p,1,(const PetscReal**)&N1);									//Derivatives of the shape functions

		AppCtx    *user  = (AppCtx *)ctx;
		PetscReal dt     = user->dt;

		//PetscReal Vs[2],dVs[2][2];															//Seems that if I create this one with size 3, next line gives problems
		//IGAPointFormValue(pVs,UVs,&Vs[0]);													//This fills the values
		//IGAPointFormGrad (pVs,UVs,&dVs[0][0]);												//This fills the values

		PetscReal fullVs[3]={0};
		PetscReal full_dVs[3][3]={0};
		fullVs[0]=Int1a/normS; fullVs[1]=-Int2a/normS;
		fullVs[0]=0.0;		 fullVs[1]=-1.0;			//THIS IS FOR DEBUGGING!!

		full_dVs[0][0]=0.0; full_dVs[0][1]=0.0;
		full_dVs[1][0]=0.0; full_dVs[1][1]=0.0;

		//PetscReal S[8];																			//Create array to receive S and grad(S) (backwards Euler)
		PetscReal S[8], dS[8][2];																//Create array to receive S and grad(S) (Trapezoidal rule)
		IGAPointFormValue(pS,US,&S[0]);															//This fills the values of S (remember that S has 8 non zero components in 2D)
		IGAPointFormGrad (pS,US,&dS[0][0]);													//This fills the values of the gradient

		PetscReal fullS[3][3][3]={0};
		fullS[0][0][0]=S[0]; fullS[0][0][1]=S[1];												//Expand S to full 3rd order form, only non-zero elements
		fullS[0][1][0]=S[2]; fullS[0][1][1]=S[3];
		fullS[1][0][0]=S[4]; fullS[1][0][1]=S[5];
		fullS[1][1][0]=S[6]; fullS[1][1][1]=S[7];

		PetscReal fulld_S[3][3][3][3]={0};													//Expand grad(S) to full tensor order form, only non-zero elements
		fulld_S[0][0][0][0]=dS[0][0]; fulld_S[0][0][0][1]=dS[0][1]; 
		fulld_S[0][0][1][0]=dS[1][0]; fulld_S[0][0][1][1]=dS[1][1]; 
		fulld_S[0][1][0][0]=dS[2][0]; fulld_S[0][1][0][1]=dS[2][1]; 
		fulld_S[0][1][1][0]=dS[3][0]; fulld_S[0][1][1][1]=dS[3][1];
		fulld_S[1][0][0][0]=dS[4][0]; fulld_S[1][0][0][1]=dS[4][1]; 
		fulld_S[1][0][1][0]=dS[5][0]; fulld_S[1][0][1][1]=dS[5][1]; 
		fulld_S[1][1][0][0]=dS[6][0]; fulld_S[1][1][0][1]=dS[6][1]; 
		fulld_S[1][1][1][0]=dS[7][0]; fulld_S[1][1][1][1]=dS[7][1];

		PetscReal St[3][3][3]={0};
		PetscReal dSt[3][3][3][3]={0};
		PetscReal v[3][3][3]={0};
		PetscReal dv[3][3][3][3]={0};

		PetscReal (*KSdot)[dof][nen][dof] = (typeof(KSdot))K;
		PetscReal (*FSdot)[dof] = (PetscReal (*)[dof])F;
		if (p->atboundary)
		{
			return 0;
		}
		else
		{
			for (a=0; a<nen; a++) 
			{
				PetscReal Na=N0[a];
				PetscReal Na_x = N1[a][0];			PetscReal Na_y = N1[a][1];

				for (i=0; i<dof; i++)
				{
					if (i==0)
					{
							v[0][0][0]=Na;		v[0][0][1]=0.0;		v[0][1][0]=0.0;	v[0][1][1]=0.0;	v[1][0][0]=0.0;	v[1][0][1]=0.0;	v[1][1][0]=0.0;	v[1][1][1]=0.0;
							dv[0][0][0][0]=Na_x;	dv[0][0][0][1]=Na_y;	dv[0][0][0][2]=0.0;
							dv[0][0][1][0]=0.0;		dv[0][0][1][1]=0.0;		dv[0][0][1][2]=0.0;
							dv[0][1][0][0]=0.0;		dv[0][1][0][1]=0.0;		dv[0][1][0][2]=0.0;
							dv[0][1][1][0]=0.0;		dv[0][1][1][1]=0.0;		dv[0][1][1][2]=0.0;
							dv[1][0][0][0]=0.0;		dv[1][0][0][1]=0.0;		dv[1][0][0][2]=0.0;
							dv[1][0][1][0]=0.0;		dv[1][0][1][1]=0.0;		dv[1][0][1][2]=0.0;
							dv[1][1][0][0]=0.0;		dv[1][1][0][1]=0.0;		dv[1][1][0][2]=0.0;
							dv[1][1][1][0]=0.0;		dv[1][1][1][1]=0.0;		dv[1][1][1][2]=0.0;
					}
					else if(i==1)
					{
							v[0][0][0]=0.0;		v[0][0][1]=Na;		v[0][1][0]=0.0;	v[0][1][1]=0.0;	v[1][0][0]=0.0;	v[1][0][1]=0.0;	v[1][1][0]=0.0;	v[1][1][1]=0.0;
							dv[0][0][0][0]=0.0;		dv[0][0][0][1]=0.0;		dv[0][0][0][2]=0.0;
							dv[0][0][1][0]=Na_x;	dv[0][0][1][1]=Na_y;	dv[0][0][1][2]=0.0;
							dv[0][1][0][0]=0.0;		dv[0][1][0][1]=0.0;		dv[0][1][0][2]=0.0;
							dv[0][1][1][0]=0.0;		dv[0][1][1][1]=0.0;		dv[0][1][1][2]=0.0;
							dv[1][0][0][0]=0.0;		dv[1][0][0][1]=0.0;		dv[1][0][0][2]=0.0;
							dv[1][0][1][0]=0.0;		dv[1][0][1][1]=0.0;		dv[1][0][1][2]=0.0;
							dv[1][1][0][0]=0.0;		dv[1][1][0][1]=0.0;		dv[1][1][0][2]=0.0;
							dv[1][1][1][0]=0.0;		dv[1][1][1][1]=0.0;		dv[1][1][1][2]=0.0;
					}
					else if(i==2)
					{
							v[0][0][0]=0.0;		v[0][0][1]=0.0;		v[0][1][0]=Na;	v[0][1][1]=0.0;	v[1][0][0]=0.0;	v[1][0][1]=0.0;	v[1][1][0]=0.0;	v[1][1][1]=0.0;
							dv[0][0][0][0]=0.0;		dv[0][0][0][1]=0.0;		dv[0][0][0][2]=0.0;
							dv[0][0][1][0]=0.0;		dv[0][0][1][1]=0.0;		dv[0][0][1][2]=0.0;
							dv[0][1][0][0]=Na_x;	dv[0][1][0][1]=Na_x;	dv[0][1][0][2]=0.0;
							dv[0][1][1][0]=0.0;		dv[0][1][1][1]=0.0;		dv[0][1][1][2]=0.0;
							dv[1][0][0][0]=0.0;		dv[1][0][0][1]=0.0;		dv[1][0][0][2]=0.0;
							dv[1][0][1][0]=0.0;		dv[1][0][1][1]=0.0;		dv[1][0][1][2]=0.0;
							dv[1][1][0][0]=0.0;		dv[1][1][0][1]=0.0;		dv[1][1][0][2]=0.0;
							dv[1][1][1][0]=0.0;		dv[1][1][1][1]=0.0;		dv[1][1][1][2]=0.0;
					}
					else if(i==3)
					{
							v[0][0][0]=0.0;		v[0][0][1]=0.0;		v[0][1][0]=0.0;	v[0][1][1]=Na;	v[1][0][0]=0.0;	v[1][0][1]=0.0;	v[1][1][0]=0.0;	v[1][1][1]=0.0;
							dv[0][0][0][0]=0.0;		dv[0][0][0][1]=0.0;		dv[0][0][0][2]=0.0;
							dv[0][0][1][0]=0.0;		dv[0][0][1][1]=0.0;		dv[0][0][1][2]=0.0;
							dv[0][1][0][0]=0.0;		dv[0][1][0][1]=0.0;		dv[0][1][0][2]=0.0;
							dv[0][1][1][0]=Na_x;	dv[0][1][1][1]=Na_y;	dv[0][1][1][2]=0.0;
							dv[1][0][0][0]=0.0;		dv[1][0][0][1]=0.0;		dv[1][0][0][2]=0.0;
							dv[1][0][1][0]=0.0;		dv[1][0][1][1]=0.0;		dv[1][0][1][2]=0.0;
							dv[1][1][0][0]=0.0;		dv[1][1][0][1]=0.0;		dv[1][1][0][2]=0.0;
							dv[1][1][1][0]=0.0;		dv[1][1][1][1]=0.0;		dv[1][1][1][2]=0.0;
					}
					else if(i==4)
					{
							v[0][0][0]=0.0;		v[0][0][1]=0.0;		v[0][1][0]=0.0;	v[0][1][1]=0.0;	v[1][0][0]=Na;	v[1][0][1]=0.0;	v[1][1][0]=0.0;	v[1][1][1]=0.0;
							dv[0][0][0][0]=0.0;		dv[0][0][0][1]=0.0;		dv[0][0][0][2]=0.0;
							dv[0][0][1][0]=0.0;		dv[0][0][1][1]=0.0;		dv[0][0][1][2]=0.0;
							dv[0][1][0][0]=0.0;		dv[0][1][0][1]=0.0;		dv[0][1][0][2]=0.0;
							dv[0][1][1][0]=0.0;		dv[0][1][1][1]=0.0;		dv[0][1][1][2]=0.0;
							dv[1][0][0][0]=Na_x;	dv[1][0][0][1]=Na_y;	dv[1][0][0][2]=0.0;
							dv[1][0][1][0]=0.0;		dv[1][0][1][1]=0.0;		dv[1][0][1][2]=0.0;
							dv[1][1][0][0]=0.0;		dv[1][1][0][1]=0.0;		dv[1][1][0][2]=0.0;
							dv[1][1][1][0]=0.0;		dv[1][1][1][1]=0.0;		dv[1][1][1][2]=0.0;
					}
					else if(i==5)
					{
							v[0][0][0]=0.0;		v[0][0][1]=0.0;		v[0][1][0]=0.0;	v[0][1][1]=0.0;	v[1][0][0]=0.0;	v[1][0][1]=Na;	v[1][1][0]=0.0;	v[1][1][1]=0.0;
							dv[0][0][0][0]=0.0;		dv[0][0][0][1]=0.0;		dv[0][0][0][2]=0.0;
							dv[0][0][1][0]=0.0;		dv[0][0][1][1]=0.0;		dv[0][0][1][2]=0.0;
							dv[0][1][0][0]=0.0;		dv[0][1][0][1]=0.0;		dv[0][1][0][2]=0.0;
							dv[0][1][1][0]=0.0;		dv[0][1][1][1]=0.0;		dv[0][1][1][2]=0.0;
							dv[1][0][0][0]=0.0;		dv[1][0][0][1]=0.0;		dv[1][0][0][2]=0.0;
							dv[1][0][1][0]=Na_x;	dv[1][0][1][1]=Na_y;	dv[1][0][1][2]=0.0;
							dv[1][1][0][0]=0.0;		dv[1][1][0][1]=0.0;		dv[1][1][0][2]=0.0;
							dv[1][1][1][0]=0.0;		dv[1][1][1][1]=0.0;		dv[1][1][1][2]=0.0;
					}
					else if(i==6)
					{
							v[0][0][0]=0.0;		v[0][0][1]=0.0;		v[0][1][0]=0.0;	v[0][1][1]=0.0;	v[1][0][0]=0.0;	v[1][0][1]=0.0;	v[1][1][0]=Na;	v[1][1][1]=0.0;
							dv[0][0][0][0]=0.0;		dv[0][0][0][1]=0.0;		dv[0][0][0][2]=0.0;
							dv[0][0][1][0]=0.0;		dv[0][0][1][1]=0.0;		dv[0][0][1][2]=0.0;
							dv[0][1][0][0]=0.0;		dv[0][1][0][1]=0.0;		dv[0][1][0][2]=0.0;
							dv[0][1][1][0]=0.0;		dv[0][1][1][1]=0.0;		dv[0][1][1][2]=0.0;
							dv[1][0][0][0]=0.0;		dv[1][0][0][1]=0.0;		dv[1][0][0][2]=0.0;
							dv[1][0][1][0]=0.0;		dv[1][0][1][1]=0.0;		dv[1][0][1][2]=0.0;
							dv[1][1][0][0]=Na_x;	dv[1][1][0][1]=Na_y;	dv[1][1][0][2]=0.0;
							dv[1][1][1][0]=0.0;		dv[1][1][1][1]=0.0;		dv[1][1][1][2]=0.0;
					}
					else if(i==7)
					{
							v[0][0][0]=0.0;		v[0][0][1]=0.0;		v[0][1][0]=0.0;	v[0][1][1]=0.0;	v[1][0][0]=0.0;	v[1][0][1]=0.0;	v[1][1][0]=0.0;	v[1][1][1]=Na;
							dv[0][0][0][0]=0.0;		dv[0][0][0][1]=0.0;		dv[0][0][0][2]=0.0;
							dv[0][0][1][0]=0.0;		dv[0][0][1][1]=0.0;		dv[0][0][1][2]=0.0;
							dv[0][1][0][0]=0.0;		dv[0][1][0][1]=0.0;		dv[0][1][0][2]=0.0;
							dv[0][1][1][0]=0.0;		dv[0][1][1][1]=0.0;		dv[0][1][1][2]=0.0;
							dv[1][0][0][0]=0.0;		dv[1][0][0][1]=0.0;		dv[1][0][0][2]=0.0;
							dv[1][0][1][0]=0.0;		dv[1][0][1][1]=0.0;		dv[1][0][1][2]=0.0;
							dv[1][1][0][0]=0.0;		dv[1][1][0][1]=0.0;		dv[1][1][0][2]=0.0;
							dv[1][1][1][0]=Na_x;	dv[1][1][1][1]=Na_y;	dv[1][1][1][2]=0.0;
					}

					for (b=0; b<nen; b++) 
					{
						PetscReal Nb=N0[b];
						PetscReal Nb_x = N1[b][0];		PetscReal Nb_y = N1[b][1];

						for (j=0; j<dof; j++)
						{
							if (j==0)
							{
								St[0][0][0]=Nb;		St[0][0][1]=0.0;		St[0][1][0]=0.0;	St[0][1][1]=0.0;	St[1][0][0]=0.0;	St[1][0][1]=0.0;	St[1][1][0]=0.0;	St[1][1][1]=0.0;
								dSt[0][0][0][0]=Nb_x;	dSt[0][0][0][1]=Nb_y;	dSt[0][0][0][2]=0.0;
								dSt[0][0][1][0]=0.0;	dSt[0][0][1][1]=0.0;	dSt[0][0][1][2]=0.0;
								dSt[0][1][0][0]=0.0;	dSt[0][1][0][1]=0.0;	dSt[0][1][0][2]=0.0;
								dSt[0][1][1][0]=0.0;	dSt[0][1][1][1]=0.0;	dSt[0][1][1][2]=0.0;
								dSt[1][0][0][0]=0.0;	dSt[1][0][0][1]=0.0;	dSt[1][0][0][2]=0.0;
								dSt[1][0][1][0]=0.0;	dSt[1][0][1][1]=0.0;	dSt[1][0][1][2]=0.0;
								dSt[1][1][0][0]=0.0;	dSt[1][1][0][1]=0.0;	dSt[1][1][0][2]=0.0;
								dSt[1][1][1][0]=0.0;	dSt[1][1][1][1]=0.0;	dSt[1][1][1][2]=0.0;
							}
							else if(j==1)
							{
								St[0][0][0]=0.0;	St[0][0][1]=Nb;			St[0][1][0]=0.0;	St[0][1][1]=0.0;	St[1][0][0]=0.0;	St[1][0][1]=0.0;	St[1][1][0]=0.0;	St[1][1][1]=0.0;
								dSt[0][0][0][0]=0.0;	dSt[0][0][0][1]=0.0;	dSt[0][0][0][2]=0.0;
								dSt[0][0][1][0]=Nb_x;	dSt[0][0][1][1]=Nb_y;	dSt[0][0][1][2]=0.0;
								dSt[0][1][0][0]=0.0;	dSt[0][1][0][1]=0.0;	dSt[0][1][0][2]=0.0;
								dSt[0][1][1][0]=0.0;	dSt[0][1][1][1]=0.0;	dSt[0][1][1][2]=0.0;
								dSt[1][0][0][0]=0.0;	dSt[1][0][0][1]=0.0;	dSt[1][0][0][2]=0.0;
								dSt[1][0][1][0]=0.0;	dSt[1][0][1][1]=0.0;	dSt[1][0][1][2]=0.0;
								dSt[1][1][0][0]=0.0;	dSt[1][1][0][1]=0.0;	dSt[1][1][0][2]=0.0;
								dSt[1][1][1][0]=0.0;	dSt[1][1][1][1]=0.0;	dSt[1][1][1][2]=0.0;
							}
							else if(j==2)
							{
								St[0][0][0]=0.0;	St[0][0][1]=0.0;		St[0][1][0]=Nb;		St[0][1][1]=0.0;	St[1][0][0]=0.0;	St[1][0][1]=0.0;	St[1][1][0]=0.0;	St[1][1][1]=0.0;
								dSt[0][0][0][0]=0.0;	dSt[0][0][0][1]=0.0;	dSt[0][0][0][2]=0.0;
								dSt[0][0][1][0]=0.0;	dSt[0][0][1][1]=0.0;	dSt[0][0][1][2]=0.0;
								dSt[0][1][0][0]=Nb_x;	dSt[0][1][0][1]=Nb_x;	dSt[0][1][0][2]=0.0;
								dSt[0][1][1][0]=0.0;	dSt[0][1][1][1]=0.0;	dSt[0][1][1][2]=0.0;
								dSt[1][0][0][0]=0.0;	dSt[1][0][0][1]=0.0;	dSt[1][0][0][2]=0.0;
								dSt[1][0][1][0]=0.0;	dSt[1][0][1][1]=0.0;	dSt[1][0][1][2]=0.0;
								dSt[1][1][0][0]=0.0;	dSt[1][1][0][1]=0.0;	dSt[1][1][0][2]=0.0;
								dSt[1][1][1][0]=0.0;	dSt[1][1][1][1]=0.0;	dSt[1][1][1][2]=0.0;
							}
							else if(j==3)
							{
								St[0][0][0]=0.0;	St[0][0][1]=0.0;		St[0][1][0]=0.0;	St[0][1][1]=Nb;		St[1][0][0]=0.0;	St[1][0][1]=0.0;	St[1][1][0]=0.0;	St[1][1][1]=0.0;
								dSt[0][0][0][0]=0.0;	dSt[0][0][0][1]=0.0;	dSt[0][0][0][2]=0.0;
								dSt[0][0][1][0]=0.0;	dSt[0][0][1][1]=0.0;	dSt[0][0][1][2]=0.0;
								dSt[0][1][0][0]=0.0;	dSt[0][1][0][1]=0.0;	dSt[0][1][0][2]=0.0;
								dSt[0][1][1][0]=Nb_x;	dSt[0][1][1][1]=Nb_y;	dSt[0][1][1][2]=0.0;
								dSt[1][0][0][0]=0.0;	dSt[1][0][0][1]=0.0;	dSt[1][0][0][2]=0.0;
								dSt[1][0][1][0]=0.0;	dSt[1][0][1][1]=0.0;	dSt[1][0][1][2]=0.0;
								dSt[1][1][0][0]=0.0;	dSt[1][1][0][1]=0.0;	dSt[1][1][0][2]=0.0;
								dSt[1][1][1][0]=0.0;	dSt[1][1][1][1]=0.0;	dSt[1][1][1][2]=0.0;
							}
							else if(j==4)
							{
								St[0][0][0]=0.0;	St[0][0][1]=0.0;		St[0][1][0]=0.0;	St[0][1][1]=0.0;	St[1][0][0]=Nb;		St[1][0][1]=0.0;	St[1][1][0]=0.0;	St[1][1][1]=0.0;
								dSt[0][0][0][0]=0.0;	dSt[0][0][0][1]=0.0;	dSt[0][0][0][2]=0.0;
								dSt[0][0][1][0]=0.0;	dSt[0][0][1][1]=0.0;	dSt[0][0][1][2]=0.0;
								dSt[0][1][0][0]=0.0;	dSt[0][1][0][1]=0.0;	dSt[0][1][0][2]=0.0;
								dSt[0][1][1][0]=0.0;	dSt[0][1][1][1]=0.0;	dSt[0][1][1][2]=0.0;
								dSt[1][0][0][0]=Nb_x;	dSt[1][0][0][1]=Nb_y;	dSt[1][0][0][2]=0.0;
								dSt[1][0][1][0]=0.0;	dSt[1][0][1][1]=0.0;	dSt[1][0][1][2]=0.0;
								dSt[1][1][0][0]=0.0;	dSt[1][1][0][1]=0.0;	dSt[1][1][0][2]=0.0;
								dSt[1][1][1][0]=0.0;	dSt[1][1][1][1]=0.0;	dSt[1][1][1][2]=0.0;
							}
							else if(j==5)
							{
								St[0][0][0]=0.0;	St[0][0][1]=0.0;		St[0][1][0]=0.0;	St[0][1][1]=0.0;	St[1][0][0]=0.0;	St[1][0][1]=Nb;		St[1][1][0]=0.0;	St[1][1][1]=0.0;
								dSt[0][0][0][0]=0.0;	dSt[0][0][0][1]=0.0;	dSt[0][0][0][2]=0.0;
								dSt[0][0][1][0]=0.0;	dSt[0][0][1][1]=0.0;	dSt[0][0][1][2]=0.0;
								dSt[0][1][0][0]=0.0;	dSt[0][1][0][1]=0.0;	dSt[0][1][0][2]=0.0;
								dSt[0][1][1][0]=0.0;	dSt[0][1][1][1]=0.0;	dSt[0][1][1][2]=0.0;
								dSt[1][0][0][0]=0.0;	dSt[1][0][0][1]=0.0;	dSt[1][0][0][2]=0.0;
								dSt[1][0][1][0]=Nb_x;	dSt[1][0][1][1]=Nb_y;	dSt[1][0][1][2]=0.0;
								dSt[1][1][0][0]=0.0;	dSt[1][1][0][1]=0.0;	dSt[1][1][0][2]=0.0;
								dSt[1][1][1][0]=0.0;	dSt[1][1][1][1]=0.0;	dSt[1][1][1][2]=0.0;
							}
							else if(j==6)
							{
								St[0][0][0]=0.0;	St[0][0][1]=0.0;		St[0][1][0]=0.0;	St[0][1][1]=0.0;	St[1][0][0]=0.0;	St[1][0][1]=0.0;	St[1][1][0]=Nb;		St[1][1][1]=0.0;
								dSt[0][0][0][0]=0.0;	dSt[0][0][0][1]=0.0;	dSt[0][0][0][2]=0.0;
								dSt[0][0][1][0]=0.0;	dSt[0][0][1][1]=0.0;	dSt[0][0][1][2]=0.0;
								dSt[0][1][0][0]=0.0;	dSt[0][1][0][1]=0.0;	dSt[0][1][0][2]=0.0;
								dSt[0][1][1][0]=0.0;	dSt[0][1][1][1]=0.0;	dSt[0][1][1][2]=0.0;
								dSt[1][0][0][0]=0.0;	dSt[1][0][0][1]=0.0;	dSt[1][0][0][2]=0.0;
								dSt[1][0][1][0]=0.0;	dSt[1][0][1][1]=0.0;	dSt[1][0][1][2]=0.0;
								dSt[1][1][0][0]=Nb_x;	dSt[1][1][0][1]=Nb_y;	dSt[1][1][0][2]=0.0;
								dSt[1][1][1][0]=0.0;	dSt[1][1][1][1]=0.0;	dSt[1][1][1][2]=0.0;
							}
							else if(j==7)
							{
								St[0][0][0]=0.0;	St[0][0][1]=0.0;		St[0][1][0]=0.0;	St[0][1][1]=0.0;	St[1][0][0]=0.0;	St[1][0][1]=0.0;	St[1][1][0]=0.0;	St[1][1][1]=Nb;
								dSt[0][0][0][0]=0.0;	dSt[0][0][0][1]=0.0;	dSt[0][0][0][2]=0.0;
								dSt[0][0][1][0]=0.0;	dSt[0][0][1][1]=0.0;	dSt[0][0][1][2]=0.0;
								dSt[0][1][0][0]=0.0;	dSt[0][1][0][1]=0.0;	dSt[0][1][0][2]=0.0;
								dSt[0][1][1][0]=0.0;	dSt[0][1][1][1]=0.0;	dSt[0][1][1][2]=0.0;
								dSt[1][0][0][0]=0.0;	dSt[1][0][0][1]=0.0;	dSt[1][0][0][2]=0.0;
								dSt[1][0][1][0]=0.0;	dSt[1][0][1][1]=0.0;	dSt[1][0][1][2]=0.0;
								dSt[1][1][0][0]=0.0;	dSt[1][1][0][1]=0.0;	dSt[1][1][0][2]=0.0;
								dSt[1][1][1][0]=Nb_x;	dSt[1][1][1][1]=Nb_y;	dSt[1][1][1][2]=0.0;
							}

							KSdot[a][i][b][j]=0.0;
							//Galerkin part for K
							for (k=0; k<3; k++)
							{
								for (l=0; l<3; l++)
								{
									for (u=0; u<3; u++)
									{
										KSdot[a][i][b][j]+=St[k][l][u]*v[k][l][u];		//This is for backwards Euler
										//KSdot[a][i][b][j]+=St[k][l][u]*v[k][l][u];			//This is for trapezoidal rule
									}
								}
							}
							for (k=0; k<3; k++)
							{
								for (l=0; l<3; l++)
								{
									for (u=0; u<3; u++)
									{
										for (w=0; w<3; w++)
										{
											KSdot[a][i][b][j]+=dt*St[k][l][w]*fullVs[w]*dv[k][l][u][u];		//This is for backwards Euler
											//KSdot[a][i][b][j]+=0.5*dt*(St[k][l][w]*fullVs[w])*dv[k][l][u][u];		//This is for trapezoidal rule

										}
									}
								}
							}
							//Least squares part for K
							for (k=0; k<3; k++)
							{
								for (l=0; l<3; l++)
								{
									for (u=0; u<3; u++)
									{
										KSdot[a][i][b][j]+=St[k][l][u]*v[k][l][u];		//This is for backwards Euler
										//KSdot[a][i][b][j]+=St[k][l][u]*v[k][l][u];			//This is for trapezoidal rule
									}
								}
							}
							for (k=0; k<3; k++)
							{
								for (l=0; l<3; l++)
								{
									for (u=0; u<3; u++)
									{
										for (w=0; w<3; w++)
										{
											//Backwards Euler
											KSdot[a][i][b][j]+=-dt*St[k][l][u]*(dv[k][l][w][u]*fullVs[w]+v[k][l][w]*full_dVs[w][u])
															   -dt*(dSt[k][l][w][u]*fullVs[w]+St[k][l][w]*full_dVs[w][u])*v[k][l][u]
															   +dt*dt*(dSt[k][l][w][u]*fullVs[w]+St[k][l][w]*full_dVs[w][u])*(dv[k][l][w][u]*fullVs[w]+v[k][l][w]*full_dVs[w][u]);			//This is for backwards Euler

											//Trapezoid Rule
											//KSdot[a][i][b][j]+=-0.5*dt*St[k][l][u]*(dv[k][l][w][u]*fullVs[w]+v[k][l][w]*full_dVs[w][u])
											//                  -0.5*dt*(dSt[k][l][w][u]*fullVs[w]+St[k][l][w]*full_dVs[w][u])*v[k][l][u]
											//                   +0.25*dt*dt*(dSt[k][l][w][u]*fullVs[w]+St[k][l][w]*full_dVs[w][u])*(dv[k][l][w][u]*fullVs[w]+v[k][l][w]*full_dVs[w][u]);			//This is for trpezoidal rule
										}
									}
								}
							}
						}
					}

					FSdot[a][i] = 0.0;
					//Galerkin part of F
					for (k=0; k<3; k++)
					{
						for (l=0; l<3; l++)
						{
							for (u=0; u<3; u++)
							{
								//This is backwards Euler
								FSdot[a][i]+=fullS[k][l][u]*v[k][l][u];		//Term with (Pi x VPi) should be included here in the future
								//This is trapezoidal rule
								//FSdot[a][i]+=fullS[k][l][u]*v[k][l][u];		//Term with (Pi x VPi) should be included here in the future (these are equal, consider deleting)
							}
						}	
					}
					for (k=0; k<3; k++)
					{
						for (l=0; l<3; l++)
						{
							for (u=0; u<3; u++)
							{
								for (w=0;w<3;w++)
								{
									//This is for trapezoidal rule, there is no equivalent term in backwards Euler
									//FSdot[a][i]+=0.5*dt*(fulld_S[k][l][w][u]*fullVs[w]+fullS[k][l][w]*full_dVs[w][u])*v[k][l][u];		//Term with (Pi x VPi) should be included here in the future (these are equal, consider deleting)
								}
							}
						}	
					}
					//Least squares part of F
					for (k=0; k<3; k++)
					{
						for (l=0; l<3; l++)
						{
							for (u=0; u<3; u++)
							{
								//This is for backwards Euler
								FSdot[a][i]+=fullS[k][l][u]*v[k][l][u];		//Term with (Pi x VPi) should be included here too in the future
								//This is for trapezoidal rule
								//FSdot[a][i]+=fullS[k][l][u]*v[k][l][u];		//Term with (Pi x VPi) should be included here too in the future (these terms are equal, consider deleting)
							}
						}	
					}
					for (k=0; k<3; k++)
					{
						for (l=0; l<3; l++)
						{
							for (u=0; u<3; u++)
							{
								for (w=0; w<3; w++)
								{
									//This term is for backwards Euler
									FSdot[a][i]+=-dt*fullS[k][l][u]*(dv[k][l][w][u]*fullVs[w]+v[k][l][w]*full_dVs[w][u]);		//Term with dt^2*(Pi x VPi) should be included here in the future
									//This term is for trapezoidal rule
									//FSdot[a][i]+=0.5*dt*(fulld_S[k][l][w][u]*fullVs[w]+fullS[k][l][w]*full_dVs[w][u])*v[k][l][u]
									//            -0.5*dt*fullS[k][l][u]*(dv[k][l][w][u]*fullVs[w]+v[k][l][w]*full_dVs[w][u])
									//            -0.25*dt*dt*(fulld_S[k][l][w][u]*fullVs[w]+fullS[k][l][w]*full_dVs[w][u])*(dv[k][l][w][u]*fullVs[w]+v[k][l][w]*full_dVs[w][u]);		//Term with dt^2*(Pi x VPi) should be included here in the future
								}
							}
						}	
					}
				}
			}
		}
		return 0;	
	}

	//This function fills just the right hand side, because the stiffness matrix does not change, at least in small deformation setting,
	//this speeds up the computation tremendously.
	#undef  __FUNCT__
	#define __FUNCT__ "SdotFuncF"
				 //SdotFuncF(point_Sdot,point_S,Fpoint_Sdot,Int1a,Int2a,S_norm,S0Sdot,&user);CHKERRQ(ierr);
	PetscErrorCode SdotFuncF(IGAPoint p,IGAPoint pS,PetscReal *F,PetscReal Int1a,PetscReal Int2a,PetscReal normS,PetscReal *US, void *ctx)
	{
		PetscInt a,i,k,l,u,w,dof=p->dof,nen=p->nen;	

		const PetscReal *N0,(*N1)[2];
		IGAPointGetShapeFuns(p,0,(const PetscReal**)&N0);									//Value of the shape functions
		IGAPointGetShapeFuns(p,1,(const PetscReal**)&N1);									//Derivatives of the shape functions

		AppCtx    *user  = (AppCtx *)ctx;
		PetscReal dt     = user->dt;

		//PetscReal Vs[2],dVs[2][2];															//Seems that if I create this one with size 3, next line gives problems
		//IGAPointFormValue(pVs,UVs,&Vs[0]);													//This fills the values
		//IGAPointFormGrad (pVs,UVs,&dVs[0][0]);												//This fills the values

		PetscReal fullVs[3]={0};
		PetscReal full_dVs[3][3]={0};

		fullVs[0]=Int1a/normS; fullVs[1]=-Int2a/normS;
		fullVs[0]=0.0;		 fullVs[1]=-1.0;			//THIS IS FOR DEBUGGING!!
		full_dVs[0][0]=0.0; full_dVs[0][1]=0.0;
		full_dVs[1][0]=0.0; full_dVs[1][1]=0.0;

		//PetscReal S[8];																			//Create array to receive S and grad(S) (backwards Euler)
		PetscReal S[8], dS[8][2];																//Create array to receive S and grad(S) (Trapezoidal rule)
		IGAPointFormValue(pS,US,&S[0]);															//This fills the values of S (remember that S has 8 non zero components in 2D)
		IGAPointFormGrad (pS,US,&dS[0][0]);													//This fills the values of the gradient

		PetscReal fullS[3][3][3]={0};
		fullS[0][0][0]=S[0]; fullS[0][0][1]=S[1];												//Expand S to full 3rd order form, only non-zero elements
		fullS[0][1][0]=S[2]; fullS[0][1][1]=S[3];
		fullS[1][0][0]=S[4]; fullS[1][0][1]=S[5];
		fullS[1][1][0]=S[6]; fullS[1][1][1]=S[7];

		PetscReal fulld_S[3][3][3][3]={0};													//Expand grad(S) to full tensor order form, only non-zero elements
		fulld_S[0][0][0][0]=dS[0][0]; fulld_S[0][0][0][1]=dS[0][1]; 
		fulld_S[0][0][1][0]=dS[1][0]; fulld_S[0][0][1][1]=dS[1][1]; 
		fulld_S[0][1][0][0]=dS[2][0]; fulld_S[0][1][0][1]=dS[2][1]; 
		fulld_S[0][1][1][0]=dS[3][0]; fulld_S[0][1][1][1]=dS[3][1];
		fulld_S[1][0][0][0]=dS[4][0]; fulld_S[1][0][0][1]=dS[4][1]; 
		fulld_S[1][0][1][0]=dS[5][0]; fulld_S[1][0][1][1]=dS[5][1]; 
		fulld_S[1][1][0][0]=dS[6][0]; fulld_S[1][1][0][1]=dS[6][1]; 
		fulld_S[1][1][1][0]=dS[7][0]; fulld_S[1][1][1][1]=dS[7][1];

		PetscReal v[3][3][3]={0};
		PetscReal dv[3][3][3][3]={0};

		PetscReal (*FSdot)[dof] = (PetscReal (*)[dof])F;
		if (p->atboundary)
		{
			return 0;
		}
		else
		{
			for(a=0;a<nen;a++)
			{
				PetscReal Na=N0[a];
				PetscReal Na_x=N1[a][0];		PetscReal Na_y=N1[a][1];

				for (i=0; i<dof; i++)
				{
					if (i==0)
						{
							v[0][0][0]=Na;		v[0][0][1]=0.0;		v[0][1][0]=0.0;	v[0][1][1]=0.0;	v[1][0][0]=0.0;	v[1][0][1]=0.0;	v[1][1][0]=0.0;	v[1][1][1]=0.0;
							dv[0][0][0][0]=Na_x;	dv[0][0][0][1]=Na_y;	dv[0][0][0][2]=0.0;
							dv[0][0][1][0]=0.0;		dv[0][0][1][1]=0.0;		dv[0][0][1][2]=0.0;
							dv[0][1][0][0]=0.0;		dv[0][1][0][1]=0.0;		dv[0][1][0][2]=0.0;
							dv[0][1][1][0]=0.0;		dv[0][1][1][1]=0.0;		dv[0][1][1][2]=0.0;
							dv[1][0][0][0]=0.0;		dv[1][0][0][1]=0.0;		dv[1][0][0][2]=0.0;
							dv[1][0][1][0]=0.0;		dv[1][0][1][1]=0.0;		dv[1][0][1][2]=0.0;
							dv[1][1][0][0]=0.0;		dv[1][1][0][1]=0.0;		dv[1][1][0][2]=0.0;
							dv[1][1][1][0]=0.0;		dv[1][1][1][1]=0.0;		dv[1][1][1][2]=0.0;
						}
						else if(i==1)
						{
							v[0][0][0]=0.0;		v[0][0][1]=Na;		v[0][1][0]=0.0;	v[0][1][1]=0.0;	v[1][0][0]=0.0;	v[1][0][1]=0.0;	v[1][1][0]=0.0;	v[1][1][1]=0.0;
							dv[0][0][0][0]=0.0;		dv[0][0][0][1]=0.0;		dv[0][0][0][2]=0.0;
							dv[0][0][1][0]=Na_x;	dv[0][0][1][1]=Na_y;	dv[0][0][1][2]=0.0;
							dv[0][1][0][0]=0.0;		dv[0][1][0][1]=0.0;		dv[0][1][0][2]=0.0;
							dv[0][1][1][0]=0.0;		dv[0][1][1][1]=0.0;		dv[0][1][1][2]=0.0;
							dv[1][0][0][0]=0.0;		dv[1][0][0][1]=0.0;		dv[1][0][0][2]=0.0;
							dv[1][0][1][0]=0.0;		dv[1][0][1][1]=0.0;		dv[1][0][1][2]=0.0;
							dv[1][1][0][0]=0.0;		dv[1][1][0][1]=0.0;		dv[1][1][0][2]=0.0;
							dv[1][1][1][0]=0.0;		dv[1][1][1][1]=0.0;		dv[1][1][1][2]=0.0;
						}
						else if(i==2)
						{
							v[0][0][0]=0.0;		v[0][0][1]=0.0;		v[0][1][0]=Na;	v[0][1][1]=0.0;	v[1][0][0]=0.0;	v[1][0][1]=0.0;	v[1][1][0]=0.0;	v[1][1][1]=0.0;
							dv[0][0][0][0]=0.0;		dv[0][0][0][1]=0.0;		dv[0][0][0][2]=0.0;
							dv[0][0][1][0]=0.0;		dv[0][0][1][1]=0.0;		dv[0][0][1][2]=0.0;
							dv[0][1][0][0]=Na_x;	dv[0][1][0][1]=Na_x;	dv[0][1][0][2]=0.0;
							dv[0][1][1][0]=0.0;		dv[0][1][1][1]=0.0;		dv[0][1][1][2]=0.0;
							dv[1][0][0][0]=0.0;		dv[1][0][0][1]=0.0;		dv[1][0][0][2]=0.0;
							dv[1][0][1][0]=0.0;		dv[1][0][1][1]=0.0;		dv[1][0][1][2]=0.0;
							dv[1][1][0][0]=0.0;		dv[1][1][0][1]=0.0;		dv[1][1][0][2]=0.0;
							dv[1][1][1][0]=0.0;		dv[1][1][1][1]=0.0;		dv[1][1][1][2]=0.0;
						}
						else if(i==3)
						{
							v[0][0][0]=0.0;		v[0][0][1]=0.0;		v[0][1][0]=0.0;	v[0][1][1]=Na;	v[1][0][0]=0.0;	v[1][0][1]=0.0;	v[1][1][0]=0.0;	v[1][1][1]=0.0;
							dv[0][0][0][0]=0.0;		dv[0][0][0][1]=0.0;		dv[0][0][0][2]=0.0;
							dv[0][0][1][0]=0.0;		dv[0][0][1][1]=0.0;		dv[0][0][1][2]=0.0;
							dv[0][1][0][0]=0.0;		dv[0][1][0][1]=0.0;		dv[0][1][0][2]=0.0;
							dv[0][1][1][0]=Na_x;	dv[0][1][1][1]=Na_y;	dv[0][1][1][2]=0.0;
							dv[1][0][0][0]=0.0;		dv[1][0][0][1]=0.0;		dv[1][0][0][2]=0.0;
							dv[1][0][1][0]=0.0;		dv[1][0][1][1]=0.0;		dv[1][0][1][2]=0.0;
							dv[1][1][0][0]=0.0;		dv[1][1][0][1]=0.0;		dv[1][1][0][2]=0.0;
							dv[1][1][1][0]=0.0;		dv[1][1][1][1]=0.0;		dv[1][1][1][2]=0.0;
						}
						else if(i==4)
						{
							v[0][0][0]=0.0;		v[0][0][1]=0.0;		v[0][1][0]=0.0;	v[0][1][1]=0.0;	v[1][0][0]=Na;	v[1][0][1]=0.0;	v[1][1][0]=0.0;	v[1][1][1]=0.0;
							dv[0][0][0][0]=0.0;		dv[0][0][0][1]=0.0;		dv[0][0][0][2]=0.0;
							dv[0][0][1][0]=0.0;		dv[0][0][1][1]=0.0;		dv[0][0][1][2]=0.0;
							dv[0][1][0][0]=0.0;		dv[0][1][0][1]=0.0;		dv[0][1][0][2]=0.0;
							dv[0][1][1][0]=0.0;		dv[0][1][1][1]=0.0;		dv[0][1][1][2]=0.0;
							dv[1][0][0][0]=Na_x;	dv[1][0][0][1]=Na_y;	dv[1][0][0][2]=0.0;
							dv[1][0][1][0]=0.0;		dv[1][0][1][1]=0.0;		dv[1][0][1][2]=0.0;
							dv[1][1][0][0]=0.0;		dv[1][1][0][1]=0.0;		dv[1][1][0][2]=0.0;
							dv[1][1][1][0]=0.0;		dv[1][1][1][1]=0.0;		dv[1][1][1][2]=0.0;
						}
						else if(i==5)
						{
							v[0][0][0]=0.0;		v[0][0][1]=0.0;		v[0][1][0]=0.0;	v[0][1][1]=0.0;	v[1][0][0]=0.0;	v[1][0][1]=Na;	v[1][1][0]=0.0;	v[1][1][1]=0.0;
							dv[0][0][0][0]=0.0;		dv[0][0][0][1]=0.0;		dv[0][0][0][2]=0.0;
							dv[0][0][1][0]=0.0;		dv[0][0][1][1]=0.0;		dv[0][0][1][2]=0.0;
							dv[0][1][0][0]=0.0;		dv[0][1][0][1]=0.0;		dv[0][1][0][2]=0.0;
							dv[0][1][1][0]=0.0;		dv[0][1][1][1]=0.0;		dv[0][1][1][2]=0.0;
							dv[1][0][0][0]=0.0;		dv[1][0][0][1]=0.0;		dv[1][0][0][2]=0.0;
							dv[1][0][1][0]=Na_x;	dv[1][0][1][1]=Na_y;	dv[1][0][1][2]=0.0;
							dv[1][1][0][0]=0.0;		dv[1][1][0][1]=0.0;		dv[1][1][0][2]=0.0;
							dv[1][1][1][0]=0.0;		dv[1][1][1][1]=0.0;		dv[1][1][1][2]=0.0;
						}
						else if(i==6)
						{
							v[0][0][0]=0.0;		v[0][0][1]=0.0;		v[0][1][0]=0.0;	v[0][1][1]=0.0;	v[1][0][0]=0.0;	v[1][0][1]=0.0;	v[1][1][0]=Na;	v[1][1][1]=0.0;
							dv[0][0][0][0]=0.0;		dv[0][0][0][1]=0.0;		dv[0][0][0][2]=0.0;
							dv[0][0][1][0]=0.0;		dv[0][0][1][1]=0.0;		dv[0][0][1][2]=0.0;
							dv[0][1][0][0]=0.0;		dv[0][1][0][1]=0.0;		dv[0][1][0][2]=0.0;
							dv[0][1][1][0]=0.0;		dv[0][1][1][1]=0.0;		dv[0][1][1][2]=0.0;
							dv[1][0][0][0]=0.0;		dv[1][0][0][1]=0.0;		dv[1][0][0][2]=0.0;
							dv[1][0][1][0]=0.0;		dv[1][0][1][1]=0.0;		dv[1][0][1][2]=0.0;
							dv[1][1][0][0]=Na_x;	dv[1][1][0][1]=Na_y;	dv[1][1][0][2]=0.0;
							dv[1][1][1][0]=0.0;		dv[1][1][1][1]=0.0;		dv[1][1][1][2]=0.0;
						}
						else if(i==7)
						{
							v[0][0][0]=0.0;		v[0][0][1]=0.0;		v[0][1][0]=0.0;	v[0][1][1]=0.0;	v[1][0][0]=0.0;	v[1][0][1]=0.0;	v[1][1][0]=0.0;	v[1][1][1]=Na;
							dv[0][0][0][0]=0.0;		dv[0][0][0][1]=0.0;		dv[0][0][0][2]=0.0;
							dv[0][0][1][0]=0.0;		dv[0][0][1][1]=0.0;		dv[0][0][1][2]=0.0;
							dv[0][1][0][0]=0.0;		dv[0][1][0][1]=0.0;		dv[0][1][0][2]=0.0;
							dv[0][1][1][0]=0.0;		dv[0][1][1][1]=0.0;		dv[0][1][1][2]=0.0;
							dv[1][0][0][0]=0.0;		dv[1][0][0][1]=0.0;		dv[1][0][0][2]=0.0;
							dv[1][0][1][0]=0.0;		dv[1][0][1][1]=0.0;		dv[1][0][1][2]=0.0;
							dv[1][1][0][0]=0.0;		dv[1][1][0][1]=0.0;		dv[1][1][0][2]=0.0;
							dv[1][1][1][0]=Na_x;	dv[1][1][1][1]=Na_y;	dv[1][1][1][2]=0.0;
						}

					FSdot[a][i] = 0.0;
					//Galerkin part of F
					for (k=0; k<3; k++)
					{
						for (l=0; l<3; l++)
						{
							for (u=0; u<3; u++)
							{
								//This is backwards Euler
								FSdot[a][i]+=fullS[k][l][u]*v[k][l][u];		//Term with (Pi x VPi) should be included here in the future
								//This is trapezoidal rule
								//FSdot[a][i]+=fullS[k][l][u]*v[k][l][u];		//Term with (Pi x VPi) should be included here in the future (these are equal, consider deleting)
							}
						}	
					}
					for (k=0; k<3; k++)
					{
						for (l=0; l<3; l++)
						{
							for (u=0; u<3; u++)
							{
								for (w=0;w<3;w++)
								{
									//This is for trapezoidal rule, there is no equivalent term in backwards Euler
									//FSdot[a][i]+=0.5*dt*(fulld_S[k][l][w][u]*fullVs[w]+fullS[k][l][w]*full_dVs[w][u])*v[k][l][u];		//Term with (Pi x VPi) should be included here in the future (these are equal, consider deleting)
								}
							}
						}	
					}
					//Least squares part of F
					for (k=0; k<3; k++)
					{
						for (l=0; l<3; l++)
						{
							for (u=0; u<3; u++)
							{
								//This is for backwards Euler
								FSdot[a][i]+=fullS[k][l][u]*v[k][l][u];		//Term with (Pi x VPi) should be included here too in the future
								//This is for trapezoidal rule
								//FSdot[a][i]+=fullS[k][l][u]*v[k][l][u];		//Term with (Pi x VPi) should be included here too in the future (these terms are equal, consider deleting)
							}
						}	
					}
					for (k=0; k<3; k++)
					{
						for (l=0; l<3; l++)
						{
							for (u=0; u<3; u++)
							{
								for (w=0; w<3; w++)
								{
									//This term is for backwards Euler
									FSdot[a][i]+=-dt*fullS[k][l][u]*(dv[k][l][w][u]*fullVs[w]+v[k][l][w]*full_dVs[w][u]);		//Term with dt^2*(Pi x VPi) should be included here in the future
									//This term is for trapezoidal rule
									//FSdot[a][i]+=0.5*dt*(fulld_S[k][l][w][u]*fullVs[w]+fullS[k][l][w]*full_dVs[w][u])*v[k][l][u]
									//            -0.5*dt*fullS[k][l][u]*(dv[k][l][w][u]*fullVs[w]+v[k][l][w]*full_dVs[w][u])
									//            -0.25*dt*dt*(fulld_S[k][l][w][u]*fullVs[w]+fullS[k][l][w]*full_dVs[w][u])*(dv[k][l][w][u]*fullVs[w]+v[k][l][w]*full_dVs[w][u]);		//Term with dt^2*(Pi x VPi) should be included here in the future
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

//System for norm of inital S
	#undef  __FUNCT__
	#define __FUNCT__ "Snorm"
	PetscErrorCode Snorm(IGAPoint pS,PetscReal *SNorm,PetscReal *Norm,PetscReal *US,void *ctx)
	{
		PetscInt dof=pS->dof,nen=pS->nen;
		PetscInt a,i;

		PetscReal S[8];																		//Create array to receive Alfa
		IGAPointFormValue(pS,US,&S[0]);														//This fills the values

		PetscReal xi=0.0;
		for (i=0;i<8;i++)
		{
			if (S[i]>0.05 || S[i]<-0.05)
			{
				xi=1.0;
			}	
		}

		PetscReal (*FS)[dof] = (PetscReal (*)[dof])SNorm;
		PetscReal (*FI)[dof] = (PetscReal (*)[dof])Norm;			//This vector will have xi(S), so we can integrate it.

		if (pS->atboundary)
		{
			for(a=0;a<nen;a++)
			{
				for(i=0;i<dof;i++)
				{
					FS[a][i]=0.0;
				}
			}
			return 0;
		}
		else
		{
			for(a=0;a<nen;a++)
			{
				for(i=0;i<dof;i++)
				{
					FS[a][i]=0.0;
				}
			}
			for(i=0;i<dof;i++)
			{
				FS[0][i]=S[i];
			}
			
			FI[0][0]=xi;

		}
		return 0;
	}
//

//System for \dot{z} (debug correct this with new Vs!!!)
	#undef  __FUNCT__
	#define __FUNCT__ "ZdotSystem"
				 //ZdotSystem(pointZdot,point_Al_hat,pointVa,point_S,KpointZdot,FpointZdot,AlZdot,VaZdot,SZdot,Int1a,Int2a,S_norm,NULL);CHKERRQ(ierr);
	PetscErrorCode ZdotSystem(IGAPoint p,IGAPoint pAl,IGAPoint pVa,IGAPoint pS,PetscReal *K,PetscReal *F,PetscReal *UAl,PetscReal *UVa,PetscReal *US,PetscReal Vs1,PetscReal Vs2,PetscReal normS, void *ctx)
	{
		const PetscReal *N0,(*N1)[2];
		IGAPointGetShapeFuns(p,0,(const PetscReal**)&N0);
		IGAPointGetShapeFuns(p,1,(const PetscReal**)&N1);

		PetscInt a,b,i,j,u,w,k,l,nen=p->nen, dof=p->dof;

		//Definition of alternating tensor
		const PetscReal e[3][3][3]=
		{
			{{0.0,0.0,0.0},{0.0,0.0,1.0},{0.0,-1.0,0.0}},
			{{0.0,0.0,-1.0},{0.0,0.0,0.0},{1.0,0.0,0.0}},
			{{0.0,1.0,0.0},{-1.0,0.0,0.0},{0.0,0.0,0.0}}
		};

		PetscReal alfa[2];																		//Create array to receive Alfa																
		IGAPointFormValue(pAl,UAl,&alfa[0]);													//This fills the values

		PetscReal fullAlfa[3][3]={0};
		fullAlfa[0][2]=alfa[0]; fullAlfa[1][2]=alfa[1];
	
		PetscReal Va[2];																	//Seems that if I create this one with size 3, next line gives problems
		IGAPointFormValue(pVa,UVa,&Va[0]);													//This fills the values

		PetscReal fullVa[3]={0};	
		fullVa[0]=Va[0]; fullVa[1]=Va[1];

		//PetscReal Vs[2],dVs[2][2];															//Seems that if I create this one with size 3, next line gives problems
		//IGAPointFormValue(pVs,UVs,&Vs[0]);													//This fills the values

		PetscReal fullVs[3]={0};
		fullVs[0]=Vs1/normS; fullVs[1]=-Vs2/normS;
		fullVs[0]=0.0;		 fullVs[1]=-1.0;			//THIS IS FOR DEBUGGING!!

		PetscReal S[8];																		//Create array to receive S
		IGAPointFormValue(pS,US,&S[0]);														//This fills the values of S (remember that S has 8 non zero components in 2D)

		PetscReal fullS[3][3][3]={0};
		fullS[0][0][0]=S[0]; fullS[0][0][1]=S[1];											//Expand S to full 3rd order form, only non-zero elements
		fullS[0][1][0]=S[2]; fullS[0][1][1]=S[3];
		fullS[1][0][0]=S[4]; fullS[1][0][1]=S[5];
		fullS[1][1][0]=S[6]; fullS[1][1][1]=S[7];

		for(i=0;i<3;i++)
		{
			for(j=0;j<3;j++)
			{
				for(k=0;k<3;k++)
				{
					for(l=0;l<3;l++)
					{
						fullAlfa[i][j]=fullAlfa[i][j]-fullS[i][k][l]*e[j][k][l];
					}
				}
			}
		}

		PetscReal dv[3][3]={0};
		PetscReal dz[3][3]={0};

		PetscReal (*Keq)[dof][nen][dof] = (typeof(Keq)) K;
		PetscReal (*Feq)[dof] = (PetscReal (*)[dof]) F;

		if (p->atboundary)
		{
			return 0;
		}
		else
		{
			for (a=0; a<nen; a++) 
			{
				PetscReal Na_x = N1[a][0];			PetscReal Na_y = N1[a][1];

				for (b=0; b<nen; b++) 
				{
					PetscReal Nb_x = N1[b][0];		PetscReal Nb_y = N1[b][1];

					for (i=0;i<dof;i++)
					{
						if (i==0)
						{
							dv[0][0]=Na_x; dv[0][1]=Na_y; dv[0][2]=0.0;
							dv[1][0]=0.0;  dv[1][1]=0.0;  dv[1][2]=0.0;
							dv[2][0]=0.0;  dv[2][1]=0.0;  dv[2][2]=0.0;
						}
						else if(i==1)
						{
							dv[0][0]=0.0;  dv[0][1]=0.0;  dv[0][2]=0.0;
							dv[1][0]=Na_x; dv[1][1]=Na_y; dv[1][2]=0.0;
							dv[2][0]=0.0;  dv[2][1]=0.0;  dv[2][2]=0.0;
						}

						for (j=0;j<dof;j++)
						{
							if (j==0)
							{
								dz[0][0]=Nb_x; dz[0][1]=Nb_y; dz[0][2]=0.0;
								dz[1][0]=0.0;  dz[1][1]=0.0;  dz[1][2]=0.0;
								dz[2][0]=0.0;  dz[2][1]=0.0;  dz[2][2]=0.0;
							}
							else if(j==1)
							{
								dz[0][0]=0.0;  dz[0][1]=0.0;  dz[0][2]=0.0;
								dz[1][0]=Nb_x; dz[1][1]=Nb_y; dz[1][2]=0.0;
								dz[2][0]=0.0;  dz[2][1]=0.0;  dz[2][2]=0.0;
							}

							Keq[a][i][b][j]=0.0;
							for (k=0; k<3; k++)
							{
								for (l=0; l<3; l++)
								{
									Keq[a][i][b][j]+=dz[k][l]*dv[k][l];
								}
							}
						}
					}
				}
			}

			for(a=0;a<nen;a++)
			{
				PetscReal Na_x=N1[a][0];		PetscReal Na_y=N1[a][1];

				for (i=0;i<dof;i++)
				{
					if (i==0)
					{
						dv[0][0]=Na_x; dv[0][1]=Na_y; dv[0][2]=0.0;
						dv[1][0]=0.0;  dv[1][1]=0.0;  dv[1][2]=0.0;
						dv[2][0]=0.0;  dv[2][1]=0.0;  dv[2][2]=0.0;
					}
					else if (i==1)
					{
						dv[0][0]=0.0;  dv[0][1]=0.0;  dv[0][2]=0.0;
						dv[1][0]=Na_x; dv[1][1]=Na_y; dv[1][2]=0.0;
						dv[2][0]=0.0;  dv[2][1]=0.0;  dv[2][2]=0.0;
					}

					Feq[a][i] = 0.0;
					//(alpha X V^a)*grad(v)
					for (k=0;k<3;k++)
					{
						for (l=0;l<3;l++)
						{
							for (u=0;u<3;u++)
							{
								for (w=0;w<3;w++)
								{
									Feq[a][i]+=(e[l][u][w]*fullAlfa[k][u]*fullVa[w])*dv[k][l];
								}
							}
						}	
					}
					//(S*V^s)*grad(v)
					for (k=0;k<3;k++)
					{
						for (l=0;l<3;l++)
						{
							for (u=0;u<3;u++)
							{
								Feq[a][i]+=(fullS[k][l][u]*fullVs[u])*dv[k][l];
							}
						}	
					}
				}
			}
		}
		return 0;
	}
//

//System for L2 projection of grad(z)
	#undef  __FUNCT__
	#define __FUNCT__ "gradz"
	//PetscErrorCode Stress(IGAPoint p,IGAPoint pU, IGAPoint pHs,IGAPoint pChi,IGAPoint pZu,PetscReal *K,PetscReal *F,PetscReal *U,PetscReal *HS, PetscReal *Chi,PetscReal *Zu,void *ctx)		//When other fields like Pi or S (etc.) come into this, declare a IGAPoint pPi or pS and a new PescReal *UPi or *US for each
	PetscErrorCode gradz(IGAPoint p,IGAPoint pZu,PetscReal *K,PetscReal *F,PetscReal *Zu,void *ctx)		//When other fields like Pi or S (etc.) come into this, declare a IGAPoint pPi or pS and a new PetscReal *UPi or *US for each
	{
		//This functions generates the L2 projection matrix and the RHS
		const PetscReal *N0;
		IGAPointGetShapeFuns(p,0,(const PetscReal**)&N0);									//Value of the shape functions
		PetscInt a,b,i,j,k,nen=p->nen, dof=p->dof;

		PetscReal d_Z0[2][2];																//Container for values of the gradient
		IGAPointFormGrad (pZu,Zu,&d_Z0[0][0]);												//Fills the values of the gradient

		//Expanding z (and derivatives) to 3 components, more convenient for sums in for loops
		PetscReal fulld_z[3][3]={0};
		fulld_z[0][0]=d_Z0[0][0]; fulld_z[0][1]=d_Z0[0][1];
		fulld_z[1][0]=d_Z0[1][0]; fulld_z[1][1]=d_Z0[1][1];

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
						v[0][0]=N0[a];	v[0][1]=0.0;
						v[1][0]=0.0;	v[1][1]=0.0;

					}
					else if (i==1)
					{
						v[0][0]=0.0;	v[0][1]=N0[a];
						v[1][0]=0.0;	v[1][1]=0.0;
					}
					else if (i==2)
					{
						v[0][0]=0.0;	v[0][1]=0.0;
						v[1][0]=N0[a];	v[1][1]=0.0;
					}
					else if (i==3)
					{
						v[0][0]=0.0;	v[0][1]=0.0;
						v[1][0]=0.0;	v[1][1]=N0[a];
					}
					
					FCS[a][i]=0.0;
					for (j=0;j<3;j++)
					{
						for (k=0;k<3;k++)
						{
							FCS[a][i]+=fulld_z[j][k]*v[j][k];
						}
					}
				}
			}
		}
		return 0;
	}

	#undef  __FUNCT__
	#define __FUNCT__ "gradzF"
	//PetscErrorCode Stress(IGAPoint p,IGAPoint pU, IGAPoint pHs,IGAPoint pChi,IGAPoint pZu,PetscReal *K,PetscReal *F,PetscReal *U,PetscReal *HS, PetscReal *Chi,PetscReal *Zu,void *ctx)		//When other fields like Pi or S (etc.) come into this, declare a IGAPoint pPi or pS and a new PescReal *UPi or *US for each
	PetscErrorCode gradzF(IGAPoint p,IGAPoint pZu,PetscReal *F,PetscReal *Zu,void *ctx)		//When other fields like Pi or S (etc.) come into this, declare a IGAPoint pPi or pS and a new PetscReal *UPi or *US for each
	{
		//This functions generates the L2 projection matrix and the RHS
		const PetscReal *N0;
		IGAPointGetShapeFuns(p,0,(const PetscReal**)&N0);									//Value of the shape functions
		PetscInt a,i,j,k,nen=p->nen, dof=p->dof;

		PetscReal x[2];																		//Vector of reals, size equal to problem's dimension
		IGAPointFormGeomMap(p,x);															//Fills x with the coordinates of p, Gauss's point

		PetscReal d_Z0[2][2];																//Same for its gradient
		IGAPointFormGrad (pZu,Zu,&d_Z0[0][0]);												//Same for the gradient

		//Expanding z (and derivatives) to 3 components, more convenient for sums in for loops
		PetscReal fulld_z[3][3]={0};
		fulld_z[0][0]=d_Z0[0][0]; fulld_z[0][1]=d_Z0[0][1];
		fulld_z[1][0]=d_Z0[1][0]; fulld_z[1][1]=d_Z0[1][1];

		PetscReal (*FCS)[dof] = (PetscReal (*)[dof])F;

		PetscReal v[3][3]={0};

		if (p->atboundary)
		{
			return 0;
		}
		else
		{
			for(a=0 ;a<nen; a++)
			{
				for (i=0; i<dof; i++)
				{
					if (i==0)
					{
						v[0][0]=N0[a];	v[0][1]=0.0; 
						v[1][0]=0.0;	v[1][1]=0.0;

					}
					else if (i==1)
					{
						v[0][0]=0.0;	v[0][1]=N0[a]; 
						v[1][0]=0.0;	v[1][1]=0.0;
					}
					else if (i==2)
					{
						v[0][0]=0.0;	v[0][1]=0.0; 
						v[1][0]=N0[a];	v[1][1]=0.0;
					}
					else if (i==3)
					{
						v[0][0]=0.0;	v[0][1]=0.0; 
						v[1][0]=0.0;	v[1][1]=N0[a];
					}
					
					FCS[a][i]=0.0;

					for (j=0;j<3;j++)
					{
						for (k=0;k<3;k++)
						{
							FCS[a][i]+=fulld_z[j][k]*v[j][k];			
						}
					}
				}
			}
		}
		return 0;
	}
//

//System for L2 projection of Up
	#undef  __FUNCT__
	#define __FUNCT__ "SysUp"
				 //SysUp(pointUp,point_z,point_chi,KpointUp,FpointUp,Z0Up,ChiUp,NULL);CHKERRQ(ierr);
	PetscErrorCode SysUp(IGAPoint p,IGAPoint pZu,IGAPoint pChi,PetscReal *K,PetscReal *F,PetscReal *UZ,PetscReal *UChi,void *ctx)		//When other fields like Pi or S (etc.) come into this, declare a IGAPoint pPi or pS and a new PetscReal *UPi or *US for each
	{
		//This functions generates the L2 projection matrix and the RHS
		const PetscReal *N0;
		IGAPointGetShapeFuns(p,0,(const PetscReal**)&N0);									//Value of the shape functions
		PetscInt a,b,i,j,k,nen=p->nen, dof=p->dof;

		PetscReal d_Z0[2][2];																//Container for values of the gradient
		IGAPointFormGrad (pZu,UZ,&d_Z0[0][0]);												//Fills the values of the gradient

		//Expanding z (and derivatives) to 3 components, more convenient for sums in for loops
		PetscReal fulld_z[3][3]={0};
		fulld_z[0][0]=d_Z0[0][0]; fulld_z[0][1]=d_Z0[0][1];
		fulld_z[1][0]=d_Z0[1][0]; fulld_z[1][1]=d_Z0[1][1];

		//The four non-zero components of Chi are stored as a vector, restore them to an array with the correct indexing for value and derivative
		PetscReal Chi0[4];																	//Array to contain the vector chi(0)
		IGAPointFormValue(pChi,UChi,&Chi0[0]);												//Assign chi to its container

		PetscReal fullChi[3][3]={0};
		fullChi[0][0]=Chi0[0]; 	fullChi[0][1]=Chi0[1];
		fullChi[1][0]=Chi0[2]; 	fullChi[1][1]=Chi0[3];

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
						v[0][0]=N0[a];	v[0][1]=0.0;
						v[1][0]=0.0;	v[1][1]=0.0;

					}
					else if (i==1)
					{
						v[0][0]=0.0;	v[0][1]=N0[a];
						v[1][0]=0.0;	v[1][1]=0.0;
					}
					else if (i==2)
					{
						v[0][0]=0.0;	v[0][1]=0.0;
						v[1][0]=N0[a];	v[1][1]=0.0;
					}
					else if (i==3)
					{
						v[0][0]=0.0;	v[0][1]=0.0;
						v[1][0]=0.0;	v[1][1]=N0[a];
					}
					
					FCS[a][i]=0.0;
					for (j=0;j<3;j++)
					{
						for (k=0;k<3;k++)
						{
							FCS[a][i]+=(fulld_z[j][k]+fullChi[j][k])*v[j][k];
						}
					}
				}
			}
		}
		return 0;
	}

	#undef  __FUNCT__
	#define __FUNCT__ "SysUpF"
	//PetscErrorCode Stress(IGAPoint p,IGAPoint pU, IGAPoint pHs,IGAPoint pChi,IGAPoint pZu,PetscReal *K,PetscReal *F,PetscReal *U,PetscReal *HS, PetscReal *Chi,PetscReal *Zu,void *ctx)		//When other fields like Pi or S (etc.) come into this, declare a IGAPoint pPi or pS and a new PescReal *UPi or *US for each
	PetscErrorCode SysUpF(IGAPoint p,IGAPoint pZu,IGAPoint pChi,PetscReal *F,PetscReal *UZ,PetscReal *UChi,void *ctx)		//When other fields like Pi or S (etc.) come into this, declare a IGAPoint pPi or pS and a new PetscReal *UPi or *US for each
	{
		//This functions generates the L2 projection matrix and the RHS
		const PetscReal *N0;
		IGAPointGetShapeFuns(p,0,(const PetscReal**)&N0);									//Value of the shape functions
		PetscInt a,i,j,k,nen=p->nen, dof=p->dof;

		PetscReal d_Z0[2][2];																//Container for values of the gradient
		IGAPointFormGrad (pZu,UZ,&d_Z0[0][0]);												//Fills the values of the gradient

		//Expanding z (and derivatives) to 3 components, more convenient for sums in for loops
		PetscReal fulld_z[3][3]={0};
		fulld_z[0][0]=d_Z0[0][0]; fulld_z[0][1]=d_Z0[0][1];
		fulld_z[1][0]=d_Z0[1][0]; fulld_z[1][1]=d_Z0[1][1];

		//The four non-zero components of Chi are stored as a vector, restore them to an array with the correct indexing for value and derivative
		PetscReal Chi0[4];																	//Array to contain the vector chi(0)
		IGAPointFormValue(pChi,UChi,&Chi0[0]);												//Assign chi to its container

		PetscReal fullChi[3][3]={0};
		fullChi[0][0]=Chi0[0]; 	fullChi[0][1]=Chi0[1];
		fullChi[1][0]=Chi0[2]; 	fullChi[1][1]=Chi0[3];

		PetscReal (*FCS)[dof] = (PetscReal (*)[dof])F;

		PetscReal v[3][3]={0};

		if (p->atboundary)
		{
			return 0;
		}
		else
		{
			for(a=0 ;a<nen; a++)
			{
				for (i=0; i<dof; i++)
				{
					if (i==0)
					{
						v[0][0]=N0[a];	v[0][1]=0.0; 
						v[1][0]=0.0;	v[1][1]=0.0;

					}
					else if (i==1)
					{
						v[0][0]=0.0;	v[0][1]=N0[a]; 
						v[1][0]=0.0;	v[1][1]=0.0;
					}
					else if (i==2)
					{
						v[0][0]=0.0;	v[0][1]=0.0; 
						v[1][0]=N0[a];	v[1][1]=0.0;
					}
					else if (i==3)
					{
						v[0][0]=0.0;	v[0][1]=0.0; 
						v[1][0]=0.0;	v[1][1]=N0[a];
					}
					
					FCS[a][i]=0.0;

					for (j=0;j<3;j++)
					{
						for (k=0;k<3;k++)
						{
							FCS[a][i]+=(fulld_z[j][k]+fullChi[j][k])*v[j][k];			
						}
					}
				}
			}
		}
		return 0;
	}
//

//System for general L2 projection
	//Build up until it takes all fields, so they can be easily combined in this function.
	#undef  __FUNCT__
	#define __FUNCT__ "proj"
	PetscErrorCode proj(IGAPoint p,IGAPoint pAl,IGAPoint pS, IGAPoint pVa,IGAPoint pZu,IGAPoint pu,IGAPoint pChi,PetscReal *K,PetscReal *F,PetscReal *UAl,PetscReal *US,PetscReal *UVa,
						PetscReal *UZu,PetscReal *UChi,PetscReal *Uu,PetscReal Vs1,PetscReal Vs2,PetscReal normS,void *ctx)
	{
		//This functions generates the L2 projection matrix and the RHS
		const PetscReal *N0;
		IGAPointGetShapeFuns(p,0,(const PetscReal**)&N0);								//Value of the shape functions
		PetscInt a,b,i,j,k,l,m,nen=p->nen, dof=p->dof;

		//Create and inflate alpha
		PetscReal alfa[2];																//Container for values of alpha
		IGAPointFormValue(pAl,UAl,&alfa[0]);											//Fills the values of the alphas
		PetscReal fullAlfa[3][3]={0};
		fullAlfa[0][2]=alfa[0]; fullAlfa[1][2]=alfa[1];									//Expand Alfa to full 2nd order form, only non-zero elements

		PetscReal Chi0[4];																//Array to contain the vector chi(0)
		IGAPointFormValue(pChi,UChi,&Chi0[0]);											//Assign chi to its container
		PetscReal full_Chi[3][3]={0};
		full_Chi[0][0]=Chi0[0]; 	full_Chi[0][1]=Chi0[1];
		full_Chi[1][0]=Chi0[2]; 	full_Chi[1][1]=Chi0[3];

		PetscReal du[2][2];
		IGAPointFormGrad (pu,Uu,&du[0][0]);
		PetscReal full_du[3][3]={0};
		full_du[0][0]=du[0][0];
		full_du[0][1]=du[0][1];
		full_du[1][0]=du[1][0];
		full_du[1][1]=du[1][1];

		PetscReal dZ0[2][2];															//Array to store 3rd order partial derivative of z
		IGAPointFormGrad (pZu,UZu,&dZ0[0][0]);											//Calculate and store in array
		PetscReal full_dz[3][3]={0};
		full_dz[0][0]=dZ0[0][0];
		full_dz[0][1]=dZ0[0][1];
		full_dz[1][0]=dZ0[1][0];
		full_dz[1][1]=dZ0[1][1];

		PetscReal S[8];																		//Create array to receive S
		IGAPointFormValue(pS,US,&S[0]);														//This fills the values of S (remember that S has 8 non zero components in 2D)
		PetscReal fullS[3][3][3]={0};
		fullS[0][0][0]=S[0]; fullS[0][0][1]=S[1];											//Expand S to full 3rd order form, only non-zero elements
		fullS[0][1][0]=S[2]; fullS[0][1][1]=S[3];
		fullS[1][0][0]=S[4]; fullS[1][0][1]=S[5];
		fullS[1][1][0]=S[6]; fullS[1][1][1]=S[7];

		PetscReal Va[2];															//Seems that if I create this one with size 3, next line gives problems
		IGAPointFormValue(pVa,UVa,&Va[0]);													//This fills the values
		PetscReal fullVa[3]={0};	
		fullVa[0]=Va[0]; fullVa[1]=Va[1];

		PetscReal fullVs[3]={0};	
		fullVs[0]=Vs1/normS; fullVs[1]=-Vs2/normS;
		fullVs[0]=0.0;		 fullVs[1]=-1.0;			//THIS IS FOR DEBUGGING!!

		//Change for G=1
		const PetscReal nu=0.33;
		const PetscReal mu=1.0;
		const PetscReal lambda=2.0*mu*nu/(1.0-2.0*nu);

		PetscReal C[3][3][3][3]={0};														//Initialization of elastic tensor
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

		const PetscReal e[3][3][3]=
		{
			{{0.0,0.0,0.0},{0.0,0.0,1.0},{0.0,-1.0,0.0}},
			{{0.0,0.0,-1.0},{0.0,0.0,0.0},{1.0,0.0,0.0}},
			{{0.0,1.0,0.0},{-1.0,0.0,0.0},{0.0,0.0,0.0}}
		};

		for (j=0;j<3;j++)
		{
			for (k=0;k<3;k++)
			{
				for (l=0;l<3;l++)
				{
					for (m=0;m<3;m++)
					{
						fullAlfa[j][k]=fullAlfa[j][k]-fullS[j][l][m]*e[k][l][m];
					}
				}
			}
		}

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
						v[0][0]=N0[a];	v[0][1]=0.0;
						v[1][0]=0.0;	v[1][1]=0.0;
					}
					else if (i==1)
					{
						v[0][0]=0.0;	v[0][1]=N0[a];
						v[1][0]=0.0;	v[1][1]=0.0;
					}
					else if (i==2)
					{
						v[0][0]=0.0;	v[0][1]=0.0;
						v[1][0]=N0[a];	v[1][1]=0.0;
					}
					else if (i==3)
					{
						v[0][0]=0.0;	v[0][1]=0.0;
						v[1][0]=0.0;	v[1][1]=N0[a];
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
									FCS[a][i]+=C[j][k][l][m]*(full_du[l][m]-full_dz[l][m]-full_Chi[l][m])*v[j][k];
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

//Creation of solution systems
	PetscErrorCode  ierr;
	ierr = PetscInitialize(&argc,&argv,0,0);CHKERRQ(ierr);										//Always initialize PETSc
	PetscInt dir,side;
	PetscPrintf(PETSC_COMM_WORLD,"Start of PruebaV4S_V2 \n");

	PetscInt commsize,rank;
	ierr = MPI_Comm_size(PETSC_COMM_WORLD,&commsize);CHKERRQ(ierr);
	ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank);CHKERRQ(ierr);
//

//App context creation and some data
	//Mesh parameters (to fix specific point_S in z0 system)
	PetscInt b=300;				//Parameter to choose size of cores, must always be odd, core will be of size 1 unit, rest of the body will be of size b-1 units in each direction
	PetscReal Lx=10.0;
	PetscReal Ly=10.0;
	PetscInt  nx=b;
	PetscInt  ny=b;
	
	//alpha has to be at MOST 0.99 for things to work, (due to dt having to be strictly smaller than 0.5h), the smaller it is, the better accuracy, but more timesteps.
	PetscReal alpha=0.9*(1.0/8.0);
	PetscReal dt=alpha*0.5*fmin(Lx/nx,Ly/ny);

	AppCtx user;
	user.dt     = dt;

	PetscInt i=0;

	ierr=PetscPrintf(PETSC_COMM_WORLD,"Timestep is %f \n",dt);CHKERRQ(ierr);

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
		source = fopen("./PruebaV4S_V2.c","r");
		dest   = fopen("../Results/PruebaV4S_V2.c","w");

		while (0 < (bytes = fread(buffer, 1, sizeof(buffer), source)))
			fwrite(buffer, 1, bytes, dest);

		fclose(source);
		fclose(dest);
	}
//

//
//Recall that, for the systems to converge, dt<h/2, h=min(Lx/nx,Ly/ny) and the "less than" is strict 
//

//Creation of types and systems for the Initialization of S0
	//This block reads S from an external file. "generateInputS.c" should take care of generating the needed file. 
	PetscPrintf(PETSC_COMM_WORLD,"\nSystem for Initialization for S starting \n\n");
	T=time(NULL);
	tm=*localtime(&T);
	PetscPrintf(PETSC_COMM_WORLD,"Current time is %02d:%02d:%02d \n",tm.tm_hour,tm.tm_min,tm.tm_sec);
	IGA iga_S;
	ierr = IGACreate(PETSC_COMM_WORLD,&iga_S);CHKERRQ(ierr);
	ierr = IGASetDim(iga_S,2);CHKERRQ(ierr);														//Spatial dimension of the problem
	ierr = IGASetDof(iga_S,8);CHKERRQ(ierr);														//Number of degrees of freedom, per node
	ierr = IGASetOrder(iga_S,1);CHKERRQ(ierr);													//Number of spatial derivatives to calculate
	ierr = IGASetFromOptions(iga_S);CHKERRQ(ierr);												//Note: The order (or degree) of the shape functions is given by the mesh!
	ierr = IGARead(iga_S,"./geometry.dat");CHKERRQ(ierr);
	
	for (dir=0; dir<2; dir++)
	{
		ierr = IGASetRuleType(iga_S,dir,IGA_RULE_LEGENDRE);CHKERRQ(ierr);
		ierr = IGASetRuleSize(iga_S,dir,6);CHKERRQ(ierr);
	}
	ierr = IGASetUp(iga_S);CHKERRQ(ierr);

	Vec S0;
	ierr = IGACreateVec(iga_S,&S0);CHKERRQ(ierr);  

	//char nameS[]="/Input-S-2d-0.dat";			//Modify this so that the "-0.dat" part changes with the loop index 
	char nameS[512];
	sprintf(nameS,"%s%d%s","/Input-S-2d-",i,".dat");

	PetscPrintf(PETSC_COMM_WORLD,"%s\n",nameS);

	char pathS[1024];
	sprintf(pathS,"%s%s",direct,nameS);
	ierr = IGAReadVec(iga_S,S0,pathS); CHKERRQ(ierr);
//

//Creation of types and systems for the Helmholtz decomposition of S, curl part for S based input
	PetscPrintf(PETSC_COMM_WORLD,"\nSystem for curl part of Helmholtz of S starting \n\n");
	T=time(NULL);
	tm=*localtime(&T);
	PetscPrintf(PETSC_COMM_WORLD,"Current time is %02d:%02d:%02d \n",tm.tm_hour,tm.tm_min,tm.tm_sec);
	IGA iga_Sp;
	ierr = IGACreate(PETSC_COMM_WORLD,&iga_Sp);CHKERRQ(ierr);
	ierr = IGASetDim(iga_Sp,2);CHKERRQ(ierr);														//Spatial dimension of the problem
	ierr = IGASetDof(iga_Sp,8);CHKERRQ(ierr);														//Number of degrees of freedom, per node
	ierr = IGASetOrder(iga_Sp,1);CHKERRQ(ierr);														//Number of spatial derivatives to calculate
	ierr = IGASetFromOptions(iga_Sp);CHKERRQ(ierr);													//Note: The order (or degree) of the shape functions is given by the mesh!
	ierr = IGARead(iga_Sp,"./geometry.dat");CHKERRQ(ierr);
	
	for (dir=0; dir<2; dir++)
	{
		ierr = IGASetRuleType(iga_Sp,dir,IGA_RULE_LEGENDRE);CHKERRQ(ierr);
		ierr = IGASetRuleSize(iga_Sp,dir,6);CHKERRQ(ierr);
	}
	ierr = IGASetUp(iga_Sp);CHKERRQ(ierr);
	for (dir=0; dir<2; dir++) 
	{
		for (side=0; side<2; side++)
		{
			//ierr = IGASetBoundaryValue(iga,dir,side,dof,0.0);CHKERRQ(ierr);						// Dirichlet boundary conditions
			//ierr = IGASetBoundaryForm(iga,dir,side,PETSC_TRUE);CHKERRQ(ierr);						// Neumann boundary conditions
			if(dir==0 && side==0)
			{
				ierr = IGASetBoundaryValue(iga_Sp,dir,side,0,0.0);CHKERRQ(ierr);						// Dirichlet boundary conditions
				ierr = IGASetBoundaryValue(iga_Sp,dir,side,2,0.0);CHKERRQ(ierr);						// Dirichlet boundary conditions
				ierr = IGASetBoundaryValue(iga_Sp,dir,side,4,0.0);CHKERRQ(ierr);						// Dirichlet boundary conditions
				ierr = IGASetBoundaryValue(iga_Sp,dir,side,6,0.0);CHKERRQ(ierr);						// Dirichlet boundary conditions
			}
			if(dir==0 && side==1)
			{
				ierr = IGASetBoundaryValue(iga_Sp,dir,side,0,0.0);CHKERRQ(ierr);    					// Dirichlet boundary conditions
				ierr = IGASetBoundaryValue(iga_Sp,dir,side,2,0.0);CHKERRQ(ierr);    					// Dirichlet boundary conditions
				ierr = IGASetBoundaryValue(iga_Sp,dir,side,4,0.0);CHKERRQ(ierr);    					// Dirichlet boundary conditions
				ierr = IGASetBoundaryValue(iga_Sp,dir,side,6,0.0);CHKERRQ(ierr);    					// Dirichlet boundary conditions
			}
			if(dir==1 && side==0)
			{
				ierr = IGASetBoundaryValue(iga_Sp,dir,side,1,0.0);CHKERRQ(ierr);						// Dirichlet boundary conditions
				ierr = IGASetBoundaryValue(iga_Sp,dir,side,3,0.0);CHKERRQ(ierr);						// Dirichlet boundary conditions
				ierr = IGASetBoundaryValue(iga_Sp,dir,side,5,0.0);CHKERRQ(ierr);						// Dirichlet boundary conditions
				ierr = IGASetBoundaryValue(iga_Sp,dir,side,7,0.0);CHKERRQ(ierr);						// Dirichlet boundary conditions
			}
			if(dir==1 && side==1)
			{
				ierr = IGASetBoundaryValue(iga_Sp,dir,side,1,0.0);CHKERRQ(ierr);						// Dirichlet boundary conditions
				ierr = IGASetBoundaryValue(iga_Sp,dir,side,3,0.0);CHKERRQ(ierr);						// Dirichlet boundary conditions
				ierr = IGASetBoundaryValue(iga_Sp,dir,side,5,0.0);CHKERRQ(ierr);						// Dirichlet boundary conditions
				ierr = IGASetBoundaryValue(iga_Sp,dir,side,7,0.0);CHKERRQ(ierr);						// Dirichlet boundary conditions
			}
		}
	}
	
	Mat KSp;
	Vec S_perp,FSp;
	ierr = IGACreateMat(iga_Sp,&KSp);CHKERRQ(ierr);
	ierr = IGACreateVec(iga_Sp,&S_perp);CHKERRQ(ierr);
	ierr = IGACreateVec(iga_Sp,&FSp);CHKERRQ(ierr);

	IGAPoint        point_Sp, point_S;
	IGAElement      elem_Sp,  elem_S;						//element
	PetscReal       *Kloc_Sp,  *Floc_Sp;					//AA y BB
	PetscReal       *Kpoint_Sp,*Fpoint_Sp;					//KKK y FFF
	const PetscReal *array_S0_Sp;							//arrayU
	Vec  			local_S0_Sp;							//localU
	PetscReal       *S0_Sp;									//U

  	IGAFormSystem  wtfSp;
 	void           *wtf2Sp;

 	KSP ksp_Sp;
	ierr = IGACreateKSP(iga_Sp,&ksp_Sp);CHKERRQ(ierr);

	//Get local vectors S0 and arrays
	ierr = IGAGetLocalVecArray(iga_S,S0,&local_S0_Sp,&array_S0_Sp);CHKERRQ(ierr);

	//Element loop
	ierr = IGABeginElement(iga_Sp,&elem_Sp);CHKERRQ(ierr);
	ierr = IGABeginElement(iga_S,&elem_S);CHKERRQ(ierr);

	while (IGANextElement(iga_Sp,elem_Sp))
	{
		IGANextElement(iga_S,elem_S);
		ierr = IGAElementGetWorkMat(elem_Sp,&Kloc_Sp);CHKERRQ(ierr);
		ierr = IGAElementGetWorkVec(elem_Sp,&Floc_Sp);CHKERRQ(ierr);
		ierr = IGAElementGetValues(elem_S,array_S0_Sp,&S0_Sp);CHKERRQ(ierr);

		//FormSystem loop
		while (IGAElementNextFormSystem(elem_Sp,&wtfSp,&wtf2Sp)) 
		{
			//Quadrature loop
			ierr = IGAElementBeginPoint(elem_Sp,&point_Sp);CHKERRQ(ierr);
			ierr = IGAElementBeginPoint(elem_S,&point_S);CHKERRQ(ierr);
			while (IGAElementNextPoint(elem_Sp,point_Sp)) 
			{
				if(point_Sp->atboundary==0 && point_S->atboundary==0)
				{
					IGAElementNextPoint(elem_S,point_S);
					ierr = IGAPointGetWorkMat(point_Sp,&Kpoint_Sp);CHKERRQ(ierr);
					ierr = IGAPointGetWorkVec(point_Sp,&Fpoint_Sp);CHKERRQ(ierr);
					ierr = curlChiSM(point_Sp,Kpoint_Sp,NULL);CHKERRQ(ierr);
					ierr = curlChiSF(point_Sp,point_S,Fpoint_Sp,S0_Sp,NULL);CHKERRQ(ierr);
					ierr = IGAPointAddMat(point_Sp,Kpoint_Sp,Kloc_Sp);CHKERRQ(ierr);
					ierr = IGAPointAddVec(point_Sp,Fpoint_Sp,Floc_Sp);CHKERRQ(ierr);
				}
			}
			while (point_S->index != -1)
			{
				IGAElementNextPoint(elem_S,point_S);
			}
			ierr = IGAElementEndPoint(elem_Sp,&point_Sp);CHKERRQ(ierr);
			ierr = IGAElementEndPoint(elem_S,&point_S);CHKERRQ(ierr);
		}
		ierr = IGAElementFixSystem(elem_Sp,Kloc_Sp,Floc_Sp);CHKERRQ(ierr);
		ierr = IGAElementAssembleMat(elem_Sp,Kloc_Sp,KSp);CHKERRQ(ierr);
		ierr = IGAElementAssembleVec(elem_Sp,Floc_Sp,FSp);CHKERRQ(ierr);

	}
	IGANextElement(iga_S,elem_S);
	ierr = IGAEndElement(iga_Sp,&elem_Sp);CHKERRQ(ierr);
	ierr = IGAEndElement(iga_S,&elem_S);CHKERRQ(ierr);

	// Restore local vectors S0 and arrays
	ierr = IGARestoreLocalVecArray(iga_S,S0,&local_S0_Sp,&array_S0_Sp);CHKERRQ(ierr);

	//Form system matrix and vector
	ierr = MatAssemblyBegin(KSp,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyEnd  (KSp,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

	ierr = VecAssemblyBegin(FSp);CHKERRQ(ierr);
	ierr = VecAssemblyEnd  (FSp);CHKERRQ(ierr);

	ierr = KSPSetOperators(ksp_Sp,KSp,KSp);CHKERRQ(ierr);
	PC pc_Sp;
	ierr = KSPGetPC(ksp_Sp,&pc_Sp); CHKERRQ(ierr);
	ierr = PCSetType(pc_Sp,PCLU); CHKERRQ(ierr);
	ierr = PCFactorSetMatSolverType(pc_Sp,MATSOLVERMUMPS); CHKERRQ(ierr);
	//ierr = KSPSetFromOptions(ksp_Sp);CHKERRQ(ierr);
	//ierr = KSPSetTolerances(ksp_Sp,1.0e-12,5.0e-24,PETSC_DEFAULT,PETSC_DEFAULT);CHKERRQ(ierr);
	ierr = KSPSolve(ksp_Sp,FSp,S_perp);CHKERRQ(ierr);

	//ierr = KSPDestroy(&ksp_Sp);CHKERRQ(ierr);
	//ierr = MatDestroy(&KchiS);CHKERRQ(ierr);
	//ierr = VecDestroy(&FSp);CHKERRQ(ierr);
	
	char nameSp[512];
	sprintf(nameSp,"%s%d%s","/Sp-2d-",i,".dat");
	char pathSp[1024];
	sprintf(pathSp,"%s%s",direct,nameSp);
	ierr = IGAWriteVec(iga_Sp,S_perp,pathSp);CHKERRQ(ierr);
//

//System for norm of initial S
	PetscPrintf(PETSC_COMM_WORLD,"\nSystem for norm of S starting \n\n");
	T=time(NULL);
	tm=*localtime(&T);
	PetscPrintf(PETSC_COMM_WORLD,"Current time is %02d:%02d:%02d \n",tm.tm_hour,tm.tm_min,tm.tm_sec);
	IGA igaSnorm;
	ierr = IGACreate(PETSC_COMM_WORLD,&igaSnorm);CHKERRQ(ierr);
	ierr = IGASetDim(igaSnorm,2);CHKERRQ(ierr);													//Spatial dimension of the problem
	ierr = IGASetDof(igaSnorm,8);CHKERRQ(ierr);													//Number of degrees of freedom, per node
	ierr = IGASetOrder(igaSnorm,1);CHKERRQ(ierr);												//Number of spatial derivatives to calculate
	ierr = IGASetFromOptions(igaSnorm);CHKERRQ(ierr);											//Note: The order (or degree) of the shape functions is given by the mesh!
	ierr = IGARead(igaSnorm,"./geometry.dat");CHKERRQ(ierr);
	
	for (dir=0; dir<2; dir++)
	{
		ierr = IGASetRuleType(igaSnorm,dir,IGA_RULE_LEGENDRE);CHKERRQ(ierr);
		ierr = IGASetRuleSize(igaSnorm,dir,6);CHKERRQ(ierr);
	}
	ierr = IGASetUp(igaSnorm);CHKERRQ(ierr);

	//for (dir=0; dir<2; dir++) 
	//{
	//	for (side=0; side<2; side++) 
	//	{
	//		//ierr = IGASetBoundaryValue(iga,dir,side,dof,0.0);CHKERRQ(ierr);					// Dirichlet boundary conditions
	//		ierr = IGASetBoundaryForm(igaSnorm,dir,side,PETSC_TRUE);CHKERRQ(ierr);  				// Neumann boundary conditions
	//	}
	//}

	Vec FSNorm,FNorm;

	ierr = IGACreateVec(igaSnorm,&FSNorm);CHKERRQ(ierr);
	ierr = IGACreateVec(igaSnorm,&FNorm);CHKERRQ(ierr);

	IGAPoint		point_Snorm;															//point
	IGAElement		elem_Snorm;															//element
	PetscReal		*FlocSnorm,*FlocNorm;												//AA y BB
	PetscReal		*Fpoint_Snorm,*FpointNorm;											//KKK y FFF
	const PetscReal	*arrayS0Snorm;														//arrayU
	Vec				localS0Snorm;														//localU
	PetscReal		*S0Snorm;																//U

  	IGAFormSystem	wtfSnorm;
 	void			*wtf2Snorm;

	// Get local vectors Z0 and Chi0 and arrays
	ierr = IGAGetLocalVecArray(iga_S,S0,&localS0Snorm,&arrayS0Snorm);CHKERRQ(ierr);

	// Element loop
	ierr = IGABeginElement(igaSnorm,&elem_Snorm);CHKERRQ(ierr);
	ierr = IGABeginElement(iga_S,&elem_S);CHKERRQ(ierr);

	while (IGANextElement(igaSnorm,elem_Snorm)) 
	{
		IGANextElement(iga_S,elem_S);

		ierr = IGAElementGetWorkVec(elem_Snorm,&FlocNorm);CHKERRQ(ierr);
		ierr = IGAElementGetWorkVec(elem_Snorm,&FlocSnorm);CHKERRQ(ierr);
		ierr = IGAElementGetValues(elem_S,arrayS0Snorm,&S0Snorm);CHKERRQ(ierr);

		// FormSystem loop
		while (IGAElementNextFormSystem(elem_Snorm,&wtfSnorm,&wtf2Snorm)) 
		{
		// Quadrature loop
			ierr = IGAElementBeginPoint(elem_Snorm,&point_Snorm);CHKERRQ(ierr);
			ierr = IGAElementBeginPoint(elem_S,&point_S);CHKERRQ(ierr);

			while (IGAElementNextPoint(elem_Snorm,point_Snorm))
			{
				if(point_Snorm->atboundary==1)
				{

				}
				if(point_Snorm->atboundary==0 && point_S->atboundary==0)
				{
					IGAElementNextPoint(elem_S,point_S);

					ierr = IGAPointGetWorkVec(point_Snorm,&Fpoint_Snorm);CHKERRQ(ierr);
					ierr = IGAPointGetWorkVec(point_Snorm,&FpointNorm);CHKERRQ(ierr);
					//	   Snorm(IGAPoint pS,PetscReal *SNorm,PetscReal *Norm,PetscReal *US,void *ctx)
					ierr = Snorm(point_Snorm,Fpoint_Snorm,FpointNorm,S0Snorm,NULL);CHKERRQ(ierr);
					ierr = IGAPointAddVec(point_Snorm,Fpoint_Snorm,FlocSnorm);CHKERRQ(ierr);
					ierr = IGAPointAddVec(point_Snorm,FpointNorm,FlocNorm);CHKERRQ(ierr);
				}
			}
			while (point_S->index != -1)
			{
				IGAElementNextPoint(elem_S,point_S);
			}
			ierr = IGAElementEndPoint(elem_Snorm,&point_Snorm);CHKERRQ(ierr);
			ierr = IGAElementEndPoint(elem_S,&point_S);CHKERRQ(ierr);
		}
		ierr = IGAElementAssembleVec(elem_Snorm,FlocSnorm,FSNorm);CHKERRQ(ierr);
		ierr = IGAElementAssembleVec(elem_Snorm,FlocNorm,FNorm);CHKERRQ(ierr);
	}
	IGANextElement(iga_S,elem_S);

	ierr = IGAEndElement(igaSnorm,&elem_Snorm);CHKERRQ(ierr);
	ierr = IGAEndElement(iga_S,&elem_S);CHKERRQ(ierr);

	// Restore local vectors u, Z0, Chi0 and arrays
	ierr = IGARestoreLocalVecArray(iga_S,S0,&localS0Snorm,&arrayS0Snorm);CHKERRQ(ierr);

	ierr = VecAssemblyBegin(FSNorm);CHKERRQ(ierr);
	ierr = VecAssemblyEnd  (FSNorm);CHKERRQ(ierr);

	ierr = VecAssemblyBegin(FNorm);CHKERRQ(ierr);
	ierr = VecAssemblyEnd  (FNorm);CHKERRQ(ierr);

	//Here integrate xi(S)
	PetscReal IntS0;
	ierr = VecSum(FNorm,&IntS0);CHKERRQ(ierr);
	ierr = PetscPrintf(PETSC_COMM_WORLD,"IntXiS = %f \n",IntS0);CHKERRQ(ierr);

	//Here integrate S(i,j,2)
	PetscInt ind1,ind2,vec_len,*indices;
	PetscReal *values,*IntSField, S_norm;
	Vec vec_work;

	ierr = VecGetSize(FSNorm,&vec_len);
	ierr = PetscCalloc1(vec_len,&indices);
	ierr = PetscCalloc1(vec_len,&values);
	ierr = PetscCalloc1(8,&IntSField);

	for(ind2=0;ind2<8;ind2++)
	{
		ierr = VecDuplicate(FSNorm,&vec_work);CHKERRQ(ierr);
		ierr = VecCopy(FSNorm,vec_work);CHKERRQ(ierr);

		ierr = VecSetOption(vec_work,VEC_IGNORE_NEGATIVE_INDICES,PETSC_TRUE);CHKERRQ(ierr);
		ierr = VecAssemblyBegin(vec_work);CHKERRQ(ierr);
		ierr = VecAssemblyEnd  (vec_work);CHKERRQ(ierr);

		for(ind1=0;ind1<vec_len;ind1++)
		{
			indices[ind1]=ind1;
			values[ind1]=0.0;
		}

		for(ind1=0;ind1<vec_len/8-1;ind1++)
		{
			indices[8*ind1+ind2]=-1;
		}

		ierr = VecSetValues(vec_work,vec_len,indices,values,INSERT_VALUES);CHKERRQ(ierr);

		ierr = VecAssemblyBegin(vec_work);CHKERRQ(ierr);
		ierr = VecAssemblyEnd  (vec_work);CHKERRQ(ierr);

		ierr = VecSum(vec_work,&IntSField[ind2]);CHKERRQ(ierr);
	}
	for(ind1=0;ind1<8;ind1++)
	{
		IntSField[ind1]=IntSField[ind1]/IntS0;
		ierr = PetscPrintf(PETSC_COMM_WORLD,"IntS[%d] = %f \n",ind1,IntSField[ind1]);CHKERRQ(ierr);
	}
	S_norm=0.0;
	for(ind1=0;ind1<8;ind1++)
	{
		S_norm=S_norm+IntSField[ind1]*IntSField[ind1];
	}
	S_norm=sqrt(S_norm);

	ierr = PetscPrintf(PETSC_COMM_WORLD,"Norm of S = %f \n",S_norm);CHKERRQ(ierr);

	ierr = VecDestroy(&FNorm);CHKERRQ(ierr);
	ierr = VecDestroy(&FSNorm);CHKERRQ(ierr);
	ierr = IGADestroy(&igaSnorm);CHKERRQ(ierr);
//

//Creation of types and systems for the Helmholtz decomposition of S, grad part
	//System for chiS
	PetscPrintf(PETSC_COMM_WORLD,"\n System for grad part of Helmholtz of S starting \n\n");
	IGA iga_ZS;
	ierr = IGACreate(PETSC_COMM_WORLD,&iga_ZS);CHKERRQ(ierr);
	ierr = IGASetDim(iga_ZS,2);CHKERRQ(ierr);														//Spatial dimension of the problem
	ierr = IGASetDof(iga_ZS,4);CHKERRQ(ierr);														//Number of degrees of freedom, per node
	ierr = IGASetOrder(iga_ZS,2);CHKERRQ(ierr);														//Number of spatial derivatives to calculate
	ierr = IGASetFromOptions(iga_ZS);CHKERRQ(ierr);													//Note: The order (or degree) of the shape functions is given by the mesh!
	ierr = IGARead(iga_ZS,"./geometry2.dat");CHKERRQ(ierr);
	ierr = IGASetUp(iga_ZS);CHKERRQ(ierr);
	
	//PetscInt dir,side;
	for (dir=0; dir<2; dir++)
	{
		for (side=0; side<2; side++)
		{
			//ierr = IGASetBoundaryValue(iga,dir,side,dof,0.0);CHKERRQ(ierr);    				// Dirichlet boundary conditions
			ierr = IGASetBoundaryForm(iga_ZS,dir,side,PETSC_TRUE);CHKERRQ(ierr);  				// Neumann boundary conditions
		}
	}
	
	Mat K_ZS;
	Vec ZS0,F_ZS;
	ierr = IGACreateMat(iga_ZS,&K_ZS);CHKERRQ(ierr);
	ierr = IGACreateVec(iga_ZS,&ZS0);CHKERRQ(ierr);
	ierr = IGACreateVec(iga_ZS,&F_ZS);CHKERRQ(ierr);

	IGAPoint        point_ZS;
	IGAElement      elem_ZS;					//element
	PetscReal       *Kloc_ZS,*Floc_ZS;			//AA y BB
	PetscReal       *Kpoint_ZS,*Fpoint_ZS;		//KKK y FFF
	const PetscReal *array_S0_ZS;				//arrayU
	Vec  			local_S0_ZS;				//localU
	PetscReal       *S0_ZS;						//U

  	IGAFormSystem  wtfZS;
 	void           *wtf2ZS;

 	KSP ksp_ZS;
	ierr = IGACreateKSP(iga_ZS,&ksp_ZS);CHKERRQ(ierr);

	//Get local vectors S0 and arrays
	ierr = IGAGetLocalVecArray(iga_S,S0,&local_S0_ZS,&array_S0_ZS);CHKERRQ(ierr);

	//Element loop
	ierr = IGABeginElement(iga_ZS,&elem_ZS);CHKERRQ(ierr);
	ierr = IGABeginElement(iga_S,&elem_S);CHKERRQ(ierr);

	while (IGANextElement(iga_ZS,elem_ZS))
	{
		IGANextElement(iga_S,elem_S);
		ierr = IGAElementGetWorkMat(elem_ZS,&Kloc_ZS);CHKERRQ(ierr);
		ierr = IGAElementGetWorkVec(elem_ZS,&Floc_ZS);CHKERRQ(ierr);
		ierr = IGAElementGetValues(elem_S,array_S0_ZS,&S0_ZS);CHKERRQ(ierr);

		//FormSystem loop
		while (IGAElementNextFormSystem(elem_ZS,&wtfZS,&wtf2ZS)) 
		{
			//Quadrature loop
			ierr = IGAElementBeginPoint(elem_ZS,&point_ZS);CHKERRQ(ierr);
			ierr = IGAElementBeginPoint(elem_S,&point_S);CHKERRQ(ierr);
			if(point_ZS->atboundary==0 && point_S->atboundary==0)
			{
				while (IGAElementNextPoint(elem_ZS,point_ZS)) 
				{
					IGAElementNextPoint(elem_S,point_S);
					ierr = IGAPointGetWorkMat(point_ZS,&Kpoint_ZS);CHKERRQ(ierr);
					ierr = IGAPointGetWorkVec(point_ZS,&Fpoint_ZS);CHKERRQ(ierr);
					ierr = gradZSM(point_ZS,Kpoint_ZS,NULL);CHKERRQ(ierr);
					ierr = gradZSF(point_ZS,point_S,Fpoint_ZS,S0_ZS,NULL);CHKERRQ(ierr);
					ierr = IGAPointAddMat(point_ZS,Kpoint_ZS,Kloc_ZS);CHKERRQ(ierr);
					ierr = IGAPointAddVec(point_ZS,Fpoint_ZS,Floc_ZS);CHKERRQ(ierr);
				}
				IGAElementNextPoint(elem_S,point_S);

				ierr = IGAElementEndPoint(elem_ZS,&point_ZS);CHKERRQ(ierr);
				ierr = IGAElementEndPoint(elem_S,&point_S);CHKERRQ(ierr);
			}
		}
		ierr = IGAElementAssembleMat(elem_ZS,Kloc_ZS,K_ZS);CHKERRQ(ierr);
		ierr = IGAElementAssembleVec(elem_ZS,Floc_ZS,F_ZS);CHKERRQ(ierr);

	}
	IGANextElement(iga_S,elem_S);
	ierr = IGAEndElement(iga_ZS,&elem_ZS);CHKERRQ(ierr);
	ierr = IGAEndElement(iga_S,&elem_S);CHKERRQ(ierr);
	
	// Restore local vectors S0 and arrays
	ierr = IGARestoreLocalVecArray(iga_S,S0,&local_S0_ZS,&array_S0_ZS);CHKERRQ(ierr);
	
	//Impose Dirichlet condition on a single point
	//First we replace rows and columns with associated rows and columns from identity matrix
	PetscInt mz,nz;
	ierr = MatGetSize(K_ZS,&nz,&mz);CHKERRQ(ierr);
	ierr = MatSetOption(K_ZS, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);CHKERRQ(ierr);
	ierr = MatAssemblyBegin(K_ZS,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyEnd  (K_ZS,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);	

	PetscInt rowsZ, *colsZ;
	PetscReal valZ, *valsZ;

	ierr = PetscCalloc1(mz,&colsZ);CHKERRQ(ierr);													//This creates the array colsZ having mz elements, memory aligned and zeroed
	ierr = PetscCalloc1(mz,&valsZ);CHKERRQ(ierr);													//the type in the arrays is picked up from the declaration of the pointer just before this

	for(int i=0;i<mz;i++)
	{
		colsZ[i]=i;
		valsZ[i]=0.0;
	}

	rowsZ=0;
	valsZ[rowsZ]=1.0e6;																				//This basically turns row 0 into the equation 1e6*u0=0, the big number makes sure the value is strongly enforced by the iterative solver 
	ierr = MatSetValues(K_ZS,1,&rowsZ,mz,colsZ,valsZ,INSERT_VALUES);CHKERRQ(ierr);
	ierr = MatSetValues(K_ZS,nz,colsZ,1,&rowsZ,valsZ,INSERT_VALUES);CHKERRQ(ierr);
	valsZ[rowsZ]=0.0;
	rowsZ=1;
	valsZ[rowsZ]=1.0e6;
	ierr = MatSetValues(K_ZS,1,&rowsZ,mz,colsZ,valsZ,INSERT_VALUES);CHKERRQ(ierr);
	ierr = MatSetValues(K_ZS,nz,colsZ,1,&rowsZ,valsZ,INSERT_VALUES);CHKERRQ(ierr);
	valsZ[rowsZ]=0.0;
	rowsZ=2;
	valsZ[rowsZ]=1.0e6;
	ierr = MatSetValues(K_ZS,1,&rowsZ,mz,colsZ,valsZ,INSERT_VALUES);CHKERRQ(ierr);
	ierr = MatSetValues(K_ZS,nz,colsZ,1,&rowsZ,valsZ,INSERT_VALUES);CHKERRQ(ierr);
	valsZ[rowsZ]=0.0;
	rowsZ=3;
	valsZ[rowsZ]=1.0e6;
	ierr = MatSetValues(K_ZS,1,&rowsZ,mz,colsZ,valsZ,INSERT_VALUES);CHKERRQ(ierr);
	ierr = MatSetValues(K_ZS,nz,colsZ,1,&rowsZ,valsZ,INSERT_VALUES);CHKERRQ(ierr);
	valsZ[rowsZ]=0.0;

	//After assigning values to matrix, collect it again
	ierr = MatAssemblyBegin(K_ZS,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyEnd  (K_ZS,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

	//Here we set values to the vector directly
	ierr = VecAssemblyBegin(F_ZS);CHKERRQ(ierr);
	ierr = VecAssemblyEnd  (F_ZS);CHKERRQ(ierr);

	rowsZ=0;
	valZ=0.0;
	ierr = VecSetValue(F_ZS,rowsZ,valZ,INSERT_VALUES);CHKERRQ(ierr);
	rowsZ=1;
	valZ=0.0;
	ierr = VecSetValue(F_ZS,rowsZ,valZ,INSERT_VALUES);CHKERRQ(ierr);
	rowsZ=2;
	valZ=0.0;
	ierr = VecSetValue(F_ZS,rowsZ,valZ,INSERT_VALUES);CHKERRQ(ierr);
	rowsZ=3;
	valZ=0.0;
	ierr = VecSetValue(F_ZS,rowsZ,valZ,INSERT_VALUES);CHKERRQ(ierr);

	ierr = VecAssemblyBegin(F_ZS);CHKERRQ(ierr);
	ierr = VecAssemblyEnd  (F_ZS);CHKERRQ(ierr);
	//

	ierr = KSPSetOperators(ksp_ZS,K_ZS,K_ZS);CHKERRQ(ierr);
	//PC pc_ZS;
	//ierr = KSPGetPC(ksp_ZS,&pc_ZS); CHKERRQ(ierr);
	//ierr = PCSetType(pc_ZS,PCLU); CHKERRQ(ierr);
	//ierr = PCFactorSetMatSolverType(pc_ZS,MATSOLVERMUMPS); CHKERRQ(ierr);
	ierr = KSPSetFromOptions(ksp_ZS);CHKERRQ(ierr);
	ierr = KSPSetTolerances(ksp_ZS,1.0e-12,1.0e-24,PETSC_DEFAULT,PETSC_DEFAULT);CHKERRQ(ierr);
	ierr = KSPSolve(ksp_ZS,F_ZS,ZS0);CHKERRQ(ierr);

	ierr = KSPDestroy(&ksp_ZS);CHKERRQ(ierr);
	ierr = MatDestroy(&K_ZS);CHKERRQ(ierr);
	ierr = VecDestroy(&F_ZS);CHKERRQ(ierr);
	
	char nameZS[]="/ZS-2d-0.dat";
	char pathZS[512];
	sprintf(pathZS,"%s%s",direct,nameZS);
	ierr = IGAWriteVec(iga_ZS,ZS0,pathZS);CHKERRQ(ierr);
//

//Creation of types and systems for the L2 projection of alpha_tilde and alpha_hat
	//System for Alfa
	PetscPrintf(PETSC_COMM_WORLD,"\nSystem for L2 projection for Alfa-Sp:X starting \n\n");
	T=time(NULL);
	tm=*localtime(&T);
	PetscPrintf(PETSC_COMM_WORLD,"Current time is %02d:%02d:%02d \n",tm.tm_hour,tm.tm_min,tm.tm_sec);
	IGA iga_Al_hat;
	ierr = IGACreate(PETSC_COMM_WORLD,&iga_Al_hat);CHKERRQ(ierr);
	ierr = IGASetDim(iga_Al_hat,2);CHKERRQ(ierr);														//Spatial dimension of the problem
	ierr = IGASetDof(iga_Al_hat,2);CHKERRQ(ierr);														//Number of degrees of freedom, per node
	ierr = IGASetOrder(iga_Al_hat,1);CHKERRQ(ierr);													//Number of spatial derivatives to calculate
	ierr = IGASetFromOptions(iga_Al_hat);CHKERRQ(ierr);												//Note: The order (or degree) of the shape functions is given by the mesh!
	ierr = IGARead(iga_Al_hat,"./geometry.dat");CHKERRQ(ierr);
	
	for (dir=0; dir<2; dir++)
	{
		ierr = IGASetRuleType(iga_Al_hat,dir,IGA_RULE_LEGENDRE);CHKERRQ(ierr);
		ierr = IGASetRuleSize(iga_Al_hat,dir,6);CHKERRQ(ierr);
	}
	ierr = IGASetUp(iga_Al_hat);CHKERRQ(ierr);
	//PetscInt dir,side;
	for (dir=0; dir<2; dir++) 
	{
		for (side=0; side<2; side++) 
		{
			//ierr = IGASetBoundaryValue(iga,dir,side,dof,0.0);CHKERRQ(ierr);    				// Dirichlet boundary conditions
			//ierr = IGASetBoundaryForm(iga_Al_hat,dir,side,PETSC_TRUE);CHKERRQ(ierr);  				// Neumann boundary conditions
			//Not required for L2 projection
		}
	}
	
	Mat KL2_2GDL;
	Vec al_hat,al_tilde,alInput,F_alpha_hat,F_alpha_tilde;
	ierr = IGACreateMat(iga_Al_hat,&KL2_2GDL);CHKERRQ(ierr);
	ierr = IGACreateVec(iga_Al_hat,&al_hat);CHKERRQ(ierr);
	ierr = IGACreateVec(iga_Al_hat,&al_tilde);CHKERRQ(ierr);
	ierr = IGACreateVec(iga_Al_hat,&alInput);CHKERRQ(ierr);
	ierr = IGACreateVec(iga_Al_hat,&F_alpha_hat);CHKERRQ(ierr);
	ierr = IGACreateVec(iga_Al_hat,&F_alpha_tilde);CHKERRQ(ierr);

	IGAPoint        point_Al_hat;
	IGAElement      elem_Al_hat;											//element
	PetscReal       *Kloc_Al_hat,*Floc_Al_hat,*Floc_Al_tilde;				//AA y BB
	PetscReal       *Kpoint_Al_hat,*Fpoint_Al_hat,*Fpoint_Al_tilde;			//KKK y FFF
	const PetscReal *array_Sperp_Al_hat,*array_S_Al_tilde;						//arrayU
	Vec  			local_Sperp_Al_hat,local_S_Al_tilde;						//localU
	PetscReal       *Sperp_Al_hat,*S_Al_tilde;									//U

  	IGAFormSystem  wtfAlhat;
 	void           *wtf2Alhat;

 	KSP ksp_Al_hat;
	ierr = IGACreateKSP(iga_Al_hat,&ksp_Al_hat);CHKERRQ(ierr);

	//Get local vectors S0 and arrays
	ierr = IGAGetLocalVecArray(iga_Sp,S_perp,&local_Sperp_Al_hat,&array_Sperp_Al_hat);CHKERRQ(ierr);
	ierr = IGAGetLocalVecArray(iga_S,S0,&local_S_Al_tilde,&array_S_Al_tilde);CHKERRQ(ierr);

	//Element loop
	ierr = IGABeginElement(iga_Al_hat,&elem_Al_hat);CHKERRQ(ierr);
	ierr = IGABeginElement(iga_Sp,&elem_Sp);CHKERRQ(ierr);
	ierr = IGABeginElement(iga_S,&elem_S);CHKERRQ(ierr);

	while (IGANextElement(iga_Al_hat,elem_Al_hat))
	{
		IGANextElement(iga_Sp,elem_Sp);
		IGANextElement(iga_S,elem_S);
		ierr = IGAElementGetWorkMat(elem_Al_hat,&Kloc_Al_hat);CHKERRQ(ierr);
		ierr = IGAElementGetWorkVec(elem_Al_hat,&Floc_Al_hat);CHKERRQ(ierr);
		ierr = IGAElementGetWorkVec(elem_Al_hat,&Floc_Al_tilde);CHKERRQ(ierr);
		ierr = IGAElementGetValues(elem_Sp,array_Sperp_Al_hat,&Sperp_Al_hat);CHKERRQ(ierr);
		ierr = IGAElementGetValues(elem_S,array_S_Al_tilde,&S_Al_tilde);CHKERRQ(ierr);

		//FormSystem loop
		while (IGAElementNextFormSystem(elem_Al_hat,&wtfAlhat,&wtf2Alhat)) 
		{
			//Quadrature loop
			ierr = IGAElementBeginPoint(elem_Al_hat,&point_Al_hat);CHKERRQ(ierr);
			ierr = IGAElementBeginPoint(elem_Sp,&point_Sp);CHKERRQ(ierr);
			ierr = IGAElementBeginPoint(elem_S,&point_S);CHKERRQ(ierr);
			
			while (IGAElementNextPoint(elem_Al_hat,point_Al_hat)) 
			{
				if(point_Al_hat->atboundary==0 && point_Sp->atboundary==0)
				{
					IGAElementNextPoint(elem_Sp,point_Sp);
					IGAElementNextPoint(elem_S,point_S);

					ierr = IGAPointGetWorkMat(point_Al_hat,&Kpoint_Al_hat);CHKERRQ(ierr);
					ierr = IGAPointGetWorkVec(point_Al_hat,&Fpoint_Al_hat);CHKERRQ(ierr);
					ierr = IGAPointGetWorkVec(point_Al_hat,&Fpoint_Al_tilde);CHKERRQ(ierr);
					ierr = L2Projection2DOF(point_Al_hat,Kpoint_Al_hat,NULL);CHKERRQ(ierr);
					ierr = Sp_X(point_Al_hat,point_Sp,Fpoint_Al_hat,Sperp_Al_hat,NULL);CHKERRQ(ierr);
					ierr = S_X(point_Al_hat,point_Sp,Fpoint_Al_tilde,S_Al_tilde,NULL);CHKERRQ(ierr);
					ierr = IGAPointAddMat(point_Al_hat,Kpoint_Al_hat,Kloc_Al_hat);CHKERRQ(ierr);
					ierr = IGAPointAddVec(point_Al_hat,Fpoint_Al_hat,Floc_Al_hat);CHKERRQ(ierr);
					ierr = IGAPointAddVec(point_Al_hat,Fpoint_Al_tilde,Floc_Al_tilde);CHKERRQ(ierr);
				}
			}
			IGAElementNextPoint(elem_Sp,point_Sp);
			IGAElementNextPoint(elem_S,point_S);

			ierr = IGAElementEndPoint(elem_Al_hat,&point_Al_hat);CHKERRQ(ierr);
			ierr = IGAElementEndPoint(elem_Sp,&point_Sp);CHKERRQ(ierr);
			ierr = IGAElementEndPoint(elem_S,&point_S);CHKERRQ(ierr);
		}
		//ierr = IGAElementFixSystem(elem_Al_hat,Kloc_Al_hat,Floc_Al_hat);CHKERRQ(ierr);					//This sets Dirichlet condition, not used for L2 proj
		ierr = IGAElementAssembleMat(elem_Al_hat,Kloc_Al_hat,KL2_2GDL);CHKERRQ(ierr);
		ierr = IGAElementAssembleVec(elem_Al_hat,Floc_Al_hat,F_alpha_hat);CHKERRQ(ierr);
		ierr = IGAElementAssembleVec(elem_Al_hat,Floc_Al_tilde,F_alpha_tilde);CHKERRQ(ierr);

	}
	IGANextElement(iga_Sp,elem_Sp);
	ierr = IGAEndElement(iga_Al_hat,&elem_Al_hat);CHKERRQ(ierr);
	ierr = IGAEndElement(iga_Sp,&elem_Sp);CHKERRQ(ierr);
	ierr = IGAEndElement(iga_S,&elem_S);CHKERRQ(ierr);

	// Restore local vectors S0 and arrays
	ierr = IGARestoreLocalVecArray(iga_Sp,S_perp,&local_Sperp_Al_hat,&array_Sperp_Al_hat);CHKERRQ(ierr);			//CHANGE S_perp to S_perp, more descriptive name
	ierr = IGARestoreLocalVecArray(iga_S,S0,&local_S_Al_tilde,&array_S_Al_tilde);CHKERRQ(ierr);			//CHANGE S_perp to S_perp, more descriptive name

	//Form system matrix and vector
	ierr = MatAssemblyBegin(KL2_2GDL,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyEnd  (KL2_2GDL,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	
	ierr = VecAssemblyBegin(F_alpha_hat);CHKERRQ(ierr);
	ierr = VecAssemblyEnd  (F_alpha_hat);CHKERRQ(ierr);

	ierr = VecAssemblyBegin(F_alpha_tilde);CHKERRQ(ierr);
	ierr = VecAssemblyEnd  (F_alpha_tilde);CHKERRQ(ierr);

	ierr = KSPSetOperators(ksp_Al_hat,KL2_2GDL,KL2_2GDL);CHKERRQ(ierr);
	ierr = KSPSetFromOptions(ksp_Al_hat);CHKERRQ(ierr);
	ierr = KSPSetTolerances(ksp_Al_hat,1.0e-25,1.0e-40,PETSC_DEFAULT,PETSC_DEFAULT);CHKERRQ(ierr);
	ierr = KSPSolve(ksp_Al_hat,F_alpha_hat,al_hat);CHKERRQ(ierr);
	ierr = KSPSolve(ksp_Al_hat,F_alpha_tilde,al_tilde);CHKERRQ(ierr);

	char nameAlInput[]="/Input-Al-2d-0.dat";
	char pathAlInput[1024];
	sprintf(pathAlInput,"%s%s",direct,nameAlInput);
	ierr = IGAReadVec(iga_Al_hat,alInput,pathAlInput); CHKERRQ(ierr);

	ierr = VecAXPY(al_hat,1.0,alInput); CHKERRQ(ierr);
	ierr = VecAXPY(al_tilde,1.0,alInput); CHKERRQ(ierr);

	//ierr = KSPDestroy(&ksp_Al_hat);CHKERRQ(ierr);
	//ierr = MatDestroy(&KL2_2GDL);CHKERRQ(ierr);
	//ierr = VecDestroy(&F_alpha_hat);CHKERRQ(ierr);
	
	char nameAlp[512];
	sprintf(nameAlp,"%s%d%s","/Al_hat-2d-",i,".dat");
	char pathAlp[1024];
	sprintf(pathAlp,"%s%s",direct,nameAlp);
	ierr = IGAWriteVec(iga_Al_hat,al_hat,pathAlp);CHKERRQ(ierr);

	char nameAlt[512];
	sprintf(nameAlt,"%s%d%s","/Al_tilde-2d-",i,".dat");
	char pathAlt[1024];
	sprintf(pathAlt,"%s%s",direct,nameAlt);
	ierr = IGAWriteVec(iga_Al_hat,al_tilde,pathAlt);CHKERRQ(ierr);

	//ierr = IGADestroy(&iga_Sp);CHKERRQ(ierr);		//This system is one of the more memory intensive. It's not required from here down, so better to destroy it.
//

//Creation of types and systems for the Helmholtz decomposition of Up (or Ue), curl part
	PetscPrintf(PETSC_COMM_WORLD,"\nSystem for curl part of Helmholtz of Up starting \n\n");
	T=time(NULL);
	tm=*localtime(&T);
	PetscPrintf(PETSC_COMM_WORLD,"Current time is %02d:%02d:%02d \n",tm.tm_hour,tm.tm_min,tm.tm_sec);
	IGA iga_chi;
	ierr = IGACreate(PETSC_COMM_WORLD,&iga_chi);CHKERRQ(ierr);
	ierr = IGASetDim(iga_chi,2);CHKERRQ(ierr);														//Spatial dimension of the problem
	ierr = IGASetDof(iga_chi,4);CHKERRQ(ierr);														//Number of degrees of freedom, per node
	ierr = IGASetOrder(iga_chi,1);CHKERRQ(ierr);													//Number of spatial derivatives to calculate
	ierr = IGASetFromOptions(iga_chi);CHKERRQ(ierr);												//Note: The order (or degree) of the shape functions is given by the mesh!
	ierr = IGARead(iga_chi,"./geometry.dat");CHKERRQ(ierr);
	
	for (dir=0; dir<2; dir++)
	{
		ierr = IGASetRuleType(iga_chi,dir,IGA_RULE_LEGENDRE);CHKERRQ(ierr);
		ierr = IGASetRuleSize(iga_chi,dir,6);CHKERRQ(ierr);
	}
	ierr = IGASetUp(iga_chi);CHKERRQ(ierr);

	//ierr = IGASetBoundaryValue(iga,dir,side,dof,val);CHKERRQ(ierr);
	//dir=0 and side=0
	ierr = IGASetBoundaryValue(iga_chi,0,0,0,0.0);CHKERRQ(ierr);					// Dirichlet boundary conditions to get chi \dot n =0
	ierr = IGASetBoundaryValue(iga_chi,0,0,2,0.0);CHKERRQ(ierr);
	//dir=0 and side=1
	ierr = IGASetBoundaryValue(iga_chi,0,1,0,0.0);CHKERRQ(ierr);					// Dirichlet boundary conditions to get chi \dot n =0
	ierr = IGASetBoundaryValue(iga_chi,0,1,2,0.0);CHKERRQ(ierr);
	//dir=1 and side=0
	ierr = IGASetBoundaryValue(iga_chi,1,0,1,0.0);CHKERRQ(ierr);					// Dirichlet boundary conditions to get chi \dot n =0
	ierr = IGASetBoundaryValue(iga_chi,1,0,3,0.0);CHKERRQ(ierr);
	//dir=1 and side=1
	ierr = IGASetBoundaryValue(iga_chi,1,1,1,0.0);CHKERRQ(ierr);					// Dirichlet boundary conditions to get chi \dot n =0
	ierr = IGASetBoundaryValue(iga_chi,1,1,3,0.0);CHKERRQ(ierr);
	
	Mat K_chi;
	Vec chi,F_chi;
	ierr = IGACreateMat(iga_chi,&K_chi);CHKERRQ(ierr);
	ierr = IGACreateVec(iga_chi,&chi);CHKERRQ(ierr);
	ierr = IGACreateVec(iga_chi,&F_chi);CHKERRQ(ierr);

	IGAPoint        point_chi;
	IGAElement      elem_chi;							//element
	PetscReal       *Kloc_chi,*Floc_chi;				//AA y BB
	PetscReal       *Kpoint_chi,*Fpoint_chi;			//KKK y FFF
	const PetscReal *array_Al_chi;		//arrayU
	Vec  			local_Al_chi;		//localU
	PetscReal       *Al_chi;							//U

  	IGAFormSystem  wtfchi;
 	void           *wtf2chi;

	//Get local vectors Al0 and arrays
	ierr = IGAGetLocalVecArray(iga_Al_hat,al_hat,&local_Al_chi,&array_Al_chi);CHKERRQ(ierr);

	//Element loop
	ierr = IGABeginElement(iga_chi,&elem_chi);CHKERRQ(ierr);
	ierr = IGABeginElement(iga_Al_hat,&elem_Al_hat);CHKERRQ(ierr);

	while (IGANextElement(iga_chi,elem_chi))
	{
		IGANextElement(iga_Al_hat,elem_Al_hat);
		ierr = IGAElementGetWorkMat(elem_chi,&Kloc_chi);CHKERRQ(ierr);
		ierr = IGAElementGetWorkVec(elem_chi,&Floc_chi);CHKERRQ(ierr);
		ierr = IGAElementGetValues(elem_Al_hat,array_Al_chi,&Al_chi);CHKERRQ(ierr);

		//FormSystem loop
		while (IGAElementNextFormSystem(elem_chi,&wtfchi,&wtf2chi)) 
		{
			//Quadrature loop
			ierr = IGAElementBeginPoint(elem_chi,&point_chi);CHKERRQ(ierr);
			ierr = IGAElementBeginPoint(elem_Al_hat,&point_Al_hat);CHKERRQ(ierr);

			while (IGAElementNextPoint(elem_chi,point_chi)) 
			{
				if(point_chi->atboundary==0 && point_Al_hat->atboundary==0)
				{
					IGAElementNextPoint(elem_Al_hat,point_Al_hat);
					ierr = IGAPointGetWorkMat(point_chi,&Kpoint_chi);CHKERRQ(ierr);
					ierr = IGAPointGetWorkVec(point_chi,&Fpoint_chi);CHKERRQ(ierr);
					ierr = curlChiUM(point_chi,Kpoint_chi,NULL);CHKERRQ(ierr);
					ierr = curlChiUF(point_chi,point_Al_hat,Fpoint_chi,Al_chi,NULL);CHKERRQ(ierr);
					ierr = IGAPointAddMat(point_chi,Kpoint_chi,Kloc_chi);CHKERRQ(ierr);
					ierr = IGAPointAddVec(point_chi,Fpoint_chi,Floc_chi);CHKERRQ(ierr);
				}
			}
			IGAElementNextPoint(elem_Al_hat,point_Al_hat);

			ierr = IGAElementEndPoint(elem_chi,&point_chi);CHKERRQ(ierr);
			ierr = IGAElementEndPoint(elem_Al_hat,&point_Al_hat);CHKERRQ(ierr);
		}
		ierr = IGAElementFixSystem(elem_chi,Kloc_chi,Floc_chi);CHKERRQ(ierr);
		ierr = IGAElementAssembleMat(elem_chi,Kloc_chi,K_chi);CHKERRQ(ierr);
		ierr = IGAElementAssembleVec(elem_chi,Floc_chi,F_chi);CHKERRQ(ierr);

	}
	IGANextElement(iga_Al_hat,elem_Al_hat);
	ierr = IGAEndElement(iga_chi,&elem_chi);CHKERRQ(ierr);
	ierr = IGAEndElement(iga_Al_hat,&elem_Al_hat);CHKERRQ(ierr);

	// Restore local vectors S0 and arrays
	ierr = IGARestoreLocalVecArray(iga_Al_hat,al_hat,&local_Al_chi,&array_Al_chi);CHKERRQ(ierr);

	//Form system matrix and vector
	ierr = MatAssemblyBegin(K_chi,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyEnd  (K_chi,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	
	ierr = VecAssemblyBegin(F_chi);CHKERRQ(ierr);
	ierr = VecAssemblyEnd  (F_chi);CHKERRQ(ierr);

	KSP kspchiUp;
	ierr = IGACreateKSP(iga_chi,&kspchiUp);CHKERRQ(ierr);
	ierr = KSPSetOperators(kspchiUp,K_chi,K_chi);CHKERRQ(ierr);
	PC pcChiUp;
	ierr = KSPGetPC(kspchiUp,&pcChiUp); CHKERRQ(ierr);
	ierr = PCSetType(pcChiUp,PCLU); CHKERRQ(ierr);
	ierr = PCFactorSetMatSolverType(pcChiUp,MATSOLVERMUMPS); CHKERRQ(ierr);
	//ierr = KSPSetFromOptions(kspchiUp);CHKERRQ(ierr);
	//ierr = KSPSetTolerances(kspchiUp,1.0e-8,1.0e-20,PETSC_DEFAULT,PETSC_DEFAULT);CHKERRQ(ierr);
	ierr = KSPSolve(kspchiUp,F_chi,chi);CHKERRQ(ierr);

	char namechiUp[512];//="/ChiUp-2d-0.dat";
	sprintf(namechiUp,"%s%d%s","/ChiUp-2d-",i,".dat");
	char pathchiUp[1024];
	sprintf(pathchiUp,"%s%s",direct,namechiUp);
	ierr = IGAWriteVec(iga_chi,chi,pathchiUp);CHKERRQ(ierr);
//

//System for initial state of z0
	PetscPrintf(PETSC_COMM_WORLD,"\nSystem for Z0 starting \n\n");
	T=time(NULL);
	tm=*localtime(&T);
	PetscPrintf(PETSC_COMM_WORLD,"Current time is %02d:%02d:%02d \n",tm.tm_hour,tm.tm_min,tm.tm_sec);
	IGA iga_z;
	ierr = IGACreate(PETSC_COMM_WORLD,&iga_z);CHKERRQ(ierr);
	ierr = IGASetDim(iga_z,2);CHKERRQ(ierr);													//Spatial dimension of the problem
	ierr = IGASetDof(iga_z,2);CHKERRQ(ierr);													//Number of degrees of freedom, per node
	ierr = IGASetOrder(iga_z,1);CHKERRQ(ierr);													//Number of spatial derivatives to calculate
	ierr = IGASetFromOptions(iga_z);CHKERRQ(ierr);												//Note: The order (or degree) of the shape functions is given by the mesh!
	ierr = IGARead(iga_z,"./geometry.dat");CHKERRQ(ierr);
	
	for (dir=0; dir<2; dir++)
	{
		ierr = IGASetRuleType(iga_z,dir,IGA_RULE_LEGENDRE);CHKERRQ(ierr);
		ierr = IGASetRuleSize(iga_z,dir,6);CHKERRQ(ierr);
	}
	ierr = IGASetUp(iga_z);CHKERRQ(ierr);

	ierr = IGASetUp(iga_z);CHKERRQ(ierr);
	ierr = IGASetUp(iga_chi);CHKERRQ(ierr);


	for (dir=0; dir<2; dir++) 
	{
		for (side=0; side<2; side++) 
		{
			//ierr = IGASetBoundaryValue(iga,dir,side,dof,0.0);CHKERRQ(ierr);					// Dirichlet boundary conditions
			ierr = IGASetBoundaryForm(iga_z,dir,side,PETSC_TRUE);CHKERRQ(ierr);  				// Neumann boundary conditions
		}
	}

	Mat K_z;
	Vec Z0,F_z;

	ierr = IGACreateMat(iga_z,&K_z);CHKERRQ(ierr);
	ierr = IGACreateVec(iga_z,&Z0);CHKERRQ(ierr);
	ierr = IGACreateVec(iga_z,&F_z);CHKERRQ(ierr);

	IGAPoint		point_z;								//point
	IGAElement		elem_z;									//element
	PetscReal		*Kloc_z,*Floc_z;						//AA y BB
	PetscReal		*Kpoint_z,*Fpoint_z;					//KKK y FFF
	const PetscReal	*array_chi_z;							//arrayU
	Vec				local_chi_z;							//localU
	PetscReal		*Chi_z;								//U

  	IGAFormSystem	wtfZ0;
 	void			*wtf2Z0;

 	KSP kspZ0;
	ierr = IGACreateKSP(iga_z,&kspZ0);CHKERRQ(ierr);

	// Get local vectors Chi0  and arrays
	ierr = IGAGetLocalVecArray(iga_chi,chi,&local_chi_z,&array_chi_z);CHKERRQ(ierr);

	// Element loop
	ierr = IGABeginElement(iga_z,&elem_z);CHKERRQ(ierr);
	ierr = IGABeginElement(iga_chi,&elem_chi);CHKERRQ(ierr);

	while (IGANextElement(iga_z,elem_z)) 
	{
		IGANextElement(iga_chi,elem_chi);

		ierr = IGAElementGetWorkMat(elem_z,&Kloc_z);CHKERRQ(ierr);
		ierr = IGAElementGetWorkVec(elem_z,&Floc_z);CHKERRQ(ierr);
		ierr = IGAElementGetValues(elem_chi,array_chi_z,&Chi_z);CHKERRQ(ierr);

		// FormSystem loop
		while (IGAElementNextFormSystem(elem_z,&wtfZ0,&wtf2Z0)) 
		{
		// Quadrature loop
			ierr = IGAElementBeginPoint(elem_z,&point_z);CHKERRQ(ierr);
			ierr = IGAElementBeginPoint(elem_chi,&point_chi);CHKERRQ(ierr);

			while (IGAElementNextPoint(elem_z,point_z))
			{
				if(point_z->atboundary==0 && point_chi->atboundary==0)
				{
					IGAElementNextPoint(elem_chi,point_chi);

					ierr = IGAPointGetWorkMat(point_z,&Kpoint_z);CHKERRQ(ierr);
					ierr = IGAPointGetWorkVec(point_z,&Fpoint_z);CHKERRQ(ierr);
					//	   Z0sys(IGAPoint p, IGAPoint pChi, IGAPoint pS,IGAPoint pZ, PetscReal *K, PetscReal *F, PetscReal *UChi, PetscReal *S,PetscReal *ZS, void *ctx)
					ierr = Z0sys(point_z,point_chi,Kpoint_z,Fpoint_z,Chi_z,NULL);CHKERRQ(ierr);
					ierr = IGAPointAddMat(point_z,Kpoint_z,Kloc_z);CHKERRQ(ierr);
					ierr = IGAPointAddVec(point_z,Fpoint_z,Floc_z);CHKERRQ(ierr);
				}
			}
			ierr = IGAElementEndPoint(elem_z,&point_z);CHKERRQ(ierr);
			while (point_chi->index != -1)
			{
				IGAElementNextPoint(elem_chi,point_chi);
			}
			ierr = IGAElementEndPoint(elem_chi,&point_chi);CHKERRQ(ierr);
		}

		ierr = IGAElementAssembleMat(elem_z,Kloc_z,K_z);CHKERRQ(ierr);
		ierr = IGAElementAssembleVec(elem_z,Floc_z,F_z);CHKERRQ(ierr);

	}
	IGANextElement(iga_chi,elem_chi);

	ierr = IGAEndElement(iga_z,&elem_z);CHKERRQ(ierr);
	ierr = IGAEndElement(iga_chi,&elem_chi);CHKERRQ(ierr);

	// Restore local vectors Chi0 and arrays
	ierr = IGARestoreLocalVecArray(iga_chi,chi,&local_chi_z,&array_chi_z);CHKERRQ(ierr);

	ierr = MatAssemblyBegin(K_z,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyEnd  (K_z,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

	PetscInt n,m;
	PetscInt rows, *cols;
	PetscReal val, *vals;

	//Here we set values to the Matrix directly, to impose Dirichlet condition in a single point.
	//Note: Lower Left corner is gdl's  0 and 1,
	//		Lower Right corner is gdl's 2*(nx+2)-2 and 2*(nx+2)-1 
	//		Upper Left corner is gdl's  2*(nx+2)*(ny+2)-2*(nx+1)-2 and 2*(nx+2)*(ny+2)-2*(nx+1)-1
	//		Upper Right corner is gdl's 2*(nx+2)*(ny+2)-2 and 2*(nx+2)*(ny+2)-1
	//All of these for when z is a 2nd order nurbs
	
	ierr = MatGetSize(K_z,&n,&m);CHKERRQ(ierr);
	ierr = MatSetOption(K_z, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);CHKERRQ(ierr);

	ierr = PetscMalloc1(m,&cols);CHKERRQ(ierr);
	ierr = PetscMalloc1(m,&vals);CHKERRQ(ierr);

	for(int i=0;i<m;i++)
	{
		cols[i]=i;
		vals[i]=0.0;
	}

	rows=0;
	vals[rows]=1.0e6;
	ierr = MatSetValues(K_z,1,&rows,m,cols,vals,INSERT_VALUES);CHKERRQ(ierr);
	ierr = MatSetValues(K_z,n,cols,1,&rows,vals,INSERT_VALUES);CHKERRQ(ierr);
	vals[rows]=0.0;

	rows=1;
	vals[rows]=1.0e6;
	ierr = MatSetValues(K_z,1,&rows,m,cols,vals,INSERT_VALUES);CHKERRQ(ierr);
	ierr = MatSetValues(K_z,n,cols,1,&rows,vals,INSERT_VALUES);CHKERRQ(ierr);
	vals[rows]=0.0;

	rows=2*(nx+1)*(ny+1)-2;														//This is the dof in x on the upper right corner for 1st order elements
	//rows=2*(nx+2)*(ny+2)-2*(nx+1)-2; 											//This is for when z is a 2nd order nurb
	//rows=2*(nx+3)*(ny+3)-2; 													//This is for when z is a 3rd order nurb
	vals[rows]=1.0e6;
	ierr = MatSetValues(K_z,1,&rows,m,cols,vals,INSERT_VALUES);CHKERRQ(ierr);
	ierr = MatSetValues(K_z,n,cols,1,&rows,vals,INSERT_VALUES);CHKERRQ(ierr);
	vals[rows]=0.0;

	ierr = MatAssemblyBegin(K_z,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyEnd  (K_z,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

	//Here we set values to the vector directly (to impose Dirichlet condition in a single point)
	//Note: Lower Left corner is gdl's  0 and 1,
	//		Lower Right corner is gdl's 2*(nx+2)-2 and 2*(nx+2)-1 
	//		Upper Left corner is gdl's  2*(nx+2)*(ny+2)-2-2*(nx+1) and 2*(nx+2)*(ny+2)-1-2*(nx+1)
	//		Upper Right corner is gdl's 2*(nx+2)*(ny+2)-2 and 2*(nx+2)*(ny+2)-1
	//All of these for when z is a 2nd order nurbs
	ierr = VecAssemblyBegin(F_z);CHKERRQ(ierr);
	ierr = VecAssemblyEnd  (F_z);CHKERRQ(ierr);

	rows=0;
	val=0.0;
	ierr = VecSetValue(F_z,rows,val,INSERT_VALUES);CHKERRQ(ierr);

	ierr = VecAssemblyBegin(F_z);CHKERRQ(ierr);
	ierr = VecAssemblyEnd  (F_z);CHKERRQ(ierr);

	rows=1;
	val=0.0;
	ierr = VecSetValue(F_z,rows,val,INSERT_VALUES);CHKERRQ(ierr);

	ierr = VecAssemblyBegin(F_z);CHKERRQ(ierr);
	ierr = VecAssemblyEnd  (F_z);CHKERRQ(ierr);
	
	rows=2*(nx+1)*(ny+1)-1;
	//rows=2*(nx+2)*(ny+2)-2*(nx+1)-2; 										//This is for when z is a 2nd order nurb
	//rows=2*(nx+3)*(ny+3)-2;
	val=0.0;
	ierr = VecSetValue(F_z,rows,val,INSERT_VALUES);CHKERRQ(ierr);

	ierr = VecAssemblyBegin(F_z);CHKERRQ(ierr);
	ierr = VecAssemblyEnd  (F_z);CHKERRQ(ierr);

	ierr = KSPSetOperators(kspZ0,K_z,K_z);CHKERRQ(ierr);
	//PC pcZ0;
	//ierr = KSPGetPC(kspZ0,&pcZ0); CHKERRQ(ierr);
	//ierr = PCSetType(pcZ0,PCLU); CHKERRQ(ierr);
	//ierr = PCFactorSetMatSolverType(pcZ0,MATSOLVERMUMPS); CHKERRQ(ierr);
	ierr = KSPSetFromOptions(kspZ0);CHKERRQ(ierr);
	ierr = KSPSetTolerances(kspZ0,1e-18,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);CHKERRQ(ierr);
	ierr = KSPSolve(kspZ0,F_z,Z0);CHKERRQ(ierr);

	ierr = KSPDestroy(&kspZ0);CHKERRQ(ierr);
	ierr = MatDestroy(&K_z);CHKERRQ(ierr);
	ierr = VecDestroy(&F_z);CHKERRQ(ierr);
	//ierr = IGAWriteVec(iga_z,z0,"./results/z0-2d-0.dat");CHKERRQ(ierr);
	char nameZ0[512];//="/Z0-2d-0.dat";
	sprintf(nameZ0,"%s%d%s","/Z0-2d-",i,".dat");
	char pathZ0[1024];
	sprintf(pathZ0,"%s%s",direct,nameZ0);
	ierr = IGAWriteVec(iga_z,Z0,pathZ0);CHKERRQ(ierr);	
//

//System for initial state of u
	PetscPrintf(PETSC_COMM_WORLD,"\nSystem for U0 starting \n\n");
	T=time(NULL);
	tm=*localtime(&T);
	PetscPrintf(PETSC_COMM_WORLD,"Current time is %02d:%02d:%02d \n",tm.tm_hour,tm.tm_min,tm.tm_sec);
	IGA iga_u;
	ierr = IGACreate(PETSC_COMM_WORLD,&iga_u);CHKERRQ(ierr);
	ierr = IGASetDim(iga_u,2);CHKERRQ(ierr);													//Spatial dimension of the problem
	ierr = IGASetDof(iga_u,2);CHKERRQ(ierr);													//Number of degrees of freedom, per node
	ierr = IGASetOrder(iga_u,1);CHKERRQ(ierr);													//Number of spatial derivatives to calculate
	ierr = IGASetFromOptions(iga_u);CHKERRQ(ierr);												//Note: The order (or degree) of the shape functions is given by the mesh!
	ierr = IGARead(iga_u,"./geometry.dat");CHKERRQ(ierr);
	
	for (dir=0; dir<2; dir++)
	{
		ierr = IGASetRuleType(iga_u,dir,IGA_RULE_LEGENDRE);CHKERRQ(ierr);
		ierr = IGASetRuleSize(iga_u,dir,6);CHKERRQ(ierr);
	}
	ierr = IGASetUp(iga_u);CHKERRQ(ierr);

	ierr = IGASetUp(iga_u);CHKERRQ(ierr);

	for (dir=0; dir<2; dir++) 
	{
		for (side=0; side<2; side++) 
		{
			//ierr = IGASetBoundaryValue(iga,dir,side,dof,0.0);CHKERRQ(ierr);					// Dirichlet boundary conditions
			ierr = IGASetBoundaryForm(iga_u,dir,side,PETSC_TRUE);CHKERRQ(ierr);  				// Neumann boundary conditions
		}
	}

	//Set Dirichlet boundary conditions here
	//ierr = IGASetBoundaryValue(iga,dir,side,dof,value);CHKERRQ(ierr);					// Dirichlet boundary conditions
	//ierr = IGASetBoundaryValue(iga_z,0,0,0,0.0);CHKERRQ(ierr);	//Left side, 1st dof = 0
	//ierr = IGASetBoundaryValue(iga_z,0,0,1,0.0);CHKERRQ(ierr);	//Left side, 2nd dof = 0

	//ierr = IGASetBoundaryValue(iga_z,0,1,0,0.0);CHKERRQ(ierr);	//Right side, 1st dof=0
	//ierr = IGASetBoundaryValue(iga_z,0,1,1,0.0);CHKERRQ(ierr);	//Right side, 2nd dof=0

	ierr = IGASetBoundaryValue(iga_u,1,0,0,0.0);CHKERRQ(ierr);	//Bottom side, 1st dof=0
	ierr = IGASetBoundaryValue(iga_u,1,0,1,0.0);CHKERRQ(ierr);	//Bottom side, 2nd dof=0
	
	//ierr = IGASetBoundaryValue(iga_z,1,1,0,0.0);CHKERRQ(ierr);	//Top side, 1st dof=0
	//ierr = IGASetBoundaryValue(iga_z,1,1,1,0.0);CHKERRQ(ierr);	//Top side, 2nd dof=0

	Mat K_u;
	Vec U,F_u;

	ierr = IGACreateMat(iga_u,&K_u);CHKERRQ(ierr);
	ierr = IGACreateVec(iga_u,&U);CHKERRQ(ierr);
	ierr = IGACreateVec(iga_u,&F_u);CHKERRQ(ierr);

	IGAPoint		point_u;								//point
	IGAElement		elem_u;								//element
	PetscReal		*Kloc_u,*Floc_u;						//AA y BB
	PetscReal		*Kpoint_u,*Fpoint_u;					//KKK y FFF
	const PetscReal	*array_chi_u, *array_z_u;				//arrayU
	Vec				local_chi_u, local_z_u;				//localU
	PetscReal		*Chi_u, *z_u;						//U

	IGAFormSystem	wtfU0;
	void			*wtf2U0;

	KSP kspU;
	ierr = IGACreateKSP(iga_u,&kspU);CHKERRQ(ierr);

	// Get local vectors Chi0  and arrays
	ierr = IGAGetLocalVecArray(iga_chi,chi,&local_chi_u,&array_chi_u);CHKERRQ(ierr);
	ierr = IGAGetLocalVecArray(iga_z,chi,&local_chi_u,&array_chi_u);CHKERRQ(ierr);
	
	// Element loop
	ierr = IGABeginElement(iga_u,&elem_u);CHKERRQ(ierr);
	ierr = IGABeginElement(iga_chi,&elem_chi);CHKERRQ(ierr);
	ierr = IGABeginElement(iga_z,&elem_z);CHKERRQ(ierr);

	while (IGANextElement(iga_u,elem_u)) 
	{
		IGANextElement(iga_chi,elem_chi);
		IGANextElement(iga_z,elem_z);

		ierr = IGAElementGetWorkMat(elem_u,&Kloc_u);CHKERRQ(ierr);
		ierr = IGAElementGetWorkVec(elem_u,&Floc_u);CHKERRQ(ierr);
		ierr = IGAElementGetValues(elem_chi,array_chi_u,&Chi_u);CHKERRQ(ierr);
		ierr = IGAElementGetValues(elem_z,array_z_u,&z_u);CHKERRQ(ierr);

		// FormSystem loop
		while (IGAElementNextFormSystem(elem_u,&wtfU0,&wtf2U0)) 
		{
		// Quadrature loop
			ierr = IGAElementBeginPoint(elem_u,&point_u);CHKERRQ(ierr);
			ierr = IGAElementBeginPoint(elem_chi,&point_chi);CHKERRQ(ierr);
			ierr = IGAElementBeginPoint(elem_z,&point_z);CHKERRQ(ierr);

			while (IGAElementNextPoint(elem_u,point_u))
			{
				if(point_u->atboundary==1)
				{
					ierr = IGAPointGetWorkVec(point_u,&Fpoint_u);CHKERRQ(ierr);
					//	   UsysF(IGAPoint p,IGAPoint pChi,IGAPoint pZ,PetscReal *F,PetscReal *UChi,PetscReal *UZ,void *ctx)
					ierr = UsysF(point_u,point_chi,point_z,Fpoint_u,Chi_u,z_u,NULL);CHKERRQ(ierr);
					ierr = IGAPointAddVec(point_u,Fpoint_u,Floc_u);CHKERRQ(ierr);
				}
				if(point_u->atboundary==0)
				{
					IGAElementNextPoint(elem_chi,point_chi);
					IGAElementNextPoint(elem_z,point_z);

					ierr = IGAPointGetWorkMat(point_u,&Kpoint_u);CHKERRQ(ierr);
					ierr = IGAPointGetWorkVec(point_u,&Fpoint_u);CHKERRQ(ierr);
					ierr = UsysF(point_u,point_chi,point_z,Fpoint_u,Chi_u,z_u,NULL);CHKERRQ(ierr);
					ierr = UsysM(point_u,Kpoint_u,NULL);CHKERRQ(ierr);
					ierr = IGAPointAddMat(point_u,Kpoint_u,Kloc_u);CHKERRQ(ierr);
					ierr = IGAPointAddVec(point_u,Fpoint_u,Floc_u);CHKERRQ(ierr);
				}
			}
			while (point_chi->index != -1)
			{
				IGAElementNextPoint(elem_chi,point_chi);
			}
			while (point_z->index != -1)
			{
				IGAElementNextPoint(elem_z,point_z);
			}

			ierr = IGAElementEndPoint(elem_u,&point_u);CHKERRQ(ierr);
			ierr = IGAElementEndPoint(elem_chi,&point_chi);CHKERRQ(ierr);
			ierr = IGAElementEndPoint(elem_z,&point_z);CHKERRQ(ierr);
		}
		ierr = IGAElementFixSystem(elem_u,Kloc_u,Floc_u);CHKERRQ(ierr);					//This sets Dirichlet condition ¿? (Yes, this applies the conditions from IGASetBoundaryValue)
		ierr = IGAElementAssembleMat(elem_u,Kloc_u,K_u);CHKERRQ(ierr);
		ierr = IGAElementAssembleVec(elem_u,Floc_u,F_u);CHKERRQ(ierr);
	}
	IGANextElement(iga_chi,elem_chi);
	IGANextElement(iga_u,elem_u);

	ierr = IGAEndElement(iga_u,&elem_u);CHKERRQ(ierr);
	ierr = IGAEndElement(iga_chi,&elem_chi);CHKERRQ(ierr);
	ierr = IGAEndElement(iga_z,&elem_z);CHKERRQ(ierr);

	ierr = MatAssemblyBegin(K_u,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyEnd  (K_u,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

	ierr = VecAssemblyBegin(F_u);CHKERRQ(ierr);
	ierr = VecAssemblyEnd  (F_u);CHKERRQ(ierr);

	ierr = KSPSetOperators(kspU,K_u,K_u);CHKERRQ(ierr);
	//PC pcU;
	//ierr = KSPGetPC(kspU,&pcU); CHKERRQ(ierr);
	//ierr = PCSetType(pcU,PCLU); CHKERRQ(ierr);
	//ierr = PCFactorSetMatSolverType(pcU,MATSOLVERMUMPS); CHKERRQ(ierr);
	ierr = KSPSetFromOptions(kspU);CHKERRQ(ierr);
	ierr = KSPSetTolerances(kspU,1e-18,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);CHKERRQ(ierr);
	ierr = KSPSolve(kspU,F_u,U);CHKERRQ(ierr);

	char nameU[512];
	sprintf(nameU,"%s%d%s","/U-2d-",i,".dat");
	char pathU[1024];
	sprintf(pathU,"%s%s",direct,nameU);
	ierr = IGAWriteVec(iga_u,U,pathU);CHKERRQ(ierr);
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
	ierr = IGARead(igaVs,"./geometry.dat");CHKERRQ(ierr);
	
	for (dir=0; dir<2; dir++)
	{
		ierr = IGASetRuleType(igaVs,dir,IGA_RULE_LEGENDRE);CHKERRQ(ierr);
		ierr = IGASetRuleSize(igaVs,dir,6);CHKERRQ(ierr);
	}
	ierr = IGASetUp(igaVs);CHKERRQ(ierr);

	//for (dir=0; dir<2; dir++) 
	//{
		//for (side=0; side<2; side++) 
		//{
			//ierr = IGASetBoundaryValue(iga,dir,side,dof,0.0);CHKERRQ(ierr);					// Dirichlet boundary conditions
			//ierr = IGASetBoundaryForm(igaVs,dir,side,PETSC_TRUE);CHKERRQ(ierr);  				// Neumann boundary conditions
		//}
	//}

	Mat KVs;
	Vec Vs0,FVs;

	ierr = IGACreateMat(igaVs,&KVs);CHKERRQ(ierr);
	ierr = IGACreateVec(igaVs,&Vs0);CHKERRQ(ierr);
	ierr = IGACreateVec(igaVs,&FVs);CHKERRQ(ierr);

	IGAPoint		pointVs;													//point
	IGAElement		elemVs;														//element
	PetscReal		*KlocVs,*FlocVs;											//AA y BB
	PetscReal		*KpointVs,*FpointVs;										//KKK y FFF
	const PetscReal	*arrayZ0Vs,*arrayChi0Vs,*arrayS0Vs,*arrayUVs;				//arrayU
	Vec				localZ0Vs,localChi0Vs,localS0Vs,localUVs;					//localU
	PetscReal		*Chi0Vs,*Z0Vs,*S0Vs,*UVs;									//U

  	IGAFormSystem	wtfVs;
 	void			*wtf2Vs;

 	KSP kspVs;
	ierr = IGACreateKSP(igaVs,&kspVs);CHKERRQ(ierr);

	// Get local vectors Z0 and Chi0 and arrays
	ierr = IGAGetLocalVecArray(iga_chi,chi,&localChi0Vs,&arrayChi0Vs);CHKERRQ(ierr);
	ierr = IGAGetLocalVecArray(iga_z,Z0,&localZ0Vs,&arrayZ0Vs);CHKERRQ(ierr);
	ierr = IGAGetLocalVecArray(iga_u,U,&localUVs,&arrayUVs);CHKERRQ(ierr);
	ierr = IGAGetLocalVecArray(iga_S,S0,&localS0Vs,&arrayS0Vs);CHKERRQ(ierr);

	// Element loop
	ierr = IGABeginElement(igaVs,&elemVs);CHKERRQ(ierr);
	ierr = IGABeginElement(iga_z,&elem_z);CHKERRQ(ierr);
	ierr = IGABeginElement(iga_chi,&elem_chi);CHKERRQ(ierr);
	ierr = IGABeginElement(iga_u,&elem_u);CHKERRQ(ierr);
	ierr = IGABeginElement(iga_S,&elem_S);CHKERRQ(ierr);

	while (IGANextElement(igaVs,elemVs)) 
	{
		IGANextElement(iga_z,elem_z);
		IGANextElement(iga_chi,elem_chi);
		IGANextElement(iga_u,elem_u);
		IGANextElement(iga_S,elem_S);

		ierr = IGAElementGetWorkMat(elemVs,&KlocVs);CHKERRQ(ierr);
		ierr = IGAElementGetWorkVec(elemVs,&FlocVs);CHKERRQ(ierr);
		ierr = IGAElementGetValues(elem_z,arrayZ0Vs,&Z0Vs);CHKERRQ(ierr);
		ierr = IGAElementGetValues(elem_chi,arrayChi0Vs,&Chi0Vs);CHKERRQ(ierr);
		ierr = IGAElementGetValues(elem_u,arrayUVs,&UVs);CHKERRQ(ierr);
		ierr = IGAElementGetValues(elem_S,arrayS0Vs,&S0Vs);CHKERRQ(ierr);

		// FormSystem loop
		while (IGAElementNextFormSystem(elemVs,&wtfVs,&wtf2Vs)) 
		{
		// Quadrature loop
			ierr = IGAElementBeginPoint(elemVs,&pointVs);CHKERRQ(ierr);
			ierr = IGAElementBeginPoint(elem_z,&point_z);CHKERRQ(ierr);
			ierr = IGAElementBeginPoint(elem_chi,&point_chi);CHKERRQ(ierr);
			ierr = IGAElementBeginPoint(elem_u,&point_u);CHKERRQ(ierr);
			ierr = IGAElementBeginPoint(elem_S,&point_S);CHKERRQ(ierr);

			while (IGAElementNextPoint(elemVs,pointVs))
			{
				if(pointVs->atboundary==1)
				{

				}
				if(pointVs->atboundary==0 && point_z->atboundary==0 && point_chi->atboundary==0 && point_S->atboundary==0 && point_u->atboundary==0)
				{
					IGAElementNextPoint(elem_z,point_z);
					IGAElementNextPoint(elem_chi,point_chi);
					IGAElementNextPoint(elem_S,point_S);

					ierr = IGAPointGetWorkMat(pointVs,&KpointVs);CHKERRQ(ierr);
					ierr = IGAPointGetWorkVec(pointVs,&FpointVs);CHKERRQ(ierr);
					//VS(IGAPoint p,IGAPoint pChi,IGAPoint pZu,IGAPoint pU,IGAPoint pS,PetscReal *K,PetscReal *F,PetscReal *UChi,PetscReal *UZu,PetscReal *Uu,PetscReal *US,void *ctx)
					ierr = VS(pointVs,point_chi,point_z,point_u,point_S,KpointVs,FpointVs,Chi0Vs,Z0Vs,UVs,S0Vs,NULL);CHKERRQ(ierr);
					ierr = IGAPointAddMat(pointVs,KpointVs,KlocVs);CHKERRQ(ierr);
					ierr = IGAPointAddVec(pointVs,FpointVs,FlocVs);CHKERRQ(ierr);
				}
				else
				{
					PetscPrintf(PETSC_COMM_SELF,"Problema en Vs \n");CHKERRQ(ierr);
				}
			}
			while (point_z->index != -1)
			{
				IGAElementNextPoint(elem_z,point_z);
			}
			while (point_chi->index != -1)
			{
				IGAElementNextPoint(elem_chi,point_chi);
			}
			while (point_u->index != -1)
			{
				IGAElementNextPoint(elem_u,point_u);
			}
			while (point_S->index != -1)
			{
				IGAElementNextPoint(elem_S,point_S);
			}
			ierr = IGAElementEndPoint(elemVs,&pointVs);CHKERRQ(ierr);
			ierr = IGAElementEndPoint(elem_z,&point_z);CHKERRQ(ierr);
			ierr = IGAElementEndPoint(elem_chi,&point_chi);CHKERRQ(ierr);
			ierr = IGAElementEndPoint(elem_u,&point_u);CHKERRQ(ierr);
			ierr = IGAElementEndPoint(elem_S,&point_S);CHKERRQ(ierr);
		}

		//ierr = IGAElementFixSystem(elem_Stress,KlocStress,FlocStress);CHKERRQ(ierr);					//This sets Dirichlet condition ¿? (Yes, this applies the conditions from IGASetBoundaryValue)
		ierr = IGAElementAssembleMat(elemVs,KlocVs,KVs);CHKERRQ(ierr);
		ierr = IGAElementAssembleVec(elemVs,FlocVs,FVs);CHKERRQ(ierr);
	}
	IGANextElement(iga_z,elem_z);
	IGANextElement(iga_chi,elem_chi);
	IGANextElement(iga_u,elem_u);
	IGANextElement(iga_S,elem_S);

	ierr = IGAEndElement(igaVs,&elemVs);CHKERRQ(ierr);
	ierr = IGAEndElement(iga_z,&elem_z);CHKERRQ(ierr);
	ierr = IGAEndElement(iga_chi,&elem_chi);CHKERRQ(ierr);
	ierr = IGAEndElement(iga_u,&elem_u);CHKERRQ(ierr);
	ierr = IGAEndElement(iga_S,&elem_S);CHKERRQ(ierr);

	// Restore local vectors u, Z0, Chi0 and arrays
	ierr = IGARestoreLocalVecArray(iga_z,Z0,&localZ0Vs,&arrayZ0Vs);CHKERRQ(ierr);
	ierr = IGARestoreLocalVecArray(iga_chi,chi,&localChi0Vs,&arrayChi0Vs);CHKERRQ(ierr);
	ierr = IGARestoreLocalVecArray(iga_u,U,&localUVs,&arrayUVs);CHKERRQ(ierr);
	ierr = IGARestoreLocalVecArray(iga_S,S0,&localS0Vs,&arrayS0Vs);CHKERRQ(ierr);
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
	ierr = KSPSetTolerances(kspVs,1.0e-24,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);CHKERRQ(ierr);
	ierr = KSPSolve(kspVs,FVs,Vs0);CHKERRQ(ierr);

	//VecChop(Vec v, PetscReal tol) Sets anything with an absolute value less than the tolerance to 0
	ierr = VecChop(Vs0,1e-11);CHKERRQ(ierr);

	char nameVs[512];//="/Vs-2d-0.dat";
	sprintf(nameVs,"%s%d%s","/Vs-2d-",i,".dat");
	char pathVs[1024];
	sprintf(pathVs,"%s%s",direct,nameVs);
	ierr = IGAWriteVec(igaVs,Vs0,pathVs);CHKERRQ(ierr);	
//

//System for L2 protection of smoothed V^{S}
	PetscPrintf(PETSC_COMM_WORLD,"\nSystem for Smooth Vs starting \n\n");
	T=time(NULL);
	tm=*localtime(&T);
	PetscPrintf(PETSC_COMM_WORLD,"Current time is %02d:%02d:%02d \n",tm.tm_hour,tm.tm_min,tm.tm_sec);

	//First integrate velocity for both defects
	Vec FVsInt, FVsXi, Int1sVec, Int2sVec;
	ierr = IGACreateVec(igaVs,&FVsInt);CHKERRQ(ierr);
	ierr = IGACreateVec(igaVs,&FVsXi);CHKERRQ(ierr);
	ierr = IGACreateVec(igaVs,&Int1sVec);CHKERRQ(ierr);
	ierr = IGACreateVec(igaVs,&Int2sVec);CHKERRQ(ierr);
																														//Fix these comments
	PetscReal		*PointInt1a,*PointInt2a,*PointInt,*locInt1a,*locInt2a;												//AA y BB
	const PetscReal	*arrayVs0,*arraySVs;																				//arrayU
	Vec				localVs0,localSVs;																					//localU
	PetscReal 		*SVs,*Vs0Values,Int1a=0.0,Int2a=0.0;

	IGAFormSystem	wtfVsInt;
	void			*wtf2VsInt;

	// Get local vectors V0 and InputAl and arrays
	ierr = IGAGetLocalVecArray(igaVs,Vs0,&localVs0,&arrayVs0);CHKERRQ(ierr);
	ierr = IGAGetLocalVecArray(iga_S,S0,&localSVs,&arraySVs);CHKERRQ(ierr);

	// Element loop
	ierr = IGABeginElement(igaVs,&elemVs);CHKERRQ(ierr);
	ierr = IGABeginElement(iga_S,&elem_S);CHKERRQ(ierr);

	while (IGANextElement(igaVs,elemVs))
	{
		IGANextElement(iga_S,elem_S);

		ierr = IGAElementGetWorkVec(elemVs,&locInt1a);CHKERRQ(ierr);
		ierr = IGAElementGetWorkVec(elemVs,&locInt2a);CHKERRQ(ierr);
		ierr = IGAElementGetValues(elemVs,arrayVs0,&Vs0Values);CHKERRQ(ierr);
		ierr = IGAElementGetValues(elem_S,arraySVs,&SVs);CHKERRQ(ierr);

		// FormSystem loop
		while (IGAElementNextFormSystem(elemVs,&wtfVsInt,&wtf2VsInt))
		{
		// Quadrature loop
			ierr = IGAElementBeginPoint(elemVs,&pointVs);CHKERRQ(ierr);
			ierr = IGAElementBeginPoint(elem_S,&point_S);CHKERRQ(ierr);
			
			while (IGAElementNextPoint(elemVs,pointVs))
			{
				if(pointVs->atboundary==1)
				{

				}
				if(pointVs->atboundary==0 && point_S->atboundary==0)
				{
					IGAElementNextPoint(elem_S,point_S);
					ierr = IGAPointGetWorkVec(pointVs,&PointInt1a);CHKERRQ(ierr);
					ierr = IGAPointGetWorkVec(pointVs,&PointInt2a);CHKERRQ(ierr);
					//     Int_Xi_Vs(IGAPoint pV,IGAPoint pS,PetscReal *FInt1a,PetscReal *FInt2a,PetscReal *VS,PetscReal *US,void *ctx)
					ierr = Int_Xi_Vs(pointVs,point_S,PointInt1a,PointInt2a,Vs0Values,SVs,NULL);CHKERRQ(ierr);
					ierr = IGAPointAddVec(pointVs,PointInt1a,locInt1a);CHKERRQ(ierr);
					ierr = IGAPointAddVec(pointVs,PointInt2a,locInt2a);CHKERRQ(ierr);
				}
			}
			while (point_S->index != -1)
			{
				IGAElementNextPoint(elem_S,point_S);
			}
			ierr = IGAElementEndPoint(elemVs,&pointVs);CHKERRQ(ierr);
			ierr = IGAElementEndPoint(elem_S,&point_S);CHKERRQ(ierr);
		}
		ierr = IGAElementAssembleVec(elemVs,locInt1a,Int1sVec);CHKERRQ(ierr);
		ierr = IGAElementAssembleVec(elemVs,locInt2a,Int2sVec);CHKERRQ(ierr);
	}
	IGANextElement(iga_S,elem_S);

	ierr = IGAEndElement(igaVs,&elemVs);CHKERRQ(ierr);
	ierr = IGAEndElement(iga_S,&elem_S);CHKERRQ(ierr);

	ierr = IGARestoreLocalVecArray(iga_S,S0,&localSVs,&arraySVs);CHKERRQ(ierr);
	ierr = IGARestoreLocalVecArray(igaVs,Vs0,&localVs0,&arrayVs0);CHKERRQ(ierr);

	ierr = VecAssemblyBegin(Int1sVec);CHKERRQ(ierr);
	ierr = VecAssemblyEnd  (Int1sVec);CHKERRQ(ierr);
	ierr = VecAssemblyBegin(Int2sVec);CHKERRQ(ierr);
	ierr = VecAssemblyEnd  (Int2sVec);CHKERRQ(ierr);

	//Here add all values at Gauss point_S
	ierr = VecSum(Int1sVec,&Int1a);CHKERRQ(ierr);
	ierr = VecSum(Int2sVec,&Int2a);CHKERRQ(ierr);

	ierr = PetscPrintf(PETSC_COMM_WORLD,"Int1a=%f, Int2a=%f \n",Int1a,Int2a);CHKERRQ(ierr);
	//Int1a and Int2a are the integrals of Vs(1) and Vs(2) ovet the domain, respectively.	
//

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

	IGAPoint		pointVa;														//point
	IGAElement		elemVa;															//element
	PetscReal		*KlocVa,*FlocVa;												//AA y BB
	PetscReal		*KpointVa,*FpointVa;											//KKK y FFF
	const PetscReal	*arrayZ0Va,*arrayChi0Va,*arrayAl0Va,*arraySVa,*arrayUVa;		//arrayU
	Vec				localZ0Va,localChi0Va,localAl0Va,localSVa,localUVa;				//localU
	PetscReal		*Chi0Va,*Z0Va,*Al0Va,*SVa,*UVa;									//U

	IGAFormSystem	wtfVa;
	void			*wtf2Va;

 	KSP kspVa;
	ierr = IGACreateKSP(igaVa,&kspVa);CHKERRQ(ierr);

	// Get local vectors Z0 and Chi0 and arrays
	ierr = IGAGetLocalVecArray(iga_chi,chi,&localChi0Va,&arrayChi0Va);CHKERRQ(ierr);
	ierr = IGAGetLocalVecArray(iga_z,Z0,&localZ0Va,&arrayZ0Va);CHKERRQ(ierr);
	ierr = IGAGetLocalVecArray(iga_u,U,&localUVa,&arrayUVa);CHKERRQ(ierr);
	ierr = IGAGetLocalVecArray(iga_Al_hat,alInput,&localAl0Va,&arrayAl0Va);CHKERRQ(ierr);
	ierr = IGAGetLocalVecArray(iga_S,S0,&localSVa,&arraySVa);CHKERRQ(ierr);

	// Element loop
	ierr = IGABeginElement(igaVa,&elemVa);CHKERRQ(ierr);
	ierr = IGABeginElement(iga_z,&elem_z);CHKERRQ(ierr);
	ierr = IGABeginElement(iga_chi,&elem_chi);CHKERRQ(ierr);
	ierr = IGABeginElement(iga_u,&elem_u);CHKERRQ(ierr);
	ierr = IGABeginElement(iga_Al_hat,&elem_Al_hat);CHKERRQ(ierr);
	ierr = IGABeginElement(iga_S,&elem_S);CHKERRQ(ierr);

	while (IGANextElement(igaVa,elemVa)) 
	{
		IGANextElement(iga_z,elem_z);
		IGANextElement(iga_chi,elem_chi);
		IGANextElement(iga_u,elem_u);
		IGANextElement(iga_Al_hat,elem_Al_hat);
		IGANextElement(iga_S,elem_S);

		ierr = IGAElementGetWorkMat(elemVa,&KlocVa);CHKERRQ(ierr);
		ierr = IGAElementGetWorkVec(elemVa,&FlocVa);CHKERRQ(ierr);
		ierr = IGAElementGetValues(elem_z,arrayZ0Va,&Z0Va);CHKERRQ(ierr);
		ierr = IGAElementGetValues(elem_chi,arrayChi0Va,&Chi0Va);CHKERRQ(ierr);
		ierr = IGAElementGetValues(elem_u,arrayUVa,&UVa);CHKERRQ(ierr);
		ierr = IGAElementGetValues(elem_Al_hat,arrayAl0Va,&Al0Va);CHKERRQ(ierr);
		ierr = IGAElementGetValues(elem_S,arraySVa,&SVa);CHKERRQ(ierr);

		// FormSystem loop
		while (IGAElementNextFormSystem(elemVa,&wtfVa,&wtf2Va)) 
		{
		// Quadrature loop
			ierr = IGAElementBeginPoint(elemVa,&pointVa);CHKERRQ(ierr);
			ierr = IGAElementBeginPoint(elem_z,&point_z);CHKERRQ(ierr);
			ierr = IGAElementBeginPoint(elem_chi,&point_chi);CHKERRQ(ierr);
			ierr = IGAElementBeginPoint(elem_u,&point_u);CHKERRQ(ierr);
			ierr = IGAElementBeginPoint(elem_Al_hat,&point_Al_hat);CHKERRQ(ierr);
			ierr = IGAElementBeginPoint(elem_S,&point_S);CHKERRQ(ierr);

			while (IGAElementNextPoint(elemVa,pointVa))
			{
				if(pointVa->atboundary==1)
				{

				}
				if(pointVa->atboundary==0 && point_z->atboundary==0 && point_chi->atboundary==0 && point_Al_hat->atboundary==0 && point_S->atboundary==0)
				{
					IGAElementNextPoint(elem_z,point_z);
					IGAElementNextPoint(elem_chi,point_chi);
					IGAElementNextPoint(elem_u,point_u);
					IGAElementNextPoint(elem_Al_hat,point_Al_hat);
					IGAElementNextPoint(elem_S,point_S);

					ierr = IGAPointGetWorkMat(pointVa,&KpointVa);CHKERRQ(ierr);
					ierr = IGAPointGetWorkVec(pointVa,&FpointVa);CHKERRQ(ierr);
					//Valpha(IGAPoint p,IGAPoint pChi,IGAPoint pZu,IGAPoint pU,IGAPoint pAl,IGAPoint pS,PetscReal *K,PetscReal *F,PetscReal *Chi,PetscReal *Zu,PetscReal *Uu,PetscReal *UAl,PetscReal *US,void *ctx)
					ierr = Valpha(pointVa,point_chi,point_z,point_u,point_Al_hat,point_S,KpointVa,FpointVa,Chi0Va,Z0Va,UVa,Al0Va,SVa,NULL);CHKERRQ(ierr);
					ierr = IGAPointAddMat(pointVa,KpointVa,KlocVa);CHKERRQ(ierr);
					ierr = IGAPointAddVec(pointVa,FpointVa,FlocVa);CHKERRQ(ierr);
				}
			}
			while (point_z->index != -1)
			{
				IGAElementNextPoint(elem_z,point_z);
			}
			while (point_chi->index != -1)
			{
				IGAElementNextPoint(elem_chi,point_chi);
			}
			while (point_u->index != -1)
			{
				IGAElementNextPoint(elem_u,point_u);
			}
			while (point_Al_hat->index != -1)
			{
				IGAElementNextPoint(elem_Al_hat,point_Al_hat);
			}
			while (point_S->index != -1)
			{
				IGAElementNextPoint(elem_S,point_S);
			}
			ierr = IGAElementEndPoint(elemVa,&pointVa);CHKERRQ(ierr);
			ierr = IGAElementEndPoint(elem_z,&point_z);CHKERRQ(ierr);
			ierr = IGAElementEndPoint(elem_chi,&point_chi);CHKERRQ(ierr);
			ierr = IGAElementEndPoint(elem_u,&point_u);CHKERRQ(ierr);
			ierr = IGAElementEndPoint(elem_Al_hat,&point_Al_hat);CHKERRQ(ierr);
			ierr = IGAElementEndPoint(elem_S,&point_S);CHKERRQ(ierr);
		}

		//ierr = IGAElementFixSystem(elem_Stress,KlocStress,FlocStress);CHKERRQ(ierr);					//This sets Dirichlet condition ¿? (Yes, this applies the conditions from IGASetBoundaryValue)
		ierr = IGAElementAssembleMat(elemVa,KlocVa,KVa);CHKERRQ(ierr);
		ierr = IGAElementAssembleVec(elemVa,FlocVa,FVa);CHKERRQ(ierr);
	}
	IGANextElement(iga_z,elem_z);
	IGANextElement(iga_chi,elem_chi);
	IGANextElement(iga_u,elem_u);
	IGANextElement(iga_Al_hat,elem_Al_hat);
	IGANextElement(iga_S,elem_S);

	ierr = IGAEndElement(igaVa,&elemVa);CHKERRQ(ierr);
	ierr = IGAEndElement(iga_z,&elem_z);CHKERRQ(ierr);
	ierr = IGAEndElement(iga_chi,&elem_chi);CHKERRQ(ierr);
	ierr = IGAEndElement(iga_u,&elem_u);
	ierr = IGAEndElement(iga_Al_hat,&elem_Al_hat);CHKERRQ(ierr);
	ierr = IGAEndElement(iga_S,&elem_S);CHKERRQ(ierr);

	// Restore local vectors u, Z0, Chi0 and arrays
	ierr = IGARestoreLocalVecArray(iga_z,Z0,&localZ0Va,&arrayZ0Va);CHKERRQ(ierr);
	ierr = IGARestoreLocalVecArray(iga_chi,chi,&localChi0Va,&arrayChi0Va);CHKERRQ(ierr);
	ierr = IGARestoreLocalVecArray(iga_u,U,&localUVa,&arrayUVa);CHKERRQ(ierr);
	ierr = IGARestoreLocalVecArray(iga_Al_hat,alInput,&localAl0Va,&arrayAl0Va);CHKERRQ(ierr);
	ierr = IGARestoreLocalVecArray(iga_S,S0,&localSVa,&arraySVa);CHKERRQ(ierr);
	ierr = MatAssemblyBegin(KVa,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyEnd  (KVa,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

	ierr = VecAssemblyBegin(FVa);CHKERRQ(ierr);
	ierr = VecAssemblyEnd  (FVa);CHKERRQ(ierr);

	ierr = KSPSetOperators(kspVa,KVa,KVa);CHKERRQ(ierr);
	//PC pcVa;
	//ierr = KSPGetPC(kspStress,&pcStress); CHKERRQ(ierr);
	//ierr = PCSetType(pcStress,PCLU); CHKERRQ(ierr);
	//ierr = PCFactorSetMatSolverType(pcStress,MATSOLVERMUMPS); CHKERRQ(ierr);
	ierr = KSPSetFromOptions(kspVa);CHKERRQ(ierr);
	ierr = KSPSetTolerances(kspVa,1.0e-24,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);CHKERRQ(ierr);
	ierr = KSPSolve(kspVa,FVa,Va0);CHKERRQ(ierr);

	//ierr = KSPDestroy(&kspVa);CHKERRQ(ierr);
	//ierr = MatDestroy(&KVa);CHKERRQ(ierr);
	//ierr = VecDestroy(&FVa);CHKERRQ(ierr);

	//VecChop(Vec v, PetscReal tol) Sets anything with an absolute value less than the tolerance to 0
	ierr = VecChop(Va0,1e-11);CHKERRQ(ierr);

	char nameVa[512];
	sprintf(nameVa,"%s%d%s","/Va-2d-",i,".dat");
	char pathVa[1024];
	sprintf(pathVa,"%s%s",direct,nameVa);
	ierr = IGAWriteVec(igaVa,Va0,pathVa);CHKERRQ(ierr);
//

//System for evolution of S (debug and ask \dot{S}*n=0 in boundaries where Vs*n=0 or other boundary conditions)
	PetscPrintf(PETSC_COMM_WORLD,"\nSystem for S_dot starting \n\n");
	T=time(NULL);
	tm=*localtime(&T);
	PetscPrintf(PETSC_COMM_WORLD,"Current time is %02d:%02d:%02d \n",tm.tm_hour,tm.tm_min,tm.tm_sec);
	IGA igaSdot;
	ierr = IGACreate(PETSC_COMM_WORLD,&igaSdot);CHKERRQ(ierr);
	ierr = IGASetDim(igaSdot,2);CHKERRQ(ierr);													//Spatial dimension of the problem
	ierr = IGASetDof(igaSdot,8);CHKERRQ(ierr);													//Number of degrees of freedom, per node
	ierr = IGASetOrder(igaSdot,1);CHKERRQ(ierr);													//Number of spatial derivatives to calculate
	ierr = IGASetFromOptions(igaSdot);CHKERRQ(ierr);												//Note: The order (or degree) of the shape functions is given by the mesh!
	ierr = IGARead(igaSdot,"./geometry.dat");CHKERRQ(ierr);
	
	for (dir=0; dir<2; dir++)
	{
		ierr = IGASetRuleType(igaSdot,dir,IGA_RULE_LEGENDRE);CHKERRQ(ierr);
		ierr = IGASetRuleSize(igaSdot,dir,6);CHKERRQ(ierr);
	}

	for (dir=0; dir<2; dir++) 
	{
		for (side=0; side<2; side++) 
		{
			//ierr = IGASetBoundaryValue(iga,dir,side,dof,0.0);CHKERRQ(ierr);					// Dirichlet boundary conditions
			ierr = IGASetBoundaryForm(igaSdot,dir,side,PETSC_TRUE);CHKERRQ(ierr);  				// Neumann boundary conditions
		}
	}

	for (dir=0; dir<2; dir++) 
	{
		for (side=0; side<2; side++) 
		{
			if(dir==0 && side==0)
			{
				//ierr = IGASetBoundaryValue(igaSdot,dir,side,0,0.0);CHKERRQ(ierr);    				// Dirichlet boundary conditions
				//ierr = IGASetBoundaryValue(igaSdot,dir,side,2,0.0);CHKERRQ(ierr);    				// Dirichlet boundary conditions
				//ierr = IGASetBoundaryValue(igaSdot,dir,side,4,0.0);CHKERRQ(ierr);    				// Dirichlet boundary conditions
				//ierr = IGASetBoundaryValue(igaSdot,dir,side,6,0.0);CHKERRQ(ierr);    				// Dirichlet boundary conditions
			}
			if(dir==0 && side==1)
			{
				//ierr = IGASetBoundaryValue(igaSdot,dir,side,0,0.0);CHKERRQ(ierr);    				// Dirichlet boundary conditions
				//ierr = IGASetBoundaryValue(igaSdot,dir,side,2,0.0);CHKERRQ(ierr);    				// Dirichlet boundary conditions
				//ierr = IGASetBoundaryValue(igaSdot,dir,side,4,0.0);CHKERRQ(ierr);    				// Dirichlet boundary conditions
				//ierr = IGASetBoundaryValue(igaSdot,dir,side,6,0.0);CHKERRQ(ierr);    				// Dirichlet boundary conditions
			}
			//ierr = IGASetBoundaryValue(iga,dir,side,dof,0.0);CHKERRQ(ierr);    					// Dirichlet boundary conditions
			//ierr = IGASetBoundaryForm(iga_Sp,dir,side,PETSC_TRUE);CHKERRQ(ierr);  				// Neumann boundary conditions
		}
	}
	ierr = IGASetUp(igaSdot);CHKERRQ(ierr);

	Mat KSdot;
	Vec Sdot,FSdot;

	ierr = IGACreateMat(igaSdot,&KSdot);CHKERRQ(ierr);
	ierr = IGACreateVec(igaSdot,&Sdot);CHKERRQ(ierr);
	ierr = IGACreateVec(igaSdot,&FSdot);CHKERRQ(ierr);

	IGAPoint		point_Sdot;													//point
	IGAElement		elem_Sdot;													//element
	PetscReal		*KlocSdot,*FlocSdot;										//AA y BB
	PetscReal		*Kpoint_Sdot,*Fpoint_Sdot;									//KKK y FFF
	//const PetscReal	*arrayVsSdot,*arrayS0Sdot;									//arrayU
	//Vec				localVsSdot,localS0Sdot;									//localU
	//PetscReal		*VsSdot,*S0Sdot;											//U
	const PetscReal	*arrayS0Sdot;									//arrayU
	Vec				localS0Sdot;									//localU
	PetscReal		*S0Sdot;											//U

  	IGAFormSystem	wtfSdot;
 	void			*wtf2Sdot;

 	KSP kspSdot;
	ierr = IGACreateKSP(igaSdot,&kspSdot);CHKERRQ(ierr);

	// Get local vectors Z0 and Chi0 and arrays
	//ierr = IGAGetLocalVecArray(igaVs,Vs0,&localVsSdot,&arrayVsSdot);CHKERRQ(ierr); Not anymore, use integral values
	ierr = IGAGetLocalVecArray(iga_S,S0,&localS0Sdot,&arrayS0Sdot);CHKERRQ(ierr);

	// Element loop
	//ierr = IGABeginElement(igaVs,&elemVs);CHKERRQ(ierr);
	ierr = IGABeginElement(iga_S,&elem_S);CHKERRQ(ierr);
	ierr = IGABeginElement(igaSdot,&elem_Sdot);CHKERRQ(ierr);

	while (IGANextElement(igaSdot,elem_Sdot)) 
	{
		//IGANextElement(igaVs,elemVs);
		IGANextElement(iga_S,elem_S);

		ierr = IGAElementGetWorkMat(elem_Sdot,&KlocSdot);CHKERRQ(ierr);
		ierr = IGAElementGetWorkVec(elem_Sdot,&FlocSdot);CHKERRQ(ierr);
		//ierr = IGAElementGetValues(elemVs,arrayVsSdot,&VsSdot);CHKERRQ(ierr);
		ierr = IGAElementGetValues(elem_S,arrayS0Sdot,&S0Sdot);CHKERRQ(ierr);

		// FormSystem loop
		while (IGAElementNextFormSystem(elem_Sdot,&wtfSdot,&wtf2Sdot)) 
		{
		// Quadrature loop
			ierr = IGAElementBeginPoint(elem_Sdot,&point_Sdot);CHKERRQ(ierr);
			//ierr = IGAElementBeginPoint(elemVs,&pointVs);CHKERRQ(ierr);
			ierr = IGAElementBeginPoint(elem_S,&point_S);CHKERRQ(ierr);

			while (IGAElementNextPoint(elem_Sdot,point_Sdot))
			{
				if(point_Sdot->atboundary==1)
				{
					//Code boundary condition, will be required for general case
				}
				if(point_Sdot->atboundary==0 && point_S->atboundary==0)
				{
					//IGAElementNextPoint(elemVs,pointVs);
					IGAElementNextPoint(elem_S,point_S);

					ierr = IGAPointGetWorkMat(point_Sdot,&Kpoint_Sdot);CHKERRQ(ierr);
					ierr = IGAPointGetWorkVec(point_Sdot,&Fpoint_Sdot);CHKERRQ(ierr);
					//ierr = SdotFunc(point_Sdot,pointVs,point_S,Kpoint_Sdot,Fpoint_Sdot,VsSdot,S0Sdot,&user);CHKERRQ(ierr);
					ierr = SdotFunc(point_Sdot,point_S,Kpoint_Sdot,Fpoint_Sdot,Int1a,Int2a,S_norm,S0Sdot,&user);CHKERRQ(ierr);
					ierr = IGAPointAddMat(point_Sdot,Kpoint_Sdot,KlocSdot);CHKERRQ(ierr);
					ierr = IGAPointAddVec(point_Sdot,Fpoint_Sdot,FlocSdot);CHKERRQ(ierr);
				}
			}
			//while (pointVs->index != -1)
			//{
			//	IGAElementNextPoint(elemVs,pointVs);
			//}
			
			while (point_S->index != -1)
			{
				IGAElementNextPoint(elem_S,point_S);
			}
			ierr = IGAElementEndPoint(elem_Sdot,&point_Sdot);CHKERRQ(ierr);
			//ierr = IGAElementEndPoint(elemVs,&pointVs);CHKERRQ(ierr);
			ierr = IGAElementEndPoint(elem_S,&point_S);CHKERRQ(ierr);
		}

		//ierr = IGAElementFixSystem(elem_Sdot,KlocSdot,FlocSdot);CHKERRQ(ierr);					//This sets Dirichlet condition ¿? (Yes, this applies the conditions from IGASetBoundaryValue)
		ierr = IGAElementAssembleMat(elem_Sdot,KlocSdot,KSdot);CHKERRQ(ierr);
		ierr = IGAElementAssembleVec(elem_Sdot,FlocSdot,FSdot);CHKERRQ(ierr);
	}
	//IGANextElement(igaVs,elemVs);
	IGANextElement(iga_S,elem_S);

	ierr = IGAEndElement(igaSdot,&elem_Sdot);CHKERRQ(ierr);
	//ierr = IGAEndElement(igaVs,&elemVs);CHKERRQ(ierr);
	ierr = IGAEndElement(iga_S,&elem_S);CHKERRQ(ierr);

	// Restore local vectors u, Z0, Chi0 and arrays
	//ierr = IGARestoreLocalVecArray(igaVs,VsSmooth,&localVsSdot,&arrayVsSdot);CHKERRQ(ierr);
	ierr = IGARestoreLocalVecArray(iga_S,S0,&localS0Sdot,&arrayS0Sdot);CHKERRQ(ierr);
	ierr = MatAssemblyBegin(KSdot,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyEnd  (KSdot,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

	ierr = VecAssemblyBegin(FSdot);CHKERRQ(ierr);
	ierr = VecAssemblyEnd  (FSdot);CHKERRQ(ierr);

	ierr = KSPSetOperators(kspSdot,KSdot,KSdot);CHKERRQ(ierr);
	PC pcSdot;
	ierr = KSPGetPC(kspSdot,&pcSdot); CHKERRQ(ierr);
	ierr = PCSetType(pcSdot,PCLU); CHKERRQ(ierr);
	ierr = PCFactorSetMatSolverType(pcSdot,MATSOLVERMUMPS); CHKERRQ(ierr);
	//ierr = KSPSetFromOptions(kspSdot);CHKERRQ(ierr);
	//ierr = KSPSetTolerances(kspSdot,1.0e-16,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);CHKERRQ(ierr);
	ierr = KSPSolve(kspSdot,FSdot,Sdot);CHKERRQ(ierr);

	//ierr = KSPDestroy(&kspSdot);CHKERRQ(ierr);
	//ierr = MatDestroy(&KSdot);CHKERRQ(ierr);
	//ierr = VecDestroy(&FSdot);CHKERRQ(ierr);

	//Kill everything in S[0],S[2],S[4].S[6], just because I haven't been able to stop it from coming out from the solution
		ierr = VecGetSize(Sdot,&n);CHKERRQ(ierr);

		PetscInt *colsS;
		PetscReal *valsS;

		ierr = PetscMalloc1(n,&colsS);CHKERRQ(ierr);
		ierr = PetscMalloc1(n,&valsS);CHKERRQ(ierr);

		for(int i=0;i<n;i=i+8)
		{
			//A negative index in colsS means that index is ignored, the value -1.0 is passed as a way to detect if the ignoring is not working
			colsS[i]=i;		colsS[i+1]=-1;		colsS[i+2]=i+2;		colsS[i+3]=-1;		colsS[i+4]=i+4;		colsS[i+5]=-1;		colsS[i+6]=i+6;		colsS[i+7]=-1;
			valsS[i]=0.0;	valsS[i+1]=-1.0;	valsS[i+2]=0.0;		valsS[i+3]=-1.0;	valsS[i+4]=0.0;		valsS[i+5]=-1.0;	valsS[i+6]=0.0;		valsS[i+7]=-1.0;
		}

		ierr = VecSetOption(Sdot, VEC_IGNORE_NEGATIVE_INDICES,PETSC_TRUE);CHKERRQ(ierr);
		ierr = VecAssemblyBegin(Sdot);CHKERRQ(ierr);
		ierr = VecAssemblyEnd  (Sdot);CHKERRQ(ierr);

		//VecSetValues(Vec x,PetscInt ni,const PetscInt ix[],const PetscScalar y[],InsertMode iora)
		ierr = VecSetValues(Sdot,n,colsS,valsS,INSERT_VALUES);CHKERRQ(ierr);

		ierr = VecAssemblyBegin(Sdot);CHKERRQ(ierr);
		ierr = VecAssemblyEnd  (Sdot);CHKERRQ(ierr);
	//

	char nameSdot[512];//="/Input-S-2d-1.dat";
	sprintf(nameSdot,"%s%d%s","/Input-S-2d-",i+1,".dat");
	char pathSdot[1024];
	sprintf(pathSdot,"%s%s",direct,nameSdot);
	ierr = IGAWriteVec(igaSdot,Sdot,pathSdot);CHKERRQ(ierr);	
//

//System for L2 projection of gradz (debug)
	PetscPrintf(PETSC_COMM_WORLD,"\nSystem for grad(z) starting \n\n");
	T=time(NULL);
	tm=*localtime(&T);
	PetscPrintf(PETSC_COMM_WORLD,"Current time is %02d:%02d:%02d \n",tm.tm_hour,tm.tm_min,tm.tm_sec);
	IGA igaGradz;
	ierr = IGACreate(PETSC_COMM_WORLD,&igaGradz);CHKERRQ(ierr);
	ierr = IGASetDim(igaGradz,2);CHKERRQ(ierr);														//Spatial dimension of the problem
	ierr = IGASetDof(igaGradz,4);CHKERRQ(ierr);														//Number of degrees of freedom, per node
	ierr = IGASetOrder(igaGradz,1);CHKERRQ(ierr);													//Number of spatial derivatives to calculate
	ierr = IGASetFromOptions(igaGradz);CHKERRQ(ierr);												//Note: The order (or degree) of the shape functions is given by the mesh!
	ierr = IGARead(igaGradz,"./geometry.dat");CHKERRQ(ierr);
	
	for (dir=0; dir<2; dir++)
	{
		ierr = IGASetRuleType(igaGradz,dir,IGA_RULE_LEGENDRE);CHKERRQ(ierr);
		ierr = IGASetRuleSize(igaGradz,dir,6);CHKERRQ(ierr);
	}
	ierr = IGASetUp(igaGradz);CHKERRQ(ierr);
	
	Mat KGradz0;
	Vec gradZ0,FGradz0;
	ierr = IGACreateMat(igaGradz,&KGradz0);CHKERRQ(ierr);
	ierr = IGACreateVec(igaGradz,&gradZ0);CHKERRQ(ierr);
	ierr = IGACreateVec(igaGradz,&FGradz0);CHKERRQ(ierr);

	IGAPoint		pointGradz;
	IGAElement		elemGradz;							//element
	PetscReal		*KlocGradz,*FlocGradz;				//AA y BB
	PetscReal		*KpointGradz,*FpointGradz;			//KKK y FFF
	const PetscReal *arrayZ0Gradz;						//arrayU
	Vec				localZ0Gradz;						//localU
	PetscReal		*Z0Gradz;								//U

  	IGAFormSystem  wtfGradz;
 	void           *wtf2Gradz;

 	KSP kspGradz;
	ierr = IGACreateKSP(igaGradz,&kspGradz);CHKERRQ(ierr);

	//Get local vectors z0 and arrays
	ierr = IGAGetLocalVecArray(iga_z,Z0,&localZ0Gradz,&arrayZ0Gradz);CHKERRQ(ierr);

	//Element loop
	ierr = IGABeginElement(igaGradz,&elemGradz);CHKERRQ(ierr);
	ierr = IGABeginElement(iga_z,&elem_z);CHKERRQ(ierr);

	while (IGANextElement(igaGradz,elemGradz))
	{
		IGANextElement(iga_z,elem_z);

		ierr = IGAElementGetWorkMat(elemGradz,&KlocGradz);CHKERRQ(ierr);
		ierr = IGAElementGetWorkVec(elemGradz,&FlocGradz);CHKERRQ(ierr);
		ierr = IGAElementGetValues(elem_z,arrayZ0Gradz,&Z0Gradz);CHKERRQ(ierr);

		//FormSystem loop
		while (IGAElementNextFormSystem(elemGradz,&wtfGradz,&wtf2Gradz)) 
		{
			//Quadrature loop
			ierr = IGAElementBeginPoint(elemGradz,&pointGradz);CHKERRQ(ierr);
			ierr = IGAElementBeginPoint(elem_z,&point_z);CHKERRQ(ierr);

			while (IGAElementNextPoint(elemGradz,pointGradz)) 
			{
				if(pointGradz->atboundary==0 && point_z->atboundary==0)
				{
					IGAElementNextPoint(elem_z,point_z);
					ierr = IGAPointGetWorkMat(pointGradz,&KpointGradz);CHKERRQ(ierr);
					ierr = IGAPointGetWorkVec(pointGradz,&FpointGradz);CHKERRQ(ierr);
					ierr = gradz(pointGradz,point_z,KpointGradz,FpointGradz,Z0Gradz,NULL);CHKERRQ(ierr);
					ierr = IGAPointAddMat(pointGradz,KpointGradz,KlocGradz);CHKERRQ(ierr);
					ierr = IGAPointAddVec(pointGradz,FpointGradz,FlocGradz);CHKERRQ(ierr);
				}
			}
			IGAElementNextPoint(elem_z,point_z);

			ierr = IGAElementEndPoint(elemGradz,&pointGradz);CHKERRQ(ierr);
			while (point_z->index != -1)
			{
				IGAElementNextPoint(elem_z,point_z);
			}
			ierr = IGAElementEndPoint(elem_z,&point_z);CHKERRQ(ierr);
		}
		ierr = IGAElementAssembleMat(elemGradz,KlocGradz,KGradz0);CHKERRQ(ierr);
		ierr = IGAElementAssembleVec(elemGradz,FlocGradz,FGradz0);CHKERRQ(ierr);

	}
	IGANextElement(iga_z,elem_z);

	ierr = IGAEndElement(igaGradz,&elemGradz);CHKERRQ(ierr);
	ierr = IGAEndElement(iga_z,&elem_z);CHKERRQ(ierr);

	// Restore local vectors S0 and arrays
	ierr = IGARestoreLocalVecArray(iga_z,Z0,&localZ0Gradz,&arrayZ0Gradz);CHKERRQ(ierr);

	//Form system matrix and vector
	ierr = MatAssemblyBegin(KGradz0,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyEnd  (KGradz0,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	
	ierr = VecAssemblyBegin(FGradz0);CHKERRQ(ierr);
	ierr = VecAssemblyEnd  (FGradz0);CHKERRQ(ierr);

	ierr = KSPSetOperators(kspGradz,KGradz0,KGradz0);CHKERRQ(ierr);
	ierr = KSPSetFromOptions(kspGradz);CHKERRQ(ierr);
	ierr = KSPSetTolerances(kspGradz,1.0e-16,1.0e-30,PETSC_DEFAULT,PETSC_DEFAULT);CHKERRQ(ierr);
	ierr = KSPSolve(kspGradz,FGradz0,gradZ0);CHKERRQ(ierr);

	char nameGradz0[512];
	sprintf(nameGradz0,"%s%d%s","/Gradz-2d-",i,".dat");
	char pathGradz0[1024];
	sprintf(pathGradz0,"%s%s",direct,nameGradz0);
	ierr = IGAWriteVec(igaGradz,gradZ0,pathGradz0);CHKERRQ(ierr);
//

//System for L2 projection of Up (Is this neccesary? I can use the general system)
	PetscPrintf(PETSC_COMM_WORLD,"\nSystem for L2 projection of U^p starting \n\n");
	T=time(NULL);
	tm=*localtime(&T);
	PetscPrintf(PETSC_COMM_WORLD,"Current time is %02d:%02d:%02d \n",tm.tm_hour,tm.tm_min,tm.tm_sec);
	IGA igaUp;
	ierr = IGACreate(PETSC_COMM_WORLD,&igaUp);CHKERRQ(ierr);
	ierr = IGASetDim(igaUp,2);CHKERRQ(ierr);														//Spatial dimension of the problem
	ierr = IGASetDof(igaUp,4);CHKERRQ(ierr);														//Number of degrees of freedom, per node
	ierr = IGASetOrder(igaUp,1);CHKERRQ(ierr);													//Number of spatial derivatives to calculate
	ierr = IGASetFromOptions(igaUp);CHKERRQ(ierr);												//Note: The order (or degree) of the shape functions is given by the mesh!
	ierr = IGARead(igaUp,"./geometry.dat");CHKERRQ(ierr);
	
	for (dir=0; dir<2; dir++)
	{
		ierr = IGASetRuleType(igaUp,dir,IGA_RULE_LEGENDRE);CHKERRQ(ierr);
		ierr = IGASetRuleSize(igaUp,dir,6);CHKERRQ(ierr);
	}
	ierr = IGASetUp(igaUp);CHKERRQ(ierr);
	
	Mat KUp0;
	Vec Up0,FUp0;
	ierr = IGACreateMat(igaUp,&KUp0);CHKERRQ(ierr);
	ierr = IGACreateVec(igaUp,&Up0);CHKERRQ(ierr);
	ierr = IGACreateVec(igaUp,&FUp0);CHKERRQ(ierr);

	IGAPoint		pointUp;
	IGAElement		elemUp;							//element
	PetscReal		*KlocUp,*FlocUp;				//AA y BB
	PetscReal		*KpointUp,*FpointUp;			//KKK y FFF
	const PetscReal *arrayZ0Up,*arrayChiUp;			//arrayU
	Vec				localZ0Up,localChiUp;		//localU
	PetscReal		*Z0Up,*ChiUp;				//U

  	IGAFormSystem  wtfUp;
 	void           *wtf2Up;

 	KSP kspUp;
	ierr = IGACreateKSP(igaUp,&kspUp);CHKERRQ(ierr);

	//Get local vectors z0 and arrays
	ierr = IGAGetLocalVecArray(iga_z,Z0,&localZ0Up,&arrayZ0Up);CHKERRQ(ierr);
	ierr = IGAGetLocalVecArray(iga_chi,chi,&localChiUp,&arrayChiUp);CHKERRQ(ierr);

	//Element loop
	ierr = IGABeginElement(igaUp,&elemUp);CHKERRQ(ierr);
	ierr = IGABeginElement(iga_z,&elem_z);CHKERRQ(ierr);
	ierr = IGABeginElement(iga_chi,&elem_chi);CHKERRQ(ierr);

	while (IGANextElement(igaUp,elemUp))
	{
		IGANextElement(iga_z,elem_z);
		IGANextElement(iga_chi,elem_chi);

		ierr = IGAElementGetWorkMat(elemUp,&KlocUp);CHKERRQ(ierr);
		ierr = IGAElementGetWorkVec(elemUp,&FlocUp);CHKERRQ(ierr);
		ierr = IGAElementGetValues(elem_z,arrayZ0Up,&Z0Up);CHKERRQ(ierr);
		ierr = IGAElementGetValues(elemUp,arrayChiUp,&ChiUp);CHKERRQ(ierr);

		//FormSystem loop
		while (IGAElementNextFormSystem(elemUp,&wtfUp,&wtf2Up)) 
		{
			//Quadrature loop
			ierr = IGAElementBeginPoint(elemUp,&pointUp);CHKERRQ(ierr);
			ierr = IGAElementBeginPoint(elem_z,&point_z);CHKERRQ(ierr);
			ierr = IGAElementBeginPoint(elem_chi,&point_chi);CHKERRQ(ierr);

			while (IGAElementNextPoint(elemUp,pointUp)) 
			{

				if(pointUp->atboundary==0 && point_z->atboundary==0 && point_chi->atboundary==0)
				{
					IGAElementNextPoint(elem_z,point_z);
					IGAElementNextPoint(elem_chi,point_chi);

					ierr = IGAPointGetWorkMat(pointUp,&KpointUp);CHKERRQ(ierr);
					ierr = IGAPointGetWorkVec(pointUp,&FpointUp);CHKERRQ(ierr);
						 //SysUp(IGAPoint p,IGAPoint pZu,IGAPoint pChi,PetscReal *K,PetscReal *F,PetscReal *UZ,PetscReal *UChi,void *ctx)
					ierr = SysUp(pointUp,point_z,point_chi,KpointUp,FpointUp,Z0Up,ChiUp,NULL);CHKERRQ(ierr);
					ierr = IGAPointAddMat(pointUp,KpointUp,KlocUp);CHKERRQ(ierr);
					ierr = IGAPointAddVec(pointUp,FpointUp,FlocUp);CHKERRQ(ierr);
				}
			}
			IGAElementNextPoint(elem_z,point_z);
			IGAElementNextPoint(elem_chi,point_chi);
		
			ierr = IGAElementEndPoint(elemUp,&pointUp);CHKERRQ(ierr);
			while (point_z->index != -1)
			{
				IGAElementNextPoint(elem_z,point_z);
			}
			ierr = IGAElementEndPoint(elem_z,&point_z);CHKERRQ(ierr);
			while (point_chi->index != -1)
			{
				IGAElementNextPoint(elem_chi,point_chi);
			}
			ierr = IGAElementEndPoint(elem_chi,&point_chi);CHKERRQ(ierr);
		}
		ierr = IGAElementAssembleMat(elemUp,KlocUp,KUp0);CHKERRQ(ierr);
		ierr = IGAElementAssembleVec(elemUp,FlocUp,FUp0);CHKERRQ(ierr);
	}
	IGANextElement(iga_z,elem_z);
	IGANextElement(iga_chi,elem_chi);

	ierr = IGAEndElement(igaUp,&elemUp);CHKERRQ(ierr);
	ierr = IGAEndElement(iga_z,&elem_z);CHKERRQ(ierr);
	ierr = IGAEndElement(iga_chi,&elem_chi);CHKERRQ(ierr);

	// Restore local vectors and arrays
	ierr = IGARestoreLocalVecArray(iga_z,Z0,&localZ0Up,&arrayZ0Up);CHKERRQ(ierr);
	ierr = IGARestoreLocalVecArray(iga_chi,chi,&localChiUp,&arrayChiUp);CHKERRQ(ierr);

	//Form system matrix and vector
	ierr = MatAssemblyBegin(KUp0,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyEnd  (KUp0,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	
	ierr = VecAssemblyBegin(FUp0);CHKERRQ(ierr);
	ierr = VecAssemblyEnd  (FUp0);CHKERRQ(ierr);

	ierr = KSPSetOperators(kspUp,KUp0,KUp0);CHKERRQ(ierr);
	ierr = KSPSetFromOptions(kspUp);CHKERRQ(ierr);
	ierr = KSPSetTolerances(kspUp,1.0e-16,1.0e-30,PETSC_DEFAULT,PETSC_DEFAULT);CHKERRQ(ierr);
	ierr = KSPSolve(kspUp,FUp0,Up0);CHKERRQ(ierr);

	char nameUp0[512];
	sprintf(nameUp0,"%s%d%s","/Up-2d-",i,".dat");
	char pathUp0[1024];
	sprintf(pathUp0,"%s%s",direct,nameUp0);
	ierr = IGAWriteVec(igaUp,Up0,pathUp0);CHKERRQ(ierr);
//

//System for general L2 projection (debug, but not really important, just for looking at fields)
	PetscPrintf(PETSC_COMM_WORLD,"\nSystem for general L2 projection starting \n\n");
	PetscPrintf(PETSC_COMM_WORLD,"\nCurrently projecting C*Ue \n\n");
	T=time(NULL);
	tm=*localtime(&T);
	PetscPrintf(PETSC_COMM_WORLD,"Current time is %02d:%02d:%02d \n",tm.tm_hour,tm.tm_min,tm.tm_sec);
	IGA igaProj;
	ierr = IGACreate(PETSC_COMM_WORLD,&igaProj);CHKERRQ(ierr);
	ierr = IGASetDim(igaProj,2);CHKERRQ(ierr);														//Spatial dimension of the problem
	ierr = IGASetDof(igaProj,4);CHKERRQ(ierr);														//Number of degrees of freedom, per node
	ierr = IGASetOrder(igaProj,1);CHKERRQ(ierr);													//Number of spatial derivatives to calculate
	ierr = IGASetFromOptions(igaProj);CHKERRQ(ierr);												//Note: The order (or degree) of the shape functions is given by the mesh!
	ierr = IGARead(igaProj,"./geometry.dat");CHKERRQ(ierr);
	
	for (dir=0; dir<2; dir++)
	{
		ierr = IGASetRuleType(igaProj,dir,IGA_RULE_LEGENDRE);CHKERRQ(ierr);
		ierr = IGASetRuleSize(igaProj,dir,6);CHKERRQ(ierr);
	}
	ierr = IGASetUp(igaProj);CHKERRQ(ierr);
	
	Mat KProj;
	Vec proj0,FProj;
	ierr = IGACreateMat(igaProj,&KProj);CHKERRQ(ierr);
	ierr = IGACreateVec(igaProj,&proj0);CHKERRQ(ierr);
	ierr = IGACreateVec(igaProj,&FProj);CHKERRQ(ierr);

	IGAPoint		pointProj;
	IGAElement		elemProj;																				//element
	PetscReal		*KlocProj,*FlocProj;																	//AA y BB
	PetscReal		*KpointProj,*FpointProj;																//KKK y FFF
	const PetscReal *arrayAlProj,*arraySProj,*arrayVaProj,*arrayZ0Proj,*arrayChiProj,*arrayUProj;			//arrayU
	Vec				localAlProj,localSProj,localVaProj,localZ0Proj,localChiProj,localUProj;					//localU
	PetscReal		*AlProj,*SProj,*VaProj,*Z0Proj,*ChiProj,*UProj;											//U

  	IGAFormSystem  wtfProj;
 	void           *wtf2Proj;

 	KSP kspProj;
	ierr = IGACreateKSP(igaProj,&kspProj);CHKERRQ(ierr);

	//Get local vectors and arrays
	ierr = IGAGetLocalVecArray(iga_Al_hat,alInput,&localAlProj,&arrayAlProj);CHKERRQ(ierr);
	ierr = IGAGetLocalVecArray(iga_S,S0,&localSProj,&arraySProj);CHKERRQ(ierr);
	ierr = IGAGetLocalVecArray(igaVa,Va0,&localVaProj,&arrayVaProj);CHKERRQ(ierr);
	ierr = IGAGetLocalVecArray(iga_z,Z0,&localZ0Proj,&arrayZ0Proj);CHKERRQ(ierr);
	ierr = IGAGetLocalVecArray(iga_chi,chi,&localChiProj,&arrayChiProj);CHKERRQ(ierr);
	ierr = IGAGetLocalVecArray(iga_u,U,&localUProj,&arrayUProj);CHKERRQ(ierr);

	//Element loop
	ierr = IGABeginElement(igaProj,&elemProj);CHKERRQ(ierr);
	ierr = IGABeginElement(iga_Al_hat,&elem_Al_hat);CHKERRQ(ierr);
	ierr = IGABeginElement(iga_S,&elem_S);CHKERRQ(ierr);
	ierr = IGABeginElement(igaVa,&elemVa);CHKERRQ(ierr);
	ierr = IGABeginElement(iga_z,&elem_z);CHKERRQ(ierr);
	ierr = IGABeginElement(iga_chi,&elem_chi);CHKERRQ(ierr);
	ierr = IGABeginElement(iga_u,&elem_u);CHKERRQ(ierr);

	while (IGANextElement(igaProj,elemProj))
	{
		IGANextElement(iga_Al_hat,elem_Al_hat);
		IGANextElement(iga_S,elem_S);
		IGANextElement(igaVa,elemVa);
		IGANextElement(iga_z,elem_z);
		IGANextElement(iga_chi,elem_chi);
		IGANextElement(iga_u,elem_u);

		ierr = IGAElementGetWorkMat(elemProj,&KlocProj);CHKERRQ(ierr);
		ierr = IGAElementGetWorkVec(elemProj,&FlocProj);CHKERRQ(ierr);
		ierr = IGAElementGetValues(elem_Al_hat,arrayAlProj,&AlProj);CHKERRQ(ierr);
		ierr = IGAElementGetValues(elem_S,arraySProj,&SProj);CHKERRQ(ierr);
		ierr = IGAElementGetValues(elemVa,arrayVaProj,&VaProj);CHKERRQ(ierr);
		ierr = IGAElementGetValues(elem_z,arrayZ0Proj,&Z0Proj);CHKERRQ(ierr);
		ierr = IGAElementGetValues(elem_chi,arrayChiProj,&ChiProj);CHKERRQ(ierr);
		ierr = IGAElementGetValues(elem_u,arrayUProj,&UProj);CHKERRQ(ierr);

		//FormSystem loop
		while (IGAElementNextFormSystem(elemProj,&wtfProj,&wtf2Proj)) 
		{
			//Quadrature loop
			ierr = IGAElementBeginPoint(elemProj,&pointProj);CHKERRQ(ierr);
			ierr = IGAElementBeginPoint(elem_Al_hat,&point_Al_hat);CHKERRQ(ierr);
			ierr = IGAElementBeginPoint(elem_S,&point_S);CHKERRQ(ierr);
			ierr = IGAElementBeginPoint(elemVa,&pointVa);CHKERRQ(ierr);
			ierr = IGAElementBeginPoint(elem_z,&point_z);CHKERRQ(ierr);
			ierr = IGAElementBeginPoint(elem_chi,&point_chi);CHKERRQ(ierr);
			ierr = IGAElementBeginPoint(elem_u,&point_u);CHKERRQ(ierr);

			while (IGAElementNextPoint(elemProj,pointProj)) 
			{
				if(pointProj->atboundary==0 && point_Al_hat->atboundary==0 && point_S->atboundary==0 && pointVa->atboundary==0)
				{
					IGAElementNextPoint(elem_Al_hat,point_Al_hat);
					IGAElementNextPoint(elem_S,point_S);
					IGAElementNextPoint(elemVa,pointVa);
					IGAElementNextPoint(elem_z,point_z);
					IGAElementNextPoint(elem_chi,point_chi);
					IGAElementNextPoint(elem_u,point_u);

					ierr = IGAPointGetWorkMat(pointProj,&KpointProj);CHKERRQ(ierr);
					ierr = IGAPointGetWorkVec(pointProj,&FpointProj);CHKERRQ(ierr);
					//proj(IGAPoint p,IGAPoint pAl,IGAPoint pS, IGAPoint pVa,IGAPoint pZu,IGAPoint pu,IGAPoint pChi,PetscReal *K,PetscReal *F,PetscReal *UAl,PetscReal *US,PetscReal *UVa,
					//	PetscReal *UZu,PetscReal *UChi,PetscReal *Uu,PetscReal Vs1,PetscReal Vs2,PetscReal normS,void *ctx)
					ierr = proj(pointProj,point_Al_hat,point_S,pointVa,point_z,point_u,point_chi,KpointProj,FpointProj,AlProj,SProj,VaProj,Z0Proj,ChiProj,UProj,Int1a,Int2a,S_norm,NULL);CHKERRQ(ierr);
					ierr = IGAPointAddMat(pointProj,KpointProj,KlocProj);CHKERRQ(ierr);
					ierr = IGAPointAddVec(pointProj,FpointProj,FlocProj);CHKERRQ(ierr);
				}
			}
			IGAElementNextPoint(elem_Al_hat,point_Al_hat);
			IGAElementNextPoint(elem_S,point_S);
			IGAElementNextPoint(elemVa,pointVa);
			IGAElementNextPoint(elem_z,point_z);
			IGAElementNextPoint(elem_chi,point_chi);
			IGAElementNextPoint(elem_u,point_u);

			ierr = IGAElementEndPoint(elemProj,&pointProj);CHKERRQ(ierr);
			ierr = IGAElementEndPoint(elem_Al_hat,&point_Al_hat);CHKERRQ(ierr);
			ierr = IGAElementEndPoint(elem_S,&point_S);CHKERRQ(ierr);
			ierr = IGAElementEndPoint(elemVa,&pointVa);CHKERRQ(ierr);
			ierr = IGAElementEndPoint(elem_z,&point_z);CHKERRQ(ierr);
			ierr = IGAElementEndPoint(elem_chi,&point_chi);CHKERRQ(ierr);
			ierr = IGAElementEndPoint(elem_u,&point_u);CHKERRQ(ierr);
		}
		ierr = IGAElementAssembleMat(elemProj,KlocProj,KProj);CHKERRQ(ierr);
		ierr = IGAElementAssembleVec(elemProj,FlocProj,FProj);CHKERRQ(ierr);
	}
	IGANextElement(iga_Al_hat,elem_Al_hat);
	IGANextElement(iga_S,elem_S);
	IGANextElement(igaVa,elemVa);
	IGANextElement(iga_z,elem_z);
	IGANextElement(iga_chi,elem_chi);
	IGANextElement(iga_u,elem_u);

	ierr = IGAEndElement(igaProj,&elemProj);CHKERRQ(ierr);
	ierr = IGAEndElement(iga_Al_hat,&elem_Al_hat);CHKERRQ(ierr);
	ierr = IGAEndElement(iga_S,&elem_S);CHKERRQ(ierr);
	ierr = IGAEndElement(igaVa,&elemVa);CHKERRQ(ierr);
	ierr = IGAEndElement(iga_z,&elem_z);CHKERRQ(ierr);
	ierr = IGAEndElement(iga_chi,&elem_chi);CHKERRQ(ierr);
	ierr = IGAEndElement(iga_u,&elem_u);CHKERRQ(ierr);

	// Restore local vectors S0 and arrays
	ierr = IGARestoreLocalVecArray(iga_Al_hat,alInput,&localAlProj,&arrayAlProj);CHKERRQ(ierr);
	ierr = IGARestoreLocalVecArray(iga_S,S0,&localSProj,&arraySProj);CHKERRQ(ierr);
	ierr = IGARestoreLocalVecArray(igaVa,Va0,&localVaProj,&arrayVaProj);CHKERRQ(ierr);
	ierr = IGARestoreLocalVecArray(iga_z,Z0,&localZ0Proj,&arrayZ0Proj);CHKERRQ(ierr);
	ierr = IGARestoreLocalVecArray(iga_chi,chi,&localChiProj,&arrayChiProj);CHKERRQ(ierr);
	ierr = IGARestoreLocalVecArray(iga_u,U,&localUProj,&arrayUProj);CHKERRQ(ierr);

	//Form system matrix and vector
	ierr = MatAssemblyBegin(KProj,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyEnd  (KProj,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	
	ierr = VecAssemblyBegin(FProj);CHKERRQ(ierr);
	ierr = VecAssemblyEnd  (FProj);CHKERRQ(ierr);

	ierr = KSPSetOperators(kspProj,KProj,KProj);CHKERRQ(ierr);
	ierr = KSPSetFromOptions(kspProj);CHKERRQ(ierr);
	ierr = KSPSetTolerances(kspProj,1.0e-16,1.0e-30,PETSC_DEFAULT,PETSC_DEFAULT);CHKERRQ(ierr);
	ierr = KSPSolve(kspProj,FProj,proj0);CHKERRQ(ierr);

	char nameProj[512];
	sprintf(nameProj,"%s%d%s","/Proj-2d-",i,".dat");
	char pathProj[1024];
	sprintf(pathProj,"%s%s",direct,nameProj);
	ierr = IGAWriteVec(igaProj,proj0,pathProj);CHKERRQ(ierr);
//

//All the necessary fields for time t=0 are calculated. Now loop for time update

//Create things for the solution loop
	//Things for \dot(z)
		IGA igaZdot;
		ierr = IGACreate(PETSC_COMM_WORLD,&igaZdot);CHKERRQ(ierr);
		ierr = IGASetDim(igaZdot,2);CHKERRQ(ierr);													//Spatial dimension of the problem
		ierr = IGASetDof(igaZdot,2);CHKERRQ(ierr);													//Number of degrees of freedom, per node
		ierr = IGASetOrder(igaZdot,1);CHKERRQ(ierr);													//Number of spatial derivatives to calculate
		ierr = IGASetFromOptions(igaZdot);CHKERRQ(ierr);												//Note: The order (or degree) of the shape functions is given by the mesh!
		ierr = IGARead(igaZdot,"./geometry.dat");CHKERRQ(ierr);
		
		for (dir=0; dir<2; dir++)
		{
			ierr = IGASetRuleType(igaZdot,dir,IGA_RULE_LEGENDRE);CHKERRQ(ierr);
			ierr = IGASetRuleSize(igaZdot,dir,6);CHKERRQ(ierr);
		}
		ierr = IGASetUp(igaZdot);CHKERRQ(ierr);

		for (dir=0; dir<2; dir++) 
		{
			for (side=0; side<2; side++) 
			{
				//ierr = IGASetBoundaryValue(iga,dir,side,dof,0.0);CHKERRQ(ierr);					// Dirichlet boundary conditions
				//ierr = IGASetBoundaryForm(igaZdot,dir,side,PETSC_TRUE);CHKERRQ(ierr);  				// Neumann boundary conditions
			}
		}

		Mat KZdot;
		Vec Zdot,FZdot;
		ierr = IGACreateMat(igaZdot,&KZdot);CHKERRQ(ierr);
		ierr = IGACreateVec(igaZdot,&Zdot);CHKERRQ(ierr);
		ierr = IGACreateVec(igaZdot,&FZdot);CHKERRQ(ierr);

		IGAPoint		pointZdot;															//point
		IGAElement		elemZdot;															//element
		PetscReal		*KlocZdot,*FlocZdot;												//AA y BB
		PetscReal		*KpointZdot,*FpointZdot;											//KKK y FFF
		const PetscReal	*arrayAlZdot,*arrayVaZdot,*arraySZdot;								//arrayU
		Vec				localAlZdot,localVaZdot,localSZdot;									//localU
		PetscReal		*AlZdot,*VaZdot,*SZdot;												//U

		IGAFormSystem	wtfZdot;
		void			*wtf2Zdot;

		//PetscInt n,m;
		//PetscInt rows, *cols;
		//PetscReal val, *vals;

		KSP kspZdot;
		PC pcZdot;
		ierr = IGACreateKSP(igaZdot,&kspZdot);CHKERRQ(ierr);
		char nameZdot[512];
		char pathZdot[1024];
		char nameZ[512];
		char pathZ[1024];
	//
//

//From here on all matrices will be reused, as in small deformation they don't change.
//Integrating only the right hand side vector speeds things up by about 2 orders of magnitude.
for (i=1;i<=2;i++)
{
	ierr = PetscPrintf(PETSC_COMM_WORLD,"\n\n Start of iteration %d \n\n",i);CHKERRQ(ierr);
//Creation of types and systems for the Initialization of S0
	PetscPrintf(PETSC_COMM_WORLD,"\nSystem for Initialization for S starting \n\n");
	T=time(NULL);
	tm=*localtime(&T);
	PetscPrintf(PETSC_COMM_WORLD,"Current time is %02d:%02d:%02d \n",tm.tm_hour,tm.tm_min,tm.tm_sec);
	
	sprintf(nameS,"%s%d%s","/Input-S-2d-",i,".dat");
	PetscPrintf(PETSC_COMM_WORLD,"Input is %s\n",nameS);

	sprintf(pathS,"%s%s",direct,nameS);
	ierr = IGAReadVec(iga_S,S0,pathS); CHKERRQ(ierr);
//

//Creation of types and systems for the Helmholtz decomposition of S, curl part for S based input
	PetscPrintf(PETSC_COMM_WORLD,"\nSystem for curl part of Helmholtz of S starting for iteration %d \n\n",i);
	T=time(NULL);
	tm=*localtime(&T);
	PetscPrintf(PETSC_COMM_WORLD,"Current time is %02d:%02d:%02d \n",tm.tm_hour,tm.tm_min,tm.tm_sec);
	
	ierr = VecZeroEntries(FSp);CHKERRQ(ierr);							//This sets all elements of the vector to 0
	//Not necessary to zero out S_perp, it's fully overwritten every time. We also don't zero the matrix, it's reused in every iteration, as it doesn't change.

	//Get local vectors S0 and arrays
	ierr = IGAGetLocalVecArray(iga_S,S0,&local_S0_Sp,&array_S0_Sp);CHKERRQ(ierr);

	//Element loop
	ierr = IGABeginElement(iga_Sp,&elem_Sp);CHKERRQ(ierr);
	ierr = IGABeginElement(iga_S,&elem_S);CHKERRQ(ierr);

	while (IGANextElement(iga_Sp,elem_Sp))
	{
		IGANextElement(iga_S,elem_S);

		ierr = IGAElementGetWorkVec(elem_Sp,&Floc_Sp);CHKERRQ(ierr);
		ierr = IGAElementGetValues(elem_S,array_S0_Sp,&S0_Sp);CHKERRQ(ierr);

		//FormSystem loop
		while (IGAElementNextFormSystem(elem_Sp,&wtfSp,&wtf2Sp)) 
		{
			//Quadrature loop
			ierr = IGAElementBeginPoint(elem_Sp,&point_Sp);CHKERRQ(ierr);
			ierr = IGAElementBeginPoint(elem_S,&point_S);CHKERRQ(ierr);
			
			while (IGAElementNextPoint(elem_Sp,point_Sp)) 
			{
				if(point_Sp->atboundary==0 && point_S->atboundary==0)
				{
					IGAElementNextPoint(elem_S,point_S);

					ierr = IGAPointGetWorkVec(point_Sp,&Fpoint_Sp);CHKERRQ(ierr);
					ierr = curlChiSF(point_Sp,point_S,Fpoint_Sp,S0_Sp,NULL);CHKERRQ(ierr);
					ierr = IGAPointAddVec(point_Sp,Fpoint_Sp,Floc_Sp);CHKERRQ(ierr);
				}
			}
			while (point_S->index != -1)
			{
				IGAElementNextPoint(elem_S,point_S);
			}
			ierr = IGAElementEndPoint(elem_Sp,&point_Sp);CHKERRQ(ierr);
			ierr = IGAElementEndPoint(elem_S,&point_S);CHKERRQ(ierr);
		}
		ierr = IGAElementFixSystem(elem_Sp,Kpoint_Sp,Floc_Sp);CHKERRQ(ierr);
		ierr = IGAElementAssembleVec(elem_Sp,Floc_Sp,FSp);CHKERRQ(ierr);

	}
	IGANextElement(iga_S,elem_S);
	ierr = IGAEndElement(iga_Sp,&elem_Sp);CHKERRQ(ierr);
	ierr = IGAEndElement(iga_S,&elem_S);CHKERRQ(ierr);

	// Restore local vectors S0 and arrays
	ierr = IGARestoreLocalVecArray(iga_S,S0,&local_S0_Sp,&array_S0_Sp);CHKERRQ(ierr);

	ierr = VecAssemblyBegin(FSp);CHKERRQ(ierr);
	ierr = VecAssemblyEnd  (FSp);CHKERRQ(ierr);

	ierr = KSPSolve(ksp_Sp,FSp,S_perp);CHKERRQ(ierr);
	
	sprintf(nameSp,"%s%d%s","/Sp-2d-",i,".dat");
	sprintf(pathSp,"%s%s",direct,nameSp);
	ierr = IGAWriteVec(iga_Sp,S_perp,pathSp);CHKERRQ(ierr);
//

//Creation of types and systems for the L2 projection of Alfa0-Sp:X
	PetscPrintf(PETSC_COMM_WORLD,"\nSystem for L2 projection for Alfa-Sp:X starting for iteration %d \n\n",i);
	T=time(NULL);
	tm=*localtime(&T);
	PetscPrintf(PETSC_COMM_WORLD,"Current time is %02d:%02d:%02d \n",tm.tm_hour,tm.tm_min,tm.tm_sec);

	ierr = VecZeroEntries(F_alpha_hat);CHKERRQ(ierr);						//This sets all elements of the vector to 0
	//Not necessary to zero out alInput and al_hat, they are fully overwritten every time.

	//Get local vectors S0 and arrays
	ierr = IGAGetLocalVecArray(iga_Sp,S_perp,&local_Sperp_Al_hat,&array_Sperp_Al_hat);CHKERRQ(ierr);

	//Element loop
	ierr = IGABeginElement(iga_Al_hat,&elem_Al_hat);CHKERRQ(ierr);
	ierr = IGABeginElement(iga_Sp,&elem_Sp);CHKERRQ(ierr);

	while (IGANextElement(iga_Al_hat,elem_Al_hat))
	{
		IGANextElement(iga_Sp,elem_Sp);
		ierr = IGAElementGetWorkVec(elem_Al_hat,&Floc_Al_hat);CHKERRQ(ierr);
		ierr = IGAElementGetValues(elem_Sp,array_Sperp_Al_hat,&Sperp_Al_hat);CHKERRQ(ierr);

		//FormSystem loop
		while (IGAElementNextFormSystem(elem_Al_hat,&wtfAlhat,&wtf2Alhat)) 
		{
			//Quadrature loop
			ierr = IGAElementBeginPoint(elem_Al_hat,&point_Al_hat);CHKERRQ(ierr);
			ierr = IGAElementBeginPoint(elem_Sp,&point_Sp);CHKERRQ(ierr);
			while (IGAElementNextPoint(elem_Al_hat,point_Al_hat)) 
			{
				if(point_Al_hat->atboundary==0 && point_Sp->atboundary==0)
				{
					IGAElementNextPoint(elem_Sp,point_Sp);
					ierr = IGAPointGetWorkVec(point_Al_hat,&Fpoint_Al_hat);CHKERRQ(ierr);
					ierr = Sp_X(point_Al_hat,point_Sp,Fpoint_Al_hat,Sperp_Al_hat,NULL);CHKERRQ(ierr);
					ierr = IGAPointAddVec(point_Al_hat,Fpoint_Al_hat,Floc_Al_hat);CHKERRQ(ierr);
				}
			}
			IGAElementNextPoint(elem_Sp,point_Sp);
			
			ierr = IGAElementEndPoint(elem_Al_hat,&point_Al_hat);CHKERRQ(ierr);
			ierr = IGAElementEndPoint(elem_Sp,&point_Sp);CHKERRQ(ierr);
		}
		//ierr = IGAElementFixSystem(elem_Al_hat,Kloc_Al_hat,Floc_Al_hat);CHKERRQ(ierr);					//This sets Dirichlet condition ¿?
		ierr = IGAElementAssembleVec(elem_Al_hat,Floc_Al_hat,F_alpha_hat);CHKERRQ(ierr);

	}
	IGANextElement(iga_Sp,elem_Sp);
	ierr = IGAEndElement(iga_Al_hat,&elem_Al_hat);CHKERRQ(ierr);
	ierr = IGAEndElement(iga_Sp,&elem_Sp);CHKERRQ(ierr);

	// Restore local vectors S0 and arrays
	ierr = IGARestoreLocalVecArray(iga_Sp,S_perp,&local_Sperp_Al_hat,&array_Sperp_Al_hat);CHKERRQ(ierr);

	//Form right hand side vector
	ierr = VecAssemblyBegin(F_alpha_hat);CHKERRQ(ierr);
	ierr = VecAssemblyEnd  (F_alpha_hat);CHKERRQ(ierr);

	ierr = KSPSetOperators(ksp_Al_hat,KL2_2GDL,KL2_2GDL);CHKERRQ(ierr);
	ierr = KSPSetFromOptions(ksp_Al_hat);CHKERRQ(ierr);
	ierr = KSPSetTolerances(ksp_Al_hat,1.0e-25,1.0e-40,PETSC_DEFAULT,PETSC_DEFAULT);CHKERRQ(ierr);
	ierr = KSPSolve(ksp_Al_hat,F_alpha_hat,al_hat);CHKERRQ(ierr);

	ierr = IGAReadVec(iga_Al_hat,alInput,pathAlInput); CHKERRQ(ierr);

	//Computes y = alpha*x + y. VecAXPY(Vec y,PetscScalar alpha,Vec x)
	ierr = VecAXPY(al_hat,1.0,alInput); CHKERRQ(ierr);

	sprintf(nameAlp,"%s%d%s","/Alp-2d-",i,".dat");
	sprintf(pathAlp,"%s%s",direct,nameAlp);
	ierr = IGAWriteVec(iga_Al_hat,al_hat,pathAlp);CHKERRQ(ierr);
//

//Creation of types and systems for the Helmholtz decomposition of Up (or Ue), curl part
	PetscPrintf(PETSC_COMM_WORLD,"\nSystem for curl part of Helmholtz of Up starting for iteration %d \n\n",i);
	T=time(NULL);
	tm=*localtime(&T);
	PetscPrintf(PETSC_COMM_WORLD,"Current time is %02d:%02d:%02d \n",tm.tm_hour,tm.tm_min,tm.tm_sec);
	
	ierr = VecZeroEntries(F_chi);CHKERRQ(ierr);						//This sets all elements of the vector to 0
	//Not necessary to zero out chi, it's fully overwritten every time, we don't zero the matrix either, as we reuse it from the first part.
	//Reusing the matrix increases speed about tenfold

	//Get local vectors S0 and arrays
	ierr = IGAGetLocalVecArray(iga_Al_hat,al_hat,&local_Al_chi,&array_Al_chi);CHKERRQ(ierr);

	//Element loop
	ierr = IGABeginElement(iga_chi,&elem_chi);CHKERRQ(ierr);
	ierr = IGABeginElement(iga_Al_hat,&elem_Al_hat);CHKERRQ(ierr);

	while (IGANextElement(iga_chi,elem_chi))
	{
		IGANextElement(iga_Al_hat,elem_Al_hat);
		ierr = IGAElementGetWorkVec(elem_chi,&Floc_chi);CHKERRQ(ierr);
		ierr = IGAElementGetValues(elem_Al_hat,array_Al_chi,&Al_chi);CHKERRQ(ierr);

		//FormSystem loop
		while (IGAElementNextFormSystem(elem_chi,&wtfchi,&wtf2chi)) 
		{
			//Quadrature loop
			ierr = IGAElementBeginPoint(elem_chi,&point_chi);CHKERRQ(ierr);
			ierr = IGAElementBeginPoint(elem_Al_hat,&point_Al_hat);CHKERRQ(ierr);

			while (IGAElementNextPoint(elem_chi,point_chi)) 
			{
				if(point_chi->atboundary==0 && point_Al_hat->atboundary==0)
				{
					IGAElementNextPoint(elem_Al_hat,point_Al_hat);

					ierr = IGAPointGetWorkVec(point_chi,&Fpoint_chi);CHKERRQ(ierr);
					ierr = curlChiUF(point_chi,point_Al_hat,Fpoint_chi,Al_chi,NULL);CHKERRQ(ierr);
					ierr = IGAPointAddVec(point_chi,Fpoint_chi,Floc_chi);CHKERRQ(ierr);
				}
			}
			IGAElementNextPoint(elem_Al_hat,point_Al_hat);

			ierr = IGAElementEndPoint(elem_chi,&point_chi);CHKERRQ(ierr);
			ierr = IGAElementEndPoint(elem_Al_hat,&point_Al_hat);CHKERRQ(ierr);
		}
		ierr = IGAElementFixSystem(elem_chi,Kloc_chi,Floc_chi);CHKERRQ(ierr);
		ierr = IGAElementAssembleVec(elem_chi,Floc_chi,F_chi);CHKERRQ(ierr);

	}
	IGANextElement(iga_Al_hat,elem_Al_hat);
	ierr = IGAEndElement(iga_chi,&elem_chi);CHKERRQ(ierr);
	ierr = IGAEndElement(iga_Al_hat,&elem_Al_hat);CHKERRQ(ierr);

	// Restore local vectors S0 and arrays
	ierr = IGARestoreLocalVecArray(iga_Al_hat,al_hat,&local_Al_chi,&array_Al_chi);CHKERRQ(ierr);
	
	ierr = VecAssemblyBegin(F_chi);CHKERRQ(ierr);
	ierr = VecAssemblyEnd  (F_chi);CHKERRQ(ierr);

	ierr = KSPSolve(kspchiUp,F_chi,chi);CHKERRQ(ierr);

	sprintf(namechiUp,"%s%d%s","/ChiUp-2d-",i,".dat");
	sprintf(pathchiUp,"%s%s",direct,namechiUp);
	ierr = IGAWriteVec(iga_chi,chi,pathchiUp);CHKERRQ(ierr);
//

//System for evolution of z0 (debug)
	PetscPrintf(PETSC_COMM_WORLD,"\nSystem for dot{z} starting for iteration %d \n\n",i);
	T=time(NULL);
	tm=*localtime(&T);
	PetscPrintf(PETSC_COMM_WORLD,"Current time is %02d:%02d:%02d \n",tm.tm_hour,tm.tm_min,tm.tm_sec);

	ierr = MatZeroEntries(KZdot);CHKERRQ(ierr);						//This makes all non-zero elements of KchiS 0.0, while keeping the sparse structure of the matrix
	ierr = VecZeroEntries(FZdot);CHKERRQ(ierr);						//This sets all elements of the vector to 0

	// Get local vectors Z0 and Chi0 and arrays
	ierr = IGAGetLocalVecArray(iga_Al_hat,alInput,&localAlZdot,&arrayAlZdot);CHKERRQ(ierr);
	ierr = IGAGetLocalVecArray(igaVa,Va0,&localVaZdot,&arrayVaZdot);CHKERRQ(ierr);
	ierr = IGAGetLocalVecArray(iga_S,S0,&localSZdot,&arraySZdot);CHKERRQ(ierr);

	// Element loop
	ierr = IGABeginElement(igaZdot,&elemZdot);CHKERRQ(ierr);
	ierr = IGABeginElement(igaVa,&elemVa);CHKERRQ(ierr);
	ierr = IGABeginElement(iga_Al_hat,&elem_Al_hat);CHKERRQ(ierr);
	ierr = IGABeginElement(iga_S,&elem_S);CHKERRQ(ierr);

	while (IGANextElement(igaZdot,elemZdot)) 
	{
		IGANextElement(igaVa,elemVa);
		IGANextElement(iga_Al_hat,elem_Al_hat);
		IGANextElement(iga_S,elem_S);
		
		ierr = IGAElementGetWorkMat(elemZdot,&KlocZdot);CHKERRQ(ierr);
		ierr = IGAElementGetWorkVec(elemZdot,&FlocZdot);CHKERRQ(ierr);
		ierr = IGAElementGetValues(elemVa,arrayVaZdot,&VaZdot);CHKERRQ(ierr);
		ierr = IGAElementGetValues(elem_Al_hat,arrayAlZdot,&AlZdot);CHKERRQ(ierr);
		ierr = IGAElementGetValues(elem_S,arraySZdot,&SZdot);CHKERRQ(ierr);

		// FormSystem loop
		while (IGAElementNextFormSystem(elemZdot,&wtfZdot,&wtf2Zdot)) 
		{
		// Quadrature loop
			ierr = IGAElementBeginPoint(elemZdot,&pointZdot);CHKERRQ(ierr);
			ierr = IGAElementBeginPoint(elemVa,&pointVa);CHKERRQ(ierr);
			ierr = IGAElementBeginPoint(elem_Al_hat,&point_Al_hat);CHKERRQ(ierr);
			ierr = IGAElementBeginPoint(elem_S,&point_S);CHKERRQ(ierr);

			while (IGAElementNextPoint(elemZdot,pointZdot))
			{
				if(pointZdot->atboundary==1)
				{

				}
				if(pointZdot->atboundary==0 && point_Al_hat->atboundary==0 && pointVa->atboundary==0 && point_S->atboundary==0)
				{
					IGAElementNextPoint(elem_Al_hat,point_Al_hat);
					IGAElementNextPoint(elemVa,pointVa);
					IGAElementNextPoint(elem_S,point_S);

					ierr = IGAPointGetWorkMat(pointZdot,&KpointZdot);CHKERRQ(ierr);
					ierr = IGAPointGetWorkVec(pointZdot,&FpointZdot);CHKERRQ(ierr);
						 //ZdotSystem(IGAPoint p,IGAPoint pAl,IGAPoint pVa,IGAPoint pS,PetscReal *K,PetscReal *F,PetscReal *UAl,PetscReal *UVa,PetscReal *US,PetscReal Vs1,PetscReal Vs2,PetscReal normS, void *ctx)
					ierr = ZdotSystem(pointZdot,point_Al_hat,pointVa,point_S,KpointZdot,FpointZdot,AlZdot,VaZdot,SZdot,Int1a,Int2a,S_norm,NULL);CHKERRQ(ierr);
					ierr = IGAPointAddMat(pointZdot,KpointZdot,KlocZdot);CHKERRQ(ierr);
					ierr = IGAPointAddVec(pointZdot,FpointZdot,FlocZdot);CHKERRQ(ierr);
				}
			}
			//while (point_Al_hat->index != -1)
			//{
				IGAElementNextPoint(elem_Al_hat,point_Al_hat);
			//}
			//while (pointVa->index != -1)
			//{
				IGAElementNextPoint(elemVa,pointVa);
			//}
			//while (point_S->index != -1)
			//{
				IGAElementNextPoint(elem_S,point_S);
			//}

			ierr = IGAElementEndPoint(elemVa,&pointVa);CHKERRQ(ierr);
			ierr = IGAElementEndPoint(elem_Al_hat,&point_Al_hat);CHKERRQ(ierr);
			ierr = IGAElementEndPoint(elem_S,&point_S);CHKERRQ(ierr);
			ierr = IGAElementEndPoint(elemZdot,&pointZdot);CHKERRQ(ierr);
		}
		ierr = IGAElementAssembleMat(elemZdot,KlocZdot,KZdot);CHKERRQ(ierr);
		ierr = IGAElementAssembleVec(elemZdot,FlocZdot,FZdot);CHKERRQ(ierr);
	}
	IGANextElement(igaVa,elemVa);
	IGANextElement(iga_Al_hat,elem_Al_hat);
	IGANextElement(iga_S,elem_S);

	ierr = IGAEndElement(igaZdot,&elemZdot);CHKERRQ(ierr);
	ierr = IGAEndElement(iga_Al_hat,&elem_Al_hat);CHKERRQ(ierr);
	ierr = IGAEndElement(igaVa,&elemVa);CHKERRQ(ierr);
	ierr = IGAEndElement(iga_S,&elem_S);CHKERRQ(ierr);

	// Restore local vectors u, Z0, Chi0 and arrays
	ierr = IGARestoreLocalVecArray(igaVa,Va0,&localVaZdot,&arrayVaZdot);CHKERRQ(ierr);
	ierr = IGARestoreLocalVecArray(iga_Al_hat,alInput,&localAlZdot,&arrayAlZdot);CHKERRQ(ierr);
	ierr = IGARestoreLocalVecArray(iga_S,S0,&localSZdot,&arraySZdot);CHKERRQ(ierr);
	
	ierr = MatAssemblyBegin(KZdot,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyEnd  (KZdot,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

	//Fix a single point to make solution unique, in this case, lower left corner.
		ierr = MatGetSize(KZdot,&n,&m);CHKERRQ(ierr);
		ierr = MatSetOption(KZdot, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);CHKERRQ(ierr);

		ierr = PetscMalloc1(m,&cols);CHKERRQ(ierr);
		ierr = PetscMalloc1(m,&vals);CHKERRQ(ierr);

		for(int i=0;i<m;i++)
		{
			cols[i]=i;
			vals[i]=0.0;
		}

		rows=0;
		vals[rows]=1.0e12;
		ierr = MatSetValues(KZdot,1,&rows,m,cols,vals,INSERT_VALUES);CHKERRQ(ierr);		//To keep matrix symmetric this zeros a row and sets diagonal value to 1e12
		ierr = MatSetValues(KZdot,n,cols,1,&rows,vals,INSERT_VALUES);CHKERRQ(ierr);		//then this zeros out the column and sets diagonal to 1e12
		vals[rows]=0.0;

		rows=1;
		vals[rows]=1.0e12;
		ierr = MatSetValues(KZdot,1,&rows,m,cols,vals,INSERT_VALUES);CHKERRQ(ierr);
		ierr = MatSetValues(KZdot,n,cols,1,&rows,vals,INSERT_VALUES);CHKERRQ(ierr);
		vals[rows]=0.0;

		ierr = MatAssemblyBegin(KZdot,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
		ierr = MatAssemblyEnd  (KZdot,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

		//Here we set values to the vector directly (to impose Dirichlet condition in a single point)
		ierr = VecAssemblyBegin(FZdot);CHKERRQ(ierr);
		ierr = VecAssemblyEnd  (FZdot);CHKERRQ(ierr);

		rows=0;
		val=0.0;
		ierr = VecSetValue(FZdot,rows,val,INSERT_VALUES);CHKERRQ(ierr);

		rows=1;
		val=0.0;
		ierr = VecSetValue(FZdot,rows,val,INSERT_VALUES);CHKERRQ(ierr);

		ierr = VecAssemblyBegin(FZdot);CHKERRQ(ierr);
		ierr = VecAssemblyEnd  (FZdot);CHKERRQ(ierr);
	//

	ierr = KSPSetOperators(kspZdot,KZdot,KZdot);CHKERRQ(ierr);
	ierr = KSPGetPC(kspZdot,&pcZdot); CHKERRQ(ierr);
	ierr = PCSetType(pcZdot,PCLU); CHKERRQ(ierr);
	ierr = PCFactorSetMatSolverType(pcZdot,MATSOLVERMUMPS); CHKERRQ(ierr);
	//ierr = KSPSetFromOptions(kspZdot);CHKERRQ(ierr);
	//ierr = KSPSetTolerances(kspZdot,1.0e-24,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);CHKERRQ(ierr);
	ierr = KSPSolve(kspZdot,FZdot,Zdot);CHKERRQ(ierr);

	sprintf(nameZdot,"%s%d%s","/Zdot-2d-",i,".dat");
	sprintf(pathZdot,"%s%s",direct,nameZdot);
	ierr = PetscPrintf(PETSC_COMM_WORLD,"%s \n",pathZdot);CHKERRQ(ierr);

	ierr = IGAWriteVec(igaZdot,Zdot,pathZdot);CHKERRQ(ierr);

	sprintf(nameZ,"%s%d%s","/Z0-2d-",i-1,".dat");
	sprintf(pathZ,"%s%s",direct,nameZ);
	ierr = IGAReadVec(iga_z,Z0,pathZ); CHKERRQ(ierr);

	//Now z_(t+1)=z_(t)+dt*\dot{z}
	//VecAXPY(Vec y,PetscScalar alpha,Vec x)
	ierr = VecAXPY(Z0,dt,Zdot);CHKERRQ(ierr);					//This does y= alpha*x + y, in this case z0=z0 + dt*\dot{z}

	sprintf(nameZ0,"%s%d%s","/Z0-2d-",i,".dat");
	sprintf(pathZ0,"%s%s",direct,nameZ0);
	ierr = IGAWriteVec(iga_z,Z0,pathZ0);CHKERRQ(ierr);			//Save updated Z to file
//

//System for u (debug)
	PetscPrintf(PETSC_COMM_WORLD,"\nSystem for u starting for iteration %d \n\n",i);
	T=time(NULL);
	tm=*localtime(&T);
	PetscPrintf(PETSC_COMM_WORLD,"Current time is %02d:%02d:%02d \n",tm.tm_hour,tm.tm_min,tm.tm_sec);
	
	ierr = MatZeroEntries(K_u);CHKERRQ(ierr);						//This makes all non-zero elements of K_u 0.0, while keeping the sparse structure of the matrix
	ierr = VecZeroEntries(F_u);CHKERRQ(ierr);						//This sets all elements of the vector to 0

	for (dir=0;dir<2;dir++) 
	{
		for (side=0;side<2;side++) 
		{
			ierr = IGASetBoundaryForm(iga_u,dir,side,PETSC_TRUE);CHKERRQ(ierr);  				// Neumann boundary conditions
		}
	}

	//If we are not fixing a single point, set Dirichlet conditions here
	//ierr = IGASetBoundaryValue(iga,dir,side,dof,value);CHKERRQ(ierr);					// Dirichlet boundary conditions
	//ierr = IGASetBoundaryValue(iga_u,0,0,0,0.0);CHKERRQ(ierr);	//Left side, 1st dof = 0
	//ierr = IGASetBoundaryValue(iga_u,0,0,1,0.0);CHKERRQ(ierr);	//Left side, 2nd dof = 0

	//ierr = IGASetBoundaryValue(iga_u,0,1,0,0.0);CHKERRQ(ierr);	//Right side, 1st dof=0
	//ierr = IGASetBoundaryValue(iga_u,0,1,1,0.0);CHKERRQ(ierr);	//Right side, 2nd dof=0

	ierr = IGASetBoundaryValue(iga_u,1,0,0,0.0);CHKERRQ(ierr);	//Bottom side, 1st dof=0
	ierr = IGASetBoundaryValue(iga_u,1,0,1,0.0);CHKERRQ(ierr);	//Bottom side, 2nd dof=0
	
	//ierr = IGASetBoundaryValue(iga_u,1,1,0,0.0);CHKERRQ(ierr);	//Top side, 1st dof=0
	//ierr = IGASetBoundaryValue(iga_u,1,1,1,0.0);CHKERRQ(ierr);	//Top side, 2nd dof=0

	// Get local vectors Chi0, Z0 and arrays
	ierr = IGAGetLocalVecArray(iga_chi,chi,&local_chi_u,&array_chi_u);CHKERRQ(ierr);
	ierr = IGAGetLocalVecArray(iga_z,Z0,&local_z_u,&array_z_u);CHKERRQ(ierr);

	// Element loop
	ierr = IGABeginElement(iga_u,&elem_u);CHKERRQ(ierr);
	ierr = IGABeginElement(iga_chi,&elem_chi);CHKERRQ(ierr);
	ierr = IGABeginElement(iga_z,&elem_z);CHKERRQ(ierr);

	while (IGANextElement(iga_u,elem_u))
	{
		IGANextElement(iga_chi,elem_chi);
		IGANextElement(iga_z,elem_z);

		ierr = IGAElementGetWorkMat(elem_u,&Kloc_u);CHKERRQ(ierr);
		ierr = IGAElementGetWorkVec(elem_u,&Floc_u);CHKERRQ(ierr);
		ierr = IGAElementGetValues(elem_chi,array_chi_u,&Chi_u);CHKERRQ(ierr);
		ierr = IGAElementGetValues(elem_z,array_z_u,&z_u);CHKERRQ(ierr);

		// FormSystem loop
		while (IGAElementNextFormSystem(elem_u,&wtfU0,&wtf2U0)) 
		{
		// Quadrature loop
			ierr = IGAElementBeginPoint(elem_u,&point_u);CHKERRQ(ierr);
			ierr = IGAElementBeginPoint(elem_chi,&point_chi);CHKERRQ(ierr);
			ierr = IGAElementBeginPoint(elem_z,&point_z);CHKERRQ(ierr);

			while (IGAElementNextPoint(elem_u,point_u))
			{
				if(point_u->atboundary==1)
				{
					ierr = IGAPointGetWorkVec(point_u,&Fpoint_u);CHKERRQ(ierr);
					//	   Usys(IGAPoint p, IGAPoint pChi, IGAPoint pZ, PetscReal *K, PetscReal *F, PetscReal *UChi, PetscReal *UZ, void *ctx)
					ierr = UsysF(point_u,point_chi,point_z,Fpoint_u,Chi_u,z_u,NULL);CHKERRQ(ierr);
					ierr = IGAPointAddVec(point_u,Fpoint_u,Floc_u);CHKERRQ(ierr);
				}
				if(point_u->atboundary==0 && point_chi->atboundary==0 && point_z->atboundary==0)
				{
					IGAElementNextPoint(elem_chi,point_chi);
					IGAElementNextPoint(elem_z,point_z);

					ierr = IGAPointGetWorkVec(point_u,&Fpoint_u);CHKERRQ(ierr);
					//	   Usys(IGAPoint p, IGAPoint pChi, IGAPoint pZ, PetscReal *K, PetscReal *F, PetscReal *UChi, PetscReal *UZ, void *ctx)
					ierr = UsysF(point_u,point_chi,point_z,Fpoint_u,Chi_u,z_u,NULL);CHKERRQ(ierr);
					ierr = IGAPointAddVec(point_u,Fpoint_u,Floc_u);CHKERRQ(ierr);
				}
			}
			while (point_chi->index != -1)
			{
				IGAElementNextPoint(elem_chi,point_chi);
			}
			while (point_z->index != -1)
			{
				IGAElementNextPoint(elem_z,point_z);
			}

			ierr = IGAElementEndPoint(elem_u,&point_u);CHKERRQ(ierr);
			ierr = IGAElementEndPoint(elem_chi,&point_chi);CHKERRQ(ierr);
			ierr = IGAElementEndPoint(elem_z,&point_z);CHKERRQ(ierr);
		}

		ierr = IGAElementFixSystem(elem_u,Kloc_u,Floc_u);CHKERRQ(ierr);					//This sets Dirichlet condition ¿? (Yes, this applies the conditions from IGASetBoundaryValue)

		ierr = IGAElementAssembleMat(elem_u,Kloc_u,K_u);CHKERRQ(ierr);
		ierr = IGAElementAssembleVec(elem_u,Floc_u,F_u);CHKERRQ(ierr);
	}
	IGANextElement(iga_chi,elem_chi);
	IGANextElement(iga_z,elem_z);

	ierr = IGAEndElement(iga_u,&elem_u);CHKERRQ(ierr);
	ierr = IGAEndElement(iga_chi,&elem_chi);CHKERRQ(ierr);
	ierr = IGAEndElement(iga_z,&elem_z);CHKERRQ(ierr);

	// Restore local vectors Chi0 and arrays
	ierr = IGARestoreLocalVecArray(iga_chi,chi,&local_chi_u,&array_chi_u);CHKERRQ(ierr);
	ierr = IGARestoreLocalVecArray(iga_z,Z0,&local_z_u,&array_z_u);CHKERRQ(ierr);

	ierr = MatAssemblyBegin(K_u,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyEnd  (K_u,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

	ierr = VecAssemblyBegin(F_u);CHKERRQ(ierr);
	ierr = VecAssemblyEnd  (F_u);CHKERRQ(ierr);

	ierr = KSPSetOperators(kspU,K_u,K_u);CHKERRQ(ierr);
	PC pcU;
	ierr = KSPGetPC(kspU,&pcU); CHKERRQ(ierr);
	ierr = PCSetType(pcU,PCLU); CHKERRQ(ierr);
	ierr = PCFactorSetMatSolverType(pcU,MATSOLVERMUMPS); CHKERRQ(ierr);
	//ierr = KSPSetFromOptions(kspU);CHKERRQ(ierr);
	//ierr = KSPSetTolerances(kspU,1e-24,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);CHKERRQ(ierr);
	ierr = KSPSolve(kspU,F_u,U);CHKERRQ(ierr);

	sprintf(nameU,"%s%d%s","/U-2d-",i,".dat");
	sprintf(pathU,"%s%s",direct,nameU);
	ierr = IGAWriteVec(iga_u,U,pathU);CHKERRQ(ierr);
//

//System for evolution of S (debug)
	//Now reuses K from the initial case, just recalculates RHS
	PetscPrintf(PETSC_COMM_WORLD,"\nSystem for S_dot starting for iteration %d \n\n",i);
	T=time(NULL);
	tm=*localtime(&T);
	PetscPrintf(PETSC_COMM_WORLD,"Current time is %02d:%02d:%02d \n",tm.tm_hour,tm.tm_min,tm.tm_sec);

	ierr = VecZeroEntries(FSdot);CHKERRQ(ierr);						//This sets all elements of the vector to 0
	//Not necessary to zero out Sdot, it's fully overwritten every time. We reuse the matrix KSdot from the initial step, it does not change.

	// Get local vectors Z0 and Chi0 and arrays
	//ierr = IGAGetLocalVecArray(igaVs,Vs0,&localVsSdot,&arrayVsSdot);CHKERRQ(ierr); Not anymore, use integral values
	ierr = IGAGetLocalVecArray(iga_S,S0,&localS0Sdot,&arrayS0Sdot);CHKERRQ(ierr);

	// Element loop
	//ierr = IGABeginElement(igaVs,&elemVs);CHKERRQ(ierr);
	ierr = IGABeginElement(iga_S,&elem_S);CHKERRQ(ierr);
	ierr = IGABeginElement(igaSdot,&elem_Sdot);CHKERRQ(ierr);

	while (IGANextElement(igaSdot,elem_Sdot)) 
	{
		//IGANextElement(igaVs,elemVs);
		IGANextElement(iga_S,elem_S);

		ierr = IGAElementGetWorkVec(elem_Sdot,&FlocSdot);CHKERRQ(ierr);
		//ierr = IGAElementGetValues(elemVs,arrayVsSdot,&VsSdot);CHKERRQ(ierr);
		ierr = IGAElementGetValues(elem_S,arrayS0Sdot,&S0Sdot);CHKERRQ(ierr);

		// FormSystem loop
		while (IGAElementNextFormSystem(elem_Sdot,&wtfSdot,&wtf2Sdot)) 
		{
		// Quadrature loop
			ierr = IGAElementBeginPoint(elem_Sdot,&point_Sdot);CHKERRQ(ierr);
			//ierr = IGAElementBeginPoint(elemVs,&pointVs);CHKERRQ(ierr);
			ierr = IGAElementBeginPoint(elem_S,&point_S);CHKERRQ(ierr);

			while (IGAElementNextPoint(elem_Sdot,point_Sdot))
			{
				if(point_Sdot->atboundary==1)
				{
					//Code boundary condition, will be required for general case
				}
				if(point_Sdot->atboundary==0 && point_S->atboundary==0)
				{
					//IGAElementNextPoint(elemVs,pointVs);
					IGAElementNextPoint(elem_S,point_S);

					ierr = IGAPointGetWorkVec(point_Sdot,&Fpoint_Sdot);CHKERRQ(ierr);
						 //SdotFuncF(IGAPoint p,IGAPoint pS,PetscReal *F,PetscReal Int1a,PetscReal Int2a,PetscReal normS,PetscReal *US,void *ctx)
					ierr = SdotFuncF(point_Sdot,point_S,Fpoint_Sdot,Int1a,Int2a,S_norm,S0Sdot,&user);CHKERRQ(ierr);
					ierr = IGAPointAddVec(point_Sdot,Fpoint_Sdot,FlocSdot);CHKERRQ(ierr);
				}
			}
			//while (pointVs->index != -1)
			//{
			//	IGAElementNextPoint(elemVs,pointVs);
			//}
			
			while (point_S->index != -1)
			{
				IGAElementNextPoint(elem_S,point_S);
			}
			ierr = IGAElementEndPoint(elem_Sdot,&point_Sdot);CHKERRQ(ierr);
			//ierr = IGAElementEndPoint(elemVs,&pointVs);CHKERRQ(ierr);
			ierr = IGAElementEndPoint(elem_S,&point_S);CHKERRQ(ierr);
		}
		//ierr = IGAElementFixSystemF(elem_Sdot,FlocSdot);CHKERRQ(ierr);					//This sets Dirichlet condition ¿? (Yes, this applies the conditions from IGASetBoundaryValue)
		ierr = IGAElementAssembleVec(elem_Sdot,FlocSdot,FSdot);CHKERRQ(ierr);
	}
	//IGANextElement(igaVs,elemVs);
	IGANextElement(iga_S,elem_S);

	ierr = IGAEndElement(igaSdot,&elem_Sdot);CHKERRQ(ierr);
	//ierr = IGAEndElement(igaVs,&elemVs);CHKERRQ(ierr);
	ierr = IGAEndElement(iga_S,&elem_S);CHKERRQ(ierr);

	// Restore local vectors u, Z0, Chi0 and arrays
	//ierr = IGARestoreLocalVecArray(igaVs,VsSmooth,&localVsSdot,&arrayVsSdot);CHKERRQ(ierr);
	ierr = IGARestoreLocalVecArray(iga_S,S0,&localS0Sdot,&arrayS0Sdot);CHKERRQ(ierr);

	ierr = VecAssemblyBegin(FSdot);CHKERRQ(ierr);
	ierr = VecAssemblyEnd  (FSdot);CHKERRQ(ierr);

	ierr = KSPSolve(kspSdot,FSdot,Sdot);CHKERRQ(ierr);

	//VecChop(Vec v, PetscReal tol) Sets anything with an absolute value less than the tolerance to 0
	//ierr = VecChop(Sdot,1e-8);CHKERRQ(ierr);

	//Kill everything in S[0],S[2],S[4].S[6], just because I haven't been able to stop it from coming out from the solution
		ierr = VecGetSize(Sdot,&n);CHKERRQ(ierr);

		for(int i=0;i<n;i=i+8)
		{
			//A negative index in colsS means that index is ignored, the value -1.0 is passed as a way to detect if the ignoring is not working
			colsS[i]=i;		colsS[i+1]=-1;		colsS[i+2]=i+2;		colsS[i+3]=-1;		colsS[i+4]=i+4;		colsS[i+5]=-1;		colsS[i+6]=i+6;		colsS[i+7]=-1;
			valsS[i]=0.0;	valsS[i+1]=-1.0;	valsS[i+2]=0.0;		valsS[i+3]=-1.0;	valsS[i+4]=0.0;		valsS[i+5]=-1.0;	valsS[i+6]=0.0;		valsS[i+7]=-1.0;
		}

		//Here we set values to the vector directly (to impose Dirichlet condition in a single point)
		ierr = VecAssemblyBegin(Sdot);CHKERRQ(ierr);
		ierr = VecAssemblyEnd  (Sdot);CHKERRQ(ierr);

		//VecSetValues(Vec x,PetscInt ni,const PetscInt ix[],const PetscScalar y[],InsertMode iora)
		ierr = VecSetValues(Sdot,n,colsS,valsS,INSERT_VALUES);CHKERRQ(ierr);

		ierr = VecAssemblyBegin(Sdot);CHKERRQ(ierr);
		ierr = VecAssemblyEnd  (Sdot);CHKERRQ(ierr);
	//

	sprintf(nameSdot,"%s%d%s","/Input-S-2d-",i+1,".dat");
	sprintf(pathSdot,"%s%s",direct,nameSdot);
	ierr = IGAWriteVec(igaSdot,Sdot,pathSdot);CHKERRQ(ierr);	
//

//System for L2 projection of V^{S} (debug)
	PetscPrintf(PETSC_COMM_WORLD,"\nSystem for V-S starting for iteration %d \n\n",i);
	T=time(NULL);
	tm=*localtime(&T);
	PetscPrintf(PETSC_COMM_WORLD,"Current time is %02d:%02d:%02d \n",tm.tm_hour,tm.tm_min,tm.tm_sec);

	ierr = VecZeroEntries(FVs);CHKERRQ(ierr);						//This sets all elements of the vector to 0

	// Get local vectors Z0 and Chi0 and arrays
	ierr = IGAGetLocalVecArray(iga_chi,chi,&localChi0Vs,&arrayChi0Vs);CHKERRQ(ierr);
	ierr = IGAGetLocalVecArray(iga_z,Z0,&localZ0Vs,&arrayZ0Vs);CHKERRQ(ierr);
	ierr = IGAGetLocalVecArray(iga_S,S0,&localS0Vs,&arrayS0Vs);CHKERRQ(ierr);
	ierr = IGAGetLocalVecArray(iga_u,U,&localUVs,&arrayUVs);CHKERRQ(ierr);

	// Element loop
	ierr = IGABeginElement(igaVs,&elemVs);CHKERRQ(ierr);
	ierr = IGABeginElement(iga_z,&elem_z);CHKERRQ(ierr);
	ierr = IGABeginElement(iga_chi,&elem_chi);CHKERRQ(ierr);
	ierr = IGABeginElement(iga_S,&elem_S);CHKERRQ(ierr);
	ierr = IGABeginElement(iga_u,&elem_u);CHKERRQ(ierr);

	while (IGANextElement(igaVs,elemVs)) 
	{
		IGANextElement(iga_z,elem_z);
		IGANextElement(iga_chi,elem_chi);
		IGANextElement(iga_S,elem_S);
		IGANextElement(iga_u,elem_u);

		//ierr = IGAElementGetWorkMat(elemVs,&KlocVs);CHKERRQ(ierr);
		ierr = IGAElementGetWorkVec(elemVs,&FlocVs);CHKERRQ(ierr);
		ierr = IGAElementGetValues(elem_z,arrayZ0Vs,&Z0Vs);CHKERRQ(ierr);
		ierr = IGAElementGetValues(elem_chi,arrayChi0Vs,&Chi0Vs);CHKERRQ(ierr);
		ierr = IGAElementGetValues(elem_S,arrayS0Vs,&S0Vs);CHKERRQ(ierr);
		ierr = IGAElementGetValues(elem_u,arrayUVs,&UVs);CHKERRQ(ierr);

		// FormSystem loop
		while (IGAElementNextFormSystem(elemVs,&wtfVs,&wtf2Vs)) 
		{
		// Quadrature loop
			ierr = IGAElementBeginPoint(elemVs,&pointVs);CHKERRQ(ierr);
			ierr = IGAElementBeginPoint(elem_z,&point_z);CHKERRQ(ierr);
			ierr = IGAElementBeginPoint(elem_chi,&point_chi);CHKERRQ(ierr);
			ierr = IGAElementBeginPoint(elem_S,&point_S);CHKERRQ(ierr);
			ierr = IGAElementBeginPoint(elem_u,&point_u);CHKERRQ(ierr);

			while (IGAElementNextPoint(elemVs,pointVs))
			{
				if(pointVs->atboundary==1)
				{

				}
				if(pointVs->atboundary==0 && point_z->atboundary==0 && point_chi->atboundary==0 && point_S->atboundary==0)
				{
					IGAElementNextPoint(elem_z,point_z);
					IGAElementNextPoint(elem_chi,point_chi);
					IGAElementNextPoint(elem_S,point_S);
					IGAElementNextPoint(elem_u,point_u);

					ierr = IGAPointGetWorkVec(pointVs,&FpointVs);CHKERRQ(ierr);
					//VSF(IGAPoint p,IGAPoint pChi,IGAPoint pZu,IGAPoint pu,IGAPoint pS,PetscReal *F, PetscReal *UChi,PetscReal *UZu,PetscReal *Uu,PetscReal *US,void *ctx)
					ierr = VSF(pointVs,point_chi,point_z,point_u,point_S,FpointVs,Chi0Vs,Z0Vs,UVs,S0Vs,NULL);CHKERRQ(ierr);
					ierr = IGAPointAddVec(pointVs,FpointVs,FlocVs);CHKERRQ(ierr);
				}
			}
			while (point_z->index != -1)
			{
				IGAElementNextPoint(elem_z,point_z);
			}
			while (point_chi->index != -1)
			{
				IGAElementNextPoint(elem_chi,point_chi);
			}
			while (point_S->index != -1)
			{
				IGAElementNextPoint(elem_S,point_S);
			}
			while (point_u->index != -1)
			{
				IGAElementNextPoint(elem_u,point_u);
			}
			ierr = IGAElementEndPoint(elemVs,&pointVs);CHKERRQ(ierr);
			ierr = IGAElementEndPoint(elem_z,&point_z);CHKERRQ(ierr);
			ierr = IGAElementEndPoint(elem_chi,&point_chi);CHKERRQ(ierr);
			ierr = IGAElementEndPoint(elem_S,&point_S);CHKERRQ(ierr);
			ierr = IGAElementEndPoint(elem_u,&point_u);CHKERRQ(ierr);
		}

		//ierr = IGAElementFixSystem(elem_Stress,KlocStress,FlocStress);CHKERRQ(ierr);					//This sets Dirichlet condition ¿? (Yes, this applies the conditions from IGASetBoundaryValue)
		ierr = IGAElementAssembleVec(elemVs,FlocVs,FVs);CHKERRQ(ierr);
	}
	IGANextElement(iga_z,elem_z);
	IGANextElement(iga_chi,elem_chi);
	IGANextElement(iga_S,elem_S);
	IGANextElement(iga_u,elem_u);

	ierr = IGAEndElement(igaVs,&elemVs);CHKERRQ(ierr);
	ierr = IGAEndElement(iga_z,&elem_z);CHKERRQ(ierr);
	ierr = IGAEndElement(iga_chi,&elem_chi);CHKERRQ(ierr);
	ierr = IGAEndElement(iga_S,&elem_S);CHKERRQ(ierr);
	ierr = IGAEndElement(iga_u,&elem_u);CHKERRQ(ierr);

	// Restore local vectors u, Z0, Chi0 and arrays
	ierr = IGARestoreLocalVecArray(iga_z,Z0,&localZ0Vs,&arrayZ0Vs);CHKERRQ(ierr);
	ierr = IGARestoreLocalVecArray(iga_chi,chi,&localChi0Vs,&arrayChi0Vs);CHKERRQ(ierr);
	ierr = IGARestoreLocalVecArray(iga_S,S0,&localS0Vs,&arrayS0Vs);CHKERRQ(ierr);
	ierr = IGARestoreLocalVecArray(iga_u,U,&localUVs,&arrayUVs);CHKERRQ(ierr);

	ierr = VecAssemblyBegin(FVs);CHKERRQ(ierr);
	ierr = VecAssemblyEnd  (FVs);CHKERRQ(ierr);

	ierr = KSPSolve(kspVs,FVs,Vs0);CHKERRQ(ierr);

	//VecChop(Vec v, PetscReal tol) Sets anything with an absolute value less than the tolerance to 0
	//ierr = VecChop(Vs0,1e-11);CHKERRQ(ierr);

	sprintf(nameVs,"%s%d%s","/Vs-2d-",i,".dat");
	sprintf(pathVs,"%s%s",direct,nameVs);
	ierr = IGAWriteVec(igaVs,Vs0,pathVs);CHKERRQ(ierr);	
//

//System for L2 projection of smoothed V^{S} (debug)
	PetscPrintf(PETSC_COMM_WORLD,"\nSystem for Smooth Vs starting for iteration %d \n\n",i);
	T=time(NULL);
	tm=*localtime(&T);
	PetscPrintf(PETSC_COMM_WORLD,"Current time is %02d:%02d:%02d \n",tm.tm_hour,tm.tm_min,tm.tm_sec);

	ierr = VecZeroEntries(FVsInt);CHKERRQ(ierr);						//This sets all elements of the vector to 0
	ierr = VecZeroEntries(FVsXi);CHKERRQ(ierr);							//This sets all elements of the vector to 0
	ierr = VecZeroEntries(Int1sVec);CHKERRQ(ierr);						//This sets all elements of the vector to 0
	ierr = VecZeroEntries(Int2sVec);CHKERRQ(ierr);						//This sets all elements of the vector to 0
	//First integrate velocity for both defects

	// Get local vectors V0 and InputAl and arrays
	ierr = IGAGetLocalVecArray(igaVs,Vs0,&localVs0,&arrayVs0);CHKERRQ(ierr);
	ierr = IGAGetLocalVecArray(iga_S,S0,&localSVs,&arraySVs);CHKERRQ(ierr);

	// Element loop
	ierr = IGABeginElement(igaVs,&elemVs);CHKERRQ(ierr);
	ierr = IGABeginElement(iga_S,&elem_S);CHKERRQ(ierr);

	while (IGANextElement(igaVs,elemVs))
	{
		IGANextElement(iga_S,elem_S);

		ierr = IGAElementGetWorkVec(elemVs,&locInt1a);CHKERRQ(ierr);
		ierr = IGAElementGetWorkVec(elemVs,&locInt2a);CHKERRQ(ierr);
		ierr = IGAElementGetValues(elemVs,arrayVs0,&Vs0Values);CHKERRQ(ierr);
		ierr = IGAElementGetValues(elem_S,arraySVs,&SVs);CHKERRQ(ierr);

		// FormSystem loop
		while (IGAElementNextFormSystem(elemVs,&wtfVsInt,&wtf2VsInt))
		{
		// Quadrature loop
			ierr = IGAElementBeginPoint(elemVs,&pointVs);CHKERRQ(ierr);
			ierr = IGAElementBeginPoint(elem_S,&point_S);CHKERRQ(ierr);
			
			while (IGAElementNextPoint(elemVs,pointVs))
			{
				if(pointVs->atboundary==1)
				{

				}
				if(pointVs->atboundary==0 && point_S->atboundary==0)
				{
					IGAElementNextPoint(elem_S,point_S);

					ierr = IGAPointGetWorkVec(pointVs,&PointInt1a);CHKERRQ(ierr);
					ierr = IGAPointGetWorkVec(pointVs,&PointInt2a);CHKERRQ(ierr);
					ierr = IGAPointGetWorkVec(pointVs,&PointInt);CHKERRQ(ierr);
					//Int_Xi_Vs(IGAPoint pV,IGAPoint pS,PetscReal *FInt1a,PetscReal *FInt2a,PetscReal *VS,PetscReal *US,void *ctx)
					ierr = Int_Xi_Vs(pointVs,point_S,PointInt1a,PointInt2a,Vs0Values,SVs,NULL);CHKERRQ(ierr);
					ierr = IGAPointAddVec(pointVs,PointInt1a,locInt1a);CHKERRQ(ierr);
					ierr = IGAPointAddVec(pointVs,PointInt2a,locInt2a);CHKERRQ(ierr);
				}
			}
			while (point_S->index != -1)
			{
				IGAElementNextPoint(elem_S,point_S);
			}
			ierr = IGAElementEndPoint(elemVs,&pointVs);CHKERRQ(ierr);
			ierr = IGAElementEndPoint(elem_S,&point_S);CHKERRQ(ierr);
		}
		ierr = IGAElementAssembleVec(elemVs,locInt1a,Int1sVec);CHKERRQ(ierr);
		ierr = IGAElementAssembleVec(elemVs,locInt2a,Int2sVec);CHKERRQ(ierr);
	}
	IGANextElement(iga_S,elem_S);

	ierr = IGAEndElement(igaVs,&elemVs);CHKERRQ(ierr);
	ierr = IGAEndElement(iga_S,&elem_S);CHKERRQ(ierr);

	ierr = IGARestoreLocalVecArray(iga_S,S0,&localSVs,&arraySVs);CHKERRQ(ierr);
	ierr = IGARestoreLocalVecArray(igaVs,Vs0,&localVs0,&arrayVs0);CHKERRQ(ierr);

	ierr = VecAssemblyBegin(Int1sVec);CHKERRQ(ierr);
	ierr = VecAssemblyEnd  (Int1sVec);CHKERRQ(ierr);
	ierr = VecAssemblyBegin(Int2sVec);CHKERRQ(ierr);
	ierr = VecAssemblyEnd  (Int2sVec);CHKERRQ(ierr);

	//Here add all values at Gauss point_S
	ierr = VecSum(Int1sVec,&Int1a);CHKERRQ(ierr);
	ierr = VecSum(Int2sVec,&Int2a);CHKERRQ(ierr);

	ierr = PetscPrintf(PETSC_COMM_WORLD,"Int1a=%f, Int2a=%f \n",Int1a,Int2a);CHKERRQ(ierr);
	//Use Int1a and Int2a as Vs(1) and Vs(2)
//

//System for L2 projection of V^{alpha} (debug)
	PetscPrintf(PETSC_COMM_WORLD,"\nSystem for V-alpha starting for iteration %d \n\n",i);
	T=time(NULL);
	tm=*localtime(&T);
	PetscPrintf(PETSC_COMM_WORLD,"Current time is %02d:%02d:%02d \n",tm.tm_hour,tm.tm_min,tm.tm_sec);
	
	ierr = VecZeroEntries(FVa);CHKERRQ(ierr);						//This sets all elements of the vector to 0

	// Get local vectors Z0 and Chi0 and arrays
	ierr = IGAGetLocalVecArray(iga_chi,chi,&localChi0Va,&arrayChi0Va);CHKERRQ(ierr);
	ierr = IGAGetLocalVecArray(iga_z,Z0,&localZ0Va,&arrayZ0Va);CHKERRQ(ierr);
	ierr = IGAGetLocalVecArray(iga_Al_hat,alInput,&localAl0Va,&arrayAl0Va);CHKERRQ(ierr);
	ierr = IGAGetLocalVecArray(iga_S,S0,&localSVa,&arraySVa);CHKERRQ(ierr);
	ierr = IGAGetLocalVecArray(iga_u,U,&localUVa,&arrayUVa);CHKERRQ(ierr);

	// Element loop
	ierr = IGABeginElement(igaVa,&elemVa);CHKERRQ(ierr);
	ierr = IGABeginElement(iga_z,&elem_z);CHKERRQ(ierr);
	ierr = IGABeginElement(iga_chi,&elem_chi);CHKERRQ(ierr);
	ierr = IGABeginElement(iga_Al_hat,&elem_Al_hat);CHKERRQ(ierr);
	ierr = IGABeginElement(iga_S,&elem_S);CHKERRQ(ierr);
	ierr = IGABeginElement(iga_u,&elem_u);CHKERRQ(ierr);

	while (IGANextElement(igaVa,elemVa)) 
	{
		IGANextElement(iga_z,elem_z);
		IGANextElement(iga_chi,elem_chi);
		IGANextElement(iga_Al_hat,elem_Al_hat);
		IGANextElement(iga_S,elem_S);
		IGANextElement(iga_u,elem_u);

		ierr = IGAElementGetWorkMat(elemVa,&KlocVa);CHKERRQ(ierr);
		ierr = IGAElementGetWorkVec(elemVa,&FlocVa);CHKERRQ(ierr);
		ierr = IGAElementGetValues(elem_z,arrayZ0Va,&Z0Va);CHKERRQ(ierr);
		ierr = IGAElementGetValues(elem_chi,arrayChi0Va,&Chi0Va);CHKERRQ(ierr);
		ierr = IGAElementGetValues(elem_Al_hat,arrayAl0Va,&Al0Va);CHKERRQ(ierr);
		ierr = IGAElementGetValues(elem_S,arraySVa,&SVa);CHKERRQ(ierr);
		ierr = IGAElementGetValues(elem_u,arrayUVa,&UVa);CHKERRQ(ierr);

		// FormSystem loop
		while (IGAElementNextFormSystem(elemVa,&wtfVa,&wtf2Va)) 
		{
		// Quadrature loop
			ierr = IGAElementBeginPoint(elemVa,&pointVa);CHKERRQ(ierr);
			ierr = IGAElementBeginPoint(elem_z,&point_z);CHKERRQ(ierr);
			ierr = IGAElementBeginPoint(elem_chi,&point_chi);CHKERRQ(ierr);
			ierr = IGAElementBeginPoint(elem_Al_hat,&point_Al_hat);CHKERRQ(ierr);
			ierr = IGAElementBeginPoint(elem_S,&point_S);CHKERRQ(ierr);
			ierr = IGAElementBeginPoint(elem_u,&point_u);CHKERRQ(ierr);

			while (IGAElementNextPoint(elemVa,pointVa))
			{
				if(pointVa->atboundary==1)
				{

				}
				if(pointVa->atboundary==0 && point_z->atboundary==0 && point_chi->atboundary==0 && point_Al_hat->atboundary==0 && point_S->atboundary==0 && point_u->atboundary==0)
				{
					IGAElementNextPoint(elem_z,point_z);
					IGAElementNextPoint(elem_chi,point_chi);
					IGAElementNextPoint(elem_Al_hat,point_Al_hat);
					IGAElementNextPoint(elem_S,point_S);
					IGAElementNextPoint(elem_u,point_u);

					ierr = IGAPointGetWorkVec(pointVa,&FpointVa);CHKERRQ(ierr);
					//ValphaF(IGAPoint p,IGAPoint pChi,IGAPoint pu,IGAPoint pZu,IGAPoint pAl,IGAPoint pS,PetscReal *F,PetscReal *Chi,PetscReal *Zu,PetscReal *US,PetscReal *Uu,PetscReal *UAl,void *ctx)
					ierr = ValphaF(pointVa,point_chi,point_u,point_z,point_Al_hat,point_S,FpointVa,Chi0Va,Z0Va,SVa,UVa,Al0Va,NULL);CHKERRQ(ierr);
					ierr = IGAPointAddVec(pointVa,FpointVa,FlocVa);CHKERRQ(ierr);
				}
			}
			while (point_z->index != -1)
			{
				IGAElementNextPoint(elem_z,point_z);
			}
			while (point_chi->index != -1)
			{
				IGAElementNextPoint(elem_chi,point_chi);
			}
			while (point_Al_hat->index != -1)
			{
				IGAElementNextPoint(elem_Al_hat,point_Al_hat);
			}
			while (point_S->index != -1)
			{
				IGAElementNextPoint(elem_S,point_S);
			}
			while (point_u->index != -1)
			{
				IGAElementNextPoint(elem_u,point_u);
			}
			ierr = IGAElementEndPoint(elemVa,&pointVa);CHKERRQ(ierr);
			ierr = IGAElementEndPoint(elem_z,&point_z);CHKERRQ(ierr);
			ierr = IGAElementEndPoint(elem_chi,&point_chi);CHKERRQ(ierr);
			ierr = IGAElementEndPoint(elem_Al_hat,&point_Al_hat);CHKERRQ(ierr);
			ierr = IGAElementEndPoint(elem_S,&point_S);CHKERRQ(ierr);
			ierr = IGAElementEndPoint(elem_u,&point_u);CHKERRQ(ierr);
		}

		//ierr = IGAElementFixSystem(elem_Stress,KlocStress,FlocStress);CHKERRQ(ierr);					//This sets Dirichlet condition ¿? (Yes, this applies the conditions from IGASetBoundaryValue)
		ierr = IGAElementAssembleVec(elemVa,FlocVa,FVa);CHKERRQ(ierr);
	}
	IGANextElement(iga_z,elem_z);
	IGANextElement(iga_chi,elem_chi);
	IGANextElement(iga_Al_hat,elem_Al_hat);
	IGANextElement(iga_S,elem_S);
	IGANextElement(iga_u,elem_u);

	ierr = IGAEndElement(igaVa,&elemVa);CHKERRQ(ierr);
	ierr = IGAEndElement(iga_z,&elem_z);CHKERRQ(ierr);
	ierr = IGAEndElement(iga_chi,&elem_chi);CHKERRQ(ierr);
	ierr = IGAEndElement(iga_Al_hat,&elem_Al_hat);CHKERRQ(ierr);
	ierr = IGAEndElement(iga_u,&elem_u);CHKERRQ(ierr);
	ierr = IGAEndElement(iga_S,&elem_S);CHKERRQ(ierr);

	// Restore local vectors u, Z0, Chi0 and arrays
	ierr = IGARestoreLocalVecArray(iga_z,Z0,&localZ0Va,&arrayZ0Va);CHKERRQ(ierr);
	ierr = IGARestoreLocalVecArray(iga_chi,chi,&localChi0Va,&arrayChi0Va);CHKERRQ(ierr);
	ierr = IGARestoreLocalVecArray(iga_Al_hat,alInput,&localAl0Va,&arrayAl0Va);CHKERRQ(ierr);
	ierr = IGARestoreLocalVecArray(iga_S,S0,&localSVa,&arraySVa);CHKERRQ(ierr);
	ierr = IGARestoreLocalVecArray(iga_u,U,&localUVa,&arrayUVa);CHKERRQ(ierr);

	ierr = VecAssemblyBegin(FVa);CHKERRQ(ierr);
	ierr = VecAssemblyEnd  (FVa);CHKERRQ(ierr);

	//ierr = KSPSetOperators(kspVa,KVa,KVa);CHKERRQ(ierr);
	//PC pcVa;
	//ierr = KSPGetPC(kspStress,&pcStress); CHKERRQ(ierr);
	//ierr = PCSetType(pcStress,PCLU); CHKERRQ(ierr);
	//ierr = PCFactorSetMatSolverType(pcStress,MATSOLVERMUMPS); CHKERRQ(ierr);
	//ierr = KSPSetFromOptions(kspVa);CHKERRQ(ierr);
	//ierr = KSPSetTolerances(kspVa,1.0e-24,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);CHKERRQ(ierr);
	ierr = KSPSolve(kspVa,FVa,Va0);CHKERRQ(ierr);

	//VecChop(Vec v, PetscReal tol) Sets anything with an absolute value less than the tolerance to 0
	ierr = VecChop(Va0,1e-11);CHKERRQ(ierr);

	sprintf(nameVa,"%s%d%s","/Va-2d-",i,".dat");
	sprintf(pathVa,"%s%s",direct,nameVa);
	ierr = IGAWriteVec(igaVa,Va0,pathVa);CHKERRQ(ierr);
//

//System for L2 projection of gradz (debug)
	PetscPrintf(PETSC_COMM_WORLD,"\nSystem for grad(z) starting for iteration %d \n\n",i);
	T=time(NULL);
	tm=*localtime(&T);
	PetscPrintf(PETSC_COMM_WORLD,"Current time is %02d:%02d:%02d \n",tm.tm_hour,tm.tm_min,tm.tm_sec);
	
	ierr = VecZeroEntries(FGradz0);CHKERRQ(ierr);						//This sets all elements of the vector to 0

	//Get local vectors z0 and arrays
	ierr = IGAGetLocalVecArray(iga_z,Z0,&localZ0Gradz,&arrayZ0Gradz);CHKERRQ(ierr);

	//Element loop
	ierr = IGABeginElement(igaGradz,&elemGradz);CHKERRQ(ierr);
	ierr = IGABeginElement(iga_z,&elem_z);CHKERRQ(ierr);

	while (IGANextElement(igaGradz,elemGradz))
	{
		IGANextElement(iga_z,elem_z);

		ierr = IGAElementGetWorkVec(elemGradz,&FlocGradz);CHKERRQ(ierr);
		ierr = IGAElementGetValues(elem_z,arrayZ0Gradz,&Z0Gradz);CHKERRQ(ierr);

		//FormSystem loop
		while (IGAElementNextFormSystem(elemGradz,&wtfGradz,&wtf2Gradz)) 
		{
			//Quadrature loop
			ierr = IGAElementBeginPoint(elemGradz,&pointGradz);CHKERRQ(ierr);
			ierr = IGAElementBeginPoint(elem_z,&point_z);CHKERRQ(ierr);

			while (IGAElementNextPoint(elemGradz,pointGradz))
			{
				if(pointGradz->atboundary==0 && point_z->atboundary==0)
				{
					IGAElementNextPoint(elem_z,point_z);

					ierr = IGAPointGetWorkVec(pointGradz,&FpointGradz);CHKERRQ(ierr);
					ierr = gradzF(pointGradz,point_z,FpointGradz,Z0Gradz,NULL);CHKERRQ(ierr);
					ierr = IGAPointAddVec(pointGradz,FpointGradz,FlocGradz);CHKERRQ(ierr);
				}
			}
			IGAElementNextPoint(elem_z,point_z);

			ierr = IGAElementEndPoint(elemGradz,&pointGradz);CHKERRQ(ierr);
			while (point_z->index != -1)
			{
				IGAElementNextPoint(elem_z,point_z);
			}
			ierr = IGAElementEndPoint(elem_z,&point_z);CHKERRQ(ierr);
		}
		ierr = IGAElementAssembleVec(elemGradz,FlocGradz,FGradz0);CHKERRQ(ierr);

	}
	IGANextElement(iga_z,elem_z);
	ierr = IGAEndElement(igaGradz,&elemGradz);CHKERRQ(ierr);
	ierr = IGAEndElement(iga_z,&elem_z);CHKERRQ(ierr);

	// Restore local vectors S0 and arrays
	ierr = IGARestoreLocalVecArray(iga_z,Z0,&localZ0Gradz,&arrayZ0Gradz);CHKERRQ(ierr);

	//Form system vector
	ierr = VecAssemblyBegin(FGradz0);CHKERRQ(ierr);
	ierr = VecAssemblyEnd  (FGradz0);CHKERRQ(ierr);

	ierr = KSPSolve(kspGradz,FGradz0,gradZ0);CHKERRQ(ierr);

	sprintf(nameGradz0,"%s%d%s","/Gradz-2d-",i,".dat");
	sprintf(pathGradz0,"%s%s",direct,nameGradz0);
	ierr = IGAWriteVec(igaGradz,gradZ0,pathGradz0);CHKERRQ(ierr);
//

/*
//System for general L2 projection (debug, but not really important, just for looking at fields)
	PetscPrintf(PETSC_COMM_WORLD,"\nSystem for general L2 projection starting \n\n");
	PetscPrintf(PETSC_COMM_WORLD,"\nCurrently projecting C*Ue \n\n");
	T=time(NULL);
	tm=*localtime(&T);
	PetscPrintf(PETSC_COMM_WORLD,"Current time is %02d:%02d:%02d \n",tm.tm_hour,tm.tm_min,tm.tm_sec);

	ierr = MatZeroEntries(KProj);CHKERRQ(ierr);						//This makes all non-zero elements of KProj 0.0, while keeping the sparse structure of the matrix
	ierr = VecZeroEntries(FProj);CHKERRQ(ierr);						//This sets all elements of the vector to 0

	//Get local vectors and arrays
	ierr = IGAGetLocalVecArray(iga_Al_hat,alInput,&localAlProj,&arrayAlProj);CHKERRQ(ierr);
	ierr = IGAGetLocalVecArray(iga_S,S0,&localSProj,&arraySProj);CHKERRQ(ierr);
	ierr = IGAGetLocalVecArray(igaVa,Va0,&localVaProj,&arrayVaProj);CHKERRQ(ierr);
	ierr = IGAGetLocalVecArray(iga_z,Z0,&localZ0Proj,&arrayZ0Proj);CHKERRQ(ierr);
	ierr = IGAGetLocalVecArray(iga_chi,chi,&localChiProj,&arrayChiProj);CHKERRQ(ierr);
	ierr = IGAGetLocalVecArray(iga_u,U,&localUProj,&arrayUProj);CHKERRQ(ierr);

	//Element loop
	ierr = IGABeginElement(igaProj,&elemProj);CHKERRQ(ierr);
	ierr = IGABeginElement(iga_Al_hat,&elem_Al_hat);CHKERRQ(ierr);
	ierr = IGABeginElement(iga_S,&elem_S);CHKERRQ(ierr);
	ierr = IGABeginElement(igaVa,&elemVa);CHKERRQ(ierr);
	ierr = IGABeginElement(iga_z,&elem_z);CHKERRQ(ierr);
	ierr = IGABeginElement(iga_chi,&elem_chi);CHKERRQ(ierr);
	ierr = IGABeginElement(iga_u,&elem_u);CHKERRQ(ierr);

	while (IGANextElement(igaProj,elemProj))
	{
		IGANextElement(iga_Al_hat,elem_Al_hat);
		IGANextElement(iga_S,elem_S);
		IGANextElement(igaVa,elemVa);
		IGANextElement(iga_z,elem_z);
		IGANextElement(iga_chi,elem_chi);
		IGANextElement(iga_u,elem_u);

		ierr = IGAElementGetWorkMat(elemProj,&KlocProj);CHKERRQ(ierr);
		ierr = IGAElementGetWorkVec(elemProj,&FlocProj);CHKERRQ(ierr);
		ierr = IGAElementGetValues(elem_Al_hat,arrayAlProj,&AlProj);CHKERRQ(ierr);
		ierr = IGAElementGetValues(elem_S,arraySProj,&SProj);CHKERRQ(ierr);
		ierr = IGAElementGetValues(elemVa,arrayVaProj,&VaProj);CHKERRQ(ierr);
		ierr = IGAElementGetValues(elem_z,arrayZ0Proj,&Z0Proj);CHKERRQ(ierr);
		ierr = IGAElementGetValues(elem_chi,arrayChiProj,&ChiProj);CHKERRQ(ierr);
		ierr = IGAElementGetValues(elem_u,arrayUProj,&UProj);CHKERRQ(ierr);

		//FormSystem loop
		while (IGAElementNextFormSystem(elemProj,&wtfProj,&wtf2Proj)) 
		{
			//Quadrature loop
			ierr = IGAElementBeginPoint(elemProj,&pointProj);CHKERRQ(ierr);
			ierr = IGAElementBeginPoint(elem_Al_hat,&point_Al_hat);CHKERRQ(ierr);
			ierr = IGAElementBeginPoint(elem_S,&point_S);CHKERRQ(ierr);
			ierr = IGAElementBeginPoint(elemVa,&pointVa);CHKERRQ(ierr);
			ierr = IGAElementBeginPoint(elem_z,&point_z);CHKERRQ(ierr);
			ierr = IGAElementBeginPoint(elem_chi,&point_chi);CHKERRQ(ierr);
			ierr = IGAElementBeginPoint(elem_u,&point_u);CHKERRQ(ierr);

			while (IGAElementNextPoint(elemProj,pointProj)) 
			{
				if(pointProj->atboundary==0 && point_Al_hat->atboundary==0 && point_S->atboundary==0 && pointVa->atboundary==0 && point_z->atboundary==0 && point_chi->atboundary==0)
				{
					IGAElementNextPoint(elem_Al_hat,point_Al_hat);
					IGAElementNextPoint(elem_S,point_S);
					IGAElementNextPoint(elemVa,pointVa);
					IGAElementNextPoint(elem_z,point_z);
					IGAElementNextPoint(elem_chi,point_chi);
					IGAElementNextPoint(elem_u,point_u);

					ierr = IGAPointGetWorkMat(pointProj,&KpointProj);CHKERRQ(ierr);
					ierr = IGAPointGetWorkVec(pointProj,&FpointProj);CHKERRQ(ierr);
					//proj(IGAPoint p,IGAPoint pAl,IGAPoint pS, IGAPoint pVa,IGAPoint pZu,IGAPoint pu,IGAPoint pChi,PetscReal *K,PetscReal *F,PetscReal *UAl,PetscReal *US,PetscReal *UVa,
					//	PetscReal *UZu,PetscReal *UChi,PetscReal *Uu,PetscReal Vs1,PetscReal Vs2,PetscReal normS,void *ctx)
					ierr = proj(pointProj,point_Al_hat,point_S,pointVa,point_z,point_u,point_chi,KpointProj,FpointProj,AlProj,SProj,VaProj,Z0Proj,ChiProj,UProj,Int1a,Int2a,S_norm,NULL);CHKERRQ(ierr);
					ierr = IGAPointAddMat(pointProj,KpointProj,KlocProj);CHKERRQ(ierr);
					ierr = IGAPointAddVec(pointProj,FpointProj,FlocProj);CHKERRQ(ierr);
				}
			}
			
			IGAElementNextPoint(elem_Al_hat,point_Al_hat);
			IGAElementNextPoint(elem_S,point_S);
			IGAElementNextPoint(elemVa,pointVa);
			IGAElementNextPoint(elem_z,point_z);
			IGAElementNextPoint(elem_chi,point_chi);
			IGAElementNextPoint(elem_u,point_u);

			ierr = IGAElementEndPoint(elemProj,&pointProj);CHKERRQ(ierr);
			ierr = IGAElementEndPoint(elem_Al_hat,&point_Al_hat);CHKERRQ(ierr);
			ierr = IGAElementEndPoint(elem_S,&point_S);CHKERRQ(ierr);
			ierr = IGAElementEndPoint(elemVa,&pointVa);CHKERRQ(ierr);
			ierr = IGAElementEndPoint(elem_z,&point_z);CHKERRQ(ierr);
			ierr = IGAElementEndPoint(elem_chi,&point_chi);CHKERRQ(ierr);
			ierr = IGAElementEndPoint(elem_u,&point_u);CHKERRQ(ierr);
		}
		ierr = IGAElementAssembleMat(elemProj,KlocProj,KProj);CHKERRQ(ierr);
		ierr = IGAElementAssembleVec(elemProj,FlocProj,FProj);CHKERRQ(ierr);
	}

	IGANextElement(iga_Al_hat,elem_Al_hat);
	IGANextElement(iga_S,elem_S);
	IGANextElement(igaVa,elemVa);
	IGANextElement(iga_z,elem_z);
	IGANextElement(iga_chi,elem_chi);
	IGANextElement(iga_u,elem_u);

	ierr = IGAEndElement(igaProj,&elemProj);CHKERRQ(ierr);
	ierr = IGAEndElement(iga_Al_hat,&elem_Al_hat);CHKERRQ(ierr);
	ierr = IGAEndElement(iga_S,&elem_S);CHKERRQ(ierr);
	ierr = IGAEndElement(igaVa,&elemVa);CHKERRQ(ierr);
	ierr = IGAEndElement(iga_z,&elem_z);CHKERRQ(ierr);
	ierr = IGAEndElement(iga_chi,&elem_chi);CHKERRQ(ierr);
	ierr = IGAEndElement(iga_u,&elem_u);CHKERRQ(ierr);

	// Restore local vectors S0 and arrays
	ierr = IGARestoreLocalVecArray(iga_Al_hat,alInput,&localAlProj,&arrayAlProj);CHKERRQ(ierr);
	ierr = IGARestoreLocalVecArray(iga_S,S0,&localSProj,&arraySProj);CHKERRQ(ierr);
	ierr = IGARestoreLocalVecArray(igaVa,Va0,&localVaProj,&arrayVaProj);CHKERRQ(ierr);
	ierr = IGARestoreLocalVecArray(iga_z,Z0,&localZ0Proj,&arrayZ0Proj);CHKERRQ(ierr);
	ierr = IGARestoreLocalVecArray(iga_chi,chi,&localChiProj,&arrayChiProj);CHKERRQ(ierr);
	ierr = IGARestoreLocalVecArray(iga_u,U,&localUProj,&arrayUProj);CHKERRQ(ierr);

	//Form system matrix and vector
	ierr = MatAssemblyBegin(KProj,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyEnd  (KProj,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	
	ierr = VecAssemblyBegin(FProj);CHKERRQ(ierr);
	ierr = VecAssemblyEnd  (FProj);CHKERRQ(ierr);

	ierr = KSPSolve(kspProj,FProj,proj0);CHKERRQ(ierr);

	sprintf(nameProj,"%s%d%s","/Proj-2d-",i,".dat");
	sprintf(pathProj,"%s%s",direct,nameProj);
	ierr = IGAWriteVec(igaProj,proj0,pathProj);CHKERRQ(ierr);
//
*/
}

/*
//Destroy all objects not needed anymore (Better to do it here in case different codes call the same IGA, move if memory is a problem)
	ierr = IGADestroy(&iga_Al_hat); CHKERRQ(ierr);
	ierr = IGADestroy(&iga_z); CHKERRQ(ierr);
	ierr = IGADestroy(&iga_chi); CHKERRQ(ierr);
	ierr = KSPDestroy(&kspZdot);CHKERRQ(ierr);
	ierr = MatDestroy(&KZdot);CHKERRQ(ierr);
	ierr = VecDestroy(&FZdot);CHKERRQ(ierr);
//
*/

ierr = PetscFinalize();CHKERRQ(ierr);

return 0;
}