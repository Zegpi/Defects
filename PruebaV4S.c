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
	#define __FUNCT__ "curlChiS"
	PetscErrorCode curlChiS(IGAPoint p,IGAPoint pS,PetscReal *K,PetscReal *F,PetscReal *U,void *ctx)
	{
		const PetscReal *N0,(*N1)[2];
		IGAPointGetShapeFuns(p,0,(const PetscReal**)&N0);
		IGAPointGetShapeFuns(p,1,(const PetscReal**)&N1);

		PetscInt a,b,u,w,i,j,k,m,nen=p->nen, dof=p->dof;

		//PetscReal S[8];																		//Create array to recieve S
		PetscReal dS[8][2];																		//Create array to recieve dS
		//IGAPointFormValue(pS,U,&S[0]);														//This fills the values
		IGAPointFormGrad (pS,U,&dS[0][0]);														//This fills the values of the derivatives

		//PetscReal fullS[3][3][3]={0};
		PetscReal fulldS[3][3][3][3]={0};

		//fullS[0][0][0]=S[0]; fullS[0][0][1]=S[1]; fullS[0][1][0]=S[2]; fullS[0][1][1]=S[3];		//Expand S to full 3rd order form, only non-zero elements
		//fullS[1][0][0]=S[4]; fullS[1][0][1]=S[5]; fullS[1][1][0]=S[6]; fullS[1][1][1]=S[7];

		//Same for gradS
		fulldS[0][0][0][0]=dS[0][0]; fulldS[0][0][0][1]=dS[0][1]; 
		fulldS[0][0][1][0]=dS[1][0]; fulldS[0][0][1][1]=dS[1][1]; 
		fulldS[0][1][0][0]=dS[2][0]; fulldS[0][1][0][1]=dS[2][1]; 
		fulldS[0][1][1][0]=dS[3][0]; fulldS[0][1][1][1]=dS[3][1];
		fulldS[1][0][0][0]=dS[4][0]; fulldS[1][0][0][1]=dS[4][1]; 
		fulldS[1][0][1][0]=dS[5][0]; fulldS[1][0][1][1]=dS[5][1]; 
		fulldS[1][1][0][0]=dS[6][0]; fulldS[1][1][0][1]=dS[6][1]; 
		fulldS[1][1][1][0]=dS[7][0]; fulldS[1][1][1][1]=dS[7][1];

		PetscReal (*KchiS)[dof][nen][dof] = (typeof(KchiS)) K;
		PetscReal (*FchiS)[dof] = (PetscReal (*)[dof]) F;

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

					for(i=0;i<3;i++)
					{
						for(j=0;j<3;j++)
						{
							for(k=0;k<3;k++)
							{
								for(m=0;m<3;m++)
								{
									FchiS[a][u]+=fulldS[i][j][k][m]*(dv[i][j][k][m]-dv[i][j][m][k]);
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
	PetscErrorCode curlChiSF(IGAPoint p,IGAPoint pS,PetscReal *F,PetscReal *U,void *ctx)
	{
		const PetscReal *N0,(*N1)[2];
		IGAPointGetShapeFuns(p,0,(const PetscReal**)&N0);
		IGAPointGetShapeFuns(p,1,(const PetscReal**)&N1);

		PetscInt a,u,i,j,k,m,nen=p->nen, dof=p->dof;

		//PetscReal S[8];																		//Create array to receive S
		PetscReal dS[8][2];																		//Create array to receive dS
		//IGAPointFormValue(pS,U,&S[0]);														//This fills the values
		IGAPointFormGrad (pS,U,&dS[0][0]);														//This fills the values of the derivatives

		//PetscReal fullS[3][3][3]={0};
		PetscReal fulldS[3][3][3][3]={0};

		//fullS[0][0][0]=S[0]; fullS[0][0][1]=S[1]; fullS[0][1][0]=S[2]; fullS[0][1][1]=S[3];		//Expand S to full 3rd order form, only non-zero elements
		//fullS[1][0][0]=S[4]; fullS[1][0][1]=S[5]; fullS[1][1][0]=S[6]; fullS[1][1][1]=S[7];

		//Same for gradS
		fulldS[0][0][0][0]=dS[0][0]; fulldS[0][0][0][1]=dS[0][1]; 
		fulldS[0][0][1][0]=dS[1][0]; fulldS[0][0][1][1]=dS[1][1]; 
		fulldS[0][1][0][0]=dS[2][0]; fulldS[0][1][0][1]=dS[2][1]; 
		fulldS[0][1][1][0]=dS[3][0]; fulldS[0][1][1][1]=dS[3][1];
		fulldS[1][0][0][0]=dS[4][0]; fulldS[1][0][0][1]=dS[4][1]; 
		fulldS[1][0][1][0]=dS[5][0]; fulldS[1][0][1][1]=dS[5][1]; 
		fulldS[1][1][0][0]=dS[6][0]; fulldS[1][1][0][1]=dS[6][1]; 
		fulldS[1][1][1][0]=dS[7][0]; fulldS[1][1][1][1]=dS[7][1];

		PetscReal (*FchiS)[dof] = (PetscReal (*)[dof]) F;

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

					for(i=0;i<3;i++)
					{
						for(j=0;j<3;j++)
						{
							for(k=0;k<3;k++)
							{
								for(m=0;m<3;m++)
								{
									FchiS[a][u]+=fulldS[i][j][k][m]*(dv[i][j][k][m]-dv[i][j][m][k]);
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

/*
//Helmholtz decomposition of S, grad part
	#undef  __FUNCT__
	#define __FUNCT__ "gradZS"
	PetscErrorCode gradZS(IGAPoint p,IGAPoint pS,PetscReal *K,PetscReal *F,PetscReal *U,void *ctx)
	{

		const PetscReal *N0,(*N1)[2];
		IGAPointGetShapeFuns(p,0,(const PetscReal**)&N0);
		IGAPointGetShapeFuns(p,1,(const PetscReal**)&N1);

		PetscInt a,b,u,w,i,j,k,nen=p->nen, dof=p->dof;

		//PetscReal S[8];																		//Create array to receive Alfa
		PetscReal dS[8][2];																	//Create array to receive dAlfa
		//IGAPointFormValue(pS,U,&S[0]);															//This fills the values
		IGAPointFormGrad (pS,U,&dS[0][0]);														//This fills the values of the derivatives

		//PetscReal fullS[3][3][3]={0};
		PetscReal fulldS[3][3][3][3]={0};

		//fullS[0][0][0]=S[0]; fullS[0][0][1]=S[1]; 											//Expand S to full 3rd order form, only non-zero elements
		//fullS[0][1][0]=S[2]; fullS[0][1][1]=S[3];		
		//fullS[1][0][0]=S[4]; fullS[1][0][1]=S[5]; 
		//fullS[1][1][0]=S[6]; fullS[1][1][1]=S[7];

		//Same for gradS
		fulldS[0][0][0][0]=dS[0][0]; fulldS[0][0][0][1]=dS[0][1]; 
		fulldS[0][0][1][0]=dS[1][0]; fulldS[0][0][1][1]=dS[1][1]; 
		fulldS[0][1][0][0]=dS[2][0]; fulldS[0][1][0][1]=dS[2][1]; 
		fulldS[0][1][1][0]=dS[3][0]; fulldS[0][1][1][1]=dS[3][1];
		fulldS[1][0][0][0]=dS[4][0]; fulldS[1][0][0][1]=dS[4][1]; 
		fulldS[1][0][1][0]=dS[5][0]; fulldS[1][0][1][1]=dS[5][1]; 
		fulldS[1][1][0][0]=dS[6][0]; fulldS[1][1][0][1]=dS[6][1]; 
		fulldS[1][1][1][0]=dS[7][0]; fulldS[1][1][1][1]=dS[7][1];

		PetscReal (*KZS)[dof][nen][dof] = (typeof(KZS)) K;
		PetscReal (*FZS)[dof] = (PetscReal (*)[dof]) F;

		if (p->atboundary)
		{
			return 0;
		}
		else
		{
			for (a=0; a<nen; a++)
			{
				PetscReal Na =N0[a];
				PetscReal Na_x = N1[a][0];
				PetscReal Na_y = N1[a][1];
				for (u=0; u<dof; u++) 
				{
					PetscReal v[3][3]={0};
					PetscReal dv[3][3][3]={0};

					if(u==0){
						v[0][0]=Na;
						dv[0][0][0]=Na_x; dv[0][0][1]=Na_y;}
					if(u==1){
						v[0][1]=Na;
						dv[0][1][0]=Na_x; dv[0][1][1]=Na_y;}
					if(u==2){
						v[1][0]=Na;
						dv[1][0][0]=Na_x; dv[1][0][1]=Na_y;}
					if(u==3){
						v[1][1]=Na;
						dv[1][1][0]=Na_x; dv[1][1][1]=Na_y;}

					for (b=0; b<nen; b++)
					{
						PetscReal Nb_x = N1[b][0];
						PetscReal Nb_y = N1[b][1];
						for (w=0; w<dof; w++)
						{
							PetscReal dZS[3][3][3]={0};

							if(w==0){
								dZS[0][0][0]=Nb_x; dZS[0][0][1]=Nb_y;}
							if(w==1){
								dZS[0][1][0]=Nb_x; dZS[0][1][1]=Nb_y;}
							if(w==2){
								dZS[1][0][0]=Nb_x; dZS[1][0][1]=Nb_y;}
							if(w==3){
								dZS[1][1][0]=Nb_x; dZS[1][1][1]=Nb_y;}

							for(i=0;i<3;i++)
							{
								for(j=0;j<3;j++)
								{
									for(k=0;k<3;k++)
									{
											KZS[a][u][b][w]+=dZS[i][j][k]*dv[i][j][k];
									}
								}
							}
						}
					}

					for(i=0;i<3;i++)
					{
						for(j=0;j<3;j++)
						{
							for(k=0;k<3;k++)
							{
								FZS[a][u]+=-fulldS[i][j][k][k] * v[i][j];
							}
						}
					}
				}
			}
		}

		return 0;
	}
//
*/

//L2 projection of Al(0)-Sp:X
	#undef  __FUNCT__
	#define __FUNCT__ "L2ProjectionAlSp"
	PetscErrorCode L2ProjectionAlSp(IGAPoint p,IGAPoint pS,PetscReal *K,PetscReal *F,PetscReal *U,void *ctx)
	{
		//This function just creates -S:X, alpha is read from file and added to it on the function call routine
		//This function should be rewritten to be more legible
		if (p->atboundary)
		{
			return 0.0;																		//Pure Neumann condition 
		}

		//Modified to add -Sp:X to alpha
		//Definition of alternating tensor
		const PetscReal e[3][3][3]=
		{
			{{0.0,0.0,0.0},{0.0,0.0,1.0},{0.0,-1.0,0.0}},
			{{0.0,0.0,-1.0},{0.0,0.0,0.0},{1.0,0.0,0.0}},
			{{0.0,1.0,0.0},{-1.0,0.0,0.0},{0.0,0.0,0.0}}
		};

		PetscReal Sp[8];
		IGAPointFormValue(pS,U,&Sp[0]);														//This fills the values of S
		PetscReal fullSp[3][3][3]={0};

		//This kills spurious gradients due to meshes, could be removed
		fullSp[0][0][0]=round(1.0e12*Sp[0])/1.0e12; fullSp[0][0][1]=round(1.0e12*Sp[1])/1.0e12; fullSp[0][1][0]=round(1.0e12*Sp[2])/1.0e12; fullSp[0][1][1]=round(1.0e12*Sp[3])/1.0e12;		//Expand S to full 3rd order form, only non-zero elements
		fullSp[1][0][0]=round(1.0e12*Sp[4])/1.0e12; fullSp[1][0][1]=round(1.0e12*Sp[5])/1.0e12; fullSp[1][1][0]=round(1.0e12*Sp[6])/1.0e12; fullSp[1][1][1]=round(1.0e12*Sp[7])/1.0e12;

		PetscInt a,b,i,j,k,l;
		PetscInt nen = p->nen;																//Number of shape functions, 9 in this case.
		PetscInt dim = p->dim;																//Number of spatial dimensions
		PetscInt dof = p->dof;

		PetscReal SpX[3][3]={0};																//Stores S:X

		for(i=0; i<3;i++)
		{
			for (j=0; j<3; j++)
			{
				for (k=0; k<3; k++)
				{
					for (l=0; l<3; l++)
					{
						SpX[i][j]=SpX[i][j]+fullSp[i][k][l]*e[j][k][l];
					}
				}
			}
		}

		PetscReal sp[2]={0.0};
		sp[0]=SpX[0][2]; sp[1]=SpX[1][2];

		PetscReal x[dim];																	//Vector of reals, size equal to problem's dimension
		IGAPointFormGeomMap(p,x);															//Fills x with the coordinates of p, Gauss's point

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
				FF[a][i] = N[a]*(-sp[i]);
			}
		}
		return 0;
	}
//

//Helmholtz decomposition of Ûp (same as Ûe), curl part
	#undef  __FUNCT__
	#define __FUNCT__ "curlChiU"
	PetscErrorCode curlChiU(IGAPoint p,IGAPoint pAl, PetscReal *K,PetscReal *F,PetscReal *U, void *ctx)
	{
		//Receiving Z but not using, check if it's necessary

		const PetscReal *N0,(*N1)[2];
		IGAPointGetShapeFuns(p,0,(const PetscReal**)&N0);
		IGAPointGetShapeFuns(p,1,(const PetscReal**)&N1);

		PetscInt a,b,u,w,i,j,k,l,nen=p->nen, dof=p->dof;

		//Definition of alternating tensor
		const PetscReal e[3][3][3]=
		{
			{{0.0,0.0,0.0},{0.0,0.0,1.0},{0.0,-1.0,0.0}},
			{{0.0,0.0,-1.0},{0.0,0.0,0.0},{1.0,0.0,0.0}},
			{{0.0,1.0,0.0},{-1.0,0.0,0.0},{0.0,0.0,0.0}}
		};

		PetscReal alfa[2];																		//Create array to recieve Alfa
		IGAPointFormValue(pAl,U,&alfa[0]);														//This fills the values

		PetscReal fullAlfa[3][3]={0};
		fullAlfa[0][2]=alfa[0]; fullAlfa[1][2]=alfa[1];												//Expand Alfa to full 2nd order form, only non-zero elements

		PetscReal (*KchiU)[dof][nen][dof] = (typeof(KchiU)) K;
		PetscReal (*FchiU)[dof] = (PetscReal (*)[dof]) F;

		if (p->atboundary)
		{
			return 0;
		}
		else
		{
			for (a=0; a<nen; a++)
			{
				//PetscReal Na =N0[a];
				PetscReal Na_x = N1[a][0];
				PetscReal Na_y = N1[a][1];
				for (u=0; u<dof; u++) 
				{
					//PetscReal v[3][3]={0};
					PetscReal dv[3][3][3]={0};

					if(u==0){
						//v[0][0]=Na;
						dv[0][0][0]=Na_x; dv[0][0][1]=Na_y;}
					if(u==1){
						//v[0][1]=Na;
						dv[0][1][0]=Na_x; dv[0][1][1]=Na_y;}
					if(u==2){
						//v[1][0]=Na;
						dv[1][0][0]=Na_x; dv[1][0][1]=Na_y;}
					if(u==3){
						//v[1][1]=Na;
						dv[1][1][0]=Na_x; dv[1][1][1]=Na_y;}

					for (b=0; b<nen; b++)
					{
						PetscReal Nb_x = N1[b][0];
						PetscReal Nb_y = N1[b][1];
						for (w=0; w<dof; w++)
						{
							PetscReal dChiU[3][3][3]={0};

							if(w==0){
								dChiU[0][0][0]=Nb_x; dChiU[0][0][1]=Nb_y;}
							if(w==1){
								dChiU[0][1][0]=Nb_x; dChiU[0][1][1]=Nb_y;}
							if(w==2){
								dChiU[1][0][0]=Nb_x; dChiU[1][0][1]=Nb_y;}
							if(w==3){
								dChiU[1][1][0]=Nb_x; dChiU[1][1][1]=Nb_y;}

							for(i=0;i<3;i++)
							{
								for(j=0;j<3;j++)
								{
									for(k=0;k<3;k++)
									{
											KchiU[a][u][b][w]+=dChiU[i][j][k]*dv[i][j][k]-dChiU[i][j][k]*dv[i][k][j]+dChiU[i][j][j]*dv[i][k][k];
									}
								}
							}
						}
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

	#undef  __FUNCT__
	#define __FUNCT__ "curlChiUF"
	PetscErrorCode curlChiUF(IGAPoint p,IGAPoint pAl,PetscReal *F,PetscReal *U, void *ctx)
	{
		//Generates right hand side only. Faster for iterations, as K is constant 

		const PetscReal *N0,(*N1)[2];
		IGAPointGetShapeFuns(p,0,(const PetscReal**)&N0);
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
				//PetscReal Na =N0[a];
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

//System for z(0) [make extra sure that sigma=C*(-grad(z)-chi) ]
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

		//////////////Delete this part later, loads should come from appCtx or be 0
		PetscReal f[3]={0.0, 0.0, 0.0};			//Distributed load in body
		PetscReal g[3]={0.0, 0.0, 0.0};			//Boundary load (applied wherever is defined in the p->atboundary block)
		/////////////////////////////////////

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

		PetscReal n[3]={0};
		PetscReal v[3]={0};
		PetscReal dv[3][3]={0};
		PetscReal dz[3][3]={0};

		if (p->atboundary)
		{
			PetscReal Sborde[3][3]={0};

			//No stress in boundary
			Sborde[0][0]=0.0;
			Sborde[0][1]=0.0;
			Sborde[1][0]=0.0;
			Sborde[1][1]=0.0;
			Sborde[0][2]=0.0; Sborde[1][2]=0.0; Sborde[2][0]=0.0; Sborde[2][1]=0.0; Sborde[2][2]=0.0;

			PetscInt dir  = p->boundary_id / 2;
			PetscInt side = p->boundary_id % 2;

			for (a=0; a<nen; a++)
			{
				PetscReal Na   = N0[a];

				for (i=0; i<dof; i++)
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
							for (j=0; j<3; j++)
							{
								for(k=0; k<3;k++)
								{
									Feq[a][i]+=Sborde[j][k]*n[k]*v[j];
								}
							}
						}
						if(side==1)
						{
							n[0]=1.0; n[1]=0.0; n[2]=0.0;
							for (j=0; j<3; j++)
							{
								for(k=0; k<3;k++)
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
							for (j=0; j<3; j++)
							{
								for(k=0; k<3;k++)
								{
									Feq[a][i]+=Sborde[j][k]*n[k]*v[j];
								}
							}
						}
						if(side==1)
						{
							n[0]=0.0; n[1]=1.0; n[2]=0.0;
							for (j=0; j<3; j++)
							{
								for(k=0; k<3;k++)
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

						for (j=0; j<dof; j++)
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
									for (u=0; u<3; u++)
									{
										for (w=0; w<3; w++)
										{
											Keq[a][i][b][j]+=0.5*(C[k][l][u][w]*(-dz[u][w])+C[l][k][u][w]*(-dz[u][w]))*0.5*(dv[k][l]+dv[l][k]);
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
				PetscReal Na=N0[a];
				PetscReal Na_x=N1[a][0];		PetscReal Na_y=N1[a][1];

				for (i=0; i<dof; i++)
				{
					if (i==0)
					{
						v[0]=Na; 	   v[1]=0.0; 	  v[2]=0.0;
						dv[0][0]=Na_x; dv[0][1]=Na_y; dv[0][2]=0.0;
						dv[1][0]=0.0;  dv[1][1]=0.0;  dv[1][2]=0.0;
						dv[2][0]=0.0;  dv[2][1]=0.0;  dv[2][2]=0.0;

					}
					else if (i==1)
					{
						v[0]=0.0; 	   v[1]=Na; 	  v[2]=0.0;
						dv[0][0]=0.0;  dv[0][1]=0.0;  dv[0][2]=0.0;
						dv[1][0]=Na_x; dv[1][1]=Na_y; dv[1][2]=0.0;
						dv[2][0]=0.0;  dv[2][1]=0.0;  dv[2][2]=0.0;

					}

					Feq[a][i] = 0.0;
					for(k=0; k<3; k++)
					{
						Feq[a][i]+=f[k]*v[k];
					}

					for (k=0; k<3; k++)
					{
						for (l=0; l<3; l++)
						{
							for (u=0; u<3; u++)
							{
								for (w=0; w<3; w++)
								{
									Feq[a][i]+=-0.5*(C[k][l][u][w]*(-fullChi[u][w])+C[l][k][u][w]*(-fullChi[u][w]))*0.5*(dv[k][l]+dv[l][k]);
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

//System for \dot{z} (debug)
	#undef  __FUNCT__
	#define __FUNCT__ "ZdotSystem"
	PetscErrorCode ZdotSystem(IGAPoint p,IGAPoint pAl,IGAPoint pVa,IGAPoint pS,PetscReal *K,PetscReal *F,PetscReal *UAl,PetscReal *UVa,PetscReal *US, void *ctx)
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
		//Restore this later
		//fullVs[0]=Vs[0]; fullVs[1]=Vs[1];
		//This is for testing
		fullVs[0]=0.0; fullVs[1]=-1.0;

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
					for (k=0; k<3; k++)
					{
						for (l=0; l<3; l++)
						{
							for (u=0; u<3; u++)
							{
								for (w=0; w<3; w++)
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

//System for u (debug)
	#undef  __FUNCT__
	#define __FUNCT__ "Usys"
	PetscErrorCode Usys(IGAPoint p, IGAPoint pChi, IGAPoint pZ, PetscReal *K, PetscReal *F, PetscReal *UChi, PetscReal *UZ, void *ctx)
	{
		const PetscReal *N0,(*N1)[2];
		IGAPointGetShapeFuns(p,0,(const PetscReal**)&N0);
		IGAPointGetShapeFuns(p,1,(const PetscReal**)&N1);

		PetscInt a,b,i,j,u,w,k,l,nen=p->nen, dof=p->dof;

		//Change for G=1
		const PetscReal nu=0.33;
		const PetscReal mu=1.0;
		const PetscReal lambda=2.0*mu*nu/(1.0-2.0*nu);

		PetscReal Chi0[4];																	//Assign chi to a vector
		IGAPointFormValue(pChi,UChi,&Chi0[0]);

		PetscReal fullChi[3][3]={0};
		fullChi[0][0]=Chi0[0]; 	fullChi[0][1]=Chi0[1];
		fullChi[1][0]=Chi0[2]; 	fullChi[1][1]=Chi0[3];
	
		PetscReal dz0[2][2];																//Same for partial derivatives
		IGAPointFormGrad (pZ,UZ,&dz0[0][0]);

		PetscReal full_dz[3][3]={0};
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
		PetscReal du[3][3]={0};
		PetscReal n[3]={0};

		PetscReal (*Keq)[dof][nen][dof] = (typeof(Keq)) K;
		PetscReal (*Feq)[dof] = (PetscReal (*)[dof])F;

		if (p->atboundary)
		{
			PetscReal Sborde[3][3]={0};

			//Stress in boundary
			Sborde[0][0]=0.0;
			Sborde[0][1]=1.0;
			Sborde[1][0]=1.0;
			Sborde[1][1]=0.0;
			Sborde[0][2]=0.0; Sborde[1][2]=0.0; Sborde[2][0]=0.0; Sborde[2][1]=0.0; Sborde[2][2]=0.0;

			PetscInt dir  = p->boundary_id / 2;
			PetscInt side = p->boundary_id % 2;

			for (a=0; a<nen; a++)
			{
				PetscReal Na   = N0[a];

				for (i=0; i<dof; i++)
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
							for (j=0; j<3; j++)
							{
								for(k=0; k<3;k++)
								{
									Feq[a][i]+=Sborde[j][k]*n[k]*v[j];
								}
							}
						}
						if(side==1)
						{
							n[0]=1.0; n[1]=0.0; n[2]=0.0;
							for (j=0; j<3; j++)
							{
								for(k=0; k<3;k++)
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
							for (j=0; j<3; j++)
							{
								for(k=0; k<3;k++)
								{
									Feq[a][i]+=Sborde[j][k]*n[k]*v[j];
								}
							}
						}
						if(side==1)
						{
							n[0]=0.0; n[1]=1.0; n[2]=0.0;
							for (j=0; j<3; j++)
							{
								for(k=0; k<3;k++)
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

						for (j=0; j<dof; j++)
						{
							if (j==0)
							{
								du[0][0]=Nb_x; du[0][1]=Nb_y; du[0][2]=0.0;
								du[1][0]=0.0;  du[1][1]=0.0;  du[1][2]=0.0;
								du[2][0]=0.0;  du[2][1]=0.0;  du[2][2]=0.0;
							}
							else if(j==1)
							{
								du[0][0]=0.0;  du[0][1]=0.0;  du[0][2]=0.0;
								du[1][0]=Nb_x; du[1][1]=Nb_y; du[1][2]=0.0;
								du[2][0]=0.0;  du[2][1]=0.0;  du[2][2]=0.0;
							}

							Keq[a][i][b][j]=0.0;
							for (k=0; k<3; k++)
							{
								for (l=0; l<3; l++)
								{
									for (u=0; u<3; u++)
									{
										for(w=0; w<3; w++)
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

			for(a=0;a<nen;a++)
			{
				//PetscReal Na=N0[a];
				PetscReal Na_x = N1[a][0];			PetscReal Na_y = N1[a][1];

				for (i=0; i<dof; i++)
				{
					if (i==0)
					{
						//v[0]=Na; 	   v[1]=0.0; 	  v[2]=0.0;
						dv[0][0]=Na_x; dv[0][1]=Na_y; dv[0][2]=0.0;
						dv[1][0]=0.0;  dv[1][1]=0.0;  dv[1][2]=0.0;
						dv[2][0]=0.0;  dv[2][1]=0.0;  dv[2][2]=0.0;
					}
					else if (i==1)
					{
						//v[0]=0.0; 	   v[1]=Na; 	  v[2]=0.0;
						dv[0][0]=0.0;  dv[0][1]=0.0;  dv[0][2]=0.0;
						dv[1][0]=Na_x; dv[1][1]=Na_y; dv[1][2]=0.0;
						dv[2][0]=0.0;  dv[2][1]=0.0;  dv[2][2]=0.0;
					}

					Feq[a][i] = 0.0;
					for (k=0; k<3; k++)
					{
						for (l=0; l<3; l++)
						{
							for (u=0; u<3; u++)
							{
								for (w=0; w<3; w++)
								{
									Feq[a][i]+=C[k][l][u][w]*(full_dz[u][w]+fullChi[u][w])*dv[k][l];
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
	PetscErrorCode VS(IGAPoint p,IGAPoint pChi,IGAPoint pZu,IGAPoint pS,PetscReal *K,PetscReal *F, PetscReal *UChi,PetscReal *UZu,PetscReal *US,void *ctx)		//When other fields like Pi or S (etc.) come into this, declare a IGAPoint pPi or pS and a new PetscReal *UPi or *US for each
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
										FVS[a][i]+=C[j][k][l][m]*(-full_dz[l][m]-fullChi[l][m])*fullS[j][k][w]*v[w];
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
										FVS[a][i]+=C[j][k][l][m]*(full_du[l][m]-full_dz[l][m]-fullChi[l][m])*fullS[j][k][w]*v[w];
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
	//IntValpha(pointVa,pointAlp,FpointVaInt,PointInt1,PointInt2,Va0Values,Al1Va,NULL);CHKERRQ(ierr);		//When other fields like Pi or S (etc.) come into this, declare a IGAPoint pPi or pS and a new PescReal *UPi or *US for each
	PetscErrorCode Int_Xi_Vs(IGAPoint pV,IGAPoint pS,PetscReal *FInt1a,PetscReal *FInt2a,PetscReal *FInt,PetscReal *VS,PetscReal *US,void *ctx)		//When other fields like Pi or S (etc.) come into this, declare a IGAPoint pPi or pS and a new PetscReal *UPi or *US for each
	{
		PetscInt dof=pV->dof;
		PetscInt i;

		PetscReal Vs[2];																	//Create array to receive Alfa
		IGAPointFormValue(pV,VS,&Vs[0]);													//This fills the values

		PetscReal S[8];																		//Create array to receive Alfa
		IGAPointFormValue(pS,US,&S[0]);														//This fills the values

		PetscReal xi=0.0;
		PetscInt ind=0;
		for (i=0;i<8;i++)
		{
			if (S[i]>0)
			{
				ind=i;
				xi=1.0;
			}	
		}
		
		if (S[0]<0 || S[1]<0 || S[2]<0 || S[3]<0 || S[4]<0 || S[5]<0 || S[6]<0 || S[7]<0)
		{
			xi=1.0;
		}

		PetscReal (*FI1a)[dof] = (PetscReal (*)[dof])FInt1a;		//This vector will have just VS[0]
		PetscReal (*FI2a)[dof] = (PetscReal (*)[dof])FInt2a;		//This vector will have just VS[1]
		PetscReal (*FI)[dof]   = (PetscReal (*)[dof])FInt;			//This vector will have xi(S), so we can integrate it.

		if (pV->atboundary)
		{
			return 0;
		}
		else
		{
			FI1a[0][0]+=xi*Vs[0];
			FI2a[0][1]+=xi*Vs[1];

			FI[0][0]+=S[ind];
			FI[0][1]+=0.0;

		}
		return 0;
	}

	#undef  __FUNCT__
	#define __FUNCT__ "ProjVS"
	//ierr = ProjVS(pointVs,pointS,VsSmoothPoint,SVs,Int1a,Int2a,IntS,absSmax,NULL);CHKERRQ(ierr);
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
			if (S[i]>1.0e-6 || S[i]<-1.0e-6)
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
	//PetscErrorCode Stress(IGAPoint p,IGAPoint pU, IGAPoint pHs,IGAPoint pChi,IGAPoint pZu,PetscReal *K,PetscReal *F,PetscReal *U,PetscReal *HS, PetscReal *Chi,PetscReal *Zu,void *ctx)		//When other fields like Pi or S (etc.) come into this, declare a IGAPoint pPi or pS and a new PescReal *UPi or *US for each
	PetscErrorCode Valpha(IGAPoint p,IGAPoint pChi,IGAPoint pZu,IGAPoint pAl,IGAPoint pS,PetscReal *K,PetscReal *F, PetscReal *Chi,PetscReal *Zu,PetscReal *US,PetscReal *UAl,void *ctx)		//When other fields like Pi or S (etc.) come into this, declare a IGAPoint pPi or pS and a new PetscReal *UPi or *US for each
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

		PetscReal alfa[2];																	//Create array to receive Alfa
		IGAPointFormValue(pAl,UAl,&alfa[0]);												//This fills the values

		PetscReal chi0[4];																	//Array to contain the vector chi(0)
		IGAPointFormValue(pChi,Chi,&chi0[0]);												//Assign chi to its container

		PetscReal d_Z0[2][2];																//Same for its gradient
		IGAPointFormGrad (pZu,Zu,&d_Z0[0][0]);												//Same for the gradient

		//Inflate stored vectors to full tensor form
		PetscReal fullAlfa[3][3]={0};
		fullAlfa[0][2]=alfa[0]; fullAlfa[1][2]=alfa[1];										//Expand Alfa to full 2nd order form, only non-zero elements

		//The four non-zero components of Chi are stored as a vector, restore them to an array with the correct indexing for value and derivative
		PetscReal fullChi[3][3]={0};
		fullChi[0][0]=chi0[0]; 	fullChi[0][1]=chi0[1];
		fullChi[1][0]=chi0[2]; 	fullChi[1][1]=chi0[3];

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
											FCS[a][i]+=C[j][k][l][m]*(-fulld_z[l][m]-fullChi[l][m])*e[k][c][d]*fullAlfa[j][c]*v[d];
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

		PetscReal chi0[4];																	//Array to contain the vector chi(0)
		IGAPointFormValue(pChi,Chi,&chi0[0]);												//Assign chi to its container

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
		fullChi[0][0]=chi0[0]; 	fullChi[0][1]=chi0[1];
		fullChi[1][0]=chi0[2]; 	fullChi[1][1]=chi0[3];

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
	//Sdot(pointSdot,pointVs,pointS,KpointSdot,FpointSdot,VsSdot,S0Sdot,NULL);CHKERRQ(ierr);
	PetscErrorCode SdotFunc(IGAPoint p, IGAPoint pVs, IGAPoint pS, PetscReal *K, PetscReal *F, PetscReal *UVs, PetscReal *US, void *ctx)
	{
		PetscInt a,b,i,j,k,l,u,w,dof=p->dof,nen=p->nen;	

		const PetscReal *N0,(*N1)[2];
		IGAPointGetShapeFuns(p,0,(const PetscReal**)&N0);									//Value of the shape functions
		IGAPointGetShapeFuns(p,1,(const PetscReal**)&N1);									//Derivatives of the shape functions

		AppCtx    *user  = (AppCtx *)ctx;
		PetscReal dt     = user->dt;

		PetscReal Vs[2],dVs[2][2];															//Seems that if I create this one with size 3, next line gives problems
		IGAPointFormValue(pVs,UVs,&Vs[0]);													//This fills the values
		IGAPointFormGrad (pVs,UVs,&dVs[0][0]);												//This fills the values

		PetscReal fullVs[3]={0};
		PetscReal full_dVs[3][3]={0};
		//Restore this later
		//fullVs[0]=Vs[0]; fullVs[1]=Vs[1];
		//full_dVs[0][0]=dVs[0][0]; full_dVs[0][1]=dVs[0][1];
		//full_dVs[1][0]=dVs[1][0]; full_dVs[1][1]=dVs[1][1];
		//This is for testing
		fullVs[0]=0.0; fullVs[1]=-1.0;
		full_dVs[0][0]=0.0; full_dVs[0][1]=0.0;
		full_dVs[1][0]=0.0; full_dVs[1][1]=0.0;

		PetscReal S[8], dS[8][2];															//Create array to receive S and grad(S)
		IGAPointFormValue(pS,US,&S[0]);														//This fills the values of S (remember that S has 8 non zero components in 2D)
		IGAPointFormGrad (pS,US,&dS[0][0]);													//This fills the values

		PetscReal fullS[3][3][3]={0};
		fullS[0][0][0]=S[0]; fullS[0][0][1]=S[1];											//Expand S to full 3rd order form, only non-zero elements
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
										//KSdot[a][i][b][j]+=St[k][l][u]*v[k][l][u];		//This is for backwards Euler
										KSdot[a][i][b][j]+=St[k][l][u]*v[k][l][u];			//This is for trapezoidal rule (these are equal, consider deleting one?)
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
											//KSdot[a][i][b][j]+=dt*St[k][l][w]*fullVs[w]*dv[k][l][u][u];		//This is for backwards Euler
											KSdot[a][i][b][j]+=0.5*dt*(St[k][l][w]*fullVs[w])*dv[k][l][u][u];		//This is for trapezoidal rule

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
										//KSdot[a][i][b][j]+=St[k][l][u]*v[k][l][u];		//This is for backwards Euler
										KSdot[a][i][b][j]+=St[k][l][u]*v[k][l][u];			//This is for trapezoidal rule
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
											//KSdot[a][i][b][j]+=-dt*St[k][l][u]*(dv[k][l][w][u]*fullVs[w]+v[k][l][w]*full_dVs[w][u])
											//				   -dt*(dSt[k][l][w][u]*fullVs[w]+St[k][l][w]*full_dVs[w][u])*v[k][l][u]
											//				   +dt*dt*(dSt[k][l][w][u]*fullVs[w]+St[k][l][w]*full_dVs[w][u])*(dv[k][l][w][u]*fullVs[w]+v[k][l][w]*full_dVs[w][u]);			//This is for backwards Euler

											KSdot[a][i][b][j]+=-0.5*dt*St[k][l][u]*(dv[k][l][w][u]*fullVs[w]+v[k][l][w]*full_dVs[w][u])
											                   -0.5*dt*(dSt[k][l][w][u]*fullVs[w]+St[k][l][w]*full_dVs[w][u])*v[k][l][u]
											                   +0.25*dt*dt*(dSt[k][l][w][u]*fullVs[w]+St[k][l][w]*full_dVs[w][u])*(dv[k][l][w][u]*fullVs[w]+v[k][l][w]*full_dVs[w][u]);			//This is for trpezoidal rule
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
								//FSdot[a][i]+=fullS[k][l][u]*v[k][l][u];		//Term with (Pi x VPi) should be included here in the future
								//This is trapezoidal rule
								FSdot[a][i]+=fullS[k][l][u]*v[k][l][u];		//Term with (Pi x VPi) should be included here in the future (these are equal, consider deleting)
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
									FSdot[a][i]+=0.5*dt*(fulld_S[k][l][w][u]*fullVs[w]+fullS[k][l][w]*full_dVs[w][u])*v[k][l][u];		//Term with (Pi x VPi) should be included here in the future (these are equal, consider deleting)
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
								//FSdot[a][i]+=fullS[k][l][u]*v[k][l][u];		//Term with (Pi x VPi) should be included here too in the future
								//This is for trapezoidal rule
								FSdot[a][i]+=fullS[k][l][u]*v[k][l][u];		//Term with (Pi x VPi) should be included here too in the future (these terms are equal, consider deleting)
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
									//FSdot[a][i]+=-dt*fullS[k][l][u]*(dv[k][l][w][u]*fullVs[w]+v[k][l][w]*full_dVs[w][u]);		//Term with dt^2*(Pi x VPi) should be included here in the future
									//This term is for trapezoidal rule
									FSdot[a][i]+=0.5*dt*(fulld_S[k][l][w][u]*fullVs[w]+fullS[k][l][w]*full_dVs[w][u])*v[k][l][u]
									            -0.5*dt*fullS[k][l][u]*(dv[k][l][w][u]*fullVs[w]+v[k][l][w]*full_dVs[w][u])
									            -0.25*dt*dt*(fulld_S[k][l][w][u]*fullVs[w]+fullS[k][l][w]*full_dVs[w][u])*(dv[k][l][w][u]*fullVs[w]+v[k][l][w]*full_dVs[w][u]);		//Term with dt^2*(Pi x VPi) should be included here in the future
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
	//Sdot(pointSdot,pointVs,pointS,KpointSdot,FpointSdot,VsSdot,S0Sdot,NULL);CHKERRQ(ierr);
	PetscErrorCode SdotFuncF(IGAPoint p, IGAPoint pVs, IGAPoint pS, PetscReal *F, PetscReal *UVs, PetscReal *US, void *ctx)
	{
		PetscInt a,i,k,l,u,w,dof=p->dof,nen=p->nen;	

		const PetscReal *N0,(*N1)[2];
		IGAPointGetShapeFuns(p,0,(const PetscReal**)&N0);									//Value of the shape functions
		IGAPointGetShapeFuns(p,1,(const PetscReal**)&N1);									//Derivatives of the shape functions

		AppCtx    *user  = (AppCtx *)ctx;
		PetscReal dt     = user->dt;

		PetscReal Vs[2],dVs[2][2];															//Seems that if I create this one with size 3, next line gives problems
		IGAPointFormValue(pVs,UVs,&Vs[0]);													//This fills the values
		IGAPointFormGrad (pVs,UVs,&dVs[0][0]);												//This fills the values

		PetscReal fullVs[3]={0};
		PetscReal full_dVs[3][3]={0};
		//Restore this later
		//fullVs[0]=Vs[0]; fullVs[1]=Vs[1];
		//full_dVs[0][0]=dVs[0][0]; full_dVs[0][1]=dVs[0][1];
		//full_dVs[1][0]=dVs[1][0]; full_dVs[1][1]=dVs[1][1];
		
		//This is for testing
		fullVs[0]=0.0; fullVs[1]=-1.0;
		full_dVs[0][0]=0.0; full_dVs[0][1]=0.0;
		full_dVs[1][0]=0.0; full_dVs[1][1]=0.0;

		PetscReal S[8];																		//Create array to receive S
		IGAPointFormValue(pS,US,&S[0]);														//This fills the values of S (remember that S has 8 non zero components in 2D)

		PetscReal fullS[3][3][3]={0};
		fullS[0][0][0]=S[0]; fullS[0][0][1]=S[1];											//Expand S to full 3rd order form, only non-zero elements
		fullS[0][1][0]=S[2]; fullS[0][1][1]=S[3];
		fullS[1][0][0]=S[4]; fullS[1][0][1]=S[5];
		fullS[1][1][0]=S[6]; fullS[1][1][1]=S[7];

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
								FSdot[a][i]+=fullS[k][l][u]*v[k][l][u];		//Term with (Pi x VPi) should be included here in the future
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
								FSdot[a][i]+=fullS[k][l][u]*v[k][l][u];		//Term with (Pi x VPi) should be included here too in the future
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
									FSdot[a][i]+=-dt*fullS[k][l][u]*(dv[k][l][w][u]*fullVs[w]+v[k][l][w]*full_dVs[w][u]);		//Term with dt^2*(Pi x VPi) should be included here in the future
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
				 //SysUp(pointUp,pointZ0,pointchiUp,KpointUp,FpointUp,Z0Up,ChiUp,NULL);CHKERRQ(ierr);
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
	PetscErrorCode proj(IGAPoint p,IGAPoint pAl,IGAPoint pS, IGAPoint pVa,PetscReal *K,PetscReal *F,PetscReal *UAl,PetscReal *US,PetscReal *UVa,void *ctx)		//When other fields like Pi or S (etc.) come into this, declare a IGAPoint pPi or pS and a new PetscReal *UPi or *US for each
	{
		//This functions generates the L2 projection matrix and the RHS
		const PetscReal *N0;
		IGAPointGetShapeFuns(p,0,(const PetscReal**)&N0);									//Value of the shape functions
		PetscInt a,b,i,j,k,l,m,nen=p->nen, dof=p->dof;

		//Create and inflate alpha
		PetscReal alfa[2];																//Container for values of alpha
		IGAPointFormValue(pAl,UAl,&alfa[0]);											//Fills the values of the alphas
		PetscReal fullAlfa[3][3]={0};
		fullAlfa[0][2]=alfa[0]; fullAlfa[1][2]=alfa[1];										//Expand Alfa to full 2nd order form, only non-zero elements

		//Create and inflate S
		PetscReal S[8];																		//Create array to receive S
		IGAPointFormValue(pS,US,&S[0]);														//This fills the values of S (remember that S has 8 non zero components in 2D)
		PetscReal fullS[3][3][3]={0};
		fullS[0][0][0]=S[0]; fullS[0][0][1]=S[1];											//Expand S to full 3rd order form, only non-zero elements
		fullS[0][1][0]=S[2]; fullS[0][1][1]=S[3];
		fullS[1][0][0]=S[4]; fullS[1][0][1]=S[5];
		fullS[1][1][0]=S[6]; fullS[1][1][1]=S[7];

		//Creation and inflation of Va
		PetscReal Va[2];															//Seems that if I create this one with size 3, next line gives problems
		IGAPointFormValue(pVa,UVa,&Va[0]);													//This fills the values
		PetscReal fullVa[3]={0};	
		fullVa[0]=Va[0]; fullVa[1]=Va[1];

		PetscReal fullVs[3]={0};	
		fullVs[0]=0.0; fullVs[1]=-1.0;

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
									FCS[a][i]+=(e[k][l][m]*fullAlfa[j][l]*fullVa[m])*v[j][k];
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
								FCS[a][i]+=(fullS[j][k][l]*fullVs[l])*v[j][k];
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
	PetscPrintf(PETSC_COMM_WORLD,"Start of PruebaV4S \n");

	PetscInt commsize,rank;
	ierr = MPI_Comm_size(PETSC_COMM_WORLD,&commsize);CHKERRQ(ierr);
	ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank);CHKERRQ(ierr);
//

//App context creation and some data
	//Mesh parameters (to fix specific points in z0 system)
	PetscInt b=300;				//Parameter to choose size of cores, must always be odd, core will be of size 1 unit, rest of the body will be of size b-1 units in each direction
	PetscReal Lx=10.0;
	PetscReal Ly=10.0;
	PetscInt  nx=b;
	PetscInt  ny=b;
	
	//alpha has to be at MOST 0.99 for things to work, the smaller it is, the better accuracy
	PetscReal alpha=0.9*(1.0/8.0);
	PetscReal dt=alpha*0.5*fmin(Lx/nx,Ly/ny);

	AppCtx user;
	//user.Lx     = Lx;
	//user.Ly     = Ly;
	//user.nx     = nx;
	//user.ny     = ny;
	user.dt     = dt;
	//void* user;

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
		source = fopen("./PruebaV4S.c","r");
		dest   = fopen("../Results/PruebaV4S.c","w");

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
	PetscPrintf(PETSC_COMM_WORLD,"\nSystem for Initialization for S starting \n\n");
	T=time(NULL);
	tm=*localtime(&T);
	PetscPrintf(PETSC_COMM_WORLD,"Current time is %02d:%02d:%02d \n",tm.tm_hour,tm.tm_min,tm.tm_sec);
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

	//char nameS[]="/Input-S-2d-0.dat";			//Modify this so that the "-0.dat" part changes with the loop index 
	char nameS[512];
	sprintf(nameS,"%s%d%s","/Input-S-2d-",i,".dat");

	PetscPrintf(PETSC_COMM_WORLD,"%s\n",nameS);

	char pathS[512];
	sprintf(pathS,"%s%s",direct,nameS);
	ierr = IGAReadVec(igaS,s0,pathS); CHKERRQ(ierr);
//

//Creation of types and systems for the Helmholtz decomposition of S, curl part for S based input
	PetscPrintf(PETSC_COMM_WORLD,"\nSystem for curl part of Helmholtz of S starting \n\n");
	T=time(NULL);
	tm=*localtime(&T);
	PetscPrintf(PETSC_COMM_WORLD,"Current time is %02d:%02d:%02d \n",tm.tm_hour,tm.tm_min,tm.tm_sec);
	IGA igachiS;
	ierr = IGACreate(PETSC_COMM_WORLD,&igachiS);CHKERRQ(ierr);
	ierr = IGASetDim(igachiS,2);CHKERRQ(ierr);														//Spatial dimension of the problem
	ierr = IGASetDof(igachiS,8);CHKERRQ(ierr);														//Number of degrees of freedom, per node
	ierr = IGASetOrder(igachiS,1);CHKERRQ(ierr);														//Number of spatial derivatives to calculate
	ierr = IGASetFromOptions(igachiS);CHKERRQ(ierr);													//Note: The order (or degree) of the shape functions is given by the mesh!
	ierr = IGARead(igachiS,"./geometry.dat");CHKERRQ(ierr);
	
	for (dir=0; dir<2; dir++)
	{
		ierr = IGASetRuleType(igachiS,dir,IGA_RULE_LEGENDRE);CHKERRQ(ierr);
		ierr = IGASetRuleSize(igachiS,dir,6);CHKERRQ(ierr);
	}
	ierr = IGASetUp(igachiS);CHKERRQ(ierr);
	//PetscInt dir,side;
	for (dir=0; dir<2; dir++) 
	{
		for (side=0; side<2; side++) 
		{
			if(dir==0 && side==0)
			{
				ierr = IGASetBoundaryValue(igachiS,dir,side,0,0.0);CHKERRQ(ierr);    				// Dirichlet boundary conditions
				ierr = IGASetBoundaryValue(igachiS,dir,side,2,0.0);CHKERRQ(ierr);    				// Dirichlet boundary conditions
				ierr = IGASetBoundaryValue(igachiS,dir,side,4,0.0);CHKERRQ(ierr);    				// Dirichlet boundary conditions
				ierr = IGASetBoundaryValue(igachiS,dir,side,6,0.0);CHKERRQ(ierr);    				// Dirichlet boundary conditions
			}
			if(dir==0 && side==1)
			{
				ierr = IGASetBoundaryValue(igachiS,dir,side,0,0.0);CHKERRQ(ierr);    				// Dirichlet boundary conditions
				ierr = IGASetBoundaryValue(igachiS,dir,side,2,0.0);CHKERRQ(ierr);    				// Dirichlet boundary conditions
				ierr = IGASetBoundaryValue(igachiS,dir,side,4,0.0);CHKERRQ(ierr);    				// Dirichlet boundary conditions
				ierr = IGASetBoundaryValue(igachiS,dir,side,6,0.0);CHKERRQ(ierr);    				// Dirichlet boundary conditions
			}
			if(dir==1 && side==0)
			{
				ierr = IGASetBoundaryValue(igachiS,dir,side,1,0.0);CHKERRQ(ierr);    				// Dirichlet boundary conditions
				ierr = IGASetBoundaryValue(igachiS,dir,side,3,0.0);CHKERRQ(ierr);    				// Dirichlet boundary conditions
				ierr = IGASetBoundaryValue(igachiS,dir,side,5,0.0);CHKERRQ(ierr);    				// Dirichlet boundary conditions
				ierr = IGASetBoundaryValue(igachiS,dir,side,7,0.0);CHKERRQ(ierr);    				// Dirichlet boundary conditions
			}
			if(dir==1 && side==1)
			{
				ierr = IGASetBoundaryValue(igachiS,dir,side,1,0.0);CHKERRQ(ierr);    				// Dirichlet boundary conditions
				ierr = IGASetBoundaryValue(igachiS,dir,side,3,0.0);CHKERRQ(ierr);    				// Dirichlet boundary conditions
				ierr = IGASetBoundaryValue(igachiS,dir,side,5,0.0);CHKERRQ(ierr);    				// Dirichlet boundary conditions
				ierr = IGASetBoundaryValue(igachiS,dir,side,7,0.0);CHKERRQ(ierr);    				// Dirichlet boundary conditions
			}
			//ierr = IGASetBoundaryValue(iga,dir,side,dof,0.0);CHKERRQ(ierr);    				// Dirichlet boundary conditions
			//ierr = IGASetBoundaryForm(igachiS,dir,side,PETSC_TRUE);CHKERRQ(ierr);  				// Neumann boundary conditions
		}
	}
	
	Mat KchiS;
	Vec chiS0,FchiS;
	ierr = IGACreateMat(igachiS,&KchiS);CHKERRQ(ierr);
	ierr = IGACreateVec(igachiS,&chiS0);CHKERRQ(ierr);
	ierr = IGACreateVec(igachiS,&FchiS);CHKERRQ(ierr);

	IGAPoint        pointchiS, pointS;
	IGAElement      elemchiS, elemS;					//element
	PetscReal       *KlocchiS,*FlocchiS;			//AA y BB
	PetscReal       *KpointchiS,*FpointchiS;	//KKK y FFF
	const PetscReal *arrayS0chiS;				//arrayU
	Vec  			localS0chiS;				//localU
	PetscReal       *S0chiS;					//U

  	IGAFormSystem  wtfchiS;
 	void           *wtf2chiS;

 	KSP kspchiS;
	ierr = IGACreateKSP(igachiS,&kspchiS);CHKERRQ(ierr);

	//Get local vectors s0 and arrays
	ierr = IGAGetLocalVecArray(igaS,s0,&localS0chiS,&arrayS0chiS);CHKERRQ(ierr);

	//Element loop
	ierr = IGABeginElement(igachiS,&elemchiS);CHKERRQ(ierr);
	ierr = IGABeginElement(igaS,&elemS);CHKERRQ(ierr);

	while (IGANextElement(igachiS,elemchiS))
	{
		IGANextElement(igaS,elemS);
		ierr = IGAElementGetWorkMat(elemchiS,&KlocchiS);CHKERRQ(ierr);
		ierr = IGAElementGetWorkVec(elemchiS,&FlocchiS);CHKERRQ(ierr);
		ierr = IGAElementGetValues(elemS,arrayS0chiS,&S0chiS);CHKERRQ(ierr);

		//FormSystem loop
		while (IGAElementNextFormSystem(elemchiS,&wtfchiS,&wtf2chiS)) 
		{
			//Quadrature loop
			ierr = IGAElementBeginPoint(elemchiS,&pointchiS);CHKERRQ(ierr);
			ierr = IGAElementBeginPoint(elemS,&pointS);CHKERRQ(ierr);
			while (IGAElementNextPoint(elemchiS,pointchiS)) 
			{
				if(pointchiS->atboundary==0 && pointS->atboundary==0)
				{
					IGAElementNextPoint(elemS,pointS);
					ierr = IGAPointGetWorkMat(pointchiS,&KpointchiS);CHKERRQ(ierr);
					ierr = IGAPointGetWorkVec(pointchiS,&FpointchiS);CHKERRQ(ierr);
					ierr = curlChiS(pointchiS,pointS,KpointchiS,FpointchiS,S0chiS,NULL);CHKERRQ(ierr);
					ierr = IGAPointAddMat(pointchiS,KpointchiS,KlocchiS);CHKERRQ(ierr);
					ierr = IGAPointAddVec(pointchiS,FpointchiS,FlocchiS);CHKERRQ(ierr);
				}
			}
			while (pointS->index != -1)
			{
				IGAElementNextPoint(elemS,pointS);
			}
			ierr = IGAElementEndPoint(elemchiS,&pointchiS);CHKERRQ(ierr);
			ierr = IGAElementEndPoint(elemS,&pointS);CHKERRQ(ierr);
		}
		ierr = IGAElementFixSystem(elemchiS,KlocchiS,FlocchiS);CHKERRQ(ierr);
		ierr = IGAElementAssembleMat(elemchiS,KlocchiS,KchiS);CHKERRQ(ierr);
		ierr = IGAElementAssembleVec(elemchiS,FlocchiS,FchiS);CHKERRQ(ierr);

	}
	IGANextElement(igaS,elemS);
	ierr = IGAEndElement(igachiS,&elemchiS);CHKERRQ(ierr);
	ierr = IGAEndElement(igaS,&elemS);CHKERRQ(ierr);

	// Restore local vectors s0 and arrays
	ierr = IGARestoreLocalVecArray(igaS,s0,&localS0chiS,&arrayS0chiS);CHKERRQ(ierr);

	//Form system matrix and vector
	ierr = MatAssemblyBegin(KchiS,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyEnd  (KchiS,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

	ierr = VecAssemblyBegin(FchiS);CHKERRQ(ierr);
	ierr = VecAssemblyEnd  (FchiS);CHKERRQ(ierr);

	ierr = KSPSetOperators(kspchiS,KchiS,KchiS);CHKERRQ(ierr);
	PC pcchiS;
	ierr = KSPGetPC(kspchiS,&pcchiS); CHKERRQ(ierr);
	ierr = PCSetType(pcchiS,PCLU); CHKERRQ(ierr);
	ierr = PCFactorSetMatSolverType(pcchiS,MATSOLVERMUMPS); CHKERRQ(ierr);
	//ierr = KSPSetFromOptions(kspchiS);CHKERRQ(ierr);
	//ierr = KSPSetTolerances(kspchiS,1.0e-12,5.0e-24,PETSC_DEFAULT,PETSC_DEFAULT);CHKERRQ(ierr);
	ierr = KSPSolve(kspchiS,FchiS,chiS0);CHKERRQ(ierr);

	//ierr = KSPDestroy(&kspchiS);CHKERRQ(ierr);
	//ierr = MatDestroy(&KchiS);CHKERRQ(ierr);
	//ierr = VecDestroy(&FchiS);CHKERRQ(ierr);
	
	char namechiS[512];//="/ChiS-2d-0.dat";
	sprintf(namechiS,"%s%d%s","/ChiS-2d-",i,".dat");
	char pathchiS[512];
	sprintf(pathchiS,"%s%s",direct,namechiS);
	ierr = IGAWriteVec(igachiS,chiS0,pathchiS);CHKERRQ(ierr);
//

/*
//Creation of types and systems for the Helmholtz decomposition of S, grad part
	//System for chiS
	PetscPrintf(PETSC_COMM_WORLD,"\n System for grad part of Helmholtz of S starting \n\n");
	IGA igaZS;
	ierr = IGACreate(PETSC_COMM_WORLD,&igaZS);CHKERRQ(ierr);
	ierr = IGASetDim(igaZS,2);CHKERRQ(ierr);														//Spatial dimension of the problem
	ierr = IGASetDof(igaZS,4);CHKERRQ(ierr);														//Number of degrees of freedom, per node
	ierr = IGASetOrder(igaZS,2);CHKERRQ(ierr);														//Number of spatial derivatives to calculate
	ierr = IGASetFromOptions(igaZS);CHKERRQ(ierr);													//Note: The order (or degree) of the shape functions is given by the mesh!
	ierr = IGARead(igaZS,"./geometry2.dat");CHKERRQ(ierr);
	ierr = IGASetUp(igaZS);CHKERRQ(ierr);
	
	//PetscInt dir,side;
	for (dir=0; dir<2; dir++) 
	{
		for (side=0; side<2; side++) 
		{
			//ierr = IGASetBoundaryValue(iga,dir,side,dof,0.0);CHKERRQ(ierr);    				// Dirichlet boundary conditions
			ierr = IGASetBoundaryForm(igaZS,dir,side,PETSC_TRUE);CHKERRQ(ierr);  				// Neumann boundary conditions
		}
	}
	
	Mat KZS;
	Vec ZS0,FZS;
	ierr = IGACreateMat(igaZS,&KZS);CHKERRQ(ierr);
	ierr = IGACreateVec(igaZS,&ZS0);CHKERRQ(ierr);
	ierr = IGACreateVec(igaZS,&FZS);CHKERRQ(ierr);

	IGAPoint        pointZS;
	IGAElement      elemZS;					//element
	PetscReal       *KlocZS,*FlocZS;			//AA y BB
	PetscReal       *KpointZS,*FpointZS;	//KKK y FFF
	const PetscReal *arrayS0ZS;				//arrayU
	Vec  			localS0ZS;				//localU
	PetscReal       *S0ZS;					//U

  	IGAFormSystem  wtfZS;
 	void           *wtf2ZS;

 	KSP kspZS;
	ierr = IGACreateKSP(igaZS,&kspZS);CHKERRQ(ierr);

	//Get local vectors s0 and arrays
	ierr = IGAGetLocalVecArray(igaS,s0,&localS0ZS,&arrayS0ZS);CHKERRQ(ierr);

	//Element loop
	ierr = IGABeginElement(igaZS,&elemZS);CHKERRQ(ierr);
	ierr = IGABeginElement(igaS,&elemS);CHKERRQ(ierr);

	while (IGANextElement(igaZS,elemZS))
	{
		IGANextElement(igaS,elemS);
		ierr = IGAElementGetWorkMat(elemZS,&KlocZS);CHKERRQ(ierr);
		ierr = IGAElementGetWorkVec(elemZS,&FlocZS);CHKERRQ(ierr);
		ierr = IGAElementGetValues(elemS,arrayS0ZS,&S0ZS);CHKERRQ(ierr);

		//FormSystem loop
		while (IGAElementNextFormSystem(elemZS,&wtfZS,&wtf2ZS)) 
		{
			//Quadrature loop
			ierr = IGAElementBeginPoint(elemZS,&pointZS);CHKERRQ(ierr);
			ierr = IGAElementBeginPoint(elemS,&pointS);CHKERRQ(ierr);
			if(pointZS->atboundary==0 && pointS->atboundary==0)
			{
				while (IGAElementNextPoint(elemZS,pointZS)) 
				{
					IGAElementNextPoint(elemS,pointS);
					ierr = IGAPointGetWorkMat(pointZS,&KpointZS);CHKERRQ(ierr);
					ierr = IGAPointGetWorkVec(pointZS,&FpointZS);CHKERRQ(ierr);
					ierr = gradZS(pointZS,pointS,KpointZS,FpointZS,S0ZS,NULL);CHKERRQ(ierr);
					ierr = IGAPointAddMat(pointZS,KpointZS,KlocZS);CHKERRQ(ierr);
					ierr = IGAPointAddVec(pointZS,FpointZS,FlocZS);CHKERRQ(ierr);
				}
				IGAElementNextPoint(elemS,pointS);

				ierr = IGAElementEndPoint(elemZS,&pointZS);CHKERRQ(ierr);
				ierr = IGAElementEndPoint(elemS,&pointS);CHKERRQ(ierr);
			}
		}
		//ierr = IGAElementFixSystem(elemZS,KlocZS,FlocZS);CHKERRQ(ierr);
		ierr = IGAElementAssembleMat(elemZS,KlocZS,KZS);CHKERRQ(ierr);
		ierr = IGAElementAssembleVec(elemZS,FlocZS,FZS);CHKERRQ(ierr);

	}
	IGANextElement(igaS,elemS);
	ierr = IGAEndElement(igaZS,&elemZS);CHKERRQ(ierr);
	ierr = IGAEndElement(igaS,&elemS);CHKERRQ(ierr);
	
	// Restore local vectors s0 and arrays
	ierr = IGARestoreLocalVecArray(igaS,s0,&localS0ZS,&arrayS0ZS);CHKERRQ(ierr);
	
	//Impose Dirichlet condition on a single point
	//First we replace rows and columns with associated rows and columns from identity matrix
	PetscInt mz,nz;
	ierr = MatGetSize(KZS,&nz,&mz);CHKERRQ(ierr);
	ierr = MatSetOption(KZS, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);CHKERRQ(ierr);
	ierr = MatAssemblyBegin(KZS,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyEnd  (KZS,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);	

	PetscInt rowsZ, *colsZ;
	PetscReal valZ, *valsZ;

	ierr = PetscMalloc1(mz,&colsZ);CHKERRQ(ierr);
	ierr = PetscMalloc1(mz,&valsZ);CHKERRQ(ierr);

	for(int i=0;i<mz;i++)
	{
		colsZ[i]=i;
		valsZ[i]=0.0;
	}

	rowsZ=0;
	valsZ[rowsZ]=1.0e6;
	ierr = MatSetValues(KZS,1,&rowsZ,mz,colsZ,valsZ,INSERT_VALUES);CHKERRQ(ierr);
	ierr = MatSetValues(KZS,nz,colsZ,1,&rowsZ,valsZ,INSERT_VALUES);CHKERRQ(ierr);
	valsZ[rowsZ]=0.0;
	rowsZ=1;
	valsZ[rowsZ]=1.0e6;
	ierr = MatSetValues(KZS,1,&rowsZ,mz,colsZ,valsZ,INSERT_VALUES);CHKERRQ(ierr);
	ierr = MatSetValues(KZS,nz,colsZ,1,&rowsZ,valsZ,INSERT_VALUES);CHKERRQ(ierr);
	valsZ[rowsZ]=0.0;
	rowsZ=2;
	valsZ[rowsZ]=1.0e6;
	ierr = MatSetValues(KZS,1,&rowsZ,mz,colsZ,valsZ,INSERT_VALUES);CHKERRQ(ierr);
	ierr = MatSetValues(KZS,nz,colsZ,1,&rowsZ,valsZ,INSERT_VALUES);CHKERRQ(ierr);
	valsZ[rowsZ]=0.0;
	rowsZ=3;
	valsZ[rowsZ]=1.0e6;
	ierr = MatSetValues(KZS,1,&rowsZ,mz,colsZ,valsZ,INSERT_VALUES);CHKERRQ(ierr);
	ierr = MatSetValues(KZS,nz,colsZ,1,&rowsZ,valsZ,INSERT_VALUES);CHKERRQ(ierr);
	valsZ[rowsZ]=0.0;

	ierr = MatAssemblyBegin(KZS,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyEnd  (KZS,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

	//Here we set values to the vector directly
	ierr = VecAssemblyBegin(FZS);CHKERRQ(ierr);
	ierr = VecAssemblyEnd  (FZS);CHKERRQ(ierr);

	rowsZ=0;
	valZ=0.0;
	ierr = VecSetValue(FZS,rowsZ,valZ,INSERT_VALUES);CHKERRQ(ierr);
	rowsZ=1;
	valZ=0.0;
	ierr = VecSetValue(FZS,rowsZ,valZ,INSERT_VALUES);CHKERRQ(ierr);
	rowsZ=2;
	valZ=0.0;
	ierr = VecSetValue(FZS,rowsZ,valZ,INSERT_VALUES);CHKERRQ(ierr);
	rowsZ=3;
	valZ=0.0;
	ierr = VecSetValue(FZS,rowsZ,valZ,INSERT_VALUES);CHKERRQ(ierr);

	ierr = VecAssemblyBegin(FZS);CHKERRQ(ierr);
	ierr = VecAssemblyEnd  (FZS);CHKERRQ(ierr);
	//

	ierr = KSPSetOperators(kspZS,KZS,KZS);CHKERRQ(ierr);
	//PC pcZs;
	//ierr = KSPGetPC(kspZS,&pcZs); CHKERRQ(ierr);
	//ierr = PCSetType(pcZs,PCLU); CHKERRQ(ierr);
	//ierr = PCFactorSetMatSolverType(pcZs,MATSOLVERMUMPS); CHKERRQ(ierr);
	ierr = KSPSetFromOptions(kspZS);CHKERRQ(ierr);
	ierr = KSPSetType(kspZS,KSPCG);																				//UNTESTED LINE, BE CAREFUL!!!
	ierr = KSPSetTolerances(kspZS,1.0e-12,1.0e-24,PETSC_DEFAULT,PETSC_DEFAULT);CHKERRQ(ierr);
	ierr = KSPSolve(kspZS,FZS,ZS0);CHKERRQ(ierr);

	ierr = KSPDestroy(&kspZS);CHKERRQ(ierr);
	ierr = MatDestroy(&KZS);CHKERRQ(ierr);
	ierr = VecDestroy(&FZS);CHKERRQ(ierr);
	
	char nameZS[]="/ZS-2d-0.dat";
	char pathZS[512];
	sprintf(pathZS,"%s%s",direct,nameZS);
	ierr = IGAWriteVec(igaZS,ZS0,pathZS);CHKERRQ(ierr);
//
*/

//Creation of types and systems for the L2 projection of Alfa0-Sp:X
	//System for Alfa
	PetscPrintf(PETSC_COMM_WORLD,"\nSystem for L2 projection for Alfa-Sp:X starting \n\n");
	T=time(NULL);
	tm=*localtime(&T);
	PetscPrintf(PETSC_COMM_WORLD,"Current time is %02d:%02d:%02d \n",tm.tm_hour,tm.tm_min,tm.tm_sec);
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
	//PetscInt dir,side;
	for (dir=0; dir<2; dir++) 
	{
		for (side=0; side<2; side++) 
		{
			//ierr = IGASetBoundaryValue(iga,dir,side,dof,0.0);CHKERRQ(ierr);    				// Dirichlet boundary conditions
			//ierr = IGASetBoundaryForm(igaAl,dir,side,PETSC_TRUE);CHKERRQ(ierr);  				// Neumann boundary conditions
		}
	}
	
	Mat KAlp;
	Vec alp0,alInput,FAlp;
	ierr = IGACreateMat(igaAl,&KAlp);CHKERRQ(ierr);
	ierr = IGACreateVec(igaAl,&alp0);CHKERRQ(ierr);
	ierr = IGACreateVec(igaAl,&alInput);CHKERRQ(ierr);
	ierr = IGACreateVec(igaAl,&FAlp);CHKERRQ(ierr);

	IGAPoint        pointAlp, pointSp;
	IGAElement      elemAlp, elemSp;				//element
	PetscReal       *KlocAlp,*FlocAlp;				//AA y BB
	PetscReal       *KpointAlp,*FpointAlp;			//KKK y FFF
	const PetscReal *arrayS0Alp;					//arrayU
	Vec  			localS0Alp;						//localU
	PetscReal       *S0Alp;							//U

  	IGAFormSystem  wtfAlp;
 	void           *wtf2Alp;

 	KSP kspAlp;
	ierr = IGACreateKSP(igaAl,&kspAlp);CHKERRQ(ierr);

	//Get local vectors s0 and arrays
	ierr = IGAGetLocalVecArray(igachiS,chiS0,&localS0Alp,&arrayS0Alp);CHKERRQ(ierr);

	//Element loop
	ierr = IGABeginElement(igaAl,&elemAlp);CHKERRQ(ierr);
	ierr = IGABeginElement(igachiS,&elemSp);CHKERRQ(ierr);

	while (IGANextElement(igaAl,elemAlp))
	{
		IGANextElement(igachiS,elemSp);
		ierr = IGAElementGetWorkMat(elemAlp,&KlocAlp);CHKERRQ(ierr);
		ierr = IGAElementGetWorkVec(elemAlp,&FlocAlp);CHKERRQ(ierr);
		ierr = IGAElementGetValues(elemSp,arrayS0Alp,&S0Alp);CHKERRQ(ierr);

		//FormSystem loop
		while (IGAElementNextFormSystem(elemAlp,&wtfAlp,&wtf2Alp)) 
		{
			//Quadrature loop
			ierr = IGAElementBeginPoint(elemAlp,&pointAlp);CHKERRQ(ierr);
			ierr = IGAElementBeginPoint(elemSp,&pointSp);CHKERRQ(ierr);
			
			while (IGAElementNextPoint(elemAlp,pointAlp)) 
			{
				if(pointAlp->atboundary==0 && pointSp->atboundary==0)
				{
					IGAElementNextPoint(elemSp,pointSp);

					ierr = IGAPointGetWorkMat(pointAlp,&KpointAlp);CHKERRQ(ierr);
					ierr = IGAPointGetWorkVec(pointAlp,&FpointAlp);CHKERRQ(ierr);
					ierr = L2ProjectionAlSp(pointAlp,pointSp,KpointAlp,FpointAlp,S0Alp,&user);CHKERRQ(ierr);
					ierr = IGAPointAddMat(pointAlp,KpointAlp,KlocAlp);CHKERRQ(ierr);
					ierr = IGAPointAddVec(pointAlp,FpointAlp,FlocAlp);CHKERRQ(ierr);
				}
			}
			IGAElementNextPoint(elemSp,pointSp);

			ierr = IGAElementEndPoint(elemAlp,&pointAlp);CHKERRQ(ierr);
			ierr = IGAElementEndPoint(elemSp,&pointSp);CHKERRQ(ierr);
		}
		//ierr = IGAElementFixSystem(elemAlp,KlocAlp,FlocAlp);CHKERRQ(ierr);					//This sets Dirichlet condition ¿?
		ierr = IGAElementAssembleMat(elemAlp,KlocAlp,KAlp);CHKERRQ(ierr);
		ierr = IGAElementAssembleVec(elemAlp,FlocAlp,FAlp);CHKERRQ(ierr);

	}
	IGANextElement(igachiS,elemSp);
	ierr = IGAEndElement(igaAl,&elemAlp);CHKERRQ(ierr);
	ierr = IGAEndElement(igachiS,&elemSp);CHKERRQ(ierr);

	// Restore local vectors s0 and arrays
	ierr = IGARestoreLocalVecArray(igachiS,chiS0,&localS0Alp,&arrayS0Alp);CHKERRQ(ierr);

	//Form system matrix and vector
	ierr = MatAssemblyBegin(KAlp,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyEnd  (KAlp,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	
	ierr = VecAssemblyBegin(FAlp);CHKERRQ(ierr);
	ierr = VecAssemblyEnd  (FAlp);CHKERRQ(ierr);

	ierr = KSPSetOperators(kspAlp,KAlp,KAlp);CHKERRQ(ierr);
	ierr = KSPSetFromOptions(kspAlp);CHKERRQ(ierr);
	ierr = KSPSetTolerances(kspAlp,1.0e-25,1.0e-40,PETSC_DEFAULT,PETSC_DEFAULT);CHKERRQ(ierr);
	ierr = KSPSolve(kspAlp,FAlp,alp0);CHKERRQ(ierr);

	char nameAlInput[]="/Input-Al-2d-0.dat";
	char pathAlInput[512];
	sprintf(pathAlInput,"%s%s",direct,nameAlInput);
	ierr = IGAReadVec(igaAl,alInput,pathAlInput); CHKERRQ(ierr);

	ierr = VecAXPY(alp0,1.0,alInput); CHKERRQ(ierr);

	//ierr = KSPDestroy(&kspAlp);CHKERRQ(ierr);
	//ierr = MatDestroy(&KAlp);CHKERRQ(ierr);
	//ierr = VecDestroy(&FAlp);CHKERRQ(ierr);
	
	char nameAlp[512];//="/Alp-2d-0.dat";
	sprintf(nameAlp,"%s%d%s","/Alp-2d-",i,".dat");
	char pathAlp[512];
	sprintf(pathAlp,"%s%s",direct,nameAlp);
	ierr = IGAWriteVec(igaAl,alp0,pathAlp);CHKERRQ(ierr);

	//ierr = IGADestroy(&igachiS);CHKERRQ(ierr);		//This system is one of the more memory intensive. It's not required from here down, so better to destroy it.
//

//Creation of types and systems for the Helmholtz decomposition of Up (or Ue), curl part
	PetscPrintf(PETSC_COMM_WORLD,"\nSystem for curl part of Helmholtz of Up starting \n\n");
	T=time(NULL);
	tm=*localtime(&T);
	PetscPrintf(PETSC_COMM_WORLD,"Current time is %02d:%02d:%02d \n",tm.tm_hour,tm.tm_min,tm.tm_sec);
	IGA igachiUp;
	ierr = IGACreate(PETSC_COMM_WORLD,&igachiUp);CHKERRQ(ierr);
	ierr = IGASetDim(igachiUp,2);CHKERRQ(ierr);														//Spatial dimension of the problem
	ierr = IGASetDof(igachiUp,4);CHKERRQ(ierr);														//Number of degrees of freedom, per node
	ierr = IGASetOrder(igachiUp,1);CHKERRQ(ierr);													//Number of spatial derivatives to calculate
	ierr = IGASetFromOptions(igachiUp);CHKERRQ(ierr);												//Note: The order (or degree) of the shape functions is given by the mesh!
	ierr = IGARead(igachiUp,"./geometry.dat");CHKERRQ(ierr);
	
	for (dir=0; dir<2; dir++)
	{
		ierr = IGASetRuleType(igachiUp,dir,IGA_RULE_LEGENDRE);CHKERRQ(ierr);
		ierr = IGASetRuleSize(igachiUp,dir,6);CHKERRQ(ierr);
	}
	ierr = IGASetUp(igachiUp);CHKERRQ(ierr);

	//ierr = IGASetBoundaryValue(iga,dir,side,dof,val);CHKERRQ(ierr);
	//dir=0 and side=0
	ierr = IGASetBoundaryValue(igachiUp,0,0,0,0.0);CHKERRQ(ierr);					// Dirichlet boundary conditions to get chi \dot n =0
	ierr = IGASetBoundaryValue(igachiUp,0,0,2,0.0);CHKERRQ(ierr);
	//dir=0 and side=1
	ierr = IGASetBoundaryValue(igachiUp,0,1,0,0.0);CHKERRQ(ierr);					// Dirichlet boundary conditions to get chi \dot n =0
	ierr = IGASetBoundaryValue(igachiUp,0,1,2,0.0);CHKERRQ(ierr);
	//dir=1 and side=0
	ierr = IGASetBoundaryValue(igachiUp,1,0,1,0.0);CHKERRQ(ierr);					// Dirichlet boundary conditions to get chi \dot n =0
	ierr = IGASetBoundaryValue(igachiUp,1,0,3,0.0);CHKERRQ(ierr);
	//dir=1 and side=1
	ierr = IGASetBoundaryValue(igachiUp,1,1,1,0.0);CHKERRQ(ierr);					// Dirichlet boundary conditions to get chi \dot n =0
	ierr = IGASetBoundaryValue(igachiUp,1,1,3,0.0);CHKERRQ(ierr);
	
	Mat KchiUp;
	Vec chiUp0,FchiUp;
	ierr = IGACreateMat(igachiUp,&KchiUp);CHKERRQ(ierr);
	ierr = IGACreateVec(igachiUp,&chiUp0);CHKERRQ(ierr);
	ierr = IGACreateVec(igachiUp,&FchiUp);CHKERRQ(ierr);

	IGAPoint        pointchiUp;
	IGAElement      elemchiUp;							//element
	PetscReal       *KlocchiUp,*FlocchiUp;				//AA y BB
	PetscReal       *KpointchiUp,*FpointchiUp;			//KKK y FFF
	const PetscReal *arrayAl0chiUp;		//arrayU
	Vec  			localAl0chiUp;		//localU
	PetscReal       *Al0chiUp;							//U

  	IGAFormSystem  wtfchiUp;
 	void           *wtf2chiUp;

 	KSP kspchiUp;
	ierr = IGACreateKSP(igachiUp,&kspchiUp);CHKERRQ(ierr);

	//Get local vectors s0 and arrays
	ierr = IGAGetLocalVecArray(igaAl,alp0,&localAl0chiUp,&arrayAl0chiUp);CHKERRQ(ierr);

	//Element loop
	ierr = IGABeginElement(igachiUp,&elemchiUp);CHKERRQ(ierr);
	ierr = IGABeginElement(igaAl,&elemAlp);CHKERRQ(ierr);

	while (IGANextElement(igachiUp,elemchiUp))
	{
		IGANextElement(igaAl,elemAlp);
		ierr = IGAElementGetWorkMat(elemchiUp,&KlocchiUp);CHKERRQ(ierr);
		ierr = IGAElementGetWorkVec(elemchiUp,&FlocchiUp);CHKERRQ(ierr);
		ierr = IGAElementGetValues(elemAlp,arrayAl0chiUp,&Al0chiUp);CHKERRQ(ierr);

		//FormSystem loop
		while (IGAElementNextFormSystem(elemchiUp,&wtfchiUp,&wtf2chiUp)) 
		{
			//Quadrature loop
			ierr = IGAElementBeginPoint(elemchiUp,&pointchiUp);CHKERRQ(ierr);
			ierr = IGAElementBeginPoint(elemAlp,&pointAlp);CHKERRQ(ierr);

			while (IGAElementNextPoint(elemchiUp,pointchiUp)) 
			{
				if(pointchiUp->atboundary==0 && pointAlp->atboundary==0)
				{
					IGAElementNextPoint(elemAlp,pointAlp);
					ierr = IGAPointGetWorkMat(pointchiUp,&KpointchiUp);CHKERRQ(ierr);
					ierr = IGAPointGetWorkVec(pointchiUp,&FpointchiUp);CHKERRQ(ierr);
					ierr = curlChiU(pointchiUp,pointAlp,KpointchiUp,FpointchiUp,Al0chiUp,NULL);CHKERRQ(ierr);
					ierr = IGAPointAddMat(pointchiUp,KpointchiUp,KlocchiUp);CHKERRQ(ierr);
					ierr = IGAPointAddVec(pointchiUp,FpointchiUp,FlocchiUp);CHKERRQ(ierr);
				}
			}
			IGAElementNextPoint(elemAlp,pointAlp);

			ierr = IGAElementEndPoint(elemchiUp,&pointchiUp);CHKERRQ(ierr);
			ierr = IGAElementEndPoint(elemAlp,&pointAlp);CHKERRQ(ierr);
		}
		ierr = IGAElementFixSystem(elemchiUp,KlocchiUp,FlocchiUp);CHKERRQ(ierr);
		ierr = IGAElementAssembleMat(elemchiUp,KlocchiUp,KchiUp);CHKERRQ(ierr);
		ierr = IGAElementAssembleVec(elemchiUp,FlocchiUp,FchiUp);CHKERRQ(ierr);

	}
	IGANextElement(igaAl,elemAlp);
	ierr = IGAEndElement(igachiUp,&elemchiUp);CHKERRQ(ierr);
	ierr = IGAEndElement(igaAl,&elemAlp);CHKERRQ(ierr);

	// Restore local vectors s0 and arrays
	ierr = IGARestoreLocalVecArray(igaAl,alp0,&localAl0chiUp,&arrayAl0chiUp);CHKERRQ(ierr);

	//Form system matrix and vector
	ierr = MatAssemblyBegin(KchiUp,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyEnd  (KchiUp,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	
	ierr = VecAssemblyBegin(FchiUp);CHKERRQ(ierr);
	ierr = VecAssemblyEnd  (FchiUp);CHKERRQ(ierr);

	ierr = KSPSetOperators(kspchiUp,KchiUp,KchiUp);CHKERRQ(ierr);
	PC pcChiUp;
	ierr = KSPGetPC(kspchiUp,&pcChiUp); CHKERRQ(ierr);
	ierr = PCSetType(pcChiUp,PCLU); CHKERRQ(ierr);
	ierr = PCFactorSetMatSolverType(pcChiUp,MATSOLVERMUMPS); CHKERRQ(ierr);
	//ierr = KSPSetFromOptions(kspchiUp);CHKERRQ(ierr);
	//ierr = KSPSetTolerances(kspchiUp,1.0e-8,1.0e-20,PETSC_DEFAULT,PETSC_DEFAULT);CHKERRQ(ierr);
	ierr = KSPSolve(kspchiUp,FchiUp,chiUp0);CHKERRQ(ierr);

	//ierr = KSPDestroy(&kspchiUp);CHKERRQ(ierr);
	//ierr = MatDestroy(&KchiUp);CHKERRQ(ierr);
	//ierr = VecDestroy(&FchiUp);CHKERRQ(ierr);

	//Si hay problemas, borrar esto
	//ierr = VecDestroy(&alp0); CHKERRQ(ierr);
	//hasta aquí

	char namechiUp[512];//="/ChiUp-2d-0.dat";
	sprintf(namechiUp,"%s%d%s","/ChiUp-2d-",i,".dat");
	char pathchiUp[512];
	sprintf(pathchiUp,"%s%s",direct,namechiUp);
	ierr = IGAWriteVec(igachiUp,chiUp0,pathchiUp);CHKERRQ(ierr);
//

//System for initial state of z0
	PetscPrintf(PETSC_COMM_WORLD,"\nSystem for Z0 starting \n\n");
	T=time(NULL);
	tm=*localtime(&T);
	PetscPrintf(PETSC_COMM_WORLD,"Current time is %02d:%02d:%02d \n",tm.tm_hour,tm.tm_min,tm.tm_sec);
	IGA igaZ0;
	ierr = IGACreate(PETSC_COMM_WORLD,&igaZ0);CHKERRQ(ierr);
	ierr = IGASetDim(igaZ0,2);CHKERRQ(ierr);													//Spatial dimension of the problem
	ierr = IGASetDof(igaZ0,2);CHKERRQ(ierr);													//Number of degrees of freedom, per node
	ierr = IGASetOrder(igaZ0,1);CHKERRQ(ierr);													//Number of spatial derivatives to calculate
	ierr = IGASetFromOptions(igaZ0);CHKERRQ(ierr);												//Note: The order (or degree) of the shape functions is given by the mesh!
	ierr = IGARead(igaZ0,"./geometry.dat");CHKERRQ(ierr);
	
	for (dir=0; dir<2; dir++)
	{
		ierr = IGASetRuleType(igaZ0,dir,IGA_RULE_LEGENDRE);CHKERRQ(ierr);
		ierr = IGASetRuleSize(igaZ0,dir,6);CHKERRQ(ierr);
	}
	ierr = IGASetUp(igaZ0);CHKERRQ(ierr);

	PetscInt fijaPunto=0;																		//Fix a single point (1) or a side (chosen in blocks below)

	ierr = IGASetUp(igaZ0);CHKERRQ(ierr);
	ierr = IGASetUp(igachiUp);CHKERRQ(ierr);


	for (dir=0; dir<2; dir++) 
	{
		for (side=0; side<2; side++) 
		{
			//ierr = IGASetBoundaryValue(iga,dir,side,dof,0.0);CHKERRQ(ierr);					// Dirichlet boundary conditions
			ierr = IGASetBoundaryForm(igaZ0,dir,side,PETSC_TRUE);CHKERRQ(ierr);  				// Neumann boundary conditions
		}
	}

	if (fijaPunto==0)
	{
		//If we are not fixing a single point, set Dirichlet conditions here
		//ierr = IGASetBoundaryValue(iga,dir,side,dof,value);CHKERRQ(ierr);					// Dirichlet boundary conditions
		//ierr = IGASetBoundaryValue(igaZ0,0,0,0,0.0);CHKERRQ(ierr);	//Left side, 1st dof = 0
		//ierr = IGASetBoundaryValue(igaZ0,0,0,1,0.0);CHKERRQ(ierr);	//Left side, 2nd dof = 0

		//ierr = IGASetBoundaryValue(igaZ0,0,1,0,0.0);CHKERRQ(ierr);	//Right side, 1st dof=0
		//ierr = IGASetBoundaryValue(igaZ0,0,1,1,0.0);CHKERRQ(ierr);	//Right side, 2nd dof=0

		ierr = IGASetBoundaryValue(igaZ0,1,0,0,0.0);CHKERRQ(ierr);	//Bottom side, 1st dof=0
		ierr = IGASetBoundaryValue(igaZ0,1,0,1,0.0);CHKERRQ(ierr);	//Bottom side, 2nd dof=0
		
		//ierr = IGASetBoundaryValue(igaZ0,1,1,0,0.0);CHKERRQ(ierr);	//Top side, 1st dof=0
		//ierr = IGASetBoundaryValue(igaZ0,1,1,1,0.0);CHKERRQ(ierr);	//Top side, 2nd dof=0
	}

	Mat KZ0;
	Vec Z0,FZ0;

	ierr = IGACreateMat(igaZ0,&KZ0);CHKERRQ(ierr);
	ierr = IGACreateVec(igaZ0,&Z0);CHKERRQ(ierr);
	ierr = IGACreateVec(igaZ0,&FZ0);CHKERRQ(ierr);

	IGAPoint		pointZ0;								//point
	IGAElement		elemZ0;									//element
	PetscReal		*KlocZ0,*FlocZ0;						//AA y BB
	PetscReal		*KpointZ0,*FpointZ0;					//KKK y FFF
	const PetscReal	*arrayChi0Z0;							//arrayU
	Vec				localChi0Z0;							//localU
	PetscReal		*Chi0Z0;								//U

  	IGAFormSystem	wtfZ0;
 	void			*wtf2Z0;

 	KSP kspZ0;
	ierr = IGACreateKSP(igaZ0,&kspZ0);CHKERRQ(ierr);

	// Get local vectors Chi0, S and arrays
	ierr = IGAGetLocalVecArray(igachiUp,chiUp0,&localChi0Z0,&arrayChi0Z0);CHKERRQ(ierr);

	// Element loop
	ierr = IGABeginElement(igaZ0,&elemZ0);CHKERRQ(ierr);
	ierr = IGABeginElement(igachiUp,&elemchiUp);CHKERRQ(ierr);

	while (IGANextElement(igaZ0,elemZ0)) 
	{
		IGANextElement(igachiUp,elemchiUp);

		ierr = IGAElementGetWorkMat(elemZ0,&KlocZ0);CHKERRQ(ierr);
		ierr = IGAElementGetWorkVec(elemZ0,&FlocZ0);CHKERRQ(ierr);
		ierr = IGAElementGetValues(elemchiUp,arrayChi0Z0,&Chi0Z0);CHKERRQ(ierr);

		// FormSystem loop
		while (IGAElementNextFormSystem(elemZ0,&wtfZ0,&wtf2Z0)) 
		{
		// Quadrature loop
			ierr = IGAElementBeginPoint(elemZ0,&pointZ0);CHKERRQ(ierr);
			ierr = IGAElementBeginPoint(elemchiUp,&pointchiUp);CHKERRQ(ierr);

			while (IGAElementNextPoint(elemZ0,pointZ0))
			{
				if(pointZ0->atboundary==1)
				{
					ierr = IGAPointGetWorkMat(pointZ0,&KpointZ0);CHKERRQ(ierr);
					ierr = IGAPointGetWorkVec(pointZ0,&FpointZ0);CHKERRQ(ierr);
					//	   Z0sys(IGAPoint p, IGAPoint pChi, IGAPoint pZ, PetscReal *K, PetscReal *F, PetscReal *UChi, PetscReal *S, void *ctx)
					ierr = Z0sys(pointZ0,pointchiUp,KpointZ0,FpointZ0,Chi0Z0,NULL);CHKERRQ(ierr);
					ierr = IGAPointAddMat(pointZ0,KpointZ0,KlocZ0);CHKERRQ(ierr);
					ierr = IGAPointAddVec(pointZ0,FpointZ0,FlocZ0);CHKERRQ(ierr);
				}
				if(pointZ0->atboundary==0 && pointchiUp->atboundary==0)
				{
					IGAElementNextPoint(elemchiUp,pointchiUp);

					ierr = IGAPointGetWorkMat(pointZ0,&KpointZ0);CHKERRQ(ierr);
					ierr = IGAPointGetWorkVec(pointZ0,&FpointZ0);CHKERRQ(ierr);
					//	   Z0sys(IGAPoint p, IGAPoint pChi, IGAPoint pS,IGAPoint pZ, PetscReal *K, PetscReal *F, PetscReal *UChi, PetscReal *S,PetscReal *ZS, void *ctx)
					ierr = Z0sys(pointZ0,pointchiUp,KpointZ0,FpointZ0,Chi0Z0,NULL);CHKERRQ(ierr);
					ierr = IGAPointAddMat(pointZ0,KpointZ0,KlocZ0);CHKERRQ(ierr);
					ierr = IGAPointAddVec(pointZ0,FpointZ0,FlocZ0);CHKERRQ(ierr);
				}
			}
			ierr = IGAElementEndPoint(elemZ0,&pointZ0);CHKERRQ(ierr);
			while (pointchiUp->index != -1)
			{
				IGAElementNextPoint(elemchiUp,pointchiUp);
			}
			ierr = IGAElementEndPoint(elemchiUp,&pointchiUp);CHKERRQ(ierr);
		}

		if(fijaPunto==0)
		{
			ierr = IGAElementFixSystem(elemZ0,KlocZ0,FlocZ0);CHKERRQ(ierr);					//This sets Dirichlet condition ¿? (Yes, this applies the conditions from IGASetBoundaryValue)
		}
		ierr = IGAElementAssembleMat(elemZ0,KlocZ0,KZ0);CHKERRQ(ierr);
		ierr = IGAElementAssembleVec(elemZ0,FlocZ0,FZ0);CHKERRQ(ierr);

	}
	IGANextElement(igachiUp,elemchiUp);

	ierr = IGAEndElement(igaZ0,&elemZ0);CHKERRQ(ierr);
	ierr = IGAEndElement(igachiUp,&elemchiUp);CHKERRQ(ierr);

	// Restore local vectors Chi0 and arrays
	ierr = IGARestoreLocalVecArray(igachiUp,chiUp0,&localChi0Z0,&arrayChi0Z0);CHKERRQ(ierr);

	ierr = MatAssemblyBegin(KZ0,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyEnd  (KZ0,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

	PetscInt n,m;
	PetscInt rows, *cols;
	PetscReal val, *vals;

	if (fijaPunto==1)
	{
		//Here we set values to the Matrix directly, to impose Dirichlet condition in a single point.
		//Note: Lower Left corner is gdl's  0 and 1,
		//		Lower Right corner is gdl's 2*(nx+2)-2 and 2*(nx+2)-1 
		//		Upper Left corner is gdl's  2*(nx+2)*(ny+2)-2*(nx+1)-2 and 2*(nx+2)*(ny+2)-2*(nx+1)-1
		//		Upper Right corner is gdl's 2*(nx+2)*(ny+2)-2 and 2*(nx+2)*(ny+2)-1
		//All of these for when z is a 2nd order nurbs
		
		ierr = MatGetSize(KZ0,&n,&m);CHKERRQ(ierr);
		ierr = MatSetOption(KZ0, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);CHKERRQ(ierr);

		ierr = PetscMalloc1(m,&cols);CHKERRQ(ierr);
		ierr = PetscMalloc1(m,&vals);CHKERRQ(ierr);

		for(int i=0;i<m;i++)
		{
			cols[i]=i;
			vals[i]=0.0;
		}

		rows=0;
		vals[rows]=1.0e6;
		ierr = MatSetValues(KZ0,1,&rows,m,cols,vals,INSERT_VALUES);CHKERRQ(ierr);
		ierr = MatSetValues(KZ0,n,cols,1,&rows,vals,INSERT_VALUES);CHKERRQ(ierr);
		vals[rows]=0.0;

		rows=1;
		vals[rows]=1.0e6;
		ierr = MatSetValues(KZ0,1,&rows,m,cols,vals,INSERT_VALUES);CHKERRQ(ierr);
		ierr = MatSetValues(KZ0,n,cols,1,&rows,vals,INSERT_VALUES);CHKERRQ(ierr);
		vals[rows]=0.0;

		rows=2*(nx+3)*(ny+3)-2; 										//This is for when z is a 3rd order nurb
		//rows=2*(nx+2)*(ny+2)-2*(nx+1)-2; 										//This is for when z is a 2nd order nurb
		//rows=2*(nx+1)*(ny+1)-1;														//This is the dof in x on the upper right corner for 1st order elements
		vals[rows]=1.0e6;
		ierr = MatSetValues(KZ0,1,&rows,m,cols,vals,INSERT_VALUES);CHKERRQ(ierr);
		ierr = MatSetValues(KZ0,n,cols,1,&rows,vals,INSERT_VALUES);CHKERRQ(ierr);
		vals[rows]=0.0;

		ierr = MatAssemblyBegin(KZ0,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
		ierr = MatAssemblyEnd  (KZ0,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

		//Here we set values to the vector directly (to impose Dirichlet condition in a single point)
		//Note: Lower Left corner is gdl's  0 and 1,
		//		Lower Right corner is gdl's 2*(nx+2)-2 and 2*(nx+2)-1 
		//		Upper Left corner is gdl's  2*(nx+2)*(ny+2)-2-2*(nx+1) and 2*(nx+2)*(ny+2)-1-2*(nx+1)
		//		Upper Right corner is gdl's 2*(nx+2)*(ny+2)-2 and 2*(nx+2)*(ny+2)-1
		//All of these for when z is a 2nd order nurbs
		ierr = VecAssemblyBegin(FZ0);CHKERRQ(ierr);
		ierr = VecAssemblyEnd  (FZ0);CHKERRQ(ierr);

		rows=0;
		val=0.0;
		ierr = VecSetValue(FZ0,rows,val,INSERT_VALUES);CHKERRQ(ierr);

		ierr = VecAssemblyBegin(FZ0);CHKERRQ(ierr);
		ierr = VecAssemblyEnd  (FZ0);CHKERRQ(ierr);

		rows=1;
		val=0.0;
		ierr = VecSetValue(FZ0,rows,val,INSERT_VALUES);CHKERRQ(ierr);

		ierr = VecAssemblyBegin(FZ0);CHKERRQ(ierr);
		ierr = VecAssemblyEnd  (FZ0);CHKERRQ(ierr);
		
		rows=2*(nx+3)*(ny+3)-2;
		//rows=2*(nx+2)*(ny+2)-2*(nx+1)-2; 										//This is for when z is a 2nd order nurb
		//rows=2*(nx+1)-1;
		val=0.0;
		ierr = VecSetValue(FZ0,rows,val,INSERT_VALUES);CHKERRQ(ierr);
	}

	ierr = VecAssemblyBegin(FZ0);CHKERRQ(ierr);
	ierr = VecAssemblyEnd  (FZ0);CHKERRQ(ierr);

	ierr = KSPSetOperators(kspZ0,KZ0,KZ0);CHKERRQ(ierr);
	PC pcZ0;
	ierr = KSPGetPC(kspZ0,&pcZ0); CHKERRQ(ierr);
	ierr = PCSetType(pcZ0,PCLU); CHKERRQ(ierr);
	ierr = PCFactorSetMatSolverType(pcZ0,MATSOLVERMUMPS); CHKERRQ(ierr);
	//ierr = KSPSetFromOptions(kspZ0);CHKERRQ(ierr);
	//ierr = KSPSetTolerances(kspZ0,1e-18,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);CHKERRQ(ierr);
	ierr = KSPSolve(kspZ0,FZ0,Z0);CHKERRQ(ierr);

	ierr = KSPDestroy(&kspZ0);CHKERRQ(ierr);
	ierr = MatDestroy(&KZ0);CHKERRQ(ierr);
	ierr = VecDestroy(&FZ0);CHKERRQ(ierr);
	//ierr = IGAWriteVec(igaZ0,z0,"./results/z0-2d-0.dat");CHKERRQ(ierr);
	char nameZ0[512];//="/Z0-2d-0.dat";
	sprintf(nameZ0,"%s%d%s","/Z0-2d-",i,".dat");
	char pathZ0[512];
	sprintf(pathZ0,"%s%s",direct,nameZ0);
	ierr = IGAWriteVec(igaZ0,Z0,pathZ0);CHKERRQ(ierr);	
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
	PetscReal		*KlocVs,*FlocVs;												//AA y BB
	PetscReal		*KpointVs,*FpointVs;											//KKK y FFF
	const PetscReal	*arrayZ0Vs,*arrayChi0Vs,*arrayS0Vs;			//arrayU
	Vec				localZ0Vs,localChi0Vs,localS0Vs;				//localU
	PetscReal		*Chi0Vs,*Z0Vs,*S0Vs;											//U

  	IGAFormSystem	wtfVs;
 	void			*wtf2Vs;

 	KSP kspVs;
	ierr = IGACreateKSP(igaVs,&kspVs);CHKERRQ(ierr);

	// Get local vectors Z0 and Chi0 and arrays
	ierr = IGAGetLocalVecArray(igachiUp,chiUp0,&localChi0Vs,&arrayChi0Vs);CHKERRQ(ierr);
	ierr = IGAGetLocalVecArray(igaZ0,Z0,&localZ0Vs,&arrayZ0Vs);CHKERRQ(ierr);
	ierr = IGAGetLocalVecArray(igaS,s0,&localS0Vs,&arrayS0Vs);CHKERRQ(ierr);

	// Element loop
	ierr = IGABeginElement(igaVs,&elemVs);CHKERRQ(ierr);
	ierr = IGABeginElement(igaZ0,&elemZ0);CHKERRQ(ierr);
	ierr = IGABeginElement(igachiUp,&elemchiUp);CHKERRQ(ierr);
	ierr = IGABeginElement(igaS,&elemS);CHKERRQ(ierr);

	while (IGANextElement(igaVs,elemVs)) 
	{
		IGANextElement(igaZ0,elemZ0);
		IGANextElement(igachiUp,elemchiUp);
		IGANextElement(igaS,elemS);

		ierr = IGAElementGetWorkMat(elemVs,&KlocVs);CHKERRQ(ierr);
		ierr = IGAElementGetWorkVec(elemVs,&FlocVs);CHKERRQ(ierr);
		ierr = IGAElementGetValues(elemZ0,arrayZ0Vs,&Z0Vs);CHKERRQ(ierr);
		ierr = IGAElementGetValues(elemchiUp,arrayChi0Vs,&Chi0Vs);CHKERRQ(ierr);
		ierr = IGAElementGetValues(elemS,arrayS0Vs,&S0Vs);CHKERRQ(ierr);

		// FormSystem loop
		while (IGAElementNextFormSystem(elemVs,&wtfVs,&wtf2Vs)) 
		{
		// Quadrature loop
			ierr = IGAElementBeginPoint(elemVs,&pointVs);CHKERRQ(ierr);
			ierr = IGAElementBeginPoint(elemZ0,&pointZ0);CHKERRQ(ierr);
			ierr = IGAElementBeginPoint(elemchiUp,&pointchiUp);CHKERRQ(ierr);
			ierr = IGAElementBeginPoint(elemS,&pointS);CHKERRQ(ierr);

			while (IGAElementNextPoint(elemVs,pointVs))
			{
				if(pointVs->atboundary==1)
				{

				}
				if(pointVs->atboundary==0 && pointZ0->atboundary==0 && pointchiUp->atboundary==0 && pointS->atboundary==0)
				{
					IGAElementNextPoint(elemZ0,pointZ0);
					IGAElementNextPoint(elemchiUp,pointchiUp);
					IGAElementNextPoint(elemS,pointS);

					ierr = IGAPointGetWorkMat(pointVs,&KpointVs);CHKERRQ(ierr);
					ierr = IGAPointGetWorkVec(pointVs,&FpointVs);CHKERRQ(ierr);
					//VS(IGAPoint p,IGAPoint pChi,IGAPoint pZu,IGAPoint pS,PetscReal *K,PetscReal *F, PetscReal *Chi,PetscReal *Zu,PetscReal *S,void *ctx)
					ierr = VS(pointVs,pointchiUp,pointZ0,pointS,KpointVs,FpointVs,Chi0Vs,Z0Vs,S0Vs,NULL);CHKERRQ(ierr);
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
			ierr = IGAElementEndPoint(elemVs,&pointVs);CHKERRQ(ierr);
			ierr = IGAElementEndPoint(elemZ0,&pointZ0);CHKERRQ(ierr);
			ierr = IGAElementEndPoint(elemchiUp,&pointchiUp);CHKERRQ(ierr);
			ierr = IGAElementEndPoint(elemS,&pointS);CHKERRQ(ierr);
		}

		//ierr = IGAElementFixSystem(elemStress,KlocStress,FlocStress);CHKERRQ(ierr);					//This sets Dirichlet condition ¿? (Yes, this applies the conditions from IGASetBoundaryValue)
		ierr = IGAElementAssembleMat(elemVs,KlocVs,KVs);CHKERRQ(ierr);
		ierr = IGAElementAssembleVec(elemVs,FlocVs,FVs);CHKERRQ(ierr);
	}
	IGANextElement(igaZ0,elemZ0);
	IGANextElement(igachiUp,elemchiUp);
	IGANextElement(igaS,elemS);

	ierr = IGAEndElement(igaVs,&elemVs);CHKERRQ(ierr);
	ierr = IGAEndElement(igaZ0,&elemZ0);CHKERRQ(ierr);
	ierr = IGAEndElement(igachiUp,&elemchiUp);CHKERRQ(ierr);
	ierr = IGAEndElement(igaS,&elemS);CHKERRQ(ierr);

	// Restore local vectors u, Z0, Chi0 and arrays
	ierr = IGARestoreLocalVecArray(igaZ0,Z0,&localZ0Vs,&arrayZ0Vs);CHKERRQ(ierr);
	ierr = IGARestoreLocalVecArray(igachiUp,chiUp0,&localChi0Vs,&arrayChi0Vs);CHKERRQ(ierr);
	ierr = IGARestoreLocalVecArray(igaS,s0,&localS0Vs,&arrayS0Vs);CHKERRQ(ierr);
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
	char pathVs[512];
	sprintf(pathVs,"%s%s",direct,nameVs);
	ierr = IGAWriteVec(igaVs,Vs0,pathVs);CHKERRQ(ierr);	
//

//System for L2 protection of smoothed V^{S}
	PetscPrintf(PETSC_COMM_WORLD,"\nSystem for Smooth Vs starting \n\n");
	T=time(NULL);
	tm=*localtime(&T);
	PetscPrintf(PETSC_COMM_WORLD,"Current time is %02d:%02d:%02d \n",tm.tm_hour,tm.tm_min,tm.tm_sec);


	//First integrate velocity for both defects
	Vec FVsInt, FVsXi, Int1sVec, Int2sVec, IntVec;
	ierr = IGACreateVec(igaVs,&FVsInt);CHKERRQ(ierr);
	ierr = IGACreateVec(igaVs,&FVsXi);CHKERRQ(ierr);
	ierr = IGACreateVec(igaVs,&Int1sVec);CHKERRQ(ierr);
	ierr = IGACreateVec(igaVs,&Int2sVec);CHKERRQ(ierr);
	ierr = IGACreateVec(igaVs,&IntVec);CHKERRQ(ierr);
																																		//Fix these comments
	PetscReal		*PointInt1a,*PointInt2a,*PointInt,*locInt1a,*locInt2a,*locInt;														//AA y BB
	const PetscReal	*arrayVs0,*arraySVs;																								//arrayU
	Vec				localVs0,localSVs;																								//localU
	PetscReal 		*SVs,*Vs0Values,Int1a=0.0,Int2a=0.0,IntS=0.0;

	IGAFormSystem	wtfVsInt;
	void			*wtf2VsInt;

	// Get local vectors V0 and InputAl and arrays
	ierr = IGAGetLocalVecArray(igaVs,Vs0,&localVs0,&arrayVs0);CHKERRQ(ierr);
	ierr = IGAGetLocalVecArray(igaS,s0,&localSVs,&arraySVs);CHKERRQ(ierr);

	// Element loop
	ierr = IGABeginElement(igaVs,&elemVs);CHKERRQ(ierr);
	ierr = IGABeginElement(igaS,&elemS);CHKERRQ(ierr);

	while (IGANextElement(igaVs,elemVs))
	{
		IGANextElement(igaS,elemS);

		ierr = IGAElementGetWorkVec(elemVs,&locInt1a);CHKERRQ(ierr);
		ierr = IGAElementGetWorkVec(elemVs,&locInt2a);CHKERRQ(ierr);
		ierr = IGAElementGetWorkVec(elemVs,&locInt);CHKERRQ(ierr);
		ierr = IGAElementGetValues(elemVs,arrayVs0,&Vs0Values);CHKERRQ(ierr);
		ierr = IGAElementGetValues(elemS,arraySVs,&SVs);CHKERRQ(ierr);

		// FormSystem loop
		while (IGAElementNextFormSystem(elemVs,&wtfVsInt,&wtf2VsInt))
		{
		// Quadrature loop
			ierr = IGAElementBeginPoint(elemVs,&pointVs);CHKERRQ(ierr);
			ierr = IGAElementBeginPoint(elemS,&pointS);CHKERRQ(ierr);
			
			while (IGAElementNextPoint(elemVs,pointVs))
			{
				if(pointVs->atboundary==1)
				{

				}
				if(pointVs->atboundary==0 && pointS->atboundary==0)
				{
					IGAElementNextPoint(elemS,pointS);
					ierr = IGAPointGetWorkVec(pointVs,&PointInt1a);CHKERRQ(ierr);
					ierr = IGAPointGetWorkVec(pointVs,&PointInt2a);CHKERRQ(ierr);
					ierr = IGAPointGetWorkVec(pointVs,&PointInt);CHKERRQ(ierr);
					//Int_Xi_Vs(IGAPoint pV,IGAPoint pS,PetscReal *FInt1a,PetscReal *FInt2a,PetscReal *FInt,PetscReal *VS,PetscReal *US,void *ctx)
					ierr = Int_Xi_Vs(pointVs,pointS,PointInt1a,PointInt2a,PointInt,Vs0Values,SVs,NULL);CHKERRQ(ierr);
					ierr = IGAPointAddVec(pointVs,PointInt1a,locInt1a);CHKERRQ(ierr);
					ierr = IGAPointAddVec(pointVs,PointInt2a,locInt2a);CHKERRQ(ierr);
					ierr = IGAPointAddVec(pointVs,PointInt,locInt);CHKERRQ(ierr);
				}
			}
			while (pointS->index != -1)
			{
				IGAElementNextPoint(elemS,pointS);
			}
			ierr = IGAElementEndPoint(elemVs,&pointVs);CHKERRQ(ierr);
			ierr = IGAElementEndPoint(elemS,&pointS);CHKERRQ(ierr);
		}
		ierr = IGAElementAssembleVec(elemVs,locInt1a,Int1sVec);CHKERRQ(ierr);
		ierr = IGAElementAssembleVec(elemVs,locInt2a,Int2sVec);CHKERRQ(ierr);
		ierr = IGAElementAssembleVec(elemVs,locInt,IntVec);CHKERRQ(ierr);
	}
	IGANextElement(igaS,elemS);

	ierr = IGAEndElement(igaVs,&elemVs);CHKERRQ(ierr);
	ierr = IGAEndElement(igaS,&elemS);CHKERRQ(ierr);

	ierr = IGARestoreLocalVecArray(igaS,s0,&localSVs,&arraySVs);CHKERRQ(ierr);
	ierr = IGARestoreLocalVecArray(igaVs,Vs0,&localVs0,&arrayVs0);CHKERRQ(ierr);

	ierr = VecAssemblyBegin(Int1sVec);CHKERRQ(ierr);
	ierr = VecAssemblyEnd  (Int1sVec);CHKERRQ(ierr);
	ierr = VecAssemblyBegin(Int2sVec);CHKERRQ(ierr);
	ierr = VecAssemblyEnd  (Int2sVec);CHKERRQ(ierr);

	ierr = VecAssemblyBegin(IntVec);CHKERRQ(ierr);
	ierr = VecAssemblyEnd  (IntVec);CHKERRQ(ierr);

	//Here add all values at Gauss points
	ierr = VecSum(Int1sVec,&Int1a);CHKERRQ(ierr);
	ierr = VecSum(Int2sVec,&Int2a);CHKERRQ(ierr);
	ierr = VecSum(IntVec,&IntS);CHKERRQ(ierr);

	ierr = PetscPrintf(PETSC_COMM_WORLD,"Int1a=%f, Int2a=%f, IntS=%f \n",Int1a,Int2a,IntS);CHKERRQ(ierr);

	//Now put those velocities on each defect location
	Vec FVsSmooth,VsSmooth;
	ierr = IGACreateVec(igaVs,&FVsSmooth);CHKERRQ(ierr);
	ierr = IGACreateVec(igaVs,&VsSmooth);CHKERRQ(ierr);

	PetscReal		*VsSmoothPoint,*locSmooth;										//AA y BB

	PetscReal Smin,Smax,absSmax;

	ierr = VecMin(s0,NULL,&Smin);
	ierr = VecMax(s0,NULL,&Smax);
	absSmax=fmax(fabs(Smin),fabs(Smax));

	ierr = PetscPrintf(PETSC_COMM_WORLD,"Smin=%f, Smax=%f, absSmax=%f \n",Smin,Smax,absSmax);CHKERRQ(ierr);

	// Get local vectors V0 and InputAl and arrays
	ierr = IGAGetLocalVecArray(igaVs,Vs0,&localVs0,&arrayVs0);CHKERRQ(ierr);
	ierr = IGAGetLocalVecArray(igaS,s0,&localSVs,&arraySVs);CHKERRQ(ierr);

	// Element loop
	ierr = IGABeginElement(igaVs,&elemVs);CHKERRQ(ierr);
	ierr = IGABeginElement(igaS,&elemS);CHKERRQ(ierr);

	while (IGANextElement(igaVs,elemVs))
	{
		IGANextElement(igaS,elemS);

		ierr = IGAElementGetWorkVec(elemVs,&locSmooth);CHKERRQ(ierr);
		ierr = IGAElementGetValues(elemVs,arrayVs0,&Vs0Values);CHKERRQ(ierr);
		ierr = IGAElementGetValues(elemS,arraySVs,&SVs);CHKERRQ(ierr);

		// FormSystem loop
		while (IGAElementNextFormSystem(elemVs,&wtfVsInt,&wtf2VsInt))
		{
		// Quadrature loop
			ierr = IGAElementBeginPoint(elemVs,&pointVs);CHKERRQ(ierr);
			ierr = IGAElementBeginPoint(elemS,&pointS);CHKERRQ(ierr);
			
			while (IGAElementNextPoint(elemVs,pointVs))
			{
				if(pointVs->atboundary==1)
				{

				}
				if(pointVs->atboundary==0 && pointS->atboundary==0)
				{
					IGAElementNextPoint(elemS,pointS);
					ierr = IGAPointGetWorkVec(pointVs,&VsSmoothPoint);CHKERRQ(ierr);
					//ProjVS(IGAPoint pV,IGAPoint pS,PetscReal *FS1,PetscReal *US,PetscReal Vs1,PetscReal Vs2,PetscReal IntS,void *ctx)
					ierr = ProjVS(pointVs,pointS,VsSmoothPoint,SVs,Int1a,Int2a,IntS,absSmax,NULL);CHKERRQ(ierr);
					ierr = IGAPointAddVec(pointVs,VsSmoothPoint,locSmooth);CHKERRQ(ierr);
				}
			}
			while (pointS->index != -1)
			{
				IGAElementNextPoint(elemS,pointS);
			}
			ierr = IGAElementEndPoint(elemVs,&pointVs);CHKERRQ(ierr);
			ierr = IGAElementEndPoint(elemS,&pointS);CHKERRQ(ierr);
		}
		ierr = IGAElementAssembleVec(elemVs,locSmooth,FVsSmooth);CHKERRQ(ierr);
	}
	IGANextElement(igaS,elemS);

	ierr = IGAEndElement(igaVs,&elemVs);CHKERRQ(ierr);
	ierr = IGAEndElement(igaS,&elemS);CHKERRQ(ierr);

	ierr = IGARestoreLocalVecArray(igaS,s0,&localSVs,&arraySVs);CHKERRQ(ierr);
	ierr = IGARestoreLocalVecArray(igaVs,Vs0,&localVs0,&arrayVs0);CHKERRQ(ierr);

	//Here we have a vector with an L2 Projection of the integrated velocity.
	ierr = VecAssemblyBegin(FVsSmooth);CHKERRQ(ierr);
	ierr = VecAssemblyEnd  (FVsSmooth);CHKERRQ(ierr);

	//Here we solved the L2 projection system
	ierr = KSPSolve(kspVs,FVsSmooth,VsSmooth);CHKERRQ(ierr);

	char nameVsS1[512];//="/VsSmooth-2d-0.dat";
	sprintf(nameVsS1,"%s%d%s","/VsSmooth-2d-",i,".dat");
	char pathVsS1[512];
	sprintf(pathVsS1,"%s%s",direct,nameVsS1);
	ierr = IGAWriteVec(igaVs,VsSmooth,pathVsS1);CHKERRQ(ierr);

	//ierr = KSPDestroy(&kspVs);CHKERRQ(ierr);
	//ierr = MatDestroy(&KVs);CHKERRQ(ierr);
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
	const PetscReal	*arrayZ0Va,*arrayChi0Va,*arrayAl0Va,*arraySVa;					//arrayU
	Vec				localZ0Va,localChi0Va,localAl0Va,localSVa;						//localU
	PetscReal		*Chi0Va,*Z0Va,*Al0Va,*SVa;										//U

	IGAFormSystem	wtfVa;
	void			*wtf2Va;

 	KSP kspVa;
	ierr = IGACreateKSP(igaVa,&kspVa);CHKERRQ(ierr);

	// Get local vectors Z0 and Chi0 and arrays
	ierr = IGAGetLocalVecArray(igachiUp,chiUp0,&localChi0Va,&arrayChi0Va);CHKERRQ(ierr);
	ierr = IGAGetLocalVecArray(igaZ0,Z0,&localZ0Va,&arrayZ0Va);CHKERRQ(ierr);
	ierr = IGAGetLocalVecArray(igaAl,alInput,&localAl0Va,&arrayAl0Va);CHKERRQ(ierr);
	ierr = IGAGetLocalVecArray(igaS,s0,&localSVa,&arraySVa);CHKERRQ(ierr);

	// Element loop
	ierr = IGABeginElement(igaVa,&elemVa);CHKERRQ(ierr);
	ierr = IGABeginElement(igaZ0,&elemZ0);CHKERRQ(ierr);
	ierr = IGABeginElement(igachiUp,&elemchiUp);CHKERRQ(ierr);
	ierr = IGABeginElement(igaAl,&elemAlp);CHKERRQ(ierr);
	ierr = IGABeginElement(igaS,&elemS);CHKERRQ(ierr);

	while (IGANextElement(igaVa,elemVa)) 
	{
		IGANextElement(igaZ0,elemZ0);
		IGANextElement(igachiUp,elemchiUp);
		IGANextElement(igaAl,elemAlp);
		IGANextElement(igaS,elemS);

		ierr = IGAElementGetWorkMat(elemVa,&KlocVa);CHKERRQ(ierr);
		ierr = IGAElementGetWorkVec(elemVa,&FlocVa);CHKERRQ(ierr);
		ierr = IGAElementGetValues(elemZ0,arrayZ0Va,&Z0Va);CHKERRQ(ierr);
		ierr = IGAElementGetValues(elemchiUp,arrayChi0Va,&Chi0Va);CHKERRQ(ierr);
		ierr = IGAElementGetValues(elemAlp,arrayAl0Va,&Al0Va);CHKERRQ(ierr);
		ierr = IGAElementGetValues(elemS,arraySVa,&SVa);CHKERRQ(ierr);

		// FormSystem loop
		while (IGAElementNextFormSystem(elemVa,&wtfVa,&wtf2Va)) 
		{
		// Quadrature loop
			ierr = IGAElementBeginPoint(elemVa,&pointVa);CHKERRQ(ierr);
			ierr = IGAElementBeginPoint(elemZ0,&pointZ0);CHKERRQ(ierr);
			ierr = IGAElementBeginPoint(elemchiUp,&pointchiUp);CHKERRQ(ierr);
			ierr = IGAElementBeginPoint(elemAlp,&pointAlp);CHKERRQ(ierr);
			ierr = IGAElementBeginPoint(elemS,&pointS);CHKERRQ(ierr);

			while (IGAElementNextPoint(elemVa,pointVa))
			{
				if(pointVa->atboundary==1)
				{

				}
				if(pointVa->atboundary==0 && pointZ0->atboundary==0 && pointchiUp->atboundary==0 && pointAlp->atboundary==0 && pointS->atboundary==0)
				{
					IGAElementNextPoint(elemZ0,pointZ0);
					IGAElementNextPoint(elemchiUp,pointchiUp);
					IGAElementNextPoint(elemAlp,pointAlp);
					IGAElementNextPoint(elemS,pointS);

					ierr = IGAPointGetWorkMat(pointVa,&KpointVa);CHKERRQ(ierr);
					ierr = IGAPointGetWorkVec(pointVa,&FpointVa);CHKERRQ(ierr);
					//Valpha(IGAPoint p,IGAPoint pChi,IGAPoint pZu,IGAPoint pAl,IGAPoint pS,PetscReal *K,PetscReal *F, PetscReal *Chi,PetscReal *Zu,PetscReal *US,PetscReal *UAl,void *ctx)
					ierr = Valpha(pointVa,pointchiUp,pointZ0,pointAlp,pointS,KpointVa,FpointVa,Chi0Va,Z0Va,SVa,Al0Va,NULL);CHKERRQ(ierr);
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
			while (pointS->index != -1)
			{
				IGAElementNextPoint(elemS,pointS);
			}
			ierr = IGAElementEndPoint(elemVa,&pointVa);CHKERRQ(ierr);
			ierr = IGAElementEndPoint(elemZ0,&pointZ0);CHKERRQ(ierr);
			ierr = IGAElementEndPoint(elemchiUp,&pointchiUp);CHKERRQ(ierr);
			ierr = IGAElementEndPoint(elemAlp,&pointAlp);CHKERRQ(ierr);
			ierr = IGAElementEndPoint(elemS,&pointS);CHKERRQ(ierr);
		}

		//ierr = IGAElementFixSystem(elemStress,KlocStress,FlocStress);CHKERRQ(ierr);					//This sets Dirichlet condition ¿? (Yes, this applies the conditions from IGASetBoundaryValue)
		ierr = IGAElementAssembleMat(elemVa,KlocVa,KVa);CHKERRQ(ierr);
		ierr = IGAElementAssembleVec(elemVa,FlocVa,FVa);CHKERRQ(ierr);
	}
	IGANextElement(igaZ0,elemZ0);
	IGANextElement(igachiUp,elemchiUp);
	IGANextElement(igaAl,elemAlp);
	IGANextElement(igaS,elemS);

	ierr = IGAEndElement(igaVa,&elemVa);CHKERRQ(ierr);
	ierr = IGAEndElement(igaZ0,&elemZ0);CHKERRQ(ierr);
	ierr = IGAEndElement(igachiUp,&elemchiUp);CHKERRQ(ierr);
	ierr = IGAEndElement(igaAl,&elemAlp);CHKERRQ(ierr);
	ierr = IGAEndElement(igaS,&elemS);CHKERRQ(ierr);

	// Restore local vectors u, Z0, Chi0 and arrays
	ierr = IGARestoreLocalVecArray(igaZ0,Z0,&localZ0Va,&arrayZ0Va);CHKERRQ(ierr);
	ierr = IGARestoreLocalVecArray(igachiUp,chiUp0,&localChi0Va,&arrayChi0Va);CHKERRQ(ierr);
	ierr = IGARestoreLocalVecArray(igaAl,alInput,&localAl0Va,&arrayAl0Va);CHKERRQ(ierr);
	ierr = IGARestoreLocalVecArray(igaS,s0,&localSVa,&arraySVa);CHKERRQ(ierr);
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

	char nameVa[]="/Va-2d-0.dat";
	char pathVa[512];
	sprintf(pathVa,"%s%s",direct,nameVa);
	ierr = IGAWriteVec(igaVa,Va0,pathVa);CHKERRQ(ierr);
//

//System for evolution of S (debug and ask \dot{S}*n=0 in boundaries where Vs*n=0)
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
			//ierr = IGASetBoundaryForm(igachiS,dir,side,PETSC_TRUE);CHKERRQ(ierr);  				// Neumann boundary conditions
		}
	}
	ierr = IGASetUp(igaSdot);CHKERRQ(ierr);

	Mat KSdot;
	Vec Sdot,FSdot;

	ierr = IGACreateMat(igaSdot,&KSdot);CHKERRQ(ierr);
	ierr = IGACreateVec(igaSdot,&Sdot);CHKERRQ(ierr);
	ierr = IGACreateVec(igaSdot,&FSdot);CHKERRQ(ierr);

	IGAPoint		pointSdot;													//point
	IGAElement		elemSdot;													//element
	PetscReal		*KlocSdot,*FlocSdot;										//AA y BB
	PetscReal		*KpointSdot,*FpointSdot;									//KKK y FFF
	const PetscReal	*arrayVsSdot,*arrayS0Sdot;									//arrayU
	Vec				localVsSdot,localS0Sdot;									//localU
	PetscReal		*VsSdot,*S0Sdot;											//U

  	IGAFormSystem	wtfSdot;
 	void			*wtf2Sdot;

 	KSP kspSdot;
	ierr = IGACreateKSP(igaSdot,&kspSdot);CHKERRQ(ierr);

	// Get local vectors Z0 and Chi0 and arrays
	ierr = IGAGetLocalVecArray(igaVs,VsSmooth,&localVsSdot,&arrayVsSdot);CHKERRQ(ierr);
	ierr = IGAGetLocalVecArray(igaS,s0,&localS0Sdot,&arrayS0Sdot);CHKERRQ(ierr);

	// Element loop
	ierr = IGABeginElement(igaVs,&elemVs);CHKERRQ(ierr);
	ierr = IGABeginElement(igaS,&elemS);CHKERRQ(ierr);
	ierr = IGABeginElement(igaSdot,&elemSdot);CHKERRQ(ierr);

	while (IGANextElement(igaSdot,elemSdot)) 
	{
		IGANextElement(igaVs,elemVs);
		IGANextElement(igaS,elemS);

		ierr = IGAElementGetWorkMat(elemSdot,&KlocSdot);CHKERRQ(ierr);
		ierr = IGAElementGetWorkVec(elemSdot,&FlocSdot);CHKERRQ(ierr);
		ierr = IGAElementGetValues(elemVs,arrayVsSdot,&VsSdot);CHKERRQ(ierr);
		ierr = IGAElementGetValues(elemS,arrayS0Sdot,&S0Sdot);CHKERRQ(ierr);

		// FormSystem loop
		while (IGAElementNextFormSystem(elemSdot,&wtfSdot,&wtf2Sdot)) 
		{
		// Quadrature loop
			ierr = IGAElementBeginPoint(elemSdot,&pointSdot);CHKERRQ(ierr);
			ierr = IGAElementBeginPoint(elemVs,&pointVs);CHKERRQ(ierr);
			ierr = IGAElementBeginPoint(elemS,&pointS);CHKERRQ(ierr);

			while (IGAElementNextPoint(elemSdot,pointSdot))
			{
				if(pointSdot->atboundary==1)
				{
					//Code boundary condition, will be required for general case
				}
				if(pointSdot->atboundary==0 && pointVs->atboundary==0 && pointS->atboundary==0)
				{
					IGAElementNextPoint(elemVs,pointVs);
					IGAElementNextPoint(elemS,pointS);

					ierr = IGAPointGetWorkMat(pointSdot,&KpointSdot);CHKERRQ(ierr);
					ierr = IGAPointGetWorkVec(pointSdot,&FpointSdot);CHKERRQ(ierr);
					//SdotFunc(IGAPoint p, IGAPoint pVs, IGAPoint pS, PetscReal *K, PetscReal *F, PetscReal *UVs, PetscReal *US, void *ctx)
					ierr = SdotFunc(pointSdot,pointVs,pointS,KpointSdot,FpointSdot,VsSdot,S0Sdot,&user);CHKERRQ(ierr);
					ierr = IGAPointAddMat(pointSdot,KpointSdot,KlocSdot);CHKERRQ(ierr);
					ierr = IGAPointAddVec(pointSdot,FpointSdot,FlocSdot);CHKERRQ(ierr);
				}
			}
			while (pointVs->index != -1)
			{
				IGAElementNextPoint(elemVs,pointVs);
			}
			
			while (pointS->index != -1)
			{
				IGAElementNextPoint(elemS,pointS);
			}
			ierr = IGAElementEndPoint(elemSdot,&pointSdot);CHKERRQ(ierr);
			ierr = IGAElementEndPoint(elemVs,&pointVs);CHKERRQ(ierr);
			ierr = IGAElementEndPoint(elemS,&pointS);CHKERRQ(ierr);
		}

		//ierr = IGAElementFixSystem(elemSdot,KlocSdot,FlocSdot);CHKERRQ(ierr);					//This sets Dirichlet condition ¿? (Yes, this applies the conditions from IGASetBoundaryValue)
		ierr = IGAElementAssembleMat(elemSdot,KlocSdot,KSdot);CHKERRQ(ierr);
		ierr = IGAElementAssembleVec(elemSdot,FlocSdot,FSdot);CHKERRQ(ierr);
	}
	IGANextElement(igaVs,elemVs);
	IGANextElement(igaS,elemS);

	ierr = IGAEndElement(igaSdot,&elemSdot);CHKERRQ(ierr);
	ierr = IGAEndElement(igaVs,&elemVs);CHKERRQ(ierr);
	ierr = IGAEndElement(igaS,&elemS);CHKERRQ(ierr);

	// Restore local vectors u, Z0, Chi0 and arrays
	ierr = IGARestoreLocalVecArray(igaVs,VsSmooth,&localVsSdot,&arrayVsSdot);CHKERRQ(ierr);
	ierr = IGARestoreLocalVecArray(igaS,s0,&localS0Sdot,&arrayS0Sdot);CHKERRQ(ierr);
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
		//PetscInt n;
		ierr = VecGetSize(Sdot,&n);CHKERRQ(ierr);

		PetscInt *colsS;
		PetscReal *valsS;

		ierr = PetscMalloc1(n,&colsS);CHKERRQ(ierr);
		ierr = PetscMalloc1(n,&valsS);CHKERRQ(ierr);

		for(int i=0;i<n;i=i+8)
		{
			colsS[i]=i;		colsS[i+1]=-1;		colsS[i+2]=i+2;		colsS[i+3]=-1;		colsS[i+4]=i+4;		colsS[i+5]=-1;		colsS[i+6]=i+6;		colsS[i+7]=-1;
			valsS[i]=0.0;	valsS[i+1]=-1.0;	valsS[i+2]=0.0;		valsS[i+3]=-1.0;	valsS[i+4]=0.0;		valsS[i+5]=-1.0;	valsS[i+6]=0.0;		valsS[i+7]=-1.0;
		}

		//Here we set values to the vector directly (to impose Dirichlet condition in a single point)
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
	char pathSdot[512];
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
	ierr = IGAGetLocalVecArray(igaZ0,Z0,&localZ0Gradz,&arrayZ0Gradz);CHKERRQ(ierr);

	//Element loop
	ierr = IGABeginElement(igaGradz,&elemGradz);CHKERRQ(ierr);
	ierr = IGABeginElement(igaZ0,&elemZ0);CHKERRQ(ierr);

	while (IGANextElement(igaGradz,elemGradz))
	{
		IGANextElement(igaZ0,elemZ0);

		ierr = IGAElementGetWorkMat(elemGradz,&KlocGradz);CHKERRQ(ierr);
		ierr = IGAElementGetWorkVec(elemGradz,&FlocGradz);CHKERRQ(ierr);
		ierr = IGAElementGetValues(elemZ0,arrayZ0Gradz,&Z0Gradz);CHKERRQ(ierr);

		//FormSystem loop
		while (IGAElementNextFormSystem(elemGradz,&wtfGradz,&wtf2Gradz)) 
		{
			//Quadrature loop
			ierr = IGAElementBeginPoint(elemGradz,&pointGradz);CHKERRQ(ierr);
			ierr = IGAElementBeginPoint(elemZ0,&pointZ0);CHKERRQ(ierr);

			while (IGAElementNextPoint(elemGradz,pointGradz)) 
			{
				if(pointGradz->atboundary==0 && pointZ0->atboundary==0)
				{
					IGAElementNextPoint(elemZ0,pointZ0);
					ierr = IGAPointGetWorkMat(pointGradz,&KpointGradz);CHKERRQ(ierr);
					ierr = IGAPointGetWorkVec(pointGradz,&FpointGradz);CHKERRQ(ierr);
					ierr = gradz(pointGradz,pointZ0,KpointGradz,FpointGradz,Z0Gradz,NULL);CHKERRQ(ierr);
					ierr = IGAPointAddMat(pointGradz,KpointGradz,KlocGradz);CHKERRQ(ierr);
					ierr = IGAPointAddVec(pointGradz,FpointGradz,FlocGradz);CHKERRQ(ierr);
				}
			}
			IGAElementNextPoint(elemZ0,pointZ0);

			ierr = IGAElementEndPoint(elemGradz,&pointGradz);CHKERRQ(ierr);
			while (pointZ0->index != -1)
			{
				IGAElementNextPoint(elemZ0,pointZ0);
			}
			ierr = IGAElementEndPoint(elemZ0,&pointZ0);CHKERRQ(ierr);
		}
		ierr = IGAElementAssembleMat(elemGradz,KlocGradz,KGradz0);CHKERRQ(ierr);
		ierr = IGAElementAssembleVec(elemGradz,FlocGradz,FGradz0);CHKERRQ(ierr);

	}
	IGANextElement(igaZ0,elemZ0);

	ierr = IGAEndElement(igaGradz,&elemGradz);CHKERRQ(ierr);
	ierr = IGAEndElement(igaZ0,&elemZ0);CHKERRQ(ierr);

	// Restore local vectors s0 and arrays
	ierr = IGARestoreLocalVecArray(igaZ0,Z0,&localZ0Gradz,&arrayZ0Gradz);CHKERRQ(ierr);

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
	char pathGradz0[512];
	sprintf(pathGradz0,"%s%s",direct,nameGradz0);
	ierr = IGAWriteVec(igaGradz,gradZ0,pathGradz0);CHKERRQ(ierr);
//

//System for L2 projection of Up
	PetscPrintf(PETSC_COMM_WORLD,"\nSystem for U^p starting \n\n");
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
	ierr = IGAGetLocalVecArray(igaZ0,Z0,&localZ0Up,&arrayZ0Up);CHKERRQ(ierr);
	ierr = IGAGetLocalVecArray(igachiUp,chiUp0,&localChiUp,&arrayChiUp);CHKERRQ(ierr);

	//Element loop
	ierr = IGABeginElement(igaUp,&elemUp);CHKERRQ(ierr);
	ierr = IGABeginElement(igaZ0,&elemZ0);CHKERRQ(ierr);
	ierr = IGABeginElement(igachiUp,&elemchiUp);CHKERRQ(ierr);

	while (IGANextElement(igaUp,elemUp))
	{
		IGANextElement(igaZ0,elemZ0);
		IGANextElement(igachiUp,elemchiUp);

		ierr = IGAElementGetWorkMat(elemUp,&KlocUp);CHKERRQ(ierr);
		ierr = IGAElementGetWorkVec(elemUp,&FlocUp);CHKERRQ(ierr);
		ierr = IGAElementGetValues(elemZ0,arrayZ0Up,&Z0Up);CHKERRQ(ierr);
		ierr = IGAElementGetValues(elemUp,arrayChiUp,&ChiUp);CHKERRQ(ierr);

		//FormSystem loop
		while (IGAElementNextFormSystem(elemUp,&wtfUp,&wtf2Up)) 
		{
			//Quadrature loop
			ierr = IGAElementBeginPoint(elemUp,&pointUp);CHKERRQ(ierr);
			ierr = IGAElementBeginPoint(elemZ0,&pointZ0);CHKERRQ(ierr);
			ierr = IGAElementBeginPoint(elemchiUp,&pointchiUp);CHKERRQ(ierr);

			while (IGAElementNextPoint(elemUp,pointUp)) 
			{
				if(pointUp->atboundary==0 && pointZ0->atboundary==0 && pointchiUp==0)
				{
					IGAElementNextPoint(elemZ0,pointZ0);
					IGAElementNextPoint(elemchiUp,pointchiUp);

					ierr = IGAPointGetWorkMat(pointUp,&KpointUp);CHKERRQ(ierr);
					ierr = IGAPointGetWorkVec(pointUp,&FpointUp);CHKERRQ(ierr);
					ierr = SysUp(pointUp,pointZ0,pointchiUp,KpointUp,FpointUp,Z0Up,ChiUp,NULL);CHKERRQ(ierr);
					ierr = IGAPointAddMat(pointUp,KpointUp,KlocUp);CHKERRQ(ierr);
					ierr = IGAPointAddVec(pointUp,FpointUp,FlocUp);CHKERRQ(ierr);
				}
			}
			IGAElementNextPoint(elemZ0,pointZ0);
			IGAElementNextPoint(elemchiUp,pointchiUp);
		
			ierr = IGAElementEndPoint(elemUp,&pointUp);CHKERRQ(ierr);
			while (pointZ0->index != -1)
			{
				IGAElementNextPoint(elemZ0,pointZ0);
			}
			ierr = IGAElementEndPoint(elemZ0,&pointZ0);CHKERRQ(ierr);
			while (pointchiUp->index != -1)
			{
				IGAElementNextPoint(elemchiUp,pointchiUp);
			}
			ierr = IGAElementEndPoint(elemchiUp,&pointchiUp);CHKERRQ(ierr);
		}
		ierr = IGAElementAssembleMat(elemUp,KlocUp,KUp0);CHKERRQ(ierr);
		ierr = IGAElementAssembleVec(elemUp,FlocUp,FUp0);CHKERRQ(ierr);

	}
	IGANextElement(igaZ0,elemZ0);
	IGANextElement(igachiUp,elemchiUp);

	ierr = IGAEndElement(igaUp,&elemUp);CHKERRQ(ierr);
	ierr = IGAEndElement(igaZ0,&elemZ0);CHKERRQ(ierr);
	ierr = IGAEndElement(igachiUp,&elemchiUp);CHKERRQ(ierr);

	// Restore local vectors and arrays
	ierr = IGARestoreLocalVecArray(igaZ0,Z0,&localZ0Up,&arrayZ0Up);CHKERRQ(ierr);
	ierr = IGARestoreLocalVecArray(igachiUp,chiUp0,&localChiUp,&arrayChiUp);CHKERRQ(ierr);

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
	char pathUp0[512];
	sprintf(pathUp0,"%s%s",direct,nameUp0);
	ierr = IGAWriteVec(igaUp,Up0,pathUp0);CHKERRQ(ierr);
//

//System for general L2 projection (debug, but not really important, just for looking at fields)
	PetscPrintf(PETSC_COMM_WORLD,"\nSystem for general L2 projection starting \n\n");
	PetscPrintf(PETSC_COMM_WORLD,"\nCurrently projecting alpha X V^a + SV^s \n\n");
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
	IGAElement		elemProj;											//element
	PetscReal		*KlocProj,*FlocProj;								//AA y BB
	PetscReal		*KpointProj,*FpointProj;							//KKK y FFF
	const PetscReal *arrayAlProj,*arraySProj,*arrayVaProj;				//arrayU
	Vec				localAlProj,localSProj,localVaProj;					//localU
	PetscReal		*AlProj,*SProj,*VaProj;								//U

  	IGAFormSystem  wtfProj;
 	void           *wtf2Proj;

 	KSP kspProj;
	ierr = IGACreateKSP(igaProj,&kspProj);CHKERRQ(ierr);

	//Get local vectors and arrays
	ierr = IGAGetLocalVecArray(igaAl,alp0,&localAlProj,&arrayAlProj);CHKERRQ(ierr);
	ierr = IGAGetLocalVecArray(igaS,s0,&localSProj,&arraySProj);CHKERRQ(ierr);
	ierr = IGAGetLocalVecArray(igaVa,Va0,&localVaProj,&arrayVaProj);CHKERRQ(ierr);

	//Element loop
	ierr = IGABeginElement(igaProj,&elemProj);CHKERRQ(ierr);
	ierr = IGABeginElement(igaAl,&elemAlp);CHKERRQ(ierr);
	ierr = IGABeginElement(igaS,&elemS);CHKERRQ(ierr);
	ierr = IGABeginElement(igaVa,&elemVa);CHKERRQ(ierr);

	while (IGANextElement(igaProj,elemProj))
	{
		IGANextElement(igaAl,elemAlp);
		IGANextElement(igaS,elemS);
		IGANextElement(igaVa,elemVa);


		ierr = IGAElementGetWorkMat(elemProj,&KlocProj);CHKERRQ(ierr);
		ierr = IGAElementGetWorkVec(elemProj,&FlocProj);CHKERRQ(ierr);
		ierr = IGAElementGetValues(elemAlp,arrayAlProj,&AlProj);CHKERRQ(ierr);
		ierr = IGAElementGetValues(elemS,arraySProj,&SProj);CHKERRQ(ierr);
		ierr = IGAElementGetValues(elemVa,arrayVaProj,&VaProj);CHKERRQ(ierr);

		//FormSystem loop
		while (IGAElementNextFormSystem(elemProj,&wtfProj,&wtf2Proj)) 
		{
			//Quadrature loop
			ierr = IGAElementBeginPoint(elemProj,&pointProj);CHKERRQ(ierr);
			ierr = IGAElementBeginPoint(elemAlp,&pointAlp);CHKERRQ(ierr);
			ierr = IGAElementBeginPoint(elemS,&pointS);CHKERRQ(ierr);
			ierr = IGAElementBeginPoint(elemVa,&pointVa);CHKERRQ(ierr);

			while (IGAElementNextPoint(elemProj,pointProj)) 
			{
				if(pointProj->atboundary==0 && pointAlp->atboundary==0 && pointS->atboundary==0 && pointVa->atboundary==0)
				{
					IGAElementNextPoint(elemAlp,pointAlp);
					IGAElementNextPoint(elemS,pointS);
					IGAElementNextPoint(elemVa,pointVa);

					ierr = IGAPointGetWorkMat(pointProj,&KpointProj);CHKERRQ(ierr);
					ierr = IGAPointGetWorkVec(pointProj,&FpointProj);CHKERRQ(ierr);
					//proj(IGAPoint p,IGAPoint pAl,IGAPoint pS, IGAPoint pVa,PetscReal *K,PetscReal *F,PetscReal *UAl,PetscReal *US,PetscReal *UVa,void *ctx)
					ierr = proj(pointProj,pointAlp,pointS,pointVa,KpointProj,FpointProj,AlProj,SProj,VaProj,NULL);CHKERRQ(ierr);
					ierr = IGAPointAddMat(pointProj,KpointProj,KlocProj);CHKERRQ(ierr);
					ierr = IGAPointAddVec(pointProj,FpointProj,FlocProj);CHKERRQ(ierr);
				}
			}
			IGAElementNextPoint(elemAlp,pointAlp);
			IGAElementNextPoint(elemS,pointS);
			IGAElementNextPoint(elemVa,pointVa);

			ierr = IGAElementEndPoint(elemProj,&pointProj);CHKERRQ(ierr);
			ierr = IGAElementEndPoint(elemAlp,&pointAlp);CHKERRQ(ierr);
			ierr = IGAElementEndPoint(elemS,&pointS);CHKERRQ(ierr);
			ierr = IGAElementEndPoint(elemVa,&pointVa);CHKERRQ(ierr);
		}
		ierr = IGAElementAssembleMat(elemProj,KlocProj,KProj);CHKERRQ(ierr);
		ierr = IGAElementAssembleVec(elemProj,FlocProj,FProj);CHKERRQ(ierr);

	}
	IGANextElement(igaAl,elemAlp);
	IGANextElement(igaS,elemS);
	IGANextElement(igaVa,elemVa);

	ierr = IGAEndElement(igaProj,&elemProj);CHKERRQ(ierr);
	ierr = IGAEndElement(igaAl,&elemAlp);CHKERRQ(ierr);
	ierr = IGAEndElement(igaS,&elemS);CHKERRQ(ierr);
	ierr = IGAEndElement(igaVa,&elemVa);CHKERRQ(ierr);

	// Restore local vectors s0 and arrays
	ierr = IGARestoreLocalVecArray(igaAl,alp0,&localAlProj,&arrayAlProj);CHKERRQ(ierr);
	ierr = IGARestoreLocalVecArray(igaS,s0,&localSProj,&arraySProj);CHKERRQ(ierr);
	ierr = IGARestoreLocalVecArray(igaVa,Va0,&localVaProj,&arrayVaProj);CHKERRQ(ierr);

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
	char pathProj[512];
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
			ierr = IGASetRuleSize(igaZdot,dir,8);CHKERRQ(ierr);
		}
		ierr = IGASetUp(igaZdot);CHKERRQ(ierr);

		for (dir=0; dir<2; dir++) 
		{
			for (side=0; side<2; side++) 
			{
				//ierr = IGASetBoundaryValue(iga,dir,side,dof,0.0);CHKERRQ(ierr);					// Dirichlet boundary conditions
				ierr = IGASetBoundaryForm(igaZdot,dir,side,PETSC_TRUE);CHKERRQ(ierr);  				// Neumann boundary conditions
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
		char pathZdot[512];
	//

	//Things for u
		IGA igaU;
		ierr = IGACreate(PETSC_COMM_WORLD,&igaU);CHKERRQ(ierr);
		ierr = IGASetDim(igaU,2);CHKERRQ(ierr);													//Spatial dimension of the problem
		ierr = IGASetDof(igaU,2);CHKERRQ(ierr);													//Number of degrees of freedom, per node
		ierr = IGASetOrder(igaU,1);CHKERRQ(ierr);													//Number of spatial derivatives to calculate
		ierr = IGASetFromOptions(igaU);CHKERRQ(ierr);												//Note: The order (or degree) of the shape functions is given by the mesh!
		ierr = IGARead(igaU,"./geometry.dat");CHKERRQ(ierr);
		
		for (dir=0; dir<2; dir++)
		{
			ierr = IGASetRuleType(igaU,dir,IGA_RULE_LEGENDRE);CHKERRQ(ierr);
			ierr = IGASetRuleSize(igaU,dir,6);CHKERRQ(ierr);
		}
		ierr = IGASetUp(igaU);CHKERRQ(ierr);

		Mat KU;
		Vec U,FU;
		ierr = IGACreateMat(igaU,&KU);CHKERRQ(ierr);
		ierr = IGACreateVec(igaU,&U);CHKERRQ(ierr);
		ierr = IGACreateVec(igaU,&FU);CHKERRQ(ierr);

		IGAPoint		pointU;									//point
		IGAElement		elemU;									//element
		PetscReal		*KlocU,*FlocU;							//AA y BB
		PetscReal		*KpointU,*FpointU;						//KKK y FFF
		const PetscReal	*arrayChiU, *arrayZ0U;		//arrayU
		Vec				localChiU, localZ0U;			//localU
		PetscReal		*ChiU, *Z0U;						//U

		IGAFormSystem	wtfU0;
		void			*wtf2U0;

		KSP kspU;
		ierr = IGACreateKSP(igaU,&kspU);CHKERRQ(ierr);

		char nameU[512];
		char pathU[512];
	//

	//Things for Vs
		const PetscReal	*arrayUVs;						//arrayU
		Vec				localUVs;						//localU
		PetscReal		*UVs;							//U
	//

	//Things for Va
		const PetscReal	*arrayUVa;						//arrayU
		Vec				localUVa;						//localU
		PetscReal		*UVa;							//U
	//
//

for (i=1;i<800;i++)
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
	ierr = IGAReadVec(igaS,s0,pathS); CHKERRQ(ierr);
//

//Creation of types and systems for the Helmholtz decomposition of S, curl part for S based input
	PetscPrintf(PETSC_COMM_WORLD,"\nSystem for curl part of Helmholtz of S starting for iteration %d \n\n",i);
	T=time(NULL);
	tm=*localtime(&T);
	PetscPrintf(PETSC_COMM_WORLD,"Current time is %02d:%02d:%02d \n",tm.tm_hour,tm.tm_min,tm.tm_sec);
	
	ierr = VecZeroEntries(FchiS);CHKERRQ(ierr);							//This sets all elements of the vector to 0
	//Not necessary to zero out chiS0, it's fully overwritten every time. We also don't zero the matrix, it's reused in every iteration, as it doesn't change.

	//Get local vectors s0 and arrays
	ierr = IGAGetLocalVecArray(igaS,s0,&localS0chiS,&arrayS0chiS);CHKERRQ(ierr);

	//Element loop
	ierr = IGABeginElement(igachiS,&elemchiS);CHKERRQ(ierr);
	ierr = IGABeginElement(igaS,&elemS);CHKERRQ(ierr);

	while (IGANextElement(igachiS,elemchiS))
	{
		IGANextElement(igaS,elemS);

		ierr = IGAElementGetWorkVec(elemchiS,&FlocchiS);CHKERRQ(ierr);
		ierr = IGAElementGetValues(elemS,arrayS0chiS,&S0chiS);CHKERRQ(ierr);

		//FormSystem loop
		while (IGAElementNextFormSystem(elemchiS,&wtfchiS,&wtf2chiS)) 
		{
			//Quadrature loop
			ierr = IGAElementBeginPoint(elemchiS,&pointchiS);CHKERRQ(ierr);
			ierr = IGAElementBeginPoint(elemS,&pointS);CHKERRQ(ierr);
			
			while (IGAElementNextPoint(elemchiS,pointchiS)) 
			{
				if(pointchiS->atboundary==0 && pointS->atboundary==0)
				{
					IGAElementNextPoint(elemS,pointS);

					ierr = IGAPointGetWorkVec(pointchiS,&FpointchiS);CHKERRQ(ierr);
					ierr = curlChiSF(pointchiS,pointS,FpointchiS,S0chiS,NULL);CHKERRQ(ierr);
					ierr = IGAPointAddVec(pointchiS,FpointchiS,FlocchiS);CHKERRQ(ierr);
				}
			}
			while (pointS->index != -1)
			{
				IGAElementNextPoint(elemS,pointS);
			}
			ierr = IGAElementEndPoint(elemchiS,&pointchiS);CHKERRQ(ierr);
			ierr = IGAElementEndPoint(elemS,&pointS);CHKERRQ(ierr);
		}
		ierr = IGAElementFixSystem(elemchiS,KlocchiS,FlocchiS);CHKERRQ(ierr);
		ierr = IGAElementAssembleVec(elemchiS,FlocchiS,FchiS);CHKERRQ(ierr);

	}
	IGANextElement(igaS,elemS);
	ierr = IGAEndElement(igachiS,&elemchiS);CHKERRQ(ierr);
	ierr = IGAEndElement(igaS,&elemS);CHKERRQ(ierr);

	// Restore local vectors s0 and arrays
	ierr = IGARestoreLocalVecArray(igaS,s0,&localS0chiS,&arrayS0chiS);CHKERRQ(ierr);

	ierr = VecAssemblyBegin(FchiS);CHKERRQ(ierr);
	ierr = VecAssemblyEnd  (FchiS);CHKERRQ(ierr);

	ierr = KSPSolve(kspchiS,FchiS,chiS0);CHKERRQ(ierr);
	
	sprintf(namechiS,"%s%d%s","/ChiS-2d-",i,".dat");
	sprintf(pathchiS,"%s%s",direct,namechiS);
	ierr = IGAWriteVec(igachiS,chiS0,pathchiS);CHKERRQ(ierr);
//

//Creation of types and systems for the L2 projection of Alfa0-Sp:X
	PetscPrintf(PETSC_COMM_WORLD,"\nSystem for L2 projection for Alfa-Sp:X starting for iteration %d \n\n",i);
	T=time(NULL);
	tm=*localtime(&T);
	PetscPrintf(PETSC_COMM_WORLD,"Current time is %02d:%02d:%02d \n",tm.tm_hour,tm.tm_min,tm.tm_sec);

	ierr = MatZeroEntries(KAlp);CHKERRQ(ierr);						//This makes all non-zero elements of KchiS 0.0, while keeping the sparse structure of the matrix
	ierr = VecZeroEntries(FAlp);CHKERRQ(ierr);						//This sets all elements of the vector to 0
	//Not necessary to zero out alInput and alp0, they are fully overwritten every time.

	//Get local vectors s0 and arrays
	ierr = IGAGetLocalVecArray(igachiS,chiS0,&localS0Alp,&arrayS0Alp);CHKERRQ(ierr);

	//Element loop
	ierr = IGABeginElement(igaAl,&elemAlp);CHKERRQ(ierr);
	ierr = IGABeginElement(igachiS,&elemSp);CHKERRQ(ierr);

	while (IGANextElement(igaAl,elemAlp))
	{
		IGANextElement(igachiS,elemSp);
		ierr = IGAElementGetWorkMat(elemAlp,&KlocAlp);CHKERRQ(ierr);
		ierr = IGAElementGetWorkVec(elemAlp,&FlocAlp);CHKERRQ(ierr);
		ierr = IGAElementGetValues(elemSp,arrayS0Alp,&S0Alp);CHKERRQ(ierr);

		//FormSystem loop
		while (IGAElementNextFormSystem(elemAlp,&wtfAlp,&wtf2Alp)) 
		{
			//Quadrature loop
			ierr = IGAElementBeginPoint(elemAlp,&pointAlp);CHKERRQ(ierr);
			ierr = IGAElementBeginPoint(elemSp,&pointSp);CHKERRQ(ierr);
			while (IGAElementNextPoint(elemAlp,pointAlp)) 
			{
				if(pointAlp->atboundary==0 && pointSp->atboundary==0)
				{
					IGAElementNextPoint(elemSp,pointSp);
					ierr = IGAPointGetWorkMat(pointAlp,&KpointAlp);CHKERRQ(ierr);
					ierr = IGAPointGetWorkVec(pointAlp,&FpointAlp);CHKERRQ(ierr);
					ierr = L2ProjectionAlSp(pointAlp,pointSp,KpointAlp,FpointAlp,S0Alp,&user);CHKERRQ(ierr);
					ierr = IGAPointAddMat(pointAlp,KpointAlp,KlocAlp);CHKERRQ(ierr);
					ierr = IGAPointAddVec(pointAlp,FpointAlp,FlocAlp);CHKERRQ(ierr);
				}
			}
			IGAElementNextPoint(elemSp,pointSp);
			
			ierr = IGAElementEndPoint(elemAlp,&pointAlp);CHKERRQ(ierr);
			ierr = IGAElementEndPoint(elemSp,&pointSp);CHKERRQ(ierr);
		}
		//ierr = IGAElementFixSystem(elemAlp,KlocAlp,FlocAlp);CHKERRQ(ierr);					//This sets Dirichlet condition ¿?
		ierr = IGAElementAssembleMat(elemAlp,KlocAlp,KAlp);CHKERRQ(ierr);
		ierr = IGAElementAssembleVec(elemAlp,FlocAlp,FAlp);CHKERRQ(ierr);

	}
	IGANextElement(igachiS,elemSp);
	ierr = IGAEndElement(igaAl,&elemAlp);CHKERRQ(ierr);
	ierr = IGAEndElement(igachiS,&elemSp);CHKERRQ(ierr);

	// Restore local vectors s0 and arrays
	ierr = IGARestoreLocalVecArray(igachiS,chiS0,&localS0Alp,&arrayS0Alp);CHKERRQ(ierr);

	//Form system matrix and vector
	ierr = MatAssemblyBegin(KAlp,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyEnd  (KAlp,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	
	ierr = VecAssemblyBegin(FAlp);CHKERRQ(ierr);
	ierr = VecAssemblyEnd  (FAlp);CHKERRQ(ierr);

	ierr = KSPSetOperators(kspAlp,KAlp,KAlp);CHKERRQ(ierr);
	ierr = KSPSetFromOptions(kspAlp);CHKERRQ(ierr);
	ierr = KSPSetTolerances(kspAlp,1.0e-25,1.0e-40,PETSC_DEFAULT,PETSC_DEFAULT);CHKERRQ(ierr);
	ierr = KSPSolve(kspAlp,FAlp,alp0);CHKERRQ(ierr);

	ierr = IGAReadVec(igaAl,alInput,pathAlInput); CHKERRQ(ierr);

	//Computes y = alpha*x + y. VecAXPY(Vec y,PetscScalar alpha,Vec x)
	ierr = VecAXPY(alp0,1.0,alInput); CHKERRQ(ierr);

	sprintf(nameAlp,"%s%d%s","/Alp-2d-",i,".dat");
	sprintf(pathAlp,"%s%s",direct,nameAlp);
	ierr = IGAWriteVec(igaAl,alp0,pathAlp);CHKERRQ(ierr);
//

//Creation of types and systems for the Helmholtz decomposition of Up (or Ue), curl part
	PetscPrintf(PETSC_COMM_WORLD,"\nSystem for curl part of Helmholtz of Up starting for iteration %d \n\n",i);
	T=time(NULL);
	tm=*localtime(&T);
	PetscPrintf(PETSC_COMM_WORLD,"Current time is %02d:%02d:%02d \n",tm.tm_hour,tm.tm_min,tm.tm_sec);
	
	ierr = VecZeroEntries(FchiUp);CHKERRQ(ierr);						//This sets all elements of the vector to 0
	//Not necessary to zero out chiUp0, it's fully overwritten every time, we don't zero the matrix either, as we reuse it from the first part.
	//Reusing the matrix increases speed about tenfold

	//Get local vectors s0 and arrays
	ierr = IGAGetLocalVecArray(igaAl,alp0,&localAl0chiUp,&arrayAl0chiUp);CHKERRQ(ierr);

	//Element loop
	ierr = IGABeginElement(igachiUp,&elemchiUp);CHKERRQ(ierr);
	ierr = IGABeginElement(igaAl,&elemAlp);CHKERRQ(ierr);

	while (IGANextElement(igachiUp,elemchiUp))
	{
		IGANextElement(igaAl,elemAlp);
		ierr = IGAElementGetWorkVec(elemchiUp,&FlocchiUp);CHKERRQ(ierr);
		ierr = IGAElementGetValues(elemAlp,arrayAl0chiUp,&Al0chiUp);CHKERRQ(ierr);

		//FormSystem loop
		while (IGAElementNextFormSystem(elemchiUp,&wtfchiUp,&wtf2chiUp)) 
		{
			//Quadrature loop
			ierr = IGAElementBeginPoint(elemchiUp,&pointchiUp);CHKERRQ(ierr);
			ierr = IGAElementBeginPoint(elemAlp,&pointAlp);CHKERRQ(ierr);

			while (IGAElementNextPoint(elemchiUp,pointchiUp)) 
			{
				if(pointchiUp->atboundary==0 && pointAlp->atboundary==0)
				{
					IGAElementNextPoint(elemAlp,pointAlp);

					ierr = IGAPointGetWorkVec(pointchiUp,&FpointchiUp);CHKERRQ(ierr);
					ierr = curlChiUF(pointchiUp,pointAlp,FpointchiUp,Al0chiUp,NULL);CHKERRQ(ierr);
					ierr = IGAPointAddVec(pointchiUp,FpointchiUp,FlocchiUp);CHKERRQ(ierr);
				}
			}
			IGAElementNextPoint(elemAlp,pointAlp);

			ierr = IGAElementEndPoint(elemchiUp,&pointchiUp);CHKERRQ(ierr);
			ierr = IGAElementEndPoint(elemAlp,&pointAlp);CHKERRQ(ierr);
		}
		ierr = IGAElementFixSystem(elemchiUp,KlocchiUp,FlocchiUp);CHKERRQ(ierr);
		ierr = IGAElementAssembleVec(elemchiUp,FlocchiUp,FchiUp);CHKERRQ(ierr);

	}
	IGANextElement(igaAl,elemAlp);
	ierr = IGAEndElement(igachiUp,&elemchiUp);CHKERRQ(ierr);
	ierr = IGAEndElement(igaAl,&elemAlp);CHKERRQ(ierr);

	// Restore local vectors s0 and arrays
	ierr = IGARestoreLocalVecArray(igaAl,alp0,&localAl0chiUp,&arrayAl0chiUp);CHKERRQ(ierr);
	
	ierr = VecAssemblyBegin(FchiUp);CHKERRQ(ierr);
	ierr = VecAssemblyEnd  (FchiUp);CHKERRQ(ierr);

	ierr = KSPSolve(kspchiUp,FchiUp,chiUp0);CHKERRQ(ierr);

	sprintf(namechiUp,"%s%d%s","/ChiUp-2d-",i,".dat");
	sprintf(pathchiUp,"%s%s",direct,namechiUp);
	ierr = IGAWriteVec(igachiUp,chiUp0,pathchiUp);CHKERRQ(ierr);
//

//System for evolution of z0 (debug)
	PetscPrintf(PETSC_COMM_WORLD,"\nSystem for dot{z} starting for iteration %d \n\n",i);
	T=time(NULL);
	tm=*localtime(&T);
	PetscPrintf(PETSC_COMM_WORLD,"Current time is %02d:%02d:%02d \n",tm.tm_hour,tm.tm_min,tm.tm_sec);

	ierr = MatZeroEntries(KZdot);CHKERRQ(ierr);						//This makes all non-zero elements of KchiS 0.0, while keeping the sparse structure of the matrix
	ierr = VecZeroEntries(FZdot);CHKERRQ(ierr);						//This sets all elements of the vector to 0

	// Get local vectors Z0 and Chi0 and arrays
	ierr = IGAGetLocalVecArray(igaAl,alp0,&localAlZdot,&arrayAlZdot);CHKERRQ(ierr);
	ierr = IGAGetLocalVecArray(igaVa,Va0,&localVaZdot,&arrayVaZdot);CHKERRQ(ierr);
	ierr = IGAGetLocalVecArray(igaS,s0,&localSZdot,&arraySZdot);CHKERRQ(ierr);

	// Element loop
	ierr = IGABeginElement(igaZdot,&elemZdot);CHKERRQ(ierr);
	ierr = IGABeginElement(igaVa,&elemVa);CHKERRQ(ierr);
	ierr = IGABeginElement(igaAl,&elemAlp);CHKERRQ(ierr);
	ierr = IGABeginElement(igaS,&elemS);CHKERRQ(ierr);

	ierr = PetscPrintf(PETSC_COMM_WORLD,"Hola 1 \n");CHKERRQ(ierr);

	while (IGANextElement(igaZdot,elemZdot)) 
	{
		IGANextElement(igaVa,elemVa);
		IGANextElement(igaAl,elemAlp);
		IGANextElement(igaS,elemS);
		
		ierr = IGAElementGetWorkMat(elemZdot,&KlocZdot);CHKERRQ(ierr);
		ierr = IGAElementGetWorkVec(elemZdot,&FlocZdot);CHKERRQ(ierr);
		ierr = IGAElementGetValues(elemVa,arrayVaZdot,&VaZdot);CHKERRQ(ierr);
		ierr = IGAElementGetValues(elemAlp,arrayAlZdot,&AlZdot);CHKERRQ(ierr);
		ierr = IGAElementGetValues(elemS,arraySZdot,&SZdot);CHKERRQ(ierr);

		// FormSystem loop
		while (IGAElementNextFormSystem(elemZdot,&wtfZdot,&wtf2Zdot)) 
		{
		// Quadrature loop
			ierr = IGAElementBeginPoint(elemZdot,&pointZdot);CHKERRQ(ierr);
			ierr = IGAElementBeginPoint(elemVa,&pointVa);CHKERRQ(ierr);
			ierr = IGAElementBeginPoint(elemAlp,&pointAlp);CHKERRQ(ierr);
			ierr = IGAElementBeginPoint(elemS,&pointS);CHKERRQ(ierr);

			while (IGAElementNextPoint(elemZdot,pointZdot))
			{
				if(pointZdot->atboundary==1)
				{

				}
				if(pointZdot->atboundary==0 && pointAlp->atboundary==0 && pointVa->atboundary==0 && pointS->atboundary==0)
				{
					IGAElementNextPoint(elemAlp,pointAlp);
					IGAElementNextPoint(elemVa,pointVa);
					IGAElementNextPoint(elemS,pointS);

					ierr = IGAPointGetWorkMat(pointZdot,&KpointZdot);CHKERRQ(ierr);
					ierr = IGAPointGetWorkVec(pointZdot,&FpointZdot);CHKERRQ(ierr);
					//     ZdotSystem(IGAPoint p,IGAPoint pAl,IGAPoint pVa,IGAPoint pS,PetscReal *K,PetscReal *F,PetscReal *UAl,PetscReal *UVa,PetscReal *US, void *ctx)
					ierr = ZdotSystem(pointZdot,pointAlp,pointVa,pointS,KpointZdot,FpointZdot,AlZdot,VaZdot,SZdot,NULL);CHKERRQ(ierr);
					ierr = IGAPointAddMat(pointZdot,KpointZdot,KlocZdot);CHKERRQ(ierr);
					ierr = IGAPointAddVec(pointZdot,FpointZdot,FlocZdot);CHKERRQ(ierr);
				}
			}
			while (pointAlp->index != -1)
			{
				IGAElementNextPoint(elemAlp,pointAlp);
			}
			while (pointVa->index != -1)
			{
				IGAElementNextPoint(elemVa,pointVa);
			}
			while (pointS->index != -1)
			{
				IGAElementNextPoint(elemS,pointS);
			}

			ierr = IGAElementEndPoint(elemVa,&pointVa);CHKERRQ(ierr);
			ierr = IGAElementEndPoint(elemAlp,&pointAlp);CHKERRQ(ierr);
			ierr = IGAElementEndPoint(elemS,&pointS);CHKERRQ(ierr);
			ierr = IGAElementEndPoint(elemZdot,&pointZdot);CHKERRQ(ierr);
		}
		ierr = IGAElementAssembleMat(elemZdot,KlocZdot,KZdot);CHKERRQ(ierr);
		ierr = IGAElementAssembleVec(elemZdot,FlocZdot,FZdot);CHKERRQ(ierr);
	}
	IGANextElement(igaVa,elemVa);
	IGANextElement(igaAl,elemAlp);
	IGANextElement(igaS,elemS);

	ierr = PetscPrintf(PETSC_COMM_WORLD,"Hola 2 \n");CHKERRQ(ierr);

	ierr = IGAEndElement(igaZdot,&elemZdot);CHKERRQ(ierr);
	ierr = IGAEndElement(igaAl,&elemAlp);CHKERRQ(ierr);
	ierr = IGAEndElement(igaVa,&elemVa);CHKERRQ(ierr);
	ierr = IGAEndElement(igaS,&elemS);CHKERRQ(ierr);

	ierr = PetscPrintf(PETSC_COMM_WORLD,"Hola 3 \n");CHKERRQ(ierr);

	// Restore local vectors u, Z0, Chi0 and arrays
	ierr = IGARestoreLocalVecArray(igaVa,Va0,&localVaZdot,&arrayVaZdot);CHKERRQ(ierr);
	ierr = IGARestoreLocalVecArray(igaAl,alp0,&localAlZdot,&arrayAlZdot);CHKERRQ(ierr);
	ierr = IGARestoreLocalVecArray(igaS,s0,&localSZdot,&arraySZdot);CHKERRQ(ierr);
	
	ierr = MatAssemblyBegin(KZdot,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyEnd  (KZdot,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

	ierr = PetscPrintf(PETSC_COMM_WORLD,"Hola 4 \n");CHKERRQ(ierr);

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
	ierr = PetscPrintf(PETSC_COMM_WORLD,"Hola 5 \n");CHKERRQ(ierr);

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

	//Now z_(t+1)=z_(t)+dt*\dot{z}
	//VecAXPY(Vec y,PetscScalar alpha,Vec x)
	ierr = VecAXPY(Z0,dt,Zdot);CHKERRQ(ierr);					//This does y= alpha*x + y, in this case z0=z0 + dt*\dot{z}

	sprintf(nameZ0,"%s%d%s","/Z0-2d-",i,".dat");
	sprintf(pathZ0,"%s%s",direct,nameZ0);
	ierr = IGAWriteVec(igaZ0,Z0,pathZ0);CHKERRQ(ierr);			//Save updated Z to file
//

//System for u (debug)
	PetscPrintf(PETSC_COMM_WORLD,"\nSystem for u starting for iteration %d \n\n",i);
	T=time(NULL);
	tm=*localtime(&T);
	PetscPrintf(PETSC_COMM_WORLD,"Current time is %02d:%02d:%02d \n",tm.tm_hour,tm.tm_min,tm.tm_sec);
	
	ierr = MatZeroEntries(KU);CHKERRQ(ierr);						//This makes all non-zero elements of KU 0.0, while keeping the sparse structure of the matrix
	ierr = VecZeroEntries(FU);CHKERRQ(ierr);						//This sets all elements of the vector to 0

	fijaPunto=0;																		//Fix a single point (fijaPunto=1) or a side (chosen in blocks below)
	for (dir=0; dir<2; dir++) 
	{
		for (side=0; side<2; side++) 
		{
			//ierr = IGASetBoundaryValue(iga,dir,side,dof,0.0);CHKERRQ(ierr);					// Dirichlet boundary conditions
			ierr = IGASetBoundaryForm(igaU,dir,side,PETSC_TRUE);CHKERRQ(ierr);  				// Neumann boundary conditions
		}
	}

	if (fijaPunto==0)
	{
		//If we are not fixing a single point, set Dirichlet conditions here
		//ierr = IGASetBoundaryValue(iga,dir,side,dof,value);CHKERRQ(ierr);					// Dirichlet boundary conditions
		//ierr = IGASetBoundaryValue(igaU,0,0,0,0.0);CHKERRQ(ierr);	//Left side, 1st dof = 0
		//ierr = IGASetBoundaryValue(igaU,0,0,1,0.0);CHKERRQ(ierr);	//Left side, 2nd dof = 0

		//ierr = IGASetBoundaryValue(igaU,0,1,0,0.0);CHKERRQ(ierr);	//Right side, 1st dof=0
		//ierr = IGASetBoundaryValue(igaU,0,1,1,0.0);CHKERRQ(ierr);	//Right side, 2nd dof=0

		ierr = IGASetBoundaryValue(igaU,1,0,0,0.0);CHKERRQ(ierr);	//Bottom side, 1st dof=0
		ierr = IGASetBoundaryValue(igaU,1,0,1,0.0);CHKERRQ(ierr);	//Bottom side, 2nd dof=0
		
		//ierr = IGASetBoundaryValue(igaU,1,1,0,0.0);CHKERRQ(ierr);	//Top side, 1st dof=0
		//ierr = IGASetBoundaryValue(igaU,1,1,1,0.0);CHKERRQ(ierr);	//Top side, 2nd dof=0
	}

	// Get local vectors Chi0, Z0 and arrays
	ierr = IGAGetLocalVecArray(igachiUp,chiUp0,&localChiU,&arrayChiU);CHKERRQ(ierr);
	ierr = IGAGetLocalVecArray(igaZ0,Z0,&localZ0U,&arrayZ0U);CHKERRQ(ierr);

	// Element loop
	ierr = IGABeginElement(igaU,&elemU);CHKERRQ(ierr);
	ierr = IGABeginElement(igachiUp,&elemchiUp);CHKERRQ(ierr);
	ierr = IGABeginElement(igaZ0,&elemZ0);CHKERRQ(ierr);

	while (IGANextElement(igaU,elemU))
	{
		IGANextElement(igachiUp,elemchiUp);
		IGANextElement(igaZ0,elemZ0);

		ierr = IGAElementGetWorkMat(elemU,&KlocU);CHKERRQ(ierr);
		ierr = IGAElementGetWorkVec(elemU,&FlocU);CHKERRQ(ierr);
		ierr = IGAElementGetValues(elemchiUp,arrayChiU,&ChiU);CHKERRQ(ierr);
		ierr = IGAElementGetValues(elemZ0,arrayZ0U,&Z0U);CHKERRQ(ierr);

		// FormSystem loop
		while (IGAElementNextFormSystem(elemU,&wtfU0,&wtf2U0)) 
		{
		// Quadrature loop
			ierr = IGAElementBeginPoint(elemU,&pointU);CHKERRQ(ierr);
			ierr = IGAElementBeginPoint(elemchiUp,&pointchiUp);CHKERRQ(ierr);
			ierr = IGAElementBeginPoint(elemZ0,&pointZ0);CHKERRQ(ierr);

			while (IGAElementNextPoint(elemU,pointU))
			{
				if(pointU->atboundary==1)
				{
					ierr = IGAPointGetWorkMat(pointU,&KpointU);CHKERRQ(ierr);
					ierr = IGAPointGetWorkVec(pointU,&FpointU);CHKERRQ(ierr);
					//	   Usys(IGAPoint p, IGAPoint pChi, IGAPoint pZ, PetscReal *K, PetscReal *F, PetscReal *UChi, PetscReal *UZ, void *ctx)
					ierr = Usys(pointU,pointchiUp,pointZ0,KpointU,FpointU,ChiU,Z0U,NULL);CHKERRQ(ierr);
					ierr = IGAPointAddMat(pointU,KpointU,KlocU);CHKERRQ(ierr);
					ierr = IGAPointAddVec(pointU,FpointU,FlocU);CHKERRQ(ierr);
				}
				if(pointU->atboundary==0 && pointchiUp->atboundary==0 && pointZ0->atboundary==0)
				{
					IGAElementNextPoint(elemchiUp,pointchiUp);
					IGAElementNextPoint(elemZ0,pointZ0);

					ierr = IGAPointGetWorkMat(pointU,&KpointU);CHKERRQ(ierr);
					ierr = IGAPointGetWorkVec(pointU,&FpointU);CHKERRQ(ierr);
					//	   Usys(IGAPoint p, IGAPoint pChi, IGAPoint pZ, PetscReal *K, PetscReal *F, PetscReal *UChi, PetscReal *UZ, void *ctx)
					ierr = Usys(pointU,pointchiUp,pointZ0,KpointU,FpointU,ChiU,Z0U,NULL);CHKERRQ(ierr);
					ierr = IGAPointAddMat(pointU,KpointU,KlocU);CHKERRQ(ierr);
					ierr = IGAPointAddVec(pointU,FpointU,FlocU);CHKERRQ(ierr);
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

			ierr = IGAElementEndPoint(elemU,&pointU);CHKERRQ(ierr);
			ierr = IGAElementEndPoint(elemchiUp,&pointchiUp);CHKERRQ(ierr);
			ierr = IGAElementEndPoint(elemZ0,&pointZ0);CHKERRQ(ierr);
		}

		if(fijaPunto==0)
		{
			ierr = IGAElementFixSystem(elemU,KlocU,FlocU);CHKERRQ(ierr);					//This sets Dirichlet condition ¿? (Yes, this applies the conditions from IGASetBoundaryValue)
		}
		ierr = IGAElementAssembleMat(elemU,KlocU,KU);CHKERRQ(ierr);
		ierr = IGAElementAssembleVec(elemU,FlocU,FU);CHKERRQ(ierr);
	}
	IGANextElement(igachiUp,elemchiUp);
	IGANextElement(igaZ0,elemZ0);

	ierr = IGAEndElement(igaU,&elemU);CHKERRQ(ierr);
	ierr = IGAEndElement(igachiUp,&elemchiUp);CHKERRQ(ierr);
	ierr = IGAEndElement(igaZ0,&elemZ0);CHKERRQ(ierr);

	// Restore local vectors Chi0 and arrays
	ierr = IGARestoreLocalVecArray(igachiUp,chiUp0,&localChiU,&arrayChiU);CHKERRQ(ierr);
	ierr = IGARestoreLocalVecArray(igaZ0,Z0,&localZ0U,&arrayZ0U);CHKERRQ(ierr);

	ierr = MatAssemblyBegin(KU,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyEnd  (KU,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

	if (fijaPunto==1)
	{
		//Here we set values to the Matrix directly, to impose Dirichlet condition in a single point.
		ierr = MatGetSize(KU,&n,&m);CHKERRQ(ierr);
		ierr = MatSetOption(KU, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);CHKERRQ(ierr);

		ierr = PetscMalloc1(m,&cols);CHKERRQ(ierr);
		ierr = PetscMalloc1(m,&vals);CHKERRQ(ierr);

		for(int i=0;i<m;i++)
		{
			cols[i]=i;
			vals[i]=0.0;
		}

		rows=0;
		vals[rows]=1.0e6;
		ierr = MatSetValues(KU,1,&rows,m,cols,vals,INSERT_VALUES);CHKERRQ(ierr);
		ierr = MatSetValues(KU,n,cols,1,&rows,vals,INSERT_VALUES);CHKERRQ(ierr);
		vals[rows]=0.0;

		rows=1;
		vals[rows]=1.0e6;
		ierr = MatSetValues(KU,1,&rows,m,cols,vals,INSERT_VALUES);CHKERRQ(ierr);
		ierr = MatSetValues(KU,n,cols,1,&rows,vals,INSERT_VALUES);CHKERRQ(ierr);
		vals[rows]=0.0;

		rows=2*(nx+1)*(ny+1)-2;					//-2 is because last point is 2*(nx+1)*(ny+1)-1 (due to 0 index) and I want to fix the x-dof
		vals[rows]=1.0e6;
		ierr = MatSetValues(KU,1,&rows,m,cols,vals,INSERT_VALUES);CHKERRQ(ierr);
		ierr = MatSetValues(KU,n,cols,1,&rows,vals,INSERT_VALUES);CHKERRQ(ierr);
		vals[rows]=0.0;

		ierr = MatAssemblyBegin(KU,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
		ierr = MatAssemblyEnd  (KU,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

		//Here we set values to the vector directly (to impose Dirichlet condition in a single point)
		ierr = VecAssemblyBegin(FU);CHKERRQ(ierr);
		ierr = VecAssemblyEnd  (FU);CHKERRQ(ierr);

		rows=0;
		val=0.0;
		ierr = VecSetValue(FU,rows,val,INSERT_VALUES);CHKERRQ(ierr);

		ierr = VecAssemblyBegin(FU);CHKERRQ(ierr);
		ierr = VecAssemblyEnd  (FU);CHKERRQ(ierr);

		rows=1;
		val=0.0;
		ierr = VecSetValue(FU,rows,val,INSERT_VALUES);CHKERRQ(ierr);

		ierr = VecAssemblyBegin(FU);CHKERRQ(ierr);
		ierr = VecAssemblyEnd  (FU);CHKERRQ(ierr);
		
		rows=2*(nx+1)*(ny+1)-2;					//-2 is because last point is 2*(nx+1)*(ny+1)-1 (due to 0 index) and I want to fix the x-dof
		val=0.0;
		ierr = VecSetValue(FU,rows,val,INSERT_VALUES);CHKERRQ(ierr);
	}

	ierr = VecAssemblyBegin(FU);CHKERRQ(ierr);
	ierr = VecAssemblyEnd  (FU);CHKERRQ(ierr);


	ierr = KSPSetOperators(kspU,KU,KU);CHKERRQ(ierr);
	PC pcU;
	ierr = KSPGetPC(kspU,&pcU); CHKERRQ(ierr);
	ierr = PCSetType(pcU,PCLU); CHKERRQ(ierr);
	ierr = PCFactorSetMatSolverType(pcU,MATSOLVERMUMPS); CHKERRQ(ierr);
	//ierr = KSPSetFromOptions(kspU);CHKERRQ(ierr);
	//ierr = KSPSetTolerances(kspU,1e-24,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);CHKERRQ(ierr);
	ierr = KSPSolve(kspU,FU,U);CHKERRQ(ierr);

	sprintf(nameU,"%s%d%s","/U-2d-",i,".dat");
	sprintf(pathU,"%s%s",direct,nameU);
	ierr = IGAWriteVec(igaU,U,pathU);CHKERRQ(ierr);
//

//System for evolution of S  (debug)
	//Now reuses K from the initial case, just recalculates RHS
	PetscPrintf(PETSC_COMM_WORLD,"\nSystem for S_dot starting for iteration %d \n\n",i);
	T=time(NULL);
	tm=*localtime(&T);
	PetscPrintf(PETSC_COMM_WORLD,"Current time is %02d:%02d:%02d \n",tm.tm_hour,tm.tm_min,tm.tm_sec);

	ierr = VecZeroEntries(FSdot);CHKERRQ(ierr);						//This sets all elements of the vector to 0
	//Not necessary to zero out Sdot, it's fully overwritten every time. We reuse the matrix KSdot from the initial step, it does not change.

	// Get local vectors Z0 and Chi0 and arrays
	ierr = IGAGetLocalVecArray(igaVs,VsSmooth,&localVsSdot,&arrayVsSdot);CHKERRQ(ierr);
	ierr = IGAGetLocalVecArray(igaS,s0,&localS0Sdot,&arrayS0Sdot);CHKERRQ(ierr);

	// Element loop
	ierr = IGABeginElement(igaVs,&elemVs);CHKERRQ(ierr);
	ierr = IGABeginElement(igaS,&elemS);CHKERRQ(ierr);
	ierr = IGABeginElement(igaSdot,&elemSdot);CHKERRQ(ierr);

	while (IGANextElement(igaSdot,elemSdot)) 
	{
		IGANextElement(igaVs,elemVs);
		IGANextElement(igaS,elemS);

		ierr = IGAElementGetWorkVec(elemSdot,&FlocSdot);CHKERRQ(ierr);
		ierr = IGAElementGetValues(elemVs,arrayVsSdot,&VsSdot);CHKERRQ(ierr);
		ierr = IGAElementGetValues(elemS,arrayS0Sdot,&S0Sdot);CHKERRQ(ierr);

		// FormSystem loop
		while (IGAElementNextFormSystem(elemSdot,&wtfSdot,&wtf2Sdot)) 
		{
		// Quadrature loop
			ierr = IGAElementBeginPoint(elemSdot,&pointSdot);CHKERRQ(ierr);
			ierr = IGAElementBeginPoint(elemVs,&pointVs);CHKERRQ(ierr);
			ierr = IGAElementBeginPoint(elemS,&pointS);CHKERRQ(ierr);

			while (IGAElementNextPoint(elemSdot,pointSdot))
			{
				if(pointSdot->atboundary==1)
				{
					//Code boundary condition, will be required for general case
				}
				if(pointSdot->atboundary==0 && pointVs->atboundary==0 && pointS->atboundary==0)
				{
					IGAElementNextPoint(elemVs,pointVs);
					IGAElementNextPoint(elemS,pointS);

					ierr = IGAPointGetWorkVec(pointSdot,&FpointSdot);CHKERRQ(ierr);
					ierr = SdotFuncF(pointSdot,pointVs,pointS,FpointSdot,VsSdot,S0Sdot,&user);CHKERRQ(ierr);
					ierr = IGAPointAddVec(pointSdot,FpointSdot,FlocSdot);CHKERRQ(ierr);
				}
			}
			while (pointVs->index != -1)
			{
				IGAElementNextPoint(elemVs,pointVs);
			}
			
			while (pointS->index != -1)
			{
				IGAElementNextPoint(elemS,pointS);
			}
			ierr = IGAElementEndPoint(elemSdot,&pointSdot);CHKERRQ(ierr);
			ierr = IGAElementEndPoint(elemVs,&pointVs);CHKERRQ(ierr);
			ierr = IGAElementEndPoint(elemS,&pointS);CHKERRQ(ierr);
		}
		//ierr = IGAElementFixSystemF(elemSdot,FlocSdot);CHKERRQ(ierr);					//This sets Dirichlet condition ¿? (Yes, this applies the conditions from IGASetBoundaryValue)
		ierr = IGAElementAssembleVec(elemSdot,FlocSdot,FSdot);CHKERRQ(ierr);
	}
	IGANextElement(igaVs,elemVs);
	IGANextElement(igaS,elemS);

	ierr = IGAEndElement(igaSdot,&elemSdot);CHKERRQ(ierr);
	ierr = IGAEndElement(igaVs,&elemVs);CHKERRQ(ierr);
	ierr = IGAEndElement(igaS,&elemS);CHKERRQ(ierr);

	// Restore local vectors u, Z0, Chi0 and arrays
	ierr = IGARestoreLocalVecArray(igaVs,VsSmooth,&localVsSdot,&arrayVsSdot);CHKERRQ(ierr);
	ierr = IGARestoreLocalVecArray(igaS,s0,&localS0Sdot,&arrayS0Sdot);CHKERRQ(ierr);

	ierr = VecAssemblyBegin(FSdot);CHKERRQ(ierr);
	ierr = VecAssemblyEnd  (FSdot);CHKERRQ(ierr);

	ierr = KSPSolve(kspSdot,FSdot,Sdot);CHKERRQ(ierr);

	//VecChop(Vec v, PetscReal tol) Sets anything with an absolute value less than the tolerance to 0
	//ierr = VecChop(Sdot,1e-8);CHKERRQ(ierr);

	//Kill everything in S[0],S[2],S[4].S[6], just because I haven't been able to stop it from coming out from the solution
		ierr = VecGetSize(Sdot,&n);CHKERRQ(ierr);

		for(int i=0;i<n;i=i+8)
		{
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
	ierr = IGAGetLocalVecArray(igachiUp,chiUp0,&localChi0Vs,&arrayChi0Vs);CHKERRQ(ierr);
	ierr = IGAGetLocalVecArray(igaZ0,Z0,&localZ0Vs,&arrayZ0Vs);CHKERRQ(ierr);
	ierr = IGAGetLocalVecArray(igaS,s0,&localS0Vs,&arrayS0Vs);CHKERRQ(ierr);
	ierr = IGAGetLocalVecArray(igaU,U,&localUVs,&arrayUVs);CHKERRQ(ierr);

	// Element loop
	ierr = IGABeginElement(igaVs,&elemVs);CHKERRQ(ierr);
	ierr = IGABeginElement(igaZ0,&elemZ0);CHKERRQ(ierr);
	ierr = IGABeginElement(igachiUp,&elemchiUp);CHKERRQ(ierr);
	ierr = IGABeginElement(igaS,&elemS);CHKERRQ(ierr);
	ierr = IGABeginElement(igaU,&elemU);CHKERRQ(ierr);

	while (IGANextElement(igaVs,elemVs)) 
	{
		IGANextElement(igaZ0,elemZ0);
		IGANextElement(igachiUp,elemchiUp);
		IGANextElement(igaS,elemS);
		IGANextElement(igaU,elemU);

		//ierr = IGAElementGetWorkMat(elemVs,&KlocVs);CHKERRQ(ierr);
		ierr = IGAElementGetWorkVec(elemVs,&FlocVs);CHKERRQ(ierr);
		ierr = IGAElementGetValues(elemZ0,arrayZ0Vs,&Z0Vs);CHKERRQ(ierr);
		ierr = IGAElementGetValues(elemchiUp,arrayChi0Vs,&Chi0Vs);CHKERRQ(ierr);
		ierr = IGAElementGetValues(elemS,arrayS0Vs,&S0Vs);CHKERRQ(ierr);
		ierr = IGAElementGetValues(elemU,arrayUVs,&UVs);CHKERRQ(ierr);

		// FormSystem loop
		while (IGAElementNextFormSystem(elemVs,&wtfVs,&wtf2Vs)) 
		{
		// Quadrature loop
			ierr = IGAElementBeginPoint(elemVs,&pointVs);CHKERRQ(ierr);
			ierr = IGAElementBeginPoint(elemZ0,&pointZ0);CHKERRQ(ierr);
			ierr = IGAElementBeginPoint(elemchiUp,&pointchiUp);CHKERRQ(ierr);
			ierr = IGAElementBeginPoint(elemS,&pointS);CHKERRQ(ierr);
			ierr = IGAElementBeginPoint(elemU,&pointU);CHKERRQ(ierr);

			while (IGAElementNextPoint(elemVs,pointVs))
			{
				if(pointVs->atboundary==1)
				{

				}
				if(pointVs->atboundary==0 && pointZ0->atboundary==0 && pointchiUp->atboundary==0 && pointS->atboundary==0)
				{
					IGAElementNextPoint(elemZ0,pointZ0);
					IGAElementNextPoint(elemchiUp,pointchiUp);
					IGAElementNextPoint(elemS,pointS);
					IGAElementNextPoint(elemU,pointU);

					ierr = IGAPointGetWorkVec(pointVs,&FpointVs);CHKERRQ(ierr);
					//VSF(IGAPoint p,IGAPoint pChi,IGAPoint pZu,IGAPoint pu,IGAPoint pS,PetscReal *F, PetscReal *UChi,PetscReal *UZu,PetscReal *Uu,PetscReal *US,void *ctx)
					ierr = VSF(pointVs,pointchiUp,pointZ0,pointU,pointS,FpointVs,Chi0Vs,Z0Vs,UVs,S0Vs,NULL);CHKERRQ(ierr);
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
			while (pointU->index != -1)
			{
				IGAElementNextPoint(elemU,pointU);
			}
			ierr = IGAElementEndPoint(elemVs,&pointVs);CHKERRQ(ierr);
			ierr = IGAElementEndPoint(elemZ0,&pointZ0);CHKERRQ(ierr);
			ierr = IGAElementEndPoint(elemchiUp,&pointchiUp);CHKERRQ(ierr);
			ierr = IGAElementEndPoint(elemS,&pointS);CHKERRQ(ierr);
			ierr = IGAElementEndPoint(elemU,&pointU);CHKERRQ(ierr);
		}

		//ierr = IGAElementFixSystem(elemStress,KlocStress,FlocStress);CHKERRQ(ierr);					//This sets Dirichlet condition ¿? (Yes, this applies the conditions from IGASetBoundaryValue)
		ierr = IGAElementAssembleVec(elemVs,FlocVs,FVs);CHKERRQ(ierr);
	}
	IGANextElement(igaZ0,elemZ0);
	IGANextElement(igachiUp,elemchiUp);
	IGANextElement(igaS,elemS);
	IGANextElement(igaU,elemU);

	ierr = IGAEndElement(igaVs,&elemVs);CHKERRQ(ierr);
	ierr = IGAEndElement(igaZ0,&elemZ0);CHKERRQ(ierr);
	ierr = IGAEndElement(igachiUp,&elemchiUp);CHKERRQ(ierr);
	ierr = IGAEndElement(igaS,&elemS);CHKERRQ(ierr);
	ierr = IGAEndElement(igaU,&elemU);CHKERRQ(ierr);

	// Restore local vectors u, Z0, Chi0 and arrays
	ierr = IGARestoreLocalVecArray(igaZ0,Z0,&localZ0Vs,&arrayZ0Vs);CHKERRQ(ierr);
	ierr = IGARestoreLocalVecArray(igachiUp,chiUp0,&localChi0Vs,&arrayChi0Vs);CHKERRQ(ierr);
	ierr = IGARestoreLocalVecArray(igaS,s0,&localS0Vs,&arrayS0Vs);CHKERRQ(ierr);
	ierr = IGARestoreLocalVecArray(igaU,U,&localUVs,&arrayUVs);CHKERRQ(ierr);

	ierr = VecAssemblyBegin(FVs);CHKERRQ(ierr);
	ierr = VecAssemblyEnd  (FVs);CHKERRQ(ierr);

	//ierr = KSPSetOperators(kspVs,KVs,KVs);CHKERRQ(ierr);
	//PC pcVs;
	//ierr = KSPGetPC(kspVs,&pcVs); CHKERRQ(ierr);
	//ierr = PCSetType(pcVs,PCLU); CHKERRQ(ierr);
	//ierr = PCFactorSetMatSolverType(pcVs,MATSOLVERMUMPS); CHKERRQ(ierr);
	//ierr = KSPSetFromOptions(kspVs);CHKERRQ(ierr);
	//ierr = KSPSetTolerances(kspVs,1.0e-24,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);CHKERRQ(ierr);
	ierr = KSPSolve(kspVs,FVs,Vs0);CHKERRQ(ierr);

	//VecChop(Vec v, PetscReal tol) Sets anything with an absolute value less than the tolerance to 0
	ierr = VecChop(Vs0,1e-11);CHKERRQ(ierr);

	//char nameVs[512];//="/Vs-2d-0.dat";
	sprintf(nameVs,"%s%d%s","/Vs-2d-",i,".dat");
	//char pathVs[512];
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
	ierr = VecZeroEntries(IntVec);CHKERRQ(ierr);						//This sets all elements of the vector to 0
	//First integrate velocity for both defects

	// Get local vectors V0 and InputAl and arrays
	ierr = IGAGetLocalVecArray(igaVs,Vs0,&localVs0,&arrayVs0);CHKERRQ(ierr);
	ierr = IGAGetLocalVecArray(igaS,s0,&localSVs,&arraySVs);CHKERRQ(ierr);

	// Element loop
	ierr = IGABeginElement(igaVs,&elemVs);CHKERRQ(ierr);
	ierr = IGABeginElement(igaS,&elemS);CHKERRQ(ierr);

	while (IGANextElement(igaVs,elemVs))
	{
		IGANextElement(igaS,elemS);

		ierr = IGAElementGetWorkVec(elemVs,&locInt1a);CHKERRQ(ierr);
		ierr = IGAElementGetWorkVec(elemVs,&locInt2a);CHKERRQ(ierr);
		ierr = IGAElementGetWorkVec(elemVs,&locInt);CHKERRQ(ierr);
		ierr = IGAElementGetValues(elemVs,arrayVs0,&Vs0Values);CHKERRQ(ierr);
		ierr = IGAElementGetValues(elemS,arraySVs,&SVs);CHKERRQ(ierr);

		// FormSystem loop
		while (IGAElementNextFormSystem(elemVs,&wtfVsInt,&wtf2VsInt))
		{
		// Quadrature loop
			ierr = IGAElementBeginPoint(elemVs,&pointVs);CHKERRQ(ierr);
			ierr = IGAElementBeginPoint(elemS,&pointS);CHKERRQ(ierr);
			
			while (IGAElementNextPoint(elemVs,pointVs))
			{
				if(pointVs->atboundary==1)
				{

				}
				if(pointVs->atboundary==0 && pointS->atboundary==0)
				{
					IGAElementNextPoint(elemS,pointS);

					ierr = IGAPointGetWorkVec(pointVs,&PointInt1a);CHKERRQ(ierr);
					ierr = IGAPointGetWorkVec(pointVs,&PointInt2a);CHKERRQ(ierr);
					ierr = IGAPointGetWorkVec(pointVs,&PointInt);CHKERRQ(ierr);
					//Int_Xi_Vs(IGAPoint pV,IGAPoint pS,PetscReal *FInt1a,PetscReal *FInt2a,PetscReal *FInt,PetscReal *VS,PetscReal *US,void *ctx)
					ierr = Int_Xi_Vs(pointVs,pointS,PointInt1a,PointInt2a,PointInt,Vs0Values,SVs,NULL);CHKERRQ(ierr);
					ierr = IGAPointAddVec(pointVs,PointInt1a,locInt1a);CHKERRQ(ierr);
					ierr = IGAPointAddVec(pointVs,PointInt2a,locInt2a);CHKERRQ(ierr);
					ierr = IGAPointAddVec(pointVs,PointInt,locInt);CHKERRQ(ierr);
				}
			}
			while (pointS->index != -1)
			{
				IGAElementNextPoint(elemS,pointS);
			}
			ierr = IGAElementEndPoint(elemVs,&pointVs);CHKERRQ(ierr);
			ierr = IGAElementEndPoint(elemS,&pointS);CHKERRQ(ierr);
		}
		ierr = IGAElementAssembleVec(elemVs,locInt1a,Int1sVec);CHKERRQ(ierr);
		ierr = IGAElementAssembleVec(elemVs,locInt2a,Int2sVec);CHKERRQ(ierr);
		ierr = IGAElementAssembleVec(elemVs,locInt,IntVec);CHKERRQ(ierr);
	}
	IGANextElement(igaS,elemS);

	ierr = IGAEndElement(igaVs,&elemVs);CHKERRQ(ierr);
	ierr = IGAEndElement(igaS,&elemS);CHKERRQ(ierr);

	ierr = IGARestoreLocalVecArray(igaS,s0,&localSVs,&arraySVs);CHKERRQ(ierr);
	ierr = IGARestoreLocalVecArray(igaVs,Vs0,&localVs0,&arrayVs0);CHKERRQ(ierr);

	ierr = VecAssemblyBegin(Int1sVec);CHKERRQ(ierr);
	ierr = VecAssemblyEnd  (Int1sVec);CHKERRQ(ierr);
	ierr = VecAssemblyBegin(Int2sVec);CHKERRQ(ierr);
	ierr = VecAssemblyEnd  (Int2sVec);CHKERRQ(ierr);

	ierr = VecAssemblyBegin(IntVec);CHKERRQ(ierr);
	ierr = VecAssemblyEnd  (IntVec);CHKERRQ(ierr);

	//Here add all values at Gauss points
	ierr = VecSum(Int1sVec,&Int1a);CHKERRQ(ierr);
	ierr = VecSum(Int2sVec,&Int2a);CHKERRQ(ierr);
	ierr = VecSum(IntVec,&IntS);CHKERRQ(ierr);

	ierr = PetscPrintf(PETSC_COMM_WORLD,"Int1a=%f, Int2a=%f, IntS=%f \n",Int1a,Int2a,IntS);CHKERRQ(ierr);

	//Now put those velocities on each defect location
	ierr = VecZeroEntries(FVsSmooth);CHKERRQ(ierr);						//This sets all elements of the vector to 0
	ierr = VecZeroEntries(VsSmooth);CHKERRQ(ierr);						//This sets all elements of the vector to 0

	//PetscReal		*VsSmoothPoint,*locSmooth;										//AA y BB

	//PetscReal Smin,Smax,absSmax;

	ierr = VecMin(s0,NULL,&Smin);
	ierr = VecMax(s0,NULL,&Smax);
	absSmax=fmax(fabs(Smin),fabs(Smax));

	ierr = PetscPrintf(PETSC_COMM_WORLD,"Smin=%f, Smax=%f, absSmax=%f \n",Smin,Smax,absSmax);CHKERRQ(ierr);

	// Get local vectors V0 and InputAl and arrays
	ierr = IGAGetLocalVecArray(igaVs,Vs0,&localVs0,&arrayVs0);CHKERRQ(ierr);
	ierr = IGAGetLocalVecArray(igaS,s0,&localSVs,&arraySVs);CHKERRQ(ierr);

	// Element loop
	ierr = IGABeginElement(igaVs,&elemVs);CHKERRQ(ierr);
	ierr = IGABeginElement(igaS,&elemS);CHKERRQ(ierr);

	while (IGANextElement(igaVs,elemVs))
	{
		IGANextElement(igaS,elemS);

		ierr = IGAElementGetWorkVec(elemVs,&locSmooth);CHKERRQ(ierr);
		ierr = IGAElementGetValues(elemVs,arrayVs0,&Vs0Values);CHKERRQ(ierr);
		ierr = IGAElementGetValues(elemS,arraySVs,&SVs);CHKERRQ(ierr);

		// FormSystem loop
		while (IGAElementNextFormSystem(elemVs,&wtfVsInt,&wtf2VsInt))
		{
		// Quadrature loop
			ierr = IGAElementBeginPoint(elemVs,&pointVs);CHKERRQ(ierr);
			ierr = IGAElementBeginPoint(elemS,&pointS);CHKERRQ(ierr);
			
			while (IGAElementNextPoint(elemVs,pointVs))
			{
				if(pointVs->atboundary==1)
				{

				}
				if(pointVs->atboundary==0 && pointS->atboundary==0)
				{
					IGAElementNextPoint(elemS,pointS);
					ierr = IGAPointGetWorkVec(pointVs,&VsSmoothPoint);CHKERRQ(ierr);
					//ProjVS(IGAPoint pV,IGAPoint pS,PetscReal *FS1,PetscReal *US,PetscReal Vs1,PetscReal Vs2,PetscReal IntS,void *ctx)
					ierr = ProjVS(pointVs,pointS,VsSmoothPoint,SVs,Int1a,Int2a,IntS,absSmax,NULL);CHKERRQ(ierr);
					ierr = IGAPointAddVec(pointVs,VsSmoothPoint,locSmooth);CHKERRQ(ierr);
				}
			}
			while (pointS->index != -1)
			{
				IGAElementNextPoint(elemS,pointS);
			}
			ierr = IGAElementEndPoint(elemVs,&pointVs);CHKERRQ(ierr);
			ierr = IGAElementEndPoint(elemS,&pointS);CHKERRQ(ierr);
		}
		ierr = IGAElementAssembleVec(elemVs,locSmooth,FVsSmooth);CHKERRQ(ierr);
	}
	IGANextElement(igaS,elemS);

	ierr = IGAEndElement(igaVs,&elemVs);CHKERRQ(ierr);
	ierr = IGAEndElement(igaS,&elemS);CHKERRQ(ierr);

	ierr = IGARestoreLocalVecArray(igaS,s0,&localSVs,&arraySVs);CHKERRQ(ierr);
	ierr = IGARestoreLocalVecArray(igaVs,Vs0,&localVs0,&arrayVs0);CHKERRQ(ierr);

	//Here we have a vector with an L2 Projection of the integrated velocity.
	ierr = VecAssemblyBegin(FVsSmooth);CHKERRQ(ierr);
	ierr = VecAssemblyEnd  (FVsSmooth);CHKERRQ(ierr);

	//Here we solved the L2 projection system
	ierr = KSPSolve(kspVs,FVsSmooth,VsSmooth);CHKERRQ(ierr);

	//char nameVsS1[512];//="/VsSmooth-2d-0.dat";
	sprintf(nameVsS1,"%s%d%s","/VsSmooth-2d-",i,".dat");
	//char pathVsS1[512];
	sprintf(pathVsS1,"%s%s",direct,nameVsS1);
	ierr = IGAWriteVec(igaVs,VsSmooth,pathVsS1);CHKERRQ(ierr);

	//ierr = KSPDestroy(&kspVs);CHKERRQ(ierr);
	//ierr = MatDestroy(&KVs);CHKERRQ(ierr);
//

//System for L2 projection of V^{alpha} (debug)
	PetscPrintf(PETSC_COMM_WORLD,"\nSystem for V-alpha starting for iteration %d \n\n",i);
	T=time(NULL);
	tm=*localtime(&T);
	PetscPrintf(PETSC_COMM_WORLD,"Current time is %02d:%02d:%02d \n",tm.tm_hour,tm.tm_min,tm.tm_sec);
	
	ierr = VecZeroEntries(FVa);CHKERRQ(ierr);						//This sets all elements of the vector to 0

	// Get local vectors Z0 and Chi0 and arrays
	ierr = IGAGetLocalVecArray(igachiUp,chiUp0,&localChi0Va,&arrayChi0Va);CHKERRQ(ierr);
	ierr = IGAGetLocalVecArray(igaZ0,Z0,&localZ0Va,&arrayZ0Va);CHKERRQ(ierr);
	ierr = IGAGetLocalVecArray(igaAl,alInput,&localAl0Va,&arrayAl0Va);CHKERRQ(ierr);
	ierr = IGAGetLocalVecArray(igaS,s0,&localSVa,&arraySVa);CHKERRQ(ierr);
	ierr = IGAGetLocalVecArray(igaU,U,&localUVa,&arrayUVa);CHKERRQ(ierr);

	// Element loop
	ierr = IGABeginElement(igaVa,&elemVa);CHKERRQ(ierr);
	ierr = IGABeginElement(igaZ0,&elemZ0);CHKERRQ(ierr);
	ierr = IGABeginElement(igachiUp,&elemchiUp);CHKERRQ(ierr);
	ierr = IGABeginElement(igaAl,&elemAlp);CHKERRQ(ierr);
	ierr = IGABeginElement(igaS,&elemS);CHKERRQ(ierr);
	ierr = IGABeginElement(igaU,&elemU);CHKERRQ(ierr);

	while (IGANextElement(igaVa,elemVa)) 
	{
		IGANextElement(igaZ0,elemZ0);
		IGANextElement(igachiUp,elemchiUp);
		IGANextElement(igaAl,elemAlp);
		IGANextElement(igaS,elemS);
		IGANextElement(igaU,elemU);

		ierr = IGAElementGetWorkMat(elemVa,&KlocVa);CHKERRQ(ierr);
		ierr = IGAElementGetWorkVec(elemVa,&FlocVa);CHKERRQ(ierr);
		ierr = IGAElementGetValues(elemZ0,arrayZ0Va,&Z0Va);CHKERRQ(ierr);
		ierr = IGAElementGetValues(elemchiUp,arrayChi0Va,&Chi0Va);CHKERRQ(ierr);
		ierr = IGAElementGetValues(elemAlp,arrayAl0Va,&Al0Va);CHKERRQ(ierr);
		ierr = IGAElementGetValues(elemS,arraySVa,&SVa);CHKERRQ(ierr);
		ierr = IGAElementGetValues(elemU,arrayUVa,&UVa);CHKERRQ(ierr);

		// FormSystem loop
		while (IGAElementNextFormSystem(elemVa,&wtfVa,&wtf2Va)) 
		{
		// Quadrature loop
			ierr = IGAElementBeginPoint(elemVa,&pointVa);CHKERRQ(ierr);
			ierr = IGAElementBeginPoint(elemZ0,&pointZ0);CHKERRQ(ierr);
			ierr = IGAElementBeginPoint(elemchiUp,&pointchiUp);CHKERRQ(ierr);
			ierr = IGAElementBeginPoint(elemAlp,&pointAlp);CHKERRQ(ierr);
			ierr = IGAElementBeginPoint(elemS,&pointS);CHKERRQ(ierr);
			ierr = IGAElementBeginPoint(elemU,&pointU);CHKERRQ(ierr);

			while (IGAElementNextPoint(elemVa,pointVa))
			{
				if(pointVa->atboundary==1)
				{

				}
				if(pointVa->atboundary==0 && pointZ0->atboundary==0 && pointchiUp->atboundary==0 && pointAlp->atboundary==0 && pointS->atboundary==0 && pointU->atboundary==0)
				{
					IGAElementNextPoint(elemZ0,pointZ0);
					IGAElementNextPoint(elemchiUp,pointchiUp);
					IGAElementNextPoint(elemAlp,pointAlp);
					IGAElementNextPoint(elemS,pointS);
					IGAElementNextPoint(elemU,pointU);

					ierr = IGAPointGetWorkVec(pointVa,&FpointVa);CHKERRQ(ierr);
					//ValphaF(IGAPoint p,IGAPoint pChi,IGAPoint pu,IGAPoint pZu,IGAPoint pAl,IGAPoint pS,PetscReal *F,PetscReal *Chi,PetscReal *Zu,PetscReal *US,PetscReal *Uu,PetscReal *UAl,void *ctx)
					ierr = ValphaF(pointVa,pointchiUp,pointU,pointZ0,pointAlp,pointS,FpointVa,Chi0Va,Z0Va,SVa,UVa,Al0Va,NULL);CHKERRQ(ierr);
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
			while (pointS->index != -1)
			{
				IGAElementNextPoint(elemS,pointS);
			}
			while (pointU->index != -1)
			{
				IGAElementNextPoint(elemU,pointU);
			}
			ierr = IGAElementEndPoint(elemVa,&pointVa);CHKERRQ(ierr);
			ierr = IGAElementEndPoint(elemZ0,&pointZ0);CHKERRQ(ierr);
			ierr = IGAElementEndPoint(elemchiUp,&pointchiUp);CHKERRQ(ierr);
			ierr = IGAElementEndPoint(elemAlp,&pointAlp);CHKERRQ(ierr);
			ierr = IGAElementEndPoint(elemS,&pointS);CHKERRQ(ierr);
			ierr = IGAElementEndPoint(elemU,&pointU);CHKERRQ(ierr);
		}

		//ierr = IGAElementFixSystem(elemStress,KlocStress,FlocStress);CHKERRQ(ierr);					//This sets Dirichlet condition ¿? (Yes, this applies the conditions from IGASetBoundaryValue)
		ierr = IGAElementAssembleVec(elemVa,FlocVa,FVa);CHKERRQ(ierr);
	}
	IGANextElement(igaZ0,elemZ0);
	IGANextElement(igachiUp,elemchiUp);
	IGANextElement(igaAl,elemAlp);
	IGANextElement(igaS,elemS);
	IGANextElement(igaU,elemU);

	ierr = IGAEndElement(igaVa,&elemVa);CHKERRQ(ierr);
	ierr = IGAEndElement(igaZ0,&elemZ0);CHKERRQ(ierr);
	ierr = IGAEndElement(igachiUp,&elemchiUp);CHKERRQ(ierr);
	ierr = IGAEndElement(igaAl,&elemAlp);CHKERRQ(ierr);
	ierr = IGAEndElement(igaU,&elemU);CHKERRQ(ierr);
	ierr = IGAEndElement(igaS,&elemS);CHKERRQ(ierr);

	// Restore local vectors u, Z0, Chi0 and arrays
	ierr = IGARestoreLocalVecArray(igaZ0,Z0,&localZ0Va,&arrayZ0Va);CHKERRQ(ierr);
	ierr = IGARestoreLocalVecArray(igachiUp,chiUp0,&localChi0Va,&arrayChi0Va);CHKERRQ(ierr);
	ierr = IGARestoreLocalVecArray(igaAl,alInput,&localAl0Va,&arrayAl0Va);CHKERRQ(ierr);
	ierr = IGARestoreLocalVecArray(igaS,s0,&localSVa,&arraySVa);CHKERRQ(ierr);
	ierr = IGARestoreLocalVecArray(igaU,U,&localUVa,&arrayUVa);CHKERRQ(ierr);

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
	ierr = IGAGetLocalVecArray(igaZ0,Z0,&localZ0Gradz,&arrayZ0Gradz);CHKERRQ(ierr);

	//Element loop
	ierr = IGABeginElement(igaGradz,&elemGradz);CHKERRQ(ierr);
	ierr = IGABeginElement(igaZ0,&elemZ0);CHKERRQ(ierr);

	while (IGANextElement(igaGradz,elemGradz))
	{
		IGANextElement(igaZ0,elemZ0);

		ierr = IGAElementGetWorkVec(elemGradz,&FlocGradz);CHKERRQ(ierr);
		ierr = IGAElementGetValues(elemZ0,arrayZ0Gradz,&Z0Gradz);CHKERRQ(ierr);

		//FormSystem loop
		while (IGAElementNextFormSystem(elemGradz,&wtfGradz,&wtf2Gradz)) 
		{
			//Quadrature loop
			ierr = IGAElementBeginPoint(elemGradz,&pointGradz);CHKERRQ(ierr);
			ierr = IGAElementBeginPoint(elemZ0,&pointZ0);CHKERRQ(ierr);

			while (IGAElementNextPoint(elemGradz,pointGradz))
			{
				if(pointGradz->atboundary==0 && pointZ0->atboundary==0)
				{
					IGAElementNextPoint(elemZ0,pointZ0);

					ierr = IGAPointGetWorkVec(pointGradz,&FpointGradz);CHKERRQ(ierr);
					ierr = gradzF(pointGradz,pointZ0,FpointGradz,Z0Gradz,NULL);CHKERRQ(ierr);
					ierr = IGAPointAddVec(pointGradz,FpointGradz,FlocGradz);CHKERRQ(ierr);
				}
			}
			IGAElementNextPoint(elemZ0,pointZ0);

			ierr = IGAElementEndPoint(elemGradz,&pointGradz);CHKERRQ(ierr);
			while (pointZ0->index != -1)
			{
				IGAElementNextPoint(elemZ0,pointZ0);
			}
			ierr = IGAElementEndPoint(elemZ0,&pointZ0);CHKERRQ(ierr);
		}
		ierr = IGAElementAssembleVec(elemGradz,FlocGradz,FGradz0);CHKERRQ(ierr);

	}
	IGANextElement(igaZ0,elemZ0);
	ierr = IGAEndElement(igaGradz,&elemGradz);CHKERRQ(ierr);
	ierr = IGAEndElement(igaZ0,&elemZ0);CHKERRQ(ierr);

	// Restore local vectors s0 and arrays
	ierr = IGARestoreLocalVecArray(igaZ0,Z0,&localZ0Gradz,&arrayZ0Gradz);CHKERRQ(ierr);

	//Form system vector
	ierr = VecAssemblyBegin(FGradz0);CHKERRQ(ierr);
	ierr = VecAssemblyEnd  (FGradz0);CHKERRQ(ierr);

	ierr = KSPSolve(kspGradz,FGradz0,gradZ0);CHKERRQ(ierr);

	sprintf(nameGradz0,"%s%d%s","/Gradz-2d-",i,".dat");
	sprintf(pathGradz0,"%s%s",direct,nameGradz0);
	ierr = IGAWriteVec(igaGradz,gradZ0,pathGradz0);CHKERRQ(ierr);
//

//System for general L2 projection (debug, but not really important, just for looking at fields)
	PetscPrintf(PETSC_COMM_WORLD,"\nSystem for general L2 projection starting for iteration %d \n\n",i);
	PetscPrintf(PETSC_COMM_WORLD,"Currently projecting alpha X V^a + SV^s \n");
	T=time(NULL);
	tm=*localtime(&T);
	PetscPrintf(PETSC_COMM_WORLD,"Current time is %02d:%02d:%02d \n",tm.tm_hour,tm.tm_min,tm.tm_sec);
	
	ierr = MatZeroEntries(KProj);CHKERRQ(ierr);						//This makes all non-zero elements of KchiS 0.0, while keeping the sparse structure of the matrix
	ierr = VecZeroEntries(FProj);CHKERRQ(ierr);	

	//Get local vectors and arrays
	ierr = IGAGetLocalVecArray(igaAl,alp0,&localAlProj,&arrayAlProj);CHKERRQ(ierr);
	ierr = IGAGetLocalVecArray(igaS,s0,&localSProj,&arraySProj);CHKERRQ(ierr);
	ierr = IGAGetLocalVecArray(igaVa,Va0,&localVaProj,&arrayVaProj);CHKERRQ(ierr);

	//Element loop
	ierr = IGABeginElement(igaProj,&elemProj);CHKERRQ(ierr);
	ierr = IGABeginElement(igaAl,&elemAlp);CHKERRQ(ierr);
	ierr = IGABeginElement(igaS,&elemS);CHKERRQ(ierr);
	ierr = IGABeginElement(igaVa,&elemVa);CHKERRQ(ierr);

	while (IGANextElement(igaProj,elemProj))
	{
		IGANextElement(igaAl,elemAlp);
		IGANextElement(igaS,elemS);
		IGANextElement(igaVa,elemVa);

		ierr = IGAElementGetWorkMat(elemProj,&KlocProj);CHKERRQ(ierr);
		ierr = IGAElementGetWorkVec(elemProj,&FlocProj);CHKERRQ(ierr);
		ierr = IGAElementGetValues(elemAlp,arrayAlProj,&AlProj);CHKERRQ(ierr);
		ierr = IGAElementGetValues(elemS,arraySProj,&SProj);CHKERRQ(ierr);
		ierr = IGAElementGetValues(elemVa,arrayVaProj,&VaProj);CHKERRQ(ierr);

		//FormSystem loop
		while (IGAElementNextFormSystem(elemProj,&wtfProj,&wtf2Proj)) 
		{
			//Quadrature loop
			ierr = IGAElementBeginPoint(elemProj,&pointProj);CHKERRQ(ierr);
			ierr = IGAElementBeginPoint(elemAlp,&pointAlp);CHKERRQ(ierr);
			ierr = IGAElementBeginPoint(elemS,&pointS);CHKERRQ(ierr);
			ierr = IGAElementBeginPoint(elemVa,&pointVa);CHKERRQ(ierr);

			while (IGAElementNextPoint(elemProj,pointProj)) 
			{
				if(pointProj->atboundary==0 && pointAlp->atboundary==0 && pointS->atboundary==0 && pointVa->atboundary==0)
				{
					IGAElementNextPoint(elemAlp,pointAlp);
					IGAElementNextPoint(elemS,pointS);
					IGAElementNextPoint(elemVa,pointVa);

					ierr = IGAPointGetWorkMat(pointProj,&KpointProj);CHKERRQ(ierr);
					ierr = IGAPointGetWorkVec(pointProj,&FpointProj);CHKERRQ(ierr);
					//proj(IGAPoint p,IGAPoint pAl,IGAPoint pS, IGAPoint pVa,PetscReal *K,PetscReal *F,PetscReal *UAl,PetscReal *US,PetscReal *UVa,void *ctx)
					ierr = proj(pointProj,pointAlp,pointS,pointVa,KpointProj,FpointProj,AlProj,SProj,VaProj,NULL);CHKERRQ(ierr);
					ierr = IGAPointAddMat(pointProj,KpointProj,KlocProj);CHKERRQ(ierr);
					ierr = IGAPointAddVec(pointProj,FpointProj,FlocProj);CHKERRQ(ierr);
				}
			}
			IGAElementNextPoint(elemAlp,pointAlp);
			IGAElementNextPoint(elemS,pointS);
			IGAElementNextPoint(elemVa,pointVa);

			ierr = IGAElementEndPoint(elemProj,&pointProj);CHKERRQ(ierr);
			ierr = IGAElementEndPoint(elemAlp,&pointAlp);CHKERRQ(ierr);
			ierr = IGAElementEndPoint(elemS,&pointS);CHKERRQ(ierr);
			ierr = IGAElementEndPoint(elemVa,&pointVa);CHKERRQ(ierr);
		}
		ierr = IGAElementAssembleMat(elemProj,KlocProj,KProj);CHKERRQ(ierr);
		ierr = IGAElementAssembleVec(elemProj,FlocProj,FProj);CHKERRQ(ierr);

	}
	IGANextElement(igaAl,elemAlp);
	IGANextElement(igaS,elemS);
	IGANextElement(igaVa,elemVa);

	ierr = IGAEndElement(igaProj,&elemProj);CHKERRQ(ierr);
	ierr = IGAEndElement(igaAl,&elemAlp);CHKERRQ(ierr);
	ierr = IGAEndElement(igaS,&elemS);CHKERRQ(ierr);
	ierr = IGAEndElement(igaVa,&elemVa);CHKERRQ(ierr);

	// Restore local vectors s0 and arrays
	ierr = IGARestoreLocalVecArray(igaAl,alp0,&localAlProj,&arrayAlProj);CHKERRQ(ierr);
	ierr = IGARestoreLocalVecArray(igaS,s0,&localSProj,&arraySProj);CHKERRQ(ierr);
	ierr = IGARestoreLocalVecArray(igaVa,Va0,&localVaProj,&arrayVaProj);CHKERRQ(ierr);

	//Form system matrix and vector
	ierr = MatAssemblyBegin(KProj,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyEnd  (KProj,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	
	ierr = VecAssemblyBegin(FProj);CHKERRQ(ierr);
	ierr = VecAssemblyEnd  (FProj);CHKERRQ(ierr);

	ierr = KSPSetOperators(kspProj,KProj,KProj);CHKERRQ(ierr);
	//ierr = KSPSetFromOptions(kspProj);CHKERRQ(ierr);
	//ierr = KSPSetTolerances(kspProj,1.0e-16,1.0e-30,PETSC_DEFAULT,PETSC_DEFAULT);CHKERRQ(ierr);
	ierr = KSPSolve(kspProj,FProj,proj0);CHKERRQ(ierr);

	//char nameProj[512];
	sprintf(nameProj,"%s%d%s","/Proj-2d-",i,".dat");
	//char pathProj[512];
	sprintf(pathProj,"%s%s",direct,nameProj);
	ierr = IGAWriteVec(igaProj,proj0,pathProj);CHKERRQ(ierr);
//
}


/*
//Destroy all objects not needed anymore (Better to do it here in case different codes call the same IGA, move if memory is a problem)
	ierr = IGADestroy(&igaAl); CHKERRQ(ierr);
	ierr = IGADestroy(&igaZ0); CHKERRQ(ierr);
	ierr = IGADestroy(&igachiUp); CHKERRQ(ierr);
	ierr = KSPDestroy(&kspZdot);CHKERRQ(ierr);
	ierr = MatDestroy(&KZdot);CHKERRQ(ierr);
	ierr = VecDestroy(&FZdot);CHKERRQ(ierr);
//
*/

ierr = PetscFinalize();CHKERRQ(ierr);

return 0;
}