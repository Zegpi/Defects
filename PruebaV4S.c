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

		//PetscReal S[8];																		//Create array to recieve Alfa
		PetscReal dS[8][2];																	//Create array to recieve dAlfa
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
		//Recieving Z but not using, check if it's neccesary

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
//

//System for z(0)
	#undef  __FUNCT__
	#define __FUNCT__ "Z0sys"
	PetscErrorCode Z0sys(IGAPoint p, IGAPoint pChi, IGAPoint pS, PetscReal *K, PetscReal *F, PetscReal *UChi, PetscReal *US, void *ctx)
	{
		PetscReal C[3][3][3][3]={0};														//Initialization of elastic tensor
		const PetscReal *N0,(*N1)[2],(*N2)[2][2];
		IGAPointGetShapeFuns(p,0,(const PetscReal**)&N0);									//Value of the shape functions
		IGAPointGetShapeFuns(p,1,(const PetscReal**)&N1);									//Derivatives of the shape functions
		IGAPointGetShapeFuns(p,2,(const PetscReal**)&N2);									//Second derivatives of the shape funcions
		//After this command Na_xx=N2[a][0][0], Na_yy=N2[a][1][1], Na_xy=N2[a][0][1], Na_yx[a][1][0] (these last two are equal) (remember a is the index of the shape function)
		PetscInt a,b,i,j,k,l,u,w,nen=p->nen, dof=p->dof;

		PetscReal x[2];																		//Vector of reals, size equal to problem's dimension
		IGAPointFormGeomMap(p,x);															//Fills x with the coordinates of p, Gauss's point

		//E and nu should come from AppCtx in the future
		//const PetscReal E=1.0; //2100000.0*9.81*10000.0;					//[Pa]
		//const PetscReal nu=0.33;
		//const PetscReal lambda=(E*nu)/((1.0+nu)*(1.0-2.0*nu));
		//const PetscReal mu=E/(2.0*(1.0+nu));
		//const PetscReal eps=0.0*E/1000.0;							//Choose later based on whatever Amit says :)

		//Change for G=1
		const PetscReal nu=0.33;
		const PetscReal mu=1.0;
		const PetscReal lambda=2.0*mu*nu/(1.0-2.0*nu);
		const PetscReal eps=mu/100.0;								//Choose later based on whatever Amit says :)

		PetscReal Chi0[4];																	//Assign chi to a vector
		PetscReal dChi0[4][2];																//Same for partial derivatives
		IGAPointFormValue(pChi,UChi,&Chi0[0]);
		IGAPointFormGrad (pChi,UChi,&dChi0[0][0]);

		PetscReal fullChi[3][3]={0};
		fullChi[0][0]=Chi0[0]; 	fullChi[0][1]=Chi0[1];
		fullChi[1][0]=Chi0[2]; 	fullChi[1][1]=Chi0[3];

		PetscReal full_dChi[3][3][3]={0};
		full_dChi[0][0][0]=dChi0[0][0]; full_dChi[0][0][1]=dChi0[0][1];
		full_dChi[0][1][0]=dChi0[1][0]; full_dChi[0][1][1]=dChi0[1][1];
		full_dChi[1][0][0]=dChi0[2][0]; full_dChi[1][0][1]=dChi0[2][1];
		full_dChi[1][1][0]=dChi0[3][0]; full_dChi[1][1][1]=dChi0[3][1];

		PetscReal S0[8];																	//Assign S to a vector
		IGAPointFormValue(pS,US,&S0[0]);

		PetscReal fullS[3][3][3]={0};
		fullS[0][0][0]=S0[0]; fullS[0][0][1]=S0[1];											//Expand S to full 3rd order form, only non-zero elements
		fullS[0][1][0]=S0[2]; fullS[0][1][1]=S0[3];
		fullS[1][0][0]=S0[4]; fullS[1][0][1]=S0[5];
		fullS[1][1][0]=S0[6]; fullS[1][1][1]=S0[7];

		PetscReal (*Keq)[dof][nen][dof] = (typeof(Keq)) K;
		PetscReal (*Feq)[dof] = (PetscReal (*)[dof])F;

		//////////////Delete this parte later, loads should come from appCtx or be 0
		PetscReal f[3]={0.0, 0.0, 0.0};			//Distributed load in body
		//These loads are for the manufactured solutions
		//PetscReal alpha=3.0;
		//PetscReal Lx=20.0; PetscReal Lx2=Lx*Lx; PetscReal Lx3=Lx*Lx*Lx; PetscReal Lx4=Lx*Lx*Lx*Lx;
		//PetscReal Ly=20.0; PetscReal Ly2=Ly*Ly; PetscReal Ly3=Ly*Ly*Ly; PetscReal Ly4=Ly*Ly*Ly*Ly;
		//f[0]=16.0*alpha/(425.0*Lx4*Ly4)*
		//	(425.0*Lx4*Ly2 -5100.0*Lx4*x[1]*x[1] +51.0*Lx4 +1675.0*Lx2*Ly4 -3400.0*Lx2*Ly2*x[0]*x[0] -13400.0*Lx2*Ly2*x[1]*x[1] +85.0*Lx2*Ly2 +40800.0*Lx2*x[0]*x[0]*x[1]*x[1]
		//    -408.0*Lx2*x[0]*x[0] +26800.0*Lx2*x[1]*x[1]*x[1]*x[1] -1020.0*Lx2*x[1]*x[1] -20100.0*Ly4*x[0]*x[0] +102.0*Ly4 +6800.0*Ly2*x[0]*x[0]*x[0]*x[0] +160800*Ly2*x[0]*x[0]*x[1]*x[1] -1020.0*Ly2*x[0]*x[0]
		//    -408.0*Ly2*x[0]*x[1] -816.0*Ly2*x[1]*x[1] -81600.0*x[0]*x[0]*x[0]*x[0]*x[1]*x[1] +816.0*x[0]*x[0]*x[0]*x[0] -321600.0*x[0]*x[0]*x[1]*x[1]*x[1]*x[1] +12240.0*x[0]*x[0]*x[1]*x[1] 
		//    +1632.0*x[0]*x[1]*x[1]*x[1] +1632.0*x[1]*x[1]*x[1]*x[1]);
		//f[1]=16.0*alpha/(425.0*Lx4*Ly4)*
		//	(-20000.0*Lx2*Ly2*x[0]*x[1] +17.0*Lx2*Ly2 +80000.0*Lx2*x[0]*x[1]*x[1]*x[1] -408.0*Lx2*x[0]*x[1] -204.0*Lx2*x[1]*x[1] +80000.0*Ly2*x[0]*x[0]*x[0]*x[1] -204.0*Ly2*x[0]*x[0] -816.0*Ly2*x[0]*x[1]
		//	-320000.0*x[0]*x[0]*x[0]*x[1]*x[1]*x[1] +1632.0*x[0]*x[0]*x[0]*x[1] +2448.0*x[0]*x[0]*x[1]*x[1] +3264.0*x[0]*x[1]*x[1]*x[1]);
		//These loads are for the manufactured solutions when S=grad[0 -y^2;y^2 0] with u(0)=(x+Lx/2)*(Lx/2-x)*(y+Ly/2)*(Ly/2-y) *16*a/(Lx^2*Ly^2), u(1)=0.
		//PetscReal ad,bd;
		//ad=20.0/101.0;
		//bd=0.01;
		//f[0]=eps*(((tanh((ad+x[1])/bd))/(bd*bd*cosh((ad+x[1])/bd)*cosh((ad+x[1])/bd)))+((tanh((ad-x[1])/bd))/(bd*bd*cosh((ad-x[1])/bd)*cosh((ad-x[1])/bd))));
		f[0]=0.0;
		f[1]=0.0;
		f[2]=0.0;
		PetscReal g[3]={0.0, 0.0, 0.0};			//Boundary load (applied wherever is defined in the p->atboundary block)
		PetscReal M[3]={0.0, 0.0, 0.0};			//Distributed moment in body
		PetscReal h[3]={0.0, 0.0, 0.0};			//Boundary moment (applied wherever is defined in the p->atboundary block)
		/////////////////////////////////////

		//Definition of alternating tensor
		const PetscReal e[3][3][3]=
		{
			{{0.0,0.0,0.0},{0.0,0.0,1.0},{0.0,-1.0,0.0}},
			{{0.0,0.0,-1.0},{0.0,0.0,0.0},{1.0,0.0,0.0}},
			{{0.0,1.0,0.0},{-1.0,0.0,0.0},{0.0,0.0,0.0}}
		};

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

		PetscReal n[3]={0};
		PetscReal v[3]={0};
		PetscReal dv[3][3]={0};
		PetscReal dz[3][3]={0};
		PetscReal d2z[3][3][3]={0};
		PetscReal d2v[3][3][3]={0};

		if (p->atboundary)
		{
			PetscReal Sborde[3][3]={0};

			//PetscReal Omega=tan(5.0/180.0*ConstPi);
			//PetscReal rho=sqrt(x[0]*x[0]+x[1]*x[1]);
			//const PetscReal burgers[2]={1.0,0.0};

			//Stress for single disclination
			//Sborde[0][0]=mu*Omega/(2.0*ConstPi*(1.0-nu))*(log(rho)+(x[1]*x[1])/(rho*rho)+nu/(1.0-2.0*nu));
			//Sborde[0][1]=-mu*Omega/(2.0*ConstPi*(1.0-nu))*x[0]*x[1]/(rho*rho);
			//Sborde[1][0]=-mu*Omega/(2.0*ConstPi*(1.0-nu))*x[0]*x[1]/(rho*rho);
			//Sborde[1][1]=mu*Omega/(2.0*ConstPi*(1.0-nu))*(log(rho)+(x[0]*x[0])/(rho*rho)+nu/(1.0-2.0*nu));
			//Sborde[0][2]=0.0; Sborde[1][2]=0.0; Sborde[2][0]=0.0; Sborde[2][1]=0.0; Sborde[2][2]=0.0;

			//Stress for single dislocation
			//Sborde[0][0]=-mu*burgers/(ConstPi*(1.0-nu))/2.0*x[1]*(x[1]*x[1]+3*x[0]*x[0])/((x[0]*x[0]+x[1]*x[1])*(x[0]*x[0]+x[1]*x[1]));
			//Sborde[0][1]=-mu*burgers/(ConstPi*(1.0-nu))/2.0*x[0]*(x[1]*x[1]-x[0]*x[0])  /((x[0]*x[0]+x[1]*x[1])*(x[0]*x[0]+x[1]*x[1]));
			//Sborde[1][0]=-mu*burgers/(ConstPi*(1.0-nu))/2.0*x[0]*(x[1]*x[1]-x[0]*x[0])  /((x[0]*x[0]+x[1]*x[1])*(x[0]*x[0]+x[1]*x[1]));
			//Sborde[1][1]=-mu*burgers/(ConstPi*(1.0-nu))/2.0*x[1]*(x[1]*x[1]-x[0]*x[0])  /((x[0]*x[0]+x[1]*x[1])*(x[0]*x[0]+x[1]*x[1]));
			//Sborde[0][2]=0.0; Sborde[1][2]=0.0; Sborde[2][0]=0.0; Sborde[2][1]=0.0; Sborde[2][2]=0.0;

			//When b in y axis
			//Sborde[0][0]=-mu*burgers/(ConstPi*(1.0-nu))/2.0*-x[0]*(x[0]*x[0]+3.0*x[1]*x[1])/((x[0]*x[0]+x[1]*x[1])*(x[0]*x[0]+x[1]*x[1]));
			//Sborde[0][1]=-mu*burgers/(ConstPi*(1.0-nu))/2.0*x[1]*(x[0]*x[0]-x[1]*x[1])/((x[0]*x[0]+x[1]*x[1])*(x[0]*x[0]+x[1]*x[1]));
			//Sborde[1][0]=-mu*burgers/(ConstPi*(1.0-nu))/2.0*x[1]*(x[0]*x[0]-x[1]*x[1])/((x[0]*x[0]+x[1]*x[1])*(x[0]*x[0]+x[1]*x[1]));
			//Sborde[1][1]=-mu*burgers/(ConstPi*(1.0-nu))/2.0*-x[0]*(x[0]*x[0]-x[1]*x[1])/((x[0]*x[0]+x[1]*x[1])*(x[0]*x[0]+x[1]*x[1]));
			//Sborde[0][2]=0.0; Sborde[1][2]=0.0; Sborde[2][0]=0.0; Sborde[2][1]=0.0; Sborde[2][2]=0.0;

			//General dislocation
			//Sborde[0][0]=-mu*burgers[0]/(ConstPi*(1.0-nu))/2.0*x[1]*(x[1]*x[1]+3*x[0]*x[0])/((x[0]*x[0]+x[1]*x[1])*(x[0]*x[0]+x[1]*x[1]))
			// 			 -mu*burgers[1]/(ConstPi*(1.0-nu))/2.0*-x[0]*(x[0]*x[0]+3.0*x[1]*x[1])/((x[0]*x[0]+x[1]*x[1])*(x[0]*x[0]+x[1]*x[1]));
			//Sborde[0][1]=-mu*burgers[0]/(ConstPi*(1.0-nu))/2.0*x[0]*(x[1]*x[1]-x[0]*x[0])/((x[0]*x[0]+x[1]*x[1])*(x[0]*x[0]+x[1]*x[1]))
			// 			 -mu*burgers[1]/(ConstPi*(1.0-nu))/2.0*x[1]*(x[0]*x[0]-x[1]*x[1])/((x[0]*x[0]+x[1]*x[1])*(x[0]*x[0]+x[1]*x[1]));
			//Sborde[1][0]=-mu*burgers[0]/(ConstPi*(1.0-nu))/2.0*x[0]*(x[1]*x[1]-x[0]*x[0])/((x[0]*x[0]+x[1]*x[1])*(x[0]*x[0]+x[1]*x[1]))
			// 			 -mu*burgers[1]/(ConstPi*(1.0-nu))/2.0*x[1]*(x[0]*x[0]-x[1]*x[1])/((x[0]*x[0]+x[1]*x[1])*(x[0]*x[0]+x[1]*x[1]));
			//Sborde[1][1]=-mu*burgers[0]/(ConstPi*(1.0-nu))/2.0*x[1]*(x[1]*x[1]-x[0]*x[0])/((x[0]*x[0]+x[1]*x[1])*(x[0]*x[0]+x[1]*x[1]))
			// 			 -mu*burgers[1]/(ConstPi*(1.0-nu))/2.0*-x[0]*(x[0]*x[0]-x[1]*x[1])/((x[0]*x[0]+x[1]*x[1])*(x[0]*x[0]+x[1]*x[1]));
			//Sborde[0][2]=0.0; Sborde[1][2]=0.0; Sborde[2][0]=0.0; Sborde[2][1]=0.0; Sborde[2][2]=0.0;

			//No stress in boundary
			Sborde[0][0]=0.0;
			Sborde[0][1]=-1.0;
			Sborde[1][0]=0.0;
			Sborde[1][1]=0.0;
			Sborde[0][2]=0.0; Sborde[1][2]=0.0; Sborde[2][0]=0.0; Sborde[2][1]=0.0; Sborde[2][2]=0.0;

			PetscInt dir  = p->boundary_id / 2;
			PetscInt side = p->boundary_id % 2;

			for (a=0; a<nen; a++)
			{
				PetscReal Na   = N0[a];
				PetscReal Na_x=N1[a][0];		PetscReal Na_y=N1[a][1];

				for (i=0; i<dof; i++)
				{
					if(i==0)
					{
						v[0]=Na; v[1]=0.0; v[2]=0.0;
						dv[0][0]=Na_x; dv[0][1]=Na_y; dv[0][2]=0.0;
						dv[1][0]=0.0;  dv[1][1]=0.0;  dv[1][2]=0.0;
						dv[2][0]=0.0;  dv[2][1]=0.0;  dv[2][2]=0.0;
					}
					if(i==1)
					{
						v[0]=0.0; v[1]=Na; v[2]=0.0;
						dv[0][0]=0.0;  dv[0][1]=0.0;  dv[0][2]=0.0;
						dv[1][0]=Na_x; dv[1][1]=Na_y; dv[1][2]=0.0;
						dv[2][0]=0.0;  dv[2][1]=0.0;  dv[2][2]=0.0;
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
								//Feq[a][i]+=g[j]*v[j];
							}

							/*
							for (r=0; r<3; r++)
							{
								for (s=0; s<3; s++)
								{
									for (m=0; m<3; m++)
									{
										Feq[a][i]+=h[r]*0.5*e[r][s][m]*dv[m][s];
									}
								}
							}
							*/
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
				PetscReal Na_xx =N2[a][0][0];		PetscReal Na_yy =N2[a][1][1];
				PetscReal Na_xy =N2[a][0][1];		PetscReal Na_yx =N2[a][1][0];

				for (b=0; b<nen; b++) 
				{
					PetscReal Nb_x = N1[b][0];		PetscReal Nb_y = N1[b][1];
					PetscReal Nb_xx =N2[b][0][0];	PetscReal Nb_yy =N2[b][1][1];
					PetscReal Nb_xy =N2[b][0][1];	PetscReal Nb_yx =N2[b][1][0];

					for (i=0; i<dof; i++)
					{
						if (i==0)
						{
							dv[0][0]=Na_x; dv[0][1]=Na_y; dv[0][2]=0.0;
							dv[1][0]=0.0;  dv[1][1]=0.0;  dv[1][2]=0.0;
							dv[2][0]=0.0;  dv[2][1]=0.0;  dv[2][2]=0.0;

							d2v[0][0][0]=Na_xx; d2v[0][0][1]=Na_xy; d2v[0][0][2]=0.0; d2v[0][1][0]=Na_yx; d2v[0][1][1]=Na_yy; d2v[0][1][2]=0.0; d2v[0][2][0]=0.0; d2v[0][2][1]=0.0; d2v[0][2][2]=0.0;
							d2v[1][0][0]=0.00;  d2v[1][0][1]=0.00;  d2v[1][0][2]=0.0; d2v[1][1][0]=0.0;   d2v[1][1][1]=0.0;   d2v[1][1][2]=0.0; d2v[1][2][0]=0.0; d2v[1][2][1]=0.0; d2v[1][2][2]=0.0;
							d2v[2][0][0]=0.00;  d2v[2][0][1]=0.00;  d2v[2][0][2]=0.0; d2v[2][1][0]=0.0;   d2v[2][1][1]=0.0;   d2v[2][1][2]=0.0; d2v[2][2][0]=0.0; d2v[2][2][1]=0.0; d2v[2][2][2]=0.0;
							
						}
						else if(i==1)
						{
							dv[0][0]=0.0;  dv[0][1]=0.0;  dv[0][2]=0.0;
							dv[1][0]=Na_x; dv[1][1]=Na_y; dv[1][2]=0.0;
							dv[2][0]=0.0;  dv[2][1]=0.0;  dv[2][2]=0.0;

							d2v[0][0][0]=0.00;  d2v[0][0][1]=0.00;  d2v[0][0][2]=0.0; d2v[0][1][0]=0.0;   d2v[0][1][1]=0.0;   d2v[0][1][2]=0.0; d2v[0][2][0]=0.0; d2v[0][2][1]=0.0; d2v[0][2][2]=0.0;
							d2v[1][0][0]=Na_xx; d2v[1][0][1]=Na_xy; d2v[1][0][2]=0.0; d2v[1][1][0]=Na_yx; d2v[1][1][1]=Na_yy; d2v[1][1][2]=0.0; d2v[1][2][0]=0.0; d2v[1][2][1]=0.0; d2v[1][2][2]=0.0;
							d2v[2][0][0]=0.00;  d2v[2][0][1]=0.00;  d2v[2][0][2]=0.0; d2v[2][1][0]=0.0;   d2v[2][1][1]=0.0;   d2v[2][1][2]=0.0; d2v[2][2][0]=0.0; d2v[2][2][1]=0.0; d2v[2][2][2]=0.0;
							
						}

						for (j=0; j<dof; j++)
						{
							if (j==0)
							{
								dz[0][0]=Nb_x; dz[0][1]=Nb_y; dz[0][2]=0.0;
								dz[1][0]=0.0;  dz[1][1]=0.0;  dz[1][2]=0.0;
								dz[2][0]=0.0;  dz[2][1]=0.0;  dz[2][2]=0.0;

								d2z[0][0][0]=Nb_xx; d2z[0][0][1]=Nb_xy; d2z[0][0][2]=0.0; d2z[0][1][0]=Nb_yx; d2z[0][1][1]=Nb_yy; d2z[0][1][2]=0.0; d2z[0][2][0]=0.0; d2z[0][2][1]=0.0; d2z[0][2][2]=0.0;
								d2z[1][0][0]=0.00;  d2z[1][0][1]=0.00;  d2z[1][0][2]=0.0; d2z[1][1][0]=0.0;   d2z[1][1][1]=0.0;   d2z[1][1][2]=0.0; d2z[1][2][0]=0.0; d2z[1][2][1]=0.0; d2z[1][2][2]=0.0;
								d2z[2][0][0]=0.00;  d2z[2][0][1]=0.00;  d2z[2][0][2]=0.0; d2z[2][1][0]=0.0;   d2z[2][1][1]=0.0;   d2z[2][1][2]=0.0; d2z[2][2][0]=0.0; d2z[2][2][1]=0.0; d2z[2][2][2]=0.0;
							}
							else if(j==1)
							{
								dz[0][0]=0.0;  dz[0][1]=0.0;  dz[0][2]=0.0;
								dz[1][0]=Nb_x; dz[1][1]=Nb_y; dz[1][2]=0.0;
								dz[2][0]=0.0;  dz[2][1]=0.0;  dz[2][2]=0.0;

								d2z[0][0][0]=0.00;  d2z[0][0][1]=0.00;  d2z[0][0][2]=0.0; d2z[0][1][0]=0.0;   d2z[0][1][1]=0.0;   d2z[0][1][2]=0.0; d2z[0][2][0]=0.0; d2z[0][2][1]=0.0; d2z[0][2][2]=0.0;
								d2z[1][0][0]=Nb_xx; d2z[1][0][1]=Nb_xy; d2z[1][0][2]=0.0; d2z[1][1][0]=Nb_yx; d2z[1][1][1]=Nb_yy; d2z[1][1][2]=0.0; d2z[1][2][0]=0.0; d2z[1][2][1]=0.0; d2z[1][2][2]=0.0;
								d2z[2][0][0]=0.00;  d2z[2][0][1]=0.00;  d2z[2][0][2]=0.0; d2z[2][1][0]=0.0;   d2z[2][1][1]=0.0;   d2z[2][1][2]=0.0; d2z[2][2][0]=0.0; d2z[2][2][1]=0.0; d2z[2][2][2]=0.0;
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

							for (k=0; k<3; k++)
							{
								for (l=0; l<3; l++)
								{
									for (u=0; u<3; u++)
									{
										Keq[a][i][b][j]+= 0.5*eps*(-d2z[k][l][u]-d2z[l][k][u])*0.5*(d2v[k][l][u]+d2v[l][k][u])
														 -0.5*eps*(-d2z[l][k][u]+d2z[k][l][u])*d2v[k][l][u];
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
				PetscReal Na_xx =N2[a][0][0];	PetscReal Na_yy =N2[a][1][1];
				PetscReal Na_xy =N2[a][0][1];	PetscReal Na_yx =N2[a][1][0];

				for (i=0; i<dof; i++)
				{
					if (i==0)
					{
						v[0]=Na; 	   v[1]=0.0; 	  v[2]=0.0;
						dv[0][0]=Na_x; dv[0][1]=Na_y; dv[0][2]=0.0;
						dv[1][0]=0.0;  dv[1][1]=0.0;  dv[1][2]=0.0;
						dv[2][0]=0.0;  dv[2][1]=0.0;  dv[2][2]=0.0;

						d2v[0][0][0]=Na_xx; d2v[0][0][1]=Na_xy; d2v[0][0][2]=0.0; d2v[0][1][0]=Na_yx; d2v[0][1][1]=Na_yy; d2v[0][1][2]=0.0; d2v[0][2][0]=0.0; d2v[0][2][1]=0.0; d2v[0][2][2]=0.0;
						d2v[1][0][0]=0.00;  d2v[1][0][1]=0.00;  d2v[1][0][2]=0.0; d2v[1][1][0]=0.0;   d2v[1][1][1]=0.0;   d2v[1][1][2]=0.0; d2v[1][2][0]=0.0; d2v[1][2][1]=0.0; d2v[1][2][2]=0.0;
						d2v[2][0][0]=0.00;  d2v[2][0][1]=0.00;  d2v[2][0][2]=0.0; d2v[2][1][0]=0.0;   d2v[2][1][1]=0.0;   d2v[2][1][2]=0.0; d2v[2][2][0]=0.0; d2v[2][2][1]=0.0; d2v[2][2][2]=0.0;
					}
					else if (i==1)
					{
						v[0]=0.0; 	   v[1]=Na; 	  v[2]=0.0;
						dv[0][0]=0.0;  dv[0][1]=0.0;  dv[0][2]=0.0;
						dv[1][0]=Na_x; dv[1][1]=Na_y; dv[1][2]=0.0;
						dv[2][0]=0.0;  dv[2][1]=0.0;  dv[2][2]=0.0;

						d2v[0][0][0]=0.00;  d2v[0][0][1]=0.00;  d2v[0][0][2]=0.0; d2v[0][1][0]=0.0;   d2v[0][1][1]=0.0;   d2v[0][1][2]=0.0; d2v[0][2][0]=0.0; d2v[0][2][1]=0.0; d2v[0][2][2]=0.0;
						d2v[1][0][0]=Na_xx; d2v[1][0][1]=Na_xy; d2v[1][0][2]=0.0; d2v[1][1][0]=Na_yx; d2v[1][1][1]=Na_yy; d2v[1][1][2]=0.0; d2v[1][2][0]=0.0; d2v[1][2][1]=0.0; d2v[1][2][2]=0.0;
						d2v[2][0][0]=0.00;  d2v[2][0][1]=0.00;  d2v[2][0][2]=0.0; d2v[2][1][0]=0.0;   d2v[2][1][1]=0.0;   d2v[2][1][2]=0.0; d2v[2][2][0]=0.0; d2v[2][2][1]=0.0; d2v[2][2][2]=0.0;
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
								Feq[a][i]+=M[k]*e[k][l][u]*dv[u][l];
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
									Feq[a][i]+=-0.5*(C[k][l][u][w]*(-fullChi[u][w])+C[l][k][u][w]*(-fullChi[u][w]))*0.5*(dv[k][l]+dv[l][k]);
								}
							}
						}
					}

					for (k=0; k<3; k++)
					{
						for (l=0; l<3; l++)
						{
							for (u=0; u<3; u++)
							{
								Feq[a][i]+=-0.25*eps*(-full_dChi[k][l][u]-full_dChi[k][u][l]-full_dChi[l][k][u]-full_dChi[l][u][k])*0.5*(d2v[k][l][u]+d2v[l][k][u])
										   +0.25*eps*(-full_dChi[l][k][u]+full_dChi[k][l][u]-full_dChi[l][u][k]+full_dChi[k][u][l])*d2v[k][l][u]
										   //Next part is the contribution from having S in the energy function
										   -0.25*eps*(-fullS[k][l][u]-fullS[k][u][l]-fullS[l][k][u]-fullS[l][u][k])*0.5*(d2v[k][l][u]+d2v[l][k][u])
										   +0.25*eps*(-fullS[l][k][u]+fullS[k][l][u]-fullS[l][u][k]+fullS[k][u][l])*d2v[k][l][u];
							}
						}
					}

				}
			}
		}
		return 0;
	}
//

//System for L2 projection of V^{S}
	#undef  __FUNCT__
	#define __FUNCT__ "VS"
	PetscErrorCode VS(IGAPoint p,IGAPoint pChi,IGAPoint pZu,IGAPoint pS,PetscReal *K,PetscReal *F, PetscReal *Chi,PetscReal *Zu,PetscReal *S,void *ctx)		//When other fields like Pi or S (etc.) come into this, declare a IGAPoint pPi or pS and a new PestcReal *UPi or *US for each
	{
		const PetscReal *N0;
		IGAPointGetShapeFuns(p,0,(const PetscReal**)&N0);									//Value of the shape functions
		PetscInt a,b,i,j,k,l,m,nen=p->nen, dof=p->dof;

		//Change to consider G=1
		//const PetscReal nu=0.33;
		const PetscReal mu=1.0;
		//const PetscReal lambda=2.0*mu*nu/(1.0-2.0*nu);
		const PetscReal eps=mu/100.0;														//Choose later based on whatever Amit says :)

		PetscReal x[2];																		//Vector of reals, size equal to problem's dimension
		IGAPointFormGeomMap(p,x);															//Fills x with the coordinates of p, Gauss's point

		PetscReal S0[8];																	//Assign S to a vector
		PetscReal dS[8][2];																	//And its derivative
		IGAPointFormValue(pS,S,&S0[0]);																
		IGAPointFormGrad (pS,S,&dS[0][0]);													//Same for the gradient

		//PetscReal chi0[4];																	//Array to contain the vector chi(0)
		PetscReal d2_Chi0[4][2][2];															//Same for its Hessian
		//IGAPointFormValue(pChi,Chi,&chi0[0]);												//Assign chi to its container
		IGAPointFormHess (pChi,Chi,&d2_Chi0[0][0][0]);										//This should be the 3-rd order tensor Chi_{i,jk} (remember that we are storing Chi_{kl} as a column vector)

		PetscReal d3_Z0[2][2][2][2];														//Array to store 3rd order partial derivative of z
		IGAPointFormDer3 (pZu,Zu,&d3_Z0[0][0][0][0]);										//Calculate and store in array

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
									FVS[a][i]+=eps*(fulld_S[j][k][l][l]+fulld3_z[j][k][l][l]+fulld2_Chi[j][k][l][l])*fullS[j][k][m]*v[m];
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

//System for L2 projection of V^{S} smooth (Integrates Vs multiplied with xi(S), where xi(S)=1 if any component of S is different than 0, and 0 in other cases)
	#undef  __FUNCT__
	#define __FUNCT__ "Int_Xi_Vs"
	//IntValpha(pointVa,pointAlp,FpointVaInt,PointInt1,PointInt2,Va0Values,Al1Va,NULL);CHKERRQ(ierr);		//When other fields like Pi or S (etc.) come into this, declare a IGAPoint pPi or pS and a new PescReal *UPi or *US for each
	PetscErrorCode Int_Xi_Vs(IGAPoint pV,IGAPoint pS,PetscReal *FInt1a,PetscReal *FInt2a,PetscReal *FInt,PetscReal *VS,PetscReal *US,void *ctx)		//When other fields like Pi or S (etc.) come into this, declare a IGAPoint pPi or pS and a new PetscReal *UPi or *US for each
	{
		PetscInt dof=pV->dof;
		PetscInt i;

		PetscReal Vs[2];																	//Create array to recieve Alfa
		IGAPointFormValue(pV,VS,&Vs[0]);													//This fills the values

		PetscReal S[8];																		//Create array to recieve Alfa
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

		PetscReal S[8];																		//Create array to recieve Alfa
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

//System for updting S using equation for Sdot (Under construction)
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
		fullVs[0]=0.0; fullVs[1]=0.5;
		full_dVs[0][0]=0.0; full_dVs[0][1]=0.0;
		full_dVs[1][0]=0.0; full_dVs[1][1]=0.0;

		PetscReal S[8];																		//Create array to recieve S
		IGAPointFormValue(pS,US,&S[0]);														//This fills the values of S (remember that S has 8 non zero components in 2D)

		PetscReal fullS[3][3][3]={0};
		fullS[0][0][0]=S[0]; fullS[0][0][1]=S[1];											//Expand S to full 3rd order form, only non-zero elements
		fullS[0][1][0]=S[2]; fullS[0][1][1]=S[3];
		fullS[1][0][0]=S[4]; fullS[1][0][1]=S[5];
		fullS[1][1][0]=S[6]; fullS[1][1][1]=S[7];

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

				for (b=0; b<nen; b++) 
				{
					PetscReal Nb=N0[b];
					PetscReal Nb_x = N1[b][0];		PetscReal Nb_y = N1[b][1];

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
							St[0][0][0]=0.0;		St[0][0][1]=Nb;		St[0][1][0]=0.0;	St[0][1][1]=0.0;	St[1][0][0]=0.0;	St[1][0][1]=0.0;	St[1][1][0]=0.0;	St[1][1][1]=0.0;
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
							St[0][0][0]=0.0;		St[0][0][1]=0.0;		St[0][1][0]=Nb;	St[0][1][1]=0.0;	St[1][0][0]=0.0;	St[1][0][1]=0.0;	St[1][1][0]=0.0;	St[1][1][1]=0.0;
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
							St[0][0][0]=0.0;		St[0][0][1]=0.0;		St[0][1][0]=0.0;	St[0][1][1]=Nb;	St[1][0][0]=0.0;	St[1][0][1]=0.0;	St[1][1][0]=0.0;	St[1][1][1]=0.0;
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
							St[0][0][0]=0.0;		St[0][0][1]=0.0;		St[0][1][0]=0.0;	St[0][1][1]=0.0;	St[1][0][0]=Nb;	St[1][0][1]=0.0;	St[1][1][0]=0.0;	St[1][1][1]=0.0;
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
							St[0][0][0]=0.0;		St[0][0][1]=0.0;		St[0][1][0]=0.0;	St[0][1][1]=0.0;	St[1][0][0]=0.0;	St[1][0][1]=Nb;	St[1][1][0]=0.0;	St[1][1][1]=0.0;
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
							St[0][0][0]=0.0;		St[0][0][1]=0.0;		St[0][1][0]=0.0;	St[0][1][1]=0.0;	St[1][0][0]=0.0;	St[1][0][1]=0.0;	St[1][1][0]=Nb;	St[1][1][1]=0.0;
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
							St[0][0][0]=0.0;		St[0][0][1]=0.0;		St[0][1][0]=0.0;	St[0][1][1]=0.0;	St[1][0][0]=0.0;	St[1][0][1]=0.0;	St[1][1][0]=0.0;	St[1][1][1]=Nb;
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
										KSdot[a][i][b][j]+=St[k][l][u]*v[k][l][u];
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
											KSdot[a][i][b][j]+=dt*St[k][l][w]*fullVs[w]*dv[k][l][u][u];
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
										KSdot[a][i][b][j]+=St[k][l][u]*v[k][l][u];
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
											KSdot[a][i][b][j]+=-dt*St[k][l][u]*(dv[k][l][w][u]*fullVs[w]+v[k][l][w]*full_dVs[w][u])-dt*(dSt[k][l][w][u]*fullVs[w]+St[k][l][w]*full_dVs[w][u])*v[k][l][u]
															   +dt*dt*(dSt[k][l][w][u]*fullVs[w]+St[k][l][w]*full_dVs[w][u])*(dv[k][l][w][u]*fullVs[w]+v[k][l][w]*full_dVs[w][u]);
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

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc, char *argv[]) {

//Creation of solution systems
	PetscErrorCode  ierr;
	ierr = PetscInitialize(&argc,&argv,0,0);CHKERRQ(ierr);										//Always initialize PETSc
	PetscInt dir,side;
	PetscPrintf(PETSC_COMM_WORLD,"Start of PruebaV5S \n");

	PetscInt commsize,rank;
	ierr = MPI_Comm_size(PETSC_COMM_WORLD,&commsize);CHKERRQ(ierr);
	ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank);CHKERRQ(ierr);
//

//App context creation and some data
	//Mesh parameters (to fix specific points in z0 system)
	PetscInt b=400;				//Parmeter to choose size of cores, must always be odd, core will be of size 1 unit, rest of the body will be of size b-1 units in each direction
	PetscReal Lx=20.0;
	PetscReal Ly=20.0;
	PetscInt  nx=b;
	PetscInt  ny=b;
	
	PetscReal alpha=0.09;
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
		dest   = fopen("../Results/PruebaV5S.c","w");

		while (0 < (bytes = fread(buffer, 1, sizeof(buffer), source)))
			fwrite(buffer, 1, bytes, dest);

		fclose(source);
		fclose(dest);
	}
//

//
//Recall that, for the systems to converge, dt<h/2, h=min(Lx/nx,Ly/ny) and the "less than" is strict 
//

for(i=0;i<200;i++)
{
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
	//System for chiS
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
			if(pointchiS->atboundary==0 && pointS->atboundary==0)
			{
				while (IGAElementNextPoint(elemchiS,pointchiS)) 
				{
					IGAElementNextPoint(elemS,pointS);
					ierr = IGAPointGetWorkMat(pointchiS,&KpointchiS);CHKERRQ(ierr);
					ierr = IGAPointGetWorkVec(pointchiS,&FpointchiS);CHKERRQ(ierr);
					ierr = curlChiS(pointchiS,pointS,KpointchiS,FpointchiS,S0chiS,NULL);CHKERRQ(ierr);
					ierr = IGAPointAddMat(pointchiS,KpointchiS,KlocchiS);CHKERRQ(ierr);
					ierr = IGAPointAddVec(pointchiS,FpointchiS,FlocchiS);CHKERRQ(ierr);
				}
				while (pointS->index != -1)
				{
					IGAElementNextPoint(elemS,pointS);
				}
				ierr = IGAElementEndPoint(elemchiS,&pointchiS);CHKERRQ(ierr);
				ierr = IGAElementEndPoint(elemS,&pointS);CHKERRQ(ierr);
			}
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

	ierr = KSPDestroy(&kspchiS);CHKERRQ(ierr);
	ierr = MatDestroy(&KchiS);CHKERRQ(ierr);
	ierr = VecDestroy(&FchiS);CHKERRQ(ierr);
	
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
	PetscPrintf(PETSC_COMM_WORLD,"\nSystem for L2 projection for Alfa+Sp:X starting \n\n");
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
			ierr = IGASetBoundaryForm(igaAl,dir,side,PETSC_TRUE);CHKERRQ(ierr);  				// Neumann boundary conditions
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
			if(pointAlp->atboundary==0 && pointSp->atboundary==0)
			{
				while (IGAElementNextPoint(elemAlp,pointAlp)) 
				{
					IGAElementNextPoint(elemSp,pointSp);
					ierr = IGAPointGetWorkMat(pointAlp,&KpointAlp);CHKERRQ(ierr);
					ierr = IGAPointGetWorkVec(pointAlp,&FpointAlp);CHKERRQ(ierr);
					ierr = L2ProjectionAlSp(pointAlp,pointSp,KpointAlp,FpointAlp,S0Alp,&user);CHKERRQ(ierr);
					ierr = IGAPointAddMat(pointAlp,KpointAlp,KlocAlp);CHKERRQ(ierr);
					ierr = IGAPointAddVec(pointAlp,FpointAlp,FlocAlp);CHKERRQ(ierr);
				}
				IGAElementNextPoint(elemSp,pointSp);
				ierr = IGAElementEndPoint(elemAlp,&pointAlp);CHKERRQ(ierr);
				ierr = IGAElementEndPoint(elemSp,&pointSp);CHKERRQ(ierr);
			}
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

	ierr = KSPDestroy(&kspAlp);CHKERRQ(ierr);
	ierr = MatDestroy(&KAlp);CHKERRQ(ierr);
	ierr = VecDestroy(&FAlp);CHKERRQ(ierr);

	//Si hay problemas, borrar esto
	ierr = VecDestroy(&chiS0); CHKERRQ(ierr);
	//hasta aquí
	
	char nameAlp[512];//="/Alp-2d-0.dat";
	sprintf(nameAlp,"%s%d%s","/Alp-2d-",i,".dat");
	char pathAlp[512];
	sprintf(pathAlp,"%s%s",direct,nameAlp);
	ierr = IGAWriteVec(igaAl,alp0,pathAlp);CHKERRQ(ierr);

	ierr = IGADestroy(&igachiS);CHKERRQ(ierr);		//This system is one of the more memory intensive. It's not required from here down, so better to destroy it.
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
	ierr = IGASetOrder(igachiUp,2);CHKERRQ(ierr);													//Number of spatial derivatives to calculate
	ierr = IGASetFromOptions(igachiUp);CHKERRQ(ierr);												//Note: The order (or degree) of the shape functions is given by the mesh!
	ierr = IGARead(igachiUp,"./geometry2.dat");CHKERRQ(ierr);
	
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

			if(pointchiUp->atboundary==0 && pointAlp->atboundary==0)
			{
				while (IGAElementNextPoint(elemchiUp,pointchiUp)) 
				{
					IGAElementNextPoint(elemAlp,pointAlp);
					ierr = IGAPointGetWorkMat(pointchiUp,&KpointchiUp);CHKERRQ(ierr);
					ierr = IGAPointGetWorkVec(pointchiUp,&FpointchiUp);CHKERRQ(ierr);
					ierr = curlChiU(pointchiUp,pointAlp,KpointchiUp,FpointchiUp,Al0chiUp,NULL);CHKERRQ(ierr);
					ierr = IGAPointAddMat(pointchiUp,KpointchiUp,KlocchiUp);CHKERRQ(ierr);
					ierr = IGAPointAddVec(pointchiUp,FpointchiUp,FlocchiUp);CHKERRQ(ierr);
				}
				IGAElementNextPoint(elemAlp,pointAlp);

				ierr = IGAElementEndPoint(elemchiUp,&pointchiUp);CHKERRQ(ierr);
				ierr = IGAElementEndPoint(elemAlp,&pointAlp);CHKERRQ(ierr);
			}
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

	ierr = KSPDestroy(&kspchiUp);CHKERRQ(ierr);
	ierr = MatDestroy(&KchiUp);CHKERRQ(ierr);
	ierr = VecDestroy(&FchiUp);CHKERRQ(ierr);

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
	ierr = IGASetOrder(igaZ0,3);CHKERRQ(ierr);													//Number of spatial derivatives to calculate
	ierr = IGASetFromOptions(igaZ0);CHKERRQ(ierr);												//Note: The order (or degree) of the shape functions is given by the mesh!
	ierr = IGARead(igaZ0,"./geometry3.dat");CHKERRQ(ierr);
	
	for (dir=0; dir<2; dir++)
	{
		ierr = IGASetRuleType(igaZ0,dir,IGA_RULE_LEGENDRE);CHKERRQ(ierr);
		ierr = IGASetRuleSize(igaZ0,dir,6);CHKERRQ(ierr);
	}
	ierr = IGASetUp(igaZ0);CHKERRQ(ierr);

	PetscInt fijaPunto=0;																		//Fix a single point (1) or a side (chosen in blocks below)

	ierr = IGASetUp(igaZ0);CHKERRQ(ierr);
	ierr = IGASetUp(igachiUp);CHKERRQ(ierr);
	ierr = IGASetUp(igaS);CHKERRQ(ierr);


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
		//If we are not fixing a single point, set dirichlet conditions here
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
	const PetscReal	*arrayChi0Z0, *arraySZ0;	//arrayU
	Vec				localChi0Z0, localSZ0;		//localU
	PetscReal		*Chi0Z0, *SZ0;					//U

  	IGAFormSystem	wtfZ0;
 	void			*wtf2Z0;

 	KSP kspZ0;
	ierr = IGACreateKSP(igaZ0,&kspZ0);CHKERRQ(ierr);

	// Get local vectors Chi0, S and arrays
	ierr = IGAGetLocalVecArray(igachiUp,chiUp0,&localChi0Z0,&arrayChi0Z0);CHKERRQ(ierr);
	ierr = IGAGetLocalVecArray(igaS,s0,&localSZ0,&arraySZ0);CHKERRQ(ierr);

	// Element loop
	ierr = IGABeginElement(igaZ0,&elemZ0);CHKERRQ(ierr);
	ierr = IGABeginElement(igachiUp,&elemchiUp);CHKERRQ(ierr);
	ierr = IGABeginElement(igaS,&elemS);CHKERRQ(ierr);

	while (IGANextElement(igaZ0,elemZ0)) 
	{
		IGANextElement(igachiUp,elemchiUp);
		IGANextElement(igaS,elemS);

		ierr = IGAElementGetWorkMat(elemZ0,&KlocZ0);CHKERRQ(ierr);
		ierr = IGAElementGetWorkVec(elemZ0,&FlocZ0);CHKERRQ(ierr);
		ierr = IGAElementGetValues(elemchiUp,arrayChi0Z0,&Chi0Z0);CHKERRQ(ierr);
		ierr = IGAElementGetValues(elemS,arraySZ0,&SZ0);CHKERRQ(ierr);

		// FormSystem loop
		while (IGAElementNextFormSystem(elemZ0,&wtfZ0,&wtf2Z0)) 
		{
		// Quadrature loop
			ierr = IGAElementBeginPoint(elemZ0,&pointZ0);CHKERRQ(ierr);
			ierr = IGAElementBeginPoint(elemchiUp,&pointchiUp);CHKERRQ(ierr);
			ierr = IGAElementBeginPoint(elemS,&pointS);CHKERRQ(ierr);

			while (IGAElementNextPoint(elemZ0,pointZ0))
			{
				if(pointZ0->atboundary==1)
				{
					ierr = IGAPointGetWorkMat(pointZ0,&KpointZ0);CHKERRQ(ierr);
					ierr = IGAPointGetWorkVec(pointZ0,&FpointZ0);CHKERRQ(ierr);
					//	   Z0sys(IGAPoint p, IGAPoint pChi, IGAPoint pS,IGAPoint pZ, PetscReal *K, PetscReal *F, PetscReal *UChi, PetscReal *S,PetscReal *ZS, void *ctx)
					ierr = Z0sys(pointZ0,pointchiUp,pointS,KpointZ0,FpointZ0,Chi0Z0,SZ0,NULL);CHKERRQ(ierr);
					ierr = IGAPointAddMat(pointZ0,KpointZ0,KlocZ0);CHKERRQ(ierr);
					ierr = IGAPointAddVec(pointZ0,FpointZ0,FlocZ0);CHKERRQ(ierr);
				}

				if(pointchiUp->atboundary==1)
				{
					PetscPrintf(PETSC_COMM_WORLD,"Hola");
				}

				if(pointZ0->atboundary==0 && pointchiUp->atboundary==0 && pointS->atboundary==0)
				{
					IGAElementNextPoint(elemchiUp,pointchiUp);
					IGAElementNextPoint(elemS,pointS);

					ierr = IGAPointGetWorkMat(pointZ0,&KpointZ0);CHKERRQ(ierr);
					ierr = IGAPointGetWorkVec(pointZ0,&FpointZ0);CHKERRQ(ierr);
					//	   Z0sys(IGAPoint p, IGAPoint pChi, IGAPoint pS,IGAPoint pZ, PetscReal *K, PetscReal *F, PetscReal *UChi, PetscReal *S,PetscReal *ZS, void *ctx)
					ierr = Z0sys(pointZ0,pointchiUp,pointS,KpointZ0,FpointZ0,Chi0Z0,SZ0,NULL);CHKERRQ(ierr);
					ierr = IGAPointAddMat(pointZ0,KpointZ0,KlocZ0);CHKERRQ(ierr);
					ierr = IGAPointAddVec(pointZ0,FpointZ0,FlocZ0);CHKERRQ(ierr);
				}
			}
			while (pointchiUp->index != -1)
			{
				IGAElementNextPoint(elemchiUp,pointchiUp);
			}
			while (pointS->index != -1)
			{
				IGAElementNextPoint(elemS,pointS);
			}

			ierr = IGAElementEndPoint(elemZ0,&pointZ0);CHKERRQ(ierr);
			ierr = IGAElementEndPoint(elemchiUp,&pointchiUp);CHKERRQ(ierr);
			ierr = IGAElementEndPoint(elemS,&pointS);CHKERRQ(ierr);
		}

		if(fijaPunto==0)
		{
			ierr = IGAElementFixSystem(elemZ0,KlocZ0,FlocZ0);CHKERRQ(ierr);					//This sets Dirichlet condition ¿? (Yes, this applies the conditions from IGASetBoundaryValue)
		}
		ierr = IGAElementAssembleMat(elemZ0,KlocZ0,KZ0);CHKERRQ(ierr);
		ierr = IGAElementAssembleVec(elemZ0,FlocZ0,FZ0);CHKERRQ(ierr);

	}
	IGANextElement(igachiUp,elemchiUp);
	IGANextElement(igaS,elemS);

	ierr = IGAEndElement(igaZ0,&elemZ0);CHKERRQ(ierr);
	ierr = IGAEndElement(igachiUp,&elemchiUp);CHKERRQ(ierr);
	ierr = IGAEndElement(igaS,&elemS);CHKERRQ(ierr);

	// Restore local vectors Chi0 and arrays
	ierr = IGARestoreLocalVecArray(igachiUp,chiUp0,&localChi0Z0,&arrayChi0Z0);CHKERRQ(ierr);
	ierr = IGARestoreLocalVecArray(igaS,s0,&localSZ0,&arraySZ0);CHKERRQ(ierr);

	ierr = MatAssemblyBegin(KZ0,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyEnd  (KZ0,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

	if (fijaPunto==1)
	{
		//Here we set values to the Matrix directly, to impose Dirichlet condition in a single point.
		//Note: Lower Left corner is gdl's  0 and 1,
		//		Lower Right corner is gdl's 2*(nx+2)-2 and 2*(nx+2)-1 
		//		Upper Left corner is gdl's  2*(nx+2)*(ny+2)-2*(nx+1)-2 and 2*(nx+2)*(ny+2)-2*(nx+1)-1
		//		Upper Right corner is gdl's 2*(nx+2)*(ny+2)-2 and 2*(nx+2)*(ny+2)-1
		//All of these for when z is a 2nd order nurbs
		PetscInt n,m;
		ierr = MatGetSize(KZ0,&n,&m);CHKERRQ(ierr);
		ierr = MatSetOption(KZ0, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);CHKERRQ(ierr);

		PetscInt rows, *cols;
		PetscReal val, *vals;

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
	ierr = IGASetOrder(igaVs,2);CHKERRQ(ierr);													//Number of spatial derivatives to calculate
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
					//IGAElementNextPoint(elemZ0,pointZ0);
					//IGAElementNextPoint(elemchiUp,pointchiUp);
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

	//ierr = KSPDestroy(&kspVs);CHKERRQ(ierr);
	//ierr = MatDestroy(&KVs);CHKERRQ(ierr);
	ierr = VecDestroy(&FVs);CHKERRQ(ierr);

	char nameVs[512];//="/Vs-2d-0.dat";
	sprintf(nameVs,"%s%d%s","/Vs-2d-",i,".dat");
	char pathVs[512];
	sprintf(pathVs,"%s%s",direct,nameVs);
	ierr = IGAWriteVec(igaVs,Vs0,pathVs);CHKERRQ(ierr);	
//

//System for L2 proyection of smoothed V^{S}
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

	ierr = PetscPrintf(PETSC_COMM_WORLD,"Hola1 \n");CHKERRQ(ierr);

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

	ierr = PetscPrintf(PETSC_COMM_WORLD,"Hola2 \n");CHKERRQ(ierr);

	ierr = IGAEndElement(igaVs,&elemVs);CHKERRQ(ierr);
	ierr = IGAEndElement(igaS,&elemS);CHKERRQ(ierr);

	ierr = IGARestoreLocalVecArray(igaS,s0,&localSVs,&arraySVs);CHKERRQ(ierr);
	ierr = IGARestoreLocalVecArray(igaVs,Vs0,&localVs0,&arrayVs0);CHKERRQ(ierr);

	//Here we have a vector with an L2 Projection of the integrated velocity.
	ierr = VecAssemblyBegin(FVsSmooth);CHKERRQ(ierr);
	ierr = VecAssemblyEnd  (FVsSmooth);CHKERRQ(ierr);

	ierr = PetscPrintf(PETSC_COMM_WORLD,"Hola3 \n");CHKERRQ(ierr);

	//Here we solved the L2 projection system
	ierr = KSPSolve(kspVs,FVsSmooth,VsSmooth);CHKERRQ(ierr);

	char nameVsS1[512];//="/VsSmooth-2d-0.dat";
	sprintf(nameVsS1,"%s%d%s","/VsSmooth-2d-",i,".dat");
	char pathVsS1[512];
	sprintf(pathVsS1,"%s%s",direct,nameVsS1);
	ierr = IGAWriteVec(igaVs,VsSmooth,pathVsS1);CHKERRQ(ierr);

	ierr = KSPDestroy(&kspVs);CHKERRQ(ierr);
	ierr = MatDestroy(&KVs);CHKERRQ(ierr);
//

//Systm for evolution of S
	PetscPrintf(PETSC_COMM_WORLD,"\nSystem for S_dot starting \n\n");
	T=time(NULL);
	tm=*localtime(&T);
	PetscPrintf(PETSC_COMM_WORLD,"Current time is %02d:%02d:%02d \n",tm.tm_hour,tm.tm_min,tm.tm_sec);
	IGA igaSdot;
	ierr = IGACreate(PETSC_COMM_WORLD,&igaSdot);CHKERRQ(ierr);
	ierr = IGASetDim(igaSdot,2);CHKERRQ(ierr);													//Spatial dimension of the problem
	ierr = IGASetDof(igaSdot,8);CHKERRQ(ierr);													//Number of degrees of freedom, per node
	ierr = IGASetOrder(igaSdot,2);CHKERRQ(ierr);													//Number of spatial derivatives to calculate
	ierr = IGASetFromOptions(igaSdot);CHKERRQ(ierr);												//Note: The order (or degree) of the shape functions is given by the mesh!
	ierr = IGARead(igaSdot,"./geometry.dat");CHKERRQ(ierr);
	
	for (dir=0; dir<2; dir++)
	{
		ierr = IGASetRuleType(igaSdot,dir,IGA_RULE_LEGENDRE);CHKERRQ(ierr);
		ierr = IGASetRuleSize(igaSdot,dir,6);CHKERRQ(ierr);
	}
	ierr = IGASetUp(igaSdot);CHKERRQ(ierr);

	for (dir=0; dir<2; dir++) 
	{
		for (side=0; side<2; side++) 
		{
			//ierr = IGASetBoundaryValue(iga,dir,side,dof,0.0);CHKERRQ(ierr);					// Dirichlet boundary conditions
			ierr = IGASetBoundaryForm(igaSdot,dir,side,PETSC_TRUE);CHKERRQ(ierr);  				// Neumann boundary conditions
		}
	}

	Mat KSdot;
	Vec Sdot,FSdot;

	ierr = IGACreateMat(igaSdot,&KSdot);CHKERRQ(ierr);
	ierr = IGACreateVec(igaSdot,&Sdot);CHKERRQ(ierr);
	ierr = IGACreateVec(igaSdot,&FSdot);CHKERRQ(ierr);

	IGAPoint		pointSdot;															//point
	IGAElement		elemSdot;																//element
	PetscReal		*KlocSdot,*FlocSdot;												//AA y BB
	PetscReal		*KpointSdot,*FpointSdot;											//KKK y FFF
	const PetscReal	*arrayVsSdot,*arrayS0Sdot;			//arrayU
	Vec				localVsSdot,localS0Sdot;				//localU
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
					//VS(IGAPoint p,IGAPoint pChi,IGAPoint pZu,IGAPoint pS,PetscReal *K,PetscReal *F, PetscReal *Chi,PetscReal *Zu,PetscReal *S,void *ctx)
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

		//ierr = IGAElementFixSystem(elemStress,KlocStress,FlocStress);CHKERRQ(ierr);					//This sets Dirichlet condition ¿? (Yes, this applies the conditions from IGASetBoundaryValue)
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
	ierr = VecDestroy(&FSdot);CHKERRQ(ierr);

	char nameSdot[512];//="/Input-S-2d-1.dat";
	sprintf(nameSdot,"%s%d%s","/Input-S-2d-",i+1,".dat");
	char pathSdot[512];
	sprintf(pathSdot,"%s%s",direct,nameSdot);
	ierr = IGAWriteVec(igaSdot,Sdot,pathSdot);CHKERRQ(ierr);	
//
}

/*
//Destroy all objects not needed anymore (Better to do it here in case different codes call the same IGA, move if memory is a problem)
	ierr = IGADestroy(&igaAl); CHKERRQ(ierr);
	ierr = IGADestroy(&igaZ0); CHKERRQ(ierr);
	ierr = IGADestroy(&igachiUp); CHKERRQ(ierr);
//
*/

ierr = PetscFinalize();CHKERRQ(ierr);

return 0;
}