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
		//PetscReal nx     = user->nx;
		//PetscReal ny     = user->ny;

		PetscInt a,b,i;
		PetscInt nen = p->nen;																//Number of shape functions
		PetscInt dim = p->dim;																//Spatial dimensions of the problem
		PetscInt dof = p->dof;

		PetscReal x[dim];																	//Vector of reals, size equal to problem's dimension
		IGAPointFormGeomMap(p,x);															//Fills x with the coordinates of p, Gauss's point

		//g is the function to L2 project
		PetscReal dx=Lx/201.0;																//The number of the denominator should come automatically
		PetscReal dy=Ly/201.0;
		PetscReal t1=0.075;
		PetscReal t2=0.075;
		PetscReal a2;

		//S has 8 components, in order S(1,1,1), S(1,1,2), S(1,2,1), S(1,2,2), S(2,1,1), S(2,1,2), S(2,2,1), S(2,2,2)

		PetscReal g[dof];

		/*
		g[0]=(x[0]>-0.25)*(x[0]<0.25)*(x[1]>-0.01)*(x[1]<0.01);;

		g[1]=(x[0]>-0.15)*(x[0]<0.35)*(x[1]>-0.01)*(x[1]<0.01);

		g[2]=(x[0]>-0.01)*(x[0]<0.01)*(x[1]>-0.2)*(x[1]<0.2);

		g[5]=(x[0]>-0.01)*(x[0]<0.01)*(x[1]>-0.2)*(x[1]<0.2);

		g[7]=(x[0]>-0.11)*(x[0]<-0.01)*(x[1]>-0.1)*(x[1]<0.3);
		*/
		
		PetscReal norm[3]={0.0,1.0,0.0};					//S is not 0 for (i,j,k)!=3, so norm(2)=0 and A(0,0),A(0,1),A(1,0),A(1,1)=!0
		PetscReal sTemp[3][3][3]={0};
		PetscReal A[3][3]={0};

		//a2=0.25*(1.0-tanh((x[0])/t2))/1.5;
		a2=1.0/5;
		A[0][0]=0.0;
		A[0][1]=-tan(5.0/180.0*ConstPi)*0.5*(tanh((x[1]+5*dy)/t1)-tanh((x[1]-5*dy)/t1))/dy*a2;
		A[1][0]=tan(5.0/180.0*ConstPi)*0.5*(tanh((x[1]+5*dy)/t1)-tanh((x[1]-5*dy)/t1))/dy*a2;
		A[1][1]=0.0;

		for (int k=0;k<2;k++)
		{
			for (int n=0;n<2;n++)
			{
				for (int m=0;m<2;m++)
				{
					sTemp[k][n][m]=A[k][n]*norm[m];
				}
			}
		}

		g[0]=sTemp[0][0][0]; g[1]=sTemp[0][0][1];
		g[2]=sTemp[0][1][0]; g[3]=sTemp[0][1][1];
		g[4]=sTemp[1][0][0]; g[5]=sTemp[1][0][1];
		g[6]=sTemp[1][1][0]; g[7]=sTemp[1][1][1];

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
		  			KK[a][i][b][i] = N[a]*N[b];				//Check this!!
		  		}
		  		FF[a][i] = N[a]*g[i];
			}
		}
		return 0;
	}
//

/*
//Helmholtz decomposition of S, curl part
	#undef  __FUNCT__
	#define __FUNCT__ "curlChiS"
	PetscErrorCode curlChiS(IGAPoint p,IGAPoint pPi,PetscReal *K,PetscReal *F,PetscReal *U,void *ctx)
	{
		const PetscReal *N0,(*N1)[2];
		IGAPointGetShapeFuns(p,0,(const PetscReal**)&N0);
		IGAPointGetShapeFuns(p,1,(const PetscReal**)&N1);

		PetscInt a,b,u,w,i,j,k,l,m,nen=p->nen, dof=p->dof;

		PetscReal Pi[4];																		//Create array to recieve S
		//PetscReal dS[8][2];																		//Create array to recieve dS
		IGAPointFormValue(pPi,U,&Pi[0]);														//This fills the values
		//IGAPointFormGrad (pS,U,&dS[0][0]);														//This fills the values of the derivatives

		PetscReal fullPi[3][3][3]={0};
		//PetscReal fulldS[3][3][3][3]={0};

		fullPi[0][0][2]=Pi[0]; fullPi[0][1][2]=Pi[1]; fullPi[1][0][2]=Pi[2]; fullPi[1][1][2]=Pi[3];		//Expand S to full 3rd order form, only non-zero elements

		//Definition of alternating tensor
		const PetscReal e[3][3][3]=
		{
			{{0.0,0.0,0.0},{0.0,0.0,1.0},{0.0,-1.0,0.0}},
			{{0.0,0.0,-1.0},{0.0,0.0,0.0},{1.0,0.0,0.0}},
			{{0.0,1.0,0.0},{-1.0,0.0,0.0},{0.0,0.0,0.0}}
		};

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
											//KchiS[a][u][b][w]+=dchiS[i][j][k][m]*(dv[i][j][k][m]-dv[i][j][m][k])+dchiS[i][j][k][k]*dv[i][j][m][m];
											KchiS[a][u][b][w]+=dchiS[i][j][k][m]*dv[i][j][k][m]-dchiS[i][j][k][m]*dv[i][j][m][k]+dchiS[i][j][k][k]*dv[i][j][m][m];

											
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
									for(l=0;l<3;l++)
									{
										FchiS[a][u]+=fullPi[i][j][k]*e[k][l][m]*dv[i][j][m][l];
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
*/

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

//System for L2 projection of stress (sym part)
	#undef  __FUNCT__
	#define __FUNCT__ "Stress"
	//PetscErrorCode Stress(IGAPoint p,IGAPoint pU, IGAPoint pHs,IGAPoint pChi,IGAPoint pZu,PetscReal *K,PetscReal *F,PetscReal *U,PetscReal *HS, PetscReal *Chi,PetscReal *Zu,void *ctx)		//When other fields like Pi or S (etc.) come into this, declare a IGAPoint pPi or pS and a new PescReal *UPi or *US for each
	PetscErrorCode Stress(IGAPoint p,IGAPoint pChi,IGAPoint pZu,IGAPoint pS, PetscReal *K,PetscReal *F, PetscReal *Chi,PetscReal *Zu,PetscReal *S, void *ctx)		//When other fields like Pi or S (etc.) come into this, declare a IGAPoint pPi or pS and a new PestcReal *UPi or *US for each
	{
		//const PetscReal *N0,(*N1)[2],(*N2)[2][2];
		//const PetscReal *N0,(*N1)[2];
		const PetscReal *N0;
		IGAPointGetShapeFuns(p,0,(const PetscReal**)&N0);									//Value of the shape functions
		//IGAPointGetShapeFuns(p,1,(const PetscReal**)&N1);									//Derivatives of the shape functions
		//IGAPointGetShapeFuns(p,2,(const PetscReal**)&N2);									//Second derivatives of the shape funcions
		PetscInt a,b,i,j,k,l,m,n,nen=p->nen, dof=p->dof;

		//E and nu should come from AppCtx in the future
		//const PetscReal E=2100000.0*9.81*10000.0;
		//const PetscReal nu=0.33;
		//const PetscReal lambda=(E*nu)/((1.0+nu)*(1.0-2.0*nu));
		//const PetscReal mu=E/(2.0*(1.0+nu));
		//const PetscReal eps=0.0*E/1000.0;														//Choose later based on whatever Amit says :)

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

						Fstress[a][i]+= -0.25*eps*(-fulld_S[m][n][k][k]-fulld_S[m][k][n][k]-fulld_S[n][m][k][k]-fulld_S[n][k][m][k])*v[m][n];
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
	PetscErrorCode ClassicStress(IGAPoint p,IGAPoint pChi,IGAPoint pZu,PetscReal *K,PetscReal *F, PetscReal *Chi,PetscReal *Zu,void *ctx)		//When other fields like Pi or S (etc.) come into this, declare a IGAPoint pPi or pS and a new PestcReal *UPi or *US for each
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
							Fstress[a][i]+=0.5*(C[m][n][k][l]*(fulld_z[k][l]-fullChi[k][l])+C[n][m][k][l]*(fulld_z[k][l]-fullChi[k][l]))*v[m][n];
						}
					}
				}
			}
		}
		return 0;
	}
//

//System for L2 projection of couple-stresses
	#undef  __FUNCT__
	#define __FUNCT__ "CoupleStress"
	//PetscErrorCode Stress(IGAPoint p,IGAPoint pU, IGAPoint pHs,IGAPoint pChi,IGAPoint pZu,PetscReal *K,PetscReal *F,PetscReal *U,PetscReal *HS, PetscReal *Chi,PetscReal *Zu,void *ctx)		//When other fields like Pi or S (etc.) come into this, declare a IGAPoint pPi or pS and a new PescReal *UPi or *US for each
	PetscErrorCode CoupleStress(IGAPoint p, IGAPoint pChi, IGAPoint pZu,IGAPoint pS, PetscReal *K, PetscReal *F, PetscReal *Chi,PetscReal *Zu,PetscReal *s0, void *ctx)		//When other fields like Pi or S (etc.) come into this, declare a IGAPoint pPi or pS and a new PestcReal *UPi or *US for each
	{
		const PetscReal *N0;
		IGAPointGetShapeFuns(p,0,(const PetscReal**)&N0);									//Value of the shape functions
		PetscInt a,b,i,k,l,m,n,nen=p->nen, dof=p->dof;

		PetscReal x[2];																		//Vector of reals, size equal to problem's dimension
		IGAPointFormGeomMap(p,x);															//Fills x with the coordinates of p, Gauss's point

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

		PetscReal S0[8];																	//Create an array to hold the values for S
		IGAPointFormValue(pS,s0,&S0[0]);													//Assign S to its container

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

		//Expand S to full 3rd order form, only non-zero elements
		PetscReal fullS[3][3][3]={0};
		fullS[0][0][0]=S0[0]; fullS[0][0][1]=S0[1];
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
							FCS[a][i]+=-0.5*eps*e[m][k][l]*(-fulld2_z[k][l][n]-fulld_Chi[k][l][n]-fulld2_z[k][n][l]-fulld_Chi[k][n][l])*v[m][n];
						}
					}
				}
			}
		}
		return 0;
	}
//

//System for L2 projection of energy density
	#undef  __FUNCT__
	#define __FUNCT__ "EnergyDensity"
	//PetscErrorCode Stress(IGAPoint p,IGAPoint pU, IGAPoint pHs,IGAPoint pChi,IGAPoint pZu,PetscReal *K,PetscReal *F,PetscReal *U,PetscReal *HS, PetscReal *Chi,PetscReal *Zu,void *ctx)		//When other fields like Pi or S (etc.) come into this, declare a IGAPoint pPi or pS and a new PescReal *UPi or *US for each
	PetscErrorCode EnergyDensity(IGAPoint p, IGAPoint pChi, IGAPoint pZu, IGAPoint pS, PetscReal *K, PetscReal *F, PetscReal *Chi, PetscReal *Zu, PetscReal *S, void *ctx)		//When other fields like Pi or S (etc.) come into this, declare a IGAPoint pPi or pS and a new PestcReal *UPi or *US for each
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

		//Definition of alternating tensor
		//const PetscReal e[3][3][3]=
		//{
		//	{{0.0,0.0,0.0},{0.0,0.0,1.0},{0.0,-1.0,0.0}},
		//	{{0.0,0.0,-1.0},{0.0,0.0,0.0},{1.0,0.0,0.0}},
		//	{{0.0,1.0,0.0},{-1.0,0.0,0.0},{0.0,0.0,0.0}}
		//};

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

		PetscReal S0[8];																		//Array to contain the vector S0
		IGAPointFormValue(pS,S,&S0[0]);														//Assign chi to its container

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
								FED[a][i]+=0.5*eps*(+fulld2_z[k][l][m]+fulld_Chi[k][l][m]+fullS[k][l][m])*(+fulld2_z[k][l][m]+fulld_Chi[k][l][m]+fullS[k][l][m])*N0[a];
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

//System for L2 projection of elastic energy density
	#undef  __FUNCT__
	#define __FUNCT__ "ElasEnergyDensity"
	//PetscErrorCode Stress(IGAPoint p,IGAPoint pU, IGAPoint pHs,IGAPoint pChi,IGAPoint pZu,PetscReal *K,PetscReal *F,PetscReal *U,PetscReal *HS, PetscReal *Chi,PetscReal *Zu,void *ctx)		//When other fields like Pi or S (etc.) come into this, declare a IGAPoint pPi or pS and a new PescReal *UPi or *US for each
	PetscErrorCode ElasEnergyDensity(IGAPoint p, IGAPoint pChi, IGAPoint pZu, IGAPoint pS, PetscReal *K, PetscReal *F, PetscReal *Chi, PetscReal *Zu, PetscReal *S, void *ctx)		//When other fields like Pi or S (etc.) come into this, declare a IGAPoint pPi or pS and a new PestcReal *UPi or *US for each
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

		//Definition of alternating tensor
		//const PetscReal e[3][3][3]=
		//{
		//	{{0.0,0.0,0.0},{0.0,0.0,1.0},{0.0,-1.0,0.0}},
		//	{{0.0,0.0,-1.0},{0.0,0.0,0.0},{1.0,0.0,0.0}},
		//	{{0.0,1.0,0.0},{-1.0,0.0,0.0},{0.0,0.0,0.0}}
		//};

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
								//FED[a][i]+=0.5*eps*(+fulld2_z[k][l][m+fulld_Chi[k][l][m]+fullS[k][l][m])*(+fulld2_z[k][l][m]+fulld_Chi[k][l][m]+fullS[k][l][m])*N0[a];
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

//System for L2 projection of gradient energy density
	#undef  __FUNCT__
	#define __FUNCT__ "GradEnergyDensity"
	//PetscErrorCode Stress(IGAPoint p,IGAPoint pU, IGAPoint pHs,IGAPoint pChi,IGAPoint pZu,PetscReal *K,PetscReal *F,PetscReal *U,PetscReal *HS, PetscReal *Chi,PetscReal *Zu,void *ctx)		//When other fields like Pi or S (etc.) come into this, declare a IGAPoint pPi or pS and a new PescReal *UPi or *US for each
	PetscErrorCode GradEnergyDensity(IGAPoint p, IGAPoint pChi, IGAPoint pZu, IGAPoint pS, PetscReal *K, PetscReal *F, PetscReal *Chi, PetscReal *Zu, PetscReal *S, void *ctx)		//When other fields like Pi or S (etc.) come into this, declare a IGAPoint pPi or pS and a new PestcReal *UPi or *US for each
	{
		const PetscReal *N0;
		IGAPointGetShapeFuns(p,0,(const PetscReal**)&N0);									//Value of the shape functions
		PetscInt a,b,i,k,l,m,n,nen=p->nen, dof=p->dof;

		PetscReal x[2];																		//Vector of reals, size equal to problem's dimension
		IGAPointFormGeomMap(p,x);

		//Change to consider G=1
		const PetscReal mu=1.0;
		const PetscReal eps=mu/100.0;															//Choose later based on whatever Amit says :)

		PetscReal d_Chi0[4][2];																//Same for its gradient
		IGAPointFormGrad(pChi,Chi,&d_Chi0[0][0]);											//Assign grad chi to its container

		PetscReal d2_Z0[2][2][2];															//Same for its 2nd order partial derivatives
		IGAPointFormHess (pZu,Zu,&d2_Z0[0][0][0]);											//Same for the hessian (second derivatives)

		PetscReal S0[8];																		//Array to contain the vector S0
		IGAPointFormValue(pS,S,&S0[0]);														//Assign chi to its container

		PetscReal fulld_Chi[3][3][3]={0};
		fulld_Chi[0][0][0]=d_Chi0[0][0]; fulld_Chi[0][0][1]=d_Chi0[0][1];
		fulld_Chi[0][1][0]=d_Chi0[1][0]; fulld_Chi[0][1][1]=d_Chi0[1][1];
		fulld_Chi[1][0][0]=d_Chi0[2][0]; fulld_Chi[1][0][1]=d_Chi0[2][1];
		fulld_Chi[1][1][0]=d_Chi0[3][0]; fulld_Chi[1][1][1]=d_Chi0[3][1];

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
						for (k=0;k<3;k++)
					{
						for(l=0;l<3;l++)
						{
							for (m=0;m<3;m++)
							{
								for (n=0;n<3;n++)
								{
									//FED[a][i]+=0.5*((-fulld_z[k][l]-fullChi[k][l])*C[k][l][m][n]*(-fulld_z[m][n]-fullChi[m][n]))*N0[a];
								}
								FED[a][i]+=0.5*eps*(+fulld2_z[k][l][m]+fulld_Chi[k][l][m]+fullS[k][l][m])*(+fulld2_z[k][l][m]+fulld_Chi[k][l][m]+fullS[k][l][m])*N0[a];
							}
						}
					}
				}
			}
		}
		return 0;
	}
//

//System for L2 projection of total stress
	#undef  __FUNCT__
	#define __FUNCT__ "FullS"
	//PetscErrorCode Stress(IGAPoint p,IGAPoint pU, IGAPoint pHs,IGAPoint pChi,IGAPoint pZu,PetscReal *K,PetscReal *F,PetscReal *U,PetscReal *HS, PetscReal *Chi,PetscReal *Zu,void *ctx)		//When other fields like Pi or S (etc.) come into this, declare a IGAPoint pPi or pS and a new PescReal *UPi or *US for each
	PetscErrorCode FullS(IGAPoint p,IGAPoint pChi,IGAPoint pZu,IGAPoint pS, PetscReal *K,PetscReal *F, PetscReal *Chi,PetscReal *Zu,PetscReal *S, void *ctx)		//When other fields like Pi or S (etc.) come into this, declare a IGAPoint pPi or pS and a new PestcReal *UPi or *US for each
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

		PetscReal S0[8];																	//Assign S to a vector
		PetscReal dS[8][2];																	//And its derivative
		IGAPointFormValue(pS,S,&S0[0]);																
		IGAPointFormGrad (pS,S,&dS[0][0]);													//Same for the gradient

		PetscReal full_dS[3][3][3][3]={0};													//Expand grad(S) to full tensor order form, only non-zero elements
		full_dS[0][0][0][0]=dS[0][0]; full_dS[0][0][0][1]=dS[0][1]; 
		full_dS[0][0][1][0]=dS[1][0]; full_dS[0][0][1][1]=dS[1][1]; 
		full_dS[0][1][0][0]=dS[2][0]; full_dS[0][1][0][1]=dS[2][1]; 
		full_dS[0][1][1][0]=dS[3][0]; full_dS[0][1][1][1]=dS[3][1];
		full_dS[1][0][0][0]=dS[4][0]; full_dS[1][0][0][1]=dS[4][1]; 
		full_dS[1][0][1][0]=dS[5][0]; full_dS[1][0][1][1]=dS[5][1]; 
		full_dS[1][1][0][0]=dS[6][0]; full_dS[1][1][0][1]=dS[6][1]; 
		full_dS[1][1][1][0]=dS[7][0]; full_dS[1][1][1][1]=dS[7][1];


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

						Fstress[a][i]+= -0.5*eps*(-fulld3_z[m][n][k][k]-fulld2_Chi[m][n][k][k]-fulld3_z[m][k][n][k]-fulld2_Chi[m][k][n][k])*v[m][n];

						Fstress[a][i]+= -0.5*eps*(-full_dS[m][n][k][k]-full_dS[m][k][n][k])*v[m][n];

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
	PetscErrorCode Valpha(IGAPoint p,IGAPoint pChi,IGAPoint pZu,IGAPoint pAl,PetscReal *K,PetscReal *F, PetscReal *Chi,PetscReal *Zu,PetscReal *U,void *ctx)		//When other fields like Pi or S (etc.) come into this, declare a IGAPoint pPi or pS and a new PestcReal *UPi or *US for each
	{
		const PetscReal *N0;
		IGAPointGetShapeFuns(p,0,(const PetscReal**)&N0);									//Value of the shape functions
		PetscInt a,b,c,d,i,j,k,l,m,nen=p->nen, dof=p->dof;

		//Change to consider G=1
		const PetscReal nu=0.33;
		const PetscReal mu=1.0;
		const PetscReal lambda=2.0*mu*nu/(1.0-2.0*nu);
		const PetscReal eps=mu/100.0;															//Choose later based on whatever Amit says :)

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
		PetscReal Omega=tan(5.0/180.0*ConstPi);
		PetscReal rho=sqrt(x[0]*x[0]+x[1]*x[1]);

		//This is for a single disclination
		g[0]=mu*Omega/(2.0*ConstPi*(1.0-nu))*(log(rho)+(x[1]*x[1])/(rho*rho)+nu/(1.0-2.0*nu));
		g[1]=-mu*Omega/(2.0*ConstPi*(1.0-nu))*x[0]*x[1]/(rho*rho);
		g[2]=-mu*Omega/(2.0*ConstPi*(1.0-nu))*x[0]*x[1]/(rho*rho);
		g[3]=mu*Omega/(2.0*ConstPi*(1.0-nu))*(log(rho)+(x[0]*x[0])/(rho*rho)+nu/(1.0-2.0*nu));

		//This for dislocation with burgers vector in x axis
		//g[0]=-mu*burgers/(ConstPi*(1.0-nu))/2.0*x[1]*(x[1]*x[1]+3*x[0]*x[0])/((x[0]*x[0]+x[1]*x[1])*(x[0]*x[0]+x[1]*x[1]));
		//g[1]=-mu*burgers/(ConstPi*(1.0-nu))/2.0*x[0]*(x[1]*x[1]-x[0]*x[0])/((x[0]*x[0]+x[1]*x[1])*(x[0]*x[0]+x[1]*x[1]));
		//g[2]=-mu*burgers/(ConstPi*(1.0-nu))/2.0*x[0]*(x[1]*x[1]-x[0]*x[0])/((x[0]*x[0]+x[1]*x[1])*(x[0]*x[0]+x[1]*x[1]));
		//g[3]=-mu*burgers/(ConstPi*(1.0-nu))/2.0*x[1]*(x[1]*x[1]-x[0]*x[0])/((x[0]*x[0]+x[1]*x[1])*(x[0]*x[0]+x[1]*x[1]));

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

//System for L2 projection of Grad(Z0)
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

//Helmholtz decomposition of Z curl part
	#undef  __FUNCT__
	#define __FUNCT__ "curlZ"
	PetscErrorCode curlZ(IGAPoint p,IGAPoint pZ, PetscReal *K,PetscReal *F,PetscReal *Z, void *ctx)
	{

		const PetscReal (*N1)[2];
		IGAPointGetShapeFuns(p,1,(const PetscReal**)&N1);

		PetscInt a,b,u,w,i,j,k,nen=p->nen, dof=p->dof;

		PetscReal dz[4][2];																				//Create array to recieve Alfa
		IGAPointFormGrad (pZ,Z,&dz[0][0]);																//This fills the values of the derivatives
		
		PetscReal fulld_Z[3][3][3]={0};
		fulld_Z[0][0][0]=dz[0][0]; fulld_Z[0][0][1]=dz[0][1];											//Expand dZ to full 3rd order form, only non-zero elements
		fulld_Z[0][1][0]=dz[1][0]; fulld_Z[0][1][1]=dz[1][1];
		fulld_Z[1][0][0]=dz[2][0]; fulld_Z[1][0][1]=dz[2][1];
		fulld_Z[1][1][0]=dz[3][0]; fulld_Z[1][1][1]=dz[3][1];

		PetscReal (*KchiZ)[dof][nen][dof] = (typeof(KchiZ)) K;
		PetscReal (*FchiZ)[dof] = (PetscReal (*)[dof]) F;

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

					if(u==0){
						dv[0][0][0]=Na_x; dv[0][0][1]=Na_y;}
					if(u==1){
						dv[0][1][0]=Na_x; dv[0][1][1]=Na_y;}
					if(u==2){
						dv[1][0][0]=Na_x; dv[1][0][1]=Na_y;}
					if(u==3){
						dv[1][1][0]=Na_x; dv[1][1][1]=Na_y;}

					for (b=0; b<nen; b++)
					{
						PetscReal Nb_x = N1[b][0];
						PetscReal Nb_y = N1[b][1];
						for (w=0; w<dof; w++)
						{
							PetscReal dChiZ[3][3][3]={0};

							if(w==0){
								dChiZ[0][0][0]=Nb_x; dChiZ[0][0][1]=Nb_y;}
							if(w==1){
								dChiZ[0][1][0]=Nb_x; dChiZ[0][1][1]=Nb_y;}
							if(w==2){
								dChiZ[1][0][0]=Nb_x; dChiZ[1][0][1]=Nb_y;}
							if(w==3){
								dChiZ[1][1][0]=Nb_x; dChiZ[1][1][1]=Nb_y;}

							KchiZ[a][i][b][j]=0.0;
							for(i=0;i<3;i++)
							{
								for(j=0;j<3;j++)
								{
									for(k=0;k<3;k++)
									{
											KchiZ[a][u][b][w]+=dChiZ[i][j][k]*dv[i][j][k]-dChiZ[i][j][k]*dv[i][k][j]+dChiZ[i][j][j]*dv[i][k][k];
									}
								}
							}
						}
					}

					FchiZ[a][u]=0.0;
					for(i=0;i<3;i++)
					{
						for(j=0;j<3;j++)
						{
							for(k=0;k<3;k++)
							{
								FchiZ[a][u]+=fulld_Z[i][j][k]*dv[i][j][k]-fulld_Z[i][j][k]*dv[i][k][j];
							}
						}
					}
				}
			}
		}
		return 0;
	}
//

//Helmholtz decomposition of Z, grad part
	#undef  __FUNCT__
	#define __FUNCT__ "gradZa"
	PetscErrorCode gradZa(IGAPoint p,IGAPoint pZ,PetscReal *K,PetscReal *F,PetscReal *Z,void *ctx)
	{

		const PetscReal *N0,(*N1)[2];
		IGAPointGetShapeFuns(p,0,(const PetscReal**)&N0);
		IGAPointGetShapeFuns(p,1,(const PetscReal**)&N1);

		PetscInt a,b,u,w,i,j,nen=p->nen, dof=p->dof;

		PetscReal dZ[4][2];																		//Create array to recieve dAlfa
		IGAPointFormGrad (pZ,Z,&dZ[0][0]);														//This fills the values of the derivatives

		PetscReal fulldZ[3][3][3]={0};
		//Same for gradS
		fulldZ[0][0][0]=dZ[0][0]; fulldZ[0][0][1]=dZ[0][1];
		fulldZ[0][1][0]=dZ[1][0]; fulldZ[0][1][1]=dZ[1][1];
		fulldZ[1][0][0]=dZ[2][0]; fulldZ[1][0][1]=dZ[2][1];
		fulldZ[1][1][0]=dZ[3][0]; fulldZ[1][1][1]=dZ[3][1];

		PetscReal (*KZa)[dof][nen][dof] = (typeof(KZa)) K;
		PetscReal (*FZa)[dof] = (PetscReal (*)[dof]) F;

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
					PetscReal v[3]={0};
					PetscReal dv[3][3]={0};
					if(u==0){
						v[0]=Na;
						dv[0][0]=Na_x; dv[0][1]=Na_y;}
					if(u==1){
						v[1]=Na;
						dv[1][0]=Na_x; dv[1][1]=Na_y;}

					for (b=0; b<nen; b++)
					{
						PetscReal Nb_x = N1[b][0];
						PetscReal Nb_y = N1[b][1];
						for (w=0; w<dof; w++)
						{
							PetscReal dZa[3][3]={0};
							if(w==0){
								dZa[0][0]=Nb_x; dZa[0][1]=Nb_y;}
							if(w==1){
								dZa[1][0]=Nb_x; dZa[1][1]=Nb_y;}

							KZa[a][u][b][w]=0.0;
							for(i=0;i<3;i++)
							{
								for(j=0;j<3;j++)
								{
									KZa[a][u][b][w]+=dZa[i][j]*dv[i][j];
								}
							}
						}
					}

					FZa[a][u]=0.0;
					for(i=0;i<3;i++)
					{
						for(j=0;j<3;j++)
						{
							FZa[a][u]+=-fulldZ[i][j][j]*v[i];
						}
					}
				}
			}
		}

		return 0;
	}
//

//System for L2 projection of grad(grad(z))
	#undef  __FUNCT__
	#define __FUNCT__ "SysGradGradz"
	//PetscErrorCode Stress(IGAPoint p,IGAPoint pU, IGAPoint pHs,IGAPoint pChi,IGAPoint pZu,PetscReal *K,PetscReal *F,PetscReal *U,PetscReal *HS, PetscReal *Chi,PetscReal *Zu,void *ctx)		//When other fields like Pi or S (etc.) come into this, declare a IGAPoint pPi or pS and a new PescReal *UPi or *US for each
	PetscErrorCode SysGradGradz(IGAPoint p,IGAPoint pZu,PetscReal *K,PetscReal *F,PetscReal *Zu,void *ctx)		//When other fields like Pi or S (etc.) come into this, declare a IGAPoint pPi or pS and a new PestcReal *UPi or *US for each
	{
		const PetscReal *N0;
		IGAPointGetShapeFuns(p,0,(const PetscReal**)&N0);									//Value of the shape functions
		//IGAPointGetShapeFuns(p,1,(const PetscReal**)&N1);									//Derivatives of the shape functions
		//IGAPointGetShapeFuns(p,2,(const PetscReal**)&N2);									//Second derivatives of the shape funcions
		PetscInt a,b,i,k,m,n,nen=p->nen, dof=p->dof;

		PetscReal d2_Z0[2][2][2];															//Same for its 2nd order partial derivatives
		IGAPointFormHess (pZu,Zu,&d2_Z0[0][0][0]);											//Same for the hessian (second derivatives)

		//Expanding z (and derivatives) to 3 components, more convenient for sums in for loops
		PetscReal fulld2_z[3][3][3]={0};
		fulld2_z[0][0][0]=d2_Z0[0][0][0]; fulld2_z[0][0][1]=d2_Z0[0][0][1];
		fulld2_z[0][1][0]=d2_Z0[0][1][0]; fulld2_z[0][1][1]=d2_Z0[0][1][1];
		fulld2_z[1][0][0]=d2_Z0[1][0][0]; fulld2_z[1][0][1]=d2_Z0[1][0][1];
		fulld2_z[1][1][0]=d2_Z0[1][1][0]; fulld2_z[1][1][1]=d2_Z0[1][1][1];

		PetscReal (*Kj)[dof][nen][dof] = (typeof(Kj)) K;
		PetscReal (*Fj)[dof] = (PetscReal (*)[dof])F;

		PetscReal v[3][3][3]={0};

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
						Kj[a][i][b][i]=N0[a]*N0[b];
					}
				}
			}

			for(a=0 ;a<nen; a++)
			{
				for (i=0; i<dof; i++)
				{
					if (i==0)
					{
						for(m=0;m<3;m++)
						{
							for (n=0;n<3;n++)
							{
								for (k=0;k<3;k++)
								{
									v[m][n][k]=0.0;
								}
							}
						}
						m=0; n=0; k=0;
						v[m][n][k]=N0[a];
					}
					else if (i==1)
					{
						for(m=0;m<3;m++)
						{
							for (n=0;n<3;n++)
							{
								for (k=0;k<3;k++)
								{
									v[m][n][k]=0.0;
								}
							}
						}
						m=0; n=0; k=1;
						v[m][n][k]=N0[a];
					}
					else if (i==2)
					{
						for(m=0;m<3;m++)
						{
							for (n=0;n<3;n++)
							{
								for (k=0;k<3;k++)
								{
									v[m][n][k]=0.0;
								}
							}
						}
						m=0; n=1; k=0;
						v[m][n][k]=N0[a];
					}
					else if (i==3)
					{
						for(m=0;m<3;m++)
						{
							for (n=0;n<3;n++)
							{
								for (k=0;k<3;k++)
								{
									v[m][n][k]=0.0;
								}
							}
						}
						m=0; n=1; k=1;
						v[m][n][k]=N0[a];
					}
					else if (i==4)
					{
						for(m=0;m<3;m++)
						{
							for (n=0;n<3;n++)
							{
								for (k=0;k<3;k++)
								{
									v[m][n][k]=0.0;
								}
							}
						}
						m=1; n=0; k=0;
						v[m][n][k]=N0[a];
					}
					else if (i==5)
					{
						for(m=0;m<3;m++)
						{
							for (n=0;n<3;n++)
							{
								for (k=0;k<3;k++)
								{
									v[m][n][k]=0.0;
								}
							}
						}
						m=1; n=0; k=1;
						v[m][n][k]=N0[a];
					}
					else if (i==6)
					{
						for(m=0;m<3;m++)
						{
							for (n=0;n<3;n++)
							{
								for (k=0;k<3;k++)
								{
									v[m][n][k]=0.0;
								}
							}
						}
						m=1; n=1; k=0;
						v[m][n][k]=N0[a];
					}
					else if (i==7)
					{
						for(m=0;m<3;m++)
						{
							for (n=0;n<3;n++)
							{
								for (k=0;k<3;k++)
								{
									v[m][n][k]=0.0;
								}
							}
						}
						m=1; n=1; k=1;
						v[m][n][k]=N0[a];
					}
				
					Fj[a][i]=(-fulld2_z[m][n][k])*v[m][n][k];
					
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
		source = fopen("/home/eazegpi/CodigosPetIga/PruebaV5S.c","r");
		dest   = fopen("/home/eazegpi/Results/PruebaV5S.c","w");

		while (0 < (bytes = fread(buffer, 1, sizeof(buffer), source)))
			fwrite(buffer, 1, bytes, dest);

		fclose(source);
		fclose(dest);
	}
//

/*
//Creation of types and systems for the L2 projection of S0
	PetscPrintf(PETSC_COMM_WORLD,"\n System for L2 projection for S starting \n\n");
	IGA igaS;
	ierr = IGACreate(PETSC_COMM_WORLD,&igaS);CHKERRQ(ierr);
	ierr = IGASetDim(igaS,2);CHKERRQ(ierr);														//Spatial dimension of the problem
	ierr = IGASetDof(igaS,8);CHKERRQ(ierr);														//Number of degrees of freedom, per node
	ierr = IGASetOrder(igaS,1);CHKERRQ(ierr);													//Number of spatial derivatives to calculate
	ierr = IGASetFromOptions(igaS);CHKERRQ(ierr);												//Note: The order (or degree) of the shape functions is given by the mesh!
	ierr = IGARead(igaS,"./geometry.dat");CHKERRQ(ierr);
	ierr = IGASetUp(igaS);CHKERRQ(ierr);
	
	for (dir=0; dir<2; dir++) 
	{
		for (side=0; side<2; side++) 
		{
			//ierr = IGASetBoundaryValue(iga,dir,side,dof,0.0);CHKERRQ(ierr);    				// Dirichlet boundary conditions
			ierr = IGASetBoundaryForm(igaS,dir,side,PETSC_TRUE);CHKERRQ(ierr);  				// Neumann boundary conditions
		}
	}

	Vec s0;
	Mat Kl2S;
	Vec Fl2S;
	ierr = IGACreateVec(igaS,&s0);CHKERRQ(ierr);  
	ierr = IGACreateMat(igaS,&Kl2S);CHKERRQ(ierr);
	ierr = IGACreateVec(igaS,&Fl2S);CHKERRQ(ierr);
	ierr = IGASetFormSystem(igaS,L2ProjectionS,&userL2);CHKERRQ(ierr);
	ierr = IGAComputeSystem(igaS,Kl2S,Fl2S);CHKERRQ(ierr);

	//This parts set and calls KSP to solve the linear system
	KSP kspl2S;
	ierr = IGACreateKSP(igaS,&kspl2S);CHKERRQ(ierr);										
	ierr = KSPSetOperators(kspl2S,Kl2S,Kl2S);CHKERRQ(ierr); 								//This function creates the matrix for the system on the second parameter and uses the 3rd parameter as a preconditioner
	ierr = KSPSetType(kspl2S,KSPCG);CHKERRQ(ierr);											//Using KSPCG (conjugated gradient) because the matrix is symmetric
	//ierr = KSPSetOptionsPrefix(kspl2S,"l2pS_");CHKERRQ(ierr);
	ierr = KSPSetFromOptions(kspl2S);CHKERRQ(ierr);
	ierr = KSPSetTolerances(kspl2S,1.0e-18,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);CHKERRQ(ierr);
	ierr = KSPSolve(kspl2S,Fl2S,s0);CHKERRQ(ierr);										//This is a simple system, so it can be solved with just this command

	ierr = KSPDestroy(&kspl2S);CHKERRQ(ierr);
	ierr = MatDestroy(&Kl2S);CHKERRQ(ierr);
	ierr = VecDestroy(&Fl2S);CHKERRQ(ierr);
	char nameS[]="/S-2d-0.dat";
	char pathS[512];
	sprintf(pathS,"%s%s",direct,nameS);
	ierr = IGAWriteVec(igaS,s0,pathS);CHKERRQ(ierr);
//
*/

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

	char nameS[]="/Input-S-2d-0.dat";
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
	//PC pcchiS;
	//ierr = KSPGetPC(kspchiS,&pcchiS); CHKERRQ(ierr);
	//ierr = PCSetType(pcchiS,PCLU); CHKERRQ(ierr);
	//ierr = PCFactorSetMatSolverType(pcchiS,MATSOLVERMUMPS); CHKERRQ(ierr);
	ierr = KSPSetFromOptions(kspchiS);CHKERRQ(ierr);
	ierr = KSPSetTolerances(kspchiS,1.0e-12,5.0e-24,PETSC_DEFAULT,PETSC_DEFAULT);CHKERRQ(ierr);
	ierr = KSPSolve(kspchiS,FchiS,chiS0);CHKERRQ(ierr);

	ierr = KSPDestroy(&kspchiS);CHKERRQ(ierr);
	ierr = MatDestroy(&KchiS);CHKERRQ(ierr);
	ierr = VecDestroy(&FchiS);CHKERRQ(ierr);
	
	char namechiS[]="/ChiS-2d-0.dat";
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
					ierr = L2ProjectionAlSp(pointAlp,pointSp,KpointAlp,FpointAlp,S0Alp,&userL2);CHKERRQ(ierr);
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
	
	char nameAlp[]="/Alp-2d-0.dat";
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

	char namechiUp[]="/ChiUp-2d-0.dat";
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
	char nameZ0[]="/Z0-2d-0.dat";
	char pathZ0[512];
	sprintf(pathZ0,"%s%s",direct,nameZ0);
	ierr = IGAWriteVec(igaZ0,Z0,pathZ0);CHKERRQ(ierr);	
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
	const PetscReal	*arrayZ0Stress,*arrayChi0Stress, *arraySStress;			//arrayU
	Vec				localZ0Stress,localChi0Stress, localSStress;				//localU
	PetscReal		*Chi0Stress,*Z0Stress, *SStress;											//U

  	IGAFormSystem	wtfStress;
 	void			*wtf2Stress;

 	KSP kspStress;
	ierr = IGACreateKSP(igaStress,&kspStress);CHKERRQ(ierr);

	// Get local vectors Z0 and Chi0 and arrays
	ierr = IGAGetLocalVecArray(igachiUp,chiUp0,&localChi0Stress,&arrayChi0Stress);CHKERRQ(ierr);
	ierr = IGAGetLocalVecArray(igaZ0,Z0,&localZ0Stress,&arrayZ0Stress);CHKERRQ(ierr);
	ierr = IGAGetLocalVecArray(igaS,s0,&localSStress,&arraySStress);CHKERRQ(ierr);

	// Element loop
	ierr = IGABeginElement(igaStress,&elemStress);CHKERRQ(ierr);
	ierr = IGABeginElement(igaZ0,&elemZ0);CHKERRQ(ierr);
	ierr = IGABeginElement(igachiUp,&elemchiUp);CHKERRQ(ierr);
	ierr = IGABeginElement(igaS,&elemS);CHKERRQ(ierr);

	while (IGANextElement(igaStress,elemStress)) 
	{
		IGANextElement(igaZ0,elemZ0);
		IGANextElement(igachiUp,elemchiUp);
		IGANextElement(igaS,elemS);

		ierr = IGAElementGetWorkMat(elemStress,&KlocStress);CHKERRQ(ierr);
		ierr = IGAElementGetWorkVec(elemStress,&FlocStress);CHKERRQ(ierr);
		ierr = IGAElementGetValues(elemZ0,arrayZ0Stress,&Z0Stress);CHKERRQ(ierr);
		ierr = IGAElementGetValues(elemchiUp,arrayChi0Stress,&Chi0Stress);CHKERRQ(ierr);
		ierr = IGAElementGetValues(elemS,arraySStress,&SStress);CHKERRQ(ierr);

		// FormSystem loop
		while (IGAElementNextFormSystem(elemStress,&wtfStress,&wtf2Stress)) 
		{
		// Quadrature loop
			ierr = IGAElementBeginPoint(elemStress,&pointStress);CHKERRQ(ierr);
			ierr = IGAElementBeginPoint(elemZ0,&pointZ0);CHKERRQ(ierr);
			ierr = IGAElementBeginPoint(elemchiUp,&pointchiUp);CHKERRQ(ierr);
			ierr = IGAElementBeginPoint(elemS,&pointS);CHKERRQ(ierr);

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

					ierr = IGAPointGetWorkMat(pointStress,&KpointStress);CHKERRQ(ierr);
					ierr = IGAPointGetWorkVec(pointStress,&FpointStress);CHKERRQ(ierr);
					//	   Stress(IGAPoint p,IGAPoint pChi,IGAPoint pZu,IGAPoint pS, PetscReal *K,PetscReal *F, PetscReal *Chi,PetscReal *Zu,PetscReal *S, void *ctx)
					ierr = Stress(pointStress,pointchiUp,pointZ0,pointS,KpointStress,FpointStress,Chi0Stress,Z0Stress,SStress,NULL);CHKERRQ(ierr);
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

			ierr = IGAElementEndPoint(elemStress,&pointStress);CHKERRQ(ierr);
			ierr = IGAElementEndPoint(elemZ0,&pointZ0);CHKERRQ(ierr);
			ierr = IGAElementEndPoint(elemchiUp,&pointchiUp);CHKERRQ(ierr);
			ierr = IGAElementEndPoint(elemS,&pointS);CHKERRQ(ierr);
		}

		//ierr = IGAElementFixSystem(elemStress,KlocStress,FlocStress);CHKERRQ(ierr);					//This sets Dirichlet condition ¿? (Yes, this applies the conditions from IGASetBoundaryValue)
		ierr = IGAElementAssembleMat(elemStress,KlocStress,KStress);CHKERRQ(ierr);
		ierr = IGAElementAssembleVec(elemStress,FlocStress,FStress);CHKERRQ(ierr);
	}
	IGANextElement(igaZ0,elemZ0);
	IGANextElement(igachiUp,elemchiUp);
	IGANextElement(igaS,elemS);

	ierr = IGAEndElement(igaStress,&elemStress);CHKERRQ(ierr);
	ierr = IGAEndElement(igaZ0,&elemZ0);CHKERRQ(ierr);
	ierr = IGAEndElement(igachiUp,&elemchiUp);CHKERRQ(ierr);
	ierr = IGAEndElement(igaS,&elemS);CHKERRQ(ierr);

	// Restore local vectors u, Z0, Chi0 and arrays
	ierr = IGARestoreLocalVecArray(igaZ0,Z0,&localZ0Stress,&arrayZ0Stress);CHKERRQ(ierr);
	ierr = IGARestoreLocalVecArray(igachiUp,chiUp0,&localChi0Stress,&arrayChi0Stress);CHKERRQ(ierr);
	ierr = IGARestoreLocalVecArray(igaS,s0,&localSStress,&arraySStress);CHKERRQ(ierr);

	ierr = MatAssemblyBegin(KStress,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyEnd  (KStress,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

	ierr = VecAssemblyBegin(FStress);CHKERRQ(ierr);
	ierr = VecAssemblyEnd  (FStress);CHKERRQ(ierr);

	ierr = KSPSetOperators(kspStress,KStress,KStress);CHKERRQ(ierr);
	//PC pcStress;
	//ierr = KSPGetPC(kspStress,&pcStress); CHKERRQ(ierr);
	//ierr = PCSetType(pcStress,PCLU); CHKERRQ(ierr);
	//ierr = PCFactorSetMatSolverType(pcStress,MATSOLVERMUMPS); CHKERRQ(ierr);
	ierr = KSPSetFromOptions(kspStress);CHKERRQ(ierr);
	ierr = KSPSetTolerances(kspStress,1e-24,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);CHKERRQ(ierr);
	ierr = KSPSolve(kspStress,FStress,sigma0);CHKERRQ(ierr);

	ierr = KSPDestroy(&kspStress);CHKERRQ(ierr);
	ierr = MatDestroy(&KStress);CHKERRQ(ierr);
	ierr = VecDestroy(&FStress);CHKERRQ(ierr);

	char nameStress[]="/sigma-2d-0.dat";
	char pathStress[512];
	sprintf(pathStress,"%s%s",direct,nameStress);
	ierr = IGAWriteVec(igaStress,sigma0,pathStress);CHKERRQ(ierr);	
//

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
	//PC pcClassicStress;
	//ierr = KSPGetPC(kspClassicStress,&pcClassicStress); CHKERRQ(ierr);
	//ierr = PCSetType(pcClassicStress,PCLU); CHKERRQ(ierr);
	//ierr = PCFactorSetMatSolverType(pcClassicStress,MATSOLVERMUMPS); CHKERRQ(ierr);
	ierr = KSPSetFromOptions(kspClassicStress);CHKERRQ(ierr);
	ierr = KSPSetTolerances(kspClassicStress,1.0e-24,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);CHKERRQ(ierr);
	ierr = KSPSolve(kspClassicStress,FClassicStress,classicSigma0);CHKERRQ(ierr);

	ierr = KSPDestroy(&kspClassicStress);CHKERRQ(ierr);
	ierr = MatDestroy(&KClassicStress);CHKERRQ(ierr);
	ierr = VecDestroy(&FClassicStress);CHKERRQ(ierr);

	char nameClassicStress[]="/classicSigma-2d-0.dat";
	char pathClassicStress[512];
	sprintf(pathClassicStress,"%s%s",direct,nameClassicStress);
	ierr = IGAWriteVec(igaClassicStress,classicSigma0,pathClassicStress);CHKERRQ(ierr);	
//

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
	const PetscReal	*arrayZ0CS,*arrayChi0CS,*arrayS0CS;								//arrayU
	Vec				localZ0CS,localChi0CS,localS0CS;								//localU
	PetscReal		*Chi0CS,*Z0CS,*S0CS;											//U

  	IGAFormSystem	wtfCS;
 	void			*wtf2CS;

 	KSP kspCS;
	ierr = IGACreateKSP(igaCS,&kspCS);CHKERRQ(ierr);

	// Get local vectors Z0 and Chi0 and arrays
	ierr = IGAGetLocalVecArray(igachiUp,chiUp0,&localChi0CS,&arrayChi0CS);CHKERRQ(ierr);
	ierr = IGAGetLocalVecArray(igaZ0,Z0,&localZ0CS,&arrayZ0CS);CHKERRQ(ierr);
	ierr = IGAGetLocalVecArray(igaS,s0,&localS0CS,&arrayS0CS);CHKERRQ(ierr);

	// Element loop
	ierr = IGABeginElement(igaCS,&elemCS);CHKERRQ(ierr);
	ierr = IGABeginElement(igaZ0,&elemZ0);CHKERRQ(ierr);
	ierr = IGABeginElement(igachiUp,&elemchiUp);CHKERRQ(ierr);
	ierr = IGABeginElement(igaS,&elemS);CHKERRQ(ierr);

	while (IGANextElement(igaCS,elemCS)) 				
	{
		IGANextElement(igaZ0,elemZ0);
		IGANextElement(igachiUp,elemchiUp);
		IGANextElement(igaS,elemS);

		ierr = IGAElementGetWorkMat(elemCS,&KlocCS);CHKERRQ(ierr);
		ierr = IGAElementGetWorkVec(elemCS,&FlocCS);CHKERRQ(ierr);
		ierr = IGAElementGetValues(elemZ0,arrayZ0CS,&Z0CS);CHKERRQ(ierr);
		ierr = IGAElementGetValues(elemchiUp,arrayChi0CS,&Chi0CS);CHKERRQ(ierr);
		ierr = IGAElementGetValues(elemS,arrayS0CS,&S0CS);CHKERRQ(ierr);

		// FormSystem loop
		while (IGAElementNextFormSystem(elemCS,&wtfCS,&wtf2CS))
		{
		// Quadrature loop
			ierr = IGAElementBeginPoint(elemCS,&pointCS);CHKERRQ(ierr);
			ierr = IGAElementBeginPoint(elemZ0,&pointZ0);CHKERRQ(ierr);
			ierr = IGAElementBeginPoint(elemchiUp,&pointchiUp);CHKERRQ(ierr);
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

				if(pointCS->atboundary==0 && pointZ0->atboundary==0 && pointchiUp->atboundary==0)
				{
					IGAElementNextPoint(elemZ0,pointZ0);
					IGAElementNextPoint(elemchiUp,pointchiUp);
					IGAElementNextPoint(elemS,pointS);
					ierr = IGAPointGetWorkMat(pointCS,&KpointCS);CHKERRQ(ierr);
					ierr = IGAPointGetWorkVec(pointCS,&FpointCS);CHKERRQ(ierr);
					//	   CoupleStress(IGAPoint p, IGAPoint pChi, IGAPoint pZu,IGAPoint pS, PetscReal *K, PetscReal *F, PetscReal *Chi,PetscReal *Zu,PetscReal *s0, void *ctx)
					ierr = CoupleStress(pointCS,pointchiUp,pointZ0,pointS,KpointCS,FpointCS,Chi0CS,Z0CS,S0CS,NULL);CHKERRQ(ierr);
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
			while (pointS->index != -1)
			{
				IGAElementNextPoint(elemS,pointS);
			}
			ierr = IGAElementEndPoint(elemCS,&pointCS);CHKERRQ(ierr);
			ierr = IGAElementEndPoint(elemZ0,&pointZ0);CHKERRQ(ierr);
			ierr = IGAElementEndPoint(elemchiUp,&pointchiUp);CHKERRQ(ierr);
			ierr = IGAElementEndPoint(elemS,&pointS);CHKERRQ(ierr);
		}

		//ierr = IGAElementFixSystem(elemStress,KlocStress,FlocStress);CHKERRQ(ierr);					//This sets Dirichlet condition ¿? (Yes, this applies the conditions from IGASetBoundaryValue)
		ierr = IGAElementAssembleMat(elemCS,KlocCS,KCStress);CHKERRQ(ierr);
		ierr = IGAElementAssembleVec(elemCS,FlocCS,FCStress);CHKERRQ(ierr);
	}
	IGANextElement(igaZ0,elemZ0);
	IGANextElement(igachiUp,elemchiUp);
	IGANextElement(igaS,elemS);

	ierr = IGAEndElement(igaCS,&elemCS);CHKERRQ(ierr);
	ierr = IGAEndElement(igaZ0,&elemZ0);CHKERRQ(ierr);
	ierr = IGAEndElement(igachiUp,&elemchiUp);CHKERRQ(ierr);													//////////Voy aquí
	ierr = IGAEndElement(igaS,&elemS);CHKERRQ(ierr);

	// Restore local vectors u, Z0, Chi0, S and arrays
	ierr = IGARestoreLocalVecArray(igaZ0,Z0,&localZ0CS,&arrayZ0CS);CHKERRQ(ierr);
	ierr = IGARestoreLocalVecArray(igachiUp,chiUp0,&localChi0CS,&arrayChi0CS);CHKERRQ(ierr);
	ierr = IGARestoreLocalVecArray(igaS,s0,&localS0CS,&arrayS0CS);CHKERRQ(ierr);
	ierr = MatAssemblyBegin(KCStress,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyEnd  (KCStress,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

	ierr = VecAssemblyBegin(FCStress);CHKERRQ(ierr);
	ierr = VecAssemblyEnd  (FCStress);CHKERRQ(ierr);

	ierr = KSPSetOperators(kspCS,KCStress,KCStress);CHKERRQ(ierr);
	//PC pcCS;
	//ierr = KSPGetPC(kspCS,&pcCS); CHKERRQ(ierr);
	//ierr = PCSetType(pcCS,PCLU); CHKERRQ(ierr);
	//ierr = PCFactorSetMatSolverType(pcCS,MATSOLVERMUMPS); CHKERRQ(ierr);
	ierr = KSPSetFromOptions(kspCS);CHKERRQ(ierr);
	ierr = KSPSetTolerances(kspCS,1e-24,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);CHKERRQ(ierr);
	ierr = KSPSolve(kspCS,FCStress,lambda0);CHKERRQ(ierr);

	ierr = KSPDestroy(&kspStress);CHKERRQ(ierr);
	ierr = MatDestroy(&KStress);CHKERRQ(ierr);
	ierr = VecDestroy(&FStress);CHKERRQ(ierr);

	char nameCStress[]="/lambda-2d-0.dat";
	char pathCStress[512];
	sprintf(pathCStress,"%s%s",direct,nameCStress);
	ierr = IGAWriteVec(igaCS,lambda0,pathCStress);CHKERRQ(ierr);	
//

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
	const PetscReal	*arrayChi0ED,*arrayZ0ED,*arrayS0ED;									//arrayU
	Vec				localChi0ED,localZ0ED,localS0ED;									//localU
	PetscReal		*chi0ED,*Z0ED,*S0ED;												//U

  	IGAFormSystem	wtfED;
 	void			*wtf2ED;

 	KSP kspED;
	ierr = IGACreateKSP(igaED,&kspED);CHKERRQ(ierr);

	// Get local vectors Z0 and Chi0 and arrays
	ierr = IGAGetLocalVecArray(igachiUp,chiUp0,&localChi0ED,&arrayChi0ED);CHKERRQ(ierr);
	ierr = IGAGetLocalVecArray(igaZ0,Z0,&localZ0ED,&arrayZ0ED);CHKERRQ(ierr);
	ierr = IGAGetLocalVecArray(igaS,s0,&localS0ED,&arrayS0ED);CHKERRQ(ierr);

	// Element loop
	ierr = IGABeginElement(igaED,&elemED);CHKERRQ(ierr);
	ierr = IGABeginElement(igaZ0,&elemZ0);CHKERRQ(ierr);
	ierr = IGABeginElement(igachiUp,&elemchiUp);CHKERRQ(ierr);
	ierr = IGABeginElement(igaS,&elemS);CHKERRQ(ierr);

	while (IGANextElement(igaED,elemED)) 				
	{
		IGANextElement(igachiUp,elemchiUp);
		IGANextElement(igaZ0,elemZ0);
		IGANextElement(igaS,elemS);

		ierr = IGAElementGetWorkMat(elemED,&KlocED);CHKERRQ(ierr);
		ierr = IGAElementGetWorkVec(elemED,&FlocED);CHKERRQ(ierr);
		ierr = IGAElementGetValues(elemchiUp,arrayChi0ED,&chi0ED);CHKERRQ(ierr);
		ierr = IGAElementGetValues(elemZ0,arrayZ0ED,&Z0ED);CHKERRQ(ierr);
		ierr = IGAElementGetValues(elemS,arrayS0ED,&S0ED);CHKERRQ(ierr);

		// FormSystem loop
		while (IGAElementNextFormSystem(elemED,&wtfED,&wtf2ED))
		{
		// Quadrature loop
			ierr = IGAElementBeginPoint(elemED,&pointED);CHKERRQ(ierr);
			ierr = IGAElementBeginPoint(elemZ0,&pointZ0);CHKERRQ(ierr);
			ierr = IGAElementBeginPoint(elemchiUp,&pointchiUp);CHKERRQ(ierr);
			ierr = IGAElementBeginPoint(elemS,&pointS);CHKERRQ(ierr);

			while (IGAElementNextPoint(elemED,pointED))
			{
				if(pointED->atboundary==0 && pointchiUp->atboundary==0 && pointZ0->atboundary==0)
				{
					IGAElementNextPoint(elemZ0,pointZ0);
					IGAElementNextPoint(elemchiUp,pointchiUp);
					IGAElementNextPoint(elemS,pointS);

					ierr = IGAPointGetWorkMat(pointED,&KpointED);CHKERRQ(ierr);
					ierr = IGAPointGetWorkVec(pointED,&FpointED);CHKERRQ(ierr);
					//	   EnergyDensity(IGAPoint p, IGAPoint pChi, IGAPoint pZu, IGAPoint pS, PetscReal *K, PetscReal *F, PetscReal *Chi, PetscReal *Zu, PetscReal *S0, void *ctx)	
					ierr = EnergyDensity(pointED,pointchiUp,pointZ0,pointS,KpointED,FpointED,chi0ED,Z0ED,S0ED,NULL);CHKERRQ(ierr);
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
			ierr = IGAElementEndPoint(elemED,&pointED);CHKERRQ(ierr);
			ierr = IGAElementEndPoint(elemchiUp,&pointchiUp);CHKERRQ(ierr);
			ierr = IGAElementEndPoint(elemZ0,&pointZ0);CHKERRQ(ierr);
			ierr = IGAElementEndPoint(elemS,&pointS);CHKERRQ(ierr);
		}

		//ierr = IGAElementFixSystem(elemStress,KlocStress,FlocStress);CHKERRQ(ierr);					//This sets Dirichlet condition ¿? (Yes, this applies the conditions from IGASetBoundaryValue)
		ierr = IGAElementAssembleMat(elemED,KlocED,KED);CHKERRQ(ierr);
		ierr = IGAElementAssembleVec(elemED,FlocED,FED);CHKERRQ(ierr);
	}
	IGANextElement(igachiUp,elemchiUp);
	IGANextElement(igaZ0,elemZ0);
	IGANextElement(igaS,elemS);

	ierr = IGAEndElement(igaED,&elemED);CHKERRQ(ierr);
	ierr = IGAEndElement(igaZ0,&elemZ0);CHKERRQ(ierr);
	ierr = IGAEndElement(igachiUp,&elemchiUp);CHKERRQ(ierr);
	ierr = IGAEndElement(igaS,&elemS);CHKERRQ(ierr);

	// Restore local vectors u, Z0, Chi0 and arrays
	ierr = IGARestoreLocalVecArray(igachiUp,chiUp0,&localChi0ED,&arrayChi0ED);CHKERRQ(ierr);
	ierr = IGARestoreLocalVecArray(igaZ0,Z0,&localZ0ED,&arrayZ0ED);CHKERRQ(ierr);
	ierr = IGARestoreLocalVecArray(igaS,s0,&localS0ED,&arrayS0ED);CHKERRQ(ierr);

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
	ierr = KSPSetTolerances(kspED,1e-24,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);CHKERRQ(ierr);
	ierr = KSPSolve(kspED,FED,ed0);CHKERRQ(ierr);

	ierr = KSPDestroy(&kspED);CHKERRQ(ierr);
	ierr = MatDestroy(&KED);CHKERRQ(ierr);
	ierr = VecDestroy(&FED);CHKERRQ(ierr);

	char nameED[]="/DensityEnergy.dat";
	char pathED[512];
	sprintf(pathED,"%s%s",direct,nameED);
	ierr = IGAWriteVec(igaED,ed0,pathED);CHKERRQ(ierr);	
//

/*
//System for L2 projection of Elastic Energy Density
	PetscPrintf(PETSC_COMM_WORLD,"\n System for Elastic Energy Density starting \n\n");
	IGA igaEL;
	ierr = IGACreate(PETSC_COMM_WORLD,&igaEL);CHKERRQ(ierr);
	ierr = IGASetDim(igaEL,2);CHKERRQ(ierr);													//Spatial dimension of the problem
	ierr = IGASetDof(igaEL,1);CHKERRQ(ierr);													//Number of degrees of freedom, per node
	ierr = IGASetOrder(igaEL,2);CHKERRQ(ierr);												//Number of spatial derivatives to calculate
	ierr = IGASetFromOptions(igaEL);CHKERRQ(ierr);											//Note: The order (or degree) of the shape functions is given by the mesh!
	ierr = IGARead(igaEL,"./geometry3.dat");CHKERRQ(ierr);
	ierr = IGASetUp(igaEL);CHKERRQ(ierr);

	for (dir=0; dir<2; dir++) 
	{
		for (side=0; side<2; side++) 
		{
			//ierr = IGASetBoundaryValue(iga,dir,side,dof,0.0);CHKERRQ(ierr);					// Dirichlet boundary conditions
			ierr = IGASetBoundaryForm(igaEL,dir,side,PETSC_TRUE);CHKERRQ(ierr);  				// Neumann boundary conditions
		}
	}

	Mat KEL;
	Vec el0,FEL;

	ierr = IGACreateMat(igaEL,&KEL);CHKERRQ(ierr);
	ierr = IGACreateVec(igaEL,&el0);CHKERRQ(ierr);
	ierr = IGACreateVec(igaEL,&FEL);CHKERRQ(ierr);

	IGAPoint		pointEL;															//point
	IGAElement		elemEL;																//element
	PetscReal		*KlocEL,*FlocEL;													//AA y BB
	PetscReal		*KpointEL,*FpointEL;												//KKK y FFF
	const PetscReal	*arrayChi0EL,*arrayZ0EL,*arrayS0EL;									//arrayU
	Vec				localChi0EL,localZ0EL,localS0EL;									//localU
	PetscReal		*chi0EL,*Z0EL,*S0EL;												//U

  	IGAFormSystem	wtfEL;
 	void			*wtf2EL;

 	KSP kspEL;
	ierr = IGACreateKSP(igaEL,&kspEL);CHKERRQ(ierr);

	// Get local vectors Z0 and Chi0 and arrays
	ierr = IGAGetLocalVecArray(igachiUp,chiUp0,&localChi0EL,&arrayChi0EL);CHKERRQ(ierr);
	ierr = IGAGetLocalVecArray(igaZ0,Z0,&localZ0EL,&arrayZ0EL);CHKERRQ(ierr);
	ierr = IGAGetLocalVecArray(igaS,s0,&localS0EL,&arrayS0EL);CHKERRQ(ierr);

	// Element loop
	ierr = IGABeginElement(igaEL,&elemEL);CHKERRQ(ierr);
	ierr = IGABeginElement(igaZ0,&elemZ0);CHKERRQ(ierr);
	ierr = IGABeginElement(igachiUp,&elemchiUp);CHKERRQ(ierr);
	ierr = IGABeginElement(igaS,&elemS);CHKERRQ(ierr);

	while (IGANextElement(igaEL,elemEL)) 				
	{
		IGANextElement(igachiUp,elemchiUp);
		IGANextElement(igaZ0,elemZ0);
		IGANextElement(igaS,elemS);

		ierr = IGAElementGetWorkMat(elemEL,&KlocEL);CHKERRQ(ierr);
		ierr = IGAElementGetWorkVec(elemEL,&FlocEL);CHKERRQ(ierr);
		ierr = IGAElementGetValues(elemchiUp,arrayChi0EL,&chi0EL);CHKERRQ(ierr);
		ierr = IGAElementGetValues(elemZ0,arrayZ0EL,&Z0EL);CHKERRQ(ierr);
		ierr = IGAElementGetValues(elemS,arrayS0EL,&S0EL);CHKERRQ(ierr);

		// FormSystem loop
		while (IGAElementNextFormSystem(elemEL,&wtfEL,&wtf2EL))
		{
		// Quadrature loop
			ierr = IGAElementBeginPoint(elemEL,&pointEL);CHKERRQ(ierr);
			ierr = IGAElementBeginPoint(elemZ0,&pointZ0);CHKERRQ(ierr);
			ierr = IGAElementBeginPoint(elemchiUp,&pointchiUp);CHKERRQ(ierr);
			ierr = IGAElementBeginPoint(elemS,&pointS);CHKERRQ(ierr);

			while (IGAElementNextPoint(elemEL,pointEL))
			{
				if(pointEL->atboundary==0 && pointchiUp->atboundary==0 && pointZ0->atboundary==0)
				{
					IGAElementNextPoint(elemZ0,pointZ0);
					IGAElementNextPoint(elemchiUp,pointchiUp);
					IGAElementNextPoint(elemS,pointS);

					ierr = IGAPointGetWorkMat(pointEL,&KpointEL);CHKERRQ(ierr);
					ierr = IGAPointGetWorkVec(pointEL,&FpointEL);CHKERRQ(ierr);
					//	   EnergyDensity(IGAPoint p, IGAPoint pChi, IGAPoint pZu, IGAPoint pS, PetscReal *K, PetscReal *F, PetscReal *Chi, PetscReal *Zu, PetscReal *S0, void *ctx)	
					ierr = ElasEnergyDensity(pointEL,pointchiUp,pointZ0,pointS,KpointEL,FpointEL,chi0EL,Z0EL,S0EL,NULL);CHKERRQ(ierr);
					ierr = IGAPointAddMat(pointEL,KpointEL,KlocEL);CHKERRQ(ierr);
					ierr = IGAPointAddVec(pointEL,FpointEL,FlocEL);CHKERRQ(ierr);
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
			ierr = IGAElementEndPoint(elemEL,&pointEL);CHKERRQ(ierr);
			ierr = IGAElementEndPoint(elemchiUp,&pointchiUp);CHKERRQ(ierr);
			ierr = IGAElementEndPoint(elemZ0,&pointZ0);CHKERRQ(ierr);
			ierr = IGAElementEndPoint(elemS,&pointS);CHKERRQ(ierr);
		}

		//ierr = IGAElementFixSystem(elemStress,KlocStress,FlocStress);CHKERRQ(ierr);					//This sets Dirichlet condition ¿? (Yes, this applies the conditions from IGASetBoundaryValue)
		ierr = IGAElementAssembleMat(elemEL,KlocEL,KEL);CHKERRQ(ierr);
		ierr = IGAElementAssembleVec(elemEL,FlocEL,FEL);CHKERRQ(ierr);
	}
	IGANextElement(igachiUp,elemchiUp);
	IGANextElement(igaZ0,elemZ0);
	IGANextElement(igaS,elemS);

	ierr = IGAEndElement(igaEL,&elemEL);CHKERRQ(ierr);
	ierr = IGAEndElement(igaZ0,&elemZ0);CHKERRQ(ierr);
	ierr = IGAEndElement(igachiUp,&elemchiUp);CHKERRQ(ierr);
	ierr = IGAEndElement(igaS,&elemS);CHKERRQ(ierr);

	// Restore local vectors u, Z0, Chi0 and arrays
	ierr = IGARestoreLocalVecArray(igachiUp,chiUp0,&localChi0EL,&arrayChi0EL);CHKERRQ(ierr);
	ierr = IGARestoreLocalVecArray(igaZ0,Z0,&localZ0EL,&arrayZ0EL);CHKERRQ(ierr);
	ierr = IGARestoreLocalVecArray(igaS,s0,&localS0EL,&arrayS0EL);CHKERRQ(ierr);

	ierr = MatAssemblyBegin(KEL,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyEnd  (KEL,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

	ierr = VecAssemblyBegin(FEL);CHKERRQ(ierr);
	ierr = VecAssemblyEnd  (FEL);CHKERRQ(ierr);

	ierr = KSPSetOperators(kspEL,KEL,KEL);CHKERRQ(ierr);
	PC pcEL;
	ierr = KSPGetPC(kspEL,&pcEL); CHKERRQ(ierr);
	ierr = PCSetType(pcEL,PCLU); CHKERRQ(ierr);
	ierr = PCFactorSetMatSolverType(pcEL,MATSOLVERMUMPS); CHKERRQ(ierr);
	//ierr = KSPSetFromOptions(kspStress);CHKERRQ(ierr);
	//ierr = KSPSetTolerances(kspStress,1e-28,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);CHKERRQ(ierr);
	ierr = KSPSolve(kspEL,FEL,el0);CHKERRQ(ierr);

	ierr = KSPDestroy(&kspEL);CHKERRQ(ierr);
	ierr = MatDestroy(&KEL);CHKERRQ(ierr);
	ierr = VecDestroy(&FEL);CHKERRQ(ierr);

	char nameEL[]="/ElasticEnergy.dat";
	char pathEL[512];
	sprintf(pathEL,"%s%s",direct,nameEL);
	ierr = IGAWriteVec(igaEL,el0,pathEL);CHKERRQ(ierr);	
//

//System for L2 projection of Gradient Energy Density
	PetscPrintf(PETSC_COMM_WORLD,"\n System for Gradient Energy Density starting \n\n");
	IGA igaEG;
	ierr = IGACreate(PETSC_COMM_WORLD,&igaEG);CHKERRQ(ierr);
	ierr = IGASetDim(igaEG,2);CHKERRQ(ierr);													//Spatial dimension of the problem
	ierr = IGASetDof(igaEG,1);CHKERRQ(ierr);													//Number of degrees of freedom, per node
	ierr = IGASetOrder(igaEG,2);CHKERRQ(ierr);												//Number of spatial derivatives to calculate
	ierr = IGASetFromOptions(igaEG);CHKERRQ(ierr);											//Note: The order (or degree) of the shape functions is given by the mesh!
	ierr = IGARead(igaEG,"./geometry3.dat");CHKERRQ(ierr);
	ierr = IGASetUp(igaEG);CHKERRQ(ierr);

	for (dir=0; dir<2; dir++) 
	{
		for (side=0; side<2; side++) 
		{
			//ierr = IGASetBoundaryValue(iga,dir,side,dof,0.0);CHKERRQ(ierr);					// Dirichlet boundary conditions
			ierr = IGASetBoundaryForm(igaEG,dir,side,PETSC_TRUE);CHKERRQ(ierr);  				// Neumann boundary conditions
		}
	}

	Mat KEG;
	Vec eg0,FEG;

	ierr = IGACreateMat(igaEG,&KEG);CHKERRQ(ierr);
	ierr = IGACreateVec(igaEG,&eg0);CHKERRQ(ierr);
	ierr = IGACreateVec(igaEG,&FEG);CHKERRQ(ierr);

	IGAPoint		pointEG;															//point
	IGAElement		elemEG;																//element
	PetscReal		*KlocEG,*FlocEG;													//AA y BB
	PetscReal		*KpointEG,*FpointEG;												//KKK y FFF
	const PetscReal	*arrayChi0EG,*arrayZ0EG,*arrayS0EG;									//arrayU
	Vec				localChi0EG,localZ0EG,localS0EG;									//localU
	PetscReal		*chi0EG,*Z0EG,*S0EG;												//U

  	IGAFormSystem	wtfEG;
 	void			*wtf2EG;

 	KSP kspEG;
	ierr = IGACreateKSP(igaEG,&kspEG);CHKERRQ(ierr);

	// Get local vectors Z0 and Chi0 and arrays
	ierr = IGAGetLocalVecArray(igachiUp,chiUp0,&localChi0EG,&arrayChi0EG);CHKERRQ(ierr);
	ierr = IGAGetLocalVecArray(igaZ0,Z0,&localZ0EG,&arrayZ0EG);CHKERRQ(ierr);
	ierr = IGAGetLocalVecArray(igaS,s0,&localS0EG,&arrayS0EG);CHKERRQ(ierr);

	// Element loop
	ierr = IGABeginElement(igaEG,&elemEG);CHKERRQ(ierr);
	ierr = IGABeginElement(igaZ0,&elemZ0);CHKERRQ(ierr);
	ierr = IGABeginElement(igachiUp,&elemchiUp);CHKERRQ(ierr);
	ierr = IGABeginElement(igaS,&elemS);CHKERRQ(ierr);

	while (IGANextElement(igaEG,elemEG)) 				
	{
		IGANextElement(igachiUp,elemchiUp);
		IGANextElement(igaZ0,elemZ0);
		IGANextElement(igaS,elemS);

		ierr = IGAElementGetWorkMat(elemEG,&KlocEG);CHKERRQ(ierr);
		ierr = IGAElementGetWorkVec(elemEG,&FlocEG);CHKERRQ(ierr);
		ierr = IGAElementGetValues(elemchiUp,arrayChi0EG,&chi0EG);CHKERRQ(ierr);
		ierr = IGAElementGetValues(elemZ0,arrayZ0EG,&Z0EG);CHKERRQ(ierr);
		ierr = IGAElementGetValues(elemS,arrayS0EG,&S0EG);CHKERRQ(ierr);

		// FormSystem loop
		while (IGAElementNextFormSystem(elemEG,&wtfEG,&wtf2EG))
		{
		// Quadrature loop
			ierr = IGAElementBeginPoint(elemEG,&pointEG);CHKERRQ(ierr);
			ierr = IGAElementBeginPoint(elemZ0,&pointZ0);CHKERRQ(ierr);
			ierr = IGAElementBeginPoint(elemchiUp,&pointchiUp);CHKERRQ(ierr);
			ierr = IGAElementBeginPoint(elemS,&pointS);CHKERRQ(ierr);

			while (IGAElementNextPoint(elemEG,pointEG))
			{
				if(pointEG->atboundary==0 && pointchiUp->atboundary==0 && pointZ0->atboundary==0)
				{
					IGAElementNextPoint(elemZ0,pointZ0);
					IGAElementNextPoint(elemchiUp,pointchiUp);
					IGAElementNextPoint(elemS,pointS);

					ierr = IGAPointGetWorkMat(pointEG,&KpointEG);CHKERRQ(ierr);
					ierr = IGAPointGetWorkVec(pointEG,&FpointEG);CHKERRQ(ierr);
					//	   EnergyDensity(IGAPoint p, IGAPoint pChi, IGAPoint pZu, IGAPoint pS, PetscReal *K, PetscReal *F, PetscReal *Chi, PetscReal *Zu, PetscReal *S0, void *ctx)	
					ierr = GradEnergyDensity(pointEG,pointchiUp,pointZ0,pointS,KpointEG,FpointEG,chi0EG,Z0EG,S0EG,NULL);CHKERRQ(ierr);
					ierr = IGAPointAddMat(pointEG,KpointEG,KlocEG);CHKERRQ(ierr);
					ierr = IGAPointAddVec(pointEG,FpointEG,FlocEG);CHKERRQ(ierr);
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
			ierr = IGAElementEndPoint(elemEG,&pointEG);CHKERRQ(ierr);
			ierr = IGAElementEndPoint(elemchiUp,&pointchiUp);CHKERRQ(ierr);
			ierr = IGAElementEndPoint(elemZ0,&pointZ0);CHKERRQ(ierr);
			ierr = IGAElementEndPoint(elemS,&pointS);CHKERRQ(ierr);
		}

		//ierr = IGAElementFixSystem(elemStress,KlocStress,FlocStress);CHKERRQ(ierr);					//This sets Dirichlet condition ¿? (Yes, this applies the conditions from IGASetBoundaryValue)
		ierr = IGAElementAssembleMat(elemEG,KlocEG,KEG);CHKERRQ(ierr);
		ierr = IGAElementAssembleVec(elemEG,FlocEG,FEG);CHKERRQ(ierr);
	}
	IGANextElement(igachiUp,elemchiUp);
	IGANextElement(igaZ0,elemZ0);
	IGANextElement(igaS,elemS);

	ierr = IGAEndElement(igaEG,&elemEG);CHKERRQ(ierr);
	ierr = IGAEndElement(igaZ0,&elemZ0);CHKERRQ(ierr);
	ierr = IGAEndElement(igachiUp,&elemchiUp);CHKERRQ(ierr);
	ierr = IGAEndElement(igaS,&elemS);CHKERRQ(ierr);

	// Restore local vectors u, Z0, Chi0 and arrays
	ierr = IGARestoreLocalVecArray(igachiUp,chiUp0,&localChi0EG,&arrayChi0EG);CHKERRQ(ierr);
	ierr = IGARestoreLocalVecArray(igaZ0,Z0,&localZ0EG,&arrayZ0EG);CHKERRQ(ierr);
	ierr = IGARestoreLocalVecArray(igaS,s0,&localS0EG,&arrayS0EG);CHKERRQ(ierr);

	ierr = MatAssemblyBegin(KEG,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyEnd  (KEG,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

	ierr = VecAssemblyBegin(FEG);CHKERRQ(ierr);
	ierr = VecAssemblyEnd  (FEG);CHKERRQ(ierr);

	ierr = KSPSetOperators(kspEG,KEG,KEG);CHKERRQ(ierr);
	PC pcEG;
	ierr = KSPGetPC(kspEG,&pcEG); CHKERRQ(ierr);
	ierr = PCSetType(pcEG,PCLU); CHKERRQ(ierr);
	ierr = PCFactorSetMatSolverType(pcEG,MATSOLVERMUMPS); CHKERRQ(ierr);
	//ierr = KSPSetFromOptions(kspStress);CHKERRQ(ierr);
	//ierr = KSPSetTolerances(kspStress,1e-28,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);CHKERRQ(ierr);
	ierr = KSPSolve(kspEG,FEG,eg0);CHKERRQ(ierr);

	ierr = KSPDestroy(&kspEG);CHKERRQ(ierr);
	ierr = MatDestroy(&KEG);CHKERRQ(ierr);
	ierr = VecDestroy(&FEG);CHKERRQ(ierr);

	char nameEG[]="/GradEnergy.dat";
	char pathEG[512];
	sprintf(pathEG,"%s%s",direct,nameEG);
	ierr = IGAWriteVec(igaEG,eg0,pathEG);CHKERRQ(ierr);	
//
*/

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
			ierr = IGASetBoundaryForm(igaFullS,dir,side,PETSC_TRUE);CHKERRQ(ierr);  				// Neumann boundary conditions
		}
	}

	Mat KFullS;
	Vec Fs0,FFullS;

	ierr = IGACreateMat(igaFullS,&KFullS);CHKERRQ(ierr);
	ierr = IGACreateVec(igaFullS,&Fs0);CHKERRQ(ierr);
	ierr = IGACreateVec(igaFullS,&FFullS);CHKERRQ(ierr);

	IGAPoint		pointFullS;															//point
	IGAElement		elemFullS;																//element
	PetscReal		*KlocFullS,*FlocFullS;												//AA y BB
	PetscReal		*KpointFullS,*FpointFullS;											//KKK y FFF
	const PetscReal	*arrayZ0FullS,*arrayChi0FullS, *arraySFullS;			//arrayU
	Vec				localZ0FullS,localChi0FullS, localSFullS;				//localU
	PetscReal		*Chi0FullS,*Z0FullS, *SFullS;											//U

  	IGAFormSystem	wtfFullS;
 	void			*wtf2FullS;

 	KSP kspFullS;
	ierr = IGACreateKSP(igaFullS,&kspFullS);CHKERRQ(ierr);

	// Get local vectors Z0 and Chi0 and arrays
	ierr = IGAGetLocalVecArray(igachiUp,chiUp0,&localChi0FullS,&arrayChi0FullS);CHKERRQ(ierr);
	ierr = IGAGetLocalVecArray(igaZ0,Z0,&localZ0FullS,&arrayZ0FullS);CHKERRQ(ierr);
	ierr = IGAGetLocalVecArray(igaS,s0,&localSFullS,&arraySFullS);CHKERRQ(ierr);

	// Element loop
	ierr = IGABeginElement(igaFullS,&elemFullS);CHKERRQ(ierr);
	ierr = IGABeginElement(igaZ0,&elemZ0);CHKERRQ(ierr);
	ierr = IGABeginElement(igachiUp,&elemchiUp);CHKERRQ(ierr);
	ierr = IGABeginElement(igaS,&elemS);CHKERRQ(ierr);

	while (IGANextElement(igaFullS,elemFullS)) 
	{
		IGANextElement(igaZ0,elemZ0);
		IGANextElement(igachiUp,elemchiUp);
		IGANextElement(igaS,elemS);

		ierr = IGAElementGetWorkMat(elemFullS,&KlocFullS);CHKERRQ(ierr);
		ierr = IGAElementGetWorkVec(elemFullS,&FlocFullS);CHKERRQ(ierr);
		ierr = IGAElementGetValues(elemZ0,arrayZ0FullS,&Z0FullS);CHKERRQ(ierr);
		ierr = IGAElementGetValues(elemchiUp,arrayChi0FullS,&Chi0FullS);CHKERRQ(ierr);
		ierr = IGAElementGetValues(elemS,arraySFullS,&SFullS);CHKERRQ(ierr);

		// FormSystem loop
		while (IGAElementNextFormSystem(elemFullS,&wtfFullS,&wtf2FullS)) 
		{
		// Quadrature loop
			ierr = IGAElementBeginPoint(elemFullS,&pointFullS);CHKERRQ(ierr);
			ierr = IGAElementBeginPoint(elemZ0,&pointZ0);CHKERRQ(ierr);
			ierr = IGAElementBeginPoint(elemchiUp,&pointchiUp);CHKERRQ(ierr);
			ierr = IGAElementBeginPoint(elemS,&pointS);CHKERRQ(ierr);

			while (IGAElementNextPoint(elemFullS,pointFullS))
			{
				if(pointFullS->atboundary==1)
				{
					//IGAElementNextPoint(elemZ0,pointZ0);
					//IGAElementNextPoint(elemchiUp,pointchiUp);
				}

				if(pointZ0->atboundary==1)
				{
					PetscPrintf(PETSC_COMM_WORLD,"Hola1 en FullS \n");
				}

				if(pointFullS->atboundary==0 && pointZ0->atboundary==0 && pointchiUp->atboundary==0)
				{
					IGAElementNextPoint(elemZ0,pointZ0);
					IGAElementNextPoint(elemchiUp,pointchiUp);

					ierr = IGAPointGetWorkMat(pointFullS,&KpointFullS);CHKERRQ(ierr);
					ierr = IGAPointGetWorkVec(pointFullS,&FpointFullS);CHKERRQ(ierr);
					//	   FullS(IGAPoint p,IGAPoint pChi,IGAPoint pZu,IGAPoint pS, PetscReal *K,PetscReal *F, PetscReal *Chi,PetscReal *Zu,PetscReal *S, void *ctx)
					ierr = FullS(pointFullS,pointchiUp,pointZ0,pointS,KpointFullS,FpointFullS,Chi0FullS,Z0FullS,SFullS,NULL);CHKERRQ(ierr);
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
			while (pointS->index != -1)
			{
				IGAElementNextPoint(elemS,pointS);
			}

			ierr = IGAElementEndPoint(elemFullS,&pointFullS);CHKERRQ(ierr);
			ierr = IGAElementEndPoint(elemZ0,&pointZ0);CHKERRQ(ierr);
			ierr = IGAElementEndPoint(elemchiUp,&pointchiUp);CHKERRQ(ierr);
			ierr = IGAElementEndPoint(elemS,&pointS);CHKERRQ(ierr);
		}

		//ierr = IGAElementFixSystem(elemFullS,KlocFullS,FlocFullS);CHKERRQ(ierr);					//This sets Dirichlet condition ¿? (Yes, this applies the conditions from IGASetBoundaryValue)
		ierr = IGAElementAssembleMat(elemFullS,KlocFullS,KFullS);CHKERRQ(ierr);
		ierr = IGAElementAssembleVec(elemFullS,FlocFullS,FFullS);CHKERRQ(ierr);
	}
	IGANextElement(igaZ0,elemZ0);
	IGANextElement(igachiUp,elemchiUp);
	IGANextElement(igaS,elemS);

	ierr = IGAEndElement(igaFullS,&elemFullS);CHKERRQ(ierr);
	ierr = IGAEndElement(igaZ0,&elemZ0);CHKERRQ(ierr);
	ierr = IGAEndElement(igachiUp,&elemchiUp);CHKERRQ(ierr);
	ierr = IGAEndElement(igaS,&elemS);CHKERRQ(ierr);

	// Restore local vectors u, Z0, Chi0 and arrays
	ierr = IGARestoreLocalVecArray(igaZ0,Z0,&localZ0FullS,&arrayZ0FullS);CHKERRQ(ierr);
	ierr = IGARestoreLocalVecArray(igachiUp,chiUp0,&localChi0FullS,&arrayChi0FullS);CHKERRQ(ierr);
	ierr = IGARestoreLocalVecArray(igaS,s0,&localSFullS,&arraySFullS);CHKERRQ(ierr);

	ierr = MatAssemblyBegin(KFullS,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyEnd  (KFullS,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

	ierr = VecAssemblyBegin(FFullS);CHKERRQ(ierr);
	ierr = VecAssemblyEnd  (FFullS);CHKERRQ(ierr);

	ierr = KSPSetOperators(kspFullS,KFullS,KFullS);CHKERRQ(ierr);
	//PC pcFullS;
	//ierr = KSPGetPC(kspFullS,&pcFullS); CHKERRQ(ierr);
	//ierr = PCSetType(pcFullS,PCLU); CHKERRQ(ierr);
	//ierr = PCFactorSetMatSolverType(pcFullS,MATSOLVERMUMPS); CHKERRQ(ierr);
	ierr = KSPSetFromOptions(kspFullS);CHKERRQ(ierr);
	ierr = KSPSetTolerances(kspFullS,1e-24,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);CHKERRQ(ierr);
	ierr = KSPSolve(kspFullS,FFullS,Fs0);CHKERRQ(ierr);

	ierr = KSPDestroy(&kspFullS);CHKERRQ(ierr);
	ierr = MatDestroy(&KFullS);CHKERRQ(ierr);
	ierr = VecDestroy(&FFullS);CHKERRQ(ierr);

	char nameFullS[]="/FullSigma-2d-0.dat";
	char pathFullS[512];
	sprintf(pathFullS,"%s%s",direct,nameFullS);
	ierr = IGAWriteVec(igaFullS,Fs0,pathFullS);CHKERRQ(ierr);	
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
	//PC pcVa;
	//ierr = KSPGetPC(kspStress,&pcStress); CHKERRQ(ierr);
	//ierr = PCSetType(pcStress,PCLU); CHKERRQ(ierr);
	//ierr = PCFactorSetMatSolverType(pcStress,MATSOLVERMUMPS); CHKERRQ(ierr);
	ierr = KSPSetFromOptions(kspVa);CHKERRQ(ierr);
	ierr = KSPSetTolerances(kspVa,1.0e-24,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);CHKERRQ(ierr);
	ierr = KSPSolve(kspVa,FVa,Va0);CHKERRQ(ierr);

	ierr = KSPDestroy(&kspVa);CHKERRQ(ierr);
	ierr = MatDestroy(&KVa);CHKERRQ(ierr);
	ierr = VecDestroy(&FVa);CHKERRQ(ierr);

	char nameVa[]="/Va-2d-0.dat";
	char pathVa[512];
	sprintf(pathVa,"%s%s",direct,nameVa);
	ierr = IGAWriteVec(igaVa,Va0,pathVa);CHKERRQ(ierr);	
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

	ierr = KSPDestroy(&kspVs);CHKERRQ(ierr);
	ierr = MatDestroy(&KVs);CHKERRQ(ierr);
	ierr = VecDestroy(&FVs);CHKERRQ(ierr);

	char nameVs[]="/Vs-2d-0.dat";
	char pathVs[512];
	sprintf(pathVs,"%s%s",direct,nameVs);
	ierr = IGAWriteVec(igaVs,Vs0,pathVs);CHKERRQ(ierr);	
//

//System for L2 projection of exact stress
	PetscPrintf(PETSC_COMM_WORLD,"\nSystem for L2 projection for exact stress starting \n\n");
	T=time(NULL);
	tm=*localtime(&T);
	PetscPrintf(PETSC_COMM_WORLD,"Current time is %02d:%02d:%02d \n",tm.tm_hour,tm.tm_min,tm.tm_sec);
	IGA igaExact;
	ierr = IGACreate(PETSC_COMM_WORLD,&igaExact);CHKERRQ(ierr);
	ierr = IGASetDim(igaExact,2);CHKERRQ(ierr);														//Spatial dimension of the problem
	ierr = IGASetDof(igaExact,4);CHKERRQ(ierr);														//Number of degrees of freedom, per node
	ierr = IGASetOrder(igaExact,2);CHKERRQ(ierr);													//Number of spatial derivatives to calculate
	ierr = IGASetFromOptions(igaExact);CHKERRQ(ierr);												//Note: The order (or degree) of the shape functions is given by the mesh!
	ierr = IGARead(igaExact,"./geometry3.dat");CHKERRQ(ierr);
	
	for (dir=0; dir<2; dir++)
	{
		ierr = IGASetRuleType(igaExact,dir,IGA_RULE_LEGENDRE);CHKERRQ(ierr);
		ierr = IGASetRuleSize(igaExact,dir,6);CHKERRQ(ierr);
	}
	ierr = IGASetUp(igaExact);CHKERRQ(ierr);
	
	for (dir=0; dir<2; dir++) 
	{
		for (side=0; side<2; side++) 
		{
			//ierr = IGASetBoundaryValue(iga,dir,side,dof,0.0);CHKERRQ(ierr);    				// Dirichlet boundary conditions
			ierr = IGASetBoundaryForm(igaExact,dir,side,PETSC_TRUE);CHKERRQ(ierr);  				// Neumann boundary conditions
		}
	}

	Vec e0;
	Mat Kl2e;
	Vec Fl2e;
	ierr = IGACreateVec(igaExact,&e0);CHKERRQ(ierr);  
	ierr = IGACreateMat(igaExact,&Kl2e);CHKERRQ(ierr);
	ierr = IGACreateVec(igaExact,&Fl2e);CHKERRQ(ierr);
	ierr = IGASetFormSystem(igaExact,L2ProjectionExactStress,&userL2);CHKERRQ(ierr);
	ierr = IGAComputeSystem(igaExact,Kl2e,Fl2e);CHKERRQ(ierr);

	//This parts set and calls KSP to solve the linear system
	KSP kspl2e;
	ierr = IGACreateKSP(igaExact,&kspl2e);CHKERRQ(ierr);										
	ierr = KSPSetOperators(kspl2e,Kl2e,Kl2e);CHKERRQ(ierr); 								//This function creates the matrix for the system on the second parameter and uses the 3rd parameter as a preconditioner
	ierr = KSPSetType(kspl2e,KSPCG);CHKERRQ(ierr);											//Using KSPCG (conjugated gradient) because the matrix is symmetric
	//ierr = KSPSetOptionsPrefix(kspl2S,"l2pS_");CHKERRQ(ierr);
	ierr = KSPSetFromOptions(kspl2e);CHKERRQ(ierr);
	ierr = KSPSetTolerances(kspl2e,1.0e-24,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);CHKERRQ(ierr);
	ierr = KSPSolve(kspl2e,Fl2e,e0);CHKERRQ(ierr);										//This is a simple system, so it can be solved with just this command

	ierr = KSPDestroy(&kspl2e);CHKERRQ(ierr);
	ierr = MatDestroy(&Kl2e);CHKERRQ(ierr);
	ierr = VecDestroy(&Fl2e);CHKERRQ(ierr);
	char nameE[]="/ExactStress-2d-0.dat";
	char pathE[512];
	sprintf(pathE,"%s%s",direct,nameE);
	ierr = IGAWriteVec(igaExact,e0,pathE);CHKERRQ(ierr);
//

//System for L2 projection of grad(Z0)
	PetscPrintf(PETSC_COMM_WORLD,"\nSystem for L2 projection for grad(Z0) \n\n");
	T=time(NULL);
	tm=*localtime(&T);
	PetscPrintf(PETSC_COMM_WORLD,"Current time is %02d:%02d:%02d \n",tm.tm_hour,tm.tm_min,tm.tm_sec);
	IGA igaGrad;
	ierr = IGACreate(PETSC_COMM_WORLD,&igaGrad);CHKERRQ(ierr);
	ierr = IGASetDim(igaGrad,2);CHKERRQ(ierr);														//Spatial dimension of the problem
	ierr = IGASetDof(igaGrad,4);CHKERRQ(ierr);														//Number of degrees of freedom, per node
	ierr = IGASetOrder(igaGrad,2);CHKERRQ(ierr);													//Number of spatial derivatives to calculate
	ierr = IGASetFromOptions(igaGrad);CHKERRQ(ierr);												//Note: The order (or degree) of the shape functions is given by the mesh!
	ierr = IGARead(igaGrad,"./geometry3.dat");CHKERRQ(ierr);
	
	for (dir=0; dir<2; dir++)
	{
		ierr = IGASetRuleType(igaGrad,dir,IGA_RULE_LEGENDRE);CHKERRQ(ierr);
		ierr = IGASetRuleSize(igaGrad,dir,6);CHKERRQ(ierr);
	}
	ierr = IGASetUp(igaGrad);CHKERRQ(ierr);
	
	for (dir=0; dir<2; dir++) 
	{
		for (side=0; side<2; side++) 
		{
			//ierr = IGASetBoundaryValue(iga,dir,side,dof,0.0);CHKERRQ(ierr);    				// Dirichlet boundary conditions
			ierr = IGASetBoundaryForm(igaGrad,dir,side,PETSC_TRUE);CHKERRQ(ierr);  				// Neumann boundary conditions
		}
	}

	Vec gradZ0;
	Mat Kl2Grad;
	Vec Fl2Grad;
	ierr = IGACreateVec(igaGrad,&gradZ0);CHKERRQ(ierr);  
	ierr = IGACreateMat(igaGrad,&Kl2Grad);CHKERRQ(ierr);
	ierr = IGACreateVec(igaGrad,&Fl2Grad);CHKERRQ(ierr);
	
	IGAPoint		pointGrad;						//point
	IGAElement		elemGrad;						//element
	PetscReal		*KlocGrad,*FlocGrad;			//AA y BB
	PetscReal		*KpointGrad,*FpointGrad;		//KKK y FFF
	const PetscReal	*arrayGradZ0;					//arrayU
	Vec				localGradZ0;					//localU
	PetscReal		*GradZ0;						//U

  	IGAFormSystem	wtfGrad;
 	void			*wtf2Grad;

	// Get local vectors Chi0 and arrays
	ierr = IGAGetLocalVecArray(igaZ0,Z0,&localGradZ0,&arrayGradZ0);CHKERRQ(ierr);

	// Element loop
	ierr = IGABeginElement(igaGrad,&elemGrad);CHKERRQ(ierr);
	ierr = IGABeginElement(igaZ0,&elemZ0);CHKERRQ(ierr);

	while (IGANextElement(igaGrad,elemGrad)) 
	{
		IGANextElement(igaZ0,elemZ0);

		ierr = IGAElementGetWorkMat(elemGrad,&KlocGrad);CHKERRQ(ierr);
		ierr = IGAElementGetWorkVec(elemGrad,&FlocGrad);CHKERRQ(ierr);
		ierr = IGAElementGetValues(elemZ0,arrayGradZ0,&GradZ0);CHKERRQ(ierr);

		// FormSystem loop
		while (IGAElementNextFormSystem(elemGrad,&wtfGrad,&wtf2Grad)) 
		{
		// Quadrature loop
			ierr = IGAElementBeginPoint(elemGrad,&pointGrad);CHKERRQ(ierr);
			ierr = IGAElementBeginPoint(elemZ0,&pointZ0);CHKERRQ(ierr);

			while (IGAElementNextPoint(elemGrad,pointGrad))
			{
				if(pointGrad->atboundary==1)
				{
					ierr = IGAPointGetWorkMat(pointGrad,&KpointGrad);CHKERRQ(ierr);
					ierr = IGAPointGetWorkVec(pointGrad,&FpointGrad);CHKERRQ(ierr);
					ierr = L2ProjectionGradZ(pointGrad,pointZ0,KpointGrad,FpointGrad,GradZ0,NULL);CHKERRQ(ierr);
					ierr = IGAPointAddMat(pointGrad,KpointGrad,KlocGrad);CHKERRQ(ierr);
					ierr = IGAPointAddVec(pointGrad,FpointGrad,FlocGrad);CHKERRQ(ierr);
				}

				if(pointZ0->atboundary==1)
				{
					PetscPrintf(PETSC_COMM_WORLD,"Hola en Grad");
				}

				if(pointGrad->atboundary==0 && pointZ0->atboundary==0)
				{
					IGAElementNextPoint(elemZ0,pointZ0);

					ierr = IGAPointGetWorkMat(pointGrad,&KpointGrad);CHKERRQ(ierr);
					ierr = IGAPointGetWorkVec(pointGrad,&FpointGrad);CHKERRQ(ierr);
					ierr = L2ProjectionGradZ(pointGrad,pointZ0,KpointGrad,FpointGrad,GradZ0,NULL);CHKERRQ(ierr);
					ierr = IGAPointAddMat(pointGrad,KpointGrad,KlocGrad);CHKERRQ(ierr);
					ierr = IGAPointAddVec(pointGrad,FpointGrad,FlocGrad);CHKERRQ(ierr);
				}
			}
			if (pointZ0->index != -1)
			{
				IGAElementNextPoint(elemZ0,pointZ0);
			}

			ierr = IGAElementEndPoint(elemGrad,&pointGrad);CHKERRQ(ierr);
			ierr = IGAElementEndPoint(elemZ0,&pointZ0);CHKERRQ(ierr);
		}

		ierr = IGAElementAssembleMat(elemGrad,KlocGrad,Kl2Grad);CHKERRQ(ierr);
		ierr = IGAElementAssembleVec(elemGrad,FlocGrad,Fl2Grad);CHKERRQ(ierr);

	}
	IGANextElement(igaZ0,elemZ0);

	ierr = IGAEndElement(igaGrad,&elemGrad);CHKERRQ(ierr);
	ierr = IGAEndElement(igaZ0,&elemZ0);CHKERRQ(ierr);

	// Restore local vectors Chi0 and arrays
	ierr = IGARestoreLocalVecArray(igaZ0,Z0,&localGradZ0,&arrayGradZ0);CHKERRQ(ierr);

	ierr = MatAssemblyBegin(Kl2Grad,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyEnd  (Kl2Grad,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);	

	ierr = VecAssemblyBegin(Fl2Grad);CHKERRQ(ierr);
	ierr = VecAssemblyEnd  (Fl2Grad);CHKERRQ(ierr);

	//This parts set and calls KSP to solve the linear system
	KSP kspl2Grad;
	ierr = IGACreateKSP(igaGrad,&kspl2Grad);CHKERRQ(ierr);										
	ierr = KSPSetOperators(kspl2Grad,Kl2Grad,Kl2Grad);CHKERRQ(ierr); 								//This function creates the matrix for the system on the second parameter and uses the 3rd parameter as a preconditioner
	ierr = KSPSetType(kspl2Grad,KSPCG);CHKERRQ(ierr);											//Using KSPCG (conjugated gradient) because the matrix is symmetric
	//ierr = KSPSetOptionsPrefix(kspl2S,"l2pS_");CHKERRQ(ierr);
	ierr = KSPSetFromOptions(kspl2Grad);CHKERRQ(ierr);
	ierr = KSPSetTolerances(kspl2Grad,1.0e-24,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);CHKERRQ(ierr);
	ierr = KSPSolve(kspl2Grad,Fl2Grad,gradZ0);CHKERRQ(ierr);										//This is a simple system, so it can be solved with just this command

	ierr = KSPDestroy(&kspl2Grad);CHKERRQ(ierr);
	ierr = MatDestroy(&Kl2Grad);CHKERRQ(ierr);
	ierr = VecDestroy(&Fl2Grad);CHKERRQ(ierr);
	char nameGrad[]="/gradZ0-2d-0.dat";
	char pathGrad[512];
	sprintf(pathGrad,"%s%s",direct,nameGrad);
	ierr = IGAWriteVec(igaGrad,gradZ0,pathGrad);CHKERRQ(ierr);
//

/*
//Creation of types and systems for the Helmholtz decomposition of Z, curl part
	PetscPrintf(PETSC_COMM_WORLD,"\n System for curl part of Helmholtz of Z starting \n\n");
	IGA igachiZ;
	ierr = IGACreate(PETSC_COMM_WORLD,&igachiZ);CHKERRQ(ierr);
	ierr = IGASetDim(igachiZ,2);CHKERRQ(ierr);														//Spatial dimension of the problem
	ierr = IGASetDof(igachiZ,4);CHKERRQ(ierr);														//Number of degrees of freedom, per node
	ierr = IGASetOrder(igachiZ,2);CHKERRQ(ierr);													//Number of spatial derivatives to calculate
	ierr = IGASetFromOptions(igachiZ);CHKERRQ(ierr);												//Note: The order (or degree) of the shape functions is given by the mesh!
	ierr = IGARead(igachiZ,"./geometry2.dat");CHKERRQ(ierr);
	ierr = IGASetUp(igachiZ);CHKERRQ(ierr);

	//ierr = IGASetBoundaryValue(iga,dir,side,dof,val);CHKERRQ(ierr);
	//dir=0 and side=0
	ierr = IGASetBoundaryValue(igachiZ,0,0,0,0.0);CHKERRQ(ierr);					// Dirichlet boundary conditions to get chi \dot n =0
	ierr = IGASetBoundaryValue(igachiZ,0,0,2,0.0);CHKERRQ(ierr);
	//dir=0 and side=1
	ierr = IGASetBoundaryValue(igachiZ,0,1,0,0.0);CHKERRQ(ierr);					// Dirichlet boundary conditions to get chi \dot n =0
	ierr = IGASetBoundaryValue(igachiZ,0,1,2,0.0);CHKERRQ(ierr);
	//dir=1 and side=0
	ierr = IGASetBoundaryValue(igachiZ,1,0,1,0.0);CHKERRQ(ierr);					// Dirichlet boundary conditions to get chi \dot n =0
	ierr = IGASetBoundaryValue(igachiZ,1,0,3,0.0);CHKERRQ(ierr);
	//dir=1 and side=1
	ierr = IGASetBoundaryValue(igachiZ,1,1,1,0.0);CHKERRQ(ierr);					// Dirichlet boundary conditions to get chi \dot n =0
	ierr = IGASetBoundaryValue(igachiZ,1,1,3,0.0);CHKERRQ(ierr);
	
	Mat KchiZp;
	Vec chiZ,FchiZ;
	ierr = IGACreateMat(igachiZ,&KchiZp);CHKERRQ(ierr);
	ierr = IGACreateVec(igachiZ,&chiZ);CHKERRQ(ierr);
	ierr = IGACreateVec(igachiZ,&FchiZ);CHKERRQ(ierr);

	IGAPoint        pointchiZ;
	IGAElement      elemchiZ;							//element
	PetscReal       *KlocchiZ,*FlocchiZ;				//AA y BB
	PetscReal       *KpointchiZ,*FpointchiZ;			//KKK y FFF
	const PetscReal *arrayZS0chiZ;		//arrayU
	Vec  			localZS0chiZ;		//localU
	PetscReal       *ZS0chiZ;							//U

  	IGAFormSystem  wtfchiZ;
 	void           *wtf2chiZ;

 	KSP kspchiZ;
	ierr = IGACreateKSP(igachiZ,&kspchiZ);CHKERRQ(ierr);

	//Get local vectors s0 and arrays
	ierr = IGAGetLocalVecArray(igaZS,ZS0,&localZS0chiZ,&arrayZS0chiZ);CHKERRQ(ierr);

	//Element loop
	ierr = IGABeginElement(igachiZ,&elemchiZ);CHKERRQ(ierr);
	ierr = IGABeginElement(igaZS,&elemZS);CHKERRQ(ierr);

	while (IGANextElement(igachiZ,elemchiZ))
	{
		IGANextElement(igaZS,elemZS);
		ierr = IGAElementGetWorkMat(elemchiZ,&KlocchiZ);CHKERRQ(ierr);
		ierr = IGAElementGetWorkVec(elemchiZ,&FlocchiZ);CHKERRQ(ierr);
		ierr = IGAElementGetValues(elemZS,arrayZS0chiZ,&ZS0chiZ);CHKERRQ(ierr);

		//FormSystem loop
		while (IGAElementNextFormSystem(elemchiZ,&wtfchiZ,&wtf2chiZ)) 
		{
			//Quadrature loop
			ierr = IGAElementBeginPoint(elemchiZ,&pointchiZ);CHKERRQ(ierr);
			ierr = IGAElementBeginPoint(elemZS,&pointZS);CHKERRQ(ierr);

			while (IGAElementNextPoint(elemchiZ,pointchiZ)) 
			{
				IGAElementNextPoint(elemZS,pointZS);
				ierr = IGAPointGetWorkMat(pointchiZ,&KpointchiZ);CHKERRQ(ierr);
				ierr = IGAPointGetWorkVec(pointchiZ,&FpointchiZ);CHKERRQ(ierr);
				//	   curlZ(IGAPoint p,IGAPoint pZ, PetscReal *K,PetscReal *F,PetscReal *Z, void *ctx)
				ierr = curlZ(pointchiZ,pointZS,KpointchiZ,FpointchiZ,ZS0chiZ,NULL);CHKERRQ(ierr);
				ierr = IGAPointAddMat(pointchiZ,KpointchiZ,KlocchiZ);CHKERRQ(ierr);
				ierr = IGAPointAddVec(pointchiZ,FpointchiZ,FlocchiZ);CHKERRQ(ierr);
			}
			IGAElementNextPoint(elemZS,pointZS);

			while (pointZS->index != -1)
			{
				IGAElementNextPoint(elemZS,pointZS);
			}

			ierr = IGAElementEndPoint(elemchiZ,&pointchiZ);CHKERRQ(ierr);
			ierr = IGAElementEndPoint(elemZS,&pointZS);CHKERRQ(ierr);
		}
		ierr = IGAElementFixSystem(elemchiZ,KlocchiZ,FlocchiZ);CHKERRQ(ierr);
		ierr = IGAElementAssembleMat(elemchiZ,KlocchiZ,KchiZp);CHKERRQ(ierr);
		ierr = IGAElementAssembleVec(elemchiZ,FlocchiZ,FchiZ);CHKERRQ(ierr);

	}
	IGANextElement(igaZS,elemZS);
	ierr = IGAEndElement(igachiZ,&elemchiZ);CHKERRQ(ierr);
	ierr = IGAEndElement(igaZS,&elemZS);CHKERRQ(ierr);

	// Restore local vectors s0 and arrays
	ierr = IGARestoreLocalVecArray(igaZS,ZS0,&localZS0chiZ,&arrayZS0chiZ);CHKERRQ(ierr);

	//Form system matrix and vector
	ierr = MatAssemblyBegin(KchiZp,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyEnd  (KchiZp,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	
	ierr = VecAssemblyBegin(FchiZ);CHKERRQ(ierr);
	ierr = VecAssemblyEnd  (FchiZ);CHKERRQ(ierr);

	ierr = KSPSetOperators(kspchiZ,KchiZp,KchiZp);CHKERRQ(ierr);
	//PC pcChiZ;
	//ierr = KSPGetPC(kspchiZ,&pcChiZ); CHKERRQ(ierr);
	//ierr = PCSetType(pcChiZ,PCLU); CHKERRQ(ierr);
	//ierr = PCFactorSetMatSolverType(pcChiZ,MATSOLVERMUMPS); CHKERRQ(ierr);
	//ierr = KSPSetFromOptions(kspchiZ);CHKERRQ(ierr);
	ierr = KSPSetType(kspchiZ,KSPFGMRES);CHKERRQ(ierr);
	ierr = KSPSetTolerances(kspchiZ,1.0e-16,1.0e-20,PETSC_DEFAULT,PETSC_DEFAULT);CHKERRQ(ierr);
	ierr = KSPSolve(kspchiZ,FchiZ,chiZ);CHKERRQ(ierr);

	ierr = KSPDestroy(&kspchiZ);CHKERRQ(ierr);
	ierr = MatDestroy(&KchiZp);CHKERRQ(ierr);
	ierr = VecDestroy(&FchiZ);CHKERRQ(ierr);

	//Si hay problemas, borrar esto
	//ierr = VecDestroy(&alp0); CHKERRQ(ierr);
	//hasta aquí

	char namechiZ[]="/ChiZ-2d-0.dat";
	char pathchiZ[512];
	sprintf(pathchiZ,"%s%s",direct,namechiZ);
	ierr = IGAWriteVec(igachiZ,chiZ,pathchiZ);CHKERRQ(ierr);
//
*/

/*
//Creation of types and systems for the Helmholtz decomposition of S, grad part
	//System for chiS
	PetscPrintf(PETSC_COMM_WORLD,"\n System for grad part of Z starting \n\n");
	IGA igaZa;
	ierr = IGACreate(PETSC_COMM_WORLD,&igaZa);CHKERRQ(ierr);
	ierr = IGASetDim(igaZa,2);CHKERRQ(ierr);														//Spatial dimension of the problem
	ierr = IGASetDof(igaZa,2);CHKERRQ(ierr);														//Number of degrees of freedom, per node
	ierr = IGASetOrder(igaZa,2);CHKERRQ(ierr);														//Number of spatial derivatives to calculate
	ierr = IGASetFromOptions(igaZa);CHKERRQ(ierr);													//Note: The order (or degree) of the shape functions is given by the mesh!
	ierr = IGARead(igaZa,"./geometry2.dat");CHKERRQ(ierr);
	ierr = IGASetUp(igaZa);CHKERRQ(ierr);
	
	//PetscInt dir,side;
	for (dir=0; dir<2; dir++) 
	{
		for (side=0; side<2; side++) 
		{
			//ierr = IGASetBoundaryValue(iga,dir,side,dof,0.0);CHKERRQ(ierr);    				// Dirichlet boundary conditions
			ierr = IGASetBoundaryForm(igaZa,dir,side,PETSC_TRUE);CHKERRQ(ierr);  				// Neumann boundary conditions
		}
	}
	
	Mat KZa;
	Vec Za0,FZa;
	ierr = IGACreateMat(igaZa,&KZa);CHKERRQ(ierr);
	ierr = IGACreateVec(igaZa,&Za0);CHKERRQ(ierr);
	ierr = IGACreateVec(igaZa,&FZa);CHKERRQ(ierr);

	IGAPoint        pointZa;
	IGAElement      elemZa;					//element
	PetscReal       *KlocZa,*FlocZa;		//AA y BB
	PetscReal       *KpointZa,*FpointZa;	//KKK y FFF
	const PetscReal *arrayZSZa;				//arrayU
	Vec  			localZSZa;				//localU
	PetscReal       *ZSZa;					//U

  	IGAFormSystem  wtfZa;
 	void           *wtf2Za;

 	KSP kspZa;
	ierr = IGACreateKSP(igaZa,&kspZa);CHKERRQ(ierr);

	//Get local vectors s0 and arrays
	ierr = IGAGetLocalVecArray(igaZS,ZS0,&localZSZa,&arrayZSZa);CHKERRQ(ierr);

	//Element loop
	ierr = IGABeginElement(igaZa,&elemZa);CHKERRQ(ierr);
	ierr = IGABeginElement(igaZS,&elemZS);CHKERRQ(ierr);

	while (IGANextElement(igaZa,elemZa))
	{
		IGANextElement(igaZS,elemZS);
		ierr = IGAElementGetWorkMat(elemZa,&KlocZa);CHKERRQ(ierr);
		ierr = IGAElementGetWorkVec(elemZa,&FlocZa);CHKERRQ(ierr);
		ierr = IGAElementGetValues(elemZS,arrayZSZa,&ZSZa);CHKERRQ(ierr);

		//FormSystem loop
		while (IGAElementNextFormSystem(elemZa,&wtfZa,&wtf2Za)) 
		{
			//Quadrature loop
			ierr = IGAElementBeginPoint(elemZa,&pointZa);CHKERRQ(ierr);
			ierr = IGAElementBeginPoint(elemZS,&pointZS);CHKERRQ(ierr);
			if(pointZa->atboundary==0)
			{
				while (IGAElementNextPoint(elemZa,pointZa)) 
				{
					IGAElementNextPoint(elemZS,pointZS);
					ierr = IGAPointGetWorkMat(pointZa,&KpointZa);CHKERRQ(ierr);
					ierr = IGAPointGetWorkVec(pointZa,&FpointZa);CHKERRQ(ierr);
					//	   gradZa(IGAPoint p,IGAPoint pZ,PetscReal *K,PetscReal *F,PetscReal *Z,void *ctx)
					ierr = gradZa(pointZa,pointZS,KpointZa,FpointZa,ZSZa,NULL);CHKERRQ(ierr);
					ierr = IGAPointAddMat(pointZa,KpointZa,KlocZa);CHKERRQ(ierr);
					ierr = IGAPointAddVec(pointZa,FpointZa,FlocZa);CHKERRQ(ierr);
				}
				IGAElementNextPoint(elemZS,pointZS);
				ierr = IGAElementEndPoint(elemZa,&pointZa);CHKERRQ(ierr);
				ierr = IGAElementEndPoint(elemZS,&pointZS);CHKERRQ(ierr);
			}
		}
		//ierr = IGAElementFixSystem(elemZa,KlocZa,FlocZa);CHKERRQ(ierr);
		ierr = IGAElementAssembleMat(elemZa,KlocZa,KZa);CHKERRQ(ierr);
		ierr = IGAElementAssembleVec(elemZa,FlocZa,FZa);CHKERRQ(ierr);

	}
	IGANextElement(igaZS,elemZS);
	ierr = IGAEndElement(igaZa,&elemZa);CHKERRQ(ierr);
	ierr = IGAEndElement(igaZS,&elemZS);CHKERRQ(ierr);

	// Restore local vectors s0 and arrays
	ierr = IGARestoreLocalVecArray(igaZS,ZS0,&localZSZa,&arrayZSZa);CHKERRQ(ierr);
	
	//Impose Dirichlet condition on a single point
	//First we replace rows and columns with associated rows and columns from identity matrix
	PetscInt ma,na;
	ierr = MatGetSize(KZa,&na,&ma);CHKERRQ(ierr);
	ierr = MatSetOption(KZa,MAT_NEW_NONZERO_ALLOCATION_ERR,PETSC_FALSE);CHKERRQ(ierr);
	ierr = MatAssemblyBegin(KZa,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyEnd  (KZa,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);	

	PetscInt rowsZa, colsZa[ma];
	PetscReal valZa, valsZa[ma];
	for(int i=0;i<ma;i++)
	{
		colsZa[i]=i;
		valsZa[i]=0.0;
	}

	rowsZa=0;
	valsZa[rowsZa]=1.0e6;
	ierr = MatSetValues(KZa,1,&rowsZa,ma,colsZa,valsZa,INSERT_VALUES);CHKERRQ(ierr);
	ierr = MatSetValues(KZa,na,colsZa,1,&rowsZa,valsZa,INSERT_VALUES);CHKERRQ(ierr);
	valsZa[rowsZa]=0.0;
	rowsZa=1;
	valsZa[rowsZa]=1.0e6;
	ierr = MatSetValues(KZa,1,&rowsZa,ma,colsZa,valsZa,INSERT_VALUES);CHKERRQ(ierr);
	ierr = MatSetValues(KZa,na,colsZa,1,&rowsZa,valsZa,INSERT_VALUES);CHKERRQ(ierr);
	valsZa[rowsZa]=0.0;

	ierr = MatAssemblyBegin(KZa,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyEnd  (KZa,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

	//Here we set values to the vector directly
	ierr = VecAssemblyBegin(FZa);CHKERRQ(ierr);
	ierr = VecAssemblyEnd  (FZa);CHKERRQ(ierr);

	rowsZa=0;
	valZa=0.0;
	ierr = VecSetValue(FZa,rowsZa,valZa,INSERT_VALUES);CHKERRQ(ierr);
	rowsZa=1;
	valZa=0.0;
	ierr = VecSetValue(FZa,rowsZa,valZa,INSERT_VALUES);CHKERRQ(ierr);

	ierr = VecAssemblyBegin(FZa);CHKERRQ(ierr);
	ierr = VecAssemblyEnd  (FZa);CHKERRQ(ierr);
	//

	ierr = KSPSetOperators(kspZa,KZa,KZa);CHKERRQ(ierr);
	PC pcZa;
	ierr = KSPGetPC(kspZa,&pcZa); CHKERRQ(ierr);
	ierr = PCSetType(pcZa,PCLU); CHKERRQ(ierr);
	ierr = PCFactorSetMatSolverType(pcZa,MATSOLVERMUMPS); CHKERRQ(ierr);
	//ierr = KSPSetFromOptions(kspZa);CHKERRQ(ierr);
	//ierr = KSPSetTolerances(kspZa,1.0e-16,1.0e-20,PETSC_DEFAULT,PETSC_DEFAULT);CHKERRQ(ierr);
	ierr = KSPSolve(kspZa,FZa,Za0);CHKERRQ(ierr);

	ierr = KSPDestroy(&kspZa);CHKERRQ(ierr);
	ierr = MatDestroy(&KZa);CHKERRQ(ierr);
	ierr = VecDestroy(&FZa);CHKERRQ(ierr);
	
	char nameZa[]="/Za-2d-0.dat";
	char pathZa[512];
	sprintf(pathZa,"%s%s",direct,nameZa);
	ierr = IGAWriteVec(igaZa,Za0,pathZa);CHKERRQ(ierr);
//
*/

/*
//System for L2 projection of GradGradz
	PetscPrintf(PETSC_COMM_WORLD,"\n System for GradGradz starting \n\n");
	IGA igaGradGradz;
	ierr = IGACreate(PETSC_COMM_WORLD,&igaGradGradz);CHKERRQ(ierr);
	ierr = IGASetDim(igaGradGradz,2);CHKERRQ(ierr);													//Spatial dimension of the problem
	ierr = IGASetDof(igaGradGradz,8);CHKERRQ(ierr);													//Number of degrees of freedom, per node
	ierr = IGASetOrder(igaGradGradz,3);CHKERRQ(ierr);													//Number of spatial derivatives to calculate
	ierr = IGASetFromOptions(igaGradGradz);CHKERRQ(ierr);												//Note: The order (or degree) of the shape functions is given by the mesh!
	ierr = IGARead(igaGradGradz,"./geometry3.dat");CHKERRQ(ierr);
	ierr = IGASetUp(igaGradGradz);CHKERRQ(ierr);

	for (dir=0; dir<2; dir++) 
	{
		for (side=0; side<2; side++) 
		{
			//ierr = IGASetBoundaryValue(iga,dir,side,dof,0.0);CHKERRQ(ierr);					// Dirichlet boundary conditions
			ierr = IGASetBoundaryForm(igaGradGradz,dir,side,PETSC_TRUE);CHKERRQ(ierr);  				// Neumann boundary conditions
		}
	}

	Mat KGradGradz;
	Vec GradGradz0,FGradGradz;

	ierr = IGACreateMat(igaGradGradz,&KGradGradz);CHKERRQ(ierr);
	ierr = IGACreateVec(igaGradGradz,&GradGradz0);CHKERRQ(ierr);
	ierr = IGACreateVec(igaGradGradz,&FGradGradz);CHKERRQ(ierr);

	IGAPoint		pointGradGradz;															//point
	IGAElement		elemGradGradz;																//element
	PetscReal		*KlocGradGradz,*FlocGradGradz;												//AA y BB
	PetscReal		*KpointGradGradz,*FpointGradGradz;											//KKK y FFF
	const PetscReal	*arrayZ0GradGradz;			//arrayU
	Vec				localZ0GradGradz;				//localU
	PetscReal		*Z0GradGradz;											//U

  	IGAFormSystem	wtfGradGradz;
 	void			*wtf2GradGradz;

 	KSP kspGradGradz;
	ierr = IGACreateKSP(igaGradGradz,&kspGradGradz);CHKERRQ(ierr);

	// Get local vectors Z0 and Chi0 and arrays
	ierr = IGAGetLocalVecArray(igaZ0,Z0,&localZ0GradGradz,&arrayZ0GradGradz);CHKERRQ(ierr);

	// Element loop
	ierr = IGABeginElement(igaGradGradz,&elemGradGradz);CHKERRQ(ierr);
	ierr = IGABeginElement(igaZ0,&elemZ0);CHKERRQ(ierr);

	while (IGANextElement(igaGradGradz,elemGradGradz)) 
	{
		IGANextElement(igaZ0,elemZ0);

		ierr = IGAElementGetWorkMat(elemGradGradz,&KlocGradGradz);CHKERRQ(ierr);
		ierr = IGAElementGetWorkVec(elemGradGradz,&FlocGradGradz);CHKERRQ(ierr);
		ierr = IGAElementGetValues(elemZ0,arrayZ0GradGradz,&Z0GradGradz);CHKERRQ(ierr);

		// FormSystem loop
		while (IGAElementNextFormSystem(elemGradGradz,&wtfGradGradz,&wtf2GradGradz)) 
		{
		// Quadrature loop
			ierr = IGAElementBeginPoint(elemGradGradz,&pointGradGradz);CHKERRQ(ierr);
			ierr = IGAElementBeginPoint(elemZ0,&pointZ0);CHKERRQ(ierr);

			while (IGAElementNextPoint(elemGradGradz,pointGradGradz))
			{
				if(pointGradGradz->atboundary==1)
				{
					//IGAElementNextPoint(elemZ0,pointZ0);
					//IGAElementNextPoint(elemchiUp,pointchiUp);
				}

				if(pointZ0->atboundary==1)
				{
					PetscPrintf(PETSC_COMM_WORLD,"Hola1 en GradGradZ \n");
				}

				if(pointGradGradz->atboundary==0 && pointZ0->atboundary==0)
				{
					IGAElementNextPoint(elemZ0,pointZ0);

					ierr = IGAPointGetWorkMat(pointGradGradz,&KpointGradGradz);CHKERRQ(ierr);
					ierr = IGAPointGetWorkVec(pointGradGradz,&FpointGradGradz);CHKERRQ(ierr);
					//	   SysGradGradz(IGAPoint p,IGAPoint pZu,PetscReal *K,PetscReal *F,PetscReal *Zu,void *ctx)
					ierr = SysGradGradz(pointGradGradz,pointZ0,KpointGradGradz,FpointGradGradz,Z0GradGradz,NULL);CHKERRQ(ierr);
					ierr = IGAPointAddMat(pointGradGradz,KpointGradGradz,KlocGradGradz);CHKERRQ(ierr);
					ierr = IGAPointAddVec(pointGradGradz,FpointGradGradz,FlocGradGradz);CHKERRQ(ierr);
				}
			}
			if (pointZ0->index != -1)
			{
				IGAElementNextPoint(elemZ0,pointZ0);
			}

			ierr = IGAElementEndPoint(elemGradGradz,&pointGradGradz);CHKERRQ(ierr);
			ierr = IGAElementEndPoint(elemZ0,&pointZ0);CHKERRQ(ierr);
		}

		//ierr = IGAElementFixSystem(elemStress,KlocStress,FlocStress);CHKERRQ(ierr);					//This sets Dirichlet condition ¿? (Yes, this applies the conditions from IGASetBoundaryValue)
		ierr = IGAElementAssembleMat(elemGradGradz,KlocGradGradz,KGradGradz);CHKERRQ(ierr);
		ierr = IGAElementAssembleVec(elemGradGradz,FlocGradGradz,FGradGradz);CHKERRQ(ierr);
	}
	IGANextElement(igaZ0,elemZ0);

	ierr = IGAEndElement(igaGradGradz,&elemGradGradz);CHKERRQ(ierr);
	ierr = IGAEndElement(igaZ0,&elemZ0);CHKERRQ(ierr);

	// Restore local vectors u, Z0, Chi0 and arrays
	ierr = IGARestoreLocalVecArray(igaZ0,Z0,&localZ0GradGradz,&arrayZ0GradGradz);CHKERRQ(ierr);

	ierr = MatAssemblyBegin(KGradGradz,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyEnd  (KGradGradz,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

	ierr = VecAssemblyBegin(FGradGradz);CHKERRQ(ierr);
	ierr = VecAssemblyEnd  (FGradGradz);CHKERRQ(ierr);

	ierr = KSPSetOperators(kspGradGradz,KGradGradz,KGradGradz);CHKERRQ(ierr);
	//PC pcStress;
	//ierr = KSPGetPC(kspStress,&pcStress); CHKERRQ(ierr);
	//ierr = PCSetType(pcStress,PCLU); CHKERRQ(ierr);
	//ierr = PCFactorSetMatSolverType(pcStress,MATSOLVERMUMPS); CHKERRQ(ierr);
	ierr = KSPSetFromOptions(kspGradGradz);CHKERRQ(ierr);
	ierr = KSPSetTolerances(kspGradGradz,1e-24,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);CHKERRQ(ierr);
	ierr = KSPSolve(kspGradGradz,FGradGradz,GradGradz0);CHKERRQ(ierr);

	ierr = KSPDestroy(&kspGradGradz);CHKERRQ(ierr);
	ierr = MatDestroy(&KGradGradz);CHKERRQ(ierr);
	ierr = VecDestroy(&FGradGradz);CHKERRQ(ierr);

	char nameGradGradz[]="/GradGradz-2d-0.dat";
	char pathGradGradz[512];
	sprintf(pathGradGradz,"%s%s",direct,nameGradGradz);
	ierr = IGAWriteVec(igaGradGradz,GradGradz0,pathGradGradz);CHKERRQ(ierr);	
//
*/

//Destroy all objects not needed anymore (Better to do it here in case different codes call the same IGA, move if memory is a problem)
	ierr = IGADestroy(&igaAl); CHKERRQ(ierr);
	ierr = IGADestroy(&igaZ0); CHKERRQ(ierr);
	ierr = IGADestroy(&igachiUp); CHKERRQ(ierr);
	ierr = IGADestroy(&igaGrad); CHKERRQ(ierr);
	//ierr = IGADestroy(&igaExact);CHKERRQ(ierr);
//

ierr = PetscFinalize();CHKERRQ(ierr);

return 0;
}