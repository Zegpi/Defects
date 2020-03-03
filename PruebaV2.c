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
		PetscReal dx=Lx/51.0;
		PetscReal dy=Ly/51.0;
		PetscReal eps=0.075;
		PetscReal eps2=0.08;
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
		
		PetscReal norm[3]={0.0,1.0,0.0};
		PetscReal sTemp[3][3][3]={0};
		PetscReal A[3][3]={0};

		/*
		if(x[0]>=-dx/2.0)
		{
			//a=1.0/(2.0*dx)*(x[0]-dx);
			a=0.5*(0.5-tanh((x[0]-dx/2.0))/eps2);
		}
		else if((x[0]<-dx) && (x[1]>-dy && x[1]<dy))
		{
			a=1.0;
		}
		else
		{
			a=0.0;
		}
		*/

		a2=0.5*(1.0-tanh((x[0]-dx)/eps2))/dy;

		A[0][0]=0.0;
		A[0][1]=0.5*(tanh((x[1]+dy)/eps)-tanh((x[1]-dy)/eps));    //This is for the through or terminating twin;
		A[0][2]=0.0;
		A[1][0]=0.0;
		A[1][1]=0.0;
		A[1][2]=0.0;
		A[2][0]=0.0;
		A[2][1]=0.0;
		A[2][2]=0.0;


		/*
		This for the single discl.
		A[0][0]=0.0;
		A[0][1]=0.0*tan(5.0/180.0*ConstPi)*0.5*(tanh((x[1]+dy)/eps)-tanh((x[1]-dy)/eps))*a2;
		A[0][2]=0.0
		A[1][0]=0.0*tan(5.0/180.0*ConstPi)*0.5*(tanh((x[1]+dy)/eps)-tanh((x[1]-dy)/eps))*a2;
		A[1][1]=0.0;
		A[1][2]=0.0;
		A[2][0]=0.0;
		A[2][1]=0.0;
		A[2][2]=0.0;
		*/


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

//Helmholtz decomposition of S, curl part
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

//Helmholtz decomposition of S, grad part
	#undef  __FUNCT__
	#define __FUNCT__ "gradZS"
	PetscErrorCode gradZS(IGAPoint p,IGAPoint pS,PetscReal *K,PetscReal *F,PetscReal *U,void *ctx)
	{

		const PetscReal *N0,(*N1)[2];
		IGAPointGetShapeFuns(p,0,(const PetscReal**)&N0);
		IGAPointGetShapeFuns(p,1,(const PetscReal**)&N1);

		PetscInt a,b,u,w,i,j,k,nen=p->nen, dof=p->dof;

		PetscReal S[8];																		//Create array to recieve Alfa
		//PetscReal dS[8][2];																	//Create array to recieve dAlfa
		IGAPointFormValue(pS,U,&S[0]);															//This fills the values
		//IGAPointFormGrad (pS,U,&dS[0][0]);														//This fills the values of the derivatives

		PetscReal fullS[3][3][3]={0};
		//PetscReal fulldS[3][3][3][3]={0};

		fullS[0][0][0]=S[0]; fullS[0][0][1]=S[1]; 											//Expand S to full 3rd order form, only non-zero elements
		fullS[0][1][0]=S[2]; fullS[0][1][1]=S[3];		
		fullS[1][0][0]=S[4]; fullS[1][0][1]=S[5]; 
		fullS[1][1][0]=S[6]; fullS[1][1][1]=S[7];

		//Same for gradS
		//fulldS[0][0][0][0]=dS[0][0]; fulldS[0][0][0][1]=dS[0][1]; 
		//fulldS[0][0][1][0]=dS[1][0]; fulldS[0][0][1][1]=dS[1][1]; 
		//fulldS[0][1][0][0]=dS[2][0]; fulldS[0][1][0][1]=dS[2][1]; 
		//fulldS[0][1][1][0]=dS[3][0]; fulldS[0][1][1][1]=dS[3][1];
		//fulldS[1][0][0][0]=dS[4][0]; fulldS[1][0][0][1]=dS[4][1]; 
		//fulldS[1][0][1][0]=dS[5][0]; fulldS[1][0][1][1]=dS[5][1]; 
		//fulldS[1][1][0][0]=dS[6][0]; fulldS[1][1][0][1]=dS[6][1]; 
		//fulldS[1][1][1][0]=dS[7][0]; fulldS[1][1][1][1]=dS[7][1];

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
								FZS[a][u]+=fullS[i][j][k]*dv[i][j][k];
							}
						}
					}
				}
			}
		}

		return 0;
	}
//

//L2 projection of Al(0)-Sp:X
	#undef  __FUNCT__
	#define __FUNCT__ "L2ProjectionAlSp"
	PetscErrorCode L2ProjectionAlSp(IGAPoint p,IGAPoint pS,PetscReal *K,PetscReal *F,PetscReal *U,void *ctx)
	{
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

		fullSp[0][0][0]=round(1.0e12*Sp[0])/1.0e12; fullSp[0][0][1]=round(1.0e12*Sp[1])/1.0e12; fullSp[0][1][0]=round(1.0e12*Sp[2])/1.0e12; fullSp[0][1][1]=round(1.0e12*Sp[3])/1.0e12;		//Expand S to full 3rd order form, only non-zero elements
		fullSp[1][0][0]=round(1.0e12*Sp[4])/1.0e12; fullSp[1][0][1]=round(1.0e12*Sp[5])/1.0e12; fullSp[1][1][0]=round(1.0e12*Sp[6])/1.0e12; fullSp[1][1][1]=round(1.0e12*Sp[7])/1.0e12;

		AppCtxL2    *user  = (AppCtxL2 *)ctx;
		PetscReal Lx     = user->Lx;
		PetscReal Ly     = user->Ly;
		PetscReal nx     = user->nx;
		PetscReal ny     = user->ny;

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

		//g is the function to L2 project
		PetscReal dx=Lx/nx;
		PetscReal dy=Ly/ny;
		PetscReal eps=0.01;
		PetscReal g[dof];
		for (i=0; i<dof; i++)
		{
			if (i==0)
			{
				//g[i]=0*1.917545*(0.5*(tanh((x[0]+dx)/eps)+1.0)-0.5*(tanh((x[0]-dx)/eps)+1.0))*(0.5*(tanh((x[1]+dy)/eps)+1.0)-0.5*(tanh((x[1]-dy)/eps)+1.0));
				g[i]=0.0*0.9585551*0.5*(tanh((x[0]+dx)/eps)-tanh((x[0]-dx)/eps))*(tanh((x[1]+dy)/eps)-tanh((x[1]-dy)/eps));
			}
			else if (i==1)
			{
				//g[i]=1.917545*(0.5*(tanh((x[0]+dx)/eps)+1.0)-0.5*(tanh((x[0]-dx)/eps)+1.0))*(0.5*(tanh((x[1]+dy)/eps)+1.0)-0.5*(tanh((x[1]-dy)/eps)+1.0));
				g[i]=0.0*0.9585551*0.5*(tanh((x[0]+dx)/eps)-tanh((x[0]-dx)/eps))*(tanh((x[1]+dy)/eps)-tanh((x[1]-dy)/eps));
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
				FF[a][i] = N[a]*(g[i]-sp[i]);
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

		/*
		PetscReal dzU[4][2];
		IGAPointFormGrad(pZs,Z,&dzU[0][0]);

		PetscReal fulldZ[3][3][3];
		fulldZ[0][0][0]=dzU[0][0]; fulldZ[0][0][1]=dzU[0][1];
		fulldZ[0][1][0]=dzU[1][0]; fulldZ[0][1][1]=dzU[1][1];
		fulldZ[1][0][0]=dzU[2][0]; fulldZ[1][0][1]=dzU[2][1];
		fulldZ[1][1][0]=dzU[3][0]; fulldZ[1][1][1]=dzU[3][1];

		//Calculo de curl(Z)
		PetscReal cZ[3][3];
		
		for(i=0;i<3;i++)
		{
			for(j=0;j<3;j++)
			{
				for(k=0;k<3;k++)
				{
					for(l=0;l<3;l++)
					{
						cZ[i][j]=e[j][k][l]*fulldZ[i][l][k];
					}
				}
			}
		}
		*/

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
									FchiU[a][u]+=-fullAlfa[i][j]*(e[j][k][l]*dv[i][l][k]);
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
	PetscErrorCode Z0sys(IGAPoint p,IGAPoint pChi,IGAPoint pHs,PetscReal *K,PetscReal *F, PetscReal *UChi,PetscReal *UHs,void *ctx)
	{
		PetscReal C[3][3][3][3]={0};														//Initialization of elastic tensor
		const PetscReal *N0,(*N1)[2],(*N2)[2][2];
		IGAPointGetShapeFuns(p,0,(const PetscReal**)&N0);									//Value of the shape functions
		IGAPointGetShapeFuns(p,1,(const PetscReal**)&N1);									//Derivatives of the shape functions
		IGAPointGetShapeFuns(p,2,(const PetscReal**)&N2);									//Second derivatives of the shape funcions
		//After this command Na_xx=N2[a][0][0], Na_yy=N2[a][1][1], Na_xy=N2[a][0][1], Na_yx[a][1][0] (these last two are equal) (remember a is the index of the shape function)
		PetscInt a,b,i,j,k,l,u,w,r,s,m,nen=p->nen, dof=p->dof;

		PetscReal Chi0[4];																	//Assign Al(t) to a vector
		PetscReal dChi0[4][2];																//Same for partial derivatives
		IGAPointFormValue(pChi,UChi,&Chi0[0]);
		IGAPointFormGrad (pChi,UChi,&dChi0[0][0]);

		PetscReal fullChi[3][3];
		fullChi[0][0]=Chi0[0]; 	fullChi[0][1]=Chi0[1]; 	fullChi[0][2]=0;
		fullChi[1][0]=Chi0[2]; 	fullChi[1][1]=Chi0[3]; 	fullChi[1][2]=0;
		fullChi[2][0]=0; 		fullChi[2][1]=0; 		fullChi[2][2]=0;

		PetscReal full_dChi[3][3][3]={0};
		full_dChi[0][0][0]=dChi0[0][0]; full_dChi[0][0][1]=dChi0[0][1];
		full_dChi[0][1][0]=dChi0[1][0]; full_dChi[0][1][1]=dChi0[1][1];
		full_dChi[1][0][0]=dChi0[2][0]; full_dChi[1][0][1]=dChi0[2][1];
		full_dChi[1][1][0]=dChi0[3][0]; full_dChi[1][1][1]=dChi0[3][1];

		PetscReal Hs0[4];
		PetscReal dHs0[4][2];
		IGAPointFormValue(pHs,UHs,&Hs0[0]);
		IGAPointFormGrad (pHs,UHs,&dHs0[0][0]);

		/*
		PetscReal fullHs[3][3];
		fullHs[0][0]=Hs0[0]; 	fullHs[0][1]=Hs0[1]; 	fullHs[0][2]=0;
		fullHs[1][0]=Hs0[2]; 	fullHs[1][1]=Hs0[3]; 	fullHs[1][2]=0;
		fullHs[2][0]=0; 		fullHs[2][1]=0; 		fullHs[2][2]=0;

		PetscReal full_dHs[3][3][3]={0};
		full_dHs[0][0][0]=dHs0[0][0]; full_dHs[0][0][1]=dHs0[0][1];
		full_dHs[0][1][0]=dHs0[1][0]; full_dHs[0][1][1]=dHs0[1][1];
		full_dHs[1][0][0]=dHs0[2][0]; full_dHs[1][0][1]=dHs0[2][1];
		full_dHs[1][1][0]=dHs0[3][0]; full_dHs[1][1][1]=dHs0[3][1];
		*/

		PetscReal (*Keq)[dof][nen][dof] = (typeof(Keq)) K;
		PetscReal (*Feq)[dof] = (PetscReal (*)[dof])F;

		//E and nu should come from AppCtx in the future
		const PetscReal E=2100000.0*9.81*10000.0;
		const PetscReal nu=0.33;
		const PetscReal lambda=(E*nu)/((1.0+nu)*(1.0-2.0*nu));
		const PetscReal mu=E/(2.0*(1.0+nu));
		const PetscReal eps=0.0*E/1000.0;							//Choose later based on whatever Amit says :)

		//////////////Delete this parte later, loads should come from appCtx or be 0
		//PetscReal f[2]={0.0, 0.0};			//Distributed load in body
		//PetscReal g[2]={0.0, 0.0};			//Boundary load (applied wherever is defined in the p->atboundary block)
		//PetscReal M[3]={0.0, 0.0, 0.0};		//Distributed load in body
		//PetscReal h[3]={0.0, 0.0, E};		//Boundary moment (applied wherever is defined in the p->atboundary block)
		/////////////////////////////////////
		/*
		const PetscReal e[3][3][3]=
		{
			{{0.0,0.0,0.0},{0.0,0.0,1.0},{0.0,-1.0,0.0}},
			{{0.0,0.0,-1.0},{0.0,0.0,0.0},{1.0,0.0,0.0}},
			{{0.0,1.0,0.0},{-1.0,0.0,0.0},{0.0,0.0,0.0}}
		};
		*/

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

			PetscReal x[2];																		//Vector of reals, size equal to problem's dimension
			IGAPointFormGeomMap(p,x);															//Fills x with the coordinates of p, Gauss's point

			PetscReal Omega=tan(5.0/180.0*ConstPi);
			PetscReal rho=sqrt(x[0]*x[0]+x[1]*x[1]);
			PetscReal Sborde[3][3]={0};

			Sborde[0][0]=0.0*mu*Omega/(2.0*ConstPi*(1.0-nu))*(log(rho)+(x[1]*x[1])/(rho*rho)+nu/(1.0-2.0*nu));
			Sborde[0][1]=0.0*-mu*Omega/(2.0*ConstPi*(1.0-nu))*x[0]*x[1]/(rho*rho);
			Sborde[1][0]=0.0*-mu*Omega/(2.0*ConstPi*(1.0-nu))*x[0]*x[1]/(rho*rho);
			Sborde[1][1]=0.0*mu*Omega/(2.0*ConstPi*(1.0-nu))*(log(rho)+(x[0]*x[0])/(rho*rho)+nu/(1.0-2.0*nu));

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
							dv[1][0]=0.0; dv[1][1]=0.0; dv[1][2]=0.0;
							dv[2][0]=0.0; dv[2][1]=0.0; dv[2][2]=0.0;

							d2v[0][0][0]=Na_xx; d2v[0][0][1]=Na_xy; d2v[0][0][2]=0.0; d2v[0][1][0]=Na_yx; d2v[0][1][1]=Na_yy; d2v[0][1][2]=0.0; d2v[0][2][0]=0.0; d2v[0][2][1]=0.0; d2v[0][2][2]=0.0;
							d2v[1][0][0]=0.00;  d2v[1][0][1]=0.00;  d2v[1][0][2]=0.0; d2v[1][1][0]=0.0;   d2v[1][1][1]=0.0;   d2v[1][1][2]=0.0; d2v[1][2][0]=0.0; d2v[1][2][1]=0.0; d2v[1][2][2]=0.0;
							d2v[2][0][0]=0.00;  d2v[2][0][1]=0.00;  d2v[2][0][2]=0.0; d2v[2][1][0]=0.0;   d2v[2][1][1]=0.0;   d2v[2][1][2]=0.0; d2v[2][2][0]=0.0; d2v[2][2][1]=0.0; d2v[2][2][2]=0.0;
							
						}
						else if(i==1)
						{
							dv[0][0]=0.0; dv[0][1]=0.0; dv[0][2]=0.0;
							dv[1][0]=Na_x; dv[1][1]=Na_y; dv[1][2]=0.0;
							dv[2][0]=0.0; dv[2][1]=0.0; dv[2][2]=0.0;

							d2v[0][0][0]=0.00;  d2v[0][0][1]=0.00;  d2v[0][0][2]=0.0; d2v[0][1][0]=0.0;   d2v[0][1][1]=0.0;   d2v[0][1][2]=0.0; d2v[0][2][0]=0.0; d2v[0][2][1]=0.0; d2v[0][2][2]=0.0;
							d2v[1][0][0]=Na_xx; d2v[1][0][1]=Na_xy; d2v[1][0][2]=0.0; d2v[1][1][0]=Na_yx; d2v[1][1][1]=Na_yy; d2v[1][1][2]=0.0; d2v[1][2][0]=0.0; d2v[1][2][1]=0.0; d2v[1][2][2]=0.0;
							d2v[2][0][0]=0.00;  d2v[2][0][1]=0.00;  d2v[2][0][2]=0.0; d2v[2][1][0]=0.0;   d2v[2][1][1]=0.0;   d2v[2][1][2]=0.0; d2v[2][2][0]=0.0; d2v[2][2][1]=0.0; d2v[2][2][2]=0.0;
							
						}

						for (j=0; j<dof; j++)
						{
							if (j==0)
							{
								dz[0][0]=Nb_x; dz[0][1]=Nb_y; dz[0][2]=0.0;
								dz[1][0]=0.0; dz[1][1]=0.0; dz[1][2]=0.0;
								dz[2][0]=0.0; dz[2][1]=0.0; dz[2][2]=0.0;

								d2z[0][0][0]=Nb_xx; d2z[0][0][1]=Nb_xy; d2z[0][0][2]=0.0; d2z[0][1][0]=Nb_yx; d2z[0][1][1]=Nb_yy; d2z[0][1][2]=0.0; d2z[0][2][0]=0.0; d2z[0][2][1]=0.0; d2z[0][2][2]=0.0;
								d2z[1][0][0]=0.00;  d2z[1][0][1]=0.00;  d2z[1][0][2]=0.0; d2z[1][1][0]=0.0;   d2z[1][1][1]=0.0;   d2z[1][1][2]=0.0; d2z[1][2][0]=0.0; d2z[1][2][1]=0.0; d2z[1][2][2]=0.0;
								d2z[2][0][0]=0.00;  d2z[2][0][1]=0.00;  d2z[2][0][2]=0.0; d2z[2][1][0]=0.0;   d2z[2][1][1]=0.0;   d2z[2][1][2]=0.0; d2z[2][2][0]=0.0; d2z[2][2][1]=0.0; d2z[2][2][2]=0.0;
							}
							else if(j==1)
							{
								dz[0][0]=0.0; dz[0][1]=0.0; dz[0][2]=0.0;
								dz[1][0]=Nb_x; dz[1][1]=Nb_y; dz[1][2]=0.0;
								dz[2][0]=0.0; dz[2][1]=0.0; dz[2][2]=0.0;

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
											Keq[a][i][b][j]+=-C[k][l][u][w]*dz[u][w]*0.5*(dv[k][l]+dv[l][k]);
										}
									}
								}
							}

							for (r=0; r<3; r++)
							{
								for (s=0; s<3; s++)
								{
									for (m=0; m<3; m++)
									{
										Keq[a][i][b][j]+=-0.25*eps*(d2z[r][s][m]+d2z[r][m][s]+d2z[s][r][m]+d2z[s][m][r])*0.5*(d2v[r][s][m]+d2v[s][r][m])
														 -0.25*eps*(d2z[r][s][m]+d2z[r][m][s]-d2z[s][r][m]-d2z[s][m][r])*d2v[r][s][m];
										//Keq[a][i][b][j]-=0.25*eps*(d2z[r][s][m]+d2z[r][m][s]+d2z[s][r][m]+d2z[s][m][r])*0.5*(d2v[r][s][m]+2.0*d2v[s][r][m]);
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
					PetscReal Na_x=N1[a][0];		PetscReal Na_y=N1[a][1];
					PetscReal Na_xx =N2[a][0][0];	PetscReal Na_yy =N2[a][1][1];
					PetscReal Na_xy =N2[a][0][1];	PetscReal Na_yx =N2[a][1][0];

					if (i==0)
					{
						dv[0][0]=Na_x; dv[0][1]=Na_y; dv[0][2]=0.0;
						dv[1][0]=0.0; dv[1][1]=0.0; dv[1][2]=0.0;
						dv[2][0]=0.0; dv[2][1]=0.0; dv[2][2]=0.0;

						d2v[0][0][0]=Na_xx; d2v[0][0][1]=Na_xy; d2v[0][0][2]=0.0; d2v[0][1][0]=Na_yx; d2v[0][1][1]=Na_yy; d2v[0][1][2]=0.0; d2v[0][2][0]=0.0; d2v[0][2][1]=0.0; d2v[0][2][2]=0.0;
						d2v[1][0][0]=0.00;  d2v[1][0][1]=0.00;  d2v[1][0][2]=0.0; d2v[1][1][0]=0.0;   d2v[1][1][1]=0.0;   d2v[1][1][2]=0.0; d2v[1][2][0]=0.0; d2v[1][2][1]=0.0; d2v[1][2][2]=0.0;
						d2v[2][0][0]=0.00;  d2v[2][0][1]=0.00;  d2v[2][0][2]=0.0; d2v[2][1][0]=0.0;   d2v[2][1][1]=0.0;   d2v[2][1][2]=0.0; d2v[2][2][0]=0.0; d2v[2][2][1]=0.0; d2v[2][2][2]=0.0;
					}
					else if (i==1)
					{
						dv[0][0]=0.0; dv[0][1]=0.0; dv[0][2]=0.0;
						dv[1][0]=Na_x; dv[1][1]=Na_y; dv[1][2]=0.0;
						dv[2][0]=0.0; dv[2][1]=0.0; dv[2][2]=0.0;

						d2v[0][0][0]=0.00;  d2v[0][0][1]=0.00;  d2v[0][0][2]=0.0; d2v[0][1][0]=0.0;   d2v[0][1][1]=0.0;   d2v[0][1][2]=0.0; d2v[0][2][0]=0.0; d2v[0][2][1]=0.0; d2v[0][2][2]=0.0;
						d2v[1][0][0]=Na_xx; d2v[1][0][1]=Na_xy; d2v[1][0][2]=0.0; d2v[1][1][0]=Na_yx; d2v[1][1][1]=Na_yy; d2v[1][1][2]=0.0; d2v[1][2][0]=0.0; d2v[1][2][1]=0.0; d2v[1][2][2]=0.0;
						d2v[2][0][0]=0.00;  d2v[2][0][1]=0.00;  d2v[2][0][2]=0.0; d2v[2][1][0]=0.0;   d2v[2][1][1]=0.0;   d2v[2][1][2]=0.0; d2v[2][2][0]=0.0; d2v[2][2][1]=0.0; d2v[2][2][2]=0.0;
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
									Feq[a][i]+=C[k][l][u][w]*fullChi[u][w]*0.5*(dv[k][l]+dv[l][k]);
								}
							}
						}
					}

					for (r=0; r<3; r++)
					{
						for (s=0; s<3; s++)
						{
							for (m=0; m<3; m++)
							{
								Feq[a][i]+=0.25*eps*(full_dChi[r][s][m]+full_dChi[r][m][s]+full_dChi[s][r][m]+full_dChi[s][m][r])*0.5*(d2v[r][s][m]+d2v[s][r][m])
								          +0.25*eps*(full_dChi[r][s][m]+full_dChi[r][m][s]-full_dChi[s][r][m]-full_dChi[s][m][r])*d2v[r][s][m];
								//Feq[a][i]+=0.25*eps*(full_dChi[r][s][m]+full_dChi[r][m][s]+full_dChi[s][r][m]+full_dChi[s][m][r])*0.5*(d2v[r][s][m]+2.0*d2v[s][r][m]);
							}
						}
					}
				}
			}
		}
		
		return 0;
	}
//

//System for stresses
	#undef  __FUNCT__
	#define __FUNCT__ "Stress"
	//PetscErrorCode Stress(IGAPoint p,IGAPoint pU, IGAPoint pHs,IGAPoint pChi,IGAPoint pZu,PetscReal *K,PetscReal *F,PetscReal *U,PetscReal *HS, PetscReal *Chi,PetscReal *Zu,void *ctx)		//When other fields like Pi or S (etc.) come into this, declare a IGAPoint pPi or pS and a new PescReal *UPi or *US for each
	PetscErrorCode Stress(IGAPoint p,IGAPoint pChi,IGAPoint pZu,PetscReal *K,PetscReal *F, PetscReal *Chi,PetscReal *Zu,void *ctx)		//When other fields like Pi or S (etc.) come into this, declare a IGAPoint pPi or pS and a new PestcReal *UPi or *US for each
	{
		const PetscReal *N0,(*N1)[2],(*N2)[2][2];
		IGAPointGetShapeFuns(p,0,(const PetscReal**)&N0);									//Value of the shape functions
		IGAPointGetShapeFuns(p,1,(const PetscReal**)&N1);									//Derivatives of the shape functions
		IGAPointGetShapeFuns(p,2,(const PetscReal**)&N2);									//Second derivatives of the shape funcions
		PetscInt a,b,i,j,k,l,m,n,nen=p->nen, dof=p->dof;

		/*
		PetscReal u[2];																	//Array to contain vector u
		PetscReal d_u[2][2];															//Same for its gradient
		PetscReal d2_u[2][2][2];													//Same for its 2nd order partial derivatives
		IGAPointFormValue(pU,U,&u[0]);														//Assign u to its container
		IGAPointFormGrad (pU,U,&d_u[0][0]);													//Same for the gradient
		IGAPointFormHess (pU,U,&d2_u[0][0][0]);												//Same for the hessian
		*/

		/*
		PetscReal Hs[4];																	//Array to contain the vector Hs(0)
		PetscReal d_Hs[4][2];																//Same for its gradient
		IGAPointFormValue(pHs,HS,&Hs[0]);													//Assign Hs to its container
		IGAPointFormGrad (pHs,HS,&d_Hs[0][0]);												//Same for the gradient
		*/

		PetscReal chi0[4];																	//Array to contain the vector chi(0)
		PetscReal d_Chi0[4][2];																//Same for its gradient
		IGAPointFormValue(pChi,Chi,&chi0[0]);												//Assign chi to its container
		IGAPointFormGrad(pChi,Chi,&d_Chi0[0][0]);											//Assign grad chi to its container
		//IGAPointFormHess (pChi,Chi,&d2_Chi0[0][0][0]);										//This should be the 3-rd order tensor Chi_{i,jk} (remember that we are storing Chi_{kl} as a column vector)

		//PetscReal Z[2];																	//Array to contain the vector z(0)
		PetscReal d_Z0[2][2];																//Same for its gradient
		PetscReal d2_Z0[2][2][2];															//Same for its 2nd order partial derivatives
		//IGAPointFormValue(pHs,Hs,&Hs[0]);													//Assign z to its container  (Z is not used, only its derivatives)
		IGAPointFormGrad (pZu,Zu,&d_Z0[0][0]);												//Same for the gradient
		IGAPointFormHess (pZu,Zu,&d2_Z0[0][0][0]);											//Same for the hessian

		//Expanding Chi so we can use for loops with correct indexing
		PetscReal fullChi[3][3]={0};
		fullChi[0][0]=chi0[0]; 	fullChi[0][1]=chi0[1];
		fullChi[1][0]=chi0[2]; 	fullChi[1][1]=chi0[3];

		PetscReal fulld_Chi[3][3][3]={0};
		//The four non-zero components of Chi are stored as a vector, restore them to an array with the correct indexing for value and derivative
		fulld_Chi[0][0][0]=d_Chi0[0][0]; fulld_Chi[0][0][1]=d_Chi0[0][1];
		fulld_Chi[0][1][0]=d_Chi0[1][0]; fulld_Chi[0][1][1]=d_Chi0[1][1];
		fulld_Chi[1][0][0]=d_Chi0[2][0]; fulld_Chi[1][0][1]=d_Chi0[2][1];
		fulld_Chi[1][1][0]=d_Chi0[3][0]; fulld_Chi[1][1][1]=d_Chi0[3][1];

		//Expanding z (and derivatives) to 3 components, more convenient for sums in for loops

		PetscReal fulld_z[3][3]={0};
		fulld_z[0][0]=d_Z0[0][0]; fulld_z[0][1]=d_Z0[0][1];
		fulld_z[1][0]=d_Z0[1][0]; fulld_z[1][1]=d_Z0[1][1];
		
		PetscReal fulld2_z[3][3][3]={0};
		for(i=0;i<2;i++)
		{
			for(j=0;j<2;j++)
			{
				for(k=0;k<2;k++)
				{
					fulld2_z[i][j][k]=d2_Z0[i][j][k];
				}
			}
		}
		//fulld2_z[0][0][0]=d2_Z0[0][0][0]; fulld2_z[0][0][1]=d2_Z0[0][0][1];
		//fulld2_z[0][1][0]=d2_Z0[0][1][0]; fulld2_z[0][1][1]=d2_Z0[0][1][1];
		//fulld2_z[1][0][0]=d2_Z0[1][0][0]; fulld2_z[1][0][1]=d2_Z0[1][0][1];
		//fulld2_z[1][1][0]=d2_Z0[1][1][0]; fulld2_z[1][1][1]=d2_Z0[1][1][1];

		PetscReal (*Kstress)[dof][nen][dof] = (typeof(Kstress)) K;
		PetscReal (*Fstress)[dof] = (PetscReal (*)[dof])F;

		//E and nu should come from AppCtx in the future
		const PetscReal E=2100000.0*9.81*10000.0;
		const PetscReal nu=0.33;
		const PetscReal lambda=(E*nu)/((1.0+nu)*(1.0-2.0*nu));
		const PetscReal mu=E/(2.0*(1.0+nu));
		const PetscReal eps=0.0*E/1000.0;														//Choose later based on whatever Amit says :)

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
						m=0; n=0;
					}
					else if (i==1)
					{
						m=0; n=1;
					}
					else if (i==2)
					{
						m=1; n=0;
					}
					else if (i==3)
					{
						m=1; n=1;
					}
					
					for (k=0;k<3;k++)
					{
						for(l=0;l<3;l++)
						{
							Fstress[a][i]+=N0[a]*C[m][n][k][l]*(-fullChi[k][l]-fulld_z[k][l]);
						}
						
						Fstress[a][i]+=N1[a][k]*0.25*eps*(-fulld_Chi[m][n][k]-fulld2_z[m][n][k]-fulld_Chi[m][k][n]-fulld2_z[m][k][n]
														  -fulld_Chi[n][m][k]-fulld2_z[n][m][k]-fulld_Chi[n][k][m]-fulld2_z[n][k][m]);
					}
				}
			}
		}
		return 0;
	}
//

//System for L2 projection of stress for single discl
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

		const PetscReal E=2100000.0*9.81*10000.0;
		const PetscReal nu=0.33;
		const PetscReal mu=E/(2.0*(1.0+nu));
		PetscReal Omega=tan(5.0/180.0*ConstPi);
		PetscReal rho=sqrt(x[0]*x[0]+x[1]*x[1]);

		g[0]=mu*Omega/(2.0*ConstPi*(1.0-nu))*(log(rho)+(x[1]*x[1])/(rho*rho)+nu/(1.0-2.0*nu));
		g[1]=-mu*Omega/(2.0*ConstPi*(1.0-nu))*x[0]*x[1]/(rho*rho);
		g[2]=-mu*Omega/(2.0*ConstPi*(1.0-nu))*x[0]*x[1]/(rho*rho);
		g[3]=mu*Omega/(2.0*ConstPi*(1.0-nu))*(log(rho)+(x[0]*x[0])/(rho*rho)+nu/(1.0-2.0*nu));

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
	PetscPrintf(PETSC_COMM_WORLD,"Start of PruebaV2 \n");

	PetscInt commsize,rank;
	ierr = MPI_Comm_size(PETSC_COMM_WORLD,&commsize);CHKERRQ(ierr);
	ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank);CHKERRQ(ierr);

//

/*
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
*/

//Generate Mesh and copy cpp to result folder to save for reproduction
	PetscReal Lx=20.0;
	PetscReal Ly=20.0;
	PetscInt  nx=51*2;
	PetscInt  ny=51*2;
	PetscReal h=fmin(Lx/nx,Ly/ny);
	
	if (rank==0)
	{
		ierr=system("cp ~/shared/CodigosPetIGA/PruebaV2.c ~/sharedResults/");

		PetscPrintf(PETSC_COMM_WORLD,"h= %f\n",h);

		char      nombreMalla[12]="geometry";
		char 	  comando[512];
		sprintf(comando,"exec python rectangle.py %s %f %f %d %d\n",nombreMalla,Lx,Ly,nx,ny);
		ierr=system(comando);

		strcpy(nombreMalla,"geometry2");
		//nombreMalla="geometry2";
		sprintf(comando,"exec python rectangle2.py %s %f %f %d %d\n",nombreMalla,Lx,Ly,nx,ny);
		ierr=system(comando);

		//Time parameters
		PetscReal T=0.05;
		PetscReal fact=0.5;
		PetscReal dt=fact*(h/2.0);																	//dt must be less than h/2
		PetscPrintf(PETSC_COMM_WORLD,"nx<= %d ny= %d\n",nx,ny);
		PetscPrintf(PETSC_COMM_WORLD,"dt<= %f     dt= %f\n",dt/fact,dt);
		PetscInt numStep=ceil(T/dt);
		PetscPrintf(PETSC_COMM_WORLD,"Nt= %d\n",numStep);
	}
//

//App context creation
	AppCtxL2 userL2;
	userL2.Lx     = Lx;
	userL2.Ly     = Ly;
	userL2.nx     = nx;
	userL2.ny     = ny;
	//void* user;
//

//Prueba creacion matrices
	PetscInt rowsElim=4, colsElim=3;
	PetscInt M,M2,N,N2, rows2elim[rowsElim], cols2elim[colsElim];
	M=10; M2=M-rowsElim;
	N=10; N2=N-colsElim;
	Mat elimRow, elimCols, matPrueba, C, D;
	ierr=MatCreate(PETSC_COMM_WORLD,&elimRow); CHKERRQ(ierr);
	ierr=MatSetSizes(elimRow, PETSC_DECIDE, PETSC_DECIDE, M2, N); CHKERRQ(ierr);
	ierr=MatSetUp(elimRow);

	ierr=MatCreate(PETSC_COMM_WORLD,&elimCols); CHKERRQ(ierr);
	ierr=MatSetSizes(elimCols, PETSC_DECIDE, PETSC_DECIDE, M, N2); CHKERRQ(ierr);
	ierr=MatSetUp(elimCols);

	ierr=MatCreate(PETSC_COMM_WORLD,&matPrueba); CHKERRQ(ierr);
	ierr=MatSetSizes(matPrueba, PETSC_DECIDE, PETSC_DECIDE, M, N); CHKERRQ(ierr);
	ierr=MatSetUp(matPrueba);

	rows2elim[0]=0; rows2elim[1]=1; rows2elim[2]=2; rows2elim[3]=3;
	cols2elim[0]=1; cols2elim[1]=2; cols2elim[2]=3;

	PetscReal value[1];
	value[0]=1.0;

	for(int j=0; j<rowsElim; j++)
	{
		rows2elim[j]=rows2elim[j]-j;
	}

	for(int j=0; j<colsElim; j++)
	{
		cols2elim[j]=cols2elim[j]-j;
	}

	PetscInt cont=0;
	for(int i=0; i<M2; i++)
	{
		cont=i;
		for(int j=0; j<rowsElim; j++)
		{
			if(i>=rows2elim[j])
			{
				cont=cont+1;
			}
		}
		ierr=MatSetValues(elimRow, 1, &i, 1, &cont, &value[0], INSERT_VALUES); CHKERRQ(ierr);
	}

	cont=0;
	for(int i=0; i<N2; i++)
	{
		cont=i;
		for(int j=0; j<colsElim; j++)
		{
			if(i>=cols2elim[j])
			{
				cont=cont+1;
			}
		}
		ierr=MatSetValues(elimCols, 1, &cont, 1, &i, &value[0], INSERT_VALUES); CHKERRQ(ierr);
	}

	PetscReal value2[N*M];
	for (int i=0; i<(N*M); i++)
	{
		value2[i]=i+1.0;
	}

	for (int i=0; i<M; i++)
	{
		for (int j=0; j<N; j++)
		{
			ierr=MatSetValues(matPrueba, 1, &i, 1, &j, &value2[i*N+j], INSERT_VALUES); CHKERRQ(ierr);
		}
	}

	ierr=MatAssemblyBegin(elimRow,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
	ierr=MatAssemblyEnd(elimRow,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);

	ierr=MatAssemblyBegin(elimCols,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
	ierr=MatAssemblyEnd(elimCols,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);

	ierr=MatAssemblyBegin(matPrueba,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
	ierr=MatAssemblyEnd(matPrueba,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);

	PetscViewerPushFormat(PETSC_VIEWER_STDOUT_WORLD,PETSC_VIEWER_DEFAULT);
	MatView(elimRow,PETSC_VIEWER_STDOUT_WORLD);

	PetscViewerPushFormat(PETSC_VIEWER_STDOUT_WORLD,PETSC_VIEWER_DEFAULT);
	MatView(matPrueba,PETSC_VIEWER_STDOUT_WORLD);

	ierr=MatMatMult(elimRow,matPrueba,MAT_INITIAL_MATRIX,1, &C);
	ierr=MatMatMult(matPrueba,elimCols,MAT_INITIAL_MATRIX,1, &D);

	PetscViewerPushFormat(PETSC_VIEWER_STDOUT_WORLD,PETSC_VIEWER_DEFAULT);
	MatView(C,PETSC_VIEWER_STDOUT_WORLD);

	PetscViewerPushFormat(PETSC_VIEWER_STDOUT_WORLD,PETSC_VIEWER_DEFAULT);
	MatView(D,PETSC_VIEWER_STDOUT_WORLD);
//

/*
//Creation of types and systems for the L2 projection of S0
	PetscPrintf(PETSC_COMM_WORLD,"System for L2 projection for S starting \n");
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
	//ierr = KSPSetFromOptions(kspl2S);CHKERRQ(ierr);
	ierr = KSPSetTolerances(kspl2S,1.0e-16,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);CHKERRQ(ierr);
	ierr = KSPSolve(kspl2S,Fl2S,s0);CHKERRQ(ierr);										//This is a simple system, so it can be solved with just this command

	ierr = KSPDestroy(&kspl2S);CHKERRQ(ierr);
	ierr = MatDestroy(&Kl2S);CHKERRQ(ierr);
	ierr = VecDestroy(&Fl2S);CHKERRQ(ierr);
	char nameS[]="/S-2d-0.dat";
	char pathS[512];
	sprintf(pathS,"%s%s",direct,nameS);
	ierr = IGAWriteVec(igaS,s0,pathS);CHKERRQ(ierr);
//

//Creation of types and systems for the Helmholtz decomposition of S, curl part
	//System for chiS
	PetscPrintf(PETSC_COMM_WORLD,"System for curl part of Helmholtz of S starting \n");
	IGA igachiS;
	ierr = IGACreate(PETSC_COMM_WORLD,&igachiS);CHKERRQ(ierr);
	ierr = IGASetDim(igachiS,2);CHKERRQ(ierr);														//Spatial dimension of the problem
	ierr = IGASetDof(igachiS,8);CHKERRQ(ierr);														//Number of degrees of freedom, per node
	ierr = IGASetOrder(igachiS,1);CHKERRQ(ierr);														//Number of spatial derivatives to calculate
	ierr = IGASetFromOptions(igachiS);CHKERRQ(ierr);													//Note: The order (or degree) of the shape functions is given by the mesh!
	ierr = IGARead(igachiS,"./geometry.dat");CHKERRQ(ierr);
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
				IGAElementNextPoint(elemS,pointS);
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
	//ierr = KSPSetTolerances(kspchiS,1.0e-8,5.0e-18,PETSC_DEFAULT,PETSC_DEFAULT);CHKERRQ(ierr);
	ierr = KSPSolve(kspchiS,FchiS,chiS0);CHKERRQ(ierr);

	ierr = KSPDestroy(&kspchiS);CHKERRQ(ierr);
	ierr = MatDestroy(&KchiS);CHKERRQ(ierr);
	ierr = VecDestroy(&FchiS);CHKERRQ(ierr);
	
	char namechiS[]="/ChiS-2d-0.dat";
	char pathchiS[512];
	sprintf(pathchiS,"%s%s",direct,namechiS);
	ierr = IGAWriteVec(igachiS,chiS0,pathchiS);CHKERRQ(ierr);
//

//Creation of types and systems for the Helmholtz decomposition of S, grad part
	//System for chiS
	PetscPrintf(PETSC_COMM_WORLD,"System for grad part of Helmholtz of S starting \n");
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

	//Form system matrix and vector
	ierr = MatAssemblyBegin(KZS,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyEnd  (KZS,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	
	//
	//PetscInt mz,nz;
	//ierr = MatGetSize(KZS,&mz,&nz);CHKERRQ(ierr);

	PetscInt rows=0;
	PetscReal vals; 

	rows=0;
	ierr = MatZeroRows(KZS,1,&rows,1.0e50,0,0);
	rows=1;
	ierr = MatZeroRows(KZS,1,&rows,1.0e50,0,0);
	rows=2;
	ierr = MatZeroRows(KZS,1,&rows,1.0e50,0,0);
	rows=3;
	ierr = MatZeroRows(KZS,1,&rows,1.0e50,0,0);

	ierr = MatAssemblyBegin(KZS,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyEnd  (KZS,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

	//Here we set values to the vector directly (to impose Dirichlet condition in a single point)
	ierr = VecAssemblyBegin(FZS);CHKERRQ(ierr);
	ierr = VecAssemblyEnd  (FZS);CHKERRQ(ierr);
	
	rows=0;
	vals=0.0;
	ierr = VecSetValue(FZS,rows,vals,INSERT_VALUES);CHKERRQ(ierr);
	rows=1;
	vals=0.0;
	ierr = VecSetValue(FZS,rows,vals,INSERT_VALUES);CHKERRQ(ierr);
	rows=2;
	vals=0.0;
	ierr = VecSetValue(FZS,rows,vals,INSERT_VALUES);CHKERRQ(ierr);
	rows=3;
	vals=0.0;
	ierr = VecSetValue(FZS,rows,vals,INSERT_VALUES);CHKERRQ(ierr);

	ierr = VecAssemblyBegin(FZS);CHKERRQ(ierr);
	ierr = VecAssemblyEnd  (FZS);CHKERRQ(ierr);
	//

	ierr = KSPSetOperators(kspZS,KZS,KZS);CHKERRQ(ierr);
	PC pcZs;
	ierr = KSPGetPC(kspZS,&pcZs); CHKERRQ(ierr);
	ierr = PCSetType(pcZs,PCLU); CHKERRQ(ierr);
	ierr = PCFactorSetMatSolverType(pcZs,MATSOLVERMUMPS); CHKERRQ(ierr);
	//ierr = KSPSetFromOptions(kspZS);CHKERRQ(ierr);
	//ierr = KSPSetTolerances(kspZS,1.0e-9,1.0e-18,PETSC_DEFAULT,PETSC_DEFAULT);CHKERRQ(ierr);
	ierr = KSPSolve(kspZS,FZS,ZS0);CHKERRQ(ierr);

	ierr = KSPDestroy(&kspZS);CHKERRQ(ierr);
	ierr = MatDestroy(&KZS);CHKERRQ(ierr);
	ierr = VecDestroy(&FZS);CHKERRQ(ierr);
	
	char nameZS[]="/ZS-2d-0.dat";
	char pathZS[512];
	sprintf(pathZS,"%s%s",direct,nameZS);
	ierr = IGAWriteVec(igaZS,ZS0,pathZS);CHKERRQ(ierr);
//

//Creation of types and systems for the L2 projection of Alfa0-Sp:X
	//System for Alfa
	PetscPrintf(PETSC_COMM_WORLD,"System for L2 projection for Alfa+Sp:X starting \n");
	IGA igaAl;
	ierr = IGACreate(PETSC_COMM_WORLD,&igaAl);CHKERRQ(ierr);
	ierr = IGASetDim(igaAl,2);CHKERRQ(ierr);														//Spatial dimension of the problem
	ierr = IGASetDof(igaAl,2);CHKERRQ(ierr);														//Number of degrees of freedom, per node
	ierr = IGASetOrder(igaAl,1);CHKERRQ(ierr);													//Number of spatial derivatives to calculate
	ierr = IGASetFromOptions(igaAl);CHKERRQ(ierr);												//Note: The order (or degree) of the shape functions is given by the mesh!
	ierr = IGARead(igaAl,"./geometry.dat");CHKERRQ(ierr);
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
	Vec alp0,FAlp;
	ierr = IGACreateMat(igaAl,&KAlp);CHKERRQ(ierr);
	ierr = IGACreateVec(igaAl,&alp0);CHKERRQ(ierr);
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
	ierr = KSPSetTolerances(kspAlp,1.0e-8,1.0e-30,PETSC_DEFAULT,PETSC_DEFAULT);CHKERRQ(ierr);
	ierr = KSPSolve(kspAlp,FAlp,alp0);CHKERRQ(ierr);

	ierr = KSPDestroy(&kspAlp);CHKERRQ(ierr);
	ierr = MatDestroy(&KAlp);CHKERRQ(ierr);
	ierr = VecDestroy(&FAlp);CHKERRQ(ierr);
	
	char nameAlp[]="/Alp-2d-0.dat";
	char pathAlp[512];
	sprintf(pathAlp,"%s%s",direct,nameAlp);
	ierr = IGAWriteVec(igaAl,alp0,pathAlp);CHKERRQ(ierr);
//

//Creation of types and systems for the Helmholtz decomposition of Up (or Ue), curl part
	//System for chiU
	PetscPrintf(PETSC_COMM_WORLD,"System for curl part of Helmholtz of Up starting \n");
	IGA igachiUp;
	ierr = IGACreate(PETSC_COMM_WORLD,&igachiUp);CHKERRQ(ierr);
	ierr = IGASetDim(igachiUp,2);CHKERRQ(ierr);														//Spatial dimension of the problem
	ierr = IGASetDof(igachiUp,4);CHKERRQ(ierr);														//Number of degrees of freedom, per node
	ierr = IGASetOrder(igachiUp,2);CHKERRQ(ierr);													//Number of spatial derivatives to calculate
	ierr = IGASetFromOptions(igachiUp);CHKERRQ(ierr);												//Note: The order (or degree) of the shape functions is given by the mesh!
	ierr = IGARead(igachiUp,"./geometry2.dat");CHKERRQ(ierr);
	ierr = IGASetUp(igachiUp);CHKERRQ(ierr);
	//PetscInt dir,side;
	for (dir=0; dir<2; dir++) 
	{
		for (side=0; side<2; side++) 
		{
			if(dir==0 && side==0)
			{
				ierr = IGASetBoundaryValue(igachiUp,dir,side,0,0.0);CHKERRQ(ierr);    				// Dirichlet boundary conditions
			}
			if(dir==0 && side==1)
			{
				ierr = IGASetBoundaryValue(igachiUp,dir,side,0,0.0);CHKERRQ(ierr);    				// Dirichlet boundary conditions
			}
			if(dir==1 && side==0)
			{
				ierr = IGASetBoundaryValue(igachiUp,dir,side,1,0.0);CHKERRQ(ierr);    				// Dirichlet boundary conditions
			}
			if(dir==1 && side==1)
			{
				ierr = IGASetBoundaryValue(igachiUp,dir,side,1,0.0);CHKERRQ(ierr);    				// Dirichlet boundary conditions
			}
			//ierr = IGASetBoundaryValue(iga,dir,side,dof,0.0);CHKERRQ(ierr);    				// Dirichlet boundary conditions
			//ierr = IGASetBoundaryForm(igachiUp,dir,side,PETSC_TRUE);CHKERRQ(ierr);  				// Neumann boundary conditions
		}
	}
	
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
	PetscReal       *Al0chiUp;				//U

  	IGAFormSystem  wtfchiUp;
 	void           *wtf2chiUp;

 	KSP kspchiUp;
	ierr = IGACreateKSP(igachiUp,&kspchiUp);CHKERRQ(ierr);

	//Get local vectors s0 and arrays
	ierr = IGAGetLocalVecArray(igaAl,alp0,&localAl0chiUp,&arrayAl0chiUp);CHKERRQ(ierr);
	//ierr = IGAGetLocalVecArray(igaZS,ZS0,&localZSchiUp,&arrayZSchiUp);CHKERRQ(ierr);

	//Element loop
	ierr = IGABeginElement(igachiUp,&elemchiUp);CHKERRQ(ierr);
	ierr = IGABeginElement(igaAl,&elemAlp);CHKERRQ(ierr);
	//ierr = IGABeginElement(igaZS,&elemZS);CHKERRQ(ierr);

	while (IGANextElement(igachiUp,elemchiUp))
	{
		IGANextElement(igaAl,elemAlp);
		//IGANextElement(igaZS,elemZS);
		ierr = IGAElementGetWorkMat(elemchiUp,&KlocchiUp);CHKERRQ(ierr);
		ierr = IGAElementGetWorkVec(elemchiUp,&FlocchiUp);CHKERRQ(ierr);
		ierr = IGAElementGetValues(elemAlp,arrayAl0chiUp,&Al0chiUp);CHKERRQ(ierr);
		//ierr = IGAElementGetValues(elemZS,arrayZSchiUp,&ZSchiUp);CHKERRQ(ierr);

		//FormSystem loop
		while (IGAElementNextFormSystem(elemchiUp,&wtfchiUp,&wtf2chiUp)) 
		{
			//Quadrature loop
			ierr = IGAElementBeginPoint(elemchiUp,&pointchiUp);CHKERRQ(ierr);
			ierr = IGAElementBeginPoint(elemAlp,&pointAlp);CHKERRQ(ierr);
			//ierr = IGAElementBeginPoint(elemZS,&pointZS);CHKERRQ(ierr);
			if(pointchiUp->atboundary==0 && pointAlp->atboundary==0)
			{
				while (IGAElementNextPoint(elemchiUp,pointchiUp)) 
				{
					IGAElementNextPoint(elemAlp,pointAlp);
					//IGAElementNextPoint(elemZS,pointZS);
					ierr = IGAPointGetWorkMat(pointchiUp,&KpointchiUp);CHKERRQ(ierr);
					ierr = IGAPointGetWorkVec(pointchiUp,&FpointchiUp);CHKERRQ(ierr);
					ierr = curlChiU(pointchiUp,pointAlp, KpointchiUp,FpointchiUp,Al0chiUp,NULL);CHKERRQ(ierr);
					ierr = IGAPointAddMat(pointchiUp,KpointchiUp,KlocchiUp);CHKERRQ(ierr);
					ierr = IGAPointAddVec(pointchiUp,FpointchiUp,FlocchiUp);CHKERRQ(ierr);
				}
				IGAElementNextPoint(elemAlp,pointAlp);
				//while (pointZS->index != -1)
				//{
				//	IGAElementNextPoint(elemZS,pointZS);
				//}
				ierr = IGAElementEndPoint(elemchiUp,&pointchiUp);CHKERRQ(ierr);
				ierr = IGAElementEndPoint(elemAlp,&pointAlp);CHKERRQ(ierr);
				//ierr = IGAElementEndPoint(elemZS,&pointZS);CHKERRQ(ierr);
			}
		}
		ierr = IGAElementFixSystem(elemchiUp,KlocchiUp,FlocchiUp);CHKERRQ(ierr);
		ierr = IGAElementAssembleMat(elemchiUp,KlocchiUp,KchiUp);CHKERRQ(ierr);
		ierr = IGAElementAssembleVec(elemchiUp,FlocchiUp,FchiUp);CHKERRQ(ierr);

	}
	IGANextElement(igaAl,elemAlp);
	//IGANextElement(igaZS,elemZS);
	ierr = IGAEndElement(igachiUp,&elemchiUp);CHKERRQ(ierr);
	ierr = IGAEndElement(igaAl,&elemAlp);CHKERRQ(ierr);
	//ierr = IGAEndElement(igaZS,&elemZS);CHKERRQ(ierr);

	// Restore local vectors s0 and arrays
	ierr = IGARestoreLocalVecArray(igaAl,alp0,&localAl0chiUp,&arrayAl0chiUp);CHKERRQ(ierr);
	//ierr = IGARestoreLocalVecArray(igaZS,ZS0,&localZSchiUp,&arrayZSchiUp);CHKERRQ(ierr);

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
	ierr = KSPSetFromOptions(kspchiUp);CHKERRQ(ierr);
	//ierr = KSPSetTolerances(kspchiUp,1.0e-8,1.0e-20,PETSC_DEFAULT,PETSC_DEFAULT);CHKERRQ(ierr);
	ierr = KSPSolve(kspchiUp,FchiUp,chiUp0);CHKERRQ(ierr);

	ierr = KSPDestroy(&kspchiUp);CHKERRQ(ierr);
	ierr = MatDestroy(&KchiUp);CHKERRQ(ierr);
	ierr = VecDestroy(&FchiUp);CHKERRQ(ierr);
	
	char namechiUp[]="/ChiUp-2d-0.dat";
	char pathchiUp[512];
	sprintf(pathchiUp,"%s%s",direct,namechiUp);
	ierr = IGAWriteVec(igachiUp,chiUp0,pathchiUp);CHKERRQ(ierr);
//

//System for initial state of z
	PetscPrintf(PETSC_COMM_WORLD,"System for Z0 starting \n");
	IGA igaZ0;
	ierr = IGACreate(PETSC_COMM_WORLD,&igaZ0);CHKERRQ(ierr);
	ierr = IGASetDim(igaZ0,2);CHKERRQ(ierr);													//Spatial dimension of the problem
	ierr = IGASetDof(igaZ0,2);CHKERRQ(ierr);													//Number of degrees of freedom, per node
	ierr = IGASetOrder(igaZ0,2);CHKERRQ(ierr);													//Number of spatial derivatives to calculate
	ierr = IGASetFromOptions(igaZ0);CHKERRQ(ierr);												//Note: The order (or degree) of the shape functions is given by the mesh!
	ierr = IGARead(igaZ0,"./geometry2.dat");CHKERRQ(ierr);
	ierr = IGASetUp(igaZ0);CHKERRQ(ierr);

	for (dir=0; dir<2; dir++) 
	{
		for (side=0; side<2; side++) 
		{
			//ierr = IGASetBoundaryValue(iga,dir,side,dof,0.0);CHKERRQ(ierr);    				// Dirichlet boundary conditions
			ierr = IGASetBoundaryForm(igaZ0,dir,side,PETSC_TRUE);CHKERRQ(ierr);  				// Neumann boundary conditions
		}
	}


	Mat KZ0;
	Vec Z0,FZ0;

	ierr = IGACreateMat(igaZ0,&KZ0);CHKERRQ(ierr);
	ierr = IGACreateVec(igaZ0,&Z0);CHKERRQ(ierr);
	ierr = IGACreateVec(igaZ0,&FZ0);CHKERRQ(ierr);

	IGAPoint		pointZ0;						//point
	IGAElement		elemZ0;							//element
	PetscReal		*KlocZ0,*FlocZ0;				//AA y BB
	PetscReal		*KpointZ0,*FpointZ0;			//KKK y FFF
	const PetscReal	*arrayChi0Z0,*arrayHsZ0;		//arrayU
	Vec				localChi0Z0,localHsZ0;			//localU
	PetscReal		*Chi0Z0,*HsZ0;					//U

  	IGAFormSystem	wtfZ0;
 	void			*wtf2Z0;

 	KSP kspZ0;
	ierr = IGACreateKSP(igaZ0,&kspZ0);CHKERRQ(ierr);

	// Get local vectors Chi0 and arrays
	ierr = IGAGetLocalVecArray(igachiUp,chiUp0,&localChi0Z0,&arrayChi0Z0);CHKERRQ(ierr);
	ierr = IGAGetLocalVecArray(igaZS,ZS0,&localHsZ0,&arrayHsZ0);CHKERRQ(ierr);

	// Element loop
	ierr = IGABeginElement(igaZ0,&elemZ0);CHKERRQ(ierr);
	ierr = IGABeginElement(igachiUp,&elemchiUp);CHKERRQ(ierr);
	ierr = IGABeginElement(igaZS,&elemZS);CHKERRQ(ierr);

	while (IGANextElement(igaZ0,elemZ0)) 
	{
		IGANextElement(igachiUp,elemchiUp);
		IGANextElement(igaZS,elemZS);

		ierr = IGAElementGetWorkMat(elemZ0,&KlocZ0);CHKERRQ(ierr);
		ierr = IGAElementGetWorkVec(elemZ0,&FlocZ0);CHKERRQ(ierr);
		ierr = IGAElementGetValues(elemchiUp,arrayChi0Z0,&Chi0Z0);CHKERRQ(ierr);
		ierr = IGAElementGetValues(elemZS,arrayHsZ0,&HsZ0);CHKERRQ(ierr);

		// FormSystem loop
		while (IGAElementNextFormSystem(elemZ0,&wtfZ0,&wtf2Z0)) 
		{
		// Quadrature loop
			ierr = IGAElementBeginPoint(elemZ0,&pointZ0);CHKERRQ(ierr);
			ierr = IGAElementBeginPoint(elemchiUp,&pointchiUp);CHKERRQ(ierr);
			ierr = IGAElementBeginPoint(elemZS,&pointZS);CHKERRQ(ierr);

			while (IGAElementNextPoint(elemZ0,pointZ0))
			{
				if(pointZ0->atboundary==1)
				{
					ierr = IGAPointGetWorkMat(pointZ0,&KpointZ0);CHKERRQ(ierr);
					ierr = IGAPointGetWorkVec(pointZ0,&FpointZ0);CHKERRQ(ierr);
					ierr = Z0sys(pointZ0,pointchiUp,pointZS,KpointZ0,FpointZ0,Chi0Z0,HsZ0,NULL);CHKERRQ(ierr);
					ierr = IGAPointAddMat(pointZ0,KpointZ0,KlocZ0);CHKERRQ(ierr);
					ierr = IGAPointAddVec(pointZ0,FpointZ0,FlocZ0);CHKERRQ(ierr);
				}

				if(pointZ0->atboundary==0 && pointchiUp->atboundary==0 && pointZS->atboundary==0)
				{
					IGAElementNextPoint(elemchiUp,pointchiUp);
					IGAElementNextPoint(elemZS,pointZS);

					ierr = IGAPointGetWorkMat(pointZ0,&KpointZ0);CHKERRQ(ierr);
					ierr = IGAPointGetWorkVec(pointZ0,&FpointZ0);CHKERRQ(ierr);
					ierr = Z0sys(pointZ0,pointchiUp,pointZS,KpointZ0,FpointZ0,Chi0Z0,HsZ0,NULL);CHKERRQ(ierr);
					ierr = IGAPointAddMat(pointZ0,KpointZ0,KlocZ0);CHKERRQ(ierr);
					ierr = IGAPointAddVec(pointZ0,FpointZ0,FlocZ0);CHKERRQ(ierr);
				}
			}
			if (pointchiUp->index != -1)
			{
				IGAElementNextPoint(elemchiUp,pointchiUp);
			}
			if (pointZS->index != -1)
			{
				IGAElementNextPoint(elemZS,pointZS);
			}
			ierr = IGAElementEndPoint(elemZ0,&pointZ0);CHKERRQ(ierr);
			ierr = IGAElementEndPoint(elemchiUp,&pointchiUp);CHKERRQ(ierr);
			ierr = IGAElementEndPoint(elemZS,&pointZS);CHKERRQ(ierr);
		}

		//ierr = IGAElementFixSystem(elemZ0,KlocZ0,FlocZ0);CHKERRQ(ierr);					//This sets Dirichlet condition ¿? (Yes, this applies the conditions from IGASetBoundaryValue)
		ierr = IGAElementAssembleMat(elemZ0,KlocZ0,KZ0);CHKERRQ(ierr);
		ierr = IGAElementAssembleVec(elemZ0,FlocZ0,FZ0);CHKERRQ(ierr);

	}
	IGANextElement(igachiUp,elemchiUp);
	IGANextElement(igaZS,elemZS);

	ierr = IGAEndElement(igaZ0,&elemZ0);CHKERRQ(ierr);
	ierr = IGAEndElement(igachiUp,&elemchiUp);CHKERRQ(ierr);
	ierr = IGAEndElement(igaZS,&elemZS);CHKERRQ(ierr);

	// Restore local vectors Chi0 and arrays
	ierr = IGARestoreLocalVecArray(igachiUp,chiUp0,&localChi0Z0,&arrayChi0Z0);CHKERRQ(ierr);
	ierr = IGARestoreLocalVecArray(igaZS,ZS0,&localHsZ0,&arrayHsZ0);CHKERRQ(ierr);

	//Here we set values to the Matrix directly (to impose Dirichlet condition in a single point)
	//Note: Lower Left corner is gdl's  0 and 1,
	//		Lower Right corner is gdl's 2*(nx+2)-2 and 2*(nx+2)-1 
	//		Upper Left corner is gdl's  2*(nx+2)*(ny+1)-2 and 2*(nx+2)*(ny+1)-1
	//		Upper Right corner is gdl's 2*(nx+2)*(ny+2)-2 and 2*(nx+2)*(ny+2)-1
	ierr = MatAssemblyBegin(KZ0,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyEnd  (KZ0,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);	
	
	rows=0;
	ierr = MatZeroRows(KZ0,1,&rows,1.0e60,0,0);CHKERRQ(ierr);
	rows=1;
	ierr = MatZeroRows(KZ0,1,&rows,1.0e60,0,0);CHKERRQ(ierr);
	rows=2*(nx+2)*(ny+2)-1;
	ierr = MatZeroRows(KZ0,1,&rows,1.0e60,0,0);CHKERRQ(ierr);

	ierr = MatAssemblyBegin(KZ0,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyEnd  (KZ0,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	
	
	//Here we set values to the vector directly (to impose Dirichlet condition in a single point)
	//Note: Lower Left corner is gdl's  0 and 1,
	//		Lower Right corner is gdl's 2*(nx+2)-2 and 2*(nx+2)-1 
	//		Upper Left corner is gdl's  2*(nx+2)*(ny+1)-2 and 2*(nx+2)*(ny+1)-1
	//		Upper Right corner is gdl's 2*(nx+2)*(ny+2)-2 and 2*(nx+2)*(ny+2)-1
	
	ierr = VecAssemblyBegin(FZ0);CHKERRQ(ierr);
	ierr = VecAssemblyEnd  (FZ0);CHKERRQ(ierr);

	rows=0;
	vals=0.0;
	ierr = VecSetValue(FZ0,rows,vals,INSERT_VALUES);CHKERRQ(ierr);
	rows=1;
	vals=0.0;
	ierr = VecSetValue(FZ0,rows,vals,INSERT_VALUES);CHKERRQ(ierr);
	rows=2*(nx+2)*(ny+2)-1;
	vals=0.0;
	ierr = VecSetValue(FZ0,rows,vals,INSERT_VALUES);CHKERRQ(ierr);

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

//System for L2 projection of stress
	PetscPrintf(PETSC_COMM_WORLD,"System for Stress starting \n");
	IGA igaStress;
	ierr = IGACreate(PETSC_COMM_WORLD,&igaStress);CHKERRQ(ierr);
	ierr = IGASetDim(igaStress,2);CHKERRQ(ierr);													//Spatial dimension of the problem
	ierr = IGASetDof(igaStress,4);CHKERRQ(ierr);													//Number of degrees of freedom, per node
	ierr = IGASetOrder(igaStress,2);CHKERRQ(ierr);													//Number of spatial derivatives to calculate
	ierr = IGASetFromOptions(igaStress);CHKERRQ(ierr);												//Note: The order (or degree) of the shape functions is given by the mesh!
	ierr = IGARead(igaStress,"./geometry2.dat");CHKERRQ(ierr);
	ierr = IGASetUp(igaStress);CHKERRQ(ierr);

	Mat KStress;
	Vec sigma0,FStress;

	ierr = IGACreateMat(igaStress,&KStress);CHKERRQ(ierr);
	ierr = IGACreateVec(igaStress,&sigma0);CHKERRQ(ierr);
	ierr = IGACreateVec(igaStress,&FStress);CHKERRQ(ierr);

	IGAPoint		pointStress;															//point
	IGAElement		elemStress;																//element
	PetscReal		*KlocStress,*FlocStress;												//AA y BB
	PetscReal		*KpointStress,*FpointStress;											//KKK y FFF
	const PetscReal	*arrayZ0Stress,*arrayChi0Stress;			//arrayU
	Vec				localZ0Stress,localChi0Stress;				//localU
	PetscReal		*Chi0Stress,*Z0Stress;											//U

  	IGAFormSystem	wtfStress;
 	void			*wtf2Stress;

 	KSP kspStress;
	ierr = IGACreateKSP(igaStress,&kspStress);CHKERRQ(ierr);

	// Get local vectors Z0 and Chi0 and arrays
	ierr = IGAGetLocalVecArray(igachiUp,chiUp0,&localChi0Stress,&arrayChi0Stress);CHKERRQ(ierr);
	ierr = IGAGetLocalVecArray(igaZ0,Z0,&localZ0Stress,&arrayZ0Stress);CHKERRQ(ierr);

	// Element loop
	ierr = IGABeginElement(igaStress,&elemStress);CHKERRQ(ierr);
	ierr = IGABeginElement(igaZ0,&elemZ0);CHKERRQ(ierr);
	ierr = IGABeginElement(igachiUp,&elemchiUp);CHKERRQ(ierr);

	while (IGANextElement(igaStress,elemStress)) 
	{
		IGANextElement(igaZ0,elemZ0);
		IGANextElement(igachiUp,elemchiUp);

		ierr = IGAElementGetWorkMat(elemStress,&KlocStress);CHKERRQ(ierr);
		ierr = IGAElementGetWorkVec(elemStress,&FlocStress);CHKERRQ(ierr);
		ierr = IGAElementGetValues(elemZ0,arrayZ0Stress,&Z0Stress);CHKERRQ(ierr);
		ierr = IGAElementGetValues(elemchiUp,arrayChi0Stress,&Chi0Stress);CHKERRQ(ierr);

		// FormSystem loop
		while (IGAElementNextFormSystem(elemStress,&wtfStress,&wtf2Stress)) 
		{
		// Quadrature loop
			ierr = IGAElementBeginPoint(elemStress,&pointStress);CHKERRQ(ierr);
			ierr = IGAElementBeginPoint(elemZ0,&pointZ0);CHKERRQ(ierr);
			ierr = IGAElementBeginPoint(elemchiUp,&pointchiUp);CHKERRQ(ierr);

			while (IGAElementNextPoint(elemStress,pointStress))
			{
				if(pointStress->atboundary==1)
				{
					IGAElementNextPoint(elemZ0,pointZ0);
					IGAElementNextPoint(elemchiUp,pointchiUp);
				}

				if(pointStress->atboundary==0 && pointZ0->atboundary==0 && pointchiUp->atboundary==0)
				{
					IGAElementNextPoint(elemZ0,pointZ0);
					IGAElementNextPoint(elemchiUp,pointchiUp);

					ierr = IGAPointGetWorkMat(pointStress,&KpointStress);CHKERRQ(ierr);
					ierr = IGAPointGetWorkVec(pointStress,&FpointStress);CHKERRQ(ierr);
					//PetscErrorCode Stress(IGAPoint p,IGAPoint pChi,IGAPoint pZu,PetscReal *K,PetscReal *F, PetscReal *Chi,PetscReal *Zu,void *ctx)
					ierr = Stress(pointStress,pointchiUp,pointZ0,KpointStress,FpointStress,Chi0Stress,Z0Stress,NULL);CHKERRQ(ierr);
					ierr = IGAPointAddMat(pointStress,KpointStress,KlocStress);CHKERRQ(ierr);
					ierr = IGAPointAddVec(pointStress,FpointStress,FlocStress);CHKERRQ(ierr);
				}
			}
			if (pointZ0->index != -1)
			{
				IGAElementNextPoint(elemZ0,pointZ0);
			}
			if (pointchiUp->index != -1)
			{
				IGAElementNextPoint(elemchiUp,pointchiUp);
			}
			ierr = IGAElementEndPoint(elemStress,&pointStress);CHKERRQ(ierr);
			ierr = IGAElementEndPoint(elemZ0,&pointZ0);CHKERRQ(ierr);
			ierr = IGAElementEndPoint(elemchiUp,&pointchiUp);CHKERRQ(ierr);
		}

		//ierr = IGAElementFixSystem(elemStress,KlocStress,FlocStress);CHKERRQ(ierr);					//This sets Dirichlet condition ¿? (Yes, this applies the conditions from IGASetBoundaryValue)
		ierr = IGAElementAssembleMat(elemStress,KlocStress,KStress);CHKERRQ(ierr);
		ierr = IGAElementAssembleVec(elemStress,FlocStress,FStress);CHKERRQ(ierr);
	}
	IGANextElement(igaZ0,elemZ0);
	IGANextElement(igachiUp,elemchiUp);

	ierr = IGAEndElement(igaStress,&elemStress);CHKERRQ(ierr);
	ierr = IGAEndElement(igaZ0,&elemZ0);CHKERRQ(ierr);
	ierr = IGAEndElement(igachiUp,&elemchiUp);CHKERRQ(ierr);

	// Restore local vectors u, Z0, Chi0 and arrays
	ierr = IGARestoreLocalVecArray(igaZ0,Z0,&localZ0Stress,&arrayZ0Stress);CHKERRQ(ierr);
	ierr = IGARestoreLocalVecArray(igachiUp,chiUp0,&localChi0Stress,&arrayChi0Stress);CHKERRQ(ierr);
	ierr = MatAssemblyBegin(KStress,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyEnd  (KStress,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

	ierr = VecAssemblyBegin(FStress);CHKERRQ(ierr);
	ierr = VecAssemblyEnd  (FStress);CHKERRQ(ierr);

	ierr = KSPSetOperators(kspStress,KStress,KStress);CHKERRQ(ierr);
	PC pcStress;
	ierr = KSPGetPC(kspStress,&pcStress); CHKERRQ(ierr);
	ierr = PCSetType(pcStress,PCLU); CHKERRQ(ierr);
	ierr = PCFactorSetMatSolverType(pcStress,MATSOLVERMUMPS); CHKERRQ(ierr);
	ierr = KSPSetFromOptions(kspStress);CHKERRQ(ierr);
	//ierr = KSPSetTolerances(kspStress,1e-18,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);CHKERRQ(ierr);
	ierr = KSPSolve(kspStress,FStress,sigma0);CHKERRQ(ierr);

	ierr = KSPDestroy(&kspStress);CHKERRQ(ierr);
	ierr = MatDestroy(&KStress);CHKERRQ(ierr);
	ierr = VecDestroy(&FStress);CHKERRQ(ierr);

	char nameStress[]="/sigma-2d-0.dat";
	char pathStress[512];
	sprintf(pathStress,"%s%s",direct,nameStress);
	ierr = IGAWriteVec(igaStress,sigma0,pathStress);CHKERRQ(ierr);	
//

//System for L2 projection of exact stress
	PetscPrintf(PETSC_COMM_WORLD,"System for L2 projection for exact stress starting \n");
	IGA igaExact;
	ierr = IGACreate(PETSC_COMM_WORLD,&igaExact);CHKERRQ(ierr);
	ierr = IGASetDim(igaExact,2);CHKERRQ(ierr);														//Spatial dimension of the problem
	ierr = IGASetDof(igaExact,4);CHKERRQ(ierr);														//Number of degrees of freedom, per node
	ierr = IGASetOrder(igaExact,1);CHKERRQ(ierr);													//Number of spatial derivatives to calculate
	ierr = IGASetFromOptions(igaExact);CHKERRQ(ierr);												//Note: The order (or degree) of the shape functions is given by the mesh!
	ierr = IGARead(igaExact,"./geometry.dat");CHKERRQ(ierr);
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

//Destroy all objects not needed anymore (Better to do it here in case different codes call the same IGA, move if memory is a problem)
	ierr = IGADestroy(&igaAl);CHKERRQ(ierr);
	ierr = IGADestroy(&igaS);CHKERRQ(ierr);
	ierr = IGADestroy(&igachiS);CHKERRQ(ierr);
	ierr = IGADestroy(&igaZS);CHKERRQ(ierr);
	ierr = IGADestroy(&igachiUp);CHKERRQ(ierr);
	ierr = IGADestroy(&igaExact);CHKERRQ(ierr);
//
*/

ierr = PetscFinalize();CHKERRQ(ierr);

return 0;
}