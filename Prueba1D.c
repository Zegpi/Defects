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

//System for z(0)
	#undef  __FUNCT__
	#define __FUNCT__ "Z0sys"
	PetscErrorCode Z0sys(IGAPoint p,PetscReal *K,PetscReal *F,void *ctx)
	{
		PetscReal C[3][3][3][3]={0};														//Initialization of elastic tensor
		const PetscReal *N0,(*N1)[2],(*N2)[2][2];
		IGAPointGetShapeFuns(p,0,(const PetscReal**)&N0);									//Value of the shape functions
		IGAPointGetShapeFuns(p,1,(const PetscReal**)&N1);									//Derivatives of the shape functions
		IGAPointGetShapeFuns(p,2,(const PetscReal**)&N2);									//Second derivatives of the shape funcions
		//After this command Na_xx=N2[a][0][0], Na_yy=N2[a][1][1], Na_xy=N2[a][0][1], Na_yx[a][1][0] (these last two are equal) (remember a is the index of the shape function)
		PetscInt a,b,i,j,k,l,u,w,r,s,m,nen=p->nen, dof=p->dof;

		PetscReal x[2];																		//Vector of reals, size equal to problem's dimension
		IGAPointFormGeomMap(p,x);															//Fills x with the coordinates of p, Gauss's point

		//Change for G=1
		const PetscReal nu=0.33;
		const PetscReal mu=1.0;
		const PetscReal lambda=2.0*mu*nu/(1.0-2.0*nu);
		const PetscReal eps=0.0*mu/1000.0;								//Choose later based on whatever Amit says :)

		PetscReal (*Keq)[dof][nen][dof] = (typeof(Keq)) K;
		PetscReal (*Feq)[dof] = (PetscReal (*)[dof])F;

		//////////////Delete this parte later, loads should come from appCtx or be 0
		PetscReal f[3]={0.0, 0.0, 0.0};			//Distributed load in body
		//PetscReal alpha=3.0;
		//PetscReal Lx=20.0; PetscReal Lx2=Lx*Lx; PetscReal Lx4=Lx*Lx*Lx*Lx;
		//PetscReal Ly=20.0; PetscReal Ly2=Ly*Ly; PetscReal Ly4=Ly*Ly*Ly*Ly;
		//f[0]=4.0*alpha/(425.0*Lx4*Ly4)*
		//	(1700.0*Lx4*Ly2 -20400.0*Lx4*x[1]*x[1] +204.0*Lx4 +6700.0*Lx2*Ly4 -13600.0*Lx2*Ly2*x[0]*x[0] -20000.0*Lx2*Ly2*x[0]*x[1] -53600.0*Lx2*Ly2*x[1]*x[1] +357.0*Lx2*Ly2 +163200.0*Lx2*x[0]*x[0]*x[1]*x[1]
		//	-1632.0*Lx2*x[0]*x[0] +107200.0*Lx2*x[1]*x[1]*x[1]*x[1] -4080.0*Lx2*x[1]*x[1] -80400.0*Ly4*x[0]*x[0] +408.0*Ly4 +27200.0*Ly2*x[0]*x[0]*x[0]*x[0] +643200.0*Ly2*x[0]*x[0]*x[1]*x[1] 
		//	-4080.0*Ly2*x[0]*x[0] -1632.0*Ly2*x[0]*x[1] -3264.0*Ly2*x[1]*x[1] -326400.0*x[0]*x[0]*x[0]*x[0]*x[1]*x[1] +3264.0*x[0]*x[0]*x[0]*x[0] -1286400.0*x[0]*x[0]*x[1]*x[1]*x[1]*x[1] 
		//	+48960.0*x[0]*x[0]*x[1]*x[1] +6582.0*x[0]*x[1]*x[1]*x[1] +6582.0*x[1]*x[1]*x[1]*x[1]);
		//f[1]=4.0*alpha/(425.0*Lx4*Ly4)*
		//	(3350.0*Lx4*Ly2 +850.0*Lx2*Ly4 -13400.0*Lx2*Ly2*x[0]*x[0] -80000.0*Lx2*Ly2*x[0]*x[1] -3400.0*Lx2*Ly2*x[1]*x[1] +153.0*Lx2*Ly2 +320000.0*Lx2*x[0]*x[1]*x[1]*x[1] -1632.0*Lx2*x[0]*x[1]
		//	-816.0*Lx2*x[1]*x[1] +320000.0*Lx2*x[0]*x[0]*x[0]*x[1] -816.0*Ly2*x[0]*x[0] -3264.0*Ly2*x[0]*x[1] -1280000.0*x[0]*x[0]*x[0]*x[1]*x[1]*x[1] +6528.0*x[0]*x[0]*x[0]*x[1] +9792.0*x[0]*x[0]*x[1]*x[1]
		//	+13056.0*x[0]*x[1]*x[1]*x[1]);
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

			PetscReal Omega=tan(45.0/180.0*ConstPi);
			PetscReal rho=sqrt(x[0]*x[0]+x[1]*x[1]);
			//const PetscReal burgers[2]={1.0,0.0};

			//Stress for single disclination
			Sborde[0][0]=mu*Omega/(2.0*ConstPi*(1.0-nu))*(log(rho)+(x[1]*x[1])/(rho*rho)+nu/(1.0-2.0*nu));
			Sborde[0][1]=-mu*Omega/(2.0*ConstPi*(1.0-nu))*x[0]*x[1]/(rho*rho);
			Sborde[1][0]=-mu*Omega/(2.0*ConstPi*(1.0-nu))*x[0]*x[1]/(rho*rho);
			Sborde[1][1]=mu*Omega/(2.0*ConstPi*(1.0-nu))*(log(rho)+(x[0]*x[0])/(rho*rho)+nu/(1.0-2.0*nu));
			Sborde[0][2]=0.0; Sborde[1][2]=0.0; Sborde[2][0]=0.0; Sborde[2][1]=0.0; Sborde[2][2]=0.0;

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
			//Sborde[0][0]=0.0;
			//Sborde[0][1]=0.0;
			//Sborde[1][0]=0.0;
			//Sborde[1][1]=0.0;
			//Sborde[0][2]=0.0; Sborde[1][2]=0.0; Sborde[2][0]=0.0; Sborde[2][1]=0.0; Sborde[2][2]=0.0;

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
											Keq[a][i][b][j]+=C[k][l][u][w]*0.5*(dz[u][w]+dz[w][u])*0.5*(dv[k][l]+dv[l][k]);
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
										Keq[a][i][b][j]+=0.5*eps*(d2z[r][s][m]+d2z[s][r][m])*0.5*(d2v[r][s][m]+d2v[s][r][m])
														-0.5*eps*(d2z[r][s][m]-d2z[s][r][m])*d2v[s][r][m];
									}
								}
							}

						}
					}
				}
			}

			for(a=0 ;a<nen; a++)
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

					for (r=0; r<3; r++)
					{
						for (s=0; s<3; s++)
						{
							for (m=0; m<3; m++)
							{
								Feq[a][i]+=M[r]*e[r][s][m]*dv[m][s];
							}
						}	
					}
				}
			}
		}
		return 0;
	}
//

//System for L2 projection of stresses
	#undef  __FUNCT__
	#define __FUNCT__ "Stress"
	//PetscErrorCode Stress(IGAPoint p,IGAPoint pU, IGAPoint pHs,IGAPoint pChi,IGAPoint pZu,PetscReal *K,PetscReal *F,PetscReal *U,PetscReal *HS, PetscReal *Chi,PetscReal *Zu,void *ctx)		//When other fields like Pi or S (etc.) come into this, declare a IGAPoint pPi or pS and a new PescReal *UPi or *US for each
	PetscErrorCode Stress(IGAPoint p,IGAPoint pZu,PetscReal *K,PetscReal *F,PetscReal *Zu,void *ctx)		//When other fields like Pi or S (etc.) come into this, declare a IGAPoint pPi or pS and a new PestcReal *UPi or *US for each
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
		const PetscReal eps=0.0*mu/1000.0;															//Choose later based on whatever Amit says :)

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
		//PetscReal dv[3][3][3]={0};

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
							Fstress[a][i]+=C[m][n][k][l]*fulld_z[k][l]*v[m][n];
						}

						Fstress[a][i]+= -0.25*eps*(fulld3_z[m][n][k][k]+fulld3_z[m][k][n][k]
							                      +fulld3_z[n][m][k][k]+fulld3_z[n][k][m][k])*v[m][n];

					}
				}
			}
		}
		return 0;
	}
//

//System for L2 projection of classic stresses (C*Ue)
	#undef  __FUNCT__
	#define __FUNCT__ "ClassicStress"
	//PetscErrorCode Stress(IGAPoint p,IGAPoint pU, IGAPoint pHs,IGAPoint pChi,IGAPoint pZu,PetscReal *K,PetscReal *F,PetscReal *U,PetscReal *HS, PetscReal *Chi,PetscReal *Zu,void *ctx)		//When other fields like Pi or S (etc.) come into this, declare a IGAPoint pPi or pS and a new PescReal *UPi or *US for each
	PetscErrorCode ClassicStress(IGAPoint p,IGAPoint pZu,PetscReal *K,PetscReal *F,PetscReal *Zu,void *ctx)		//When other fields like Pi or S (etc.) come into this, declare a IGAPoint pPi or pS and a new PestcReal *UPi or *US for each
	{
		const PetscReal *N0;
		IGAPointGetShapeFuns(p,0,(const PetscReal**)&N0);									//Value of the shape functions
		PetscInt a,b,i,j,k,l,m,n,nen=p->nen, dof=p->dof;

		//Change to consider G=1
		const PetscReal nu=0.33;
		const PetscReal mu=1.0;
		const PetscReal lambda=2.0*mu*nu/(1.0-2.0*nu);

		PetscReal d_Z0[2][2];																//Same for its gradient
		IGAPointFormGrad (pZu,Zu,&d_Z0[0][0]);												//Same for the gradient

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
							Fstress[a][i]+=0.5*(C[m][n][k][l]*(fulld_z[k][l])+C[n][m][k][l]*(fulld_z[k][l]))*v[m][n];
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
	PetscErrorCode CoupleStress(IGAPoint p,IGAPoint pZu,PetscReal *K,PetscReal *F,PetscReal *Zu,void *ctx)		//When other fields like Pi or S (etc.) come into this, declare a IGAPoint pPi or pS and a new PestcReal *UPi or *US for each
	{
		const PetscReal *N0;
		IGAPointGetShapeFuns(p,0,(const PetscReal**)&N0);									//Value of the shape functions
		PetscInt a,b,i,k,l,m,n,nen=p->nen, dof=p->dof;

		//Change to consider G=1
		//const PetscReal nu=0.33;
		const PetscReal mu=1.0;
		//const PetscReal lambda=2.0*mu*nu/(1.0-2.0*nu);
		const PetscReal eps=0.0*mu/1000.0;															//Choose later based on whatever Amit says :)

		//Definition of alternating tensor
		const PetscReal e[3][3][3]=
		{
			{{0.0,0.0,0.0},{0.0,0.0,1.0},{0.0,-1.0,0.0}},
			{{0.0,0.0,-1.0},{0.0,0.0,0.0},{1.0,0.0,0.0}},
			{{0.0,1.0,0.0},{-1.0,0.0,0.0},{0.0,0.0,0.0}}
		};
		
		PetscReal d2_Z0[2][2][2];															//Same for its 2nd order partial derivatives
		IGAPointFormHess (pZu,Zu,&d2_Z0[0][0][0]);											//Same for the hessian (second derivatives)

		//Expanding z (and derivatives) to 3 components, more convenient for sums in for loops
		PetscReal fulld2_z[3][3][3]={0};
		fulld2_z[0][0][0]=d2_Z0[0][0][0]; fulld2_z[0][0][1]=d2_Z0[0][0][1];
		fulld2_z[0][1][0]=d2_Z0[0][1][0]; fulld2_z[0][1][1]=d2_Z0[0][1][1];
		fulld2_z[1][0][0]=d2_Z0[1][0][0]; fulld2_z[1][0][1]=d2_Z0[1][0][1];
		fulld2_z[1][1][0]=d2_Z0[1][1][0]; fulld2_z[1][1][1]=d2_Z0[1][1][1];

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
							FCS[a][i]+=-0.5*eps*e[m][k][l]*(fulld2_z[k][l][n]+fulld2_z[k][n][l])*v[m][n];
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
	PetscErrorCode EnergyDensity(IGAPoint p,IGAPoint pZu,PetscReal *K,PetscReal *F,PetscReal *Zu,void *ctx)		//When other fields like Pi or S (etc.) come into this, declare a IGAPoint pPi or pS and a new PestcReal *UPi or *US for each
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
		const PetscReal eps=0.0*mu/1000.0;															//Choose later based on whatever Amit says :)

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

		PetscReal d_Z0[2][2];																//Same for its gradient
		PetscReal d2_Z0[2][2][2];															//Same for its 2nd order partial derivatives
		IGAPointFormGrad (pZu,Zu,&d_Z0[0][0]);												//Same for the gradient
		IGAPointFormHess (pZu,Zu,&d2_Z0[0][0][0]);											//Same for the hessian (second derivatives)

		//Expanding z (and derivatives) to 3 components, more convenient for sums in for loops
		PetscReal fulld_z[3][3]={0};
		fulld_z[0][0]=d_Z0[0][0]; fulld_z[0][1]=d_Z0[0][1];
		fulld_z[1][0]=d_Z0[1][0]; fulld_z[1][1]=d_Z0[1][1];

		PetscReal fulld2_z[3][3][3]={0};
		fulld2_z[0][0][0]=d2_Z0[0][0][0]; fulld2_z[0][0][1]=d2_Z0[0][0][1];
		fulld2_z[0][1][0]=d2_Z0[0][1][0]; fulld2_z[0][1][1]=d2_Z0[0][1][1];
		fulld2_z[1][0][0]=d2_Z0[1][0][0]; fulld2_z[1][0][1]=d2_Z0[1][0][1];
		fulld2_z[1][1][0]=d2_Z0[1][1][0]; fulld2_z[1][1][1]=d2_Z0[1][1][1];

		PetscReal (*KED)[dof][nen][dof] = (typeof(KED)) K;
		PetscReal (*FED)[dof] = (PetscReal (*)[dof])F;

		//PetscReal v[3][3]={0};

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
									FED[a][i]+=0.5*((fulld_z[k][l])*C[k][l][m][n]*(fulld_z[m][n]))*N0[a];
								}
								FED[a][i]+=0.5*eps*(fulld2_z[k][l][m])*(fulld2_z[k][l][m])*N0[a];
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

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc, char *argv[]) {

//Creation of solution systems
	PetscErrorCode  ierr;
	ierr = PetscInitialize(&argc,&argv,0,0);CHKERRQ(ierr);										//Always initialize PETSc
	PetscInt dir,side;
	PetscPrintf(PETSC_COMM_WORLD,"Start of PruebaV5 \n");

	PetscInt commsize,rank;
	ierr = MPI_Comm_size(PETSC_COMM_WORLD,&commsize);CHKERRQ(ierr);
	ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank);CHKERRQ(ierr);
//

//App context creation and some data
	//Mesh parameters (to fix specific points in z0 system)
	PetscInt b=301;				//Parmeter to choose size of cores, must always be odd, core will be of size 1 unit, rest of the body will be of size b-1 units in each direction
	PetscReal Lx=20.0;
	PetscReal Ly=20.0;
	PetscInt  nx=b;
	PetscInt  ny=b;

	//Directory to write files to
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

	time_t T=time(NULL);
	struct tm tm=*localtime(&T);
	PetscPrintf(PETSC_COMM_WORLD,"Current time is %02d:%02d:%02d \n",tm.tm_hour,tm.tm_min,tm.tm_sec);

	if(rank==0)
	{
		ierr=system("cp ~/shared/CodigosPetIGA/PruebaV5.c ~/sharedResults/");
	}
//

//System for initial state of z0
	PetscPrintf(PETSC_COMM_WORLD,"System for Z0 starting \n");
	IGA igaZ0;
	ierr = IGACreate(PETSC_COMM_WORLD,&igaZ0);CHKERRQ(ierr);
	ierr = IGASetDim(igaZ0,1);CHKERRQ(ierr);													//Spatial dimension of the problem
	ierr = IGASetDof(igaZ0,1);CHKERRQ(ierr);													//Number of degrees of freedom, per node
	ierr = IGASetOrder(igaZ0,3);CHKERRQ(ierr);													//Number of spatial derivatives to calculate
	ierr = IGASetFromOptions(igaZ0);CHKERRQ(ierr);												//Note: The order (or degree) of the shape functions is given by the mesh!
	ierr = IGARead(igaZ0,"./geometry31D.dat");CHKERRQ(ierr);
	ierr = IGASetUp(igaZ0);CHKERRQ(ierr);

	PetscInt fijaPunto=1;																		//Fix a single point (1) or a side (chosen in blocks below)

	for (dir=0; dir<1; dir++) 
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
		ierr = IGASetBoundaryValue(igaZ0,0,0,0,0.0);CHKERRQ(ierr);	//Left side, 1st dof = 0
		ierr = IGASetBoundaryValue(igaZ0,0,1,0,0.0);CHKERRQ(ierr);	//Right side, 1st dof=0
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

  	IGAFormSystem	wtfZ0;
 	void			*wtf2Z0;

 	KSP kspZ0;
	ierr = IGACreateKSP(igaZ0,&kspZ0);CHKERRQ(ierr);

	// Element loop
	ierr = IGABeginElement(igaZ0,&elemZ0);CHKERRQ(ierr);

	while (IGANextElement(igaZ0,elemZ0)) 
	{
		ierr = IGAElementGetWorkMat(elemZ0,&KlocZ0);CHKERRQ(ierr);
		ierr = IGAElementGetWorkVec(elemZ0,&FlocZ0);CHKERRQ(ierr);

		// FormSystem loop
		while (IGAElementNextFormSystem(elemZ0,&wtfZ0,&wtf2Z0)) 
		{
		// Quadrature loop
			ierr = IGAElementBeginPoint(elemZ0,&pointZ0);CHKERRQ(ierr);

			while (IGAElementNextPoint(elemZ0,pointZ0))
			{
				if(pointZ0->atboundary==1)
				{
					ierr = IGAPointGetWorkMat(pointZ0,&KpointZ0);CHKERRQ(ierr);
					ierr = IGAPointGetWorkVec(pointZ0,&FpointZ0);CHKERRQ(ierr);
					ierr = Z0sys(pointZ0,KpointZ0,FpointZ0,NULL);CHKERRQ(ierr);
					ierr = IGAPointAddMat(pointZ0,KpointZ0,KlocZ0);CHKERRQ(ierr);
					ierr = IGAPointAddVec(pointZ0,FpointZ0,FlocZ0);CHKERRQ(ierr);
				}

				if(pointZ0->atboundary==0)
				{
					ierr = IGAPointGetWorkMat(pointZ0,&KpointZ0);CHKERRQ(ierr);
					ierr = IGAPointGetWorkVec(pointZ0,&FpointZ0);CHKERRQ(ierr);
					ierr = Z0sys(pointZ0,KpointZ0,FpointZ0,NULL);CHKERRQ(ierr);
					ierr = IGAPointAddMat(pointZ0,KpointZ0,KlocZ0);CHKERRQ(ierr);
					ierr = IGAPointAddVec(pointZ0,FpointZ0,FlocZ0);CHKERRQ(ierr);
				}
			}
			ierr = IGAElementEndPoint(elemZ0,&pointZ0);CHKERRQ(ierr);
		}

		if(fijaPunto==0)
		{
			ierr = IGAElementFixSystem(elemZ0,KlocZ0,FlocZ0);CHKERRQ(ierr);					//This sets Dirichlet condition ¿? (Yes, this applies the conditions from IGASetBoundaryValue)
		}
		ierr = IGAElementAssembleMat(elemZ0,KlocZ0,KZ0);CHKERRQ(ierr);
		ierr = IGAElementAssembleVec(elemZ0,FlocZ0,FZ0);CHKERRQ(ierr);

	}

	ierr = IGAEndElement(igaZ0,&elemZ0);CHKERRQ(ierr);

	ierr = MatAssemblyBegin(KZ0,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyEnd  (KZ0,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);	

	if (fijaPunto==0)
	{
		/*
		PetscInt rows=0;
		ierr = MatZeroRows(KZ0,1,&rows,1.0e90,0,0);CHKERRQ(ierr);
		ierr = MatAssemblyBegin(KZ0,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
		ierr = MatAssemblyEnd  (KZ0,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

		ierr = VecAssemblyBegin(FZ0);CHKERRQ(ierr);
		ierr = VecAssemblyEnd  (FZ0);CHKERRQ(ierr);
		
		rows=0;
		PetscReal vals=0.0;
		ierr = VecSetValue(FZ0,rows,vals,INSERT_VALUES);CHKERRQ(ierr);
		*/
	}

	if (fijaPunto==1)
	{
		//Here we set values to the Matrix directly, to impose Dirichlet condition in a single point.
		//Note: Lower Left corner is gdl's  0 and 1,
		//		Lower Right corner is gdl's 2*(nx+2)-2 and 2*(nx+2)-1 
		//		Upper Left corner is gdl's  2*(nx+2)*(ny+2)-2*(nx+1)-2 and 2*(nx+2)*(ny+2)-2*(nx+1)-1
		//		Upper Right corner is gdl's 2*(nx+2)*(ny+2)-2 and 2*(nx+2)*(ny+2)-1
		//All of these for when z is a 2nd order nurbs

		PetscInt rows=0;
		ierr = MatZeroRows(KZ0,1,&rows,1.0e90,0,0);CHKERRQ(ierr);
		rows=1;
		ierr = MatZeroRows(KZ0,1,&rows,1.0e90,0,0);CHKERRQ(ierr);
		rows=2*(nx+3)*(ny+3)-2*(nx+2)-2; 										//This is for when z is a 3rd order nurb
		//rows=2*(nx+2)*(ny+2)-2*(nx+1)-2; 										//This is for when z is a 2nd order nurb
		//rows=2*(nx+1)-1;														//This is the dof in x on the upper right corner for 1st order elements
		ierr = MatZeroRows(KZ0,1,&rows,1.0e90,0,0);CHKERRQ(ierr);

		ierr = MatAssemblyBegin(KZ0,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
		ierr = MatAssemblyEnd  (KZ0,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

		PetscPrintf(PETSC_COMM_WORLD, "\n \n info sobre la matriz \n");
		//PetscViewer viewer;
		PetscViewerPushFormat(PETSC_VIEWER_STDOUT_WORLD,PETSC_VIEWER_ASCII_INFO_DETAIL);
		MatView(KZ0,PETSC_VIEWER_STDOUT_WORLD);
		//PetscViewerDestroy(&viewer);
		
		PetscPrintf(PETSC_COMM_WORLD, "\n \n info sobre la matriz termina aqui \n");

		//Here we set values to the vector directly (to impose Dirichlet condition in a single point)
		//Note: Lower Left corner is gdl's  0 and 1,
		//		Lower Right corner is gdl's 2*(nx+2)-2 and 2*(nx+2)-1 
		//		Upper Left corner is gdl's  2*(nx+2)*(ny+2)-2-2*(nx+1) and 2*(nx+2)*(ny+2)-1-2*(nx+1)
		//		Upper Right corner is gdl's 2*(nx+2)*(ny+2)-2 and 2*(nx+2)*(ny+2)-1
		//All of these for when z is a 2nd order nurbs
		ierr = VecAssemblyBegin(FZ0);CHKERRQ(ierr);
		ierr = VecAssemblyEnd  (FZ0);CHKERRQ(ierr);

		rows=0;
		PetscReal vals=0.0;
		ierr = VecSetValue(FZ0,rows,vals,INSERT_VALUES);CHKERRQ(ierr);
		rows=1;
		vals=0.0;
		ierr = VecSetValue(FZ0,rows,vals,INSERT_VALUES);CHKERRQ(ierr);
		rows=2*(nx+3)*(ny+3)-2*(nx+2)-2;
		//rows=2*(nx+2)*(ny+2)-2*(nx+1)-2; 										//This is for when z is a 2nd order nurb
		//rows=2*(nx+1)-1;
		vals=0.0;
		ierr = VecSetValue(FZ0,rows,vals,INSERT_VALUES);CHKERRQ(ierr);
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

//System for L2 projection of stress
	PetscPrintf(PETSC_COMM_WORLD,"System for Stress starting \n");
	IGA igaStress;
	ierr = IGACreate(PETSC_COMM_WORLD,&igaStress);CHKERRQ(ierr);
	ierr = IGASetDim(igaStress,1);CHKERRQ(ierr);													//Spatial dimension of the problem
	ierr = IGASetDof(igaStress,1);CHKERRQ(ierr);													//Number of degrees of freedom, per node
	ierr = IGASetOrder(igaStress,2);CHKERRQ(ierr);													//Number of spatial derivatives to calculate
	ierr = IGASetFromOptions(igaStress);CHKERRQ(ierr);												//Note: The order (or degree) of the shape functions is given by the mesh!
	ierr = IGARead(igaStress,"./geometry31D.dat");CHKERRQ(ierr);
	ierr = IGASetUp(igaStress);CHKERRQ(ierr);

	for (dir=0; dir<2; dir++) 
	{
		for (side=0; side<1; side++) 
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
	const PetscReal	*arrayZ0Stress;															//arrayU
	Vec				localZ0Stress;															//localU
	PetscReal		*Z0Stress;																//U

  	IGAFormSystem	wtfStress;
 	void			*wtf2Stress;

 	KSP kspStress;
	ierr = IGACreateKSP(igaStress,&kspStress);CHKERRQ(ierr);

	// Get local vectors Z0 arrays
	ierr = IGAGetLocalVecArray(igaZ0,Z0,&localZ0Stress,&arrayZ0Stress);CHKERRQ(ierr);

	// Element loop
	ierr = IGABeginElement(igaStress,&elemStress);CHKERRQ(ierr);
	ierr = IGABeginElement(igaZ0,&elemZ0);CHKERRQ(ierr);

	while (IGANextElement(igaStress,elemStress)) 
	{
		IGANextElement(igaZ0,elemZ0);

		ierr = IGAElementGetWorkMat(elemStress,&KlocStress);CHKERRQ(ierr);
		ierr = IGAElementGetWorkVec(elemStress,&FlocStress);CHKERRQ(ierr);
		ierr = IGAElementGetValues(elemZ0,arrayZ0Stress,&Z0Stress);CHKERRQ(ierr);

		// FormSystem loop
		while (IGAElementNextFormSystem(elemStress,&wtfStress,&wtf2Stress)) 
		{
		// Quadrature loop
			ierr = IGAElementBeginPoint(elemStress,&pointStress);CHKERRQ(ierr);
			ierr = IGAElementBeginPoint(elemZ0,&pointZ0);CHKERRQ(ierr);

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

				if(pointStress->atboundary==0 && pointZ0->atboundary==0)
				{
					IGAElementNextPoint(elemZ0,pointZ0);

					ierr = IGAPointGetWorkMat(pointStress,&KpointStress);CHKERRQ(ierr);
					ierr = IGAPointGetWorkVec(pointStress,&FpointStress);CHKERRQ(ierr);
					//PetscErrorCode Stress(IGAPoint p,IGAPoint pChi,IGAPoint pZu,PetscReal *K,PetscReal *F, PetscReal *Chi,PetscReal *Zu,void *ctx)
					ierr = Stress(pointStress,pointZ0,KpointStress,FpointStress,Z0Stress,NULL);CHKERRQ(ierr);
					ierr = IGAPointAddMat(pointStress,KpointStress,KlocStress);CHKERRQ(ierr);
					ierr = IGAPointAddVec(pointStress,FpointStress,FlocStress);CHKERRQ(ierr);
				}
			}
			if (pointZ0->index != -1)
			{
				IGAElementNextPoint(elemZ0,pointZ0);
			}

			ierr = IGAElementEndPoint(elemStress,&pointStress);CHKERRQ(ierr);
			ierr = IGAElementEndPoint(elemZ0,&pointZ0);CHKERRQ(ierr);
		}

		//ierr = IGAElementFixSystem(elemStress,KlocStress,FlocStress);CHKERRQ(ierr);					//This sets Dirichlet condition ¿? (Yes, this applies the conditions from IGASetBoundaryValue)
		ierr = IGAElementAssembleMat(elemStress,KlocStress,KStress);CHKERRQ(ierr);
		ierr = IGAElementAssembleVec(elemStress,FlocStress,FStress);CHKERRQ(ierr);
	}
	IGANextElement(igaZ0,elemZ0);

	ierr = IGAEndElement(igaStress,&elemStress);CHKERRQ(ierr);
	ierr = IGAEndElement(igaZ0,&elemZ0);CHKERRQ(ierr);

	// Restore local vectors u, Z0, Chi0 and arrays
	ierr = IGARestoreLocalVecArray(igaZ0,Z0,&localZ0Stress,&arrayZ0Stress);CHKERRQ(ierr);
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
	ierr = KSPSetTolerances(kspStress,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);CHKERRQ(ierr);
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
	PetscPrintf(PETSC_COMM_WORLD,"System for Classic Stress starting \n");
	IGA igaClassicStress;
	ierr = IGACreate(PETSC_COMM_WORLD,&igaClassicStress);CHKERRQ(ierr);
	ierr = IGASetDim(igaClassicStress,1);CHKERRQ(ierr);													//Spatial dimension of the problem
	ierr = IGASetDof(igaClassicStress,1);CHKERRQ(ierr);													//Number of degrees of freedom, per node
	ierr = IGASetOrder(igaClassicStress,2);CHKERRQ(ierr);													//Number of spatial derivatives to calculate
	ierr = IGASetFromOptions(igaClassicStress);CHKERRQ(ierr);												//Note: The order (or degree) of the shape functions is given by the mesh!
	ierr = IGARead(igaClassicStress,"./geometry31D.dat");CHKERRQ(ierr);
	ierr = IGASetUp(igaClassicStress);CHKERRQ(ierr);

	for (dir=0; dir<1; dir++) 
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
	IGAElement		elemClassicStress;															//element
	PetscReal		*KlocClassicStress,*FlocClassicStress;										//AA y BB
	PetscReal		*KpointClassicStress,*FpointClassicStress;									//KKK y FFF
	const PetscReal	*arrayZ0ClassicStress;														//arrayU
	Vec				localZ0ClassicStress;														//localU
	PetscReal		*Z0ClassicStress;															//U

  	IGAFormSystem	wtfClassicStress;
 	void			*wtf2ClassicStress;

 	KSP kspClassicStress;
	ierr = IGACreateKSP(igaClassicStress,&kspClassicStress);CHKERRQ(ierr);

	// Get local vectors Z0 and Chi0 and arrays
	ierr = IGAGetLocalVecArray(igaZ0,Z0,&localZ0ClassicStress,&arrayZ0ClassicStress);CHKERRQ(ierr);

	// Element loop
	ierr = IGABeginElement(igaClassicStress,&elemClassicStress);CHKERRQ(ierr);
	ierr = IGABeginElement(igaZ0,&elemZ0);CHKERRQ(ierr);

	while (IGANextElement(igaClassicStress,elemClassicStress)) 
	{
		IGANextElement(igaZ0,elemZ0);

		ierr = IGAElementGetWorkMat(elemClassicStress,&KlocClassicStress);CHKERRQ(ierr);
		ierr = IGAElementGetWorkVec(elemClassicStress,&FlocClassicStress);CHKERRQ(ierr);
		ierr = IGAElementGetValues(elemZ0,arrayZ0ClassicStress,&Z0ClassicStress);CHKERRQ(ierr);

		// FormSystem loop
		while (IGAElementNextFormSystem(elemClassicStress,&wtfClassicStress,&wtf2ClassicStress)) 
		{
		// Quadrature loop
			ierr = IGAElementBeginPoint(elemClassicStress,&pointClassicStress);CHKERRQ(ierr);
			ierr = IGAElementBeginPoint(elemZ0,&pointZ0);CHKERRQ(ierr);

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

				if(pointClassicStress->atboundary==0 && pointZ0->atboundary==0)
				{
					IGAElementNextPoint(elemZ0,pointZ0);

					ierr = IGAPointGetWorkMat(pointClassicStress,&KpointClassicStress);CHKERRQ(ierr);
					ierr = IGAPointGetWorkVec(pointClassicStress,&FpointClassicStress);CHKERRQ(ierr);
					//PetscErrorCode Stress(IGAPoint p,IGAPoint pChi,IGAPoint pZu,PetscReal *K,PetscReal *F, PetscReal *Chi,PetscReal *Zu,void *ctx)
					ierr = ClassicStress(pointClassicStress,pointZ0,KpointClassicStress,FpointClassicStress,Z0ClassicStress,NULL);CHKERRQ(ierr);
					ierr = IGAPointAddMat(pointClassicStress,KpointClassicStress,KlocClassicStress);CHKERRQ(ierr);
					ierr = IGAPointAddVec(pointClassicStress,FpointClassicStress,FlocClassicStress);CHKERRQ(ierr);
				}
			}
			if (pointZ0->index != -1)
			{
				IGAElementNextPoint(elemZ0,pointZ0);
			}
			ierr = IGAElementEndPoint(elemClassicStress,&pointClassicStress);CHKERRQ(ierr);
			ierr = IGAElementEndPoint(elemZ0,&pointZ0);CHKERRQ(ierr);
		}

		//ierr = IGAElementFixSystem(elemStress,KlocStress,FlocStress);CHKERRQ(ierr);					//This sets Dirichlet condition ¿? (Yes, this applies the conditions from IGASetBoundaryValue)
		ierr = IGAElementAssembleMat(elemClassicStress,KlocClassicStress,KClassicStress);CHKERRQ(ierr);
		ierr = IGAElementAssembleVec(elemClassicStress,FlocClassicStress,FClassicStress);CHKERRQ(ierr);
	}
	IGANextElement(igaZ0,elemZ0);

	ierr = IGAEndElement(igaClassicStress,&elemClassicStress);CHKERRQ(ierr);
	ierr = IGAEndElement(igaZ0,&elemZ0);CHKERRQ(ierr);

	// Restore local vectors u, Z0, Chi0 and arrays
	ierr = IGARestoreLocalVecArray(igaZ0,Z0,&localZ0ClassicStress,&arrayZ0ClassicStress);CHKERRQ(ierr);
	ierr = MatAssemblyBegin(KClassicStress,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyEnd  (KClassicStress,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

	ierr = VecAssemblyBegin(FClassicStress);CHKERRQ(ierr);
	ierr = VecAssemblyEnd  (FClassicStress);CHKERRQ(ierr);

	ierr = KSPSetOperators(kspClassicStress,KClassicStress,KClassicStress);CHKERRQ(ierr);
	//PC pcStress;
	//ierr = KSPGetPC(kspStress,&pcStress); CHKERRQ(ierr);
	//ierr = PCSetType(pcStress,PCLU); CHKERRQ(ierr);
	//ierr = PCFactorSetMatSolverType(pcStress,MATSOLVERMUMPS); CHKERRQ(ierr);
	ierr = KSPSetFromOptions(kspClassicStress);CHKERRQ(ierr);
	ierr = KSPSetTolerances(kspClassicStress,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);CHKERRQ(ierr);
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
	PetscPrintf(PETSC_COMM_WORLD,"System for CoupleStress starting \n");
	IGA igaCS;
	ierr = IGACreate(PETSC_COMM_WORLD,&igaCS);CHKERRQ(ierr);
	ierr = IGASetDim(igaCS,1);CHKERRQ(ierr);													//Spatial dimension of the problem
	ierr = IGASetDof(igaCS,1);CHKERRQ(ierr);													//Number of degrees of freedom, per node
	ierr = IGASetOrder(igaCS,2);CHKERRQ(ierr);													//Number of spatial derivatives to calculate
	ierr = IGASetFromOptions(igaCS);CHKERRQ(ierr);												//Note: The order (or degree) of the shape functions is given by the mesh!
	ierr = IGARead(igaCS,"./geometry31D.dat");CHKERRQ(ierr);
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

	IGAPoint		pointCS;															//point
	IGAElement		elemCS;																//element
	PetscReal		*KlocCS,*FlocCS;												//AA y BB
	PetscReal		*KpointCS,*FpointCS;											//KKK y FFF
	const PetscReal	*arrayZ0CS;			//arrayU
	Vec				localZ0CS;				//localU
	PetscReal		*Z0CS;											//U

  	IGAFormSystem	wtfCS;
 	void			*wtf2CS;

 	KSP kspCS;
	ierr = IGACreateKSP(igaCS,&kspCS);CHKERRQ(ierr);

	// Get local vectors Z0 and Chi0 and arrays
	ierr = IGAGetLocalVecArray(igaZ0,Z0,&localZ0CS,&arrayZ0CS);CHKERRQ(ierr);

	// Element loop
	ierr = IGABeginElement(igaCS,&elemCS);CHKERRQ(ierr);
	ierr = IGABeginElement(igaZ0,&elemZ0);CHKERRQ(ierr);

	while (IGANextElement(igaCS,elemCS)) 				
	{
		IGANextElement(igaZ0,elemZ0);

		ierr = IGAElementGetWorkMat(elemCS,&KlocCS);CHKERRQ(ierr);
		ierr = IGAElementGetWorkVec(elemCS,&FlocCS);CHKERRQ(ierr);
		ierr = IGAElementGetValues(elemZ0,arrayZ0CS,&Z0CS);CHKERRQ(ierr);

		// FormSystem loop
		while (IGAElementNextFormSystem(elemCS,&wtfCS,&wtf2CS))
		{
		// Quadrature loop
			ierr = IGAElementBeginPoint(elemCS,&pointCS);CHKERRQ(ierr);
			ierr = IGAElementBeginPoint(elemZ0,&pointZ0);CHKERRQ(ierr);

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

				if(pointCS->atboundary==0 && pointZ0->atboundary==0)
				{
					IGAElementNextPoint(elemZ0,pointZ0);
					ierr = IGAPointGetWorkMat(pointCS,&KpointCS);CHKERRQ(ierr);
					ierr = IGAPointGetWorkVec(pointCS,&FpointCS);CHKERRQ(ierr);
					ierr = CoupleStress(pointCS,pointZ0,KpointCS,FpointCS,Z0CS,NULL);CHKERRQ(ierr);
					ierr = IGAPointAddMat(pointCS,KpointCS,KlocCS);CHKERRQ(ierr);
					ierr = IGAPointAddVec(pointCS,FpointCS,FlocCS);CHKERRQ(ierr);
				}
			}
			while (pointZ0->index != -1)
			{
				IGAElementNextPoint(elemZ0,pointZ0);
			}

			ierr = IGAElementEndPoint(elemCS,&pointCS);CHKERRQ(ierr);
			ierr = IGAElementEndPoint(elemZ0,&pointZ0);CHKERRQ(ierr);
		}

		//ierr = IGAElementFixSystem(elemStress,KlocStress,FlocStress);CHKERRQ(ierr);					//This sets Dirichlet condition ¿? (Yes, this applies the conditions from IGASetBoundaryValue)
		ierr = IGAElementAssembleMat(elemCS,KlocCS,KCStress);CHKERRQ(ierr);
		ierr = IGAElementAssembleVec(elemCS,FlocCS,FCStress);CHKERRQ(ierr);
	}
	IGANextElement(igaZ0,elemZ0);

	ierr = IGAEndElement(igaCS,&elemCS);CHKERRQ(ierr);
	ierr = IGAEndElement(igaZ0,&elemZ0);CHKERRQ(ierr);

	// Restore local vectors u, Z0, Chi0 and arrays
	ierr = IGARestoreLocalVecArray(igaZ0,Z0,&localZ0CS,&arrayZ0CS);CHKERRQ(ierr);
	ierr = MatAssemblyBegin(KCStress,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyEnd  (KCStress,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

	ierr = VecAssemblyBegin(FCStress);CHKERRQ(ierr);
	ierr = VecAssemblyEnd  (FCStress);CHKERRQ(ierr);

	ierr = KSPSetOperators(kspCS,KCStress,KCStress);CHKERRQ(ierr);
	PC pcCS;
	ierr = KSPGetPC(kspCS,&pcCS); CHKERRQ(ierr);
	ierr = PCSetType(pcCS,PCLU); CHKERRQ(ierr);
	ierr = PCFactorSetMatSolverType(pcCS,MATSOLVERMUMPS); CHKERRQ(ierr);
	//ierr = KSPSetFromOptions(kspStress);CHKERRQ(ierr);
	//ierr = KSPSetTolerances(kspStress,1e-28,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);CHKERRQ(ierr);
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
	PetscPrintf(PETSC_COMM_WORLD,"System for Energy Density starting \n");
	IGA igaED;
	ierr = IGACreate(PETSC_COMM_WORLD,&igaED);CHKERRQ(ierr);
	ierr = IGASetDim(igaED,1);CHKERRQ(ierr);													//Spatial dimension of the problem
	ierr = IGASetDof(igaED,1);CHKERRQ(ierr);													//Number of degrees of freedom, per node
	ierr = IGASetOrder(igaED,2);CHKERRQ(ierr);												//Number of spatial derivatives to calculate
	ierr = IGASetFromOptions(igaED);CHKERRQ(ierr);											//Note: The order (or degree) of the shape functions is given by the mesh!
	ierr = IGARead(igaED,"./geometry31D.dat");CHKERRQ(ierr);
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

	IGAPoint		pointED;												//point
	IGAElement		elemED;													//element
	PetscReal		*KlocED,*FlocED;										//AA y BB
	PetscReal		*KpointED,*FpointED;									//KKK y FFF
	const PetscReal	*arrayZ0ED;												//arrayU
	Vec				localZ0ED;												//localU
	PetscReal		*Z0ED;													//U

  	IGAFormSystem	wtfED;
 	void			*wtf2ED;

 	KSP kspED;
	ierr = IGACreateKSP(igaED,&kspED);CHKERRQ(ierr);

	// Get local vectors Z0 and Chi0 and arrays
	ierr = IGAGetLocalVecArray(igaZ0,Z0,&localZ0ED,&arrayZ0ED);CHKERRQ(ierr);

	// Element loop
	ierr = IGABeginElement(igaED,&elemED);CHKERRQ(ierr);
	ierr = IGABeginElement(igaZ0,&elemZ0);CHKERRQ(ierr);

	while (IGANextElement(igaED,elemED)) 				
	{
		IGANextElement(igaZ0,elemZ0);

		ierr = IGAElementGetWorkMat(elemED,&KlocED);CHKERRQ(ierr);
		ierr = IGAElementGetWorkVec(elemED,&FlocED);CHKERRQ(ierr);
		ierr = IGAElementGetValues(elemZ0,arrayZ0ED,&Z0ED);CHKERRQ(ierr);

		// FormSystem loop
		while (IGAElementNextFormSystem(elemED,&wtfED,&wtf2ED))
		{
		// Quadrature loop
			ierr = IGAElementBeginPoint(elemED,&pointED);CHKERRQ(ierr);
			ierr = IGAElementBeginPoint(elemZ0,&pointZ0);CHKERRQ(ierr);

			while (IGAElementNextPoint(elemED,pointED))
			{
				if(pointED->atboundary==0 && pointZ0->atboundary==0)
				{
					IGAElementNextPoint(elemZ0,pointZ0);

					ierr = IGAPointGetWorkMat(pointED,&KpointED);CHKERRQ(ierr);
					ierr = IGAPointGetWorkVec(pointED,&FpointED);CHKERRQ(ierr);
					//PetscErrorCode EnergyDensity(IGAPoint p,IGAPoint pChi,IGAPoint pZu,PetscReal *K,PetscReal *F, PetscReal *Chi,PetscReal *Zu,void *ctx)		//When other fields like Pi or S (etc.) come into this, declare a IGAPoint pPi or pS and a new PestcReal *UPi or *US for each
					ierr = EnergyDensity(pointED,pointZ0,KpointED,FpointED,Z0ED,NULL);CHKERRQ(ierr);
					ierr = IGAPointAddMat(pointED,KpointED,KlocED);CHKERRQ(ierr);
					ierr = IGAPointAddVec(pointED,FpointED,FlocED);CHKERRQ(ierr);
				}
			}
			while (pointZ0->index != -1)
			{
				IGAElementNextPoint(elemZ0,pointZ0);
			}
			ierr = IGAElementEndPoint(elemED,&pointED);CHKERRQ(ierr);
			ierr = IGAElementEndPoint(elemZ0,&pointZ0);CHKERRQ(ierr);
		}

		//ierr = IGAElementFixSystem(elemStress,KlocStress,FlocStress);CHKERRQ(ierr);					//This sets Dirichlet condition ¿? (Yes, this applies the conditions from IGASetBoundaryValue)
		ierr = IGAElementAssembleMat(elemED,KlocED,KED);CHKERRQ(ierr);
		ierr = IGAElementAssembleVec(elemED,FlocED,FED);CHKERRQ(ierr);
	}
	IGANextElement(igaZ0,elemZ0);

	ierr = IGAEndElement(igaED,&elemED);CHKERRQ(ierr);
	ierr = IGAEndElement(igaZ0,&elemZ0);CHKERRQ(ierr);

	// Restore local vectors u, Z0, Chi0 and arrays
	ierr = IGARestoreLocalVecArray(igaZ0,Z0,&localZ0ED,&arrayZ0ED);CHKERRQ(ierr);

	ierr = MatAssemblyBegin(KED,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyEnd  (KED,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

	ierr = VecAssemblyBegin(FED);CHKERRQ(ierr);
	ierr = VecAssemblyEnd  (FED);CHKERRQ(ierr);

	ierr = KSPSetOperators(kspED,KED,KED);CHKERRQ(ierr);
	PC pcED;
	ierr = KSPGetPC(kspED,&pcED); CHKERRQ(ierr);
	ierr = PCSetType(pcED,PCLU); CHKERRQ(ierr);
	ierr = PCFactorSetMatSolverType(pcED,MATSOLVERMUMPS); CHKERRQ(ierr);
	//ierr = KSPSetFromOptions(kspStress);CHKERRQ(ierr);
	//ierr = KSPSetTolerances(kspStress,1e-28,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);CHKERRQ(ierr);
	ierr = KSPSolve(kspED,FED,ed0);CHKERRQ(ierr);

	ierr = KSPDestroy(&kspED);CHKERRQ(ierr);
	ierr = MatDestroy(&KED);CHKERRQ(ierr);
	ierr = VecDestroy(&FED);CHKERRQ(ierr);

	char nameED[]="/DensityEnergy.dat";
	char pathED[512];
	sprintf(pathED,"%s%s",direct,nameED);
	ierr = IGAWriteVec(igaED,ed0,pathED);CHKERRQ(ierr);	
//

//Destroy all objects not needed anymore (Better to do it here in case different codes call the same IGA, move if memory is a problem)

//

ierr = PetscFinalize();CHKERRQ(ierr);

return 0;
}