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

//System for L2 projection of Ue
	#undef  __FUNCT__
	#define __FUNCT__ "L2_Ue"
	PetscErrorCode L2_Ue(IGAPoint p,PetscReal *K,PetscReal *F,void *ctx)
	{
		const PetscReal *N0;
		IGAPointGetShapeFuns(p,0,(const PetscReal**)&N0);									//Value of the shape functions
		PetscInt a,b,i,j,r,s,nen=p->nen, dof=p->dof;

		PetscReal (*Keq)[dof][nen][dof] = (typeof(Keq)) K;
		PetscReal (*Feq)[dof] = (PetscReal (*)[dof])F;

		//Definition of alternating tensor
		const PetscReal e[3][3][3]=
		{
			{{0.0,0.0,0.0},{0.0,0.0,1.0},{0.0,-1.0,0.0}},
			{{0.0,0.0,-1.0},{0.0,0.0,0.0},{1.0,0.0,0.0}},
			{{0.0,1.0,0.0},{-1.0,0.0,0.0},{0.0,0.0,0.0}}
		};

		PetscReal eps=1.0/10.0;
		PetscReal v0[3]={0};
		PetscReal v1[3]={0};

		v0[2]=1.0; v1[2]=-1.0;

		PetscReal x[2];																		//Vector of reals, size equal to problem's dimension
		IGAPointFormGeomMap(p,x);															//Fills x with the coordinates of p, Gauss's point

		PetscReal Ue[3][3]={0};

		for(int i=0; i<3; i++)
		{
			for(int j=0; j<3; j++)
			{
				for(int k=0; k<3; k++)
				{
					Ue[i][j]+=e[i][j][k]*(v0[k]+0.5*(1.0+tanh(x[0]/eps))*(v1[k]-v0[k]));
				}
				//PetscPrintf(PETSC_COMM_WORLD,"Ue[%d][%d]=%f \n",i,j,Ue[i][j]);
			}	
		}

		PetscReal u[3][3]={0};
		PetscReal v[3][3]={0};

		if (p->atboundary)
		{
			return 0;
		}
		else
		{
			for (a=0; a<nen; a++) 
			{
				PetscReal Na= N0[a];

				for (b=0; b<nen; b++) 
				{
					PetscReal Nb= N0[b];

					for (i=0; i<dof; i++)
					{
						if (i==0)
						{
							u[0][0]=Na;  u[0][1]=0.0; u[0][2]=0.0;
							u[1][0]=0.0; u[1][1]=0.0; u[1][2]=0.0;
							u[2][0]=0.0; u[2][1]=0.0; u[2][2]=0.0;
						}
						else if(i==1)
						{
							u[0][0]=0.0; u[0][1]=Na;  u[0][2]=0.0;
							u[1][0]=0.0; u[1][1]=0.0; u[1][2]=0.0;
							u[2][0]=0.0; u[2][1]=0.0; u[2][2]=0.0;
						}
						else if(i==2)
						{
							u[0][0]=0.0; u[0][1]=0.0; u[0][2]=0.0;
							u[1][0]=Na;  u[1][1]=0.0; u[1][2]=0.0;
							u[2][0]=0.0; u[2][1]=0.0; u[2][2]=0.0;
						}
						else if(i==3)
						{
							u[0][0]=0.0; u[0][1]=0.0; u[0][2]=0.0;
							u[1][0]=0.0; u[1][1]=Na;  u[1][2]=0.0;
							u[2][0]=0.0; u[2][1]=0.0; u[2][2]=0.0;
						}

						for (j=0; j<dof; j++)
						{
							if (j==0)
							{
								v[0][0]=Nb;  v[0][1]=0.0; v[0][2]=0.0;
								v[1][0]=0.0; v[1][1]=0.0; v[1][2]=0.0;
								v[2][0]=0.0; v[2][1]=0.0; v[2][2]=0.0;
							}
							else if(j==1)
							{
								v[0][0]=0.0; v[0][1]=Nb;  v[0][2]=0.0;
								v[1][0]=0.0; v[1][1]=0.0; v[1][2]=0.0;
								v[2][0]=0.0; v[2][1]=0.0; v[2][2]=0.0;
							}
							else if(j==2)
							{
								v[0][0]=0.0; v[0][1]=0.0; v[0][2]=0.0;
								v[1][0]=Nb;  v[1][1]=0.0; v[1][2]=0.0;
								v[2][0]=0.0; v[2][1]=0.0; v[2][2]=0.0;
							}
							else if(j==3)
							{
								v[0][0]=0.0; v[0][1]=0.0; v[0][2]=0.0;
								v[1][0]=0.0; v[1][1]=Nb;  v[1][2]=0.0;
								v[2][0]=0.0; v[2][1]=0.0; v[2][2]=0.0;
							}

							Keq[a][i][b][j]=0.0;

							for (r=0; r<3; r++)
							{
								for (s=0; s<3; s++)
								{
									Keq[a][i][b][j]+=u[r][s]*v[r][s];
								}
							}
						}
					}
				}
			}

			for(a=0 ;a<nen; a++)
			{
				PetscReal Na=N0[a];

				for (i=0; i<dof; i++)
				{
					if (i==0)
						{
							u[0][0]=Na;  u[0][1]=0.0; u[0][2]=0.0;
							u[1][0]=0.0; u[1][1]=0.0; u[1][2]=0.0;
							u[2][0]=0.0; u[2][1]=0.0; u[2][2]=0.0;
						}
						else if(i==1)
						{
							u[0][0]=0.0; u[0][1]=Na;  u[0][2]=0.0;
							u[1][0]=0.0; u[1][1]=0.0; u[1][2]=0.0;
							u[2][0]=0.0; u[2][1]=0.0; u[2][2]=0.0;
						}
						else if(i==2)
						{
							u[0][0]=0.0; u[0][1]=0.0; u[0][2]=0.0;
							u[1][0]=Na;  u[1][1]=0.0; u[1][2]=0.0;
							u[2][0]=0.0; u[2][1]=0.0; u[2][2]=0.0;
						}
						else if(i==3)
						{
							u[0][0]=0.0; u[0][1]=0.0; u[0][2]=0.0;
							u[1][0]=0.0; u[1][1]=Na;  u[1][2]=0.0;
							u[2][0]=0.0; u[2][1]=0.0; u[2][2]=0.0;
						}

					Feq[a][i] = 0.0;

					for (r=0; r<3; r++)
					{
						for (s=0; s<3; s++)
						{
							Feq[a][i]+=Ue[r][s]*u[r][s];
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
	PetscErrorCode Z0sys(IGAPoint p,PetscReal *K,PetscReal *F,void *ctx)
	{
		PetscReal C[3][3][3][3]={0};														//Initialization of elastic tensor
		const PetscReal *N0,(*N1)[2],(*N2)[2][2];
		IGAPointGetShapeFuns(p,0,(const PetscReal**)&N0);									//Value of the shape functions
		IGAPointGetShapeFuns(p,1,(const PetscReal**)&N1);									//Derivatives of the shape functions
		IGAPointGetShapeFuns(p,2,(const PetscReal**)&N2);									//Second derivatives of the shape funcions
		//After this command Na_xx=N2[a][0][0], Na_yy=N2[a][1][1], Na_xy=N2[a][0][1], Na_yx[a][1][0] (these last two are equal) (remember a is the index of the shape function)
		PetscInt a,b,i,j,k,l,u,w,r,s,m,nen=p->nen, dof=p->dof;

		PetscReal (*Keq)[dof][nen][dof] = (typeof(Keq)) K;
		PetscReal (*Feq)[dof] = (PetscReal (*)[dof])F;

		//E and nu should come from AppCtx in the future
		//const PetscReal E=2100000.0*9.81*10000.0;					//[Pa]
		//const PetscReal nu=0.33;
		//const PetscReal lambda=(E*nu)/((1.0+nu)*(1.0-2.0*nu));
		//const PetscReal mu=E/(2.0*(1.0+nu));
		//const PetscReal eps=0.0*E/1000.0;							//Choose later based on whatever Amit says :)

		//Change for G=1
		const PetscReal nu=0.33;
		const PetscReal mu=1.0;
		const PetscReal lambda=2.0*mu*nu/(1.0-2.0*nu);
		const PetscReal eps=mu/100.0;								//Choose later based on whatever Amit says :)

		//////////////Delete this parte later, loads should come from appCtx or be 0
		PetscReal f[3]={0.0, 0.0, 0.0};			//Distributed load in body
		//PetscReal g[3]={0.0, 0.0, 0.0};			//Boundary load (applied wherever is defined in the p->atboundary block)
		PetscReal M[3]={0.0, 0.0, 0.0};			//Distributed moment in body
		PetscReal h[3]={0.0, 0.0, 1.0};			//Boundary moment (applied wherever is defined in the p->atboundary block)
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
			PetscReal x[2];																		//Vector of reals, size equal to problem's dimension
			IGAPointFormGeomMap(p,x);															//Fills x with the coordinates of p, Gauss's point

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
														-0.5*eps*(d2z[r][s][m]-d2z[s][r][m])*d2v[s][r][m];		//////////Revisar signo de esta linea!!!!!!
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
								Feq[a][i]+=M[r]*0.5*e[r][s][m]*dv[m][s];
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
	PetscErrorCode Stress(IGAPoint p,IGAPoint pUe,PetscReal *K,PetscReal *F,PetscReal *Uex,void *ctx)		//When other fields like Pi or S (etc.) come into this, declare a IGAPoint pPi or pS and a new PestcReal *UPi or *US for each
	{
		//const PetscReal *N0,(*N1)[2],(*N2)[2][2];
		//const PetscReal *N0,(*N1)[2];
		const PetscReal *N0;
		IGAPointGetShapeFuns(p,0,(const PetscReal**)&N0);									//Value of the shape functions
		//IGAPointGetShapeFuns(p,1,(const PetscReal**)&N1);									//Derivatives of the shape functions
		//IGAPointGetShapeFuns(p,2,(const PetscReal**)&N2);									//Second derivatives of the shape funcions
		PetscInt a,b,i,j,k,l,m,n,nen=p->nen, dof=p->dof;

		//Change to consider G=1
		const PetscReal nu=0.33;
		const PetscReal mu=1.0;
		const PetscReal lambda=2.0*mu*nu/(1.0-2.0*nu);
		const PetscReal eps=mu/100.0;															//Choose later based on whatever Amit says :)

		PetscReal Ue[4];															//Same for its 2nd order partial derivatives
		PetscReal d2_Ue[4][2][2];															//Same for its 2nd order partial derivatives
		IGAPointFormValue(pUe,Uex,&Ue[0]);												//Assign chi to its container
		IGAPointFormHess (pUe,Uex,&d2_Ue[0][0][0]);											//Same for the hessian (second derivatives)

		//Expanding z (and derivatives) to 3 components, more convenient for sums in for loops
		PetscReal full_Ue[3][3];
		full_Ue[0][0]=Ue[0];	full_Ue[0][1]=Ue[1];
		full_Ue[1][0]=Ue[2];	full_Ue[1][1]=Ue[3];

		PetscReal fulld2_Ue[3][3][3][3]={0};
		fulld2_Ue[0][0][0][0]=d2_Ue[0][0][0]; fulld2_Ue[0][0][0][1]=d2_Ue[0][0][1]; fulld2_Ue[0][0][1][0]=d2_Ue[0][1][0]; fulld2_Ue[0][0][1][1]=d2_Ue[0][1][1];
		fulld2_Ue[0][1][0][0]=d2_Ue[1][0][0]; fulld2_Ue[0][1][0][1]=d2_Ue[1][0][1]; fulld2_Ue[0][1][1][0]=d2_Ue[1][1][0]; fulld2_Ue[0][1][1][1]=d2_Ue[1][1][1];
		fulld2_Ue[1][0][0][0]=d2_Ue[2][0][0]; fulld2_Ue[1][0][0][1]=d2_Ue[2][0][1]; fulld2_Ue[1][0][1][0]=d2_Ue[2][1][0]; fulld2_Ue[1][0][1][1]=d2_Ue[2][1][1];
		fulld2_Ue[1][1][0][0]=d2_Ue[3][0][0]; fulld2_Ue[1][1][0][1]=d2_Ue[3][0][1]; fulld2_Ue[1][1][1][0]=d2_Ue[3][1][0]; fulld2_Ue[1][1][1][1]=d2_Ue[3][1][1];

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
							Fstress[a][i]+=C[m][n][k][l]*(full_Ue[k][l])*v[m][n];
						}

						Fstress[a][i]+= -0.25*eps*(fulld2_Ue[m][n][k][k]+fulld2_Ue[m][k][n][k]
							                      +fulld2_Ue[n][m][k][k]+fulld2_Ue[n][k][m][k])*v[m][n];
					}
				}
			}
		}
		return 0;
	}
//

//System for couple-stresses
	#undef  __FUNCT__
	#define __FUNCT__ "CoupleStress"
	//PetscErrorCode Stress(IGAPoint p,IGAPoint pU, IGAPoint pHs,IGAPoint pChi,IGAPoint pZu,PetscReal *K,PetscReal *F,PetscReal *U,PetscReal *HS, PetscReal *Chi,PetscReal *Zu,void *ctx)		//When other fields like Pi or S (etc.) come into this, declare a IGAPoint pPi or pS and a new PescReal *UPi or *US for each
	PetscErrorCode CoupleStress(IGAPoint p,IGAPoint pUe,PetscReal *K,PetscReal *F,PetscReal *U,void *ctx)		//When other fields like Pi or S (etc.) come into this, declare a IGAPoint pPi or pS and a new PestcReal *UPi or *US for each
	{
		const PetscReal *N0;
		IGAPointGetShapeFuns(p,0,(const PetscReal**)&N0);									//Value of the shape functions
		PetscInt a,b,i,j,k,l,m,n,nen=p->nen, dof=p->dof;

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
		
		PetscReal d_Ue[2][2][2];															//Same for its 2nd order partial derivatives
		IGAPointFormGrad (pUe,U,&d_Ue[0][0][0]);											//Same for the hessian (second derivatives)

		//Expanding Ue (and derivatives) to 3 components, more convenient for sums in for loops
		PetscReal fulld_Ue[3][3][3]={0};

		fulld_Ue[0][0][0]=d_Ue[0][0][0]; fulld_Ue[0][0][1]=d_Ue[0][0][1];
		fulld_Ue[0][1][0]=d_Ue[0][1][0]; fulld_Ue[0][1][1]=d_Ue[0][1][1];
		fulld_Ue[1][0][0]=d_Ue[1][0][0]; fulld_Ue[1][0][1]=d_Ue[1][0][1];
		fulld_Ue[1][1][0]=d_Ue[1][1][0]; fulld_Ue[1][1][1]=d_Ue[1][1][1];

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
							FCS[a][i]+=-0.5*eps*e[m][k][l]*(d_Ue[k][l][n]+d_Ue[k][n][l])*v[m][n];
						}
					}
				}
			}
		}
		return 0;
	}
//

//System for Energy
	#undef  __FUNCT__
	#define __FUNCT__ "Energy"
	//PetscErrorCode Stress(IGAPoint p,IGAPoint pU, IGAPoint pHs,IGAPoint pChi,IGAPoint pZu,PetscReal *K,PetscReal *F,PetscReal *U,PetscReal *HS, PetscReal *Chi,PetscReal *Zu,void *ctx)		//When other fields like Pi or S (etc.) come into this, declare a IGAPoint pPi or pS and a new PescReal *UPi or *US for each
	PetscErrorCode Energy(IGAPoint p,IGAPoint pUe,PetscReal *K,PetscReal *F,PetscReal *U,void *ctx)		//When other fields like Pi or S (etc.) come into this, declare a IGAPoint pPi or pS and a new PestcReal *UPi or *US for each
	{
		const PetscReal *N0;
		IGAPointGetShapeFuns(p,0,(const PetscReal**)&N0);									//Value of the shape functions
		PetscInt a,b,i,j,k,l,m,n,nen=p->nen, dof=p->dof;

		//Change to consider G=1
		const PetscReal nu=0.33;
		const PetscReal mu=1.0;
		const PetscReal lambda=2.0*mu*nu/(1.0-2.0*nu);
		const PetscReal eps=mu/100.0;															//Choose later based on whatever Amit says :)

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

		//Definition of alternating tensor
		const PetscReal e[3][3][3]=
		{
			{{0.0,0.0,0.0},{0.0,0.0,1.0},{0.0,-1.0,0.0}},
			{{0.0,0.0,-1.0},{0.0,0.0,0.0},{1.0,0.0,0.0}},
			{{0.0,1.0,0.0},{-1.0,0.0,0.0},{0.0,0.0,0.0}}
		};
		
		PetscReal Ue[2][2], d_Ue[2][2][2];
		IGAPointFormValue(pUe,U,&Ue[0][0]);
		IGAPointFormGrad (pUe,U,&d_Ue[0][0][0]);

		//Expanding Ue (and derivatives) to 3 components, more convenient for sums in for loops
		PetscReal full_Ue[3][3]={0};

		full_Ue[0][0]=Ue[0][0]; full_Ue[0][0]=Ue[0][0];
		full_Ue[0][1]=Ue[0][1]; full_Ue[0][1]=Ue[0][1];
		full_Ue[1][0]=Ue[1][0]; full_Ue[1][0]=Ue[1][0];
		full_Ue[1][1]=Ue[1][1]; full_Ue[1][1]=Ue[1][1];

		PetscReal fulld_Ue[3][3][3]={0};

		fulld_Ue[0][0][0]=d_Ue[0][0][0]; fulld_Ue[0][0][1]=d_Ue[0][0][1];
		fulld_Ue[0][1][0]=d_Ue[0][1][0]; fulld_Ue[0][1][1]=d_Ue[0][1][1];
		fulld_Ue[1][0][0]=d_Ue[1][0][0]; fulld_Ue[1][0][1]=d_Ue[1][0][1];
		fulld_Ue[1][1][0]=d_Ue[1][1][0]; fulld_Ue[1][1][1]=d_Ue[1][1][1];

		PetscReal (*KE)[dof][nen][dof] = (typeof(KE)) K;
		PetscReal (*FE)[dof] = (PetscReal (*)[dof]) F;

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
						KE[a][i][b][i]=N0[a]*N0[b];
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
					
					FE[a][i]=0.0;
					for (k=0;k<3;k++)
					{
						for(l=0;l<3;l++)
						{
							for(m=0;m<3;m++)
							{
								for(n=0;n<3;n++)
								{
									FE[a][i]+=0.5*full_Ue[k][l]*C[k][l][m][n]*full_Ue[m][n]*N0[a];
									//FCS[a][i]+=-0.5*eps*e[m][k][l]*(d_Ue[k][l][n]+d_Ue[k][n][l])*v[m][n];
								}
							}
						}
					}
					for (k=0;k<3;k++)
					{
						for(l=0;l<3;l++)
						{
							for(m=0;m<3;m++)
							{
								FE[a][i]+=0.5*eps*fulld_Ue[k][l][m]*fulld_Ue[k][l][m]*N0[a];
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
	PetscPrintf(PETSC_COMM_WORLD,"Start of PruebaTest \n");

	PetscInt commsize,rank;
	ierr = MPI_Comm_size(PETSC_COMM_WORLD,&commsize);CHKERRQ(ierr);
	ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank);CHKERRQ(ierr);
//

//App context creation and some data
	//Mesh parameters (to fix specific points in z0 system)
	PetscInt b=201;				//Parmeter to choose size of cores, must always be odd, core will be of size 1 unit, rest of the body will be of size b-1 units in each direction
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
		ierr=system("cp ~/shared/CodigosPetIGA/PruebaTest.c ~/sharedResults/");
	}
//

//Creation of types and systems for the L2 projection of Ue
	PetscPrintf(PETSC_COMM_WORLD,"System for Ue starting \n");
	IGA igaUe;
	ierr = IGACreate(PETSC_COMM_WORLD,&igaUe);CHKERRQ(ierr);
	ierr = IGASetDim(igaUe,2);CHKERRQ(ierr);														//Spatial dimension of the problem
	ierr = IGASetDof(igaUe,4);CHKERRQ(ierr);														//Number of degrees of freedom, per node
	ierr = IGASetOrder(igaUe,3);CHKERRQ(ierr);													//Number of spatial derivatives to calculate
	ierr = IGASetFromOptions(igaUe);CHKERRQ(ierr);												//Note: The order (or degree) of the shape functions is given by the mesh!
	ierr = IGARead(igaUe,"./geometry3.dat");CHKERRQ(ierr);
	ierr = IGASetUp(igaUe);CHKERRQ(ierr);
	
	for (dir=0; dir<2; dir++) 
	{
		for (side=0; side<2; side++) 
		{
			//ierr = IGASetBoundaryValue(igaS,dir,side,dof,0.0);CHKERRQ(ierr);    				// Dirichlet boundary conditions
			ierr = IGASetBoundaryForm(igaUe,dir,side,PETSC_TRUE);CHKERRQ(ierr);  				// Neumann boundary conditions
		}
	}

	Vec ue0;
	Mat Kl2Ue;
	Vec Fl2Ue;
	ierr = IGACreateVec(igaUe,&ue0);CHKERRQ(ierr);  
	ierr = IGACreateMat(igaUe,&Kl2Ue);CHKERRQ(ierr);
	ierr = IGACreateVec(igaUe,&Fl2Ue);CHKERRQ(ierr);
	ierr = IGASetFormSystem(igaUe,L2_Ue,&userL2);CHKERRQ(ierr);
	ierr = IGAComputeSystem(igaUe,Kl2Ue,Fl2Ue);CHKERRQ(ierr);

	//This parts set and calls KSP to solve the linear system
	KSP kspl2Ue;
	ierr = IGACreateKSP(igaUe,&kspl2Ue);CHKERRQ(ierr);										
	ierr = KSPSetOperators(kspl2Ue,Kl2Ue,Kl2Ue);CHKERRQ(ierr); 								//This function creates the matrix for the system on the second parameter and uses the 3rd parameter as a preconditioner
	ierr = KSPSetType(kspl2Ue,KSPCG);CHKERRQ(ierr);											//Using KSPCG (conjugated gradient) because the matrix is symmetric
	ierr = KSPSetOptionsPrefix(kspl2Ue,"l2pUe_");CHKERRQ(ierr);
	ierr = KSPSetFromOptions(kspl2Ue);CHKERRQ(ierr);
	ierr = KSPSetTolerances(kspl2Ue,1e-30,1e-15,PETSC_DEFAULT,PETSC_DEFAULT);CHKERRQ(ierr);
	ierr = KSPSolve(kspl2Ue,Fl2Ue,ue0);CHKERRQ(ierr);											//This is a simple system, so it can be solved with just this command

	ierr = KSPDestroy(&kspl2Ue);CHKERRQ(ierr);
	ierr = MatDestroy(&Kl2Ue);CHKERRQ(ierr);
	ierr = VecDestroy(&Fl2Ue);CHKERRQ(ierr);
	char nameUe[]="/Ue-2d-0.dat";
	char pathUe[512];
	sprintf(pathUe,"%s%s",direct,nameUe);
	ierr = IGAWriteVec(igaUe,ue0,pathUe);CHKERRQ(ierr);
//

/*
//System for initial state of z0
	PetscPrintf(PETSC_COMM_WORLD,"System for Z0 starting \n");
	IGA igaZ0;
	ierr = IGACreate(PETSC_COMM_WORLD,&igaZ0);CHKERRQ(ierr);
	ierr = IGASetDim(igaZ0,2);CHKERRQ(ierr);													//Spatial dimension of the problem
	ierr = IGASetDof(igaZ0,2);CHKERRQ(ierr);													//Number of degrees of freedom, per node
	ierr = IGASetOrder(igaZ0,2);CHKERRQ(ierr);													//Number of spatial derivatives to calculate
	ierr = IGASetFromOptions(igaZ0);CHKERRQ(ierr);												//Note: The order (or degree) of the shape functions is given by the mesh!
	ierr = IGARead(igaZ0,"./geometry2.dat");CHKERRQ(ierr);
	ierr = IGASetUp(igaZ0);CHKERRQ(ierr);

	PetscInt fijaPunto=0;

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
		//ierr = IGASetBoundaryValue(iga,dir,side,dof,0.0);CHKERRQ(ierr);					// Dirichlet boundary conditions
		ierr = IGASetBoundaryValue(igaZ0,1,0,0,0.0);CHKERRQ(ierr);
		ierr = IGASetBoundaryValue(igaZ0,1,0,1,0.0);CHKERRQ(ierr);
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

	if (fijaPunto==1)
	{
		//Here we set values to the Matrix directly, to impose Dirichlet condition in a single point. This is valid for C1 nurbs
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
*/

//System for L2 projection of stress
	PetscPrintf(PETSC_COMM_WORLD,"System for Stress starting \n");
	IGA igaStress;
	ierr = IGACreate(PETSC_COMM_WORLD,&igaStress);CHKERRQ(ierr);
	ierr = IGASetDim(igaStress,2);CHKERRQ(ierr);													//Spatial dimension of the problem
	ierr = IGASetDof(igaStress,4);CHKERRQ(ierr);													//Number of degrees of freedom, per node
	ierr = IGASetOrder(igaStress,2);CHKERRQ(ierr);													//Number of spatial derivatives to calculate
	ierr = IGASetFromOptions(igaStress);CHKERRQ(ierr);												//Note: The order (or degree) of the shape functions is given by the mesh!
	ierr = IGARead(igaStress,"./geometry3.dat");CHKERRQ(ierr);
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

	IGAPoint		pointStress, pointUe;															//point
	IGAElement		elemStress, elemUe;																//element
	PetscReal		*KlocStress,*FlocStress;												//AA y BB
	PetscReal		*KpointStress,*FpointStress;											//KKK y FFF
	const PetscReal	*arrayUeStress;															//arrayU
	Vec				localUeStress;															//localU
	PetscReal		*UeStress;																//U

  	IGAFormSystem	wtfStress;
 	void			*wtf2Stress;

 	KSP kspStress;
	ierr = IGACreateKSP(igaStress,&kspStress);CHKERRQ(ierr);

	// Get local vectors Z0
	ierr = IGAGetLocalVecArray(igaUe,ue0,&localUeStress,&arrayUeStress);CHKERRQ(ierr);

	// Element loop
	ierr = IGABeginElement(igaStress,&elemStress);CHKERRQ(ierr);
	ierr = IGABeginElement(igaUe,&elemUe);CHKERRQ(ierr);

	while (IGANextElement(igaStress,elemStress)) 
	{
		IGANextElement(igaUe,elemUe);

		ierr = IGAElementGetWorkMat(elemStress,&KlocStress);CHKERRQ(ierr);
		ierr = IGAElementGetWorkVec(elemStress,&FlocStress);CHKERRQ(ierr);
		ierr = IGAElementGetValues(elemUe,arrayUeStress,&UeStress);CHKERRQ(ierr);

		// FormSystem loop
		while (IGAElementNextFormSystem(elemStress,&wtfStress,&wtf2Stress)) 
		{
		// Quadrature loop
			ierr = IGAElementBeginPoint(elemStress,&pointStress);CHKERRQ(ierr);
			ierr = IGAElementBeginPoint(elemUe,&pointUe);CHKERRQ(ierr);

			while (IGAElementNextPoint(elemStress,pointStress))
			{
				if(pointStress->atboundary==1)
				{
					//IGAElementNextPoint(elemZ0,pointZ0);
				}

				if(pointUe->atboundary==1)
				{
					PetscPrintf(PETSC_COMM_WORLD,"Hola1 en stress \n");
				}

				if(pointStress->atboundary==0 && pointUe->atboundary==0)
				{
					IGAElementNextPoint(elemUe,pointUe);

					ierr = IGAPointGetWorkMat(pointStress,&KpointStress);CHKERRQ(ierr);
					ierr = IGAPointGetWorkVec(pointStress,&FpointStress);CHKERRQ(ierr);
					//PetscErrorCode Stress(IGAPoint p,IGAPoint pChi,IGAPoint pZu,PetscReal *K,PetscReal *F, PetscReal *Chi,PetscReal *Zu,void *ctx)
					ierr = Stress(pointStress,pointUe,KpointStress,FpointStress,UeStress,NULL);CHKERRQ(ierr);
					ierr = IGAPointAddMat(pointStress,KpointStress,KlocStress);CHKERRQ(ierr);
					ierr = IGAPointAddVec(pointStress,FpointStress,FlocStress);CHKERRQ(ierr);
				}
			}
			while (pointUe->index != -1)
			{
				IGAElementNextPoint(elemUe,pointUe);
			}

			ierr = IGAElementEndPoint(elemStress,&pointStress);CHKERRQ(ierr);
			ierr = IGAElementEndPoint(elemUe,&pointUe);CHKERRQ(ierr);
		}

		//ierr = IGAElementFixSystem(elemStress,KlocStress,FlocStress);CHKERRQ(ierr);					//This sets Dirichlet condition ¿? (Yes, this applies the conditions from IGASetBoundaryValue)
		ierr = IGAElementAssembleMat(elemStress,KlocStress,KStress);CHKERRQ(ierr);
		ierr = IGAElementAssembleVec(elemStress,FlocStress,FStress);CHKERRQ(ierr);
	}
	IGANextElement(igaUe,elemUe);

	ierr = IGAEndElement(igaStress,&elemStress);CHKERRQ(ierr);
	ierr = IGAEndElement(igaUe,&elemUe);CHKERRQ(ierr);
	// Restore local vectors u, Z0, Chi0 and arrays
	ierr = IGARestoreLocalVecArray(igaUe,ue0,&localUeStress,&arrayUeStress);CHKERRQ(ierr);
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

//System for L2 projection of couple stress
	PetscPrintf(PETSC_COMM_WORLD,"System for CoupleStress starting \n");
	IGA igaCS;
	ierr = IGACreate(PETSC_COMM_WORLD,&igaCS);CHKERRQ(ierr);
	ierr = IGASetDim(igaCS,2);CHKERRQ(ierr);													//Spatial dimension of the problem
	ierr = IGASetDof(igaCS,2);CHKERRQ(ierr);													//Number of degrees of freedom, per node
	ierr = IGASetOrder(igaCS,2);CHKERRQ(ierr);													//Number of spatial derivatives to calculate
	ierr = IGASetFromOptions(igaCS);CHKERRQ(ierr);												//Note: The order (or degree) of the shape functions is given by the mesh!
	ierr = IGARead(igaCS,"./geometry3.dat");CHKERRQ(ierr);
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
	const PetscReal	*arrayUeCS;			//arrayU
	Vec				localUeCS;				//localU
	PetscReal		*UeCS;											//U

  	IGAFormSystem	wtfCS;
 	void			*wtf2CS;

 	KSP kspCS;
	ierr = IGACreateKSP(igaCS,&kspCS);CHKERRQ(ierr);

	// Get local vectors Z0 and Chi0 and arrays
	ierr = IGAGetLocalVecArray(igaUe,ue0,&localUeCS,&arrayUeCS);CHKERRQ(ierr);

	// Element loop
	ierr = IGABeginElement(igaCS,&elemCS);CHKERRQ(ierr);
	ierr = IGABeginElement(igaUe,&elemUe);CHKERRQ(ierr);

	while (IGANextElement(igaCS,elemCS)) 				
	{
		IGANextElement(igaUe,elemUe);

		ierr = IGAElementGetWorkMat(elemCS,&KlocCS);CHKERRQ(ierr);
		ierr = IGAElementGetWorkVec(elemCS,&FlocCS);CHKERRQ(ierr);
		ierr = IGAElementGetValues(elemUe,arrayUeCS,&UeCS);CHKERRQ(ierr);

		// FormSystem loop
		while (IGAElementNextFormSystem(elemCS,&wtfCS,&wtf2CS))
		{
		// Quadrature loop
			ierr = IGAElementBeginPoint(elemCS,&pointCS);CHKERRQ(ierr);
			ierr = IGAElementBeginPoint(elemUe,&pointUe);CHKERRQ(ierr);

			while (IGAElementNextPoint(elemCS,pointCS))
			{
				if(pointCS->atboundary==1)
				{
					//delete this block??
				}

				if(pointUe->atboundary==1)
				{
					//delete this block??
				}

				if(pointCS->atboundary==0 && pointUe->atboundary==0)
				{
					IGAElementNextPoint(elemUe,pointUe);
					ierr = IGAPointGetWorkMat(pointCS,&KpointCS);CHKERRQ(ierr);
					ierr = IGAPointGetWorkVec(pointCS,&FpointCS);CHKERRQ(ierr);
					ierr = CoupleStress(pointCS,pointUe,KpointCS,FpointCS,UeCS,NULL);CHKERRQ(ierr);
					ierr = IGAPointAddMat(pointCS,KpointCS,KlocCS);CHKERRQ(ierr);
					ierr = IGAPointAddVec(pointCS,FpointCS,FlocCS);CHKERRQ(ierr);
				}
			}
			while (pointUe->index != -1)
			{
				IGAElementNextPoint(elemUe,pointUe);
			}
			ierr = IGAElementEndPoint(elemCS,&pointCS);CHKERRQ(ierr);
			ierr = IGAElementEndPoint(elemUe,&pointUe);CHKERRQ(ierr);
		}

		//ierr = IGAElementFixSystem(elemStress,KlocStress,FlocStress);CHKERRQ(ierr);					//This sets Dirichlet condition ¿? (Yes, this applies the conditions from IGASetBoundaryValue)
		ierr = IGAElementAssembleMat(elemCS,KlocCS,KCStress);CHKERRQ(ierr);
		ierr = IGAElementAssembleVec(elemCS,FlocCS,FCStress);CHKERRQ(ierr);
	}
	IGANextElement(igaUe,elemUe);

	ierr = IGAEndElement(igaCS,&elemCS);CHKERRQ(ierr);
	ierr = IGAEndElement(igaUe,&elemUe);CHKERRQ(ierr);

	// Restore local vectors u, Z0, Chi0 and arrays
	ierr = IGARestoreLocalVecArray(igaUe,ue0,&localUeCS,&arrayUeCS);CHKERRQ(ierr);
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

	ierr = KSPDestroy(&kspCS);CHKERRQ(ierr);
	ierr = MatDestroy(&KCStress);CHKERRQ(ierr);
	ierr = VecDestroy(&FCStress);CHKERRQ(ierr);

	char nameCStress[]="/lambda-2d-0.dat";
	char pathCStress[512];
	sprintf(pathCStress,"%s%s",direct,nameCStress);
	ierr = IGAWriteVec(igaCS,lambda0,pathCStress);CHKERRQ(ierr);	
//

//System for L2 projection of Energy
	PetscPrintf(PETSC_COMM_WORLD,"System for Energy starting \n");
	IGA igaE;
	ierr = IGACreate(PETSC_COMM_WORLD,&igaE);CHKERRQ(ierr);
	ierr = IGASetDim(igaE,2);CHKERRQ(ierr);													//Spatial dimension of the problem
	ierr = IGASetDof(igaE,1);CHKERRQ(ierr);													//Number of degrees of freedom, per node
	ierr = IGASetOrder(igaE,2);CHKERRQ(ierr);													//Number of spatial derivatives to calculate
	ierr = IGASetFromOptions(igaE);CHKERRQ(ierr);												//Note: The order (or degree) of the shape functions is given by the mesh!
	ierr = IGARead(igaE,"./geometry3.dat");CHKERRQ(ierr);
	ierr = IGASetUp(igaE);CHKERRQ(ierr);

	for (dir=0; dir<2; dir++) 
	{
		for (side=0; side<2; side++) 
		{
			//ierr = IGASetBoundaryValue(iga,dir,side,dof,0.0);CHKERRQ(ierr);					// Dirichlet boundary conditions
			ierr = IGASetBoundaryForm(igaE,dir,side,PETSC_TRUE);CHKERRQ(ierr);  				// Neumann boundary conditions
		}
	}

	Mat KE;
	Vec e0,FE;

	ierr = IGACreateMat(igaE,&KE);CHKERRQ(ierr);
	ierr = IGACreateVec(igaE,&e0);CHKERRQ(ierr);
	ierr = IGACreateVec(igaE,&FE);CHKERRQ(ierr);

	IGAPoint		pointE;															//point
	IGAElement		elemE;																//element
	PetscReal		*KlocE,*FlocE;												//AA y BB
	PetscReal		*KpointE,*FpointE;											//KKK y FFF
	const PetscReal	*arrayUeE;			//arrayU
	Vec				localUeE;				//localU
	PetscReal		*UeE;											//U

  	IGAFormSystem	wtfE;
 	void			*wtf2E;

 	KSP kspE;
	ierr = IGACreateKSP(igaE,&kspE);CHKERRQ(ierr);

	// Get local vectors Z0 and Chi0 and arrays
	ierr = IGAGetLocalVecArray(igaUe,ue0,&localUeE,&arrayUeE);CHKERRQ(ierr);

	// Element loop
	ierr = IGABeginElement(igaE,&elemE);CHKERRQ(ierr);
	ierr = IGABeginElement(igaUe,&elemUe);CHKERRQ(ierr);

	while (IGANextElement(igaE,elemE)) 				
	{
		IGANextElement(igaUe,elemUe);

		ierr = IGAElementGetWorkMat(elemE,&KlocE);CHKERRQ(ierr);
		ierr = IGAElementGetWorkVec(elemE,&FlocE);CHKERRQ(ierr);
		ierr = IGAElementGetValues(elemUe,arrayUeE,&UeE);CHKERRQ(ierr);

		// FormSystem loop
		while (IGAElementNextFormSystem(elemE,&wtfE,&wtf2E))
		{
		// Quadrature loop
			ierr = IGAElementBeginPoint(elemE,&pointE);CHKERRQ(ierr);
			ierr = IGAElementBeginPoint(elemUe,&pointUe);CHKERRQ(ierr);

			while (IGAElementNextPoint(elemE,pointE))
			{
				if(pointE->atboundary==1)
				{
					//delete this block??
				}

				if(pointUe->atboundary==1)
				{
					//delete this block??
				}

				if(pointE->atboundary==0 && pointUe->atboundary==0)
				{
					IGAElementNextPoint(elemUe,pointUe);
					ierr = IGAPointGetWorkMat(pointE,&KpointE);CHKERRQ(ierr);
					ierr = IGAPointGetWorkVec(pointE,&FpointE);CHKERRQ(ierr);
					ierr = Energy(pointE,pointUe,KpointE,FpointE,UeE,NULL);CHKERRQ(ierr);
					ierr = IGAPointAddMat(pointE,KpointE,KlocE);CHKERRQ(ierr);
					ierr = IGAPointAddVec(pointE,FpointE,FlocE);CHKERRQ(ierr);
				}
			}
			while (pointUe->index != -1)
			{
				IGAElementNextPoint(elemUe,pointUe);
			}
			ierr = IGAElementEndPoint(elemE,&pointE);CHKERRQ(ierr);
			ierr = IGAElementEndPoint(elemUe,&pointUe);CHKERRQ(ierr);
		}

		//ierr = IGAElementFixSystem(elemStress,KlocStress,FlocStress);CHKERRQ(ierr);					//This sets Dirichlet condition ¿? (Yes, this applies the conditions from IGASetBoundaryValue)
		ierr = IGAElementAssembleMat(elemE,KlocE,KE);CHKERRQ(ierr);
		ierr = IGAElementAssembleVec(elemE,FlocE,FE);CHKERRQ(ierr);
	}
	IGANextElement(igaUe,elemUe);

	ierr = IGAEndElement(igaE,&elemE);CHKERRQ(ierr);
	ierr = IGAEndElement(igaUe,&elemUe);CHKERRQ(ierr);

	// Restore local vectors u, Z0, Chi0 and arrays
	ierr = IGARestoreLocalVecArray(igaUe,ue0,&localUeE,&arrayUeE);CHKERRQ(ierr);
	ierr = MatAssemblyBegin(KE,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyEnd  (KE,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

	ierr = VecAssemblyBegin(FE);CHKERRQ(ierr);
	ierr = VecAssemblyEnd  (FE);CHKERRQ(ierr);

	ierr = KSPSetOperators(kspE,KE,KE);CHKERRQ(ierr);
	PC pcE;
	ierr = KSPGetPC(kspE,&pcE); CHKERRQ(ierr);
	ierr = PCSetType(pcE,PCLU); CHKERRQ(ierr);
	ierr = PCFactorSetMatSolverType(pcE,MATSOLVERMUMPS); CHKERRQ(ierr);
	//ierr = KSPSetFromOptions(kspStress);CHKERRQ(ierr);
	//ierr = KSPSetTolerances(kspStress,1e-28,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);CHKERRQ(ierr);
	ierr = KSPSolve(kspE,FE,e0);CHKERRQ(ierr);

	ierr = KSPDestroy(&kspE);CHKERRQ(ierr);
	ierr = MatDestroy(&KE);CHKERRQ(ierr);
	ierr = VecDestroy(&FE);CHKERRQ(ierr);

	char nameE[]="/Energia.dat";
	char pathE[512];
	sprintf(pathE,"%s%s",direct,nameE);
	ierr = IGAWriteVec(igaE,e0,pathE);CHKERRQ(ierr);	
//

//Destroy all objects not needed anymore (Better to do it here in case different codes call the same IGA, move if memory is a problem)
	ierr = IGADestroy(&igaUe);CHKERRQ(ierr);
	//ierr = IGADestroy(&igaStress);CHKERRQ(ierr);
	//ierr = IGADestroy(&igaZ0);CHKERRQ(ierr);
//

T=time(NULL);			//quitar si da problemas
tm=*localtime(&T);
PetscPrintf(PETSC_COMM_WORLD,"Current time is %02d:%02d:%02d \n",tm.tm_hour,tm.tm_min,tm.tm_sec);
ierr = PetscFinalize();CHKERRQ(ierr);

return 0;
}