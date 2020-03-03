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


//L2 projection of Pi(0)
	#undef  __FUNCT__
	#define __FUNCT__ "L2ProjectionPi"
	PetscErrorCode L2ProjectionPi(IGAPoint p,PetscReal *K,PetscReal *F,void *ctx)
	{
		if (p->atboundary)
		{
			return 0;																		//Zero Neumann boundary condition
		}

		AppCtxL2    *user  = (AppCtxL2 *)ctx;
		PetscReal Lx     = user->Lx;
		PetscReal Ly     = user->Ly;
		PetscReal nx     = user->nx;
		PetscReal ny     = user->ny;

		PetscInt a,b,i;
		PetscInt nen = p->nen;																//Number of shape functions
		PetscInt dim = p->dim;																//Spatial dimensions of the problem
		PetscInt dof = p->dof;

		PetscReal x[dim];																	//Vector of reals, size equal to problem's dimension
		IGAPointFormGeomMap(p,x);															//Fills x with the coordinates of p, Gauss's point

		//g is the function to L2 project
		PetscReal dx=Lx/nx;
		PetscReal dy=Ly/ny;
		PetscReal eps=0.005;
		PetscReal g[dof];
		for (i=0; i<dof; i++)
		{
			if (i==0)
			{
				g[i]=1.917545*(0.5*(tanh((x[0]+dx)/eps)+1.0)-0.5*(tanh((x[0]-dx)/eps)+1.0))*(0.5*(tanh((x[1]+dy)/eps)+1.0)-0.5*(tanh((x[1]-dy)/eps)+1.0));
			}
			else if (i==1)
			{
				g[i]=0;
			}
			else if (i==2)
			{
				g[i]=1.917545*(0.5*(tanh((x[0]+10*dx)/eps)+1.0)-0.5*(tanh((x[0]+8*dx)/eps)+1.0))*(0.5*(tanh((x[1]+6*dy)/eps)+1.0)-0.5*(tanh((x[1]+4*dy)/eps)+1.0));
			}
			else
			{
				g[i]=1.917545*(0.5*(tanh((x[0]-dx)/eps)+1.0)-0.5*(tanh((x[0]+dx)/eps)+1.0))*(0.5*(tanh((x[1]-7*dy)/eps)+1.0)-0.5*(tanh((x[1]-5*dy)/eps)+1.0));
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
		  		FF[a][i] = N[a]*g[i];
			}
		}
		return 0;
	}
//

//L2 projection of Al(0)
	#undef  __FUNCT__
	#define __FUNCT__ "L2ProjectionAl"
	PetscErrorCode L2ProjectionAl(IGAPoint p,PetscReal *K,PetscReal *F,void *ctx)
	{
		if (p->atboundary)
		{
			return 0;																		//Pure Neumann condition 
		}

		AppCtxL2    *user  = (AppCtxL2 *)ctx;
		PetscReal Lx     = user->Lx;
		PetscReal Ly     = user->Ly;
		PetscReal nx     = user->nx;
		PetscReal ny     = user->ny;

		PetscInt a,b,i;
		PetscInt nen = p->nen;																//Number of shape functions, 9 in this case.
		PetscInt dim = p->dim;																//Number of spatial dimensions
		PetscInt dof = p->dof;

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
				g[i]=0*0.9585551*0.5*(tanh((x[0]+dx)/eps)-tanh((x[0]-dx)/eps))*(tanh((x[1]+dy)/eps)-tanh((x[1]-dy)/eps));
			}
			else if (i==1)
			{
				//g[i]=1.917545*(0.5*(tanh((x[0]+dx)/eps)+1.0)-0.5*(tanh((x[0]-dx)/eps)+1.0))*(0.5*(tanh((x[1]+dy)/eps)+1.0)-0.5*(tanh((x[1]-dy)/eps)+1.0));
				g[i]=0*0.9585551*0.5*(tanh((x[0]+dx)/eps)-tanh((x[0]-dx)/eps))*(tanh((x[1]+dy)/eps)-tanh((x[1]-dy)/eps));
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
		  		FF[a][i] = N[a]*g[i];
			}
		}
		return 0;
	}
//

//System for L2 projection of S(0)
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
		PetscReal nx     = user->nx;
		PetscReal ny     = user->ny;

		PetscInt a,b,i;
		PetscInt nen = p->nen;																//Number of shape functions
		PetscInt dim = p->dim;																//Spatial dimensions of the problem
		PetscInt dof = p->dof;

		PetscReal x[dim];																	//Vector of reals, size equal to problem's dimension
		IGAPointFormGeomMap(p,x);															//Fills x with the coordinates of p, Gauss's point

		//g is the function to L2 project
		PetscReal dx=Lx/nx;
		PetscReal dy=Ly/ny;
		PetscReal eps=0.005;
		PetscReal g[dof];
		//g[i]=1.917545*(0.5*(tanh((x[0]+dx)/eps)+1.0)-0.5*(tanh((x[0]-dx)/eps)+1.0))*(0.5*(tanh((x[1]+dy)/eps)+1.0)-0.5*(tanh((x[1]-dy)/eps)+1.0));
		g[0]=0;
		g[1]=-1.917545*(0.5*(tanh((x[0]-10*dx)/eps)+1.0)-0.5*(tanh((x[0]+8*dx)/eps)+1.0))*(0.5*(tanh((x[1]+6*dy)/eps)+1.0)-0.5*(tanh((x[1]+4*dy)/eps)+1.0));
		g[2]=0;
		g[3]=0;
		g[4]=0;
		g[5]=0;
		g[6]=0;
		g[7]=0;

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

//System for dPi/dt=curl(Pi)
	#undef  __FUNCT__
	#define __FUNCT__ "Pi"
	PetscErrorCode Pi(IGAPoint p,PetscReal *K,PetscReal *F, PetscReal *U,void *ctx)		//En este caso se necesita solo un IGAPoint, pues Pi(t) y Pi(t-1) son iguales en dimensión y ordenamiento
	{
		AppCtx    *user  = (AppCtx *)ctx;
		PetscReal dt     = user->dt;

		const PetscReal *N0,(*N1)[2];
		IGAPointGetShapeFuns(p,0,(const PetscReal**)&N0);
		IGAPointGetShapeFuns(p,1,(const PetscReal**)&N1);
		PetscInt a,b,i,nen=p->nen, dof=p->dof;

		PetscReal Pi0[dof];																	//Assign Pi(t) to a vector
		PetscReal d_Pi0[dof][2];															//Same for the partial derivatives
		IGAPointFormValue(p,U,&Pi0[0]);														
		IGAPointFormGrad (p,U,&d_Pi0[0][0]);

		PetscReal (*Kpi)[dof][nen][dof] = (typeof(Kpi)) K;
		PetscReal (*Fpi)[dof] = (PetscReal (*)[dof])F;

		PetscReal V1=1.0;
		PetscReal V2=0.0;
		PetscReal V1_x=0.0;
		PetscReal V2_y=0.0;

		PetscInt dim = p->dim;																//Dimension of the problem
		PetscReal x[dim];
		IGAPointFormGeomMap(p,x);

		if (p->atboundary)
		{
			PetscInt dir  = p->boundary_id / 2;
			PetscInt side = p->boundary_id % 2;

			if(dir==0 && side==1)					//REMEMBER to update this if to take into account the actual velocities from AppCtx in the final version
			{
				for (a=0; a<nen; a++)
				{
					PetscReal Na   = N0[a];
					for(b=0; b<nen; b++)
					{
						PetscReal Nb   = N0[b];
						for (i=0; i<dof; i++)
						{
							Kpi[a][i][b][i] = dt*Na*Nb*V1;
						}
					}
				}
			}
		}
		else
		{
			for (a=0; a<nen; a++) 
			{
				PetscReal Na   = N0[a];
				PetscReal Na_x = N1[a][0];
				PetscReal Na_y = N1[a][1];
				for (b=0; b<nen; b++) 
				{
					PetscReal Nb   = N0[b];
					PetscReal Nb_x = N1[b][0];
					PetscReal Nb_y = N1[b][1];

					for (i=0; i<dof; i++)
					{
						Kpi[a][i][b][i] = Na*Nb												//(LS1)
										+ dt*(Nb*(Na_x*V1+Na_y*V2))							//(LS2)
										+ dt*(Nb*Na*(V1_x+V2_y))							//(LS3)
										+ dt*(Na*(Nb_x*V1+Nb_y*V2))							//(LS6)
										+ dt*(Nb*Na*(V1_x+V2_y))							//(LS7)  equal to (LS3)
										+ dt*dt*((Nb_x*V1+Nb_y*V2)*(Na_x*V1+Na_x*V2))		//(LS10)
										+ dt*dt*((Nb_x*V1+Nb_y*V2)*Na*(V1_x+V2_y))			//(LS11)
										+ dt*dt*(Nb*(V1_x+V2_y)*(Na_x*V1+Na_y*V2))			//(LS14)
										+ dt*dt*(Nb*(V1_x+V2_y)*Na*(V1_x+V2_y))				//(LS15)
										+ Na*Nb												//(GL1) 
										- dt*(Nb*(Na_x*V1+Na_y*V2))							//(GL3)  equal to the negative of (LS2)
										;		//Los demás terminos son 0 en 2d
					}
				}

				for (i=0; i<dof; i++)
				{
					Fpi[a][i] = Pi0[i]*Na 													//(LS26)
								+dt*(Pi0[i]*(Na_x*V1+Na_y*V2))								//(LS27)
								+dt*(Pi0[i]*Na*(V1_x+V2_y))									//(LS28)
								+Pi0[i]*Na 													//(GL2)
								;
				}
			}
		}
		//PetscPrintf(PETSC_COMM_WORLD,"p->nen: %d \n", p->nen);
		
		return 0;
	}
//

//System for dAlpha/dt=curl(Alpha)
	#undef  __FUNCT__
	#define __FUNCT__ "Alfa"
	PetscErrorCode Alfa(IGAPoint p,PetscReal *K,PetscReal *F, PetscReal *U,void *ctx)
	{
		AppCtx    *user  = (AppCtx *)ctx;
		PetscReal dt     = user->dt;

		const PetscReal *N0,(*N1)[2];
		IGAPointGetShapeFuns(p,0,(const PetscReal**)&N0);
		IGAPointGetShapeFuns(p,1,(const PetscReal**)&N1);
		PetscInt a,b,i,nen=p->nen, dof=p->dof;

		PetscReal Al0[dof];																	//Assign Al(t) to a vector
		PetscReal d_Al0[dof][2];															//Same for partial derivatives
		IGAPointFormValue(p,U,&Al0[0]);
		IGAPointFormGrad (p,U,&d_Al0[0][0]);

		PetscReal (*Kal)[dof][nen][dof] = (typeof(Kal)) K;
		PetscReal (*Fal)[dof] = (PetscReal (*)[dof])F;

		PetscReal V1=0.0;
		PetscReal V2=1.0;
		PetscReal V1_x=0.0;
		PetscReal V2_y=0.0;

		PetscInt dim = p->dim;																//Dimension of the problem
		PetscReal x[dim];
		IGAPointFormGeomMap(p,x);

		if (p->atboundary)
		{
			PetscInt dir  = p->boundary_id / 2;
			PetscInt side = p->boundary_id % 2;

										//REMEMBER to update this if to take into account the actual velocities from AppCtx in the final version
			if(dir==1 && side==1)		//Dir=0 is the x axis, Dir=1 is the y axis. Side=0 is the "negative" side (left or down), Side=1 is the "positive" side (right or up)
			{
				for (a=0; a<nen; a++)
				{
					PetscReal Na   = N0[a];
					for(b=0; b<nen; b++)
					{
						PetscReal Nb   = N0[b];
						for (i=0; i<dof; i++)
						{
							Kal[a][i][b][i] = dt*Na*Nb*V2;										//(GL7)
						}
					}
				}
			}
		}
		else
		{
			for (a=0; a<nen; a++) 
			{
				PetscReal Na   = N0[a];
				PetscReal Na_x = N1[a][0];
				PetscReal Na_y = N1[a][1];
				for (b=0; b<nen; b++) 
				{
					PetscReal Nb   = N0[b];
					//PetscReal Nb_x = N1[b][0];
					//PetscReal Nb_y = N1[b][1];

					for (i=0; i<dof; i++)
					{
						Kal[a][i][b][i] =	Na*Nb												//(LS1)
										+dt*Nb*(Na_x*V1+Na_y*V2)								//(LS2)
										+dt*Nb*Na*(V1_x+V2_y)									//(LS3)
										+Na*Nb													//(GL1)
										-dt*Nb*(Na_x*V1+Na_y*V2)								//(GL3) (negative of (LS2))
										;		//Los demás terminos son 0 en 2d
					}
				}

				for (i=0; i<dof; i++)
				{
					Fal[a][i] = Al0[i]*Na														//(LS6)
								-dt*Na*(d_Al0[i][0]*V1+d_Al0[i][1]*V2)							//(LS7)
								-dt*Al0[i]*Na*(V1_x+V2_y)										//(LS8)
								+dt*Al0[i]*(Na_x*V1+Na_y*V2)									//(LS11)
								+dt*Al0[i]*Na*(V1_x+V2_y)										//(LS12) (negative of (LS8))
								-dt*dt*(d_Al0[i][0]*V1+d_Al0[i][1]*V2)*(Na_x*V1+Na_y*V2)		//(LS15)
								-dt*dt*Na*(d_Al0[i][0]*V1+d_Al0[i][1]*V2)*(V1_x+V2_y)			//(LS16)
								-dt*dt*Al0[i]*(V1_x+V2_y)*(Na_x*V1+Na_y*V2)						//(LS19)
								-dt*dt*Al0[i]*Na*(V1_x+V2_y)*(V1_x+V2_y)						//(LS20)
								+Al0[i]*Na														//(GL2)
								;
				}
			}
		}
		
		return 0;
	}
//

//System for curl(chi)=alpha & div(chi)=0
	#undef  __FUNCT__
	#define __FUNCT__ "Chi"
	PetscErrorCode Chi(IGAPoint p,IGAPoint pAl,PetscReal *K,PetscReal *F, PetscReal *U,void *ctx)
	{
		//Definition of alternating tensor
		const PetscReal e[3][3][3]=
		{
			{{0.0,0.0,0.0},{0.0,0.0,1.0},{0.0,-1.0,0.0}},
			{{0.0,0.0,-1.0},{0.0,0.0,0.0},{1.0,0.0,0.0}},
			{{0.0,1.0,0.0},{-1.0,0.0,0.0},{0.0,0.0,0.0}}
		};

		const PetscReal *N0,(*N1)[2];
		IGAPointGetShapeFuns(p,0,(const PetscReal**)&N0);
		IGAPointGetShapeFuns(p,1,(const PetscReal**)&N1);
		PetscInt a,b,u,w,i,j,k,l,nen=p->nen, dof=p->dof;

		PetscReal x[2];																			//Vector of reals, size equal to problem's dimension
		IGAPointFormGeomMap(p,x);																//Fills x with the coordinates of p, Gauss's point

		PetscReal Al0[2];																		//Create array to recieve Alfa
		//PetscReal dAl0[2][2];																	//Create array to recieve dAlfa
		IGAPointFormValue(pAl,U,&Al0[0]);														//This fills the values
		//IGAPointFormGrad (pAl,U,&dAl0[0][0]);													//This fills the values of the derivatives

		PetscReal alfa[3][3]={0};
		alfa[0][2]=Al0[0];	alfa[1][2]=Al0[1];													//Assign Al0 to the full alfa expression

		PetscReal (*Kchi)[dof][nen][dof] = (typeof(Kchi)) K;
		PetscReal (*Fchi)[dof] = (PetscReal (*)[dof])F;

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

					PetscReal v[3][3][3]={0};

					if(u==0){
						v[0][0][0]=Na_x; v[0][0][1]=Na_y;}
					if(u==1){
						v[0][1][0]=Na_x; v[0][1][1]=Na_y;}
					if(u==2){
						v[1][0][0]=Na_x; v[1][0][1]=Na_y;}
					if(u==3){
						v[1][1][0]=Na_x; v[1][1][1]=Na_y;}

					for (b=0; b<nen; b++)
					{
						PetscReal Nb_x = N1[b][0];
						PetscReal Nb_y = N1[b][1];
						for (w=0; w<dof; w++)
						{
							PetscReal chi[3][3][3]={0};

							if(w==0){
								chi[0][0][0]=Nb_x; chi[0][0][1]=Nb_y;}
							if(w==1){
								chi[0][1][0]=Nb_x; chi[0][1][1]=Nb_y;}
							if(w==2){
								chi[1][0][0]=Nb_x; chi[1][0][1]=Nb_y;}
							if(w==3){
								chi[1][1][0]=Nb_x; chi[1][1][1]=Nb_y;}

							for(i=0;i<3;i++)
							{
								for(j=0;j<3;j++)
								{
									for(k=0;k<3;k++)
									{
										Kchi[a][u][b][w]+=chi[i][j][k]*v[i][j][k]-chi[i][j][k]*v[i][k][j]+chi[i][j][j]*v[i][k][k];
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
								for(l=0;l<3;l++)
								{
									Fchi[a][u]+=alfa[i][j]*e[j][k][l]*v[i][l][k];
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

//System for div(grad(z\dot))=curl(alpha x V)
	#undef  __FUNCT__
	#define __FUNCT__ "divGradZ"
	PetscErrorCode divGradZ(IGAPoint p,IGAPoint pAl,PetscReal *K,PetscReal *F, PetscReal *U,void *ctx)
	{
		//AppCtx    *user  = (AppCtx *)ctx;
		//PetscReal dt     = user->dt;

		//Definition of alternating tensor
		const PetscReal e[3][3][3]=
		{
			{{0.0,0.0,0.0},{0.0,0.0,1.0},{0.0,-1.0,0.0}},
			{{0.0,0.0,-1.0},{0.0,0.0,0.0},{1.0,0.0,0.0}},
			{{0.0,1.0,0.0},{-1.0,0.0,0.0},{0.0,0.0,0.0}}
		};
		
		PetscReal alfa[3][3]={0};
		PetscReal dAlfa[3][3][3]={0};
		PetscReal N[3]={0};

		const PetscReal *N0,(*N1)[2];
		IGAPointGetShapeFuns(p,0,(const PetscReal**)&N0);
		IGAPointGetShapeFuns(p,1,(const PetscReal**)&N1);
		PetscInt a,b,i,j,k,l,m,n,nen=p->nen, dof=p->dof;

		PetscReal Al0[2];																		//Create array to recieve Alfa
		PetscReal dAl0[dof][2];																	//Create array to recieve dAlfa
		IGAPointFormValue(pAl,U,&Al0[0]);														//This fills the values
		IGAPointFormGrad (pAl,U,&dAl0[0][0]);													//This fills the values of the derivatives

		PetscReal x[2];																			//Vector of reals, size equal to problem's dimension
		IGAPointFormGeomMap(p,x);																//Fills x with the coordinates of p, Gauss's point
		PetscReal V[3]={0};		V[0]=x[0]; V[1]=x[1];
		PetscReal dV[3][3]={0}; dV[0][0]=1.0; dV[0][1]=1.0;

		alfa[0][2]=Al0[0];	alfa[1][2]=Al0[1];													//Assign Al0 to the full alfa expression
		dAlfa[0][2][0]=dAl0[0][0]; dAlfa[0][2][1]=dAl0[0][1];
		dAlfa[1][2][0]=dAl0[1][0]; dAlfa[1][2][1]=dAl0[1][1];

		PetscReal (*Kz)[dof][nen][dof] = (typeof(Kz)) K;
		PetscReal (*Fz)[dof] = (PetscReal (*)[dof])F;

		if (p->atboundary)
		{
			PetscInt dir  = p->boundary_id / 2;
			PetscInt side = p->boundary_id % 2;

			if(dir==0)
			{
				if(side==0)
					N[0]=-1.0;	
				else
					N[0]=1.0;
			}
			else
			{
				if(side==0)
					N[1]=-1.0;	
				else
					N[1]=1.0;
			}

			for (a=0;a<nen;a++)
			{
				PetscReal Na   = N0[a];
				for(i=0;i<dof;i++)
				{
					PetscReal v[3]={0};
					v[i]=Na;

					for(k=0;k<3;k++)
					{
						for(l=0;l<3;l++)
						{
							for(m=0;m<3;m++)
							{
								for(n=0;n<3;n++)
								{
									Fz[a][i]+=e[l][n][m]*alfa[k][n]*V[m]*N[l]*v[k];
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
				PetscReal Na   = N0[a];
				PetscReal Na_x = N1[a][0];
				PetscReal Na_y = N1[a][1];

				for (i=0; i<dof; i++)
				{
					PetscReal v[3]={0};
					v[i]=Na;

					PetscReal dv[3][3]={0};
					dv[i][0]=Na_x; dv[i][1]=Na_y;

					for (b=0; b<nen; b++)
					{
						PetscReal Nb_x = N1[b][0];
						PetscReal Nb_y = N1[b][1];

						for(j=0; j<dof; j++)
						{
							PetscReal dz[3][3]={0};
							dz[j][0]=Nb_x; dz[j][1]=Nb_y;

							for(k=0;k<3;k++)
							{	
								for(l=0;l<3;l++)
								{
									Kz[a][i][b][j]+=dz[k][l]*dv[k][l];
								}
							}
						}
					}

					for(k=0;k<3;k++)
					{
						for(l=0;l<3;l++)
						{
							for(m=0;m<3;m++)
							{
								for(n=0;n<3;n++)
								{
									Fz[a][i]+=e[l][n][m]*dAlfa[k][n][l]*V[m]*v[k]+e[l][n][m]*alfa[k][n]*dV[m][l]*v[k];
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

//System for Pi=curl(S)
	#undef  __FUNCT__
	#define __FUNCT__ "curlSPi"
	PetscErrorCode curlSPi(IGAPoint p,IGAPoint pS,PetscReal *K,PetscReal *F, PetscReal *U,void *ctx)
	{
		//Definition of alternating tensor
		const PetscReal e[3][3][3]=
		{
			{{0.0,0.0,0.0},{0.0,0.0,1.0},{0.0,-1.0,0.0}},
			{{0.0,0.0,-1.0},{0.0,0.0,0.0},{1.0,0.0,0.0}},
			{{0.0,1.0,0.0},{-1.0,0.0,0.0},{0.0,0.0,0.0}}
		};

		//PetscReal S[3][3][3]={0};
		PetscReal dS[3][3][3][3]={0};

		const PetscReal *N0,(*N1)[2];
		IGAPointGetShapeFuns(p,0,(const PetscReal**)&N0);
		IGAPointGetShapeFuns(p,1,(const PetscReal**)&N1);
		PetscInt a,b,i,j,k,nen=p->nen, dof=p->dof;

		//PetscReal S0[8];																//Assign Pi(t) to a vector
		PetscReal  dS0[8][2];															//Same for the partial derivatives
		//IGAPointFormValue(pS,U,&S0[0]);													//Fills the vector with the values from outside
		IGAPointFormGrad (pS,U,&dS0[0][0]);												//Fills the matrix with the partial derivatives

		dS[0][0][0][0]=dS0[0][0];		//01	dxS(1,1,1)
		dS[0][0][0][1]=dS0[0][1];		//02	dyS(1,1,1)
		dS[0][1][0][0]=dS0[1][0];		//03	dxS(1,2,1)
		dS[0][1][0][1]=dS0[1][1];		//04	dyS(1,2,1)
		dS[1][0][0][0]=dS0[2][0];		//05	dxS(2,1,1)
		dS[1][0][0][1]=dS0[2][1];		//06	dyS(2,1,1)
		dS[1][1][0][0]=dS0[3][0];		//07	dxS(2,2,1)
		dS[1][1][0][1]=dS0[3][1];		//08	dyS(2,2,1)
		dS[0][0][1][0]=dS0[4][0];		//09	dxS(1,1,2)
		dS[0][0][1][1]=dS0[4][1];		//10	etc.
		dS[0][1][1][0]=dS0[5][0];		//11
		dS[0][1][1][1]=dS0[5][1];		//12
		dS[1][0][1][0]=dS0[6][0];		//13
		dS[1][0][1][1]=dS0[6][1];		//14
		dS[1][1][1][0]=dS0[7][0];		//15
		dS[1][1][1][1]=dS0[7][1];		//16

		PetscReal curlS[3][3][3]={0};
		for(i=0; i<3; i++)
		{
			for(j=0; j<3; j++)
			{
				for(k=0; k<3; k++)
				{
					for (a=0; a<2; a++)
					{
						for(b=0; b<3; b++)
						{
							curlS[i][j][k]+=e[k][a][b]*dS[i][j][b][a];
						}
					}
				}
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
				if (i==0)
				{
					FF[a][i] = N[a]*curlS[0][0][2];
				}
				else if (i==1)
				{
					FF[a][i] = N[a]*curlS[0][1][2];
				}
				else if (i==2)
				{
					FF[a][i] = N[a]*curlS[1][0][2];
				}
				else
				{
					FF[a][i] = N[a]*curlS[1][1][2];
				}
			}
		}
		return 0;
	}
//

//System for equilibrium
	#undef  __FUNCT__
	#define __FUNCT__ "Eq"
	PetscErrorCode Eq(IGAPoint p,IGAPoint pZ,IGAPoint pChi,PetscReal *K,PetscReal *F,PetscReal *UZ, PetscReal *UChi,void *ctx)		//When other fields like Pi or S (etc.) come into this, declare a IGAPoint pPi or pS and a new PescReal *UPi or *US for each
	{
		PetscReal C[3][3][3][3]={0};														//Initialization of elastic tensor
		const PetscReal *N0,(*N1)[2],(*N2)[2][2];
		IGAPointGetShapeFuns(p,0,(const PetscReal**)&N0);									//Value of the shape functions
		IGAPointGetShapeFuns(p,1,(const PetscReal**)&N1);									//Derivatives of the shape functions
		IGAPointGetShapeFuns(p,2,(const PetscReal**)&N2);									//Second derivatives of the shape funcions
		//After this command Na_xx=N2[0][0], Na_yy=N2[1][1], Na_xy=N2[0][1], Na_yx[1][0] (these last two are equal)
		PetscInt a,b,i,j,k,l,u,w,r,s,m,nen=p->nen, dof=p->dof;

		PetscReal Z0[dof];																	//Assign Z(0) to a vector
		PetscReal d_Z0[dof][2];																//Same for its partial derivatives
		IGAPointFormValue(pZ,UZ,&Z0[0]);
		IGAPointFormGrad (pZ,UZ,&d_Z0[0][0]);

		PetscReal (*Keq)[dof][nen][dof] = (typeof(Keq)) K;
		PetscReal (*Feq)[dof] = (PetscReal (*)[dof])F;

		PetscInt dim = p->dim;																//Dimension of the problem
		PetscReal x[dim];
		IGAPointFormGeomMap(p,x);

		//E and nu should come from AppCtx in the future
		const PetscReal E=2100000.0*9.81*10000.0;
		const PetscReal nu=0.33;
		const PetscReal lambda=(E*nu)/((1.0+nu)*(1.0-2.0*nu));
		const PetscReal mu=E/(2.0*(1.0+nu));
		const PetscReal eps=0.0*1.0/1000.0;							//Choose later based on whatever Amit says :)

		//////////////Delete this parte later, loads should come from appCtx or be 0
		PetscReal f[2]={0.0, 0.0};			//Distributed load in body
		PetscReal g[2]={0.0, -E};			//Boundary load (applied wherever is defined in the p->atboundary block)
		PetscReal M[3]={0.0, 0.0, 0.0};		//Distributed load in body (only 3rd component applies in 2d)				
		PetscReal h[3]={0.0, 0.0, 0.0};		//Boundary moment (applied wherever is defined in the p->atboundary block) (only 3rd component applies in 2d)
		/////////////////////////////////////
		const PetscReal e[3][3][3]=
		{
			{{0.0,0.0,0.0},{0.0,0.0,1.0},{0.0,-1.0,0.0}},
			{{0.0,0.0,-1.0},{0.0,0.0,0.0},{1.0,0.0,0.0}},
			{{0.0,1.0,0.0},{-1.0,0.0,0.0},{0.0,0.0,0.0}}
		};

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

		PetscReal v[2]={0};
		PetscReal eu[3][3]={0};
		PetscReal ev[3][3]={0};
		//PetscReal deu[3][3][3]={0};
		PetscReal dev[3][3][3]={0};
		PetscReal du[3][3][3]={0};
		PetscReal dv[3][3][3]={0};
		PetscReal dvp[3][3]={0};

		if (p->atboundary)
		{
			if(x[0]>-1.0/100.0 && x[0]<1.0/100.0)
			{
				PetscInt dir  = p->boundary_id / 2;
				PetscInt side = p->boundary_id % 2;

				for (a=0; a<nen; a++)
				{
					PetscReal Na   = N0[a];
					PetscReal Na_x = N1[a][0];
					PetscReal Na_y = N1[a][1];

					for (i=0; i<dof; i++)
					{
						if (i==0)
						{
							dvp[0][0]=Na_x; dvp[0][1]=Na_y; dvp[0][2]=0.0;
							dvp[1][0]=0.0;  dvp[1][1]=0.0;  dvp[1][2]=0.0;
							dvp[2][0]=0.0;  dvp[2][1]=0.0;  dvp[2][2]=0.0;
						}
						else if (i==1)
						{
							dvp[0][0]=0.0;  dvp[0][1]=0.0;  dvp[0][2]=0.0;
							dvp[1][0]=Na_x; dvp[1][1]=Na_y; dvp[1][2]=0.0;
							dvp[2][0]=0.0;  dvp[2][1]=0.0;  dvp[2][2]=0.0;
						}
						if (dir==1)
						{
							if(side==1)
							{
								Feq[a][i]+=Na*g[i];

								for (m=0; m<3; m++)
								{
									for (s=0; s<3; s++)
									{
										for (r=0; r<3;r++)
										{
											Feq[a][i]+=0.5*h[m]*e[m][r][s]*dvp[s][r];
										}
									}
								}
							}
						}
					}
				}
			}
		}
		else
		{
			for (a=0; a<nen; a++) 
			{
				PetscReal Na_x = N1[a][0];
				PetscReal Na_y = N1[a][1];
				PetscReal Na_xx =N2[a][0][0];
				PetscReal Na_yy =N2[a][1][1];
				PetscReal Na_xy =N2[a][0][1];
				PetscReal Na_yx =N2[a][1][0];

				for (b=0; b<nen; b++) 
				{
					PetscReal Nb_x = N1[b][0];
					PetscReal Nb_y = N1[b][1];
					PetscReal Nb_xx =N2[b][0][0];
					PetscReal Nb_yy =N2[b][1][1];
					PetscReal Nb_xy =N2[b][0][1];
					PetscReal Nb_yx =N2[b][1][0];

					for (i=0; i<dof; i++)
					{
						if (i==0)
						{
							ev[0][0]=Na_x; 	ev[0][1]=0.5*Na_y; ev[0][2]=0.0; ev[1][0]=0.5*Na_y; ev[1][1]=0.0; ev[1][2]=0.0;	ev[2][0]=0.0; ev[2][1]=0.0; ev[2][2]=0.0;  
							dev[0][0][0]=Na_xx; dev[0][0][1]=Na_xy; dev[0][0][2]=0; dev[0][1][0]=0.5*Na_yx; dev[0][1][1]=0.5*Na_yy; dev[0][1][2]=0; dev[0][2][0]=0; dev[0][2][1]=0; dev[0][2][2]=0;
							dev[1][0][0]=0.5*Na_yx; dev[1][0][1]=0.5*Na_yy; dev[1][0][2]=0; dev[1][1][0]=0; dev[1][1][1]=0; dev[1][1][2]=0; dev[1][2][0]=0; dev[1][2][1]=0; dev[1][2][2]=0;
							dev[2][0][0]=0; dev[2][0][1]=0; dev[2][0][2]=0; dev[2][1][0]=0; dev[2][1][1]=0; dev[2][1][2]=0; dev[2][2][0]=0; dev[2][2][1]=0; dev[2][2][2]=0;
							
							dv[0][0][0]=Na_xx; dv[0][0][1]=Na_xy; dv[0][0][2]=0.0; dv[0][1][0]=Na_yx; dv[0][1][1]=Na_yy; dv[0][1][2]=0.0; dv[0][2][0]=0.0; dv[0][2][1]=0.0; dv[0][2][2]=0.0;
							dv[1][0][0]=0.00;  dv[1][0][1]=0.00;  dv[1][0][2]=0.0; dv[1][1][0]=0.0;   dv[1][1][1]=0.0;   dv[1][1][2]=0.0; dv[1][2][0]=0.0; dv[1][2][1]=0.0; dv[1][2][2]=0.0;
							dv[2][0][0]=0.00;  dv[2][0][1]=0.00;  dv[2][0][2]=0.0; dv[2][1][0]=0.0;   dv[2][1][1]=0.0;   dv[2][1][2]=0.0; dv[2][2][0]=0.0; dv[2][2][1]=0.0; dv[2][2][2]=0.0;
							
						}
						else if(i==1)
						{
							ev[0][0]=0.0; ev[0][1]=0.5*Na_x; ev[0][2]=0.0; ev[1][0]=0.5*Na_x; ev[1][1]=Na_y; ev[1][2]=0.0; ev[2][0]=0.0; ev[2][1]=0.0; ev[2][2]=0.0;  
							dev[0][0][0]=0.0; dev[0][0][1]=0.0; dev[0][0][2]=0; dev[0][1][0]=0.5*Na_xx; dev[0][1][1]=0.5*Na_xy; dev[0][1][2]=0; dev[0][2][0]=0; dev[0][2][1]=0; dev[0][2][2]=0;
							dev[1][0][0]=0.5*Na_xx; dev[1][0][1]=0.5*Na_xy; dev[1][0][2]=0; dev[1][1][0]=Na_yx; dev[1][1][1]=Na_yy; dev[1][1][2]=0; dev[1][2][0]=0; dev[1][2][1]=0; dev[1][2][2]=0;
							dev[2][0][0]=0; dev[2][0][1]=0; dev[2][0][2]=0; dev[2][1][0]=0; dev[2][1][1]=0; dev[2][1][2]=0; dev[2][2][0]=0; dev[2][2][1]=0; dev[2][2][2]=0;

							dv[0][0][0]=0.00;  dv[0][0][1]=0.00;  dv[0][0][2]=0.0; dv[0][1][0]=0.0;   dv[0][1][1]=0.0;   dv[0][1][2]=0.0; dv[0][2][0]=0.0; dv[0][2][1]=0.0; dv[0][2][2]=0.0;
							dv[1][0][0]=Na_xx; dv[1][0][1]=Na_xy; dv[1][0][2]=0.0; dv[1][1][0]=Na_yx; dv[1][1][1]=Na_yy; dv[1][1][2]=0.0; dv[1][2][0]=0.0; dv[1][2][1]=0.0; dv[1][2][2]=0.0;
							dv[2][0][0]=0.00;  dv[2][0][1]=0.00;  dv[2][0][2]=0.0; dv[2][1][0]=0.0;   dv[2][1][1]=0.0;   dv[2][1][2]=0.0; dv[2][2][0]=0.0; dv[2][2][1]=0.0; dv[2][2][2]=0.0;
							
						}

						for (j=0; j<dof; j++)
						{
							if (j==0)
							{
								eu[0][0]=Nb_x; eu[0][1]=0.5*Nb_y; eu[1][0]=0.5*Nb_y; eu[1][1]=0.0; eu[2][0]=0.0; eu[2][1]=0.0; eu[2][2]=0.0; eu[0][2]=0.0; eu[1][2]=0.0;
								du[0][0][0]=Nb_xx; du[0][0][1]=Nb_xy; du[0][0][2]=0.0; du[0][1][0]=Nb_yx; du[0][1][1]=Nb_yy; du[0][1][2]=0.0; du[0][2][0]=0.0; du[0][2][1]=0.0; du[0][2][2]=0.0;
								du[1][0][0]=0.00;  du[1][0][1]=0.00;  du[1][0][2]=0.0; du[1][1][0]=0.0;   du[1][1][1]=0.0;   du[1][1][2]=0.0; du[1][2][0]=0.0; du[1][2][1]=0.0; du[1][2][2]=0.0;
								du[2][0][0]=0.00;  du[2][0][1]=0.00;  du[2][0][2]=0.0; du[2][1][0]=0.0;   du[2][1][1]=0.0;   du[2][1][2]=0.0; du[2][2][0]=0.0; du[2][2][1]=0.0; du[2][2][2]=0.0;
							}
							else if(j==1)
							{
								eu[0][0]=0.0; eu[0][1]=0.5*Nb_x; eu[0][2]=0.0; eu[1][0]=0.5*Nb_x; eu[1][1]=Nb_y; eu[1][2]=0.0; eu[2][0]=0.0; eu[2][1]=0.0; eu[2][2]=0.0;
								du[0][0][0]=0.00;  du[0][0][1]=0.00;  du[0][0][2]=0.0; du[0][1][0]=0.0;   du[0][1][1]=0.0;   du[0][1][2]=0.0; du[0][2][0]=0.0; du[0][2][1]=0.0; du[0][2][2]=0.0;
								du[1][0][0]=Nb_xx; du[1][0][1]=Nb_xy; du[1][0][2]=0.0; du[1][1][0]=Nb_yx; du[1][1][1]=Nb_yy; du[1][1][2]=0.0; du[1][2][0]=0.0; du[1][2][1]=0.0; du[1][2][2]=0.0;
								du[2][0][0]=0.00;  du[2][0][1]=0.00;  du[2][0][2]=0.0; du[2][1][0]=0.0;   du[2][1][1]=0.0;   du[2][1][2]=0.0; du[2][2][0]=0.0; du[2][2][1]=0.0; du[2][2][2]=0.0;
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
											Keq[a][i][b][j]+=C[k][l][u][w]*eu[u][w]*ev[k][l];		//This is C_ijkl*grad(u)^sym*grad(v)^sym, which is wrong if (alpha, S and pi are not 0)
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
										Keq[a][i][b][j]+=0.5*eps*(du[r][s][m]+du[s][r][m])*dev[r][s][m]+0.5*eps*(du[s][r][m]*dv[s][r][m]-du[r][s][m]*dv[s][r][m]);
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
					PetscReal Na_x=N1[a][0];
					PetscReal Na_y=N1[a][1];

					if (i==0)
					{
						dvp[0][0]=Na_x; dvp[0][1]=Na_y; dvp[0][2]=0.0;
						dvp[1][0]=0.0;  dvp[1][1]=0.0;  dvp[1][2]=0.0;
						dvp[2][0]=0.0;  dvp[2][1]=0.0;  dvp[2][2]=0.0;
					}
					else if (i==1)
					{
						dvp[0][0]=0.0;  dvp[0][1]=0.0;  dvp[0][2]=0.0;
						dvp[1][0]=Na_x; dvp[1][1]=Na_y; dvp[1][2]=0.0;
						dvp[2][0]=0.0;  dvp[2][1]=0.0;  dvp[2][2]=0.0;
					}

					Feq[a][i] = v[i]*f[i];

					for (m=0; m<3; m++)
					{
						for (r=0; r<3; r++)
						{
							for (s=0; s<3; s++)
							{
								Feq[a][i]+=0.5*M[m]*e[m][r][s]*dvp[s][r];
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
	PetscErrorCode Z0sys(IGAPoint p,IGAPoint pChi,PetscReal *K,PetscReal *F, PetscReal *UChi,void *ctx)
	{
		PetscReal C[3][3][3][3]={0};														//Initialization of elastic tensor
		const PetscReal *N0,(*N1)[2],(*N2)[2][2];
		IGAPointGetShapeFuns(p,0,(const PetscReal**)&N0);									//Value of the shape functions
		IGAPointGetShapeFuns(p,1,(const PetscReal**)&N1);									//Derivatives of the shape functions
		IGAPointGetShapeFuns(p,2,(const PetscReal**)&N2);									//Second derivatives of the shape funcions
		//After this command Na_xx=N2[a][0][0], Na_yy=N2[a][1][1], Na_xy=N2[a][0][1], Na_yx[a][1][0] (these last two are equal, or should be at least) (remember a is the index of the shape function)
		PetscInt a,b,i,j,k,l,u,w,d,r,s,m,nen=p->nen, dof=p->dof, dofChi=pChi->dof;

		PetscReal Chi0[dofChi];																//Assign Al(t) to a vector
		PetscReal d_Chi0[dofChi][2];															//Same for partial derivatives
		IGAPointFormValue(pChi,UChi,&Chi0[0]);
		IGAPointFormGrad (pChi,UChi,&d_Chi0[0][0]);

		PetscReal fullChi[3][3];
		fullChi[0][0]=Chi0[0]; 	fullChi[0][1]=Chi0[1]; 	fullChi[0][2]=0;
		fullChi[1][0]=Chi0[2]; 	fullChi[1][1]=Chi0[3]; 	fullChi[1][2]=0;
		fullChi[2][0]=0; 		fullChi[2][1]=0; 		fullChi[2][2]=0;

		PetscReal fulld_Chi[3][3][3]={0};
		
		fulld_Chi[0][0][0]=d_Chi0[0][0]; fulld_Chi[0][0][1]=d_Chi0[0][1];
		fulld_Chi[0][1][0]=d_Chi0[1][0]; fulld_Chi[0][1][1]=d_Chi0[1][1];
		fulld_Chi[1][0][0]=d_Chi0[2][0]; fulld_Chi[1][0][1]=d_Chi0[2][1];
		fulld_Chi[1][1][0]=d_Chi0[3][0]; fulld_Chi[1][1][1]=d_Chi0[3][1];

		PetscReal (*Keq)[dof][nen][dof] = (typeof(Keq)) K;
		PetscReal (*Feq)[dof] = (PetscReal (*)[dof])F;

		PetscInt dim = p->dim;																//Dimension of the problem
		PetscReal x[dim];
		IGAPointFormGeomMap(p,x);

		//E and nu should come from AppCtx in the future
		const PetscReal E=2100000.0*9.81*10000.0;
		const PetscReal nu=0.33;
		const PetscReal lambda=(E*nu)/((1.0+nu)*(1.0-2.0*nu));
		const PetscReal mu=E/(2.0*(1.0+nu));
		const PetscReal eps=E/1000.0;							//Choose later based on whatever Amit says :)

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

		PetscReal ez[3][3]={0};
		PetscReal ev[3][3]={0};
		//PetscReal dez[3][3][3]={0};
		PetscReal dev[3][3][3]={0};
		PetscReal dz[3][3][3]={0};
		PetscReal dv[3][3][3]={0};

		//Consider deleting the "atboundary" block for speed, no boundary loads on this problem
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
				PetscReal Na_xx =N2[a][0][0];
				PetscReal Na_yy =N2[a][1][1];
				PetscReal Na_xy =N2[a][0][1];
				PetscReal Na_yx =N2[a][1][0];

				for (b=0; b<nen; b++) 
				{
					PetscReal Nb_x = N1[b][0];
					PetscReal Nb_y = N1[b][1];
					PetscReal Nb_xx =N2[b][0][0];
					PetscReal Nb_yy =N2[b][1][1];
					PetscReal Nb_xy =N2[b][0][1];
					PetscReal Nb_yx =N2[b][1][0];

					for (i=0; i<dof; i++)
					{
						if (i==0)
						{
							ev[0][0]=Na_x; 	ev[0][1]=0.5*Na_y; ev[1][0]=0.5*Na_y; ev[1][1]=0.0; ev[2][0]=0.0; ev[2][1]=0.0; ev[2][2]=0.0; ev[0][2]=0.0; ev[1][2]=0.0;
							dev[0][0][0]=Na_xx; dev[0][1][0]=0.5*Na_yx; dev[1][0][0]=0.5*Na_yx; dev[1][1][0]=0; dev[2][0][0]=0; dev[2][1][0]=0; dev[2][2][0]=0; dev[0][2][0]=0; dev[1][2][0]=0;
							dev[0][0][1]=Na_xy; dev[0][1][1]=0.5*Na_yy; dev[1][0][1]=0.5*Na_yy; dev[1][1][1]=0; dev[2][0][1]=0; dev[2][1][1]=0; dev[2][2][1]=0; dev[0][2][1]=0; dev[1][2][1]=0;
							dev[0][0][2]=0; 	dev[0][1][2]=0; 		dev[1][0][2]=0; 		dev[1][1][2]=0; dev[2][0][2]=0; dev[2][1][2]=0; dev[2][2][2]=0; dev[0][2][2]=0; dev[1][2][2]=0;
							dv[0][0][0]=Na_xx; dv[0][0][1]=Na_xy; dv[0][0][2]=0.0; dv[0][1][0]=Na_yx; dv[0][1][1]=Na_yy; dv[0][1][2]=0.0; dv[0][2][0]=0.0; dv[0][2][1]=0.0; dv[0][2][2]=0.0;
							dv[1][0][0]=0.00;  dv[1][0][1]=0.00;  dv[1][0][2]=0.0; dv[1][1][0]=0.0;   dv[1][1][1]=0.0;   dv[1][1][2]=0.0; dv[1][2][0]=0.0; dv[1][2][1]=0.0; dv[1][2][2]=0.0;
							dv[2][0][0]=0.00;  dv[2][0][1]=0.00;  dv[2][0][2]=0.0; dv[2][1][0]=0.0;   dv[2][1][1]=0.0;   dv[2][1][2]=0.0; dv[2][2][0]=0.0; dv[2][2][1]=0.0; dv[2][2][2]=0.0;
							
						}
						else if(i==1)
						{
							ev[0][0]=0.0; ev[0][1]=0.5*Na_x; ev[1][0]=0.5*Na_x; ev[1][1]=Na_y; ev[2][0]=0.0; ev[2][1]=0.0; ev[2][2]=0.0; ev[0][2]=0.0; ev[1][2]=0.0;
							dev[0][0][0]=0; dev[0][1][0]=0.5*Na_xx; dev[1][0][0]=0.5*Na_xx; dev[1][1][0]=Na_yx; dev[2][0][0]=0; dev[2][1][0]=0; dev[2][2][0]=0; dev[0][2][0]=0; dev[1][2][0]=0;
							dev[0][0][1]=0; dev[0][1][1]=0.5*Na_xy; dev[1][0][1]=0.5*Na_xy; dev[1][1][1]=Na_yy; dev[2][0][1]=0; dev[2][1][1]=0; dev[2][2][1]=0; dev[0][2][1]=0; dev[1][2][1]=0;
							dev[0][0][2]=0; dev[0][1][2]=0; 		dev[1][0][2]=0; 		dev[1][1][2]=0; 	dev[2][0][2]=0; dev[2][1][2]=0; dev[2][2][2]=0; dev[0][2][2]=0; dev[1][2][2]=0;
							dv[0][0][0]=0.00;  dv[0][0][1]=0.00;  dv[0][0][2]=0.0; dv[0][1][0]=0.0;   dv[0][1][1]=0.0;   dv[0][1][2]=0.0; dv[0][2][0]=0.0; dv[0][2][1]=0.0; dv[0][2][2]=0.0;
							dv[1][0][0]=Na_xx; dv[1][0][1]=Na_xy; dv[1][0][2]=0.0; dv[1][1][0]=Na_yx; dv[1][1][1]=Na_yy; dv[1][1][2]=0.0; dv[1][2][0]=0.0; dv[1][2][1]=0.0; dv[1][2][2]=0.0;
							dv[2][0][0]=0.00;  dv[2][0][1]=0.00;  dv[2][0][2]=0.0; dv[2][1][0]=0.0;   dv[2][1][1]=0.0;   dv[2][1][2]=0.0; dv[2][2][0]=0.0; dv[2][2][1]=0.0; dv[2][2][2]=0.0;
							
						}

						for (j=0; j<dof; j++)
						{
							if (j==0)
							{
								ez[0][0]=Nb_x; ez[0][1]=0.5*Nb_y; ez[1][0]=0.5*Nb_y; ez[1][1]=0.0; ez[2][0]=0.0; ez[2][1]=0.0; ez[2][2]=0.0; ez[0][2]=0.0; ez[1][2]=0.0;
								//dez[0][0][0]=Nb_xx; dez[0][1][0]=0.5*Nb_yx; dez[1][0][0]=0.5*Nb_yx; dez[1][1][0]=0; dez[2][0][0]=0; dez[2][1][0]=0; dez[2][2][0]=0; dez[0][2][0]=0; dez[1][2][0]=0;
								//dez[0][0][1]=Nb_xy; dez[0][1][1]=0.5*Nb_yy; dez[1][0][1]=0.5*Nb_yy; dez[1][1][1]=0; dez[2][0][1]=0; dez[2][1][1]=0; dez[2][2][1]=0; dez[0][2][1]=0; dez[1][2][1]=0;
								//dez[0][0][2]=0; 	dez[0][1][2]=0; 		dez[1][0][2]=0; 		dez[1][1][2]=0; dez[2][0][2]=0; dez[2][1][2]=0; dez[2][2][2]=0; dez[0][2][2]=0; dez[1][2][0]=0;
								dz[0][0][0]=Nb_xx; dz[0][0][1]=Nb_xy; dz[0][0][2]=0.0; dz[0][1][0]=Nb_yx; dz[0][1][1]=Nb_yy; dz[0][1][2]=0.0; dz[0][2][0]=0.0; dz[0][2][1]=0.0; dz[0][2][2]=0.0;
								dz[1][0][0]=0.00;  dz[1][0][1]=0.00;  dz[1][0][2]=0.0; dz[1][1][0]=0.0;   dz[1][1][1]=0.0;   dz[1][1][2]=0.0; dz[1][2][0]=0.0; dz[1][2][1]=0.0; dz[1][2][2]=0.0;
								dz[2][0][0]=0.00;  dz[2][0][1]=0.00;  dz[2][0][2]=0.0; dz[2][1][0]=0.0;   dz[2][1][1]=0.0;   dz[2][1][2]=0.0; dz[2][2][0]=0.0; dz[2][2][1]=0.0; dz[2][2][2]=0.0;
							}
							else if(j==1)
							{
								ez[0][0]=0.0; ez[0][1]=0.5*Nb_x; ez[1][0]=0.5*Nb_x; ez[1][1]=Nb_y; ez[2][0]=0.0; ez[2][1]=0.0; ez[2][2]=0.0; ez[0][2]=0.0; ez[1][2]=0.0;
								//dez[0][0][0]=0; dez[0][1][0]=0.5*Nb_xx; dez[1][0][0]=0.5*Nb_xx; dez[1][1][0]=Nb_yx; dez[2][0][0]=0; dez[2][1][0]=0; dez[2][2][0]=0; dez[0][2][0]=0; dez[1][2][0]=0;
								//dez[0][0][1]=0; dez[0][1][1]=0.5*Nb_xy; dez[1][0][1]=0.5*Nb_xy; dez[1][1][1]=Nb_yy; dez[2][0][1]=0; dez[2][1][1]=0; dez[2][2][1]=0; dez[0][2][1]=0; dez[1][2][1]=0;
								//dez[0][0][2]=0; dez[0][1][2]=0; 		dez[1][0][2]=0; 		dez[1][1][2]=0; 	dez[2][0][2]=0; dez[2][1][2]=0; dez[2][2][2]=0; dez[0][2][2]=0; dez[1][2][2]=0;
								dz[0][0][0]=0.00;  dz[0][0][1]=0.00;  dz[0][0][2]=0.0; dz[0][1][0]=0.0;   dz[0][1][1]=0.0;   dz[0][1][2]=0.0; dz[0][2][0]=0.0; dz[0][2][1]=0.0; dz[0][2][2]=0.0;
								dz[1][0][0]=Nb_xx; dz[1][0][1]=Nb_xy; dz[1][0][2]=0.0; dz[1][1][0]=Nb_yx; dz[1][1][1]=Nb_yy; dz[1][1][2]=0.0; dz[1][2][0]=0.0; dz[1][2][1]=0.0; dz[1][2][2]=0.0;
								dz[2][0][0]=0.00;  dz[2][0][1]=0.00;  dz[2][0][2]=0.0; dz[2][1][0]=0.0;   dz[2][1][1]=0.0;   dz[2][1][2]=0.0; dz[2][2][0]=0.0; dz[2][2][1]=0.0; dz[2][2][2]=0.0;
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
											//Keq[a][i][b][j]+=C[k][l][u][w]*ez[u][w]*ev[k][l];		//This is C_ijkl*grad(u)^sym*grad(v)^sym, which is wrong
											Keq[a][i][b][j]+=-C[k][l][u][w]*ez[u][w]*ev[k][l];										
										}
									}
								}
							}
							//This keeps the run-time at 1:20, it's just expanding the product e[d][c][n]*e[d][r][s]
							for (r=0; r<3; r++)
							{
								for (s=0; s<3; s++)
								{
									for (m=0; m<3; m++)
									{
										Keq[a][i][b][j]+=-0.5*eps*(dz[r][s][m]+dz[s][r][m])*dev[r][s][m]-0.5*eps*(dz[s][r][m]-dz[r][s][m])*dv[s][r][m];
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
					//v[0]=0; v[1]=0;
					//v[i]=N0[a];
					PetscReal Na_x=N1[a][0];
					PetscReal Na_y=N1[a][1];
					PetscReal Na_xx =N2[a][0][0];
					PetscReal Na_yy =N2[a][1][1];
					PetscReal Na_xy =N2[a][0][1];
					PetscReal Na_yx =N2[a][1][0];

					if (i==0)
					{
						//dvp[0][0]=Na_x; dvp[0][1]=Na_y; dvp[0][2]=0.0;
						//dvp[1][0]=0.0;  dvp[1][1]=0.0;  dvp[1][2]=0.0;
						//dvp[2][0]=0.0;  dvp[2][1]=0.0;  dvp[2][2]=0.0;
						ev[0][0]=Na_x; ev[0][1]=0.5*Na_y; ev[1][0]=0.5*Na_y; ev[1][1]=0.0; ev[2][0]=0.0; ev[2][1]=0.0; ev[2][2]=0.0; ev[0][2]=0.0; ev[1][2]=0.0;
						dev[0][0][0]=Na_xx; dev[0][1][0]=0.5*Na_yx; dev[1][0][0]=0.5*Na_yx; dev[1][1][0]=0; dev[2][0][0]=0; dev[2][1][0]=0; dev[2][2][0]=0; dev[0][2][0]=0; dev[1][2][0]=0;
						dev[0][0][1]=Na_xy; dev[0][1][1]=0.5*Na_yy; dev[1][0][1]=0.5*Na_yy; dev[1][1][1]=0; dev[2][0][1]=0; dev[2][1][1]=0; dev[2][2][1]=0; dev[0][2][1]=0; dev[1][2][1]=0;
						dev[0][0][2]=0; 	dev[0][1][2]=0; 		dev[1][0][2]=0; 		dev[1][1][2]=0; dev[2][0][2]=0; dev[2][1][2]=0; dev[2][2][2]=0; dev[0][2][2]=0; dev[1][2][2]=0;
						dv[0][0][0]=Na_xx; dv[0][0][1]=Na_xy; dv[0][0][2]=0.0; dv[0][1][0]=Na_yx; dv[0][1][1]=Na_yy; dv[0][1][2]=0.0; dv[0][2][0]=0.0; dv[0][2][1]=0.0; dv[0][2][2]=0.0;
						dv[1][0][0]=0.00;  dv[1][0][1]=0.00;  dv[1][0][2]=0.0; dv[1][1][0]=0.0;   dv[1][1][1]=0.0;   dv[1][1][2]=0.0; dv[1][2][0]=0.0; dv[1][2][1]=0.0; dv[1][2][2]=0.0;
						dv[2][0][0]=0.00;  dv[2][0][1]=0.00;  dv[2][0][2]=0.0; dv[2][1][0]=0.0;   dv[2][1][1]=0.0;   dv[2][1][2]=0.0; dv[2][2][0]=0.0; dv[2][2][1]=0.0; dv[2][2][2]=0.0;
					}
					else if (i==1)
					{
						//dvp[0][0]=0.0;  dvp[0][1]=0.0;  dvp[0][2]=0.0;
						//dvp[1][0]=Na_x; dvp[1][1]=Na_y; dvp[1][2]=0.0;
						//dvp[2][0]=0.0;  dvp[2][1]=0.0;  dvp[2][2]=0.0;
						ev[0][0]=0.0; ev[0][1]=0.5*Na_x; ev[1][0]=0.5*Na_x; ev[1][1]=Na_y; ev[2][0]=0.0; ev[2][1]=0.0; ev[2][2]=0.0; ev[0][2]=0.0; ev[1][2]=0.0;
						dev[0][0][0]=0; dev[0][1][0]=0.5*Na_xx; dev[1][0][0]=0.5*Na_xx; dev[1][1][0]=Na_yx; dev[2][0][0]=0; dev[2][1][0]=0; dev[2][2][0]=0; dev[0][2][0]=0; dev[1][2][0]=0;
						dev[0][0][1]=0; dev[0][1][1]=0.5*Na_xy; dev[1][0][1]=0.5*Na_xy; dev[1][1][1]=Na_yy; dev[2][0][1]=0; dev[2][1][1]=0; dev[2][2][1]=0; dev[0][2][1]=0; dev[1][2][1]=0;
						dev[0][0][2]=0; dev[0][1][2]=0; 		dev[1][0][2]=0; 		dev[1][1][2]=0; 	dev[2][0][2]=0; dev[2][1][2]=0; dev[2][2][2]=0; dev[0][2][2]=0; dev[1][2][2]=0;
						dv[0][0][0]=0.00;  dv[0][0][1]=0.00;  dv[0][0][2]=0.0; dv[0][1][0]=0.0;   dv[0][1][1]=0.0;   dv[0][1][2]=0.0; dv[0][2][0]=0.0; dv[0][2][1]=0.0; dv[0][2][2]=0.0;
						dv[1][0][0]=Na_xx; dv[1][0][1]=Na_xy; dv[1][0][2]=0.0; dv[1][1][0]=Na_yx; dv[1][1][1]=Na_yy; dv[1][1][2]=0.0; dv[1][2][0]=0.0; dv[1][2][1]=0.0; dv[1][2][2]=0.0;
						dv[2][0][0]=0.00;  dv[2][0][1]=0.00;  dv[2][0][2]=0.0; dv[2][1][0]=0.0;   dv[2][1][1]=0.0;   dv[2][1][2]=0.0; dv[2][2][0]=0.0; dv[2][2][1]=0.0; dv[2][2][2]=0.0;
					}

					Feq[a][i] = 0.0;

					for (m=0; m<3; m++)
					{
						for (r=0; r<3; r++)
						{
							for (s=0; s<3; s++)
							{
								//Feq[a][i]+=0.5*M[m]*e[m][r][s]*dvp[s][r];
								Feq[a][i]+=0.25*eps*(fulld_Chi[m][r][s]+fulld_Chi[m][s][r]+fulld_Chi[r][m][s]+fulld_Chi[r][s][m])*dev[m][r][s]
										  +0.25*eps*(fulld_Chi[m][r][s]-fulld_Chi[r][m][s]+fulld_Chi[m][s][r]-fulld_Chi[r][s][m])*dv[m][r][s];
								for (d=0; d<3; d++)
								{
									Feq[a][i]+=C[m][r][s][d]*fullChi[s][d]*ev[m][r];
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

//System for stresses
	#undef  __FUNCT__
	#define __FUNCT__ "Stress"
	PetscErrorCode Stress(IGAPoint p,IGAPoint pU, IGAPoint pZ,IGAPoint pChi,PetscReal *K,PetscReal *F,PetscReal *U,PetscReal *Z, PetscReal *Chi,void *ctx)		//When other fields like Pi or S (etc.) come into this, declare a IGAPoint pPi or pS and a new PescReal *UPi or *US for each
	{
		const PetscReal *N0,(*N1)[2],(*N2)[2][2];
		IGAPointGetShapeFuns(p,0,(const PetscReal**)&N0);									//Value of the shape functions
		IGAPointGetShapeFuns(p,1,(const PetscReal**)&N1);									//Derivatives of the shape functions
		IGAPointGetShapeFuns(p,2,(const PetscReal**)&N2);									//Second derivatives of the shape funcions
		PetscInt a,b,i,j,k,l,m,n,nen=p->nen, dof=p->dof, dofZ=pZ->dof, dimZ=pZ->dim, dofU=pU->dof, dimU=pU->dim, dofChi=pChi->dof;

		PetscReal u[dofU];																	//Array to contain vector u
		PetscReal d_u[dofU][dimU];															//Same for its gradient
		PetscReal d2_u[dofU][dimU][dimU];													//Same for its 2nd order partial derivatives
		//PetscReal d3_u[dofU][dimU][dimU][dimU];												//Same for its 3rd order partial derivatives
		IGAPointFormValue(pU,U,&u[0]);														//Assign u to its container
		IGAPointFormGrad (pU,U,&d_u[0][0]);													//Same for the gradient
		IGAPointFormHess (pU,U,&d2_u[0][0][0]);												//Same for the hessian
		//IGAPointFormDer3 (pU,U,&d3_u[0][0][0][0]);											//This should be the 4-th order tensor u_{i,jkl}

		PetscReal Z0[dofZ];																	//Array to contain the vector z(0)
		PetscReal d_Z0[dofZ][dimZ];															//Same for its gradient
		PetscReal d2_Z0[dof][dimZ][dimZ];													//Same for its 2nd order partial derivatives
		//PetscReal d3_Z0[dofZ][dimZ][dimZ][dimZ];											//Same for its 3rd order partial derivatives
		IGAPointFormValue(pZ,Z,&Z0[0]);														//Assign z to its container
		IGAPointFormGrad (pZ,Z,&d_Z0[0][0]);												//Same for the gradient
		IGAPointFormHess (pZ,Z,&d2_Z0[0][0][0]);											//Same for the hessian
		//IGAPointFormDer3 (pZ,Z,&d3_Z0[0][0][0][0]);											//This should be the 4-th order tensor z_{i,jkl}

		PetscReal chi0[dofChi];																//Array to contain the vector chi(0)
		PetscReal d_Chi0[4][2];																//Same for its gradient
		//PetscReal d2_Chi0[dofChi][dimChi][dimChi];											//Same for its 2nd order partial derivatives
		IGAPointFormValue(pChi,Chi,&chi0[0]);												//Assign chi to its container
		IGAPointFormGrad(pChi,Chi,&d_Chi0[0][0]);											//Assign grad chi to its container
		//IGAPointFormHess (pChi,Chi,&d2_Chi0[0][0][0]);										//This should be the 3-rd order tensor Chi_{i,jk} (remember that we are storing Chi_{kl} as a column vector)

		//Expanding Chi so we can use for loops with correct indexing
		PetscReal fullChi[3][3]={0};
		fullChi[0][0]=chi0[0]; 	fullChi[0][1]=chi0[1];
		fullChi[1][0]=chi0[2]; 	fullChi[1][1]=chi0[3];

		PetscReal fulld_Chi[3][3][3]={0};
		//The four non-zero components of Chi are stored as a vector, restore them to an array with the correct indexing for value and derivative
		fulld_Chi[0][0][0]=d_Chi0[0][0]; fulld_Chi[0][0][1]=d_Chi0[0][1]; fulld_Chi[0][0][2]=0.0;
		fulld_Chi[0][1][0]=d_Chi0[1][0]; fulld_Chi[0][1][1]=d_Chi0[1][1]; fulld_Chi[0][1][2]=0.0;
		fulld_Chi[1][0][0]=d_Chi0[2][0]; fulld_Chi[1][0][1]=d_Chi0[2][1]; fulld_Chi[1][0][2]=0.0;
		fulld_Chi[1][1][0]=d_Chi0[3][0]; fulld_Chi[1][1][1]=d_Chi0[3][1]; fulld_Chi[1][1][2]=0.0;

		//Expanding u and z (and derivatives) to 3 components, more convenient for sums in for loops
		PetscReal fulld_u[3][3]={0};
		fulld_u[0][0]=d_u[0][0]; fulld_u[0][1]=d_u[0][1];
		fulld_u[1][0]=d_u[1][0]; fulld_u[1][1]=d_u[1][1];
		
		PetscReal fulld2_u[3][3][3]={0};
		fulld2_u[0][0][0]=d2_u[0][0][0]; fulld2_u[0][0][1]=d2_u[0][0][1];
		fulld2_u[0][1][0]=d2_u[0][1][0]; fulld2_u[0][1][1]=d2_u[0][1][1];
		fulld2_u[1][0][0]=d2_u[1][0][0]; fulld2_u[1][0][1]=d2_u[1][0][1];
		fulld2_u[1][1][0]=d2_u[1][1][0]; fulld2_u[1][1][1]=d2_u[1][1][1];

		PetscReal fulld_z[3][3]={0};
		fulld_z[0][0]=d_Z0[0][0]; fulld_z[0][1]=d_Z0[0][1];
		fulld_z[1][0]=d_Z0[1][0]; fulld_z[1][1]=d_Z0[1][1];
		
		PetscReal fulld2_z[3][3][3]={0};
		fulld2_z[0][0][0]=d2_Z0[0][0][0]; fulld2_z[0][0][1]=d2_Z0[0][0][1];
		fulld2_z[0][1][0]=d2_Z0[0][1][0]; fulld2_z[0][1][1]=d2_Z0[0][1][1];
		fulld2_z[1][0][0]=d2_Z0[1][0][0]; fulld2_z[1][0][1]=d2_Z0[1][0][1];
		fulld2_z[1][1][0]=d2_Z0[1][1][0]; fulld2_z[1][1][1]=d2_Z0[1][1][1];

		PetscReal (*Kstress)[dof][nen][dof] = (typeof(Kstress)) K;
		PetscReal (*Fstress)[dof] = (PetscReal (*)[dof])F;

		//E and nu should come from AppCtx in the future
		const PetscReal E=2100000.0*9.81*10000.0;
		const PetscReal nu=0.33;
		const PetscReal lambda=(E*nu)/((1.0+nu)*(1.0-2.0*nu));
		const PetscReal mu=E/(2.0*(1.0+nu));
		const PetscReal eps=0.0*1.0/1000.0;														//Choose later based on whatever Amit says :)

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
							Fstress[a][i]+=N0[a]*C[m][n][k][l]*(fulld_u[k][l]-fulld_z[k][l]-fullChi[k][l]);
						}
						Fstress[a][i]+=N1[a][k]*0.25*eps*(fulld2_u[m][n][k]-fulld2_z[m][n][k]-fulld_Chi[m][n][k] +fulld2_u[m][k][n]-fulld2_z[m][k][n]-fulld_Chi[m][k][n]
														 +fulld2_u[n][m][k]-fulld2_z[n][m][k]-fulld_Chi[n][m][k] +fulld2_u[n][k][m]-fulld2_z[n][k][m]-fulld_Chi[n][k][m]);
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
	PetscReal Lx=2.0;
	PetscReal Ly=1.0;
	PetscInt  nx=401;
	PetscInt  ny=401;
	PetscReal h=fmin(Lx/nx,Ly/ny);

	printf("h= %f\n",h);

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

//Creation of types and systems for the L2 projection of S0
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
			//ierr = IGASetBoundaryValue(igaS,dir,side,dof,0.0);CHKERRQ(ierr);    				// Dirichlet boundary conditions
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
	ierr = KSPSetOptionsPrefix(kspl2S,"l2pS_");CHKERRQ(ierr);
	ierr = KSPSetFromOptions(kspl2S);CHKERRQ(ierr);
	//ierr = KSPSetTolerances(kspl2,1e-10,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);CHKERRQ(ierr);
	ierr = KSPSolve(kspl2S,Fl2S,s0);CHKERRQ(ierr);											//This is a simple system, so it can be solved with just this command

	ierr = KSPDestroy(&kspl2S);CHKERRQ(ierr);
	ierr = MatDestroy(&Kl2S);CHKERRQ(ierr);
	ierr = VecDestroy(&Fl2S);CHKERRQ(ierr);
	char nameS[]="/S-2d-0.dat";
	char pathS[512];
	sprintf(pathS,"%s%s",direct,nameS);
	ierr = IGAWriteVec(igaS,s0,pathS);CHKERRQ(ierr);
//

//Creation of types and systems for the L2 projection of Pi=curl(S)
	IGA igaCS;
	ierr = IGACreate(PETSC_COMM_WORLD,&igaCS);CHKERRQ(ierr);
	ierr = IGASetDim(igaCS,2);CHKERRQ(ierr);													//Spatial dimension of the problem
	ierr = IGASetDof(igaCS,4);CHKERRQ(ierr);													//Number of degrees of freedom, per node
	ierr = IGASetOrder(igaCS,1);CHKERRQ(ierr);													//Number of spatial derivatives to calculate
	ierr = IGASetFromOptions(igaCS);CHKERRQ(ierr);												//Note: The order (or degree) of the shape functions is given by the mesh!
	ierr = IGARead(igaCS,"./geometry.dat");CHKERRQ(ierr);
	ierr = IGASetUp(igaCS);CHKERRQ(ierr);
	for (dir=0; dir<2; dir++) 
	{
		for (side=0; side<2; side++) 
		{
			//ierr = IGASetBoundaryValue(igaCS,dir,side,numDof,0.0);CHKERRQ(ierr);    				// Dirichlet boundary conditions
			ierr = IGASetBoundaryForm(igaCS,dir,side,PETSC_TRUE);CHKERRQ(ierr);  				// Neumann boundary conditions
		}
	}

	Mat Kcs;
	Vec cs0,Fcs;

	ierr = IGACreateMat(igaCS,&Kcs);CHKERRQ(ierr);
	ierr = IGACreateVec(igaCS,&cs0);CHKERRQ(ierr);
	ierr = IGACreateVec(igaCS,&Fcs);CHKERRQ(ierr);

	IGAPoint        pointCS,pointS;				//point
	IGAElement      elemCS,elemS;				//element
	PetscReal       *KlocCS,*FlocCS;			//AA y BB
	PetscReal       *KpointCS,*FpointCS;		//KKK y FFF
	const PetscReal *arrayS0CS;					//arrayU
	Vec  			localS0CS;					//localU
	PetscReal       *S0CS;						//U

  	IGAFormSystem  wtfCS;
 	void           *wtf2CS;

 	KSP kspCS;
	ierr = IGACreateKSP(igaCS,&kspCS);CHKERRQ(ierr);

	// Get local vectors S0 and arrays
	ierr = IGAGetLocalVecArray(igaS,s0,&localS0CS,&arrayS0CS);CHKERRQ(ierr);

	// Element loop
	ierr = IGABeginElement(igaCS,&elemCS);CHKERRQ(ierr);
	ierr = IGABeginElement(igaS,&elemS);CHKERRQ(ierr);
	while (IGANextElement(igaCS,elemCS)) 
	{
		IGANextElement(igaS,elemS);
		ierr = IGAElementGetWorkMat(elemCS,&KlocCS);CHKERRQ(ierr);
		ierr = IGAElementGetWorkVec(elemCS,&FlocCS);CHKERRQ(ierr);
		ierr = IGAElementGetValues(elemS,arrayS0CS,&S0CS);CHKERRQ(ierr);

		// FormSystem loop
		while (IGAElementNextFormSystem(elemCS,&wtfCS,&wtf2CS)) 
		{
		// Quadrature loop
			ierr = IGAElementBeginPoint(elemCS,&pointCS);CHKERRQ(ierr);
			ierr = IGAElementBeginPoint(elemS,&pointS);CHKERRQ(ierr);
			if(pointCS->atboundary==0 && pointS->atboundary==0)
			{
				while (IGAElementNextPoint(elemCS,pointCS))
				{
					IGAElementNextPoint(elemS,pointS);
					ierr = IGAPointGetWorkMat(pointCS,&KpointCS);CHKERRQ(ierr);
					ierr = IGAPointGetWorkVec(pointCS,&FpointCS);CHKERRQ(ierr);
					ierr = curlSPi(pointCS,pointS,KpointCS,FpointCS,S0CS,NULL);CHKERRQ(ierr);
					ierr = IGAPointAddMat(pointCS,KpointCS,KlocCS);CHKERRQ(ierr);
					ierr = IGAPointAddVec(pointCS,FpointCS,FlocCS);CHKERRQ(ierr);
				}
				IGAElementNextPoint(elemS,pointS);
				ierr = IGAElementEndPoint(elemCS,&pointCS);CHKERRQ(ierr);
				ierr = IGAElementEndPoint(elemS,&pointS);CHKERRQ(ierr);
			}
		}
		ierr = IGAElementFixSystem(elemCS,KlocCS,FlocCS);CHKERRQ(ierr);					//This sets Dirichlet condition ¿?
		ierr = IGAElementAssembleMat(elemCS,KlocCS,Kcs);CHKERRQ(ierr);
		ierr = IGAElementAssembleVec(elemCS,FlocCS,Fcs);CHKERRQ(ierr);

	}
	IGANextElement(igaS,elemS);
	ierr = IGAEndElement(igaCS,&elemCS);CHKERRQ(ierr);
	ierr = IGAEndElement(igaS,&elemS);CHKERRQ(ierr);

	// Restore local vectors S0 and arrays
	ierr = IGARestoreLocalVecArray(igaS,s0,&localS0CS,&arrayS0CS);CHKERRQ(ierr);

	ierr = MatAssemblyBegin(Kcs,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyEnd  (Kcs,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = VecAssemblyBegin(Fcs);CHKERRQ(ierr);
	ierr = VecAssemblyEnd  (Fcs);CHKERRQ(ierr);

	ierr = KSPSetOperators(kspCS,Kcs,Kcs);CHKERRQ(ierr);
	ierr = KSPSetFromOptions(kspCS);CHKERRQ(ierr);
	ierr = KSPSolve(kspCS,Fcs,cs0);CHKERRQ(ierr);

	ierr = KSPDestroy(&kspCS);CHKERRQ(ierr);
	ierr = MatDestroy(&Kcs);CHKERRQ(ierr);
	ierr = VecDestroy(&Fcs);CHKERRQ(ierr);
	//ierr = IGAWriteVec(igaCS,cs0,"./results/CS-2d-0.dat");CHKERRQ(ierr);
	char nameCS[]="/CS-2d-0.dat";
	char pathCS[512];
	sprintf(pathCS,"%s%s",direct,nameCS);
	ierr = IGAWriteVec(igaCS,cs0,pathCS);CHKERRQ(ierr);
	
//

//Creation of types and systems for the L2 projection of Pi0
	//System for Pi
	printf("System for L2 projection for Pi starting \n");
	IGA igaPi;
	ierr = IGACreate(PETSC_COMM_WORLD,&igaPi);CHKERRQ(ierr);
	ierr = IGASetDim(igaPi,2);CHKERRQ(ierr);													//Spatial dimension of the problem
	ierr = IGASetDof(igaPi,4);CHKERRQ(ierr);													//Number of degrees of freedom, per node
	ierr = IGASetOrder(igaPi,1);CHKERRQ(ierr);													//Number of spatial derivatives to calculate
	ierr = IGASetFromOptions(igaPi);CHKERRQ(ierr);												//Note: The order (or degree) of the shape functions is given by the mesh!
	ierr = IGARead(igaPi,"./geometry.dat");CHKERRQ(ierr);
	ierr = IGASetUp(igaPi);CHKERRQ(ierr);
	
	for (dir=0; dir<2; dir++) 
	{
		for (side=0; side<2; side++) 
		{
			//ierr = IGASetBoundaryValue(iga,dir,side,dof,0.0);CHKERRQ(ierr);    				// Dirichlet boundary conditions
			ierr = IGASetBoundaryForm(igaPi,dir,side,PETSC_TRUE);CHKERRQ(ierr);  				// Neumann boundary conditions
		}
	}

	Vec pi0;
	Mat Kl2Pi;
	Vec Fl2Pi;
	ierr = IGACreateVec(igaPi,&pi0);CHKERRQ(ierr);  
	ierr = IGACreateMat(igaPi,&Kl2Pi);CHKERRQ(ierr);
	ierr = IGACreateVec(igaPi,&Fl2Pi);CHKERRQ(ierr);
	ierr = IGASetFormSystem(igaPi,L2ProjectionPi,&userL2);CHKERRQ(ierr);
	ierr = IGAComputeSystem(igaPi,Kl2Pi,Fl2Pi);CHKERRQ(ierr);

	//This parts set and calls KSP to solve the linear system
	KSP kspl2Pi;
	ierr = IGACreateKSP(igaPi,&kspl2Pi);CHKERRQ(ierr);										
	ierr = KSPSetOperators(kspl2Pi,Kl2Pi,Kl2Pi);CHKERRQ(ierr); 								//This function creates the matrix for the system on the second parameter and uses the 3rd parameter as a preconditioner
	ierr = KSPSetType(kspl2Pi,KSPCG);CHKERRQ(ierr);											//Using KSPCG (conjugated gradient) because the matrix is symmetric
	ierr = KSPSetOptionsPrefix(kspl2Pi,"l2pPi_");CHKERRQ(ierr);
	ierr = KSPSetFromOptions(kspl2Pi);CHKERRQ(ierr);
	//ierr = KSPSetTolerances(kspl2,1e-10,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);CHKERRQ(ierr);
	ierr = KSPSolve(kspl2Pi,Fl2Pi,pi0);CHKERRQ(ierr);										//This is a simple system, so it can be solved with just this command

	ierr = KSPDestroy(&kspl2Pi);CHKERRQ(ierr);
	ierr = MatDestroy(&Kl2Pi);CHKERRQ(ierr);
	ierr = VecDestroy(&Fl2Pi);CHKERRQ(ierr);
	char namePi[]="/Pi-2d-0.dat";
	char pathPi[512];
	sprintf(pathPi,"%s%s",direct,namePi);
	ierr = IGAWriteVec(igaPi,pi0,pathPi);CHKERRQ(ierr);
//

//Creation of types and systems for the L2 projection of Alfa0
	//System for Alfa
	printf("System for L2 projection for Pi starting \n");
	IGA igaAl;
	ierr = IGACreate(PETSC_COMM_WORLD,&igaAl);CHKERRQ(ierr);
	ierr = IGASetDim(igaAl,2);CHKERRQ(ierr);													//Spatial dimension of the problem
	ierr = IGASetDof(igaAl,2);CHKERRQ(ierr);													//Number of degrees of freedom, per node
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

	Vec al0;
	Mat Kl2Al;
	Vec Fl2Al;
	ierr = IGACreateVec(igaAl,&al0);CHKERRQ(ierr);  
	ierr = IGACreateMat(igaAl,&Kl2Al);CHKERRQ(ierr);
	ierr = IGACreateVec(igaAl,&Fl2Al);CHKERRQ(ierr);
	ierr = IGASetFormSystem(igaAl,L2ProjectionAl,&userL2);CHKERRQ(ierr);
	ierr = IGAComputeSystem(igaAl,Kl2Al,Fl2Al);CHKERRQ(ierr);

	//This parts set and calls KSP to solve the linear system
	KSP kspl2Al;
	ierr = IGACreateKSP(igaAl,&kspl2Al);CHKERRQ(ierr);
	ierr = KSPSetOperators(kspl2Al,Kl2Al,Kl2Al);CHKERRQ(ierr); 								//This function creates the matrix for the system on the second parameter and uses the 3rd parameter as a preconditioner
	ierr = KSPSetType(kspl2Al,KSPCG);CHKERRQ(ierr);											//Using KSPCG (conjugated gradient) because the matrix is symmetric
	ierr = KSPSetOptionsPrefix(kspl2Al,"l2pAl_");CHKERRQ(ierr);
	ierr = KSPSetFromOptions(kspl2Al);CHKERRQ(ierr);
	//ierr = KSPSetTolerances(kspl2,1e-10,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);CHKERRQ(ierr);
	ierr = KSPSolve(kspl2Al,Fl2Al,al0);CHKERRQ(ierr);										//This is a simple system, so it can be solved with just this command

	ierr = KSPDestroy(&kspl2Al);CHKERRQ(ierr);
	ierr = MatDestroy(&Kl2Al);CHKERRQ(ierr);
	ierr = VecDestroy(&Fl2Al);CHKERRQ(ierr);
	//ierr = IGAWriteVec(igaAl,al0,"./results/Al-2d-0.dat");CHKERRQ(ierr);
	char nameAl[]="/Al-2d-0.dat";
	char pathAl[512];
	sprintf(pathAl,"%s%s",direct,nameAl);
	ierr = IGAWriteVec(igaAl,al0,pathAl);CHKERRQ(ierr);
//

//Creation of types and systems for the initial state of Chi
	//System for Chi
	printf("System for Chi starting \n");
	IGA igaChi;
	ierr = IGACreate(PETSC_COMM_WORLD,&igaChi);CHKERRQ(ierr);
	ierr = IGASetDim(igaChi,2);CHKERRQ(ierr);													//Spatial dimension of the problem
	ierr = IGASetDof(igaChi,4);CHKERRQ(ierr);													//Number of degrees of freedom, per node
	ierr = IGASetOrder(igaChi,3);CHKERRQ(ierr);													//Number of spatial derivatives to calculate
	ierr = IGASetFromOptions(igaChi);CHKERRQ(ierr);												//Note: The order (or degree) of the shape functions is given by the mesh!
	//ierr = IGARead(igaChi,"./geometry.dat");CHKERRQ(ierr);
	ierr = IGARead(igaChi,"./geometry2.dat");CHKERRQ(ierr);
	ierr = IGASetUp(igaChi);CHKERRQ(ierr);

	//ierr = IGASetBoundaryValue(iga,dir,side,dof,value);CHKERRQ(ierr);							//Syntax
	ierr = IGASetBoundaryValue(igaChi,0 ,0 ,0 ,0.0);CHKERRQ(ierr);    							// Dirichlet boundary conditions x-dir, left side, chi(1,1)=0
	ierr = IGASetBoundaryValue(igaChi,0, 1, 0, 0.0);CHKERRQ(ierr);    							// Dirichlet boundary conditions x-dir, right side, chi(1,1)=0
	ierr = IGASetBoundaryValue(igaChi,0 ,0 ,2 ,0.0);CHKERRQ(ierr);    							// Dirichlet boundary conditions x-dir, left side, chi(2,1)=0
	ierr = IGASetBoundaryValue(igaChi,0, 1, 2, 0.0);CHKERRQ(ierr);    							// Dirichlet boundary conditions x-dir, right side, chi(2,1)=0
	ierr = IGASetBoundaryValue(igaChi,1, 0, 1, 0.0);CHKERRQ(ierr);    							// Dirichlet boundary conditions y-dir, lower side, chi(1,2)=0
	ierr = IGASetBoundaryValue(igaChi,1, 1, 1, 0.0);CHKERRQ(ierr);    							// Dirichlet boundary conditions y-dir, upper side, chi(1,2)=0
	ierr = IGASetBoundaryValue(igaChi,1, 0, 3, 0.0);CHKERRQ(ierr);    							// Dirichlet boundary conditions y-dir, lower side, chi(2,2)=0
	ierr = IGASetBoundaryValue(igaChi,1, 1, 3, 0.0);CHKERRQ(ierr);    							// Dirichlet boundary conditions y-dir, upper side, chi(2,2)=0
	//ierr = IGASetBoundaryForm(igaChi,dir,side,PETSC_TRUE);CHKERRQ(ierr);						// Neumann boundary conditions

	Mat Kchi;
	Vec chi0,Fchi;

	ierr = IGACreateMat(igaChi,&Kchi);CHKERRQ(ierr);
	ierr = IGACreateVec(igaChi,&chi0);CHKERRQ(ierr);
	ierr = IGACreateVec(igaChi,&Fchi);CHKERRQ(ierr);

	IGAPoint        pointChi,pointAl;			//point
	IGAElement      elemChi,elemAl;				//element
	PetscReal       *KlocChi,*FlocChi;			//AA y BB
	PetscReal       *KpointChi,*FpointChi;		//KKK y FFF
	const PetscReal *arrayAl0Chi;				//arrayU
	Vec  			localAl0Chi;				//localU
	PetscReal       *Al0Chi;					//U

  	IGAFormSystem  wtfChi;
 	void           *wtf2Chi;

 	KSP kspChi;
	ierr = IGACreateKSP(igaChi,&kspChi);CHKERRQ(ierr);

	// Get local vectors Al0 and arrays
	ierr = IGAGetLocalVecArray(igaAl,al0,&localAl0Chi,&arrayAl0Chi);CHKERRQ(ierr);

	// Element loop
	ierr = IGABeginElement(igaChi,&elemChi);CHKERRQ(ierr);
	ierr = IGABeginElement(igaAl,&elemAl);CHKERRQ(ierr);
	while (IGANextElement(igaChi,elemChi)) 
	{
		IGANextElement(igaAl,elemAl);
		ierr = IGAElementGetWorkMat(elemChi,&KlocChi);CHKERRQ(ierr);
		ierr = IGAElementGetWorkVec(elemChi,&FlocChi);CHKERRQ(ierr);
		ierr = IGAElementGetValues(elemAl,arrayAl0Chi,&Al0Chi);CHKERRQ(ierr);

		// FormSystem loop
		while (IGAElementNextFormSystem(elemChi,&wtfChi,&wtf2Chi)) 
		{
		// Quadrature loop
			ierr = IGAElementBeginPoint(elemChi,&pointChi);CHKERRQ(ierr);
			ierr = IGAElementBeginPoint(elemAl,&pointAl);CHKERRQ(ierr);
			while (IGAElementNextPoint(elemChi,pointChi)) 
			{

				IGAElementNextPoint(elemAl,pointAl);
				ierr = IGAPointGetWorkMat(pointChi,&KpointChi);CHKERRQ(ierr);
				ierr = IGAPointGetWorkVec(pointChi,&FpointChi);CHKERRQ(ierr);
				ierr = Chi(pointChi,pointAl,KpointChi,FpointChi,Al0Chi,NULL);CHKERRQ(ierr);
				ierr = IGAPointAddMat(pointChi,KpointChi,KlocChi);CHKERRQ(ierr);
				ierr = IGAPointAddVec(pointChi,FpointChi,FlocChi);CHKERRQ(ierr);
			}
			IGAElementNextPoint(elemAl,pointAl);
			ierr = IGAElementEndPoint(elemChi,&pointChi);CHKERRQ(ierr);
			
			//printf("pointAl,index=%d\n",pointAl->index);
			while (pointAl->index != -1)
			{
				IGAElementNextPoint(elemAl,pointAl);
			}
			//printf("pointAl,index=%d\n",pointAl->index);

			ierr = IGAElementEndPoint(elemAl,&pointAl);CHKERRQ(ierr);
		}
		ierr = IGAElementFixSystem(elemChi,KlocChi,FlocChi);CHKERRQ(ierr);					//This sets Dirichlet condition ¿?
		ierr = IGAElementAssembleMat(elemChi,KlocChi,Kchi);CHKERRQ(ierr);
		ierr = IGAElementAssembleVec(elemChi,FlocChi,Fchi);CHKERRQ(ierr);

	}
	IGANextElement(igaAl,elemAl);
	ierr = IGAEndElement(igaChi,&elemChi);CHKERRQ(ierr);
	ierr = IGAEndElement(igaAl,&elemAl);CHKERRQ(ierr);

	// Restore local vectors Al0 and arrays
	ierr = IGARestoreLocalVecArray(igaAl,al0,&localAl0Chi,&arrayAl0Chi);CHKERRQ(ierr);

	ierr = MatAssemblyBegin(Kchi,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyEnd  (Kchi,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = VecAssemblyBegin(Fchi);CHKERRQ(ierr);
	ierr = VecAssemblyEnd  (Fchi);CHKERRQ(ierr);

	ierr = KSPSetOperators(kspChi,Kchi,Kchi);CHKERRQ(ierr);
	ierr = KSPSetFromOptions(kspChi);CHKERRQ(ierr);
	ierr = KSPSolve(kspChi,Fchi,chi0);CHKERRQ(ierr);

	ierr = KSPDestroy(&kspChi);CHKERRQ(ierr);
	ierr = MatDestroy(&Kchi);CHKERRQ(ierr);
	ierr = VecDestroy(&Fchi);CHKERRQ(ierr);
	//ierr = IGAWriteVec(igaChi,chi0,"./results/Chi-2d-0.dat");CHKERRQ(ierr);
	char nameChi[]="/Chi-2d-0.dat";
	char pathChi[512];
	sprintf(pathChi,"%s%s",direct,nameChi);
	ierr = IGAWriteVec(igaChi,chi0,pathChi);CHKERRQ(ierr);
	
//

//Creation of types and systems for z
	//System for z\dot (div(grad(z))=curl(alphaXV_alpha))
	printf("System for z(dot) starting \n");
	IGA igaZ;
	ierr = IGACreate(PETSC_COMM_WORLD,&igaZ);CHKERRQ(ierr);
	ierr = IGASetDim(igaZ,2);CHKERRQ(ierr);														//Spatial dimension of the problem
	ierr = IGASetDof(igaZ,2);CHKERRQ(ierr);														//Number of degrees of freedom, per node
	ierr = IGASetOrder(igaZ,1);CHKERRQ(ierr);													//Number of spatial derivatives to calculate
	ierr = IGASetFromOptions(igaZ);CHKERRQ(ierr);												//Note: The order (or degree) of the shape functions is given by the mesh!
	ierr = IGARead(igaZ,"./geometry.dat");CHKERRQ(ierr);
	ierr = IGASetUp(igaZ);CHKERRQ(ierr);
	//PetscInt dir,side;
	for (dir=0; dir<2; dir++) 
	{
		for (side=0; side<2; side++) 
		{
			//ierr = IGASetBoundaryValue(iga,dir,side,dof,0.0);CHKERRQ(ierr);    				// Dirichlet boundary conditions
			ierr = IGASetBoundaryForm(igaZ,dir,side,PETSC_TRUE);CHKERRQ(ierr);  				// Neumann boundary conditions
		}
	}
	
	Mat Kz;
	Vec z0,Fz;
	ierr = IGACreateMat(igaZ,&Kz);CHKERRQ(ierr);
	ierr = IGACreateVec(igaZ,&z0);CHKERRQ(ierr);
	ierr = IGACreateVec(igaZ,&Fz);CHKERRQ(ierr);

	IGAPoint        pointZ;
	IGAElement      elemZ;					//element
	PetscReal       *KlocZ,*FlocZ;			//AA y BB
	PetscReal       *KpointZ,*FpointZ;		//KKK y FFF
	const PetscReal *arrayAl0Z;				//arrayU
	Vec  			localAl0Z;				//localU
	PetscReal       *Al0Z;					//U

  	IGAFormSystem  wtfZ;
 	void           *wtf2Z;

 	KSP kspZ;
	ierr = IGACreateKSP(igaZ,&kspZ);CHKERRQ(ierr);

	//Get local vectors Al0 and arrays
	ierr = IGAGetLocalVecArray(igaAl,al0,&localAl0Z,&arrayAl0Z);CHKERRQ(ierr);

	//Element loop
	ierr = IGABeginElement(igaZ,&elemZ);CHKERRQ(ierr);
	ierr = IGABeginElement(igaAl,&elemAl);CHKERRQ(ierr);

	while (IGANextElement(igaZ,elemZ))
	{
		IGANextElement(igaAl,elemAl);
		ierr = IGAElementGetWorkMat(elemZ,&KlocZ);CHKERRQ(ierr);
		ierr = IGAElementGetWorkVec(elemZ,&FlocZ);CHKERRQ(ierr);
		ierr = IGAElementGetValues(elemAl,arrayAl0Z,&Al0Z);CHKERRQ(ierr);

		//FormSystem loop
		while (IGAElementNextFormSystem(elemZ,&wtfZ,&wtf2Z)) 
		{
			//Quadrature loop
			ierr = IGAElementBeginPoint(elemZ,&pointZ);CHKERRQ(ierr);
			ierr = IGAElementBeginPoint(elemAl,&pointAl);CHKERRQ(ierr);
			if(pointZ->atboundary==0 && pointAl->atboundary==0)
			{
				while (IGAElementNextPoint(elemZ,pointZ)) 
				{
					IGAElementNextPoint(elemAl,pointAl);
					ierr = IGAPointGetWorkMat(pointZ,&KpointZ);CHKERRQ(ierr);
					ierr = IGAPointGetWorkVec(pointZ,&FpointZ);CHKERRQ(ierr);
					ierr = divGradZ(pointZ,pointAl,KpointZ,FpointZ,Al0Z,NULL);CHKERRQ(ierr);
					ierr = IGAPointAddMat(pointZ,KpointZ,KlocZ);CHKERRQ(ierr);
					ierr = IGAPointAddVec(pointZ,FpointZ,FlocZ);CHKERRQ(ierr);
				}
				IGAElementNextPoint(elemAl,pointAl);
				ierr = IGAElementEndPoint(elemZ,&pointZ);CHKERRQ(ierr);
				ierr = IGAElementEndPoint(elemAl,&pointAl);CHKERRQ(ierr);
			}
		}
		ierr = IGAElementFixSystem(elemZ,KlocZ,FlocZ);CHKERRQ(ierr);					//This sets Dirichlet condition ¿?
		ierr = IGAElementAssembleMat(elemZ,KlocZ,Kz);CHKERRQ(ierr);
		ierr = IGAElementAssembleVec(elemZ,FlocZ,Fz);CHKERRQ(ierr);

	}
	IGANextElement(igaAl,elemAl);
	ierr = IGAEndElement(igaZ,&elemZ);CHKERRQ(ierr);
	ierr = IGAEndElement(igaAl,&elemAl);CHKERRQ(ierr);

	// Restore local vectors Al0 and arrays
	ierr = IGARestoreLocalVecArray(igaAl,al0,&localAl0Z,&arrayAl0Z);CHKERRQ(ierr);

	//Here we set values to the matrix directly (to impose Dirichlet condition in a single point)
	ierr = MatAssemblyBegin(Kz,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyEnd  (Kz,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	
	PetscInt mz,nz;
	ierr = MatGetSize(Kz,&mz,&nz);CHKERRQ(ierr);

	PetscInt rows=0;
	PetscReal vals; 

	rows=0;
	ierr = MatZeroRows(Kz,1,&rows,1.0e16,0,0);
	rows=1;
	ierr = MatZeroRows(Kz,1,&rows,1.0e16,0,0);

	ierr = MatAssemblyBegin(Kz,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyEnd  (Kz,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

	//Here we set values to the vector directly (to impose Dirichlet condition in a single point)
	ierr = VecAssemblyBegin(Fz);CHKERRQ(ierr);
	ierr = VecAssemblyEnd  (Fz);CHKERRQ(ierr);
	rows=0;
	vals=0.0;
	ierr = VecSetValue(Fz,rows,vals,INSERT_VALUES);CHKERRQ(ierr);
	rows=1;
	vals=0.0;
	ierr = VecSetValue(Fz,rows,vals,INSERT_VALUES);CHKERRQ(ierr);
	ierr = VecAssemblyBegin(Fz);CHKERRQ(ierr);
	ierr = VecAssemblyEnd  (Fz);CHKERRQ(ierr);

	//Now finish forming system and solve
	ierr = KSPSetOperators(kspZ,Kz,Kz);CHKERRQ(ierr);
	ierr = KSPSetFromOptions(kspZ);CHKERRQ(ierr);
	ierr = KSPSolve(kspZ,Fz,z0);CHKERRQ(ierr);

	ierr = KSPDestroy(&kspZ);CHKERRQ(ierr);
	ierr = MatDestroy(&Kz);CHKERRQ(ierr);
	ierr = VecDestroy(&Fz);CHKERRQ(ierr);
	//ierr = IGAWriteVec(igaZ,z0,"./results/Z-2d-0.dat");CHKERRQ(ierr);
	char nameZ[]="/Zdot-2d-0.dat";
	char pathZ[512];
	sprintf(pathZ,"%s%s",direct,nameZ);
	ierr = IGAWriteVec(igaZ,z0,pathZ);CHKERRQ(ierr);	
//

//Creation of types and system for Pi\dot=-curl(Pi X V_pi)
	Mat KPi;
	Vec xPi,FPi;
	ierr = IGACreateMat(igaPi,&KPi);CHKERRQ(ierr);
	ierr = IGACreateVec(igaPi,&xPi);CHKERRQ(ierr);
	ierr = IGACreateVec(igaPi,&FPi);CHKERRQ(ierr);
	//ierr = IGASetFormSystem(igaPi,System,NULL);CHKERRQ(ierr);
	//ierr = IGAComputeSystem(igaPi,KPi,FPi);CHKERRQ(ierr);

	AppCtx userPi;
	userPi.dt     = dt;

	IGAPoint        pointPi;
	IGAElement      elemPi;      			//element
	PetscReal       *KlocPi,*FlocPi;     	//AA y BB
	PetscReal       *KpointPi,*FpointPi;	//KKK y FFF
	const PetscReal *arrayPi0;				//arrayU
	Vec  localPi0;							//localU
	PetscReal       *Pi0;					//U

	IGAFormSystem  wtfPi;
	void           *wtf2Pi;

	KSP kspPi;
	ierr = IGACreateKSP(igaPi,&kspPi);CHKERRQ(ierr);

	char nombrePi[24];
//

//Creation of types and system for Alpha\dot=-curl(Alpha X V_alpha)
	Mat KAl;
	Vec xAl,FAl;
	ierr = IGACreateMat(igaAl,&KAl);CHKERRQ(ierr);
	ierr = IGACreateVec(igaAl,&xAl);CHKERRQ(ierr);
	ierr = IGACreateVec(igaAl,&FAl);CHKERRQ(ierr);
	//ierr = IGASetFormSystem(igaAl,System,NULL);CHKERRQ(ierr);
	//ierr = IGAComputeSystem(igaAl,KAl,FAl);CHKERRQ(ierr);

	AppCtx userAl;
	userAl.dt     = dt;

	//IGAPoint        pointAl;
	//IGAElement      elemAl;
	PetscReal       *KlocAl,*FlocAl;
	PetscReal       *KpointAl,*FpointAl;
	const PetscReal *arrayAl0;
	Vec  localAl0;
	PetscReal       *Al0;

	IGAFormSystem  wtfAl;
	void           *wtf2Al;

	KSP kspAl;
	ierr = IGACreateKSP(igaAl,&kspAl);CHKERRQ(ierr);

	char nombreAl[24];
//

//Solution loop
	PetscInt i;
	for (i=0; i<numStep; i++) 
	{
		PetscPrintf(PETSC_COMM_WORLD,"Iter= %d \n",(i+1));
		//Reset variables after each loop
		ierr = MatZeroEntries(KPi);CHKERRQ(ierr);
		ierr = VecZeroEntries(FPi);CHKERRQ(ierr);

		ierr = MatZeroEntries(KAl);CHKERRQ(ierr);
		ierr = VecZeroEntries(FAl);CHKERRQ(ierr);

		// Get local vectors Pi0 and arrays
		ierr = IGAGetLocalVecArray(igaPi,pi0,&localPi0,&arrayPi0);CHKERRQ(ierr);

		// Get local vectors Al0 and arrays
		ierr = IGAGetLocalVecArray(igaAl,al0,&localAl0,&arrayAl0);CHKERRQ(ierr);

		// Solution cycle for Pi
		// Element loop
		ierr = IGABeginElement(igaPi,&elemPi);CHKERRQ(ierr);
		while (IGANextElement(igaPi,elemPi)) 
		{
			ierr = IGAElementGetWorkMat(elemPi,&KlocPi);CHKERRQ(ierr);
			ierr = IGAElementGetWorkVec(elemPi,&FlocPi);CHKERRQ(ierr);
			ierr = IGAElementGetValues(elemPi,arrayPi0,&Pi0);CHKERRQ(ierr);

			// FormSystem loop
			while (IGAElementNextFormSystem(elemPi,&wtfPi,&wtf2Pi)) 
			{
			// Quadrature loop
				ierr = IGAElementBeginPoint(elemPi,&pointPi);CHKERRQ(ierr);
				while (IGAElementNextPoint(elemPi,pointPi)) 
				{
					ierr = IGAPointGetWorkMat(pointPi,&KpointPi);CHKERRQ(ierr);
					ierr = IGAPointGetWorkVec(pointPi,&FpointPi);CHKERRQ(ierr);
					ierr = Pi(pointPi,KpointPi,FpointPi,Pi0,&userPi);CHKERRQ(ierr);
					ierr = IGAPointAddMat(pointPi,KpointPi,KlocPi);CHKERRQ(ierr);
					ierr = IGAPointAddVec(pointPi,FpointPi,FlocPi);CHKERRQ(ierr);
				}
				ierr = IGAElementEndPoint(elemPi,&pointPi);CHKERRQ(ierr);
			}
			ierr = IGAElementFixSystem(elemPi,KlocPi,FlocPi);CHKERRQ(ierr);					//This sets Dirichlet condition ¿?
			ierr = IGAElementAssembleMat(elemPi,KlocPi,KPi);CHKERRQ(ierr);
			ierr = IGAElementAssembleVec(elemPi,FlocPi,FPi);CHKERRQ(ierr);

		}
		ierr = IGAEndElement(igaPi,&elemPi);CHKERRQ(ierr);

		// Restore local vectors Pi0 and arrays
		ierr = IGARestoreLocalVecArray(igaPi,pi0,&localPi0,&arrayPi0);CHKERRQ(ierr);

		ierr = MatAssemblyBegin(KPi,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
		ierr = MatAssemblyEnd  (KPi,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
		ierr = VecAssemblyBegin(FPi);CHKERRQ(ierr);
		ierr = VecAssemblyEnd  (FPi);CHKERRQ(ierr);

		ierr = KSPSetOperators(kspPi,KPi,KPi);CHKERRQ(ierr);
		ierr = KSPSetFromOptions(kspPi);CHKERRQ(ierr);
		ierr = KSPSolve(kspPi,FPi,xPi);CHKERRQ(ierr);

		// Solution cycle for Alpha
		// Element loop
		ierr = IGABeginElement(igaAl,&elemAl);CHKERRQ(ierr);
		while (IGANextElement(igaAl,elemAl)) 
		{
			ierr = IGAElementGetWorkMat(elemAl,&KlocAl);CHKERRQ(ierr);
			ierr = IGAElementGetWorkVec(elemAl,&FlocAl);CHKERRQ(ierr);
			ierr = IGAElementGetValues(elemAl,arrayAl0,&Al0);CHKERRQ(ierr);

			// FormSystem loop
			while (IGAElementNextFormSystem(elemAl,&wtfAl,&wtf2Al)) 
			{
				// Quadrature loop
				ierr = IGAElementBeginPoint(elemAl,&pointAl);CHKERRQ(ierr);
				while (IGAElementNextPoint(elemAl,pointAl)) 
				{
					ierr = IGAPointGetWorkMat(pointAl,&KpointAl);CHKERRQ(ierr);
					ierr = IGAPointGetWorkVec(pointAl,&FpointAl);CHKERRQ(ierr);
					ierr = Alfa(pointAl,KpointAl,FpointAl,Al0,&userAl);CHKERRQ(ierr);
					ierr = IGAPointAddMat(pointAl,KpointAl,KlocAl);CHKERRQ(ierr);
					ierr = IGAPointAddVec(pointAl,FpointAl,FlocAl);CHKERRQ(ierr);
				}
				ierr = IGAElementEndPoint(elemAl,&pointAl);CHKERRQ(ierr);
			}
			ierr = IGAElementFixSystem(elemAl,KlocAl,FlocAl);CHKERRQ(ierr);					//This sets Dirichlet condition ¿?
			ierr = IGAElementAssembleMat(elemAl,KlocAl,KAl);CHKERRQ(ierr);
			ierr = IGAElementAssembleVec(elemAl,FlocAl,FAl);CHKERRQ(ierr);

		}
		ierr = IGAEndElement(igaAl,&elemAl);CHKERRQ(ierr);

		// Restore local vectors Pi0 and arrays
		ierr = IGARestoreLocalVecArray(igaAl,al0,&localAl0,&arrayAl0);CHKERRQ(ierr);

		ierr = MatAssemblyBegin(KAl,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
		ierr = MatAssemblyEnd  (KAl,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
		ierr = VecAssemblyBegin(FAl);CHKERRQ(ierr);
		ierr = VecAssemblyEnd  (FAl);CHKERRQ(ierr);

		ierr = KSPSetOperators(kspAl,KAl,KAl);CHKERRQ(ierr);
		ierr = KSPSetFromOptions(kspAl);CHKERRQ(ierr);
		ierr = KSPSolve(kspAl,FAl,xAl);CHKERRQ(ierr);

		if ( ((i+1)%4==0) || ((i+1)==numStep) )
		{
			sprintf(nombrePi, "/Pi-2d-%d.dat", (i+1));											//Adds iteration number to string
			char pathPi2[512];
			sprintf(pathPi2,"%s%s",direct,nombrePi);
			ierr = IGAWriteVec(igaPi,xPi,pathPi2);CHKERRQ(ierr);								//Saves result to disk
			
			sprintf(nombreAl, "/Al-2d-%d.dat", (i+1));											//Adds iteration number to string
			char pathAl2[512];
			sprintf(pathAl2,"%s%s",direct,nombreAl);
			ierr = IGAWriteVec(igaAl,xAl,pathAl2);CHKERRQ(ierr);								//Saves result to disk
		}
		ierr=VecCopy(xPi,pi0); CHKERRQ(ierr);													//Copy Pi(t+1) to Pi(t)
		ierr=VecCopy(xAl,al0); CHKERRQ(ierr);													//Copy Al(t+1) to Al(t)
	}

	//ierr = VecView(x,PETSC_VIEWER_DRAW_WORLD);CHKERRQ(ierr);

	ierr = KSPDestroy(&kspPi);CHKERRQ(ierr);
	ierr = MatDestroy(&KPi);CHKERRQ(ierr);
	ierr = VecDestroy(&xPi);CHKERRQ(ierr);
	ierr = VecDestroy(&FPi);CHKERRQ(ierr);
	

	ierr = KSPDestroy(&kspAl);CHKERRQ(ierr);
	ierr = MatDestroy(&KAl);CHKERRQ(ierr);
	ierr = VecDestroy(&xAl);CHKERRQ(ierr);
	ierr = VecDestroy(&FAl);CHKERRQ(ierr);

//

//System for initial state of z
	printf("System for Z0 starting \n");
	IGA igaZ0;
	ierr = IGACreate(PETSC_COMM_WORLD,&igaZ0);CHKERRQ(ierr);
	ierr = IGASetDim(igaZ0,2);CHKERRQ(ierr);													//Spatial dimension of the problem
	ierr = IGASetDof(igaZ0,2);CHKERRQ(ierr);													//Number of degrees of freedom, per node
	ierr = IGASetOrder(igaZ0,2);CHKERRQ(ierr);													//Number of spatial derivatives to calculate
	ierr = IGASetFromOptions(igaZ0);CHKERRQ(ierr);												//Note: The order (or degree) of the shape functions is given by the mesh!
	ierr = IGARead(igaZ0,"./geometry2.dat");CHKERRQ(ierr);
	ierr = IGASetUp(igaZ0);CHKERRQ(ierr);

	Mat KZ0;
	Vec Z0,FZ0;

	ierr = IGACreateMat(igaZ0,&KZ0);CHKERRQ(ierr);
	ierr = IGACreateVec(igaZ0,&Z0);CHKERRQ(ierr);
	ierr = IGACreateVec(igaZ0,&FZ0);CHKERRQ(ierr);

	IGAPoint		pointZ0;					//point
	IGAElement		elemZ0;						//element
	PetscReal		*KlocZ0,*FlocZ0;			//AA y BB
	PetscReal		*KpointZ0,*FpointZ0;		//KKK y FFF
	const PetscReal	*arrayChi0Z0;				//arrayU
	Vec				localChi0Z0;				//localU
	PetscReal		*Chi0Z0;					//U

  	IGAFormSystem	wtfZ0;
 	void			*wtf2Z0;

 	KSP kspZ0;
	ierr = IGACreateKSP(igaZ0,&kspZ0);CHKERRQ(ierr);

	//ierr = IGASetBoundaryValue(iga,dir,side,dof,value);CHKERRQ(ierr);							//Syntax
	//ierr = IGASetBoundaryValue(igaZ0,1 ,0 ,0 ,0.0);CHKERRQ(ierr);    							// Dirichlet boundary conditions x-dir, left side, u(0)=0
	//ierr = IGASetBoundaryValue(igaZ0,1, 0, 1, 0.0);CHKERRQ(ierr);    							// Dirichlet boundary conditions x-dir, left side, u(1)=0 //This means clamped in left side

	//ierr = IGASetBoundaryForm(iga,dir,side,PETSC_TRUE);CHKERRQ(ierr);								// Syntax	
	//ierr = IGASetBoundaryForm(igaZ0,1,1,PETSC_TRUE);CHKERRQ(ierr);								// Neumann boundary conditions right side

	// Get local vectors Chi0 and arrays
	ierr = IGAGetLocalVecArray(igaChi,chi0,&localChi0Z0,&arrayChi0Z0);CHKERRQ(ierr);

	// Element loop
	ierr = IGABeginElement(igaZ0,&elemZ0);CHKERRQ(ierr);
	ierr = IGABeginElement(igaChi,&elemChi);CHKERRQ(ierr);
	while (IGANextElement(igaZ0,elemZ0)) 
	{
		IGANextElement(igaChi,elemChi);
		ierr = IGAElementGetWorkMat(elemZ0,&KlocZ0);CHKERRQ(ierr);
		ierr = IGAElementGetWorkVec(elemZ0,&FlocZ0);CHKERRQ(ierr);
		ierr = IGAElementGetValues(elemChi,arrayChi0Z0,&Chi0Z0);CHKERRQ(ierr);

		// FormSystem loop
		while (IGAElementNextFormSystem(elemZ0,&wtfZ0,&wtf2Z0)) 
		{
		// Quadrature loop
			ierr = IGAElementBeginPoint(elemZ0,&pointZ0);CHKERRQ(ierr);
			ierr = IGAElementBeginPoint(elemChi,&pointChi);CHKERRQ(ierr);

			while (IGAElementNextPoint(elemZ0,pointZ0))
			{
				if(pointZ0->atboundary==1)
				{
					ierr = IGAPointGetWorkMat(pointZ0,&KpointZ0);CHKERRQ(ierr);
					ierr = IGAPointGetWorkVec(pointZ0,&FpointZ0);CHKERRQ(ierr);
					ierr = Z0sys(pointZ0,pointChi,KpointZ0,FpointZ0,Chi0Z0,NULL);CHKERRQ(ierr);
					ierr = IGAPointAddMat(pointZ0,KpointZ0,KlocZ0);CHKERRQ(ierr);
					ierr = IGAPointAddVec(pointZ0,FpointZ0,FlocZ0);CHKERRQ(ierr);
				}

				if(pointZ0->atboundary==0 && pointChi->atboundary==0)
				{
					IGAElementNextPoint(elemChi,pointChi);
					ierr = IGAPointGetWorkMat(pointZ0,&KpointZ0);CHKERRQ(ierr);
					ierr = IGAPointGetWorkVec(pointZ0,&FpointZ0);CHKERRQ(ierr);
					ierr = Z0sys(pointZ0,pointChi,KpointZ0,FpointZ0,Chi0Z0,NULL);CHKERRQ(ierr);
					ierr = IGAPointAddMat(pointZ0,KpointZ0,KlocZ0);CHKERRQ(ierr);
					ierr = IGAPointAddVec(pointZ0,FpointZ0,FlocZ0);CHKERRQ(ierr);
				}
			}
			if (pointChi->index != -1)
			{
				IGAElementNextPoint(elemChi,pointChi);
			}
			ierr = IGAElementEndPoint(elemZ0,&pointZ0);CHKERRQ(ierr);
			ierr = IGAElementEndPoint(elemChi,&pointChi);CHKERRQ(ierr);
		}

		ierr = IGAElementFixSystem(elemZ0,KlocZ0,FlocZ0);CHKERRQ(ierr);					//This sets Dirichlet condition ¿? (Yes, this applies the conditions from IGASetBoundaryValue)
		ierr = IGAElementAssembleMat(elemZ0,KlocZ0,KZ0);CHKERRQ(ierr);
		ierr = IGAElementAssembleVec(elemZ0,FlocZ0,FZ0);CHKERRQ(ierr);

	}
	IGANextElement(igaChi,elemChi);
	ierr = IGAEndElement(igaZ0,&elemZ0);CHKERRQ(ierr);
	ierr = IGAEndElement(igaChi,&elemChi);CHKERRQ(ierr);

	// Restore local vectors Chi0 and arrays
	ierr = IGARestoreLocalVecArray(igaChi,chi0,&localChi0Z0,&arrayChi0Z0);CHKERRQ(ierr);

	ierr = MatAssemblyBegin(KZ0,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyEnd  (KZ0,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

	ierr = VecAssemblyBegin(FZ0);CHKERRQ(ierr);
	ierr = VecAssemblyEnd  (FZ0);CHKERRQ(ierr);


	//Prueba cond borde derivada (Cambiar a un punto fijo despues)
	//for (i=2*(nx+1); i<4*(nx+1); i=i+1)
	//{
		//rows=i;
		//ierr = MatZeroRows(KZ0,1,&rows,1.0e16,0,0);
		//ierr = VecSetValue(FZ0,rows,0.0,INSERT_VALUES);CHKERRQ(ierr);
	//}
	rows=0;
	ierr = MatZeroRows(KZ0,1,&rows,1.0e30,0,0);CHKERRQ(ierr);
	rows=1;
	ierr = MatZeroRows(KZ0,1,&rows,1.0e30,0,0);CHKERRQ(ierr);
	rows=2*(nx+2);
	ierr = MatZeroRows(KZ0,1,&rows,1.0e30,0,0);CHKERRQ(ierr);

	ierr = MatAssemblyBegin(KZ0,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyEnd  (KZ0,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

	//Here we set values to the vector directly (to impose Dirichlet condition in a single point)
	rows=0;
	vals=0.0;
	ierr = VecSetValue(FZ0,rows,vals,INSERT_VALUES);CHKERRQ(ierr);
	rows=1;
	vals=0.0;
	ierr = VecSetValue(FZ0,rows,vals,INSERT_VALUES);CHKERRQ(ierr);
	rows=2*(nx+1);
	vals=0.0;
	ierr = VecSetValue(FZ0,rows,vals,INSERT_VALUES);CHKERRQ(ierr);
	ierr = VecAssemblyBegin(FZ0);CHKERRQ(ierr);
	ierr = VecAssemblyEnd  (FZ0);CHKERRQ(ierr);

	//Hasta aqui

	ierr = KSPSetOperators(kspZ0,KZ0,KZ0);CHKERRQ(ierr);
	ierr = KSPSetFromOptions(kspZ0);CHKERRQ(ierr);
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

//System for equilibrium
	IGA igaEq;
	ierr = IGACreate(PETSC_COMM_WORLD,&igaEq);CHKERRQ(ierr);
	ierr = IGASetDim(igaEq,2);CHKERRQ(ierr);													//Spatial dimension of the problem
	ierr = IGASetDof(igaEq,2);CHKERRQ(ierr);													//Number of degrees of freedom, per node
	ierr = IGASetOrder(igaEq,2);CHKERRQ(ierr);													//Number of spatial derivatives to calculate
	ierr = IGASetFromOptions(igaEq);CHKERRQ(ierr);												//Note: The order (or degree) of the shape functions is given by the mesh!
	ierr = IGARead(igaEq,"./geometry2.dat");CHKERRQ(ierr);
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
	const PetscReal	*arrayZ0Eq, *arrayChi0Eq;	//arrayU
	Vec				localZ0Eq,localChi0Eq;		//localU
	PetscReal		*Z0Eq,*Chi0Eq;				//U

  	IGAFormSystem	wtfEq;
 	void			*wtf2Eq;

 	KSP kspEq;
	ierr = IGACreateKSP(igaEq,&kspEq);CHKERRQ(ierr);

	//ierr = IGASetBoundaryValue(iga,dir,side,dof,value);CHKERRQ(ierr);							//Syntax
	ierr = IGASetBoundaryValue(igaEq,1 ,0 ,0 ,0.0);CHKERRQ(ierr);    							// Dirichlet boundary conditions x-dir, left side, u(0)=0
	ierr = IGASetBoundaryValue(igaEq,1, 0, 1, 0.0);CHKERRQ(ierr);    							// Dirichlet boundary conditions x-dir, left side, u(1)=0 //This means clamped in left side

	//ierr = IGASetBoundaryForm(iga,dir,side,PETSC_TRUE);CHKERRQ(ierr);							// Syntax	
	ierr = IGASetBoundaryForm(igaEq,1,1,PETSC_TRUE);CHKERRQ(ierr);								// Neumann boundary conditions top side

	// Get local vectors Z0  and Chi0 and arrays
	ierr = IGAGetLocalVecArray(igaChi,chi0,&localChi0Eq,&arrayChi0Eq);CHKERRQ(ierr);
	ierr = IGAGetLocalVecArray(igaZ0,Z0,&localZ0Eq,&arrayZ0Eq);CHKERRQ(ierr);

	// Element loop
	ierr = IGABeginElement(igaEq,&elemEq);CHKERRQ(ierr);
	ierr = IGABeginElement(igaZ0,&elemZ0);CHKERRQ(ierr);
	ierr = IGABeginElement(igaChi,&elemChi);CHKERRQ(ierr);
	while (IGANextElement(igaEq,elemEq)) 
	{
		IGANextElement(igaZ0,elemZ0);
		IGANextElement(igaChi,elemChi);
		ierr = IGAElementGetWorkMat(elemEq,&KlocEq);CHKERRQ(ierr);
		ierr = IGAElementGetWorkVec(elemEq,&FlocEq);CHKERRQ(ierr);
		ierr = IGAElementGetValues(elemZ0,arrayZ0Eq,&Z0Eq);CHKERRQ(ierr);
		ierr = IGAElementGetValues(elemChi,arrayChi0Eq,&Chi0Eq);CHKERRQ(ierr);

		// FormSystem loop
		while (IGAElementNextFormSystem(elemEq,&wtfEq,&wtf2Eq)) 
		{
		// Quadrature loop
			ierr = IGAElementBeginPoint(elemEq,&pointEq);CHKERRQ(ierr);
			ierr = IGAElementBeginPoint(elemZ0,&pointZ0);CHKERRQ(ierr);
			ierr = IGAElementBeginPoint(elemChi,&pointChi);CHKERRQ(ierr);

			while (IGAElementNextPoint(elemEq,pointEq))
			{
				if(pointEq->atboundary==1)
				{
					ierr = IGAPointGetWorkMat(pointEq,&KpointEq);CHKERRQ(ierr);
					ierr = IGAPointGetWorkVec(pointEq,&FpointEq);CHKERRQ(ierr);
					ierr = Eq(pointEq,pointZ0,pointChi,KpointEq,FpointEq,Z0Eq,Chi0Eq,NULL);CHKERRQ(ierr);
					ierr = IGAPointAddMat(pointEq,KpointEq,KlocEq);CHKERRQ(ierr);
					ierr = IGAPointAddVec(pointEq,FpointEq,FlocEq);CHKERRQ(ierr);
				}

				if(pointEq->atboundary==0 && pointZ0->atboundary==0)
				{
					IGAElementNextPoint(elemZ0,pointZ0);
					IGAElementNextPoint(elemChi,pointChi);
					ierr = IGAPointGetWorkMat(pointEq,&KpointEq);CHKERRQ(ierr);
					ierr = IGAPointGetWorkVec(pointEq,&FpointEq);CHKERRQ(ierr);
					ierr = Eq(pointEq,pointZ0,pointChi,KpointEq,FpointEq,Z0Eq,Chi0Eq,NULL);CHKERRQ(ierr);
					ierr = IGAPointAddMat(pointEq,KpointEq,KlocEq);CHKERRQ(ierr);
					ierr = IGAPointAddVec(pointEq,FpointEq,FlocEq);CHKERRQ(ierr);
				}
			}
			if (pointZ0->index != -1)
			{
				IGAElementNextPoint(elemZ0,pointZ0);
			}
			if (pointChi->index != -1)
			{
				IGAElementNextPoint(elemChi,pointChi);
			}
			ierr = IGAElementEndPoint(elemEq,&pointEq);CHKERRQ(ierr);
			ierr = IGAElementEndPoint(elemZ0,&pointZ0);CHKERRQ(ierr);
			ierr = IGAElementEndPoint(elemChi,&pointChi);CHKERRQ(ierr);
		}

		ierr = IGAElementFixSystem(elemEq,KlocEq,FlocEq);CHKERRQ(ierr);					//This sets Dirichlet condition ¿? (Yes, this applies the conditions from IGASetBoundaryValue)
		ierr = IGAElementAssembleMat(elemEq,KlocEq,KEq);CHKERRQ(ierr);
		ierr = IGAElementAssembleVec(elemEq,FlocEq,FEq);CHKERRQ(ierr);

	}
	IGANextElement(igaZ0,elemZ0);
	IGANextElement(igaChi,elemChi);
	ierr = IGAEndElement(igaEq,&elemEq);CHKERRQ(ierr);
	ierr = IGAEndElement(igaZ0,&elemZ0);CHKERRQ(ierr);
	ierr = IGAEndElement(igaChi,&elemChi);CHKERRQ(ierr);

	// Restore local vectors Z0, Chi0 and arrays
	ierr = IGARestoreLocalVecArray(igaZ0,Z0,&localZ0Eq,&arrayZ0Eq);CHKERRQ(ierr);
	ierr = IGARestoreLocalVecArray(igaChi,chi0,&localChi0Eq,&arrayChi0Eq);CHKERRQ(ierr);

	ierr = MatAssemblyBegin(KEq,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyEnd  (KEq,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

	ierr = VecAssemblyBegin(FEq);CHKERRQ(ierr);
	ierr = VecAssemblyEnd  (FEq);CHKERRQ(ierr);


	//Bondary condition for derivative
	for (i=2*(nx+2); i<4*(nx+2); i=i+2)
	{
		rows=i;
		ierr = MatZeroRows(KEq,1,&rows,1.0e30,0,0);
		ierr = VecSetValue(FEq,rows,0.0,INSERT_VALUES);CHKERRQ(ierr);
	}

	ierr = MatAssemblyBegin(KEq,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyEnd  (KEq,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = VecAssemblyBegin(FEq);CHKERRQ(ierr);
	ierr = VecAssemblyEnd  (FEq);CHKERRQ(ierr);

	//Hasta aqui

	ierr = KSPSetOperators(kspEq,KEq,KEq);CHKERRQ(ierr);
	ierr = KSPSetFromOptions(kspEq);CHKERRQ(ierr);
	ierr = KSPSolve(kspEq,FEq,u0);CHKERRQ(ierr);

	ierr = KSPDestroy(&kspEq);CHKERRQ(ierr);
	ierr = MatDestroy(&KEq);CHKERRQ(ierr);
	ierr = VecDestroy(&FEq);CHKERRQ(ierr);
	//ierr = IGAWriteVec(igaEq,u0,"./results/u-2d-0.dat");CHKERRQ(ierr);
	char nameEq[]="/u-2d-0.dat";
	char pathEq[512];
	sprintf(pathEq,"%s%s",direct,nameEq);
	ierr = IGAWriteVec(igaEq,u0,pathEq);CHKERRQ(ierr);
	
//

//System for L2 projection of stress
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

	IGAPoint		pointStress;										//point
	IGAElement		elemStress;											//element
	PetscReal		*KlocStress,*FlocStress;							//AA y BB
	PetscReal		*KpointStress,*FpointStress;						//KKK y FFF
	const PetscReal	*arrayUStress,*arrayZ0Stress,*arrayChi0Stress;		//arrayU
	Vec				localUStress,localZ0Stress,localChi0Stress;			//localU
	PetscReal		*UStress,*Z0Stress,*Chi0Stress;						//U

  	IGAFormSystem	wtfStress;
 	void			*wtf2Stress;

 	KSP kspStress;
	ierr = IGACreateKSP(igaStress,&kspStress);CHKERRQ(ierr);

	// Get local vectors u, Z0 and Chi0 and arrays
	ierr = IGAGetLocalVecArray(igaEq,u0,&localUStress,&arrayUStress);CHKERRQ(ierr);
	ierr = IGAGetLocalVecArray(igaChi,chi0,&localChi0Stress,&arrayChi0Stress);CHKERRQ(ierr);
	ierr = IGAGetLocalVecArray(igaZ0,Z0,&localZ0Stress,&arrayZ0Stress);CHKERRQ(ierr);

	// Element loop
	ierr = IGABeginElement(igaStress,&elemStress);CHKERRQ(ierr);
	ierr = IGABeginElement(igaEq,&elemEq);CHKERRQ(ierr);
	ierr = IGABeginElement(igaZ0,&elemZ0);CHKERRQ(ierr);
	ierr = IGABeginElement(igaChi,&elemChi);CHKERRQ(ierr);
	while (IGANextElement(igaStress,elemStress)) 
	{
		IGANextElement(igaEq,elemEq);
		IGANextElement(igaZ0,elemZ0);
		IGANextElement(igaChi,elemChi);
		ierr = IGAElementGetWorkMat(elemStress,&KlocStress);CHKERRQ(ierr);
		ierr = IGAElementGetWorkVec(elemStress,&FlocStress);CHKERRQ(ierr);
		ierr = IGAElementGetValues(elemEq,arrayUStress,&UStress);CHKERRQ(ierr);
		ierr = IGAElementGetValues(elemZ0,arrayZ0Stress,&Z0Stress);CHKERRQ(ierr);
		ierr = IGAElementGetValues(elemChi,arrayChi0Stress,&Chi0Stress);CHKERRQ(ierr);

		// FormSystem loop
		while (IGAElementNextFormSystem(elemStress,&wtfStress,&wtf2Stress)) 
		{
		// Quadrature loop
			ierr = IGAElementBeginPoint(elemStress,&pointStress);CHKERRQ(ierr);
			ierr = IGAElementBeginPoint(elemEq,&pointEq);CHKERRQ(ierr);
			ierr = IGAElementBeginPoint(elemZ0,&pointZ0);CHKERRQ(ierr);
			ierr = IGAElementBeginPoint(elemChi,&pointChi);CHKERRQ(ierr);

			while (IGAElementNextPoint(elemStress,pointStress))
			{
				if(pointStress->atboundary==1)
				{
					ierr = IGAPointGetWorkMat(pointStress,&KpointStress);CHKERRQ(ierr);
					ierr = IGAPointGetWorkVec(pointStress,&FpointStress);CHKERRQ(ierr);
					//Stress(IGAPoint p,IGAPoint pU, IGAPoint pZ,IGAPoint pChi,PetscReal *K,PetscReal *F,PetscReal *U,PetscReal *Z, PetscReal *Chi,void *ctx)
					ierr = Stress(pointStress,pointEq,pointZ0,pointChi,KpointStress,FpointStress,UStress,Z0Stress,Chi0Stress,NULL);CHKERRQ(ierr);
					ierr = IGAPointAddMat(pointStress,KpointStress,KlocStress);CHKERRQ(ierr);
					ierr = IGAPointAddVec(pointStress,FpointStress,FlocStress);CHKERRQ(ierr);
				}

				if(pointStress->atboundary==0 && pointEq->atboundary==0 && pointZ0->atboundary==0 && pointChi->atboundary==0)
				{
					IGAElementNextPoint(elemEq,pointEq);
					IGAElementNextPoint(elemZ0,pointZ0);
					IGAElementNextPoint(elemChi,pointChi);
					ierr = IGAPointGetWorkMat(pointStress,&KpointStress);CHKERRQ(ierr);
					ierr = IGAPointGetWorkVec(pointStress,&FpointStress);CHKERRQ(ierr);
					//Stress(IGAPoint p,IGAPoint pU, IGAPoint pZ,IGAPoint pChi,PetscReal *K,PetscReal *F,PetscReal *U,PetscReal *Z, PetscReal *Chi,void *ctx)
					ierr = Stress(pointStress,pointEq,pointZ0,pointChi,KpointStress,FpointStress,UStress,Z0Stress,Chi0Stress,NULL);CHKERRQ(ierr);
					ierr = IGAPointAddMat(pointStress,KpointStress,KlocStress);CHKERRQ(ierr);
					ierr = IGAPointAddVec(pointStress,FpointStress,FlocStress);CHKERRQ(ierr);
				}
			}
			if (pointEq->index != -1)
			{
				IGAElementNextPoint(elemEq,pointEq);
			}
			if (pointZ0->index != -1)
			{
				IGAElementNextPoint(elemZ0,pointZ0);
			}
			if (pointChi->index != -1)
			{
				IGAElementNextPoint(elemChi,pointChi);
			}
			ierr = IGAElementEndPoint(elemStress,&pointStress);CHKERRQ(ierr);
			ierr = IGAElementEndPoint(elemEq,&pointEq);CHKERRQ(ierr);
			ierr = IGAElementEndPoint(elemZ0,&pointZ0);CHKERRQ(ierr);
			ierr = IGAElementEndPoint(elemChi,&pointChi);CHKERRQ(ierr);
		}

		ierr = IGAElementFixSystem(elemStress,KlocStress,FlocStress);CHKERRQ(ierr);					//This sets Dirichlet condition ¿? (Yes, this applies the conditions from IGASetBoundaryValue)
		ierr = IGAElementAssembleMat(elemStress,KlocStress,KStress);CHKERRQ(ierr);
		ierr = IGAElementAssembleVec(elemStress,FlocStress,FStress);CHKERRQ(ierr);
	}
	IGANextElement(igaEq,elemEq);
	IGANextElement(igaZ0,elemZ0);
	IGANextElement(igaChi,elemChi);

	ierr = IGAEndElement(igaStress,&elemStress);CHKERRQ(ierr);
	ierr = IGAEndElement(igaEq,&elemEq);CHKERRQ(ierr);
	ierr = IGAEndElement(igaZ0,&elemZ0);CHKERRQ(ierr);
	ierr = IGAEndElement(igaChi,&elemChi);CHKERRQ(ierr);

	// Restore local vectors u, Z0, Chi0 and arrays
	ierr = IGARestoreLocalVecArray(igaEq,u0,&localUStress,&arrayUStress);CHKERRQ(ierr);
	ierr = IGARestoreLocalVecArray(igaZ0,Z0,&localZ0Stress,&arrayZ0Stress);CHKERRQ(ierr);
	ierr = IGARestoreLocalVecArray(igaChi,chi0,&localChi0Stress,&arrayChi0Stress);CHKERRQ(ierr);

	ierr = MatAssemblyBegin(KStress,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyEnd  (KStress,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

	ierr = VecAssemblyBegin(FStress);CHKERRQ(ierr);
	ierr = VecAssemblyEnd  (FStress);CHKERRQ(ierr);

	ierr = KSPSetOperators(kspStress,KStress,KStress);CHKERRQ(ierr);
	ierr = KSPSetFromOptions(kspStress);CHKERRQ(ierr);
	ierr = KSPSolve(kspStress,FStress,sigma0);CHKERRQ(ierr);

	ierr = KSPDestroy(&kspStress);CHKERRQ(ierr);
	ierr = MatDestroy(&KStress);CHKERRQ(ierr);
	ierr = VecDestroy(&FStress);CHKERRQ(ierr);

	char nameStress[]="/sigma-2d-0.dat";
	char pathStress[512];
	sprintf(pathStress,"%s%s",direct,nameStress);
	ierr = IGAWriteVec(igaStress,sigma0,pathStress);CHKERRQ(ierr);
	
//

//Destroy all objects not needed anymore (Better to do it here in case different codes call the same IGA, move if memory is a problem)
	ierr = IGADestroy(&igaAl);CHKERRQ(ierr);
	ierr = IGADestroy(&igaPi);CHKERRQ(ierr);
	ierr = IGADestroy(&igaZ);CHKERRQ(ierr);
	ierr = IGADestroy(&igaCS);CHKERRQ(ierr);
	ierr = IGADestroy(&igaChi);CHKERRQ(ierr);
	ierr = IGADestroy(&igaAl);CHKERRQ(ierr);
	ierr = IGADestroy(&igaZ0);CHKERRQ(ierr);
	ierr = IGADestroy(&igaEq);CHKERRQ(ierr);
	ierr = IGADestroy(&igaStress);CHKERRQ(ierr);

//

ierr = PetscFinalize();CHKERRQ(ierr);

return 0;
}