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

		PetscReal x[dim];																	//Vector de reales, tamaño igual a dimensión del problema
		IGAPointFormGeomMap(p,x);															//llena x con las coordenadas de p, punto de Gauss

		//g es la función a proyectar
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
				g[i]=1.917545*(0.5*(tanh((x[0]+10*dx)/eps)+1.0)-0.5*(tanh((x[0]+8*dx)/eps)+1.0))*(0.5*(tanh((x[1]+6*dy)/eps)+1.0)-0.5*(tanh((x[1]+4*dy)/eps)+1.0));;
			}
			else
			{
				g[i]=1.917545*(0.5*(tanh((x[0]-dx)/eps)+1.0)-0.5*(tanh((x[0]+dx)/eps)+1.0))*(0.5*(tanh((x[1]-7*dy)/eps)+1.0)-0.5*(tanh((x[1]-5*dy)/eps)+1.0));;;
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
			return 0;														//Condición de borde Neumann nula
		}

		AppCtxL2    *user  = (AppCtxL2 *)ctx;
		PetscReal Lx     = user->Lx;
		PetscReal Ly     = user->Ly;
		PetscReal nx     = user->nx;
		PetscReal ny     = user->ny;

		PetscInt a,b,i;
		PetscInt nen = p->nen;																//Numero de funciones de forma, en este caso es 9
		PetscInt dim = p->dim;																//Dimension del problema
		PetscInt dof = p->dof;

		PetscReal x[dim];																	//Vector de reales, tamaño igual a dimensión del problema
		IGAPointFormGeomMap(p,x);															//llena x con las coordenadas de p, punto de Gauss

		//g es la función a proyectar
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
				g[i]=1.917545*(0.5*(tanh((x[0]+10*dx)/eps)+1.0)-0.5*(tanh((x[0]+8*dx)/eps)+1.0))*(0.5*(tanh((x[1]+dy)/eps)+1.0)-0.5*(tanh((x[1]-dy)/eps)+1.0));;
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

//System for dPi/dt=curl(Pi)
	#undef  __FUNCT__
	#define __FUNCT__ "Pi"
	PetscErrorCode Pi(IGAPoint p,PetscReal *K,PetscReal *F, PetscReal *U,void *ctx)
	{
		AppCtx    *user  = (AppCtx *)ctx;
		PetscReal dt     = user->dt;

		const PetscReal *N0,(*N1)[2];
		IGAPointGetShapeFuns(p,0,(const PetscReal**)&N0);
		IGAPointGetShapeFuns(p,1,(const PetscReal**)&N1);
		PetscInt a,b,i,nen=p->nen, dof=p->dof;

		PetscReal Pi0[dof];																	//Aquí asignamos Pi(t) a un vector
		PetscReal d_Pi0[dof][2];															//Idem para las derivadas parciales
		IGAPointFormValue(p,U,&Pi0[0]);														
		IGAPointFormGrad (p,U,&d_Pi0[0][0]);

		PetscReal (*Kpi)[dof][nen][dof] = (typeof(Kpi)) K;
		PetscReal (*Fpi)[dof] = (PetscReal (*)[dof])F;

		PetscReal V1=1.0;
		PetscReal V2=0.0;
		PetscReal V1_x=0.0;
		PetscReal V2_y=0.0;

		PetscInt dim = p->dim;																//Dimension del problema
		PetscReal x[dim];
		IGAPointFormGeomMap(p,x);

		if (p->atboundary)
		{
			PetscInt dir  = p->boundary_id / 2;
			PetscInt side = p->boundary_id % 2;

			if(dir==0 && side==1)
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
										+ dt*(Nb*Na*(V1_x+V2_y))							//(LS7)  es igual a (LS3)
										+ dt*dt*((Nb_x*V1+Nb_y*V2)*(Na_x*V1+Na_x*V2))		//(LS10)
										+ dt*dt*((Nb_x*V1+Nb_y*V2)*Na*(V1_x+V2_y))			//(LS11)
										+ dt*dt*(Nb*(V1_x+V2_y)*(Na_x*V1+Na_y*V2))			//(LS14)
										+ dt*dt*(Nb*(V1_x+V2_y)*Na*(V1_x+V2_y))				//(LS15)
										+ Na*Nb												//(GL1) 
										- dt*(Nb*(Na_x*V1+Na_y*V2))							//(GL3)  igual pero negativo a (LS2)
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

		PetscReal Al0[dof];																	//Assign Pi(t) to a vector
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

			if(dir==1 && side==1)
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
										-dt*Nb*(Na_x*V1+Na_y*V2)								//(GL3) (negativo de (LS2))
										;		//Los demás terminos son 0 en 2d
					}
				}

				for (i=0; i<dof; i++)
				{
					Fal[a][i] = Al0[i]*Na														//(LS6)
								-dt*Na*(d_Al0[i][0]*V1+d_Al0[i][1]*V2)							//(LS7)
								-dt*Al0[i]*Na*(V1_x+V2_y)										//(LS8)
								+dt*Al0[i]*(Na_x*V1+Na_y*V2)									//(LS11)
								+dt*Al0[i]*Na*(V1_x+V2_y)										//(LS12) (negativo de (LS8))
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
		
		PetscReal alfa[3][3]={0};

		const PetscReal *N0,(*N1)[2];
		IGAPointGetShapeFuns(p,0,(const PetscReal**)&N0);
		IGAPointGetShapeFuns(p,1,(const PetscReal**)&N1);
		PetscInt a,b,u,w,i,j,k,l,nen=p->nen, dof=p->dof;

		PetscReal x[2];																			//Vector de reales, tamaño igual a dimensión del problema
		IGAPointFormGeomMap(p,x);																//llena x con las coordenadas de p, punto de Gauss

		PetscReal Al0[2];																		//Create array to recieve Alfa
		//PetscReal dAl0[2][2];																	//Create array to recieve dAlfa
		IGAPointFormValue(pAl,U,&Al0[0]);														//This fills the values
		//IGAPointFormGrad (pAl,U,&dAl0[0][0]);													//This fills the values of the derivatives

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

		PetscReal x[2];																			//Vector de reales, tamaño igual a dimensión del problema
		IGAPointFormGeomMap(p,x);																//llena x con las coordenadas de p, punto de Gauss

		PetscReal V[3]={0};		V[0]=x[0]; V[1]=x[1];
		PetscReal dV[3][3]={0}; dV[0][0]=1.0; dV[0][1]=1.0;

		//if(x[0]>=-0.25/101.0 && x[0]<=0.25/101.0 && x[1]>=-0.25/101.0 && x[1]<=0.25/101.0)
		//{
		//	Al0[0]=1.0; Al0[1]=1.0;
		//}
		//else
		//{
		//	Al0[0]=0.0; Al0[1]=0.0;
		//}

		alfa[0][2]=Al0[0];	alfa[1][2]=Al0[1];												//Assign Al0 to the full alfa expression
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

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc, char *argv[]) {
	PetscErrorCode  ierr;

//Delete files in result folder
	char direct1[]="./results";
	//char direct2[]="/home/eazegpi/sharedResult";
	char direct2[]="../../sharedResults";
	char *direct;

	DIR *pDir1=opendir(direct1);

	if (pDir1 != NULL)
	{
		printf("Habemus carpeta %s \n", direct1);
		closedir(pDir1);
		direct=direct1;
	}
	else
	{
		printf("No habemus carpeta %s \n", direct1);	
	}

	DIR *pDir2=opendir(direct2);

	if (pDir2 != NULL)
	{
		printf("Habemus carpeta %s \n", direct2);
		closedir(pDir2);
		direct=direct2;
	}
	else
	{
		printf("No habemus carpeta %s \n", direct2);	
	}

	printf("directorio %s\n", direct);

	DIR *folder = opendir(direct);
	struct dirent *next_file;
	char filepath[512];

	while ( (next_file = readdir(folder)) != NULL )
	{
		//sprintf(filepath, "%s/%s", "./results", next_file->d_name);
		sprintf(filepath, "%s/%s", direct, next_file->d_name);
		remove(filepath);
	}
	closedir(folder);
//

//Generate Mesh
	PetscReal Lx=0.5;
	PetscReal Ly=0.5;
	PetscInt  nx=101;
	PetscInt  ny=101;
	PetscReal h=fmin(Lx/nx,Ly/ny);

	printf("h= %f\n",h);

	char      nombreMalla[12]="geometry";
	char 	  comando[512];
	sprintf(comando,"exec python rectangle.py %s %f %f %d %d\n",nombreMalla,Lx,Ly,nx,ny);
	ierr=system(comando);

	//Time parameters
	PetscReal T=0.1;
	PetscReal fact=0.5;
	PetscReal dt=fact*(h/2.0);																	//dt must be less than h/2
	printf("dt<= %f",dt/fact);	printf("     dt= %f\n",dt);
	PetscInt numStep=ceil(T/dt);
	printf("Nt= %d\n",numStep);

//

//Creation of solution systems
	ierr = PetscInitialize(&argc,&argv,0,0);CHKERRQ(ierr);										//Always initialize PETSc

	//System for Pi 
	IGA igaPi;
	ierr = IGACreate(PETSC_COMM_WORLD,&igaPi);CHKERRQ(ierr);
	ierr = IGASetDim(igaPi,2);CHKERRQ(ierr);													//Spatial dimension of the problem
	ierr = IGASetDof(igaPi,4);CHKERRQ(ierr);													//Number of degrees of freedom, per node
	ierr = IGASetOrder(igaPi,1);CHKERRQ(ierr);													//Number of spatial derivatives to calculate
	ierr = IGASetFromOptions(igaPi);CHKERRQ(ierr);												//Nota: El orden (o grado) de las funciones de forma viene dado por la malla!
	ierr = IGARead(igaPi,"./geometry.dat");CHKERRQ(ierr);
	ierr = IGASetUp(igaPi);CHKERRQ(ierr);

	PetscInt dir,side;
	for (dir=0; dir<2; dir++) 
	{
		for (side=0; side<2; side++) 
		{
			//ierr = IGASetBoundaryValue(iga,dir,side,dof,0.0);CHKERRQ(ierr);    				// Dirichlet boundary conditions
			ierr = IGASetBoundaryForm(igaPi,dir,side,PETSC_TRUE);CHKERRQ(ierr);  				// Neumann boundary conditions
		}
	}

	//System for Alfa
	IGA igaAl;
	ierr = IGACreate(PETSC_COMM_WORLD,&igaAl);CHKERRQ(ierr);
	ierr = IGASetDim(igaAl,2);CHKERRQ(ierr);													//Dimnsion espacial del problema
	ierr = IGASetDof(igaAl,2);CHKERRQ(ierr);													//Number of degrees of freedom, per node
	ierr = IGASetOrder(igaAl,1);CHKERRQ(ierr);													//Number of spatial derivatives to calculate
	ierr = IGASetFromOptions(igaAl);CHKERRQ(ierr);												//Nota: El orden (o grado) de las funciones de forma viene dado por la malla!
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

	//System for Chi
	IGA igaChi;
	ierr = IGACreate(PETSC_COMM_WORLD,&igaChi);CHKERRQ(ierr);
	ierr = IGASetDim(igaChi,2);CHKERRQ(ierr);													//Spatial dimension of the problem
	ierr = IGASetDof(igaChi,4);CHKERRQ(ierr);													//Numero de grados de libertad
	ierr = IGASetOrder(igaChi,1);CHKERRQ(ierr);													//Numero máximo de derivadas a calcular
	ierr = IGASetFromOptions(igaChi);CHKERRQ(ierr);												//Nota: El orden (o grado) de las funciones de forma viene dado por la malla!
	ierr = IGARead(igaChi,"./geometry.dat");CHKERRQ(ierr);
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
	
	//System for z\dot
	IGA igaZ;
	ierr = IGACreate(PETSC_COMM_WORLD,&igaZ);CHKERRQ(ierr);
	ierr = IGASetDim(igaZ,2);CHKERRQ(ierr);														//Spatial dimension of the problem
	ierr = IGASetDof(igaZ,2);CHKERRQ(ierr);														//Number of degrees of freedom, per node
	ierr = IGASetOrder(igaZ,1);CHKERRQ(ierr);													//Numero máximo de derivadas a calcular
	ierr = IGASetFromOptions(igaZ);CHKERRQ(ierr);												//Nota: El orden (o grado) de las funciones de forma viene dado por la malla!
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
//

//Creation of types and systems for the L2 projection of Pi0
	AppCtxL2 userL2;
	userL2.Lx     = Lx;
	userL2.Ly     = Ly;
	userL2.nx     = nx;
	userL2.ny     = ny;
	//void* user;

	Vec pi0;
	Mat Kl2Pi;
	Vec Fl2Pi;
	ierr = IGACreateVec(igaPi,&pi0);CHKERRQ(ierr);  
	ierr = IGACreateMat(igaPi,&Kl2Pi);CHKERRQ(ierr);
	ierr = IGACreateVec(igaPi,&Fl2Pi);CHKERRQ(ierr);
	ierr = IGASetFormSystem(igaPi,L2ProjectionPi,&userL2);CHKERRQ(ierr);
	ierr = IGAComputeSystem(igaPi,Kl2Pi,Fl2Pi);CHKERRQ(ierr);

	//Esta parte llama a KSP para resolver el sistema lineal
	KSP kspl2;
	ierr = IGACreateKSP(igaPi,&kspl2);CHKERRQ(ierr);
	ierr = KSPSetOperators(kspl2,Kl2Pi,Kl2Pi);CHKERRQ(ierr); 								//Esta función crea la matriz a resolver en el segundo parámetro y usa el 3er parámetro como preacondicionador
	ierr = KSPSetType(kspl2,KSPCG);CHKERRQ(ierr);										//Usando KSPCG (gradiente conjugado) porque la matriz es simétrica
	ierr = KSPSetOptionsPrefix(kspl2,"l2pPi_");CHKERRQ(ierr);
	ierr = KSPSetFromOptions(kspl2);CHKERRQ(ierr);
	//ierr = KSPSetTolerances(kspl2,1e-10,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);CHKERRQ(ierr);
	ierr = KSPSolve(kspl2,Fl2Pi,pi0);CHKERRQ(ierr);										//This is a simple system, so it can be solved with just this command

	//ierr = KSPDestroy(&kspl2);CHKERRQ(ierr);
	ierr = MatDestroy(&Kl2Pi);CHKERRQ(ierr);
	ierr = VecDestroy(&Fl2Pi);CHKERRQ(ierr);
	//ierr = IGAWriteVec(igaPi,pi0,"./results/Pi-2d-0.dat");CHKERRQ(ierr);
	char namePi[]="/Pi-2d-0.dat";
	char pathPi[512];
	sprintf(pathPi,"%s%s",direct,namePi);
	ierr = IGAWriteVec(igaPi,pi0,pathPi);CHKERRQ(ierr);
//

//Creation of types and systems for the L2 projection of Alfa0
	Vec al0;
	Mat Kl2Al;
	Vec Fl2Al;
	ierr = IGACreateVec(igaAl,&al0);CHKERRQ(ierr);  
	ierr = IGACreateMat(igaAl,&Kl2Al);CHKERRQ(ierr);
	ierr = IGACreateVec(igaAl,&Fl2Al);CHKERRQ(ierr);
	ierr = IGASetFormSystem(igaAl,L2ProjectionAl,&userL2);CHKERRQ(ierr);
	ierr = IGAComputeSystem(igaAl,Kl2Al,Fl2Al);CHKERRQ(ierr);

	//Esta parte llama a KSP para resolver el sistema lineal
	ierr = IGACreateKSP(igaAl,&kspl2);CHKERRQ(ierr);
	ierr = KSPSetOperators(kspl2,Kl2Al,Kl2Al);CHKERRQ(ierr); 								//Esta función crea la matriz a resolver en el segundo parámetro y usa el 3er parámetro como preacondicionador
	ierr = KSPSetType(kspl2,KSPCG);CHKERRQ(ierr);										//Usando KSPCG (gradiente conjugado) porque la matriz es simétrica
	ierr = KSPSetOptionsPrefix(kspl2,"l2pAl_");CHKERRQ(ierr);
	ierr = KSPSetFromOptions(kspl2);CHKERRQ(ierr);
	//ierr = KSPSetTolerances(kspl2,1e-10,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);CHKERRQ(ierr);
	ierr = KSPSolve(kspl2,Fl2Al,al0);CHKERRQ(ierr);										//This is a simple system, so it can be solved with just this command

	ierr = KSPDestroy(&kspl2);CHKERRQ(ierr);
	ierr = MatDestroy(&Kl2Al);CHKERRQ(ierr);
	ierr = VecDestroy(&Fl2Al);CHKERRQ(ierr);
	//ierr = IGAWriteVec(igaAl,al0,"./results/Al-2d-0.dat");CHKERRQ(ierr);
	char nameAl[]="/Al-2d-0.dat";
	char pathAl[512];
	sprintf(pathAl,"%s%s",direct,nameAl);
	ierr = IGAWriteVec(igaAl,al0,pathAl);CHKERRQ(ierr);
//

//Creation of types and systems for the initial state of chi
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
			ierr = IGAElementEndPoint(elemAl,&pointAl);CHKERRQ(ierr);
		}
		ierr = IGAElementFixSystem(elemChi,KlocChi,FlocChi);CHKERRQ(ierr);					//Esto pone condicion de Dirichlet ¿?
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
	ierr = IGADestroy(&igaChi);CHKERRQ(ierr);

//Creation of types and systems for z
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
		ierr = IGAElementFixSystem(elemZ,KlocZ,FlocZ);CHKERRQ(ierr);					//Esto pone condicion de Dirichlet ¿?
		ierr = IGAElementAssembleMat(elemZ,KlocZ,Kz);CHKERRQ(ierr);
		ierr = IGAElementAssembleVec(elemZ,FlocZ,Fz);CHKERRQ(ierr);

	}
	IGANextElement(igaAl,elemAl);
	ierr = IGAEndElement(igaZ,&elemZ);CHKERRQ(ierr);
	ierr = IGAEndElement(igaAl,&elemAl);CHKERRQ(ierr);

	// Restore local vectors Al0 and arrays
	ierr = IGARestoreLocalVecArray(igaAl,al0,&localAl0Z,&arrayAl0Z);CHKERRQ(ierr);

	ierr = MatAssemblyBegin(Kz,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyEnd  (Kz,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	//Aquí meto lo que me mandó Hugo
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

	//Hasta aquí
	ierr = VecAssemblyBegin(Fz);CHKERRQ(ierr);
	ierr = VecAssemblyEnd  (Fz);CHKERRQ(ierr);
	//Aquí hay que meter lo que mandó hugo

	rows=0;
	vals=0.0;
	ierr = VecSetValue(Fz,rows,vals,INSERT_VALUES);CHKERRQ(ierr);
	rows=1;
	vals=0.0;
	ierr = VecSetValue(Fz,rows,vals,INSERT_VALUES);CHKERRQ(ierr);
	ierr = VecAssemblyBegin(Fz);CHKERRQ(ierr);
	ierr = VecAssemblyEnd  (Fz);CHKERRQ(ierr);

	//Hasta aquí
	ierr = KSPSetOperators(kspZ,Kz,Kz);CHKERRQ(ierr);
	ierr = KSPSetFromOptions(kspZ);CHKERRQ(ierr);
	ierr = KSPSolve(kspZ,Fz,z0);CHKERRQ(ierr);

	ierr = KSPDestroy(&kspZ);CHKERRQ(ierr);
	ierr = MatDestroy(&Kz);CHKERRQ(ierr);
	ierr = VecDestroy(&Fz);CHKERRQ(ierr);
	//ierr = IGAWriteVec(igaZ,z0,"./results/Z-2d-0.dat");CHKERRQ(ierr);
	char nameZ[]="/Z-2d-0.dat";
	char pathZ[512];
	sprintf(pathZ,"%s%s",direct,nameZ);
	ierr = IGAWriteVec(igaZ,z0,pathZ);CHKERRQ(ierr);
	ierr = IGADestroy(&igaZ);CHKERRQ(ierr);
//

//Creación de tipos y sistemas para Pi
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

//Creación de tipos y sistemas para Alfa
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


//Ciclo de solucion
	PetscInt i;
	for (i=0; i<numStep; i++) 
	{
		PetscPrintf(PETSC_COMM_WORLD,"Iter= %d \n",(i+1));
		//Reinicializar las varables en cada ciclo
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
			ierr = IGAElementFixSystem(elemPi,KlocPi,FlocPi);CHKERRQ(ierr);					//Esto pone condicion de Dirichlet ¿?
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

		// Solution cycle for Alfa
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
			ierr = IGAElementFixSystem(elemAl,KlocAl,FlocAl);CHKERRQ(ierr);					//Esto pone condicion de Dirichlet ¿?
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
			sprintf(nombrePi, "/Pi-2d-%d.dat", (i+1));											//Mete el numero al string
			char pathPi2[512];
			sprintf(pathPi2,"%s%s",direct,nombrePi);
			ierr = IGAWriteVec(igaPi,xPi,pathPi2);CHKERRQ(ierr);								//Graba el resultado al disco
			
			sprintf(nombreAl, "/Al-2d-%d.dat", (i+1));										//Mete el numero al string
			char pathAl2[512];
			sprintf(pathAl2,"%s%s",direct,nombreAl);
			ierr = IGAWriteVec(igaAl,xAl,pathAl2);CHKERRQ(ierr);
		}
		ierr=VecCopy(xPi,pi0); CHKERRQ(ierr);													//Aquí copiamos Pi(t+1) a Pi(t)
		ierr=VecCopy(xAl,al0); CHKERRQ(ierr);													//Aquí copiamos Al(t+1) a Al(t)
	}

	//ierr = VecView(x,PETSC_VIEWER_DRAW_WORLD);CHKERRQ(ierr);

	ierr = KSPDestroy(&kspPi);CHKERRQ(ierr);
	ierr = MatDestroy(&KPi);CHKERRQ(ierr);
	ierr = VecDestroy(&xPi);CHKERRQ(ierr);
	ierr = VecDestroy(&FPi);CHKERRQ(ierr);
	ierr = IGADestroy(&igaPi);CHKERRQ(ierr);

	ierr = KSPDestroy(&kspAl);CHKERRQ(ierr);
	ierr = MatDestroy(&KAl);CHKERRQ(ierr);
	ierr = VecDestroy(&xAl);CHKERRQ(ierr);
	ierr = VecDestroy(&FAl);CHKERRQ(ierr);
	ierr = IGADestroy(&igaAl);CHKERRQ(ierr);

	ierr = PetscFinalize();CHKERRQ(ierr);
//

return 0;
}