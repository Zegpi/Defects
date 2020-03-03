#include "petiga.h"
#include <stdio.h>
#include <dirent.h>
#include <math.h>

#if PETSC_VERSION_LT(3,5,0)
#define KSPSetOperators(ksp,A,B) KSPSetOperators(ksp,A,B,SAME_NONZERO_PATTERN)
#endif

//Definicion de AppCtx
typedef struct {
  PetscReal Lx;
  PetscReal Ly;
  PetscReal nx;
  PetscReal ny;
} AppCtxL2;

typedef struct {
  PetscReal dt;
} AppCtx;

PETSC_STATIC_INLINE
PetscBool IGAElementNextFormSystem(IGAElement element,IGAFormSystem *sys,void **ctx)
{
  IGAForm form = element->parent->form;
  if (!IGAElementNextForm(element,form->visit)) return PETSC_FALSE;
  *sys = form->ops->System;
  *ctx = form->ops->SysCtx;
  return PETSC_TRUE;
}


// Proyección L2 de Pi(0)
#undef  __FUNCT__
#define __FUNCT__ "L2Projection"
PetscErrorCode L2Projection(IGAPoint p,PetscReal *K,PetscReal *F,void *ctx)
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
			g[i]=0;
		}
		else if (i==2)
		{
			g[i]=1.917545*(0.5*(tanh((x[0]+10*dx)/eps)+1.0)-0.5*(tanh((x[0]+8*dx)/eps)+1.0))*(0.5*(tanh((x[1]+dy)/eps)+1.0)-0.5*(tanh((x[1]-dy)/eps)+1.0));;
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
	PetscReal d_Pi0[dof][2];																//Idem para las derivadas parciales
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

		PetscPrintf(PETSC_COMM_WORLD,"En el borde p->nen: %d \n", p->nen);

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
						Kpi[a][i][b][i] = Na*Nb*V1;
					}
				}
			}
		}
	}
	else
	{

		PetscPrintf(PETSC_COMM_WORLD,"En el cuerpo p->nen: %d \n", p->nen);

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
				/*
				Kpi[a][0][b][0] = Na*Nb														//(LS1)
								+ dt*(Nb*(Na_x*V1+Na_y*V2))									//(LS2)
								+ dt*(Nb*Na*(V1_x+V2_y))									//(LS3)
								+ dt*(Na*(Nb_x*V1+Nb_y*V2))									//(LS6)
								+ dt*(Nb*Na*(V1_x+V2_y))									//(LS7)  es igual a (LS3)
								+ dt*dt*((Nb_x*V1+Nb_y*V2)*(Na_x*V1+Na_x*V2))				//(LS10)
								+ dt*dt*((Nb_x*V1+Nb_y*V2)*Na*(V1_x+V2_y))					//(LS11)
								+ dt*dt*(Nb*(V1_x+V2_y)*(Na_x*V1+Na_y*V2))					//(LS14)
								+ dt*dt*(Nb*(V1_x+V2_y)*Na*(V1_x+V2_y))						//(LS15)
								+ Na*Nb														//(GL1) 
								- dt*(Nb*(Na_x*V1+Na_y*V2))									//(GL3)  igual pero negativo a (LS2)
								;		//Los demás terminos son 0 en 2d
				*/
				for (i=0; i<dof; i++)
				{
					Kpi[a][i][b][i] = Na*Nb													//(LS1)
									+ dt*(Nb*(Na_x*V1+Na_y*V2))								//(LS2)
									+ dt*(Nb*Na*(V1_x+V2_y))								//(LS3)
									+ dt*(Na*(Nb_x*V1+Nb_y*V2))								//(LS6)
									+ dt*(Nb*Na*(V1_x+V2_y))								//(LS7)  es igual a (LS3)
									+ dt*dt*((Nb_x*V1+Nb_y*V2)*(Na_x*V1+Na_x*V2))			//(LS10)
									+ dt*dt*((Nb_x*V1+Nb_y*V2)*Na*(V1_x+V2_y))				//(LS11)
									+ dt*dt*(Nb*(V1_x+V2_y)*(Na_x*V1+Na_y*V2))				//(LS14)
									+ dt*dt*(Nb*(V1_x+V2_y)*Na*(V1_x+V2_y))					//(LS15)
									+ Na*Nb													//(GL1) 
									- dt*(Nb*(Na_x*V1+Na_y*V2))								//(GL3)  igual pero negativo a (LS2)
									;		//Los demás terminos son 0 en 2d
				}
			}
			/*
			F[a] = 	Pi0[0]*Na 																//(LS26)
					+dt*(Pi0[0]*(Na_x*V1+Na_y*V2))											//(LS27)
					+dt*(Pi0[0]*Na*(V1_x+V2_y))												//(LS28)
					+Pi0[0]*Na 																//(GL2)
					;
			*/
			for (i=0; i<dof; i++)
			{
				Fpi[a][i] = Pi0[i]*Na 														//(LS26)
							+dt*(Pi0[i]*(Na_x*V1+Na_y*V2))									//(LS27)
							+dt*(Pi0[i]*Na*(V1_x+V2_y))										//(LS28)
							+Pi0[i]*Na 														//(GL2)
							;
			}
		}
	}
	//PetscPrintf(PETSC_COMM_WORLD,"p->nen: %d \n", p->nen);
	
	return 0;
}
//

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc, char *argv[]) {

//Borra resultados anteriores
	DIR *folder = opendir("/home/eazegpi/CodigosPetIGA/poisson2dmod/results");
	struct dirent *next_file;
	char filepath[256];

	while ( (next_file = readdir(folder)) != NULL )
	{
		sprintf(filepath, "%s/%s", "/home/eazegpi/CodigosPetIGA/poisson2dmod/results", next_file->d_name);
		remove(filepath);
	}
	closedir(folder);
//

//Generación de Malla
	PetscReal Lx=0.5;
	PetscReal Ly=0.5;
	PetscInt  nx=101;
	PetscInt  ny=101;
	char      nombreMalla[12]="geometry";
	char 	  comando[512];

	sprintf(comando,"exec python rectangle.py %s %f %f %d %d\n",nombreMalla,Lx,Ly,nx,ny);
	system(comando);
//


//Creación de tipos y sistemas para proyección L2
	AppCtxL2 userL2;
	userL2.Lx     = Lx;
	userL2.Ly     = Ly;
	userL2.nx     = nx;
	userL2.ny     = ny;
	//void* user;
	PetscErrorCode  ierr;
	ierr = PetscInitialize(&argc,&argv,0,0);CHKERRQ(ierr);								//Siempre inicializar petsc

	IGA iga;	//Aqui se define los parámetros para IGA
	ierr = IGACreate(PETSC_COMM_WORLD,&iga);CHKERRQ(ierr);
	ierr = IGASetDim(iga,2);CHKERRQ(ierr);												//Dimnsion espacial del problema
	ierr = IGASetDof(iga,4);CHKERRQ(ierr);												//Numero de grados de libertad
	ierr = IGASetOrder(iga,1);CHKERRQ(ierr);											//Numero máximo de derivadas a calcular
	ierr = IGASetFromOptions(iga);CHKERRQ(ierr);										//Nota: El orden (o grado) de las funciones de forma viene dado por la malla!
	ierr = IGARead(iga,"./geometry.dat");CHKERRQ(ierr);
	ierr = IGASetUp(iga);CHKERRQ(ierr);

    PetscInt dir,side;
  	for (dir=0; dir<2; dir++) 
  	{
    	for (side=0; side<2; side++) 
    	{
      		//ierr = IGASetBoundaryValue(iga,dir,side,dof,0.0);CHKERRQ(ierr);    				// Dirichlet boundary conditions
      		ierr = IGASetBoundaryForm(iga,dir,side,PETSC_TRUE);CHKERRQ(ierr);  				// Neumann boundary conditions
    	}
  	}

	//Aqui se crean la matriz y los vectores de nuestro sistema lineal
	Vec pi0;
	Mat Kl2;
	Vec Fl2;
	ierr = IGACreateVec(iga,&pi0);CHKERRQ(ierr);  
	ierr = IGACreateMat(iga,&Kl2);CHKERRQ(ierr);
	ierr = IGACreateVec(iga,&Fl2);CHKERRQ(ierr);
	ierr = IGASetFormSystem(iga,L2Projection,&userL2);CHKERRQ(ierr);
	ierr = IGAComputeSystem(iga,Kl2,Fl2);CHKERRQ(ierr);

	KSP kspl2;																			//Esta parte llama a KSP para resolver el sistema lineal
	ierr = IGACreateKSP(iga,&kspl2);CHKERRQ(ierr);
	ierr = KSPSetOperators(kspl2,Kl2,Kl2);CHKERRQ(ierr); 								//Esta función crea la matriz a resolver en el segundo parámetro y usa el 3er parámetro como preacondicionador
	ierr = KSPSetType(kspl2,KSPCG);CHKERRQ(ierr);										//Usando KSPCG (gradiente conjugado) porque la matriz es simétrica
	ierr = KSPSetOptionsPrefix(kspl2,"l2p_");CHKERRQ(ierr);
	ierr = KSPSetFromOptions(kspl2);CHKERRQ(ierr);
	//ierr = KSPSetTolerances(kspl2,1e-10,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);CHKERRQ(ierr);

	ierr = KSPSolve(kspl2,Fl2,pi0);CHKERRQ(ierr);

	ierr = KSPDestroy(&kspl2);CHKERRQ(ierr);
	ierr = MatDestroy(&Kl2);CHKERRQ(ierr);
	ierr = VecDestroy(&Fl2);CHKERRQ(ierr);
	ierr = IGAWriteVec(iga,pi0,"./results/Pl2-2d.dat");CHKERRQ(ierr);
	ierr = IGAWriteVec(iga,pi0,"./results/Pi-2d-0.dat");CHKERRQ(ierr);
//

//Creación de tipos y sistemas para Pi0
	Mat KPi;
	Vec x,FPi;
	ierr = IGACreateMat(iga,&KPi);CHKERRQ(ierr);
	ierr = IGACreateVec(iga,&x);CHKERRQ(ierr);
	ierr = IGACreateVec(iga,&FPi);CHKERRQ(ierr);
	//ierr = IGASetFormSystem(iga,System,NULL);CHKERRQ(ierr);
	//ierr = IGAComputeSystem(iga,KPi,FPi);CHKERRQ(ierr);

	AppCtx userPi;
	userPi.dt     = 0.0005;

	IGAPoint        point;
	IGAElement      elem;      			//element
	PetscReal       *Kloc,*Floc;     	//AA y BB
	PetscReal       *Kpoint,*Fpoint;	//KKK y FFF
	const PetscReal *arrayPi0;			//arrayU
	Vec  localPi0;						//localU
	PetscReal       *Pi0;				//U

  	IGAFormSystem  wtf;
 	void           *wtf2;

 	KSP ksp;
 	ierr = IGACreateKSP(iga,&ksp);CHKERRQ(ierr);

	char nombre[24];

	PetscInt i;
	for (i=0; i<900; i++) 
	{
		ierr = MatZeroEntries(KPi);CHKERRQ(ierr);
		ierr = VecZeroEntries(FPi);CHKERRQ(ierr);

		/* Get local vectors Pi0 and arrays */
		ierr = IGAGetLocalVecArray(iga,pi0,&localPi0,&arrayPi0);CHKERRQ(ierr);

		/* Element loop */
		ierr = IGABeginElement(iga,&elem);CHKERRQ(ierr);
		while (IGANextElement(iga,elem)) 
		{
			ierr = IGAElementGetWorkMat(elem,&Kloc);CHKERRQ(ierr);
			ierr = IGAElementGetWorkVec(elem,&Floc);CHKERRQ(ierr);
			ierr = IGAElementGetValues(elem,arrayPi0,&Pi0);CHKERRQ(ierr);

			/* FormSystem loop */
			while (IGAElementNextFormSystem(elem,&wtf,&wtf2)) {
			/* Quadrature loop */
			ierr = IGAElementBeginPoint(elem,&point);CHKERRQ(ierr);
			while (IGAElementNextPoint(elem,point)) 
			{
			ierr = IGAPointGetWorkMat(point,&Kpoint);CHKERRQ(ierr);
			ierr = IGAPointGetWorkVec(point,&Fpoint);CHKERRQ(ierr);
			ierr = Pi(point,Kpoint,Fpoint,Pi0,&userPi);CHKERRQ(ierr);
			ierr = IGAPointAddMat(point,Kpoint,Kloc);CHKERRQ(ierr);
			ierr = IGAPointAddVec(point,Fpoint,Floc);CHKERRQ(ierr);
			}
			ierr = IGAElementEndPoint(elem,&point);CHKERRQ(ierr);
		}

			ierr = IGAElementFixSystem(elem,Kloc,Floc);CHKERRQ(ierr);					//Esto pone condicion de Dirichlet ¿?
			ierr = IGAElementAssembleMat(elem,Kloc,KPi);CHKERRQ(ierr);
			ierr = IGAElementAssembleVec(elem,Floc,FPi);CHKERRQ(ierr);

		}
		ierr = IGAEndElement(iga,&elem);CHKERRQ(ierr);

		/* Restore local vectors Pi0 and arrays */
		ierr = IGARestoreLocalVecArray(iga,pi0,&localPi0,&arrayPi0);CHKERRQ(ierr);

		ierr = MatAssemblyBegin(KPi,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
		ierr = MatAssemblyEnd  (KPi,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
		ierr = VecAssemblyBegin(FPi);CHKERRQ(ierr);
		ierr = VecAssemblyEnd  (FPi);CHKERRQ(ierr);

		ierr = KSPSetOperators(ksp,KPi,KPi);CHKERRQ(ierr);
		ierr = KSPSetFromOptions(ksp);CHKERRQ(ierr);
		ierr = KSPSolve(ksp,FPi,x);CHKERRQ(ierr);

		if ((i+1)%4==0)
		{
			sprintf(nombre, "./results/Pi-2d-%d.dat", (i+1));								//Mete el numero al string
			ierr = IGAWriteVec(iga,x,nombre);CHKERRQ(ierr);									//Graba el resultado al disco
		}
		ierr=VecCopy(x,pi0); CHKERRQ(ierr);												//Aquí copiamos Pi(t+1) a Pi(t)
		PetscPrintf(PETSC_COMM_WORLD,"Iter= %d \n",(i+1));
	}

	//ierr = VecView(x,PETSC_VIEWER_DRAW_WORLD);CHKERRQ(ierr);

	ierr = KSPDestroy(&ksp);CHKERRQ(ierr);
	ierr = MatDestroy(&KPi);CHKERRQ(ierr);
	ierr = VecDestroy(&x);CHKERRQ(ierr);
	ierr = VecDestroy(&FPi);CHKERRQ(ierr);
	ierr = IGADestroy(&iga);CHKERRQ(ierr);

	ierr = PetscFinalize();CHKERRQ(ierr);
//

return 0;
}
