#include "petiga.h"

#if PETSC_VERSION_LT(3,5,0)
#define KSPSetOperators(ksp,A,B) KSPSetOperators(ksp,A,B,SAME_NONZERO_PATTERN)
#endif

//Definicion de AppCtx
typedef struct {
  PetscReal dt;
} AppCtx;


// Proyección L2 de Pi(0)
#undef  __FUNCT__
#define __FUNCT__ "L2Projection"
PetscErrorCode L2Projection(IGAPoint p,PetscReal *K,PetscReal *F,void *ctx)
{
	if (p->atboundary) return 0;														//Condición de borde Neumann nula

	//AppCtx *user = (AppCtx *)ctx;

	PetscInt nen = p->nen;																//Numero de funciones de forma, en este caso es 9
	PetscInt dim = p->dim;																//Dimension del problema

	PetscReal x[dim];																	//Vector de reales, tamaño igual a dimensión del problema
	IGAPointFormGeomMap(p,x);															//llena x con las coordenadas de p, punto de Gauss

  //PetscReal g;
  //PetscReal disinter = sqrt(x[0]*x[0] + x[1]*x[1]) - 64.;
  //g = 0.5 + 0.5*tanh(disinter/(sqrt(2.)*epsilon));

	//g es la función a proyectar
	PetscReal g;
	if (x[0]>-1.0/100.0 && x[0]<1.0/100.0 && x[1]>-1.0/100.0 && x[1]<1.0/100.0)
	{
		g=1.0;
	}
	else
	{
		g=0.0;
	}

	const PetscReal *N = (typeof(N)) p->shape[0];
	PetscReal (*FF)[1] = (PetscReal (*)[1])F;
	PetscReal (*KK)[1][nen][1] = (PetscReal (*)[1][nen][1])K;
	PetscInt a,b;
	for(a=0; a<nen; a++)
	{
		for(b=0; b<nen; b++) 
		{
	  		KK[a][0][b][0] = N[a]*N[b];
		}
		FF[a][0] = N[a]*g;
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
	PetscInt a,b,nen=p->nen;

	PetscReal Pi0[1];																	//Aquí asignamos Pi(t) a un vector
	PetscReal d_Pi0[1][2];																//Idem para las derivadas parciales
	IGAPointFormValue(p,U,&Pi0[0]);														
	IGAPointFormGrad (p,U,&d_Pi0[0][0]);

	PetscScalar (*Kpi)[1][nen][1] = (typeof(Kpi)) K;

	PetscReal V1=1.0;
	PetscReal V2=0.0;
	PetscReal V1_x=0.0;
	//PetscReal V2_x=0.0;
	//PetscReal V1_y=0.0;
	PetscReal V2_y=0.0;

	PetscInt dim = p->dim;																//Dimension del problema
	PetscReal x[dim];
	IGAPointFormGeomMap(p,x);


	PetscPrintf(PETSC_COMM_WORLD,"At boundary %d \n",p->atboundary);

        if (p->atboundary) goto boundary;

	//PetscPrintf(PETSC_COMM_WORLD,"p->nen: %d \n", p->nen);
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
		}
		F[a] = 	Pi0[0]*Na 																//(LS26)
				+dt*(Pi0[0]*(Na_x*V1+Na_y*V2))											//(LS27)
				+dt*(Pi0[0]*Na*(V1_x+V2_y))												//(LS28)
				+Pi0[0]*Na 																//(GL2)
				;
	}
	return 0;


 	boundary:;

  	PetscInt dir  = p->boundary_id / 2;
  	PetscInt side = p->boundary_id % 2;

	for (a=0; a<nen; a++) {
	
 	F[a] = 0.;

        }

  	return 0;
}
//

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc, char *argv[]) {

//Creación de tipos y sistemas para proyección L2
	void* user;
	PetscErrorCode  ierr;
	ierr = PetscInitialize(&argc,&argv,0,0);CHKERRQ(ierr);								//Siempre inicializar petsc

	IGA iga;	//Aqui se define los parámetros para IGA
	ierr = IGACreate(PETSC_COMM_WORLD,&iga);CHKERRQ(ierr);
	ierr = IGASetDim(iga,2);CHKERRQ(ierr);												//Dimnsion espacial del problema
	ierr = IGASetDof(iga,1);CHKERRQ(ierr);												//Numero de grados de libertad
	ierr = IGASetOrder(iga,1);CHKERRQ(ierr);											//Numero máximo de derivadas a calcular
	ierr = IGASetFromOptions(iga);CHKERRQ(ierr);										//Nota: El orden (o grado) de las funciones de forma viene dado por la malla!
	ierr = IGARead(iga,"./geometry.dat");CHKERRQ(ierr);
	ierr = IGASetUp(iga);CHKERRQ(ierr);

    	PetscInt dir,side,dof=0;
  	for (dir=0; dir<2; dir++) {
    		for (side=0; side<2; side++) {
      		//ierr = IGASetBoundaryValue(iga,dir,side,dof,0.0);CHKERRQ(ierr);    // Dirichlet boundary conditions
      		ierr = IGASetBoundaryForm(iga,dir,side,PETSC_TRUE);CHKERRQ(ierr);  // Neumann boundary conditions
    		}
  	}

	//Aqui se crean la matriz y los vectores de nuestro sistema lineal
	Vec pi0;
	Mat Kl2;
	Vec Fl2;
	ierr = IGACreateVec(iga,&pi0);CHKERRQ(ierr);  
	ierr = IGACreateMat(iga,&Kl2);CHKERRQ(ierr);
	ierr = IGACreateVec(iga,&Fl2);CHKERRQ(ierr);
	ierr = IGASetFormSystem(iga,L2Projection,&user);CHKERRQ(ierr);
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
	userPi.dt     = 0.001;

	KSP ksp;
	IGAPoint        point;
	IGAElement      elem;      			//element
	PetscReal       *Kloc,*Floc;     	//AA y BB
	PetscReal       *Kpoint,*Fpoint;	//KKK y FFF
	const PetscReal *arrayPi0;			//arrayU
	Vec  localPi0;						//localU
	PetscReal       *Pi0;				//U

	char nombre[24];

	PetscInt i;
	for (i=0; i<5; i++) 
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

		ierr = IGACreateKSP(iga,&ksp);CHKERRQ(ierr);
		ierr = KSPSetOperators(ksp,KPi,KPi);CHKERRQ(ierr);
		ierr = KSPSetFromOptions(ksp);CHKERRQ(ierr);
		ierr = KSPSolve(ksp,FPi,x);CHKERRQ(ierr);

		sprintf(nombre, "./results/Pi-2d-%d.dat", (i+1));								//Mete el numero al string
		ierr = IGAWriteVec(iga,x,nombre);CHKERRQ(ierr);

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
