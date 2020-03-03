#include <stdio.h>
#include <dirent.h>
#include <math.h>
#include "petiga.h"

#if PETSC_VERSION_LT(3,5,0)
#define KSPSetOperators(ksp,A,B) KSPSetOperators(ksp,A,B,SAME_NONZERO_PATTERN)
#endif

//Definicion de AppCtxL2
typedef struct {
  PetscReal Lx;
  PetscReal Ly;
  PetscReal nx;
  PetscReal ny;
} AppCtxL2;

//Definicion de AppCtx
typedef struct {
  PetscReal dt;
} AppCtxPi;


// Proyección L2 de Pi(0)
#undef  __FUNCT__
#define __FUNCT__ "L2Projection"
PetscErrorCode L2Projection(IGAPoint p,PetscReal *K,PetscReal *F,void *ctx)
{
	
	AppCtxL2    *user  = (AppCtxL2 *)ctx;
	PetscReal Lx     = user->Lx;
	PetscReal Ly     = user->Ly;
	PetscReal nx     = user->nx;
	PetscReal ny     = user->ny;

	if (p->atboundary) return 0;														//Condición de borde Neumann nula

	PetscInt a,b,i,j;

	PetscInt nen = p->nen;																//Numero de funciones de forma, en este caso es 9
	PetscInt dim = p->dim;																//Dimension del problema
	PetscInt dof = p->dof;																//Dimension del problema

	PetscReal x[dim];																	//Vector de reales, tamaño igual a dimensión del problema
	IGAPointFormGeomMap(p,x);															//llena x con las coordenadas de p, punto de Gauss

	//g es la función a proyectar
	PetscReal g[dof];

	PetscReal dx=Lx/nx;
	PetscReal dy=Ly/ny;
	PetscReal eps=0.005;

	g[0]=1.917545*(0.5*(tanh((x[0]+dx)/eps)+1.0)-0.5*(tanh((x[0]-dx)/eps)+1.0))*(0.5*(tanh((x[1]+dy)/eps)+1.0)-0.5*(tanh((x[1]-dy)/eps)+1.0));
	g[1]=0;
	g[2]=0;
	g[3]=0;

	const PetscReal (*N) = (typeof(N)) p->shape[0];

	PetscReal (*FF)[dof] = (PetscReal (*)[dof])F;
	PetscReal (*KK)[dof][nen][dof] = (PetscReal (*)[dof][nen][dof])K;
	
	for(a=0; a<nen; a++)
	{
		for(i=0; i<dof; i++)
		{
			for(b=0; b<nen; b++) 
			{
				for(j=0; j<dof; j++)
				{
					if(i==j)
					KK[a][i][b][j] = N[a]*N[b];
				}
			}
			FF[a][i] = N[a]*g[i];
		}
	}
	return 0;
}

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

//Creacion de la malla
	PetscReal Lx=0.5;
	PetscReal Ly=0.5;
	PetscInt  nx=101;
	PetscInt  ny=101;
	char      nombreMalla[12]="geometry";
	char 	  comando[512];

	sprintf(comando,"exec python rectangle.py %s %f %f %d %d\n",nombreMalla,Lx,Ly,nx,ny);
	system(comando);


//Creación de tipos y sistemas para proyección L2 y solución del problema
	AppCtxL2 userL2;
	userL2.Lx     = Lx;
	userL2.Ly     = Ly;
	userL2.nx     = nx;
	userL2.ny     = ny;

	PetscErrorCode  ierr;
	ierr = PetscInitialize(&argc,&argv,0,0);CHKERRQ(ierr);								//Siempre inicializar petsc

	IGA iga;	//Aqui se define los parámetros para IGA
	ierr = IGACreate(PETSC_COMM_WORLD,&iga);CHKERRQ(ierr);
	ierr = IGASetDim(iga,2);CHKERRQ(ierr);												//Dimnsion espacial del problema
	ierr = IGASetDof(iga,4);CHKERRQ(ierr);												//Numero de grados de libertad
	ierr = IGASetOrder(iga,1);CHKERRQ(ierr);											//Numero máximo de derivadas a calcular (Nota: El orden (o grado) de las funciones de forma viene dado por la malla!)
	ierr = IGASetFromOptions(iga);CHKERRQ(ierr);
	ierr = IGARead(iga,"./geometry.dat");CHKERRQ(ierr);
	ierr = IGASetUp(iga);CHKERRQ(ierr);

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

// Finalizacion
	ierr = IGADestroy(&iga);CHKERRQ(ierr);
	ierr = PetscFinalize();CHKERRQ(ierr);
//

//Post proceso
	//system("exec python pospro2.py Pi-2d*");											// Post proceso para visualizar en VTK

return 0;
}