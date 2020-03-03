#include "petiga.h"

#if PETSC_VERSION_LT(3,5,0)
#define KSPSetOperators(ksp,A,B) KSPSetOperators(ksp,A,B,SAME_NONZERO_PATTERN)
#endif


typedef struct {
  PetscReal alpha, beta, kappa;
  Vec vecG;
  PetscScalar *G;

  PetscInt nsd;
  PetscInt dim;
  PetscInt ndof;
  PetscInt nen;
  PetscInt nqp;
} AppCtx;

#undef  __FUNCT__
#define __FUNCT__ "System"
PetscErrorCode System(IGAPoint p,PetscScalar *K,PetscScalar *F, PetscScalar *U,void *ctx)
{
  const PetscReal *N0,(*N1)[2];
  IGAPointGetShapeFuns(p,0,(const PetscReal**)&N0);
  IGAPointGetShapeFuns(p,1,(const PetscReal**)&N1);
  PetscInt a,b,nen=p->nen;

  PetscScalar u[1];
  PetscScalar grad_u[1][2];
  IGAPointFormValue(p,U,&u[0]);
  IGAPointFormGrad (p,U,&grad_u[0][0]);

  //PetscPrintf(PETSC_COMM_WORLD,"p->nen: %d \n", p->nen);
  for (a=0; a<nen; a++) {
    PetscReal Na   = N0[a];
    PetscReal Na_x = N1[a][0];
    PetscReal Na_y = N1[a][1];
    for (b=0; b<nen; b++) {
      PetscReal Nb_x = N1[b][0];
      PetscReal Nb_y = N1[b][1];
      K[a*nen+b] = Na_x*Nb_x + Na_y*Nb_y;
    }
    F[a] = Na * 1.0;
  }
  return 0;
}


#undef  __FUNCT__
#define __FUNCT__ "L2Projection"
PetscErrorCode L2Projection(IGAPoint p,PetscScalar *KK,PetscScalar *FF,void *ctx)
{
  if (p->atboundary) return 0;

  AppCtx *user = (AppCtx *)ctx;

  PetscInt nen = p->nen;
  PetscInt dim = p->dim;

  PetscReal x[dim];
  IGAPointFormGeomMap(p,x);

  // g is the function that we want to project !!
  PetscReal g;
  g = sqrt(x[0]*x[0] + x[1]*x[1]);

  const PetscReal *N = (typeof(N)) p->shape[0];

  PetscScalar (*F)[1] = (PetscScalar (*)[1])FF;
  PetscScalar (*K)[1][nen][1] = (PetscScalar (*)[1][nen][1])KK;

  PetscInt a,b;
  for(a=0; a<nen; a++){

  F[a][0] = N[a]*g;

    for(b=0; b<nen; b++) {

          K[a][0][b][0] = N[a]*N[b];

    }

  }
  return 0;
}

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc, char *argv[]) {

  PetscErrorCode  ierr;
  ierr = PetscInitialize(&argc,&argv,0,0);CHKERRQ(ierr);

  /* Define simulation specific parameters */
  AppCtx user;
  user.alpha = 0.02; /*  */           // 0.005     0.02      2.0
  user.beta = 0.5;    /*  */          // 0.02      0.5      2.0
  user.kappa  = 0.02;   /*  */        // 0.005      0.2      0.

  user.nsd = 2;
  user.dim = 2;
  user.ndof = 2;
  user.nen = 9;  // quadratics 9  cubics 16   quartics 25
  user.nqp = 9;


  IGA iga;
  ierr = IGACreate(PETSC_COMM_WORLD,&iga);CHKERRQ(ierr);
  ierr = IGASetDim(iga,2);CHKERRQ(ierr);
  ierr = IGASetDof(iga,1);CHKERRQ(ierr);
  ierr = IGASetOrder(iga,1);CHKERRQ(ierr);  // Number of spatial derivatives that are needed !!
  ierr = IGASetFromOptions(iga);CHKERRQ(ierr);
  ierr = IGARead(iga,"./geometry.dat");CHKERRQ(ierr);
  ierr = IGASetUp(iga);CHKERRQ(ierr);
  
  PetscInt dir,side;
  PetscScalar value = 0.0;
  for (dir=0; dir<2; dir++) {
    for (side=0; side<2; side++) {
      ierr = IGASetBoundaryValue(iga,dir,side,0,value);CHKERRQ(ierr);
    }
  }

  Vec pi0;
  ierr = IGACreateVec(iga,&pi0);CHKERRQ(ierr);
  Mat Al2;
  Vec bl2;
  ierr = IGACreateMat(iga,&Al2);CHKERRQ(ierr);
  ierr = IGACreateVec(iga,&bl2);CHKERRQ(ierr);
  ierr = IGASetFormSystem(iga,L2Projection,&user);CHKERRQ(ierr);
  ierr = IGAComputeSystem(iga,Al2,bl2);CHKERRQ(ierr);
  KSP kspl2;
  ierr = IGACreateKSP(iga,&kspl2);CHKERRQ(ierr);
  ierr = KSPSetOperators(kspl2,Al2,Al2);CHKERRQ(ierr);
  ierr = KSPSetType(kspl2,KSPCG);CHKERRQ(ierr);    // Do always this when the matrix is symmetric !!

  ierr = KSPSetOptionsPrefix(kspl2,"l2p_");CHKERRQ(ierr);
  ierr = KSPSetFromOptions(kspl2);CHKERRQ(ierr);
  //ierr = KSPSetTolerances(kspl2,1e-10,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);CHKERRQ(ierr);

  ierr = KSPSolve(kspl2,bl2,pi0);CHKERRQ(ierr);

  ierr = KSPDestroy(&kspl2);CHKERRQ(ierr);
  ierr = MatDestroy(&Al2);CHKERRQ(ierr);
  ierr = VecDestroy(&bl2);CHKERRQ(ierr);

  //Aqui empieza PI

  Mat A;
  Vec x,b;
  ierr = IGACreateMat(iga,&A);CHKERRQ(ierr);
  ierr = IGACreateVec(iga,&x);CHKERRQ(ierr);
  ierr = IGACreateVec(iga,&b);CHKERRQ(ierr);
  //ierr = IGASetFormSystem(iga,System,NULL);CHKERRQ(ierr);
  //ierr = IGAComputeSystem(iga,A,b);CHKERRQ(ierr);


  KSP ksp;
  IGAElement     element;
  IGAPoint       point;
  PetscScalar    *AA,*BB;
  PetscScalar    *KKK,*FFF;
  Vec  localU;
  const PetscScalar *arrayU;
  PetscScalar  *U;

  PetscInt i;
  for (i=0; i<100; i++) {

  ierr = MatZeroEntries(A);CHKERRQ(ierr);
  ierr = VecZeroEntries(b);CHKERRQ(ierr);

  /* Get local vectors U and arrays */
  ierr = IGAGetLocalVecArray(iga,pi0,&localU,&arrayU);CHKERRQ(ierr);

  /* Element loop */
  ierr = IGABeginElement(iga,&element);CHKERRQ(ierr);
  while (IGANextElement(iga,element)) {
    ierr = IGAElementGetWorkMat(element,&AA);CHKERRQ(ierr);
    ierr = IGAElementGetWorkVec(element,&BB);CHKERRQ(ierr);

    ierr = IGAElementGetValues(element,arrayU,&U);CHKERRQ(ierr);
      /* Quadrature loop */
      ierr = IGAElementBeginPoint(element,&point);CHKERRQ(ierr);
      while (IGAElementNextPoint(element,point)) {
        ierr = IGAPointGetWorkMat(point,&KKK);CHKERRQ(ierr);
        ierr = IGAPointGetWorkVec(point,&FFF);CHKERRQ(ierr);
        ierr = System(point,KKK,FFF,U,&user);CHKERRQ(ierr);
        ierr = IGAPointAddMat(point,KKK,AA);CHKERRQ(ierr);
        ierr = IGAPointAddVec(point,FFF,BB);CHKERRQ(ierr);
      }
      ierr = IGAElementEndPoint(element,&point);CHKERRQ(ierr);

    ierr = IGAElementFixSystem(element,AA,BB);CHKERRQ(ierr);
    ierr = IGAElementAssembleMat(element,AA,A);CHKERRQ(ierr);
    ierr = IGAElementAssembleVec(element,BB,b);CHKERRQ(ierr);
  }
  ierr = IGAEndElement(iga,&element);CHKERRQ(ierr);

  /* Restore local vectors U and arrays */
  ierr = IGARestoreLocalVecArray(iga,pi0,&localU,&arrayU);CHKERRQ(ierr);

  ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd  (A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = VecAssemblyBegin(b);CHKERRQ(ierr);
  ierr = VecAssemblyEnd  (b);CHKERRQ(ierr);
  
  ierr = IGACreateKSP(iga,&ksp);CHKERRQ(ierr);
  ierr = KSPSetOperators(ksp,A,A);CHKERRQ(ierr);
  ierr = KSPSetFromOptions(ksp);CHKERRQ(ierr);
  ierr = KSPSolve(ksp,b,x);CHKERRQ(ierr);

  ierr = IGAWriteVec(iga,x,"./results/poisson2d.dat");CHKERRQ(ierr);

  }

  //ierr = VecView(x,PETSC_VIEWER_DRAW_WORLD);CHKERRQ(ierr);
  
  ierr = KSPDestroy(&ksp);CHKERRQ(ierr);
  ierr = MatDestroy(&A);CHKERRQ(ierr);
  ierr = VecDestroy(&x);CHKERRQ(ierr);
  ierr = VecDestroy(&b);CHKERRQ(ierr);
  ierr = IGADestroy(&iga);CHKERRQ(ierr);

  ierr = PetscFinalize();CHKERRQ(ierr);
  return 0;
}
