

//System for Alfa
	IGA igaAl;
	ierr = IGACreate(PETSC_COMM_WORLD,&igaAl);CHKERRQ(ierr);
	ierr = IGASetDim(igaAl,2);CHKERRQ(ierr);													//Spatial dimension of the problem
	ierr = IGASetDof(igaAl,2);CHKERRQ(ierr);													//Number of dof per node
	ierr = IGASetOrder(igaAl,1);CHKERRQ(ierr);													//Number of derivatives to calculate
	ierr = IGASetFromOptions(igaAl);CHKERRQ(ierr);												
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

//System for Z
	IGA igaZ;
	ierr = IGACreate(PETSC_COMM_WORLD,&igaZ);CHKERRQ(ierr);
	ierr = IGASetDim(igaZ,2);CHKERRQ(ierr);														//Spatial dimension of the problem
	ierr = IGASetDof(igaZ,2);CHKERRQ(ierr);														//Number of dof per node
	ierr = IGASetOrder(igaZ,1);CHKERRQ(ierr);													//Number of derivatives to calculate
	ierr = IGASetFromOptions(igaZ);CHKERRQ(ierr);												
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



//Solution loops
	ierr = IGABeginElement(igaZ,&elemZ);CHKERRQ(ierr);
	ierr = IGABeginElement(igaAl,&elemAl);CHKERRQ(ierr);

	PetscPrintf(PETSC_COMM_WORLD,"elemZ count= %d \n",elemZ->count);
	PetscPrintf(PETSC_COMM_WORLD,"elemAl count= %d \n",elemAl->count);

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
			PetscPrintf(PETSC_COMM_WORLD,"pointZ count= %d \n",pointZ->count);				//This returns 2
			PetscPrintf(PETSC_COMM_WORLD,"pointAl count= %d \n",pointAl->count);			//This returns 4
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
			ierr = IGAElementEndPoint(elemAl,&pointAl);CHKERRQ(ierr);		//Due to the mentioned above, program crashes at this line
		}
		ierr = IGAElementFixSystem(elemZ,KlocZ,FlocZ);CHKERRQ(ierr);
		ierr = IGAElementAssembleMat(elemZ,KlocZ,Kz);CHKERRQ(ierr);
		ierr = IGAElementAssembleVec(elemZ,FlocZ,Fz);CHKERRQ(ierr);

	}
	IGANextElement(igaAl,elemAl);
	ierr = IGAEndElement(igaZ,&elemZ);CHKERRQ(ierr);
	ierr = IGAEndElement(igaAl,&elemAl);CHKERRQ(ierr);
