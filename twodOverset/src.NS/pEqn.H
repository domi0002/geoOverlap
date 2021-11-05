// pEqn.H for the composite case
{

// RHS
PtrList<volScalarField> RHSComposite(nMeshes);

// pressure eqn
PtrList<fvScalarMatrix> pEqnComposite(nMeshes);

// Required variables
PtrList<volScalarField> rAUComposite(nMeshes);

// Required variables
PtrList<surfaceScalarField> rAUfComposite(nMeshes);

// Required variables
PtrList<volVectorField> HbyAComposite(nMeshes);

// Required variables
PtrList<surfaceScalarField> phiHbyAComposite(nMeshes);


forAll( RHSComposite, i)
{
    pComposite[i].storePrevIter();

    rAUComposite.set(i, new volScalarField(1.0/uEqnComposite[i].A()));
    HbyAComposite.set(i, new volVectorField(constrainHbyA(rAUComposite[i]*uEqnComposite[i].H(), uComposite[i],pComposite[i])));
    rAUfComposite.set(i, new surfaceScalarField("rAUf", fvc::interpolate(rAUComposite[i])));
}

//+++ Not reequired - experimental
//oSolver.explicitInterpolate(HbyAComposite);

forAll( RHSComposite, i)
{
 
    phiHbyAComposite.set(i, new surfaceScalarField("phiHbyA", fvc::flux(HbyAComposite[i])));
    MRFComposite[i].makeRelative(phiHbyAComposite[i]);
    adjustPhi(phiHbyAComposite[i],uComposite[i],pComposite[i]);


    constrainPressure(pComposite[i], uComposite[i], phiHbyAComposite[i], rAUComposite[i], MRFComposite[i]);

    // The Source is constructed at a later stage
    pEqnComposite.set
    (
	i,
	new fvScalarMatrix
	(
	      fvm::laplacian(rAUfComposite[i],pComposite[i])
	)
    );


    // If one wants to solve on a single mesh, then this can be un-commented
    /*
    pEqnComposite.set
    (
	i,
	new fvScalarMatrix
	(
	      fvm::laplacian(rAUfComposite[i],pComposite[i]) == fvc::div(phiHbyAComposite[i])
	
        )
    );


    pEqnComposite[i].setReference(pRefCell, pRefValue);

    pEqnComposite[i].solve();
    */


    // This is the RHS of the UEqn
    RHSComposite.set(i, new volScalarField("RHS",fvc::div(phiHbyAComposite[i])));

    

}

#if true
// Now comes all the extra stuff
// Frame the Global Matrix of all meshes. Now we are going to use a custom sparse-matrix representation
// The custom sparse-matrix is from https://github.com/uestla/Sparse-Matrix 
    
    // First need to create local mesh to global matrix mapping
    labelList cumulativeSum(nMeshes);
    cumulativeSum[0] = 0;
    
    forAll(cumulativeSum,j)
    {
        if (j > 0)
        {
            cumulativeSum[j]=cumulativeSum[j-1]+componentMesh[j-1].nCells();
        }
    }
    
    SM::SparseMatrix<double> & Amat = oSolver.Pmat;
    scalarField & source            = oSolver.source;
    scalarField & solution	    = oSolver.solution;
    
    
    for(label i = 0 ; i < nMeshes ; i++)
    {
        const fvMesh & mesh = componentMesh[i];
        const fvScalarMatrix & pEqn = pEqnComposite[i];
        const surfaceScalarField & deltaCoeffs = mesh.surfaceInterpolation::deltaCoeffs();
        
        // Set the diagonal and the RHS
        forAll(pComposite[i],rowI)
        {
            label gRow   = rowI + cumulativeSum[i];
            source[gRow] = mesh.V()[rowI]*RHSComposite[i][rowI];
            Amat.set(gRow,gRow,pEqn.diag()[rowI]);
        }
        
        
        // Upper and Lower Coefficients (i.e. off-diagonal coefficients)
        for(label faceI= 0 ; faceI < mesh.nInternalFaces() ; faceI++)
        {
            label own = pEqn.lduAddr().lowerAddr()[faceI] + cumulativeSum[i];
            label nei = pEqn.lduAddr().upperAddr()[faceI] + cumulativeSum[i];
        
            Amat.set(own,nei,pEqn.upper()[faceI]);
            Amat.set(nei,own,pEqn.upper()[faceI]);
        }
        
        
        // Boundary Conditions
        forAll(pComposite[i].boundaryFieldRef(),patchI)
        {
            const fvPatch & pp = pComposite[i].boundaryFieldRef()[patchI].patch();
            forAll(pp,faceI)
            {
                label cellI = pp.faceCells()[faceI] + cumulativeSum[i];
                
		// Might throw vector assignment error?
                scalar existingInternalValue = Amat.get(cellI,cellI);
                scalar newInternalValue = existingInternalValue + pEqn.internalCoeffs()[patchI][faceI];
                
                Amat.set(cellI,cellI,newInternalValue);
                source[cellI]+= pEqn.boundaryCoeffs()[patchI][faceI];
            }
        }
        
	//Info <<"pRefCell =" << pRefCell << endl;
	
	// Set the Neumann BC fix on mesh zero. This is identical to pEqn.setReference().
	if ( i==0 )
	{
	    label cell2fix=pRefCell[i];

	    for (label jc = 0 ; jc < oSolver.nCellsTotal ; jc++)
            {
                Amat.set(cell2fix,jc,0);
            }

	    // Is this a good way? While this is correct there could be other ways.
	    // Diagonal is one
	    Amat.set(cell2fix,cell2fix,1.0);
	  
	    source[cell2fix]=pRefValue[i];
	    
	}
        
        // Interpolation Coefficients (Overset)
        forAll(pComposite[i],rowI)
        {   
            // On interpolation cells only
            if ( oSolver.connect[i]->donorCells[rowI]>=0 )
            {
                label gRow = rowI + cumulativeSum[i];
                
                // Make zero on all columns
                for (label jc = 0 ; jc < oSolver.nCellsTotal ; jc++)
                {
                    Amat.set(gRow,jc,0);
                }
                
                // Diagonal is one
                Amat.set(gRow,gRow,1);
                
                // No RHS
                source[gRow] = 0;
                
                // Insert weights as coefficients
                forAll(oSolver.connect[i]->donorCellCells[rowI], in_)
                {
                    // in a general case identify the donor mesh ID and based on donorMeshID and donorCell get the global column index(to do)
                    label gCol = oSolver.connect[i]->donorCellCells[rowI][in_]+cumulativeSum[i==0?1:0];
                    Amat.set(gRow,gCol,-oSolver.connect[i]->iweights[rowI][in_]);
                }

            }
            
            
        }

    // Conservative correction if-required
    #include "cons.H"

    }// i=0 to Meshes

// Count number of non-zero elements
label nnz = 0;
for ( label i = 0 ;i < oSolver.nCellsTotal ; i++ )
{
    for (label j = 0 ; j < oSolver.nCellsTotal ; j++)
    {
        if ( mag(Amat.get(i,j)) > 0 )
        {
            nnz++;
        }
    }

    if ( mag(Amat.get(i,i)) < 1E-20 )
    {
        Info << " Zero diagonal on " << i << endl;
    }
   
}

// Solve here
 {
        ITS_CooMat  A;
        ITS_SOLVER  s;

        bzero(&A, sizeof(A));
        
        // Size of the matrix
        A.n  = oSolver.nCellsTotal;
        
        // Number of non-zeros
        A.nnz= nnz;
        
        // Allocate non-zero entry 
        A.ma = (double *)itsol_malloc(nnz * sizeof(double), "read_coo:3");
        
        // Allocate columns
        A.ja = (int *)itsol_malloc(nnz * sizeof(int), "read_coo:4");
        
        // Allocate rows
        A.ia = (int *)itsol_malloc(nnz * sizeof(int), "read_coo:5");

        
        // Insert into matrix 'A' at non-zero locations based on the current Sparse-Matrix format from 'Amat'
        label m=0;
        for( label iRow = 0; iRow < oSolver.nCellsTotal; iRow++)
        {
            for(label jCol = 0 ; jCol < oSolver.nCellsTotal ; jCol++)
            {
                scalar entry = Amat.get(iRow,jCol);
                if ( mag(entry) > 0 )
                {
                    A.ia[m] = iRow;
                    A.ja[m] = jCol;
                    A.ma[m] = entry;
                    m++;
                }
            }
        }
       
        // Set the solver type -'BiCGStab' along with convergence conditions
        itsol_solver_initialize(&s, ITS_SOLVER_BICGSTAB, ITS_PC_ILUK, &A);
        s.pars.tol=1E-15;
	itsol_solver_assemble(&s);

	scalarField Ax(solution);
	Ax=0.0;

	// Get the initial residual
	    
	s.smat.matvec(&s.smat,solution.begin(),Ax.begin());
	scalarField res(source-Ax);
	scalar iresidual = itsol_norm(res.begin(),oSolver.nCellsTotal); 


        itsol_solver_solve(&s, solution.begin(), source.begin());
	Info << " Solving for p    :: (ILU) BiCGStab Solution converged from "<<iresidual<<" to " << s.res << " in " << s.nits << " iterations " << endl;
            
	// Plug back solution to pressure
	label solutionIndex_=0;
	for ( label i = 0 ; i < nMeshes ; i++ )
	{
	    forAll(pComposite[i], cellI)
	    {
		pComposite[i].primitiveFieldRef()[cellI] = solution[solutionIndex_];
	 	solutionIndex_++;
	    }	
	    
	}

        // Clean-up the allocated arrays
        itsol_solver_finalize(&s);
        itsol_cleanCOO(&A);

	if ( iresidual < 1E-6 && convergenceSatisfied==true)
	convergenceSatisfied=true;

    }
#endif
    
    // Do the rest of the stuff after solving for pressure. 
    for ( label i = 0 ; i < nMeshes ; i++ )
    {
	pComposite[i].correctBoundaryConditions();
	phiComposite[i] = phiHbyAComposite[i] - pEqnComposite[i].flux();

	#include "fieldContinuityErrs.H"

	pComposite[i].relax();

	// Recompute U
	uComposite[i] = HbyAComposite[i] - rAUComposite[i]*fvc::grad(pComposite[i]);


	uComposite[i].correctBoundaryConditions();
	fvOptionsComposite[i].correct(uComposite[i]);
	
    }

    // We have to interpolate U explicitly here since it has been corrected.
    oSolver.explicitInterpolate(uComposite);

    for ( label i = 0 ; i < nMeshes ; i++ )
    {
	uComposite[i].correctBoundaryConditions();
    }

    	#include "massContinuityErrs.H"

}


 
