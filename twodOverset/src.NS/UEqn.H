// UEqn.H for the composite case

// RHS
PtrList<volVectorField> RHSComposite(nMeshes);

// Momentum eqn
PtrList<fvVectorMatrix> uEqnComposite(nMeshes);


forAll( RHSComposite, i)
{
    MRFComposite[i].correctBoundaryVelocity(uComposite[i]);

    uEqnComposite.set
    (
	i,
	new fvVectorMatrix
	(
	      fvm::div(phiComposite[i],uComposite[i])
	    + MRFComposite[i].DDt(uComposite[i])
	    + turbulenceComposite[i].divDevReff(uComposite[i])
	    == fvOptionsComposite[i](uComposite[i])
	)

    );

    uEqnComposite[i].relax();
    fvOptionsComposite[i].constrain(uEqnComposite[i]);

    // This is the RHS of the UEqn
    RHSComposite.set(i, new volVectorField("RHS",-fvc::grad(pComposite[i])));

    // This can be uncommented if a single mesh solution is sought
    //solve ( uEqnComposite[i] == RHSComposite[i] );

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
    
    SM::SparseMatrix<double>* Amat = oSolver.Umat;
    vectorField & source            = oSolver.sourceU;
    vectorField & solutionU	    = oSolver.solutionU;
    

    
    for(label i = 0 ; i < nMeshes ; i++)
    {
        const fvMesh & mesh = componentMesh[i];
        const fvVectorMatrix & uEqn = uEqnComposite[i];
        const surfaceScalarField & deltaCoeffs = mesh.surfaceInterpolation::deltaCoeffs();
        
        // Set the diagonal and the RHS
        // Note the addition of source() here. This is because UEqn() was relaxed and part of it has gone into the source term.
        forAll(uComposite[i],rowI)
        {
            label gRow = rowI + cumulativeSum[i];
            source[gRow].x() = uEqnComposite[i].source()[rowI].x() + mesh.V()[rowI]*RHSComposite[i][rowI].x();
 	    source[gRow].y() = uEqnComposite[i].source()[rowI].y() + mesh.V()[rowI]*RHSComposite[i][rowI].y();
	    source[gRow].z() = uEqnComposite[i].source()[rowI].z() + mesh.V()[rowI]*RHSComposite[i][rowI].z();
            Amat[0].set(gRow,gRow,uEqn.diag()[rowI]);
            Amat[1].set(gRow,gRow,uEqn.diag()[rowI]);
        }
        
        
        // Upper and Lower Coefficients (i.e. off-diagonal coefficients)
        for(label faceI= 0 ; faceI < mesh.nInternalFaces() ; faceI++)
        {
            label own = uEqn.lduAddr().lowerAddr()[faceI] + cumulativeSum[i];
            label nei = uEqn.lduAddr().upperAddr()[faceI] + cumulativeSum[i];
        
            Amat[0].set(own,nei,uEqn.upper()[faceI]);
            Amat[0].set(nei,own,uEqn.lower()[faceI]);

	    Amat[1].set(own,nei,uEqn.upper()[faceI]);
            Amat[1].set(nei,own,uEqn.lower()[faceI]);
        }
        
        
        // Boundary Conditions
        forAll(uComposite[i].boundaryFieldRef(),patchI)
        {
            const fvPatch & pp = uComposite[i].boundaryFieldRef()[patchI].patch();
            forAll(pp,faceI)
            {
                label cellI = pp.faceCells()[faceI] + cumulativeSum[i];
                
                scalar existingInternalValue = Amat[0].get(cellI,cellI);
                scalar newInternalValue = existingInternalValue + uEqn.internalCoeffs()[patchI][faceI].component(0);
                Amat[0].set(cellI,cellI,newInternalValue);

                existingInternalValue = Amat[1].get(cellI,cellI);
                newInternalValue = existingInternalValue + uEqn.internalCoeffs()[patchI][faceI].component(1);
                Amat[1].set(cellI,cellI,newInternalValue);

                source[cellI]+= uEqn.boundaryCoeffs()[patchI][faceI];
            }
        }
        
        
        // Interpolation Coefficients (Overset)
        forAll(uComposite[i],rowI)
        {   
            // On interpolation cells only
            if ( oSolver.connect[i]->donorCells[rowI]>=0 )
            {
                label gRow = rowI + cumulativeSum[i];
                
                // Make zero on all columns
                for (label jc = 0 ; jc < oSolver.nCellsTotal ; jc++)
                {
                    Amat[0].set(gRow,jc,0);
                    Amat[1].set(gRow,jc,0);
                }
                
                // Diagonal is one
                Amat[0].set(gRow,gRow,1);
                Amat[1].set(gRow,gRow,1);
                
                // No RHS
                source[gRow] = vector(0,0,0);
                
                // Insert weights as coefficients
                forAll(oSolver.connect[i]->donorCellCells[rowI], in_)
                {
                    // in a general case identify the donor mesh ID and based on donorMeshID and donorCell get the global column index(to do)
                    label gCol = oSolver.connect[i]->donorCellCells[rowI][in_]+cumulativeSum[i==0?1:0];
                    Amat[0].set(gRow,gCol,-oSolver.connect[i]->iweights[rowI][in_]);
		    Amat[1].set(gRow,gCol,-oSolver.connect[i]->iweights[rowI][in_]);
                }

            }
            
            
        }
    }

// Count number of non-zero elements
label nnz = 0;
for ( label i = 0 ;i < oSolver.nCellsTotal ; i++ )
{
    for (label j = 0 ; j < oSolver.nCellsTotal ; j++)
    {
        if ( mag(Amat[0].get(i,j)) > 0 )
        {
            nnz++;
        }
    }

    if ( mag(Amat[0].get(i,i)) < 1E-20 )
    {
        Info << " Zero diagonal on " << i << endl;
    }
   
}

// Solve here
 scalar iresidual=-1;
 for (label cmpt = 0 ; cmpt < 2 ; cmpt++)
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
                scalar entry = Amat[cmpt].get(iRow,jCol);
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

        // Solve the linear system (<>.begin() is the raw pointer to the scalarField from OpenFOAM)
	// Solve all equations segregated way
	
	scalarField solScalar(solutionU.component(cmpt));
	scalarField souScalar(source.component(cmpt));
	scalarField Ax(solScalar);
	Ax=0.0;

	// Get the initial residual
	    
	s.smat.matvec(&s.smat,solScalar.begin(),Ax.begin());
	scalarField res(souScalar-Ax);
	iresidual = max(itsol_norm(res.begin(),oSolver.nCellsTotal),iresidual); 

        itsol_solver_solve(&s, solScalar.begin(), souScalar.begin());
	Info << " Solving for U["<<cmpt<<"] :: (ILU) BiCGStab Solution converged from " << iresidual <<"  to " << s.res << " in " << s.nits << " iterations " << endl;
            
	// Replace back to original vector field
	solutionU.replace(cmpt,solScalar);


        // Clean-up the allocated arrays
        itsol_solver_finalize(&s);
        itsol_cleanCOO(&A);
	
	if ( iresidual < 1E-10 && cmpt > 0)
	convergenceSatisfied=true;

    }

    // Plug back solution to uComposite(OpenFOAM variable)
	label solutionIndex_=0;
	for ( label i = 0 ; i < nMeshes ; i++ )
	{
	    forAll(uComposite[i], cellI)
	    {
		uComposite[i].primitiveFieldRef()[cellI] = solutionU[solutionIndex_];
	        solutionIndex_++;
	    }
	    uComposite[i].correctBoundaryConditions();	
	}

#endif







 
