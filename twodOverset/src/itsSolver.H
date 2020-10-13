// ------------------------------------------------------------------------------
// IT Solver from Saad, Y. https://www-users.cs.umn.edu/~saad/software/ITSOL
// Also: https://github.com/huiscliu/itsol
// Current version links the original IT Solver and OpenFOAM.
//
// Blame:
// --------------
// D. Chandar 2020, Queen's University of Belfast
// d.chandar@qub.ac.uk, dominic.chandar@gmail.com 
// ------------------------------------------------------------------------------

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

        // Solve the linear system (<>.begin() is the raw pointer to the scalarField from OpenFOAM)
        itsol_solver_solve(&s, oSolver.solution.begin(), oSolver.source.begin());

        Info << " (ILU) BiCGStab Solution converged to " << s.res << " in " << s.nits << " iterations " << endl;

        // Clean-up the allocated arrays
        itsol_solver_finalize(&s);
        itsol_cleanCOO(&A);

    }
