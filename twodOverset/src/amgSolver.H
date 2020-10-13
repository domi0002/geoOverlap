 bool amgclSolver=false;
    
    if (amgclSolver)
    {
        // CRS format
        std::vector<double> val;
        std::vector<int>    col_ind;
        std::vector<int>    row_ptr;
        std::vector<double> rhs;
        label m = 0;
        
        val.reserve(nnz);
        col_ind.reserve(nnz);
        row_ptr.reserve(oSolver.nCellsTotal+1);
        rhs.reserve(oSolver.nCellsTotal);
        
        
        val.clear();
        col_ind.clear();
        row_ptr.clear();
        rhs.clear();
        
        // First element *always* is non-zero
        row_ptr.push_back(0);
        
        for( label iRow = 0; iRow < oSolver.nCellsTotal; iRow++)
        {
            for(label jCol = 0 ; jCol < oSolver.nCellsTotal ; jCol++)
            {
                scalar entry = Amat.get(iRow,jCol);
                if ( mag(entry) > 0 )
                {
                    val.push_back(entry);
                    col_ind.push_back(jCol);
                    m++;
                }
                // Clear out Amat for next call(time)
                //Amat.set(iRow,jCol,0.0);
                    
            }
            rhs.push_back(source[iRow]);
            row_ptr.push_back(m);
            
            // Clear out source for next call(time)
            //source[iRow]=0.0;
        }
        Info << "begin3" << endl;
        row_ptr.push_back(nnz);
        
        // Define the backend: Use th builtin backend for now
        typedef amgcl::backend::builtin<double> Backend;
        
        // bicgstab with Multigrid as preconditioner
        typedef amgcl::make_solver<
                                    amgcl::amg<
                                    Backend,
                                    amgcl::coarsening::ruge_stuben, 
                                    amgcl::relaxation::ilu0
                                  >
                                  ,amgcl::solver::bicgstab<Backend>
                                  > Solver;
        // Solver tolerance
        Solver::params prm;
        prm.solver.tol = 1e-15;
        prm.solver.maxiter=800;
         
        // Define the Solver object 'solve'
        Solver solve( std::tie(oSolver.nCellsTotal, row_ptr, col_ind, val) , prm);  
        
        
        //std::vector<double> stlSolution(nCellsTotal, 0.0);
        int    iters=0;
        double error=0;
        std::tie(iters, error) = solve(rhs, oSolver.stlSolution);

        std::cout << "Iterations = " << iters << "  error = " << error << std::endl;
        
        // Clear-out for next call
        val.clear();
        col_ind.clear();
        row_ptr.clear();
        rhs.clear();
        
        
    }
    else
    {
        ITS_CooMat  A;
        ITS_SOLVER  s;

        bzero(&A, sizeof(A));
        A.n  = oSolver.nCellsTotal;
        A.nnz= nnz;
        
        A.ma = (double *)itsol_malloc(nnz * sizeof(double), "read_coo:3");
        A.ja = (int *)itsol_malloc(nnz * sizeof(int), "read_coo:4");
        A.ia = (int *)itsol_malloc(nnz * sizeof(int), "read_coo:5");

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
       
        //A.ia = oSolver.Amat.ia();
        //A.ja = oSolver.Amat.ja();
        
        itsol_solver_initialize(&s, ITS_SOLVER_BICGSTAB, ITS_PC_ILUK, &A);
        s.pars.tol=1E-15;

        itsol_solver_solve(&s, oSolver.solution.begin(), oSolver.source.begin());

        Info << " (ILU) BiCGStab Solution converged to " << s.res << " in " << s.nits << " iterations " << endl;

        itsol_solver_finalize(&s);
        itsol_cleanCOO(&A);

    }
