    // Set the overset interface if enforcing conservation
    label oint_=0;
    bool interfaceSpecified=false;    

    if (args.optionReadIfPresent("ointerface",oint_))
    {
        Info << "Using overset interface = "<<oint_<<" for conservation" << endl;
	interfaceSpecified=true;
    }


    label ocell_=-1;
    if (args.optionReadIfPresent("ocell",ocell_))
    {
	if (interfaceSpecified)
	{
	    Info <<" Cannot specify both interface and acceptor cell for conservation, Falling back to interface" <<  endl;
	    ocell_=-1;
	}
	else
	{
            Info << "Using acceptor(interpolation) cell = "<<ocell_<<" for conservation" << endl;
	}
    }

    
    // Create the new overset solver
    // Note that '0' in the last but one entry - That's the mesh ID where conservation takes place.
    overlapSolver oSolver(componentMesh,cellManager,cellManagerTmp,cellBlank,0,oint_,ocell_);

    // Hardcode for now (I need 9 donor cells for a given receiver in 2D)
    if ( args.found("rbf") )
    {
	oSolver.RBFPoints = 9;
    }
    else
    {
	oSolver.RBFPoints = 0;
    }
		
    // Calculate overlap (needs to be called again for moving meshes)
    oSolver.preProcessAndCompute();

    // Calculate interpolation weights
    if ( args.found("inv") )
    {
        oSolver.distanceWeightedWeights();
    }
    else if ( args.found("rbf") )
    {
        oSolver.rbfWeights();
    }
    else
    {
        oSolver.polynomialWeights();
    }

