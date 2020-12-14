/*---------------------------------------------------------------------------*\
 A CONSISTENT FLUX CORRECTION APPROACH FOR THE POISSON EQUATION
 ON OVERSET MESHES
                    
  The base code has the following dependencies
  
  1. OpenFOAM-v1812+ (ESI-OpenCFD) 
        (www.openfoam.com)
  
  2. Octree search for overset 
        (https://github.com/jbehley/octree):
            @conference{behley2015icra,
            author = {Jens Behley and Volker Steinhage and Armin B. Cremers},
             title = {{Efficient Radius Neighbor Seach in Three-dimensional Point Clouds}},
         booktitle = {Proc. of the IEEE International Conference on Robotics and Automation (ICRA)},
            year = {2015}
            
  3. The sparse matrix representation classes 
        (https://github.com/uestla/Sparse-Matrix)
  
  4. Saad's Sparse Matrix solvers (BiCGStab) 
        (https://www-users.cs.umn.edu/~saad/software/ITSOL)

        
  Dominic D. J Chandar
  Queen's University of Belfast
  School of Mechanical and Aerospace Engineering
  
  d.chandar@qub.ac.uk, dominic.chandar@gmail.com
  https://github.com/domi0002/geoOverlap
  
  For issues refer to the above GitHub page.

  Updated: July 2020
  
\*---------------------------------------------------------------------------*/

#include "twoDOverset.H"


int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Poisson equation solver on an overset mesh"
    );
    
    argList::addBoolOption
    (
        "cons",
        "Conservative correction"
    );
    
    argList::addBoolOption
    (
        "inv",
        "Inverse Distance Interpolation"
    );
    
	argList::addBoolOption
    (
        "rbf",
        "Inverse Distance Interpolation"
    );
	
    argList::addBoolOption
    (
        "neumannrhs",
        "RHS consistent with Neumann BCs"
    );

    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMeshes.H"
    #include "createFields.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    
    // Create the new overset solver
    overlapSolver oSolver(componentMesh,cellManager,cellBlank,0);
    
    // Create file for dummy forces
    OFstream forceOutput("forces.dat");
    
    while(runTime.loop())
    {
        
        Info<< "Time = " << runTime.timeName() << nl << endl;
        
        // Update moving mesh (moving regions activated in constant/region<x>/dynamicMeshDict)
        oSolver.updateMovingMesh();
        
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
   
        // Create the matrix and the associated source for the linear solver
        #include "createMatrixAndRHS.H"
    
        // Solve using Saad's library
        #include "itsSolver.H"
 
        // Post-process :: write solution and compute errors
        #include "postProcessOverset.H"    
        
        runTime.write();
    }
    
    
   
    
    return 0;
}


// ************************************************************************* //
