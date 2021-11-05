/*-------------------------------------------------------------------------------*\
 A FLUX CORRECTION APPROACH FOR THE POISSON EQUATION IN INCOMPRESSIBLE FLOWS 
 ON OVERSET MESHES IN OpenFOAM
 --------------------------------------------------------------------------------
 Copyright Dominic D. J Chandar 2021
 --------------------------------------------------------------------------------
                    
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

        
  Dominic Chandar
  Queen's University of Belfast
  School of Mechanical and Aerospace Engineering
  
  d.chandar@qub.ac.uk, dominic.chandar@gmail.com
  https://github.com/domi0002/geoOverlap
  
  For issues refer to the above GitHub page.
	
License

    This code is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, version 3, i.e. GNU GPL v3

    This code are distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    For a copy of the GNU General Public License, please visit
    <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "oversetUtils.H"


int main(int argc, char *argv[])
{
    argList::addNote
    (
        "simpleFoam solver on an overset mesh with conservative correction"
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
    
    argList::addOption
    (
        "ointerface",
        "label",
        "Set the overset interface for conservation"
    );

    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMeshes.H"
    #include "createControl.H"
    #include "createFields.H"
    #include "initContinuityErrs.H"
    #include "initOverset.H"
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    
    
    Info<< "\nStarting time loop\n" << endl;

    while(runTime.loop())
    {
        
        Info<< "Time = " << runTime.timeName() << nl << endl;
        
	{        
		#include "UEqn.H"	
    		#include "pEqn.H"
        }
   
	// Note:: Turbulence fields do not yet take part in interpolation in this approach
	forAll(componentMesh, i)
	{
	    laminarTransport[i].correct();
	    turbulenceComposite[i].correct();
	}
    
        runTime.write();
	runTime.printExecutionTime(Info);

	if (convergenceSatisfied)
	{
	    runTime.writeAndEnd();
	    break;
	}
    }
    
    
    Info << "End\n" <<endl;
    
    return 0;
}


// ************************************************************************* //
