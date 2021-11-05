if ( args.found("cons" ) )
    {
	surfaceScalarField & rA = rAUfComposite[i];

        // overwrite coefficients for conservation only on the alpha mesh
        if (i == oSolver.alphaMesh)
        {
	    label lRow=-1;

	    // If user specified overset cell, return the interface else specify the interface and return the cell.
	    if (oSolver.ocell > -1 )
	    {
	 	lRow=oSolver.ocell;
		oSolver.oint = oSolver.connect[i]->locateSingleOversetCellInterface(oSolver.cellBlank[i],lRow);
	    }
	    else
	    {
            	lRow   = oSolver.connect[i]->locateSingleOversetInterfaceCell(oSolver.cellBlank[i],oSolver.oint);
 	    }

	    // +++ For Checking only Old Way+++
            // label lRow = 4673;
	    // oSolver.oint = oSolver.connect[i]->locateSingleOversetCellInterface(oSolver.cellBlank[i],lRow);
	    // +++
	    

	    assert(lRow > -1 );
            label gRow = lRow + cumulativeSum[i];//oSolver.connect[i]->locateSingleOversetInterfaceCell(oSolver.cellBlank[i],oSolver.oint);

	    // The alpha mesh is the mesh where an overset interface/cell is found to enforce conservation
            Info << "alpha Mesh = " << oSolver.alphaMesh << "  alpha cell = " << lRow << " global cell = " << gRow << " Interface = " << oSolver.oint<< endl;
           
            for ( label colIndex = 0; colIndex < oSolver.nCellsTotal ; colIndex++ )
            {
                    Amat.set(gRow,colIndex,0);
            }
            
            source[gRow] = 0;
           
            // All interfaces on the current mesh
            forAll(oSolver.connect[i]->oversetInterface,of_)
            {
                label faceID = oSolver.connect[i]->oversetInterface[of_];
                label own    = componentMesh[i].owner()[faceID];
                label nei    = componentMesh[i].neighbour()[faceID];
            
                if ( of_ == oSolver.oint )
                {
                    // Owner is c0. So flip sign of nei and gRow (usually nei has +ve coeff on the snGrad) - refer notes
                    // Later change mesh to componentMesh[alphaMesh]
                    if (own == lRow)//gRow
                    {
                        Amat.set(gRow,nei,  Amat.get(gRow,nei) +  -deltaCoeffs[faceID]*rA[faceID]*mesh.magSf()[faceID]);
                        Amat.set(gRow,gRow, Amat.get(gRow,gRow)+   deltaCoeffs[faceID]*rA[faceID]*mesh.magSf()[faceID]);
             
                    }
                    else
                    {
                        Amat.set(gRow,own,  Amat.get(gRow,own) +  -deltaCoeffs[faceID]*rA[faceID]*mesh.magSf()[faceID]);
                        Amat.set(gRow,gRow, Amat.get(gRow,gRow)+  deltaCoeffs[faceID]*rA[faceID]*mesh.magSf()[faceID]);//this is the nei-so no flipping
      
                    }
                    
                    // Owner is c0 (flip)-refer notes
                    if ( own == lRow )//gRow
                        source[gRow] -= phiHbyAComposite[i][faceID];
                    else
                        source[gRow] += phiHbyAComposite[i][faceID];
                    
                }
                else
                {
                    if ( cellBlank[i][own] == 0 )
                    {
                        Amat.set(gRow,nei, Amat.get(gRow,nei) + -deltaCoeffs[faceID]*rA[faceID]*mesh.magSf()[faceID]);
                        Amat.set(gRow,own, Amat.get(gRow,own) +  deltaCoeffs[faceID]*rA[faceID]*mesh.magSf()[faceID]);
                        source[gRow] -= phiHbyAComposite[i][faceID];
                    }
                    else
                    {
                        Amat.set(gRow,nei, Amat.get(gRow,nei)  +  deltaCoeffs[faceID]*rA[faceID]*mesh.magSf()[faceID]);
                        Amat.set(gRow,own, Amat.get(gRow,own)  + -deltaCoeffs[faceID]*rA[faceID]*mesh.magSf()[faceID]);
                        source[gRow] += phiHbyAComposite[i][faceID];
                    }
                     
                }
                
            
            }
            
            // All interfaces on other meshes
            for (int jMesh = 0 ; jMesh < nMeshes ; jMesh++ )
            {
                
                if (jMesh!=oSolver.alphaMesh)
                {
                    const fvMesh & meshj 	    	    = componentMesh[jMesh];
                    const surfaceScalarField & deltaCoeffsj = componentMesh[jMesh].surfaceInterpolation::deltaCoeffs();
		    surfaceScalarField & rAJ                = rAUfComposite[jMesh];    

                    for (int oi=0; oi < oSolver.connect[jMesh]->oversetInterface.size() ; oi++ )
                    {
                        label faceID = oSolver.connect[jMesh]->oversetInterface[oi];
                        label own = componentMesh[jMesh].owner()[faceID];
                        label nei = componentMesh[jMesh].neighbour()[faceID];
                        
                        label ownCol = own + cumulativeSum[i==0?1:0];
                        label neiCol = nei + cumulativeSum[i==0?1:0];
                        
                        if (oSolver.cellBlank[jMesh][own] == 1)
                        {
                            Amat.set(gRow,neiCol, Amat.get(gRow,neiCol) +    deltaCoeffsj[faceID]*rAJ[faceID]*meshj.magSf()[faceID]);
                            Amat.set(gRow,ownCol, Amat.get(gRow,ownCol) +   -deltaCoeffsj[faceID]*rAJ[faceID]*meshj.magSf()[faceID]);
                        
                            source[gRow] += phiHbyAComposite[jMesh][faceID];
                        }
                        else
                        {
                            Amat.set(gRow,neiCol, Amat.get(gRow,neiCol) +  -deltaCoeffsj[faceID]*rAJ[faceID]*meshj.magSf()[faceID]);
                            Amat.set(gRow,ownCol, Amat.get(gRow,ownCol) +   deltaCoeffsj[faceID]*rAJ[faceID]*meshj.magSf()[faceID]);
                        
                            source[gRow] -= phiHbyAComposite[jMesh][faceID];
                        }
                    }
                }
                    
            }
            
            
        }// end i==alphamesh
    }  // end if conservation needs to be enforced
