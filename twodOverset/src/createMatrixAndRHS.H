PtrList<volScalarField> RHSComposite(nMeshes);
PtrList<fvScalarMatrix> TEqnComposite(nMeshes);
PtrList<surfaceScalarField> phiComposite(nMeshes);
    
scalar kconst;

// If running Neumann BCs, we use 3.5 periods so that the upper and left boundaries automatically satisfy zero gradients
if ( args.found("neumannrhs"))
{
    kconst=3.5;
}
else
{
    kconst=3.0;
}

for( label mi = 0 ; mi < nMeshes ; mi++)
    {
        const dynamicFvMesh & mesh = componentMesh[mi];
        
        // Define the RHS (Note: always use scalar(1) here since mesh,Cf() is taken as a reference here!)
        // If you fail to use scalar(1) and modify rhsFace, then mesh.Cf() is also altered.
        surfaceVectorField rhsFace("cf",mesh.Cf()*scalar(1));
    
        if ( args.found("neumannrhs"))
        {    
            // Internal faces
            forAll( mesh.faces(), i)
            {
                if  (mesh.isInternalFace(i))
                {
                    rhsFace[i].x()=-kconst*M_PI*Foam::sin(kconst*M_PI*mesh.Cf()[i].x())*Foam::sin(kconst*M_PI*mesh.Cf()[i].y());
                    rhsFace[i].y()= kconst*M_PI*Foam::cos(kconst*M_PI*mesh.Cf()[i].x())*Foam::cos(kconst*M_PI*mesh.Cf()[i].y());
                    rhsFace[i].z()=0;
                }
            }
        }
        else
        {
            // Internal faces
            forAll( mesh.faces(), i)
            {
                if  (mesh.isInternalFace(i))
                {
                    rhsFace[i].x()=kconst*M_PI*Foam::cos(kconst*M_PI*mesh.Cf()[i].x())*Foam::sin(kconst*M_PI*mesh.Cf()[i].y());
                    rhsFace[i].y()=kconst*M_PI*Foam::cos(kconst*M_PI*mesh.Cf()[i].y())*Foam::sin(kconst*M_PI*mesh.Cf()[i].x());
                    rhsFace[i].z()=0;
                }
            }
            
        }
    
        surfaceVectorField::Boundary & brhs = rhsFace.boundaryFieldRef(); 
    
        // Boundary faces
        forAll(mesh.boundary(),patchI)
        {
        
            fvsPatchVectorField  & bTestVec = brhs[patchI];
            vectorField assign(bTestVec.size());
        
            const fvPatch & pp = bTestVec.patch();
            vectorField fc(pp.Cf());
        
            if ( args.found("neumannrhs"))
            {
                forAll(mesh.boundary()[patchI],faceI)
                {
                    assign[faceI].x()=-kconst*M_PI*Foam::sin(kconst*M_PI*fc[faceI].x())*Foam::sin(kconst*M_PI*fc[faceI].y());
                    assign[faceI].y()= kconst*M_PI*Foam::cos(kconst*M_PI*fc[faceI].x())*Foam::cos(kconst*M_PI*fc[faceI].y());
                    assign[faceI].z()= 0;
                }
            }
            else
            {
                forAll(mesh.boundary()[patchI],faceI)
                {
                    assign[faceI].x()=kconst*M_PI*Foam::cos(kconst*M_PI*fc[faceI].x())*Foam::sin(kconst*M_PI*fc[faceI].y());
                    assign[faceI].y()=kconst*M_PI*Foam::cos(kconst*M_PI*fc[faceI].y())*Foam::sin(kconst*M_PI*fc[faceI].x());
                    assign[faceI].z()=0;
                }
                
            }
        
            bTestVec==assign;
        }
    

        surfaceScalarField phi(rhsFace & mesh.Sf() );
    
        // This is the RHS of the Laplacian
        RHSComposite.set(mi, new volScalarField("RHS",fvc::div(phi)));
        
        // Frame the laplacian everwhere based on OpenFOAM
        TEqnComposite.set(mi, new fvScalarMatrix(fvm::laplacian(DT,TComposite[mi])));
        
        phiComposite.set(mi, new surfaceScalarField("phi",phi));
        
    }
    
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
    
    SM::SparseMatrix<double> & Amat = oSolver.Amat;
    scalarField & source            = oSolver.source;
    
    
    for(label i = 0 ; i < nMeshes ; i++)
    {
        const dynamicFvMesh & mesh = componentMesh[i];
        const fvScalarMatrix & TEqn = TEqnComposite[i];
        const surfaceScalarField & deltaCoeffs = mesh.surfaceInterpolation::deltaCoeffs();
    
        // Set the diagonal and the RHS
        forAll(TComposite[i],rowI)
        {
            label gRow = rowI + cumulativeSum[i];
            source[gRow] = mesh.V()[rowI]*RHSComposite[i][rowI];
            Amat.set(gRow,gRow,TEqn.diag()[rowI]);
        }
        
        
        // Upper and Lower Coefficients (i.e. off-diagonal coefficients)
        for(label faceI= 0 ; faceI < mesh.nInternalFaces() ; faceI++)
        {
            label nei = TEqn.lduAddr().lowerAddr()[faceI] + cumulativeSum[i];
            label own = TEqn.lduAddr().upperAddr()[faceI] + cumulativeSum[i];
        
            Amat.set(own,nei,TEqn.upper()[faceI]);
            Amat.set(nei,own,TEqn.upper()[faceI]);
        }
        
        
        // Boundary Conditions
        forAll(TComposite[i].boundaryFieldRef(),patchI)
        {
            const fvPatch & pp = TComposite[i].boundaryFieldRef()[patchI].patch();
            forAll(pp,faceI)
            {
                label cellI = pp.faceCells()[faceI] + cumulativeSum[i];
                
                scalar existingInternalValue = Amat.get(cellI,cellI);
                scalar newInternalValue = existingInternalValue + TEqn.internalCoeffs()[patchI][faceI];
                
                Amat.set(cellI,cellI,newInternalValue);
                source[cellI]+= TEqn.boundaryCoeffs()[patchI][faceI];
            }
        }
        
        
        // Interpolation Coefficients (Overset)
        forAll(TComposite[i],rowI)
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
        
    // --------------------------------------------------------
    // The consistent flux correction approach by D.Chandar
    // --------------------------------------------------------
    if ( args.found("cons" ) )
    {
        // overwrite coefficients for conservation only on the alpha mesh
        if (i == oSolver.alphaMesh)
        {
            label gRow = oSolver.connect[i]->locateSingleOversetInterfaceCell(oSolver.cellBlank[i]);
            Info << "alpha Mesh = " << oSolver.alphaMesh << "  alpha cell = " << gRow << endl;
           
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
            
                if ( of_ == 0 )
                {
                    // Owner is c0. So flip sign of nei and gRow (usually nei has +ve coeff on the snGrad) - refer notes
                    // Later change mesh to componentMesh[alphaMesh]
                    if (own == gRow)
                    {
                        Amat.set(gRow,nei,  Amat.get(gRow,nei) +  -deltaCoeffs[faceID]*mesh.magSf()[faceID]);
                        Amat.set(gRow,gRow, Amat.get(gRow,gRow)+   deltaCoeffs[faceID]*mesh.magSf()[faceID]);
             
                    }
                    else
                    {
                        Amat.set(gRow,own,  Amat.get(gRow,own) +  -deltaCoeffs[faceID]*mesh.magSf()[faceID]);
                        Amat.set(gRow,gRow, Amat.get(gRow,gRow)+  deltaCoeffs[faceID]*mesh.magSf()[faceID]);//this is the nei-so no flipping
      
                    }
                    
                    // Owner is c0 (flip)-refer notes
                    if ( own == gRow )
                        source[gRow] -= phiComposite[i][faceID];
                    else
                        source[gRow] += phiComposite[i][faceID];
                    
                }
                else
                {
                    if ( cellBlank[i][own] == 0 )
                    {
                        Amat.set(gRow,nei, Amat.get(gRow,nei) + -deltaCoeffs[faceID]*mesh.magSf()[faceID]);
                        Amat.set(gRow,own, Amat.get(gRow,own) +  deltaCoeffs[faceID]*mesh.magSf()[faceID]);
                        source[gRow] -= phiComposite[i][faceID];
                    }
                    else
                    {
                        Amat.set(gRow,nei, Amat.get(gRow,nei)  +  deltaCoeffs[faceID]*mesh.magSf()[faceID]);
                        Amat.set(gRow,own, Amat.get(gRow,own)  + -deltaCoeffs[faceID]*mesh.magSf()[faceID]);
                        source[gRow] += phiComposite[i][faceID];
                    }
                     
                }
                
            
            }
            
            // All interfaces on other meshes
            for (int jMesh = 0 ; jMesh < nMeshes ; jMesh++ )
            {
                
                if (jMesh!=oSolver.alphaMesh)
                {
                    const dynamicFvMesh & meshj = componentMesh[jMesh];
                    const surfaceScalarField & deltaCoeffsj = componentMesh[jMesh].surfaceInterpolation::deltaCoeffs();
                    
                    for (int oi=0; oi < oSolver.connect[jMesh]->oversetInterface.size() ; oi++ )
                    {
                        label faceID = oSolver.connect[jMesh]->oversetInterface[oi];
                        label own = componentMesh[jMesh].owner()[faceID];
                        label nei = componentMesh[jMesh].neighbour()[faceID];
                        
                        label ownCol = own + cumulativeSum[i==0?1:0];
                        label neiCol = nei + cumulativeSum[i==0?1:0];
                        
                        if (oSolver.cellBlank[jMesh][own] == 1)
                        {
                            Amat.set(gRow,neiCol, Amat.get(gRow,neiCol) +   deltaCoeffsj[faceID]*meshj.magSf()[faceID]);
                            Amat.set(gRow,ownCol, Amat.get(gRow,ownCol) +   -deltaCoeffsj[faceID]*meshj.magSf()[faceID]);
                        
                            source[gRow] += phiComposite[jMesh][faceID];
                        }
                        else
                        {
                            Amat.set(gRow,neiCol, Amat.get(gRow,neiCol) +  -deltaCoeffsj[faceID]*meshj.magSf()[faceID]);
                            Amat.set(gRow,ownCol, Amat.get(gRow,ownCol) +   deltaCoeffsj[faceID]*meshj.magSf()[faceID]);
                        
                            source[gRow] -= phiComposite[jMesh][faceID];
                        }
                    }
                }
                    
            }
            
            
        }// end i==alphamesh
    }  // end if conservation needs to be enforced
    
}// end loop over meshes

 
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
    
