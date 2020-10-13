
// Error Initialize
label m=0;
scalar l2Error=0;
    
label nInterpolationCellsTotal=0;
    
for(label i = 0 ; i < nMeshes ; i++)
{
    
    if (args.found("neumannrhs"))
    {
        forAll(TComposite[i],j)
        {
            TComposite[i].primitiveFieldRef()[j] = oSolver.solution[m];
            TExact[i].primitiveFieldRef()[j]=Foam::cos(kconst*M_PI*componentMesh[i].C()[j].x())*Foam::sin(kconst*M_PI*componentMesh[i].C()[j].y());
            scalar localError=cellBlank[i][j]*(TExact[i].primitiveFieldRef()[j]-TComposite[i].primitiveFieldRef()[j]);
            TError[i].primitiveFieldRef()[j]=localError;
            
            l2Error   += (localError*localError);
            m++;
        } 
                
    }
    else
    {
        forAll(TComposite[i],j)
        {
            TComposite[i].primitiveFieldRef()[j] = oSolver.solution[m];
            TExact[i].primitiveFieldRef()[j]=Foam::sin(kconst*M_PI*componentMesh[i].C()[j].x())*Foam::sin(kconst*M_PI*componentMesh[i].C()[j].y());
            scalar localError=cellBlank[i][j]*(TExact[i].primitiveFieldRef()[j]-TComposite[i].primitiveFieldRef()[j]);
            TError[i].primitiveFieldRef()[j]=localError;
            
            l2Error   += (localError*localError);
            m++;
        }
                
    }
            
    nInterpolationCellsTotal+=oSolver.connect[i]->nInterpolationCells;
        
    TError[i].write();
    TExact[i].write();
    cellBlank[i].write();
}
    
l2Error = Foam::sqrt(l2Error/(oSolver.nCellsTotal-nInterpolationCellsTotal));
    
// 1. Check for Conservation Error(unbalanced sum): (sum(gradp)-sum(rhs)) on all boundaries except mandatory overset boundaries
// 2. Also Lets compute a force based on Sum(TndS) along the "upper" boundary as an analog to pressure.
// 3. Compute the divergence free velocity field q* = q - grad(p) on the boundaries and sum it up as 'mass'

scalar unbalancedSum=0;
scalar mass=0;
vector forcesExact(0,0,0), forces(0,0,0);
for (label i = 0 ; i < nMeshes ; i++)
{
    const dynamicFvMesh & mesh = componentMesh[i];
        
    surfaceScalarField gradTDotS(fvc::snGrad(TComposite[i])*mesh.magSf());
    //const surfaceScalarField & deltaCoeffs = mesh.surfaceInterpolation::deltaCoeffs();    
    //const surfaceScalarField::Boundary & deltaBound = deltaCoeffs.boundaryField();
    
    volVectorField gradT("gradT",fvc::grad(TComposite[i]));
    const volVectorField::Boundary & gradTBoundary = gradT.boundaryFieldRef();

    forAll(mesh.boundary(),patchI)
    {
        if ( mesh.boundaryMesh()[patchI].physicalType()!="Overset")
        {
            unbalancedSum+= sum(gradTDotS.boundaryFieldRef()[patchI]-phiComposite[i].boundaryFieldRef()[patchI]);
            forAll(mesh.boundary()[patchI],faceI)
            {		    
            	mass+=phiComposite[i].boundaryFieldRef()[patchI][faceI] 
                     - (gradTBoundary[patchI][faceI] & mesh.Sf().boundaryField()[patchI][faceI]);
                    
                    /*(TComposite[i].boundaryFieldRef()[patchI][faceI]
                    -TComposite[i].primitiveFieldRef()[mesh.boundary()[patchI].faceCells()[faceI]])
                    *deltaBound[patchI][faceI];
                    */
            }

            // Sum T.ndS on the upper boundary
            if (mesh.boundaryMesh()[patchI].name()=="upper")
            {
                forAll(mesh.boundary()[patchI],fI)
                {
                    forces+=TComposite[i].primitiveFieldRef()[mesh.boundary()[patchI].faceCells()[fI]]*mesh.boundary()[patchI].Sf()[fI];
                    forcesExact+=TExact[i].primitiveFieldRef()[mesh.boundary()[patchI].faceCells()[fI]]*mesh.boundary()[patchI].Sf()[fI];
                    
                }
            }
        }
    }
        
}
    
forceOutput << runTime.timeName() <<"\t"<<mag(forcesExact-forces)<<"\t"<<nInterpolationCellsTotal<<endl;    
    
Info << "Mesh.size = " << oSolver.nCellsTotal << nl
    <<" L2 error  = " << l2Error  << nl
    <<" C error   = " << unbalancedSum<<nl
    <<" Actual F  = " << forcesExact <<nl
    <<" Compu  F  = " << forces << nl
    <<" Force Err = " << mag(forcesExact-forces)<<nl
    <<" mass Err  = " << mass << endl;
