
scalar mass=0;

for (label i = 0 ; i < nMeshes ; i++)
{
    const fvMesh & mesh = componentMesh[i];
        
    forAll(mesh.boundary(),patchI)
    {
        if ( mesh.boundaryMesh()[patchI].physicalType()!="Overset")
        {
            mass += sum(phiComposite[i].boundaryFieldRef()[patchI]);
        }
    }
        
}

Info << " Total Mass flux across external boundaries --->>>" << mass << endl;
