
regionProperties rp(runTime);

const wordList fluidNames(rp["componentMeshes"]);

//- Dynamic Generalization
PtrList<dynamicFvMesh> componentMesh(fluidNames.size());

forAll(fluidNames, i)
{
        Info<< "Create mesh for region " << fluidNames[i]
            << " for time = " << runTime.timeName() << nl << endl;

        componentMesh.set
        (
            i,
            dynamicFvMesh::New
            (
                IOobject
                (
                    fluidNames[i],
                    runTime.timeName(),
                    runTime,
                    IOobject::MUST_READ
                )
            )
        );
}    

    
const label nMeshes = componentMesh.size();
