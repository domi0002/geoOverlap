
regionProperties rp(runTime);

const wordList fluidNames(rp["componentMeshes"]);
  
PtrList<fvMesh> componentMesh(fluidNames.size());

    forAll(fluidNames, i)
    {
        Info<< "Create mesh for region " << fluidNames[i]
            << " for time = " << runTime.timeName() << nl << endl;

        componentMesh.set
        (
            i,
            new fvMesh
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
