Info<< "Reading field T\n" << endl;

// Composite solution to the variable 'T' (inclusive of all meshes)
PtrList<volScalarField> TComposite(componentMesh.size());

// Flag to determine whether the cell is a mandatory interpolation cell or not.
PtrList<volScalarField> cellManager(componentMesh.size());

// Exact solution to the Poisson equation
PtrList<volScalarField> TExact(componentMesh.size());

// Flag to determine whether the cell is a solution cell(1) or interpolation cell(0)
PtrList<volScalarField> cellBlank(componentMesh.size());

// Error between exact and numerical solution
PtrList<volScalarField> TError(componentMesh.size());

forAll(componentMesh,i)
{

    TComposite.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "T",
                runTime.timeName(),
                componentMesh[i],
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            componentMesh[i]
        )
    );
    TComposite[i]=0.0;
    
    cellManager.set(i, new volScalarField("cellManager",scalar(0)*TComposite[i]));
    cellManager[i]=0.0;
    TExact.set(i, new volScalarField("TExact",scalar(0)*TComposite[i]));
    TError.set(i, new volScalarField("TError",scalar(0)*TComposite[i]));
    cellBlank.set(i, new volScalarField("cellBlank",scalar(0)*TComposite[i]));
    cellBlank[i]=1.0;
}

Info<< "Reading transportProperties\n" << endl;

IOdictionary transportProperties
(
    IOobject
    (
        "transportProperties",
        runTime.constant(),
        componentMesh[0],
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);


Info<< "Reading diffusivity DT\n" << endl;

dimensionedScalar DT("DT", dimViscosity, transportProperties);


