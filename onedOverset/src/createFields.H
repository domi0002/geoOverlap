Info<< "Reading field T\n" << endl;

PtrList<volScalarField> TComposite(componentMesh.size());
PtrList<volScalarField> cellManager(componentMesh.size());
PtrList<volScalarField> TExact(componentMesh.size());
PtrList<volScalarField> cellBlank(componentMesh.size());
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


