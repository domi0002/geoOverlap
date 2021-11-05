

// Composite solution to the variable 'p' (inclusive of all meshes)
PtrList<volScalarField> pComposite(componentMesh.size());

// Composite solution to the variable 'U' (inclusive of all meshes)
PtrList<volVectorField> uComposite(componentMesh.size());

// Composite solution to phi (U.nds)
PtrList<surfaceScalarField> phiComposite(nMeshes);

// Flag to determine whether the cell is a mandatory interpolation cell or not.
PtrList<volScalarField> cellManager(componentMesh.size());

// Flag to determine whether the cell is a mandatory interpolation cell or not.
PtrList<volScalarField> cellManagerTmp(componentMesh.size());


// Flag to determine whether the cell is a solution cell(1) or interpolation cell(0)
PtrList<volScalarField> cellBlank(componentMesh.size());


// Transport
PtrList<singlePhaseTransportModel> laminarTransport(componentMesh.size());

// Turbulence
PtrList<incompressible::turbulenceModel> turbulenceComposite(componentMesh.size());

// Misc.
PtrList<IOMRFZoneList> MRFComposite(componentMesh.size());

PtrList<fv::options> fvOptionsComposite(componentMesh.size());


// Ref Cell and Ref Value also in a List
// Note: While each mesh has its own ref/value, only one mesh  will use this value

List<label> pRefCell(componentMesh.size());

List<scalar> pRefValue(componentMesh.size());

//label  pRefCell  = 226;
//scalar pRefValue = 0.0;


Info<< "Reading all fields\n" << endl;

forAll(componentMesh,i)
{
    // Assign pressure
    pComposite.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "p",
                runTime.timeName(),
                componentMesh[i],
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            componentMesh[i]
        )
    );
   
    // Assign velocity
    uComposite.set
    (
        i,
        new volVectorField
        (
            IOobject
            (
                "U",
                runTime.timeName(),
                componentMesh[i],
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            componentMesh[i]
        )
    );

    // Assign phi
    phiComposite.set
    (
        i,
        new surfaceScalarField
        (
            IOobject
            (
                "phi",
                runTime.timeName(),
                componentMesh[i],
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            fvc::flux(uComposite[i])
        )
    );

    wordList boundaryTypes
    (
                componentMesh[i].boundary().size(),
                zeroGradientFvPatchScalarField::typeName
    );

    // Assign cellManager as a copy of pressure and rewrite it.
    cellManager.set(i, new volScalarField("cellManager",scalar(0)*pComposite[i]));
    cellManager[i]=0.0;

    cellManagerTmp.set(i, new volScalarField("cellManagerTmp",scalar(0)*pComposite[i]));
    cellManagerTmp[i]=0.0;

    // Assign cellBlank in the same way
    //cellBlank.set(i, new volScalarField("cellBlank",scalar(0)*pComposite[i]));
    //cellBlank[i]=1.0;

    cellBlank.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "cellBlank",
                runTime.timeName(),
                componentMesh[i],
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            componentMesh[i],
	    dimensionedScalar("FOAM_IBLANK", dimensionSet(0,0,0,0,0,0,0), 1.0),
            boundaryTypes
        )
    );

    componentMesh[i].setFluxRequired(pComposite[i].name());

    // Assign transport
    laminarTransport.set
    (
	i,
	new singlePhaseTransportModel(uComposite[i],phiComposite[i])
    );

    //Assign turbulence
    turbulenceComposite.set
    (
	i,
	incompressible::turbulenceModel::New
	(
	    uComposite[i],
	    phiComposite[i],
	    laminarTransport[i]
	).ptr()
    );

    MRFComposite.set
    (
        i,
        new IOMRFZoneList(componentMesh[i])
    );

    fvOptionsComposite.set
    (
        i,
        new fv::options(componentMesh[i])
    );

    turbulenceComposite[i].validate();


    
    setRefCell(pComposite[i], simple.dict(), pRefCell[i], pRefValue[i]);
    mesh.setFluxRequired(pComposite[i].name());

    Info << "Mesh " << i << " : pRef Cell = " << pRefCell[i] <<"  pRefValue = " << pRefValue[i] << endl;


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

bool convergenceSatisfied=false;





