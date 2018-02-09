
    Info<< "Reading field p\n" << endl;
    volScalarField p
    (
        IOobject
        (
            "p",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

     
    Info<< "Reading field U\n" << endl;
    volVectorField U
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );


    #include "createPhi.H"


    label pRefCell = 0;
    scalar pRefValue = 0.0;
    setRefCell(p, mesh.solutionDict().subDict("SIMPLE"), pRefCell, pRefValue);

    // Create constitutive equation
    
    constitutiveModel constEq(U, phi);
    
    IOdictionary cttProperties
    (
        IOobject
        (
            "constitutiveProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    );

    int    nInIter = mesh.solutionDict().subDict("SIMPLE").lookupOrDefault<int>("nInIter", 1);
    bool   sPS = cttProperties.subDict("passiveScalarProperties").lookupOrDefault<Switch>("solvePassiveScalar", false); 

    // Read C field if present, otherwise the field is (always) initialized as 0
    volScalarField C
    (
        IOobject
        (
            "C",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            sPS == false ? (IOobject::NO_WRITE) : (IOobject::AUTO_WRITE)
        ),
        mesh,
        dimensionedScalar
        (
                "0",
                dimless,
                pTraits<scalar>::zero
        ),
        zeroGradientFvPatchField<scalar>::typeName
    );
