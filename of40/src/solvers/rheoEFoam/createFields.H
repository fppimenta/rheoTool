
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
   setRefCell(p, simple.dict(), pRefCell, pRefValue);
   mesh.setFluxRequired(p.name());

    // Create constitutive equation
    constitutiveModel constEq(U, phi);
    
    // Create the electric model    
    EDFModel elecM(phi);

    // Read C field if present, otherwise the field is (always) initialized as 0
    volScalarField C
    (
        IOobject
        (
            "C",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar
        (
                "0",
                dimless,
                pTraits<scalar>::zero
        ),
        extrapolatedCalculatedFvPatchField<scalar>::typeName
    );
  
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
    
