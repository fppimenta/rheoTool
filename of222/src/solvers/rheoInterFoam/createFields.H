
    Info<< "Reading field p_rgh\n" << endl;
    volScalarField p_rgh
    (
        IOobject
        (
            "p_rgh",
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

 

    // Create the mixture
    
    constitutiveTwoPhaseMixture twoPhaseProperties(U, phi);
   
    volScalarField& alpha1(twoPhaseProperties.alpha1());
    volScalarField& alpha2(twoPhaseProperties.alpha2());

    const dimensionedScalar& rho1 = twoPhaseProperties.rho1();
    const dimensionedScalar& rho2 = twoPhaseProperties.rho2();


    // Need to store rho for ddt(rho, U)
    volScalarField rho
    (
        IOobject
        (
            "rho",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT
        ),
        alpha1*rho1 + alpha2*rho2,
        alpha1.boundaryField().types()
    );
    rho.oldTime();


    // Mass flux
    // Initialisation does not matter because rhoPhi is reset after the
    // alpha1 solution before it is used in the U equation.
    surfaceScalarField rhoPhi
    (
        IOobject
        (
            "rho*phi",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        rho1*phi
    );


    // Construct interface from alpha1 distribution
    interfaceProperties interface(alpha1, U, twoPhaseProperties);

    #include "readGravitationalAcceleration.H"

    Info<< "Calculating field g.h\n" << endl;
    volScalarField gh("gh", g & mesh.C());
    surfaceScalarField ghf("ghf", g & mesh.Cf());

    volScalarField p
    (
        IOobject
        (
            "p",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        p_rgh + rho*gh
    );

    label pRefCell = 0;
    scalar pRefValue = 0.0;
    setRefCell
    (
        p,
        p_rgh,
        mesh.solutionDict().subDict("PIMPLE"),
        pRefCell,
        pRefValue
    );

    if (p_rgh.needReference())
    {
        p += dimensionedScalar
        (
            "p",
            p.dimensions(),
            pRefValue - getRefCellValue(p, pRefCell)
        );
        p_rgh = p - rho*gh;
    }

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
 
    int    nInIter = mesh.solutionDict().subDict("PIMPLE").lookupOrDefault<int>("nInIter", 1);
    bool   simplec = mesh.solutionDict().subDict("PIMPLE").lookupOrDefault<Switch>("SIMPLEC", true);
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
    
     
