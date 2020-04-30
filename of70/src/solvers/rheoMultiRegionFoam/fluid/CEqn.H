   // Transport of passive scalar
   
    IOdictionary cttProperties
    (
      IOobject
      (
        "constitutiveProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE,
        false
      )
    );  

    dimensionedScalar  D_ = cttProperties.subDict("passiveScalarProperties").lookup("D");
    
    volScalarField& C = CFluid[i];
     
    fvScalarMatrix CEqn
    (
        fvm::ddt(C) 
      + fvm::div(phi, C)
      ==
        fvm::laplacian(D_, C)    
    );

    CEqn.relax();
    CEqn.solve();
 
