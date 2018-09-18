   // Transport of passive scalar

    dimensionedScalar  D_ = cttProperties.subDict("passiveScalarProperties").lookup("D");
     
    fvScalarMatrix CEqn
    (
        fvm::ddt(C()) 
      + fvm::div(phi, C())
      ==
        fvc::laplacian(D_, C())    
    );

    CEqn.relax();
    CEqn.solve();
 
