   // Transport of passive scalar

    dimensionedScalar  D_ = cttProperties.subDict("passiveScalarProperties").lookup("D");
     
    fvScalarMatrix CEqn
    (
        fvm::ddt(C()) 
      + fvm::div(phi, C())
      ==
        fvm::laplacian(D_, C())    
      + fvModels.source(C())
    );

    CEqn.relax();
    fvConstraints.constrain(CEqn);
    CEqn.solve();
    fvConstraints.constrain(C());
    
     
 
