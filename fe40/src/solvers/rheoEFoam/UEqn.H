   // Momentum predictor
 
    tmp<fvVectorMatrix> UEqn
    (
        fvm::ddt(U) 
      + fvm::div(phi, U)
      == 
        constEq.divTau(U) 
      + elecM.Fe()/constEq.rho()    
    );

    UEqn().relax();

    solve(UEqn() == -fvc::grad(p));
 
