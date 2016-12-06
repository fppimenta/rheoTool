   // Momentum predictor

    tmp<fvVectorMatrix> UEqn
    (
        fvm::ddt(U) 
      + fvm::div(phi, U)
      ==
        fvOptions(U)  
      + constEq.divTau(U)    
    );

    UEqn().relax();

    fvOptions.constrain(UEqn());

    solve(UEqn() == -fvc::grad(p));

    fvOptions.correct(U);

