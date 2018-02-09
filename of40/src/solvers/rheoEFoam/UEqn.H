   // Momentum predictor

    tmp<fvVectorMatrix> tUEqn
    (
        fvm::ddt(U) 
      + fvm::div(phi, U)
      ==
        fvOptions(U)  
      + constEq.divTau(U)  
      + elecM.Fe()/constEq.rho()   
    );

    fvVectorMatrix& UEqn = tUEqn.ref();

    UEqn.relax();

    fvOptions.constrain(UEqn);

    solve(UEqn == -fvc::grad(p));

    fvOptions.correct(U);

